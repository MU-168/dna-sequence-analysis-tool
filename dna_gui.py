#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
DNA 序列分析教学工具 · Tkinter 桌面前端
=====================================================================

启动方式
--------
    py -3 dna_gui.py
    py -3 dna_sequence_analysis.py        # 无参数运行也会打开此图形界面

程序启动后会先弹出一个 DNA 输入窗口。可粘贴 DNA 序列、填入内置示例，
或从 FASTA 文件载入。点击“开始分析”后，主窗口按标签页显示：

    概览 · ORF 候选 · 剪接与翻译 · UniProt / 相似蛋白 · 原始序列

分析内核全部来自 dna_sequence_analysis.py（Biopython、六读码框 ORF、
TATAAA/AATAAA 标记、GU-AG 简化剪接、标准遗传密码翻译、UniProt REST API）。
本前端只负责输入、编排、进度显示和结果可视化；核心分析、UniProt 精确检索和
基于 GenomeNet BLASTP 的相似性回退都在 dna_sequence_analysis.py 中实现。
GenomeNet BLAST 是外部公开网页服务，仅适合低频教学查询，请勿批量提交。

依赖
----
    pip install biopython
    tkinter 为 Python 标准库自带组件，无需 pip 安装。

声明
----
本工具是教学模拟器，不是基因预测或临床/科研软件。详见内核脚本头部声明。
所有结果不应直接用于正式科研、临床诊断、药物设计或生产决策。

数据与网络
----------
- 原始 DNA 仅在本机分析，不上传。
- 蛋白序列会提交给 UniProt KB REST 进行 MD5 精确检索。
- 如果精确检索无命中且开启了相似性搜索，蛋白序列还会提交给 GenomeNet 的公开
  BLASTP 服务（Swiss-Prot 数据库），随后从 UniProtKB 获取命中条目的注释。
  相似性结果按“一致度 × 查询覆盖度”的综合匹配度排序，只是可能相关的功能线索。
"""

from __future__ import annotations

import ctypes
import queue
import sys
import threading
import time
import traceback
import webbrowser
from dataclasses import dataclass, field
from typing import Callable, Optional

import tkinter as tk
import tkinter.font as tkfont
from tkinter import filedialog, messagebox, simpledialog, ttk
from Bio import SeqIO

from dna_sequence_analysis import (
    DEMO_DNA_SEQUENCE,
    DEMO_SOURCE_LABEL,
    MAX_CANDIDATES_TO_ANALYZE,
    MIN_ORF_AA,
    POLY_A_SIGNAL_MOTIF,
    REMOVE_REPEAT_REGIONS,
    REQUIRE_TATA_BOX,
    SCREEN_TANDEM_REPEATS,
    TATA_BOX_MOTIF,
    UNIPROT_MIN_INTERVAL_SECONDS,
    UNIPROT_TIMEOUT_SECONDS,
    UNIPROT_SIMILARITY_ENABLED,
    CandidateGene,
    GeneSearchResult,
    RuntimeOptions,
    SplicingResult,
    TandemRepeatRegion,
    TranslationResult,
    UniProtClient,
    UniProtResult,
    find_tandem_repeats,
    interval_1based,
    load_fasta_sequence,
    remove_marked_repeats,
    reverse_complement_dna,
    sanitize_dna,
    search_gene_candidates,
    splice_pre_mrna,
    translate_mature_mrna,
)


# ====================================================================
# 一、配色与数据模型
# ====================================================================

COLORS = {
    "bg": "#eef3f8",
    "card": "#ffffff",
    "border": "#dbe4f0",
    "text": "#172033",
    "muted": "#64748b",
    "accent": "#2563eb",
    "accent_dark": "#1d4ed8",
    "header": "#0f172a",
    "header_muted": "#aebbd0",
    "green": "#047857",
    "green_bg": "#e9f9f1",
    "amber": "#b45309",
    "amber_bg": "#fff7e6",
    "red": "#b91c1c",
    "red_bg": "#fff1f1",
    "blue_bg": "#edf4ff",
    "soft": "#f7faff",
    "sequence_bg": "#f8fafc",
}

APP_TITLE = "DNA 序列分析教学工具"


@dataclass
class GuiRequest:
    """一次 GUI 分析请求。"""
    dna: str
    source: str
    min_aa: int = MIN_ORF_AA
    require_tata: bool = REQUIRE_TATA_BOX
    remove_repeats: bool = REMOVE_REPEAT_REGIONS
    use_uniprot: bool = True
    use_similarity: bool = UNIPROT_SIMILARITY_ENABLED
    warnings: list[str] = field(default_factory=list)


@dataclass
class CandidateAnalysis:
    """一个候选在 GUI 层的完整下游结果。"""
    candidate: CandidateGene
    repeats: list[TandemRepeatRegion]
    splice_input: str
    removed_intervals: list[tuple[int, int]]
    splicing: SplicingResult
    translation: TranslationResult
    uniprot: Optional[UniProtResult] = None


@dataclass
class GuiAnalysis:
    """GUI 一次运行的完整结果。"""
    request: GuiRequest
    search: GeneSearchResult
    details: list[CandidateAnalysis]
    elapsed_seconds: float

    @property
    def gc_percent(self) -> float:
        canonical = sum(self.request.dna.count(base) for base in "ACGT")
        gc = self.request.dna.count("G") + self.request.dna.count("C")
        return (gc / canonical * 100.0) if canonical else 0.0


ProgressCallback = Callable[[str], None]


# ====================================================================
# 二、后台分析与报告生成（不接触任何 Tk 控件）
# ====================================================================

def run_gui_analysis(request: GuiRequest, progress: Optional[ProgressCallback] = None) -> GuiAnalysis:
    """
    在后台线程中执行完整分析。只操作纯 Python 对象，不访问 Tk 控件。
    """
    report = progress or (lambda _message: None)
    started = time.perf_counter()
    options = RuntimeOptions(
        min_orf_aa=request.min_aa,
        require_tata_box=request.require_tata,
        remove_repeat_regions=request.remove_repeats,
        max_candidates=MAX_CANDIDATES_TO_ANALYZE,
    )

    report("正在扫描正链与反向互补链的六个读码框…")
    search_result = search_gene_candidates(request.dna, options)

    details: list[CandidateAnalysis] = []
    if search_result.selected_candidates:
        uniprot_client = (
            UniProtClient(enable_similarity=request.use_similarity)
            if request.use_uniprot
            else None
        )
        total = len(search_result.selected_candidates)
        for index, candidate in enumerate(search_result.selected_candidates, start=1):
            report(
                f"正在处理候选 {candidate.candidate_id}（{index}/{total}）："
                "转录、重复粗筛、剪接与翻译…"
            )
            repeats = find_tandem_repeats(candidate.pre_mrna) if SCREEN_TANDEM_REPEATS else []
            removed_intervals: list[tuple[int, int]] = []
            if request.remove_repeats and repeats:
                splice_input, removed_intervals = remove_marked_repeats(
                    candidate.pre_mrna, repeats
                )
            else:
                splice_input = candidate.pre_mrna

            splicing = splice_pre_mrna(splice_input)
            translation = translate_mature_mrna(splicing.mature_mrna)
            uniprot_result: Optional[UniProtResult] = None
            if request.use_uniprot and translation.protein:
                similarity_text = (
                    "；精确未命中时进行相似蛋白 BLAST，可能需要十几秒"
                    if request.use_similarity else "；相似性搜索已关闭"
                )
                report(
                    f"候选 {candidate.candidate_id}：正在查询 UniProt"
                    f"（限流约 {UNIPROT_MIN_INTERVAL_SECONDS:.0f} 秒/请求，"
                    f"超时 {UNIPROT_TIMEOUT_SECONDS:.0f} 秒{similarity_text}）…"
                )
                if uniprot_client is None:
                    # 理论上只在 request.use_uniprot=False 时发生；保留保护。
                    uniprot_client = UniProtClient(
                        enable_similarity=request.use_similarity
                    )
                uniprot_result = uniprot_client.lookup_protein(translation.protein)

            details.append(
                CandidateAnalysis(
                    candidate=candidate,
                    repeats=repeats,
                    splice_input=splice_input,
                    removed_intervals=removed_intervals,
                    splicing=splicing,
                    translation=translation,
                    uniprot=uniprot_result,
                )
            )

    elapsed = time.perf_counter() - started
    report("分析完成，正在整理结果…")
    return GuiAnalysis(
        request=request,
        search=search_result,
        details=details,
        elapsed_seconds=elapsed,
    )


def numbered_sequence(sequence: str, width: int = 80, unit: str = "nt") -> str:
    """生成带行号、按固定宽度折行的序列文本。"""
    if not sequence:
        return "（空序列）\n"
    lines = []
    for start in range(0, len(sequence), width):
        chunk = sequence[start:start + width]
        lines.append(f"{start + 1:>6}  {chunk}")
    return "\n".join(lines) + "\n"


def short_query(query: str, limit: int = 180) -> str:
    if len(query) <= limit:
        return query
    return query[:limit - 25] + "…" + query[-18:]


def api_status_text(detail: CandidateAnalysis, use_uniprot: bool) -> str:
    if not use_uniprot:
        return "已关闭 UniProt 查询"
    if not detail.translation.protein:
        return "无蛋白产物，未查询"
    result = detail.uniprot
    if result is None:
        return "未查询"
    if result.status == "matched" and result.hits:
        hit = result.hits[0]
        if hit.match_type == "similar" or result.method == "genomenet_blast":
            combined = hit.combined_match_percent
            combined_text = f" · 综合匹配度 {combined:.1f}%" if combined is not None else ""
            return f"相似: {hit.accession} · {hit.protein_name}{combined_text}"
        return f"精确命中: {hit.accession} · {hit.protein_name}"
    if result.status == "no_match":
        return "无精确或相似匹配"
    if result.status == "error":
        return "查询失败（见 UniProt 标签页）"
    return result.status


def build_text_report(analysis: GuiAnalysis) -> str:
    """把 GUI 结果导出为结构化、可存档的纯文本报告。"""
    req = analysis.request
    search = analysis.search
    lines: list[str] = [
        "=" * 78,
        "DNA 序列分析教学报告",
        "=" * 78,
        f"输入来源: {req.source}",
        f"输入长度: {len(req.dna)} bp | GC: {analysis.gc_percent:.2f}%",
        f"最小 ORF: {req.min_aa} aa | TATA 过滤: {'开' if req.require_tata else '关'}"
        f" | 剔除重复: {'开' if req.remove_repeats else '关'}"
        f" | UniProt: {'开' if req.use_uniprot else '关'}"
        f" | 相似性 BLAST: {'开' if req.use_similarity else '关'}",
        f"完整 ORF 总数: {len(search.all_complete_orfs)}"
        f" | 通过长度阈值: {len(search.length_passing_candidates)}"
        f" | 详细分析: {len(analysis.details)}",
        f"分析用时: {analysis.elapsed_seconds:.2f} 秒",
        "",
        "重要声明: 本报告来自教学简化模拟器，不是专业基因预测或功能实验，",
        "不能用于正式科研、临床诊断、药物设计或生产决策。",
    ]
    if req.warnings:
        lines.append("")
        lines.append("[输入提示]")
        lines.extend(f"- {warning}" for warning in req.warnings)

    lines.extend(["", "=" * 78, "[1] 原始 DNA 序列", "=" * 78])
    lines.append(numbered_sequence(req.dna, width=80).rstrip())
    lines.extend(["", "[反向互补链 5'→3']"])
    lines.append(numbered_sequence(reverse_complement_dna(req.dna), width=80).rstrip())

    lines.extend(["", "=" * 78, "[2] 六读码框 ORF 与候选", "=" * 78])
    lines.append("完整 ORF（ATG 起始且含同框终止密码子）")
    if not search.length_passing_candidates:
        lines.append("无 ORF 通过最小氨基酸长度阈值。")
    for idx, candidate in enumerate(search.length_passing_candidates, start=1):
        orf = candidate.orf
        lines.append(
            f"{idx:>2}. {orf.orf_id} | {orf.strand} 链 | 读码框 {orf.strand}{orf.frame}"
            f" | 输入正向 {interval_1based(orf.genomic_start, orf.genomic_end)}"
            f" | {orf.aa_length} aa"
            f" | TATA {len(candidate.tata_hits)} 个"
            f" | AATAAA {len(candidate.polya_hits)} 个"
        )

    if analysis.details:
        lines.extend(["", "=" * 78, "[3] 候选详细结果", "=" * 78])
        for detail in analysis.details:
            cand = detail.candidate
            orf = cand.orf
            lines.extend(
                [
                    "",
                    "-" * 78,
                    f"候选 {cand.candidate_id} | {orf.orf_id} | "
                    f"{'正链' if orf.strand == '+' else '负链'} | "
                    f"读码框 {orf.strand}{orf.frame}",
                    f"ORF: {interval_1based(orf.genomic_start, orf.genomic_end)}（输入正向）"
                    f" | 候选区域: {interval_1based(cand.region_start, cand.region_end)}",
                    f"蛋白长度: {orf.aa_length} aa",
                    f"TATA 上游命中: {len(cand.tata_hits)} 个 | "
                    f"AATAAA 下游命中: {len(cand.polya_hits)} 个",
                    "",
                    "[pre-mRNA]",
                    numbered_sequence(cand.pre_mrna, width=80).rstrip(),
                    "",
                    "[重复粗筛]",
                ]
            )
            if detail.repeats:
                for repeat in detail.repeats:
                    lines.append(
                        f"- pre-mRNA {interval_1based(repeat.start, repeat.end)}"
                        f"（{repeat.length} nt），{repeat.unit} × 约 {repeat.copies} 次"
                    )
                if detail.removed_intervals:
                    lines.append("已按配置剔除以上重复标记区域。")
                else:
                    lines.append("仅标记，未剔除。")
            else:
                lines.append("未发现满足参数的串联重复。")

            lines.extend(["", "[剪接]"])
            if detail.splicing.introns:
                for number, intron in enumerate(detail.splicing.introns, start=1):
                    lines.append(
                        f"- 内含子 {number}: pre-mRNA "
                        f"{interval_1based(intron.start, intron.end)}"
                        f"（{intron.length} nt），{intron.donor}…{intron.acceptor}"
                    )
                lines.append(
                    "外显子: "
                    + " | ".join(
                        f"E{number} {interval_1based(start, end)}（{end - start} nt）"
                        for number, (start, end) in enumerate(detail.splicing.exons, start=1)
                    )
                )
            else:
                lines.append("未发现满足长度约束的 GU…AG 内含子。")
            for note in detail.splicing.notes:
                lines.append(f"提示: {note}")

            lines.extend(
                [
                    "",
                    "[成熟 mRNA]",
                    numbered_sequence(detail.splicing.mature_mrna, width=80).rstrip(),
                    "",
                    "[翻译]",
                ]
            )
            translation = detail.translation
            if translation.start_index is None:
                lines.append("未找到 AUG，无蛋白产物。")
            else:
                lines.append(
                    f"第一个 AUG 起点: {translation.start_index + 1}（成熟 mRNA 1-based）"
                )
                if translation.stop_codon and translation.stop_codon_start is not None:
                    lines.append(
                        f"终止密码子: {translation.stop_codon}，"
                        f"位置 {translation.stop_codon_start + 1}"
                    )
                else:
                    lines.append("未找到同框终止密码子。")
                lines.append(f"蛋白长度: {len(translation.protein)} aa")
                for warning in translation.warnings:
                    lines.append(f"提示: {warning}")
                lines.append("")
                lines.append(numbered_sequence(translation.protein, width=60, unit="aa").rstrip())

            lines.extend(["", "[UniProt / 相似蛋白]"])
            if not req.use_uniprot:
                lines.append("已关闭 UniProt 查询。")
            elif detail.uniprot is None:
                lines.append("无蛋白产物或未发起查询。")
            else:
                result = detail.uniprot
                lines.append(f"状态: {result.status} | 方式: {result.method or '—'}")
                if result.message:
                    lines.append(f"说明: {result.message}")
                if result.method == "genomenet_blast" and result.hits:
                    top_hit = result.hits[0]
                    lines.append(
                        f"推测参考（最相似条目）: {top_hit.accession} | "
                        f"{top_hit.protein_name} | {top_hit.organism} | "
                        f"一致度 {top_hit.identity_percent:.1f}% | "
                        f"覆盖度 {top_hit.coverage_percent:.1f}% | "
                        f"综合匹配度 {top_hit.combined_match_percent:.1f}%"
                    )
                lines.append(f"查询式: {short_query(result.query)}")
                for hit_index, hit in enumerate(result.hits, start=1):
                    lines.append(
                        f"{hit_index}. {hit.accession} | {hit.protein_name} | "
                        f"{hit.organism} | "
                        f"{'reviewed/Swiss-Prot' if hit.reviewed else 'unreviewed/TrEMBL'}"
                    )
                    if hit.match_type == "similar" or result.method == "genomenet_blast":
                        evalue = f"{hit.evalue:.2g}" if hit.evalue is not None else "—"
                        bitscore = f"{hit.bitscore:.1f}" if hit.bitscore is not None else "—"
                        lines.append(
                            f"   匹配: 一致度 {hit.identity_percent:.1f}% | "
                            f"覆盖度 {hit.coverage_percent:.1f}% | "
                            f"综合匹配度 {hit.combined_match_percent:.1f}% | "
                            f"E-value {evalue} | bit score {bitscore}"
                        )
                    lines.append(
                        f"   UniProt: https://www.uniprot.org/uniprotkb/{hit.accession}/entry"
                    )
                    lines.append(f"   功能: {hit.function or '（无 FUNCTION 注释）'}")

    lines.extend(
        [
            "",
            "=" * 78,
            "局限性",
            "=" * 78,
            "1. ORF 仅按 ATG 起始、同框终止扫描，未建模启动子、Kozak、内含子、",
            "   RNA-seq 证据或物种参数；不是 AUGUSTUS/BRAKER 等专业基因预测器。",
            "2. TATAAA/AATAAA 只是六碱基字面标记；不代表真实启动子或 poly(A) 位点。",
            "3. GU-AG 剪接只做教学简化，忽略分支点、多嘧啶区、增强子和可变剪接。",
            "4. 串联重复粗筛不等于转座子鉴定，默认只标记、不剔除。",
            "5. 翻译使用标准遗传密码，未考虑 Kozak、uORF、RNA 编辑或翻译后修饰。",
            "6. UniProt 精确命中只表示数据库中存在相同序列条目；相似性 BLAST 命中",
            "   只是通过序列相似性推断可能相关的蛋白，不能证明本 DNA 表达该蛋白。",
            "7. 相似性排名按“一致度 × 查询覆盖度”的综合匹配度计算；不是实验验证。",
            "8. 结果不得用于正式科研、临床诊断、药物设计、育种或生产决策。",
            "",
        ]
    )
    return "\n".join(lines)


# ====================================================================
# 三、通用 GUI 工具
# ====================================================================

def choose_first_font(root: tk.Misc, candidates: list[str], fallback: str) -> str:
    available = set(tkfont.families(root))
    for name in candidates:
        if name in available:
            return name
    return fallback


def center_window(window: tk.Toplevel | tk.Tk, width: Optional[int] = None, height: Optional[int] = None) -> None:
    window.update_idletasks()
    if width is not None and height is not None:
        window.geometry(f"{width}x{height}")
    win_width = width or window.winfo_width()
    win_height = height or window.winfo_height()
    screen_width = window.winfo_screenwidth()
    screen_height = window.winfo_screenheight()
    x = max(0, (screen_width - win_width) // 2)
    y = max(0, (screen_height - win_height) // 2 - 18)
    window.geometry(f"{win_width}x{win_height}+{x}+{y}")


def make_button(
    parent: tk.Misc,
    text: str,
    command: Callable[[], None],
    *,
    ui_font: tuple,
    primary: bool = False,
    dark: bool = False,
    small: bool = False,
) -> tk.Button:
    if primary:
        bg, fg, active_bg = COLORS["accent"], "#ffffff", COLORS["accent_dark"]
    elif dark:
        bg, fg, active_bg = "#1e293b", "#e2e8f0", "#334155"
    else:
        bg, fg, active_bg = "#e8eef8", COLORS["text"], "#dbe7f8"
    padx, pady = (16, 9) if not small else (10, 5)
    return tk.Button(
        parent,
        text=text,
        command=command,
        font=(ui_font[0], 10 if not small else 9, "bold"),
        bg=bg,
        fg=fg,
        activebackground=active_bg,
        activeforeground=fg if not primary else "#ffffff",
        relief="flat",
        bd=0,
        cursor="hand2",
        padx=padx,
        pady=pady,
        disabledforeground="#94a3b8",
    )


def make_text_panel(
    parent: tk.Misc,
    *,
    ui_font_name: str,
    mono_font_name: str,
    bg: str = COLORS["card"],
    height: int = 20,
) -> tuple[tk.Frame, tk.Text]:
    frame = tk.Frame(
        parent,
        bg=bg,
        highlightbackground=COLORS["border"],
        highlightthickness=1,
    )
    text = tk.Text(
        frame,
        wrap="char",
        font=(ui_font_name, 10),
        bg=bg,
        fg=COLORS["text"],
        insertbackground=COLORS["text"],
        relief="flat",
        bd=0,
        padx=14,
        pady=12,
        height=height,
        spacing1=0,
        spacing3=2,
        selectbackground="#bfdbfe",
    )
    scrollbar = ttk.Scrollbar(frame, orient="vertical", command=text.yview)
    text.configure(yscrollcommand=scrollbar.set)
    text.grid(row=0, column=0, sticky="nsew")
    scrollbar.grid(row=0, column=1, sticky="ns")
    frame.grid_rowconfigure(0, weight=1)
    frame.grid_columnconfigure(0, weight=1)
    configure_text_tags(text, ui_font_name=ui_font_name, mono_font_name=mono_font_name)
    text.configure(state="disabled")
    return frame, text


def configure_text_tags(text: tk.Text, *, ui_font_name: str, mono_font_name: str) -> None:
    """为只读文本面板统一配置标题、正文、序列、状态和超链接样式。"""
    text.tag_configure(
        "h1",
        font=(ui_font_name, 16, "bold"),
        foreground=COLORS["text"],
        spacing1=4,
        spacing3=7,
    )
    text.tag_configure(
        "h2",
        font=(ui_font_name, 12, "bold"),
        foreground=COLORS["accent"],
        spacing1=14,
        spacing3=5,
    )
    text.tag_configure(
        "h3",
        font=(ui_font_name, 10, "bold"),
        foreground=COLORS["text"],
        spacing1=8,
        spacing3=3,
    )
    text.tag_configure(
        "body",
        font=(ui_font_name, 10),
        foreground=COLORS["text"],
        spacing1=1,
        spacing3=3,
    )
    text.tag_configure(
        "muted",
        font=(ui_font_name, 9),
        foreground=COLORS["muted"],
        spacing1=1,
        spacing3=2,
    )
    text.tag_configure(
        "mono",
        font=(mono_font_name, 9),
        foreground="#0f172a",
        background=COLORS["sequence_bg"],
        lmargin1=10,
        lmargin2=10,
        rmargin=8,
        spacing1=0,
        spacing3=0,
    )
    text.tag_configure("good", foreground=COLORS["green"], font=(ui_font_name, 10, "bold"))
    text.tag_configure("warn", foreground=COLORS["amber"], font=(ui_font_name, 10, "bold"))
    text.tag_configure("bad", foreground=COLORS["red"], font=(ui_font_name, 10, "bold"))
    text.tag_configure("link", foreground=COLORS["accent"], underline=True)
    text.tag_configure(
        "chip_blue",
        background=COLORS["blue_bg"],
        foreground=COLORS["accent"],
        font=(ui_font_name, 9, "bold"),
    )
    text.tag_configure(
        "chip_green",
        background=COLORS["green_bg"],
        foreground=COLORS["green"],
        font=(ui_font_name, 9, "bold"),
    )
    text.tag_configure(
        "chip_warn",
        background=COLORS["amber_bg"],
        foreground=COLORS["amber"],
        font=(ui_font_name, 9, "bold"),
    )


def clear_children(widget: tk.Misc) -> None:
    for child in widget.winfo_children():
        child.destroy()


def set_rich_text(text_widget: tk.Text, builder: Callable[[tk.Text], None]) -> None:
    """以只读状态安全替换一个 Text 控件的全部内容。"""
    text_widget.configure(state="normal")
    text_widget.delete("1.0", "end")
    builder(text_widget)
    text_widget.configure(state="disabled")
    text_widget.yview_moveto(0.0)


def insert_sequence_block(text: tk.Text, sequence: str, *, width: int = 64, unit: str = "nt") -> None:
    if not sequence:
        text.insert("end", "（空序列）\n", "muted")
        return
    for start in range(0, len(sequence), width):
        chunk = sequence[start:start + width]
        text.insert("end", f"{start + 1:>6}  {chunk}\n", "mono")


def card_frame(parent: tk.Misc, *, bg: str = COLORS["card"]) -> tk.Frame:
    return tk.Frame(
        parent,
        bg=bg,
        highlightbackground=COLORS["border"],
        highlightthickness=1,
    )


def metric_card(
    parent: tk.Misc,
    *,
    title: str,
    value: str,
    hint: str,
    accent: str,
    ui_font_name: str,
) -> tk.Frame:
    card = card_frame(parent)
    card.grid_propagate(True)
    tk.Frame(card, bg=accent, height=4).pack(fill="x")
    body = tk.Frame(card, bg=COLORS["card"], padx=14, pady=12)
    body.pack(fill="both", expand=True)
    tk.Label(
        body,
        text=title,
        bg=COLORS["card"],
        fg=COLORS["muted"],
        font=(ui_font_name, 9),
        anchor="w",
    ).pack(fill="x")
    tk.Label(
        body,
        text=value,
        bg=COLORS["card"],
        fg=COLORS["text"],
        font=(ui_font_name, 18, "bold"),
        anchor="w",
    ).pack(fill="x", pady=(2, 1))
    tk.Label(
        body,
        text=hint,
        bg=COLORS["card"],
        fg=COLORS["muted"],
        font=(ui_font_name, 9),
        anchor="w",
        justify="left",
        wraplength=235,
    ).pack(fill="x")
    return card


class SequenceInputDialog(tk.Toplevel):
    """启动/新建分析时弹出的 DNA 输入与参数窗口。"""

    def __init__(
        self,
        master: tk.Misc,
        *,
        ui_font_name: str,
        mono_font_name: str,
        on_submit: Callable[[GuiRequest], None],
        on_cancel: Callable[[], None],
        initial: bool = False,
    ) -> None:
        super().__init__(master)
        self.ui_font_name = ui_font_name
        self.mono_font_name = mono_font_name
        self._on_submit = on_submit
        self._on_cancel = on_cancel
        self._initial = initial
        self._loaded_sequence: Optional[str] = None
        self._loaded_source: Optional[str] = None
        self._load_warnings: list[str] = []
        self._closed = False
        self._acted = False

        self.title("输入 DNA 序列 · DNA 序列分析教学工具")
        self.configure(bg=COLORS["bg"])
        self.minsize(860, 660)
        self.resizable(True, True)
        screen_width = self.winfo_screenwidth()
        screen_height = self.winfo_screenheight()
        initial_width = max(860, min(950, screen_width - 100))
        initial_height = max(660, min(750, screen_height - 100))
        center_window(self, initial_width, initial_height)
        if not initial:
            try:
                self.transient(master)
            except tk.TclError:
                pass
        self.protocol("WM_DELETE_WINDOW", self._cancel)
        self._build_ui()

        self.bind("<Escape>", lambda _event: self._cancel())
        self.bind("<Control-Return>", lambda _event: self._submit())
        self.after(80, self._make_modal_and_focus)

    # ---------- 布局 ----------
    def _build_ui(self) -> None:
        self.grid_rowconfigure(1, weight=1)
        self.grid_columnconfigure(0, weight=1)

        header = tk.Frame(self, bg=COLORS["header"], padx=24, pady=10)
        header.grid(row=0, column=0, sticky="ew")
        tk.Label(
            header,
            text="输入 DNA 序列",
            bg=COLORS["header"],
            fg="#ffffff",
            font=(self.ui_font_name, 16, "bold"),
            anchor="w",
        ).pack(fill="x")
        tk.Label(
            header,
            text="粘贴 DNA 或载入 FASTA；点击“开始分析”后按模块显示结果。",
            bg=COLORS["header"],
            fg=COLORS["header_muted"],
            font=(self.ui_font_name, 9),
            anchor="w",
            justify="left",
            wraplength=900,
        ).pack(fill="x", pady=(2, 0))

        body = tk.Frame(self, bg=COLORS["bg"], padx=20, pady=10)
        body.grid(row=1, column=0, sticky="nsew")
        body.grid_columnconfigure(0, weight=1)
        body.grid_rowconfigure(2, weight=1)

        hint = tk.Label(
            body,
            text=(
                "可粘贴原始 DNA，或从 FASTA 载入。原始 DNA 不会上传；"
                "仅翻译后的蛋白序列会用于 UniProt 查询。若开启相似性 BLAST 且精确未命中，"
                "蛋白序列还会提交给 GenomeNet BLAST。"
            ),
            bg=COLORS["bg"],
            fg=COLORS["muted"],
            font=(self.ui_font_name, 9),
            justify="left",
            anchor="w",
            wraplength=890,
        )
        hint.grid(row=0, column=0, sticky="ew", pady=(0, 5))

        action_row = tk.Frame(body, bg=COLORS["bg"])
        action_row.grid(row=1, column=0, sticky="ew", pady=(0, 8))
        make_button(
            action_row,
            "填入内置示例",
            self._fill_demo,
            ui_font=(self.ui_font_name, 10),
            small=True,
        ).pack(side="left")
        make_button(
            action_row,
            "从 FASTA 载入…",
            self._load_fasta_dialog,
            ui_font=(self.ui_font_name, 10),
            small=True,
        ).pack(side="left", padx=(8, 0))
        make_button(
            action_row,
            "清空",
            self._clear_input,
            ui_font=(self.ui_font_name, 10),
            small=True,
        ).pack(side="left", padx=(8, 0))
        self.source_var = tk.StringVar(value="来源: 手动输入")
        tk.Label(
            action_row,
            textvariable=self.source_var,
            bg=COLORS["bg"],
            fg=COLORS["muted"],
            font=(self.ui_font_name, 9),
            anchor="e",
            justify="right",
            wraplength=400,
        ).pack(side="right")

        editor_frame = tk.Frame(
            body,
            bg=COLORS["card"],
            highlightbackground=COLORS["border"],
            highlightthickness=1,
        )
        editor_frame.grid(row=2, column=0, sticky="nsew")
        editor_frame.grid_rowconfigure(0, weight=1)
        editor_frame.grid_columnconfigure(0, weight=1)
        self.sequence_text = tk.Text(
            editor_frame,
            wrap="char",
            font=(self.mono_font_name, 10),
            bg="#ffffff",
            fg=COLORS["text"],
            insertbackground=COLORS["accent"],
            relief="flat",
            bd=0,
            padx=14,
            pady=10,
            undo=True,
            height=5,
            spacing3=1,
            selectbackground="#bfdbfe",
        )
        editor_scroll = ttk.Scrollbar(editor_frame, orient="vertical", command=self.sequence_text.yview)
        self.sequence_text.configure(yscrollcommand=editor_scroll.set)
        self.sequence_text.grid(row=0, column=0, sticky="nsew")
        editor_scroll.grid(row=0, column=1, sticky="ns")
        self.sequence_text.bind("<<Modified>>", self._on_modified)
        self.sequence_text.edit_modified(False)

        counter_row = tk.Frame(body, bg=COLORS["bg"])
        counter_row.grid(row=3, column=0, sticky="ew", pady=(5, 8))
        self.counter_var = tk.StringVar(value="当前输入: 0 个碱基")
        tk.Label(
            counter_row,
            textvariable=self.counter_var,
            bg=COLORS["bg"],
            fg=COLORS["muted"],
            font=(self.ui_font_name, 9),
            anchor="w",
        ).pack(side="left")

        settings = tk.Frame(
            body,
            bg=COLORS["card"],
            highlightbackground=COLORS["border"],
            highlightthickness=1,
            padx=16,
            pady=8,
        )
        settings.grid(row=4, column=0, sticky="ew")
        settings.grid_columnconfigure(0, weight=1)
        params_row = tk.Frame(settings, bg=COLORS["card"])
        params_row.grid(row=0, column=0, sticky="ew")
        tk.Label(
            params_row,
            text="分析参数",
            bg=COLORS["card"],
            fg=COLORS["text"],
            font=(self.ui_font_name, 10, "bold"),
        ).pack(side="left")
        tk.Label(
            params_row,
            text="最小 ORF 长度（aa）",
            bg=COLORS["card"],
            fg=COLORS["muted"],
            font=(self.ui_font_name, 9),
        ).pack(side="left", padx=(18, 6))
        self.min_aa_var = tk.StringVar(value=str(MIN_ORF_AA))
        spin = ttk.Spinbox(
            params_row,
            from_=0,
            to=100000,
            width=7,
            textvariable=self.min_aa_var,
            font=(self.ui_font_name, 10),
        )
        spin.pack(side="left")

        self.require_tata_var = tk.BooleanVar(value=REQUIRE_TATA_BOX)
        self.remove_repeats_var = tk.BooleanVar(value=REMOVE_REPEAT_REGIONS)
        self.use_uniprot_var = tk.BooleanVar(value=True)
        self.use_similarity_var = tk.BooleanVar(value=UNIPROT_SIMILARITY_ENABLED)

        option_style = {
            "bg": COLORS["card"],
            "activebackground": COLORS["card"],
            "font": (self.ui_font_name, 9),
            "fg": COLORS["text"],
            "selectcolor": "#ffffff",
            "highlightthickness": 0,
            "bd": 0,
            "cursor": "hand2",
            "anchor": "w",
            "justify": "left",
            "wraplength": 200,
        }
        options_row = tk.Frame(settings, bg=COLORS["card"])
        options_row.grid(row=1, column=0, sticky="ew", pady=(4, 0))
        for column in range(4):
            options_row.grid_columnconfigure(column, weight=1)
        tk.Checkbutton(
            options_row,
            text="只保留 TATA 盒候选",
            variable=self.require_tata_var,
            **option_style,
        ).grid(row=0, column=0, sticky="w")
        tk.Checkbutton(
            options_row,
            text="剔除重复标记区域",
            variable=self.remove_repeats_var,
            **option_style,
        ).grid(row=0, column=1, sticky="w", padx=(8, 0))
        tk.Checkbutton(
            options_row,
            text="查询 UniProt",
            variable=self.use_uniprot_var,
            **option_style,
        ).grid(row=0, column=2, sticky="w", padx=(8, 0))
        tk.Checkbutton(
            options_row,
            text="相似性 BLAST（较慢）",
            variable=self.use_similarity_var,
            **option_style,
        ).grid(row=0, column=3, sticky="w", padx=(8, 0))

        tk.Label(
            settings,
            text=(
                f"最多分析排序前 {MAX_CANDIDATES_TO_ANALYZE} 个候选。"
                "相似性 BLAST 只在 UniProt 精确未命中时触发，通常需要 10–30 秒；"
                "原始 DNA 不会上传。TATA/AATAAA 与重复筛查均为教学标记。"
            ),
            bg=COLORS["card"],
            fg=COLORS["muted"],
            font=(self.ui_font_name, 8),
            anchor="w",
            justify="left",
            wraplength=880,
        ).grid(row=2, column=0, sticky="w", pady=(4, 0))

        self.error_var = tk.StringVar(value="")
        self.error_label = tk.Label(
            body,
            textvariable=self.error_var,
            bg=COLORS["bg"],
            fg=COLORS["red"],
            font=(self.ui_font_name, 9, "bold"),
            anchor="w",
            justify="left",
            wraplength=880,
        )
        self.error_label.grid(row=5, column=0, sticky="ew", pady=(6, 0))

        footer = tk.Frame(self, bg=COLORS["bg"], padx=20, pady=10)
        footer.grid(row=2, column=0, sticky="ew")
        footer.grid_columnconfigure(0, weight=1)
        tk.Label(
            footer,
            text="Ctrl+Enter 开始分析 · Esc 取消",
            bg=COLORS["bg"],
            fg=COLORS["muted"],
            font=(self.ui_font_name, 9),
            anchor="w",
        ).grid(row=0, column=0, sticky="w")
        make_button(
            footer,
            "取消",
            self._cancel,
            ui_font=(self.ui_font_name, 10),
        ).grid(row=0, column=1, sticky="e", padx=(0, 8))
        make_button(
            footer,
            "开始分析",
            self._submit,
            ui_font=(self.ui_font_name, 10),
            primary=True,
        ).grid(row=0, column=2, sticky="e")

    # ---------- 窗口生命周期 ----------
    def _make_modal_and_focus(self) -> None:
        if self._closed:
            return
        try:
            self.grab_set()
        except tk.TclError:
            pass
        self.lift()
        self.focus_force()
        self.sequence_text.focus_set()

    def _close(self) -> None:
        if self._closed:
            return
        self._closed = True
        try:
            self.grab_release()
        except tk.TclError:
            pass
        self.destroy()

    def _cancel(self) -> None:
        if self._acted:
            return
        self._acted = True
        self._close()
        self._on_cancel()

    # ---------- 输入辅助 ----------
    def _set_sequence(self, dna: str, source: str, warnings: Optional[list[str]] = None) -> None:
        self.sequence_text.delete("1.0", "end")
        self.sequence_text.insert("1.0", dna)
        self.sequence_text.edit_modified(False)
        self._loaded_sequence = dna
        self._loaded_source = source
        self._load_warnings = list(warnings or [])
        self.source_var.set(f"来源: {source}")
        self.counter_var.set(f"当前输入: {len(dna)} 个碱基")
        self.error_var.set("")

    def _fill_demo(self) -> None:
        dna, warnings = sanitize_dna(DEMO_DNA_SEQUENCE)
        self._set_sequence(dna, DEMO_SOURCE_LABEL, warnings)

    def _clear_input(self) -> None:
        self.sequence_text.delete("1.0", "end")
        self.sequence_text.edit_modified(False)
        self._loaded_sequence = None
        self._loaded_source = None
        self._load_warnings = []
        self.source_var.set("来源: 手动输入")
        self.counter_var.set("当前输入: 0 个碱基")
        self.error_var.set("")
        self.sequence_text.focus_set()

    def _on_modified(self, _event: tk.Event | None = None) -> None:
        if not self.sequence_text.edit_modified():
            return
        self.sequence_text.edit_modified(False)
        raw = self.sequence_text.get("1.0", "end")
        compact = "".join(raw.split()).upper()
        self.counter_var.set(f"当前输入: {len(compact)} 个碱基")
        if self._loaded_sequence is not None and compact == self._loaded_sequence:
            assert self._loaded_source is not None
            self.source_var.set(f"来源: {self._loaded_source}")
        elif compact:
            self.source_var.set("来源: 手动输入（已编辑）")
        else:
            self.source_var.set("来源: 手动输入")

    def _load_fasta_dialog(self) -> None:
        path = filedialog.askopenfilename(
            parent=self,
            title="载入 FASTA 文件",
            filetypes=[("FASTA", "*.fa *.fasta *.fna *.fas *.ffn *.faa"), ("所有文件", "*.*")],
        )
        if not path:
            return
        try:
            records = list(SeqIO.parse(path, "fasta"))
            if not records:
                raise ValueError("该 FASTA 中没有记录。")
            if len(records) == 1:
                record_id = records[0].id
            else:
                ids = [record.id for record in records]
                preview = ", ".join(ids[:10]) + (" …" if len(ids) > 10 else "")
                record_id = simpledialog.askstring(
                    "选择 FASTA 记录",
                    f"文件包含 {len(records)} 条记录：\n{preview}\n\n请输入要载入的 record ID：",
                    initialvalue=ids[0],
                    parent=self,
                )
                if not record_id:
                    return
            dna, source, warnings = load_fasta_sequence(path, record_id)
            self._set_sequence(dna, source, warnings)
        except Exception as exc:
            self.error_var.set(f"FASTA 载入失败: {exc}")

    def _submit(self) -> None:
        if self._acted:
            return
        raw = self.sequence_text.get("1.0", "end")
        try:
            dna, seq_warnings = sanitize_dna(raw)
        except ValueError as exc:
            self.error_var.set(f"输入无效: {exc}")
            self.sequence_text.focus_set()
            return
        try:
            min_aa = int(self.min_aa_var.get().strip())
            if min_aa < 0:
                raise ValueError
        except ValueError:
            self.error_var.set("最小 ORF 长度必须是非负整数。")
            return

        compact = "".join(raw.split()).upper()
        if self._loaded_sequence is not None and compact == self._loaded_sequence and self._loaded_source:
            source = self._loaded_source
            warnings = list(dict.fromkeys(self._load_warnings + seq_warnings))
        else:
            source = "手动输入的 DNA 序列"
            warnings = seq_warnings

        request = GuiRequest(
            dna=dna,
            source=source,
            min_aa=min_aa,
            require_tata=bool(self.require_tata_var.get()),
            remove_repeats=bool(self.remove_repeats_var.get()),
            use_uniprot=bool(self.use_uniprot_var.get()),
            use_similarity=bool(self.use_similarity_var.get()),
            warnings=warnings,
        )
        self._acted = True
        self._close()
        self._on_submit(request)


# ====================================================================
# 四、分析结果标签页
# ====================================================================

class ScrollableFrame(tk.Frame):
    """简单的可滚动容器，用于概览页。"""

    def __init__(self, parent: tk.Misc, *, bg: str = COLORS["bg"]) -> None:
        super().__init__(parent, bg=bg)
        self.canvas = tk.Canvas(self, bg=bg, highlightthickness=0, bd=0)
        self.scrollbar = ttk.Scrollbar(self, orient="vertical", command=self.canvas.yview)
        self.canvas.configure(yscrollcommand=self.scrollbar.set)
        self.canvas.pack(side="left", fill="both", expand=True)
        self.scrollbar.pack(side="right", fill="y")
        self.inner = tk.Frame(self.canvas, bg=bg)
        self._window = self.canvas.create_window((0, 0), window=self.inner, anchor="nw")
        self.inner.bind("<Configure>", self._on_inner_configure)
        self.canvas.bind("<Configure>", self._on_canvas_configure)
        self.canvas.bind("<MouseWheel>", self._on_mousewheel)
        self.inner.bind("<MouseWheel>", self._on_mousewheel)

    def _on_inner_configure(self, _event: tk.Event) -> None:
        self.canvas.configure(scrollregion=self.canvas.bbox("all"))

    def _on_canvas_configure(self, event: tk.Event) -> None:
        self.canvas.itemconfigure(self._window, width=event.width)

    def _on_mousewheel(self, event: tk.Event) -> None:
        self.canvas.yview_scroll(int(-event.delta / 120), "units")


class OverviewTab(tk.Frame):
    def __init__(self, parent: tk.Misc, app: "DNAApp") -> None:
        super().__init__(parent, bg=COLORS["bg"])
        self.app = app
        self.scrollable = ScrollableFrame(self, bg=COLORS["bg"])
        self.scrollable.pack(fill="both", expand=True)
        self.inner = self.scrollable.inner
        self.inner.grid_columnconfigure(0, weight=1)

    def _clear(self) -> None:
        clear_children(self.inner)
        for col in range(4):
            self.inner.grid_columnconfigure(col, weight=1)

    def show_loading(self, message: str = "正在分析，请稍候…") -> None:
        self._clear()
        box = card_frame(self.inner)
        box.grid(row=0, column=0, sticky="ew", padx=18, pady=18)
        tk.Label(
            box,
            text="正在分析",
            bg=COLORS["card"],
            fg=COLORS["text"],
            font=(self.app.ui_font_name, 18, "bold"),
            anchor="w",
        ).pack(fill="x", padx=22, pady=(22, 6))
        tk.Label(
            box,
            text=message,
            bg=COLORS["card"],
            fg=COLORS["muted"],
            font=(self.app.ui_font_name, 10),
            anchor="w",
            justify="left",
        ).pack(fill="x", padx=22, pady=(0, 22))

    def show_error(self, message: str, details: str) -> None:
        self._clear()
        box = card_frame(self.inner)
        box.grid(row=0, column=0, sticky="ew", padx=18, pady=18)
        tk.Label(
            box,
            text="分析失败",
            bg=COLORS["card"],
            fg=COLORS["red"],
            font=(self.app.ui_font_name, 18, "bold"),
            anchor="w",
        ).pack(fill="x", padx=22, pady=(20, 4))
        tk.Label(
            box,
            text=message,
            bg=COLORS["card"],
            fg=COLORS["text"],
            font=(self.app.ui_font_name, 10),
            anchor="w",
            justify="left",
            wraplength=1000,
        ).pack(fill="x", padx=22, pady=(0, 10))
        tk.Label(
            box,
            text=details,
            bg=COLORS["card"],
            fg=COLORS["muted"],
            font=(self.app.mono_font_name, 8),
            anchor="w",
            justify="left",
            wraplength=1000,
        ).pack(fill="x", padx=22, pady=(0, 20))

    def render(self, analysis: GuiAnalysis) -> None:
        self._clear()
        req = analysis.request
        search = analysis.search
        details = analysis.details
        first = details[0] if details else None

        heading = tk.Frame(self.inner, bg=COLORS["bg"])
        heading.grid(row=0, column=0, columnspan=4, sticky="ew", padx=18, pady=(16, 10))
        heading.grid_columnconfigure(0, weight=1)
        tk.Label(
            heading,
            text="分析概览",
            bg=COLORS["bg"],
            fg=COLORS["text"],
            font=(self.app.ui_font_name, 20, "bold"),
            anchor="w",
        ).grid(row=0, column=0, sticky="w")
        tk.Label(
            heading,
            text=f"用时 {analysis.elapsed_seconds:.2f} 秒 · {req.source}",
            bg=COLORS["bg"],
            fg=COLORS["muted"],
            font=(self.app.ui_font_name, 9),
            anchor="e",
        ).grid(row=1, column=0, sticky="w", pady=(3, 0))

        cards = tk.Frame(self.inner, bg=COLORS["bg"])
        cards.grid(row=1, column=0, columnspan=4, sticky="ew", padx=18)
        for col in range(4):
            cards.grid_columnconfigure(col, weight=1, uniform="metric")
        orf_value = str(len(search.length_passing_candidates))
        details_value = str(len(details))
        if first:
            protein_value = f"{len(first.translation.protein)} aa"
            protein_hint = api_status_text(first, req.use_uniprot)
            introns = sum(len(detail.splicing.introns) for detail in details)
            metrics = [
                ("输入序列", f"{len(req.dna):,} bp", f"GC {analysis.gc_percent:.2f}%", COLORS["accent"]),
                ("通过长度阈值", orf_value, f"完整 ORF {len(search.all_complete_orfs)} 个", "#7c3aed"),
                ("详细分析候选", details_value, f"内含子合计 {introns} 个", COLORS["green"]),
                ("首选候选蛋白", protein_value, protein_hint, "#d97706"),
            ]
        else:
            metrics = [
                ("输入序列", f"{len(req.dna):,} bp", f"GC {analysis.gc_percent:.2f}%", COLORS["accent"]),
                ("通过长度阈值", orf_value, f"完整 ORF {len(search.all_complete_orfs)} 个", "#7c3aed"),
                ("详细分析候选", details_value, "没有候选进入下游", COLORS["muted"]),
                ("UniProt", "—", "无蛋白可查询", COLORS["muted"]),
            ]
        for index, (title, value, hint, accent) in enumerate(metrics):
            card = metric_card(
                cards,
                title=title,
                value=value,
                hint=hint,
                accent=accent,
                ui_font_name=self.app.ui_font_name,
            )
            card.grid(row=0, column=index, sticky="nsew", padx=(0 if index == 0 else 8, 0))

        # 参数标签
        params = (
            f"最小 ORF {req.min_aa} aa   ·   "
            f"TATA 过滤 {'开' if req.require_tata else '关'}   ·   "
            f"剔除重复 {'开' if req.remove_repeats else '关'}   ·   "
            f"UniProt {'开' if req.use_uniprot else '关'}   ·   "
            f"相似性 BLAST {'开' if req.use_similarity and req.use_uniprot else '关'}"
        )
        tk.Label(
            self.inner,
            text=params,
            bg=COLORS["bg"],
            fg=COLORS["muted"],
            font=(self.app.ui_font_name, 9),
            anchor="w",
        ).grid(row=2, column=0, columnspan=4, sticky="w", padx=20, pady=(10, 2))

        # 候选摘要
        summary = card_frame(self.inner)
        summary.grid(row=3, column=0, columnspan=4, sticky="ew", padx=18, pady=(8, 12))
        summary.grid_columnconfigure(0, weight=1)
        summary_header = tk.Frame(summary, bg=COLORS["card"])
        summary_header.grid(row=0, column=0, sticky="ew", padx=18, pady=(16, 6))
        summary_header.grid_columnconfigure(0, weight=1)
        tk.Label(
            summary_header,
            text="候选与核心结果",
            bg=COLORS["card"],
            fg=COLORS["text"],
            font=(self.app.ui_font_name, 13, "bold"),
            anchor="w",
        ).grid(row=0, column=0, sticky="w")
        if details:
            make_button(
                summary_header,
                "查看剪接与翻译",
                self.app.show_first_detail,
                ui_font=(self.app.ui_font_name, 10),
                small=True,
            ).grid(row=0, column=1, sticky="e")

        if not details:
            if not search.length_passing_candidates:
                message = "没有 ORF 通过当前最小长度阈值。"
                suggestion = "可以降低最小 ORF 长度，或检查输入序列是否为完整 DNA。"
            elif req.require_tata and search.tata_pass_count == 0:
                message = "没有候选同时通过长度阈值和 TATA 过滤。"
                suggestion = "关闭“只保留带 TATA 盒的候选”后再试，或检查上游 TATAAA。"
            else:
                message = "没有候选进入详细分析。"
                suggestion = "请调整最小 ORF 长度、TATA 过滤开关或输入序列。"
            empty = tk.Frame(summary, bg=COLORS["card"])
            empty.grid(row=1, column=0, sticky="ew", padx=18, pady=(5, 18))
            tk.Label(
                empty,
                text=message,
                bg=COLORS["card"],
                fg=COLORS["text"],
                font=(self.app.ui_font_name, 11, "bold"),
                anchor="w",
                justify="left",
            ).pack(fill="x")
            tk.Label(
                empty,
                text=suggestion,
                bg=COLORS["card"],
                fg=COLORS["muted"],
                font=(self.app.ui_font_name, 9),
                anchor="w",
                justify="left",
            ).pack(fill="x", pady=(3, 0))
        else:
            body = tk.Frame(summary, bg=COLORS["card"])
            body.grid(row=1, column=0, sticky="ew", padx=18, pady=(2, 16))
            body.grid_columnconfigure(0, weight=1)
            for row, detail in enumerate(details):
                if row:
                    tk.Frame(body, bg=COLORS["border"], height=1).grid(
                        row=row * 2 - 1, column=0, sticky="ew", pady=8
                    )
                block = tk.Frame(body, bg=COLORS["card"])
                block.grid(row=row * 2, column=0, sticky="ew")
                block.grid_columnconfigure(0, weight=1)
                cand = detail.candidate
                orf = cand.orf
                head = tk.Frame(block, bg=COLORS["card"])
                head.grid(row=0, column=0, sticky="ew")
                tk.Label(
                    head,
                    text=f"{cand.candidate_id}  {orf.orf_id}",
                    bg=COLORS["card"],
                    fg=COLORS["accent"],
                    font=(self.app.ui_font_name, 10, "bold"),
                ).pack(side="left")
                tk.Label(
                    head,
                    text=(
                        f"  {'正链' if orf.strand == '+' else '负链'} · 读码框 "
                        f"{orf.strand}{orf.frame} · {orf.aa_length} aa"
                    ),
                    bg=COLORS["card"],
                    fg=COLORS["muted"],
                    font=(self.app.ui_font_name, 9),
                ).pack(side="left")
                api_text = api_status_text(detail, req.use_uniprot)
                tk.Label(
                    head,
                    text=api_text,
                    bg=COLORS["card"],
                    fg=COLORS["green"] if detail.uniprot and detail.uniprot.status == "matched" else COLORS["muted"],
                    font=(self.app.ui_font_name, 9, "bold"),
                ).pack(side="right")

                info = (
                    f"输入正向坐标 {interval_1based(orf.genomic_start, orf.genomic_end)}"
                    f"   ·   TATA {len(cand.tata_hits)} 个   ·   "
                    f"AATAAA {len(cand.polya_hits)} 个   ·   "
                    f"内含子 {len(detail.splicing.introns)} 个"
                )
                tk.Label(
                    block,
                    text=info,
                    bg=COLORS["card"],
                    fg=COLORS["text"],
                    font=(self.app.ui_font_name, 9),
                    anchor="w",
                    justify="left",
                    wraplength=1020,
                ).grid(row=1, column=0, sticky="w", pady=(3, 4))

                protein = detail.translation.protein or "（无蛋白产物）"
                sequence_preview = protein if len(protein) <= 90 else protein[:90] + "…"
                tk.Label(
                    block,
                    text=f"蛋白: {sequence_preview}",
                    bg=COLORS["card"],
                    fg=COLORS["muted"],
                    font=(self.app.mono_font_name, 9),
                    anchor="w",
                    justify="left",
                    wraplength=1020,
                ).grid(row=2, column=0, sticky="w")

        # 提示/警告
        notes = list(dict.fromkeys(req.warnings))
        if notes:
            warning_box = card_frame(self.inner, bg=COLORS["amber_bg"])
            warning_box.grid(row=4, column=0, columnspan=4, sticky="ew", padx=18, pady=(0, 14))
            tk.Label(
                warning_box,
                text="输入与结果提示",
                bg=COLORS["amber_bg"],
                fg=COLORS["amber"],
                font=(self.app.ui_font_name, 10, "bold"),
                anchor="w",
            ).pack(fill="x", padx=16, pady=(12, 3))
            for note in notes:
                tk.Label(
                    warning_box,
                    text=f"•  {note}",
                    bg=COLORS["amber_bg"],
                    fg=COLORS["text"],
                    font=(self.app.ui_font_name, 9),
                    anchor="w",
                    justify="left",
                    wraplength=1020,
                ).pack(fill="x", padx=16, pady=1)
            tk.Frame(warning_box, bg=COLORS["amber_bg"], height=10).pack(fill="x")

        # TATA 说明
        note = (
            "TATAAA / AATAAA 只是可供教学观察的六碱基字面标记；"
            "本程序不构建启动子模型或 poly(A) 位点模型。"
        )
        tk.Label(
            self.inner,
            text=note,
            bg=COLORS["bg"],
            fg=COLORS["muted"],
            font=(self.app.ui_font_name, 9),
            anchor="w",
            justify="left",
            wraplength=1040,
        ).grid(row=5, column=0, columnspan=4, sticky="w", padx=20, pady=(0, 18))


class CandidatesTab(tk.Frame):
    def __init__(self, parent: tk.Misc, app: "DNAApp") -> None:
        super().__init__(parent, bg=COLORS["bg"])
        self.app = app
        self.analysis: Optional[GuiAnalysis] = None
        self._candidates: list[CandidateGene] = []
        self._detail_by_orf: dict[str, CandidateAnalysis] = {}
        self._candidate_by_iid: dict[str, CandidateGene] = {}

        self.grid_columnconfigure(0, weight=1)
        self.grid_rowconfigure(2, weight=1)
        self.grid_rowconfigure(5, weight=2)

        self.heading_var = tk.StringVar(value="所有通过最小长度阈值的候选 ORF")
        tk.Label(
            self,
            textvariable=self.heading_var,
            bg=COLORS["bg"],
            fg=COLORS["text"],
            font=(app.ui_font_name, 15, "bold"),
            anchor="w",
        ).grid(row=0, column=0, sticky="ew", padx=18, pady=(14, 3))

        self.help_var = tk.StringVar(
            value="表格按氨基酸长度排序；选中一行可在下方查看该候选的坐标与标记。双击已分析候选可跳转到剪接详情。"
        )
        tk.Label(
            self,
            textvariable=self.help_var,
            bg=COLORS["bg"],
            fg=COLORS["muted"],
            font=(app.ui_font_name, 9),
            anchor="w",
            justify="left",
            wraplength=1120,
        ).grid(row=1, column=0, sticky="ew", padx=18, pady=(0, 10))

        table_frame = tk.Frame(
            self,
            bg=COLORS["card"],
            highlightbackground=COLORS["border"],
            highlightthickness=1,
        )
        table_frame.grid(row=2, column=0, sticky="nsew", padx=18)
        table_frame.grid_rowconfigure(0, weight=1)
        table_frame.grid_columnconfigure(0, weight=1)
        columns = ("idx", "orf", "strand", "frame", "coords", "aa", "tata", "polya", "analysis")
        self.tree = ttk.Treeview(table_frame, columns=columns, show="headings", height=8)
        headings = {
            "idx": ("序号", 55),
            "orf": ("ORF", 100),
            "strand": ("链", 65),
            "frame": ("框", 60),
            "coords": ("输入正向坐标", 160),
            "aa": ("aa", 60),
            "tata": ("TATA", 75),
            "polya": ("AATAAA", 80),
            "analysis": ("详细分析", 95),
        }
        for key, (label, width) in headings.items():
            self.tree.heading(key, text=label)
            self.tree.column(key, width=width, anchor="center", stretch=(key in {"orf", "coords"}))
        self.tree.tag_configure("analyzed", background=COLORS["green_bg"])
        vsb = ttk.Scrollbar(table_frame, orient="vertical", command=self.tree.yview)
        self.tree.configure(yscrollcommand=vsb.set)
        self.tree.grid(row=0, column=0, sticky="nsew")
        vsb.grid(row=0, column=1, sticky="ns")
        self.tree.bind("<<TreeviewSelect>>", self._on_select)
        self.tree.bind("<Double-1>", self._on_double_click)

        detail_header = tk.Frame(self, bg=COLORS["bg"])
        detail_header.grid(row=4, column=0, sticky="ew", padx=18, pady=(12, 5))
        detail_header.grid_columnconfigure(0, weight=1)
        tk.Label(
            detail_header,
            text="候选详情",
            bg=COLORS["bg"],
            fg=COLORS["text"],
            font=(app.ui_font_name, 12, "bold"),
            anchor="w",
        ).grid(row=0, column=0, sticky="w")
        self.detail_button = make_button(
            detail_header,
            "打开剪接与翻译",
            self._open_selected_detail,
            ui_font=(app.ui_font_name, 10),
            small=True,
        )
        self.detail_button.grid(row=0, column=1, sticky="e")
        self.detail_button.configure(state="disabled")

        detail_container = tk.Frame(self, bg=COLORS["bg"])
        detail_container.grid(row=5, column=0, sticky="nsew", padx=18, pady=(0, 16))
        detail_container.grid_rowconfigure(0, weight=1)
        detail_container.grid_columnconfigure(0, weight=1)
        panel, self.detail_text = make_text_panel(
            detail_container,
            ui_font_name=app.ui_font_name,
            mono_font_name=app.mono_font_name,
            height=10,
        )
        panel.grid(row=0, column=0, sticky="nsew")
        set_rich_text(
            self.detail_text,
            lambda text: text.insert(
                "end",
                "还没有分析结果。输入 DNA 并开始分析后，这里会显示候选信息。\n",
                "muted",
            ),
        )

    def show_loading(self) -> None:
        self.analysis = None
        self._candidates = []
        self._detail_by_orf = {}
        self._candidate_by_iid = {}
        self.tree.delete(*self.tree.get_children())
        self.heading_var.set("正在检索候选 ORF…")
        self.detail_button.configure(state="disabled")
        set_rich_text(self.detail_text, lambda text: text.insert("end", "分析进行中…\n", "muted"))

    def render(self, analysis: GuiAnalysis) -> None:
        self.analysis = analysis
        self._candidates = analysis.search.length_passing_candidates
        self._detail_by_orf = {detail.candidate.orf.orf_id: detail for detail in analysis.details}
        self._candidate_by_iid = {}
        self.tree.delete(*self.tree.get_children())

        analyzed_count = len(analysis.details)
        total = len(self._candidates)
        self.heading_var.set(
            f"候选 ORF 汇总 · 通过长度阈值 {total} 个 · 详细分析 {analyzed_count} 个"
        )
        self.help_var.set(
            "绿色行表示进入了剪接、翻译和 UniProt 查询；点击行查看基本信息和标记坐标，"
            "双击绿色行进入该候选的剪接与翻译标签页。"
        )
        if not self._candidates:
            self._on_empty("没有 ORF 通过当前的最小长度阈值或 TATA 过滤。")
            return

        for index, candidate in enumerate(self._candidates, start=1):
            orf = candidate.orf
            iid = orf.orf_id
            self._candidate_by_iid[iid] = candidate
            detail = self._detail_by_orf.get(iid)
            tags = ("analyzed",) if detail else ()
            self.tree.insert(
                "",
                "end",
                iid=iid,
                values=(
                    index,
                    orf.orf_id,
                    "+" if orf.strand == "+" else "-",
                    f"+{orf.frame}" if orf.strand == "+" else f"-{orf.frame}",
                    interval_1based(orf.genomic_start, orf.genomic_end),
                    orf.aa_length,
                    f"是 ({len(candidate.tata_hits)})" if candidate.has_tata else "否",
                    f"是 ({len(candidate.polya_hits)})" if candidate.has_polya_signal else "否",
                    "已分析" if detail else "未分析",
                ),
                tags=tags,
            )
        children = self.tree.get_children()
        if children:
            self.tree.selection_set(children[0])
            self.tree.focus(children[0])
            self._on_select()

    def _on_empty(self, message: str) -> None:
        self.heading_var.set("候选 ORF 汇总")
        self.detail_button.configure(state="disabled")
        set_rich_text(
            self.detail_text,
            lambda text: (
                text.insert("end", "没有候选\n", "h1"),
                text.insert("end", message + "\n", "warn"),
                text.insert(
                    "end",
                    "\n建议：降低最小 ORF 长度、关闭 TATA 过滤，或检查输入是否为有效的 DNA/FASTA。\n",
                    "muted",
                ),
            ),
        )

    def _on_select(self, _event: tk.Event | None = None) -> None:
        selection = self.tree.selection()
        if not selection:
            self.detail_button.configure(state="disabled")
            return
        iid = selection[0]
        candidate = self._candidate_by_iid.get(iid)
        if candidate is None:
            return
        detail = self._detail_by_orf.get(iid)
        self.detail_button.configure(state="normal" if detail else "disabled")
        self._render_detail(candidate, detail)

    def _on_double_click(self, _event: tk.Event | None = None) -> None:
        selection = self.tree.selection()
        if not selection:
            return
        iid = selection[0]
        if iid in self._detail_by_orf:
            self.app.show_candidate_detail(iid)

    def _open_selected_detail(self) -> None:
        selection = self.tree.selection()
        if selection and selection[0] in self._detail_by_orf:
            self.app.show_candidate_detail(selection[0])

    def _render_detail(self, candidate: CandidateGene, detail: Optional[CandidateAnalysis]) -> None:
        orf = candidate.orf

        def build(text: tk.Text) -> None:
            text.insert("end", f"{orf.orf_id} · 候选区域详情\n", "h1")
            text.insert(
                "end",
                f"链: {'正链(+)' if orf.strand == '+' else '负链(-)'}    "
                f"读码框: {orf.strand}{orf.frame}    "
                f"ORF 长度: {orf.aa_length} aa\n",
                "body",
            )
            text.insert(
                "end",
                f"ORF 输入正向坐标（1-based，含终止子）: "
                f"{interval_1based(orf.genomic_start, orf.genomic_end)}\n",
                "body",
            )
            text.insert(
                "end",
                f"候选区域输入正向坐标: "
                f"{interval_1based(candidate.region_start, candidate.region_end)}    "
                f"区域长度: {candidate.region_end - candidate.region_start} bp\n",
                "body",
            )
            text.insert(
                "end",
                f"ORF 在 pre-mRNA 中的坐标: "
                f"{interval_1based(candidate.orf_pre_start, candidate.orf_pre_end)}    "
                f"起始/终止: {orf.start_codon} / {orf.stop_codon}\n",
                "body",
            )

            text.insert("end", "\nTATA 上游命中\n", "h2")
            if candidate.tata_hits:
                for index, hit in enumerate(candidate.tata_hits, start=1):
                    text.insert(
                        "end",
                        f"  T{index}: pre-mRNA {interval_1based(hit.start, hit.end)}"
                        f"  |  输入正向 {interval_1based(hit.genomic_start, hit.genomic_end)}"
                        f"  |  距 ATG 上游 {hit.distance_bp} bp\n",
                        "body",
                    )
            else:
                text.insert("end", "  未在指定上游窗口内发现 TATAAA。\n", "muted")

            text.insert("end", "\nAATAAA 下游命中\n", "h2")
            if candidate.polya_hits:
                for index, hit in enumerate(candidate.polya_hits, start=1):
                    text.insert(
                        "end",
                        f"  P{index}: pre-mRNA {interval_1based(hit.start, hit.end)}"
                        f"  |  输入正向 {interval_1based(hit.genomic_start, hit.genomic_end)}"
                        f"  |  距终止子下游 {hit.distance_bp} bp\n",
                        "body",
                    )
            else:
                text.insert("end", "  未在指定下游窗口内发现 AATAAA。\n", "muted")

            text.insert("end", "\nORF 原始翻译产物（未剪接）\n", "h2")
            text.insert("end", orf.protein + "\n", "mono")

            if detail is None:
                text.insert(
                    "end",
                    "\n该候选未进入详细分析。GUI 默认只对排序靠前的候选执行剪接、翻译和 API 查询。\n",
                    "warn",
                )
            else:
                text.insert("end", "\n下游结果摘要\n", "h2")
                text.insert(
                    "end",
                    f"  正式翻译蛋白长度: {len(detail.translation.protein)} aa\n"
                    f"  识别内含子: {len(detail.splicing.introns)} 个    "
                    f"成熟 mRNA: {len(detail.splicing.mature_mrna)} nt\n"
                    f"  UniProt: {api_status_text(detail, self.app.current_analysis.request.use_uniprot if self.app.current_analysis else True)}\n",
                    "body",
                )
                if detail.translation.warnings:
                    for warning in detail.translation.warnings:
                        text.insert("end", f"  [警告] {warning}\n", "warn")

            text.insert(
                "end",
                "\n说明：以上坐标相对本次输入序列，并非参考基因组坐标；标记仅为字面匹配。\n",
                "muted",
            )

        set_rich_text(self.detail_text, build)


class SplicingTab(tk.Frame):
    def __init__(self, parent: tk.Misc, app: "DNAApp") -> None:
        super().__init__(parent, bg=COLORS["bg"])
        self.app = app
        self._details: list[CandidateAnalysis] = []
        self._labels: list[str] = []
        self._current: Optional[CandidateAnalysis] = None

        self.grid_columnconfigure(0, weight=1)
        self.grid_rowconfigure(1, weight=1)

        header = tk.Frame(self, bg=COLORS["bg"])
        header.grid(row=0, column=0, sticky="ew", padx=18, pady=(14, 8))
        header.grid_columnconfigure(2, weight=1)
        tk.Label(
            header,
            text="剪接与翻译",
            bg=COLORS["bg"],
            fg=COLORS["text"],
            font=(app.ui_font_name, 15, "bold"),
            anchor="w",
        ).grid(row=0, column=0, sticky="w")
        tk.Label(
            header,
            text="选择候选：",
            bg=COLORS["bg"],
            fg=COLORS["muted"],
            font=(app.ui_font_name, 9),
        ).grid(row=0, column=1, sticky="w", padx=(18, 6))
        self.candidate_var = tk.StringVar(value="")
        self.combo = ttk.Combobox(
            header,
            textvariable=self.candidate_var,
            state="readonly",
            width=48,
            font=(app.ui_font_name, 9),
        )
        self.combo.grid(row=0, column=2, sticky="w")
        self.combo.bind("<<ComboboxSelected>>", self._on_combo)
        self.copy_mrna_button = make_button(
            header,
            "复制成熟 mRNA",
            self._copy_mrna,
            ui_font=(app.ui_font_name, 10),
            small=True,
        )
        self.copy_mrna_button.grid(row=0, column=3, sticky="e", padx=(8, 0))
        self.copy_protein_button = make_button(
            header,
            "复制蛋白序列",
            self._copy_protein,
            ui_font=(app.ui_font_name, 10),
            small=True,
        )
        self.copy_protein_button.grid(row=0, column=4, sticky="e", padx=(8, 0))

        self.paned = ttk.Panedwindow(self, orient="horizontal")
        self.paned.grid(row=1, column=0, sticky="nsew", padx=18, pady=(0, 16))
        left_frame = tk.Frame(self.paned, bg=COLORS["bg"])
        right_frame = tk.Frame(self.paned, bg=COLORS["bg"])
        self.paned.add(left_frame, weight=3)
        self.paned.add(right_frame, weight=2)

        left_frame.grid_rowconfigure(1, weight=1)
        left_frame.grid_columnconfigure(0, weight=1)
        tk.Label(
            left_frame,
            text="转录与剪接",
            bg=COLORS["bg"],
            fg=COLORS["muted"],
            font=(app.ui_font_name, 10, "bold"),
            anchor="w",
        ).grid(row=0, column=0, sticky="w", pady=(0, 4))
        left_panel, self.left_text = make_text_panel(
            left_frame,
            ui_font_name=app.ui_font_name,
            mono_font_name=app.mono_font_name,
            height=24,
        )
        left_panel.grid(row=1, column=0, sticky="nsew")

        right_frame.grid_rowconfigure(1, weight=1)
        right_frame.grid_columnconfigure(0, weight=1)
        tk.Label(
            right_frame,
            text="翻译结果",
            bg=COLORS["bg"],
            fg=COLORS["muted"],
            font=(app.ui_font_name, 10, "bold"),
            anchor="w",
        ).grid(row=0, column=0, sticky="w", pady=(0, 4))
        right_panel, self.right_text = make_text_panel(
            right_frame,
            ui_font_name=app.ui_font_name,
            mono_font_name=app.mono_font_name,
            height=24,
        )
        right_panel.grid(row=1, column=0, sticky="nsew")

        self._set_empty()

    def _set_empty(self, message: str = "还没有剪接结果。") -> None:
        self.combo.configure(values=[], state="disabled")
        self.candidate_var.set("")
        self.copy_mrna_button.configure(state="disabled")
        self.copy_protein_button.configure(state="disabled")
        set_rich_text(self.left_text, lambda text: text.insert("end", message + "\n", "muted"))
        set_rich_text(self.right_text, lambda text: text.insert("end", message + "\n", "muted"))

    def show_loading(self) -> None:
        self._details = []
        self._labels = []
        self._current = None
        self._set_empty("正在计算转录、剪接和翻译…")

    def render(self, analysis: GuiAnalysis) -> None:
        self._details = analysis.details
        self._labels = [
            f"{detail.candidate.candidate_id} · "
            f"{'正链' if detail.candidate.orf.strand == '+' else '负链'} "
            f"{detail.candidate.orf.strand}{detail.candidate.orf.frame} · "
            f"{detail.candidate.orf.aa_length} aa · {detail.candidate.orf.orf_id}"
            for detail in self._details
        ]
        if not self._details:
            self._set_empty("没有候选进入详细分析。")
            return
        self.combo.configure(values=self._labels, state="readonly")
        self.combo.current(0)
        self._show_detail(0)

    def _on_combo(self, _event: tk.Event | None = None) -> None:
        index = self.combo.current()
        if index >= 0:
            self._show_detail(index)

    def select_orf(self, orf_id: str) -> None:
        for index, detail in enumerate(self._details):
            if detail.candidate.orf.orf_id == orf_id:
                self.combo.current(index)
                self._show_detail(index)
                return

    def _show_detail(self, index: int) -> None:
        if index < 0 or index >= len(self._details):
            return
        detail = self._details[index]
        self._current = detail
        self.copy_mrna_button.configure(state="normal" if detail.splicing.mature_mrna else "disabled")
        self.copy_protein_button.configure(state="normal" if detail.translation.protein else "disabled")
        self._render_left(detail)
        self._render_right(detail)

    def _render_left(self, detail: CandidateAnalysis) -> None:
        cand = detail.candidate
        orf = cand.orf
        markers_tata = cand.tata_hits
        markers_polya = cand.polya_hits

        def build(text: tk.Text) -> None:
            text.insert("end", f"{cand.candidate_id} · 转录与剪接\n", "h1")
            text.insert(
                "end",
                f"ORF: {orf.orf_id}    "
                f"{'正链(+)' if orf.strand == '+' else '负链(-)'}    读码框 {orf.strand}{orf.frame}\n",
                "body",
            )
            text.insert(
                "end",
                f"输入正向坐标: {interval_1based(orf.genomic_start, orf.genomic_end)}    "
                f"候选区域: {interval_1based(cand.region_start, cand.region_end)}    "
                f"蛋白: {orf.aa_length} aa\n",
                "body",
            )
            text.insert(
                "end",
                "坐标相对于本次输入序列；负链已映射回输入正向坐标。\n",
                "muted",
            )

            text.insert("end", "简易调控标记\n", "h2")
            text.insert(
                "end",
                f"TATA 盒字面基序 {TATA_BOX_MOTIF}（RNA 中为 UAUAAA）: "
                f"{len(markers_tata)} 个上游命中\n",
                "body",
            )
            for idx, hit in enumerate(markers_tata[:20], start=1):
                text.insert(
                    "end",
                    f"  T{idx}  pre-mRNA {interval_1based(hit.start, hit.end)}"
                    f"  |  {interval_1based(hit.genomic_start, hit.genomic_end)}"
                    f"  |  距 ATG {hit.distance_bp} bp\n",
                    "muted",
                )
            if len(markers_tata) > 20:
                text.insert("end", f"  … 另有 {len(markers_tata) - 20} 个命中\n", "muted")
            text.insert(
                "end",
                f"AATAAA（RNA 中为 AAUAAA）: {len(markers_polya)} 个下游命中\n",
                "body",
            )
            for idx, hit in enumerate(markers_polya[:20], start=1):
                text.insert(
                    "end",
                    f"  P{idx}  pre-mRNA {interval_1based(hit.start, hit.end)}"
                    f"  |  {interval_1based(hit.genomic_start, hit.genomic_end)}"
                    f"  |  距终止子 {hit.distance_bp} bp\n",
                    "muted",
                )
            if len(markers_polya) > 20:
                text.insert("end", f"  … 另有 {len(markers_polya) - 20} 个命中\n", "muted")

            text.insert("end", "pre-mRNA（候选区域转录；含上游侧翼）\n", "h2")
            text.insert(
                "end",
                f"长度 {len(cand.pre_mrna)} nt。教学简化：整段候选区域一起转录，"
                "真实 pre-mRNA 从转录起始位点开始。\n",
                "muted",
            )
            insert_sequence_block(text, cand.pre_mrna, width=64)

            text.insert("end", "串联重复 / 疑似转座子粗筛\n", "h2")
            if not SCREEN_TANDEM_REPEATS:
                text.insert("end", "筛查开关已关闭。\n", "muted")
            elif detail.repeats:
                for idx, repeat in enumerate(detail.repeats, start=1):
                    text.insert(
                        "end",
                        f"  R{idx}  pre-mRNA {interval_1based(repeat.start, repeat.end)}"
                        f"（{repeat.length} nt），{repeat.unit} × 约 {repeat.copies} 次\n",
                        "body",
                    )
                if detail.removed_intervals:
                    text.insert("end", "  已按本次参数剔除标记区域后再剪接。\n", "warn")
                else:
                    text.insert("end", "  只标记，不剔除（默认设置）。\n", "muted")
            else:
                text.insert("end", "  未发现满足参数的串联重复。\n", "muted")

            text.insert("end", "剪接结果（pre-mRNA 层面）\n", "h2")
            text.insert(
                "end",
                f"用于剪接的序列长度: {len(detail.splice_input)} nt；"
                f"识别内含子: {len(detail.splicing.introns)} 个\n",
                "body",
            )
            if detail.splicing.introns:
                for idx, intron in enumerate(detail.splicing.introns, start=1):
                    text.insert(
                        "end",
                        f"  内含子 {idx}: {interval_1based(intron.start, intron.end)}"
                        f"（{intron.length} nt），{intron.donor}…{intron.acceptor}\n",
                        "body",
                    )
                text.insert(
                    "end",
                    "  外显子: "
                    + " | ".join(
                        f"E{idx} {interval_1based(start, end)}（{end - start} nt）"
                        for idx, (start, end) in enumerate(detail.splicing.exons, start=1)
                    )
                    + "\n",
                    "body",
                )
            else:
                text.insert("end", "  未发现满足长度约束的 GU…AG，pre-mRNA 原样保留。\n", "muted")
            for note in detail.splicing.notes:
                text.insert("end", f"  [提示] {note}\n", "warn")
            text.insert(
                "end",
                "  注意：GU-AG 只是简化规则，未考虑分支点、多嘧啶区、增强子和可变剪接。\n",
                "muted",
            )

            text.insert("end", "成熟 mRNA\n", "h2")
            text.insert("end", f"长度: {len(detail.splicing.mature_mrna)} nt\n", "body")
            insert_sequence_block(text, detail.splicing.mature_mrna, width=64)

        set_rich_text(self.left_text, build)

    def _render_right(self, detail: CandidateAnalysis) -> None:
        translation = detail.translation
        orf_protein = detail.candidate.orf.protein

        def build(text: tk.Text) -> None:
            text.insert("end", "标准遗传密码翻译\n", "h1")
            if translation.start_index is None:
                text.insert("end", "成熟 mRNA 中未找到 AUG，无法翻译。\n", "warn")
                for warning in translation.warnings:
                    text.insert("end", warning + "\n", "muted")
                return
            text.insert(
                "end",
                f"第一个 AUG 起点: {translation.start_index + 1}（成熟 mRNA 1-based）\n",
                "body",
            )
            if translation.stop_codon and translation.stop_codon_start is not None:
                text.insert(
                    "end",
                    f"终止密码子: {translation.stop_codon}，位置 "
                    f"{translation.stop_codon_start + 1}（1-based）\n",
                    "body",
                )
            else:
                text.insert("end", "终止密码子: 未找到（翻译可能不完整）\n", "warn")
            text.insert("end", f"蛋白长度: {len(translation.protein)} aa\n", "body")
            text.insert("end", "蛋白序列\n", "h2")
            insert_sequence_block(text, translation.protein, width=50, unit="aa")

            text.insert("end", "质量检查\n", "h2")
            if translation.protein and translation.protein != orf_protein:
                text.insert(
                    "end",
                    "成熟 mRNA 的翻译产物与原始基因组 ORF 不同；可能由剪接/重复剔除改变阅读框，"
                    "也可能原始 ORF 跨越内含子。\n",
                    "warn",
                )
            else:
                text.insert("end", "成熟 mRNA 翻译产物与原始 ORF 产物一致。\n", "good")
            if translation.complete:
                text.insert("end", "已遇到同框终止密码子并停止翻译。\n", "good")
            for warning in translation.warnings:
                text.insert("end", f"[警告] {warning}\n", "warn")

            text.insert("end", "阅读提示\n", "h2")
            text.insert(
                "end",
                "从成熟 mRNA 的第一个 AUG 开始翻译；未考虑 Kozak、uORF、IRES、"
                "稀有氨基酸、线粒体密码表或翻译后修饰。\n",
                "muted",
            )

        set_rich_text(self.right_text, build)

    def _copy_mrna(self) -> None:
        if self._current:
            self.app.copy_to_clipboard(self._current.splicing.mature_mrna, "成熟 mRNA")

    def _copy_protein(self) -> None:
        if self._current:
            self.app.copy_to_clipboard(self._current.translation.protein, "蛋白序列")


class UniProtTab(tk.Frame):
    def __init__(self, parent: tk.Misc, app: "DNAApp") -> None:
        super().__init__(parent, bg=COLORS["bg"])
        self.app = app
        self._details: list[CandidateAnalysis] = []
        self._labels: list[str] = []
        self._link_tags: list[str] = []
        self._current: Optional[CandidateAnalysis] = None

        self.grid_columnconfigure(0, weight=1)
        self.grid_rowconfigure(2, weight=1)
        header = tk.Frame(self, bg=COLORS["bg"])
        header.grid(row=0, column=0, sticky="ew", padx=18, pady=(14, 6))
        header.grid_columnconfigure(2, weight=1)
        tk.Label(
            header,
            text="UniProt 注释与相似蛋白",
            bg=COLORS["bg"],
            fg=COLORS["text"],
            font=(app.ui_font_name, 15, "bold"),
            anchor="w",
        ).grid(row=0, column=0, sticky="w")
        tk.Label(
            header,
            text="选择候选：",
            bg=COLORS["bg"],
            fg=COLORS["muted"],
            font=(app.ui_font_name, 9),
        ).grid(row=0, column=1, sticky="w", padx=(18, 6))
        self.candidate_var = tk.StringVar(value="")
        self.combo = ttk.Combobox(
            header,
            textvariable=self.candidate_var,
            state="readonly",
            width=52,
            font=(app.ui_font_name, 9),
        )
        self.combo.grid(row=0, column=2, sticky="w")
        self.combo.bind("<<ComboboxSelected>>", self._on_combo)

        note = (
            "查询流程：先做 UniProt MD5 checksum 全长精确匹配；无精确命中且开启相似性搜索时，"
            "用 GenomeNet BLASTP（Swiss-Prot）找相似蛋白，再由 UniProtKB 提供名称、物种和功能。"
            "相似结果按“一致度 × 查询覆盖度”的综合匹配度排序。点击 Accession 可打开 UniProt。"
        )
        tk.Label(
            self,
            text=note,
            bg=COLORS["bg"],
            fg=COLORS["muted"],
            font=(app.ui_font_name, 9),
            anchor="w",
            justify="left",
            wraplength=1120,
        ).grid(row=1, column=0, sticky="ew", padx=18, pady=(0, 8))

        container = tk.Frame(self, bg=COLORS["bg"])
        container.grid(row=2, column=0, sticky="nsew", padx=18, pady=(0, 16))
        container.grid_rowconfigure(0, weight=1)
        container.grid_columnconfigure(0, weight=1)
        panel, self.text = make_text_panel(
            container,
            ui_font_name=app.ui_font_name,
            mono_font_name=app.mono_font_name,
            height=24,
        )
        panel.grid(row=0, column=0, sticky="nsew")
        self._set_empty("还没有 UniProt 注释。")

    def _set_empty(self, message: str) -> None:
        self.combo.configure(values=[], state="disabled")
        self.candidate_var.set("")
        self._current = None
        set_rich_text(self.text, lambda text: text.insert("end", message + "\n", "muted"))
        self._link_tags = []

    def show_loading(self) -> None:
        self._details = []
        self._labels = []
        self._set_empty("正在查询 UniProt 或等待分析完成…")

    def render(self, analysis: GuiAnalysis) -> None:
        self._details = analysis.details
        self._labels = [
            f"{detail.candidate.candidate_id} · "
            f"{len(detail.translation.protein)} aa · "
            f"{detail.candidate.orf.orf_id}"
            for detail in self._details
        ]
        if not self._details:
            self._set_empty("没有候选进入详细分析。")
            return
        self.combo.configure(values=self._labels, state="readonly")
        self.combo.current(0)
        self._show_detail(0)

    def _on_combo(self, _event: tk.Event | None = None) -> None:
        index = self.combo.current()
        if index >= 0:
            self._show_detail(index)

    def select_orf(self, orf_id: str) -> None:
        for index, detail in enumerate(self._details):
            if detail.candidate.orf.orf_id == orf_id:
                self.combo.current(index)
                self._show_detail(index)
                return

    def _show_detail(self, index: int) -> None:
        if index < 0 or index >= len(self._details):
            return
        detail = self._details[index]
        self._current = detail
        req = self.app.current_analysis.request if self.app.current_analysis else None

        self._link_tags = []
        for tag in list(self.text.tag_names()):
            if tag.startswith("uniprot_link_"):
                try:
                    self.text.tag_delete(tag)
                except tk.TclError:
                    pass

        def build(text: tk.Text) -> None:
            text.insert(
                "end",
                f"{detail.candidate.candidate_id} · {detail.candidate.orf.orf_id} · UniProt\n",
                "h1",
            )
            if req is not None and not req.use_uniprot:
                text.insert("end", "本次分析已在输入窗口关闭 UniProt 查询。\n", "warn")
                return
            if not detail.translation.protein:
                text.insert("end", "没有蛋白产物，未发起 UniProt 查询。\n", "warn")
                return
            result = detail.uniprot
            if result is None:
                text.insert("end", "没有 UniProt 查询结果。\n", "warn")
                return

            method_text = {
                "sequence_field": "UniProt sequence 字段检索",
                "md5_checksum": "MD5 checksum 精确检索（/uniprotkb/search）",
                "genomenet_blast": "GenomeNet BLASTP 相似性搜索 + UniProtKB 注释",
            }.get(result.method, result.method or "未知")
            text.insert("end", f"检索方式: {method_text}\n", "body")
            is_similarity = result.method == "genomenet_blast"
            if result.status == "matched":
                if is_similarity:
                    text.insert(
                        "end",
                        "查询状态: 找到相似蛋白（非精确全长命中，结果为功能线索）\n",
                        "warn",
                    )
                else:
                    text.insert("end", "查询状态: 成功（全长精确匹配）\n", "good")
            text.insert("end", f"查询式: {short_query(result.query)}\n", "muted")
            if result.message:
                tag = "warn" if result.status != "matched" else "muted"
                text.insert("end", f"说明: {result.message}\n", tag)

            if result.status == "error":
                text.insert(
                    "end",
                    "\n查询失败。请检查网络、超时设置，稍后重试；也可以在输入窗口中关闭相似性 BLAST。\n",
                    "bad",
                )
                return
            if result.status == "no_match":
                text.insert(
                    "end",
                    "\n未找到精确匹配，也未找到满足阈值的相似蛋白。全新蛋白、短肽、突变体或"
                    "非标准翻译产物没有命中是正常现象；这不代表程序出错。\n",
                    "warn",
                )
                return
            if not result.hits:
                text.insert("end", "\n服务器返回 matched 状态但没有可显示条目。\n", "warn")
                return

            if is_similarity:
                text.insert(
                    "end",
                    "\n相似蛋白匹配（按综合匹配度从高到低）\n",
                    "h2",
                )
                text.insert(
                    "end",
                    "排名公式：综合匹配度 = 一致度 × 查询覆盖度。"
                    "一致度来自本地 BLOSUM62 Smith-Waterman 比对；E-value/bit score 来自 BLAST。\n",
                    "muted",
                )
                top_hit = result.hits[0]
                top_evalue = f"{top_hit.evalue:.2g}" if top_hit.evalue is not None else "—"
                text.insert(
                    "end",
                    f"推测参考: 最相似条目为 {top_hit.accession} · {top_hit.protein_name}"
                    f"（{top_hit.organism}），一致度 {top_hit.identity_percent:.1f}%、"
                    f"覆盖度 {top_hit.coverage_percent:.1f}%、"
                    f"综合匹配度 {top_hit.combined_match_percent:.1f}%、"
                    f"E-value {top_evalue}。\n",
                    "warn",
                )
            else:
                text.insert("end", f"\n精确命中 {len(result.hits)} 条\n", "h2")
            for index, hit in enumerate(result.hits, start=1):
                badge = "reviewed / Swiss-Prot" if hit.reviewed else "unreviewed / TrEMBL"
                text.insert("end", f"{index}. ", "h3")
                link_tag = f"uniprot_link_{index}"
                self._link_tags.append(link_tag)
                text.insert("end", hit.accession, ("link", link_tag))
                text.insert("end", f"   {hit.protein_name}\n", "h3")
                text.insert(
                    "end",
                    f"   物种: {hit.organism}\n",
                    "body",
                )
                text.insert(
                    "end",
                    f"   状态: {badge}"
                    + (f"   |   长度: {hit.length} aa" if hit.length is not None else "")
                    + "\n",
                    "muted",
                )
                if is_similarity:
                    evalue_text = f"{hit.evalue:.2g}" if hit.evalue is not None else "—"
                    bitscore_text = (
                        f"{hit.bitscore:.1f}" if hit.bitscore is not None else "—"
                    )
                    text.insert(
                        "end",
                        f"   一致度: {hit.identity_percent:.1f}%   |   "
                        f"查询覆盖度: {hit.coverage_percent:.1f}%   |   "
                        f"综合匹配度: {hit.combined_match_percent:.1f}%\n",
                        "good",
                    )
                    text.insert(
                        "end",
                        f"   BLAST 证据: E-value {evalue_text}   |   "
                        f"bit score {bitscore_text}"
                        + (
                            f"   |   比对长度 {hit.alignment_length} aa"
                            if hit.alignment_length is not None else ""
                        )
                        + "\n",
                        "muted",
                    )
                text.insert(
                    "end",
                    f"   UniProt: https://www.uniprot.org/uniprotkb/{hit.accession}/entry\n",
                    "muted",
                )
                function = hit.function.strip() if hit.function else "（数据库未提供 FUNCTION 注释）"
                text.insert("end", "   功能: ", "body")
                text.insert("end", function + "\n", "body")
                text.insert("end", "\n")

                url = f"https://www.uniprot.org/uniprotkb/{hit.accession}/entry"
                text.tag_bind(link_tag, "<Button-1>", lambda _event, target=url: webbrowser.open(target))
                text.tag_bind(link_tag, "<Enter>", lambda _event: text.configure(cursor="hand2"))
                text.tag_bind(link_tag, "<Leave>", lambda _event: text.configure(cursor="xterm"))

            if is_similarity:
                text.insert(
                    "end",
                    "注意：相似性命中是由 BLAST 和本地比对推断出来的可能相关蛋白，不是全长精确匹配；"
                    "其 FUNCTION 注释只能作为线索，不能证明本 DNA 在细胞中表达该蛋白。\n",
                    "warn",
                )
            else:
                text.insert(
                    "end",
                    "注意：UniProt 精确命中只说明数据库中存在相同序列，不能证明本 DNA 在细胞中表达该蛋白。\n",
                    "warn",
                )

        set_rich_text(self.text, build)


class RawTab(tk.Frame):
    def __init__(self, parent: tk.Misc, app: "DNAApp") -> None:
        super().__init__(parent, bg=COLORS["bg"])
        self.app = app
        self.grid_columnconfigure(0, weight=1)
        self.grid_rowconfigure(2, weight=1)
        header = tk.Frame(self, bg=COLORS["bg"])
        header.grid(row=0, column=0, sticky="ew", padx=18, pady=(14, 6))
        tk.Label(
            header,
            text="原始序列与输入信息",
            bg=COLORS["bg"],
            fg=COLORS["text"],
            font=(app.ui_font_name, 15, "bold"),
            anchor="w",
        ).pack(side="left")
        self.copy_button = make_button(
            header,
            "复制原始序列",
            self._copy_raw,
            ui_font=(app.ui_font_name, 10),
            small=True,
        )
        self.copy_button.pack(side="right")
        self.summary_var = tk.StringVar(value="等待输入…")
        tk.Label(
            self,
            textvariable=self.summary_var,
            bg=COLORS["bg"],
            fg=COLORS["muted"],
            font=(app.ui_font_name, 9),
            anchor="w",
            justify="left",
            wraplength=1120,
        ).grid(row=1, column=0, sticky="ew", padx=18, pady=(0, 8))
        container = tk.Frame(self, bg=COLORS["bg"])
        container.grid(row=2, column=0, sticky="nsew", padx=18, pady=(0, 16))
        container.grid_rowconfigure(0, weight=1)
        container.grid_columnconfigure(0, weight=1)
        panel, self.text = make_text_panel(
            container,
            ui_font_name=app.ui_font_name,
            mono_font_name=app.mono_font_name,
            height=24,
        )
        panel.grid(row=0, column=0, sticky="nsew")
        self._dna = ""
        self.show_loading()

    def show_loading(self) -> None:
        self._dna = ""
        self.summary_var.set("等待分析完成…")
        self.copy_button.configure(state="disabled")
        set_rich_text(self.text, lambda text: text.insert("end", "分析完成后显示原始 DNA 和反向互补链。\n", "muted"))

    def render(self, analysis: GuiAnalysis) -> None:
        req = analysis.request
        self._dna = req.dna
        self.summary_var.set(
            f"来源: {req.source}   ·   长度: {len(req.dna)} bp   ·   "
            f"GC: {analysis.gc_percent:.2f}%   ·   "
            f"最小 ORF: {req.min_aa} aa"
        )
        self.copy_button.configure(state="normal")

        def build(text: tk.Text) -> None:
            text.insert("end", "原始 DNA（5'→3'）\n", "h1")
            insert_sequence_block(text, req.dna, width=80)
            text.insert("end", "\n反向互补链（5'→3'）\n", "h2")
            insert_sequence_block(text, reverse_complement_dna(req.dna), width=80)
            if req.warnings:
                text.insert("end", "\n输入提示\n", "h2")
                for warning in req.warnings:
                    text.insert("end", f"• {warning}\n", "warn")
            text.insert(
                "end",
                "\n坐标均为相对本次输入序列的位置，不是参考基因组坐标。\n",
                "muted",
            )

        set_rich_text(self.text, build)

    def _copy_raw(self) -> None:
        if self._dna:
            self.app.copy_to_clipboard(self._dna, "原始 DNA")


# ====================================================================
# 五、主窗口和线程编排
# ====================================================================

class DNAApp:
    def __init__(self, root: tk.Tk) -> None:
        self.root = root
        self.ui_font_name = choose_first_font(
            root,
            ["Microsoft YaHei UI", "Microsoft YaHei", "PingFang SC", "Noto Sans CJK SC", "Segoe UI"],
            "TkDefaultFont",
        )
        self.mono_font_name = choose_first_font(
            root,
            ["Cascadia Mono", "Consolas", "Menlo", "DejaVu Sans Mono", "Courier New"],
            "TkFixedFont",
        )
        self.current_analysis: Optional[GuiAnalysis] = None
        self._busy = False
        self._queue: queue.Queue = queue.Queue()
        self._worker_thread: Optional[threading.Thread] = None
        self._input_dialog: Optional[SequenceInputDialog] = None
        self._initial_open = True
        self._job_after_id: Optional[str] = None

        self.root.title(APP_TITLE)
        self.root.configure(bg=COLORS["bg"])
        self.root.minsize(1020, 660)
        screen_width = self.root.winfo_screenwidth()
        screen_height = self.root.winfo_screenheight()
        initial_width = max(1020, min(1200, screen_width - 120))
        initial_height = max(660, min(780, screen_height - 120))
        self._configure_styles()
        self._build_main_window()
        center_window(self.root, initial_width, initial_height)
        self.root.protocol("WM_DELETE_WINDOW", self._on_close)
        self.root.bind("<Control-n>", lambda _event: self.open_input_dialog())
        self.root.bind("<Control-e>", lambda _event: self.export_report())
        self.root.report_callback_exception = self._report_tk_exception

    # ---------- 样式 ----------
    def _configure_styles(self) -> None:
        style = ttk.Style(self.root)
        try:
            style.theme_use("clam")
        except tk.TclError:
            pass
        style.configure(
            ".",
            font=(self.ui_font_name, 10),
            background=COLORS["bg"],
            foreground=COLORS["text"],
        )
        style.configure("TFrame", background=COLORS["bg"])
        style.configure("TNotebook", background=COLORS["bg"], borderwidth=0)
        style.configure(
            "TNotebook.Tab",
            background="#dce6f3",
            foreground=COLORS["muted"],
            padding=(18, 9),
            font=(self.ui_font_name, 10, "bold"),
        )
        style.map(
            "TNotebook.Tab",
            background=[("selected", COLORS["card"])],
            foreground=[("selected", COLORS["accent"])],
        )
        style.configure(
            "Treeview",
            background=COLORS["card"],
            fieldbackground=COLORS["card"],
            foreground=COLORS["text"],
            font=(self.ui_font_name, 9),
            rowheight=28,
            borderwidth=0,
        )
        style.configure(
            "Treeview.Heading",
            font=(self.ui_font_name, 9, "bold"),
            background="#e8eef8",
            foreground=COLORS["text"],
            padding=(6, 7),
            relief="flat",
        )
        style.map("Treeview", background=[("selected", "#dbeafe")], foreground=[("selected", COLORS["text"])])
        style.configure(
            "Vertical.TScrollbar",
            background="#cbd5e1",
            troughcolor=COLORS["bg"],
            borderwidth=0,
            arrowsize=14,
        )
        style.configure(
            "TProgressbar",
            troughcolor="#dce6f3",
            background=COLORS["accent"],
            borderwidth=0,
            thickness=8,
        )
        style.configure("TSpinbox", padding=4)

    # ---------- 主窗口 ----------
    def _build_main_window(self) -> None:
        self.root.grid_rowconfigure(3, weight=1)
        self.root.grid_columnconfigure(0, weight=1)

        header = tk.Frame(self.root, bg=COLORS["header"], padx=24, pady=16)
        header.grid(row=0, column=0, sticky="ew")
        header.grid_columnconfigure(0, weight=1)
        tk.Label(
            header,
            text=APP_TITLE,
            bg=COLORS["header"],
            fg="#ffffff",
            font=(self.ui_font_name, 20, "bold"),
            anchor="w",
        ).grid(row=0, column=0, sticky="w")
        tk.Label(
            header,
            text="ORF 检索 · pre-mRNA 剪接 · 翻译 · UniProt 注释 · 相似蛋白 BLAST   |   教学模拟，不用于正式科研",
            bg=COLORS["header"],
            fg=COLORS["header_muted"],
            font=(self.ui_font_name, 9),
            anchor="w",
        ).grid(row=1, column=0, sticky="w", pady=(3, 0))

        actions = tk.Frame(header, bg=COLORS["header"])
        actions.grid(row=0, column=1, rowspan=2, sticky="e")
        self.new_button = make_button(
            actions,
            "＋ 新分析",
            self.open_input_dialog,
            ui_font=(self.ui_font_name, 10),
            dark=True,
        )
        self.new_button.pack(side="left", padx=(0, 8))
        self.export_button = make_button(
            actions,
            "导出报告",
            self.export_report,
            ui_font=(self.ui_font_name, 10),
            primary=True,
        )
        self.export_button.pack(side="left")
        self.export_button.configure(state="disabled")

        meta = tk.Frame(self.root, bg=COLORS["card"], padx=22, pady=10)
        meta.grid(row=1, column=0, sticky="ew")
        meta.grid_columnconfigure(0, weight=1)
        self.source_var = tk.StringVar(value="等待输入 DNA 序列…")
        self.params_var = tk.StringVar(value="先点击“＋ 新分析”或等待启动输入窗口。")
        self.status_var = tk.StringVar(value="空闲")
        tk.Label(
            meta,
            textvariable=self.source_var,
            bg=COLORS["card"],
            fg=COLORS["text"],
            font=(self.ui_font_name, 10, "bold"),
            anchor="w",
            justify="left",
            wraplength=930,
        ).grid(row=0, column=0, sticky="w")
        tk.Label(
            meta,
            textvariable=self.params_var,
            bg=COLORS["card"],
            fg=COLORS["muted"],
            font=(self.ui_font_name, 9),
            anchor="w",
            justify="left",
            wraplength=930,
        ).grid(row=1, column=0, sticky="w", pady=(2, 0))
        tk.Label(
            meta,
            textvariable=self.status_var,
            bg=COLORS["card"],
            fg=COLORS["accent"],
            font=(self.ui_font_name, 9, "bold"),
            anchor="e",
        ).grid(row=0, column=1, rowspan=2, sticky="e", padx=(18, 0))

        self.progress = ttk.Progressbar(self.root, mode="indeterminate", style="TProgressbar")
        self.progress.grid(row=2, column=0, sticky="new")
        self.progress.grid_remove()

        self.notebook = ttk.Notebook(self.root)
        self.notebook.grid(row=3, column=0, sticky="nsew", padx=14, pady=(10, 8))
        self.root.grid_rowconfigure(3, weight=1)

        self.overview_tab = OverviewTab(self.notebook, self)
        self.candidates_tab = CandidatesTab(self.notebook, self)
        self.splicing_tab = SplicingTab(self.notebook, self)
        self.uniprot_tab = UniProtTab(self.notebook, self)
        self.raw_tab = RawTab(self.notebook, self)
        self.notebook.add(self.overview_tab, text="概览")
        self.notebook.add(self.candidates_tab, text="ORF 候选")
        self.notebook.add(self.splicing_tab, text="剪接与翻译")
        self.notebook.add(self.uniprot_tab, text="UniProt / 相似蛋白")
        self.notebook.add(self.raw_tab, text="原始序列")

        footer = tk.Frame(self.root, bg="#e3ebf5", padx=20, pady=8)
        footer.grid(row=4, column=0, sticky="ew")
        footer.grid_columnconfigure(0, weight=1)
        tk.Label(
            footer,
            text=(
                "本程序为教学简化模拟器；ORF、TATA/AATAAA、剪接、UniProt 精确检索和"
                "相似蛋白 BLAST 推断结果均不能作为正式科研/临床/生产结论。"
            ),
            bg="#e3ebf5",
            fg=COLORS["muted"],
            font=(self.ui_font_name, 8),
            anchor="w",
            justify="left",
            wraplength=1160,
        ).grid(row=0, column=0, sticky="w")

        self.overview_tab.show_loading("等你在输入窗口粘贴 DNA 并点击“开始分析”。")
        self.candidates_tab.show_loading()
        self.splicing_tab.show_loading()
        self.uniprot_tab.show_loading()
        self.raw_tab.show_loading()

    # ---------- 输入与开始 ----------
    def open_input_dialog(self, initial: bool = False) -> None:
        if self._busy:
            return
        if self._input_dialog is not None and self._input_dialog.winfo_exists():
            self._input_dialog.lift()
            self._input_dialog.focus_force()
            return
        self._input_dialog = SequenceInputDialog(
            self.root,
            ui_font_name=self.ui_font_name,
            mono_font_name=self.mono_font_name,
            on_submit=self._start_analysis,
            on_cancel=self._input_cancelled,
            initial=initial,
        )

    def _input_cancelled(self) -> None:
        self._input_dialog = None
        if self._initial_open and self.current_analysis is None:
            self.root.destroy()
        else:
            self.status_var.set("已取消新分析")

    def _start_analysis(self, request: GuiRequest) -> None:
        self._input_dialog = None
        self._initial_open = False
        self.root.deiconify()
        self.root.lift()
        self.root.focus_force()
        self.current_analysis = None
        self._busy = True
        self.new_button.configure(state="disabled")
        self.export_button.configure(state="disabled")
        self.source_var.set(f"来源: {request.source}")
        self.params_var.set(
            f"最小 ORF {request.min_aa} aa   ·   "
            f"TATA 过滤 {'开' if request.require_tata else '关'}   ·   "
            f"剔除重复 {'开' if request.remove_repeats else '关'}   ·   "
            f"UniProt {'开' if request.use_uniprot else '关'}   ·   "
            f"相似性 BLAST {'开' if request.use_similarity and request.use_uniprot else '关'}"
        )
        self.status_var.set("正在分析…")
        self.progress.grid()
        self.progress.start(12)
        self.notebook.select(self.overview_tab)
        self.overview_tab.show_loading("正在扫描 ORF 并准备转录…")
        self.candidates_tab.show_loading()
        self.splicing_tab.show_loading()
        self.uniprot_tab.show_loading()
        self.raw_tab.show_loading()

        while True:
            try:
                self._queue.get_nowait()
            except queue.Empty:
                break

        self._worker_thread = threading.Thread(
            target=self._worker_main,
            args=(request,),
            daemon=True,
            name="DNAAnalysisWorker",
        )
        self._worker_thread.start()
        self._schedule_poll()

    def _worker_main(self, request: GuiRequest) -> None:
        try:
            analysis = run_gui_analysis(
                request,
                progress=lambda message: self._queue.put(("progress", message)),
            )
            self._queue.put(("done", analysis))
        except Exception as exc:
            self._queue.put(
                (
                    "error",
                    {
                        "message": f"{type(exc).__name__}: {exc}",
                        "traceback": traceback.format_exc(),
                    },
                )
            )

    def _schedule_poll(self) -> None:
        self._job_after_id = self.root.after(100, self._poll_worker)

    def _poll_worker(self) -> None:
        done = False
        try:
            while True:
                kind, payload = self._queue.get_nowait()
                if kind == "progress":
                    self.status_var.set(str(payload))
                elif kind == "done":
                    done = True
                    self._finish_analysis(payload)
                    break
                elif kind == "error":
                    done = True
                    self._fail_analysis(payload["message"], payload["traceback"])
                    break
        except queue.Empty:
            pass
        if self._busy and not done:
            self._schedule_poll()

    def _stop_busy(self) -> None:
        self._busy = False
        try:
            self.progress.stop()
        except tk.TclError:
            pass
        self.progress.grid_remove()
        self.new_button.configure(state="normal")
        self.export_button.configure(state="normal" if self.current_analysis else "disabled")

    def _finish_analysis(self, analysis: GuiAnalysis) -> None:
        self.current_analysis = analysis
        self._busy = False
        self._stop_busy()
        status = (
            f"分析完成 · 用时 {analysis.elapsed_seconds:.2f} 秒"
            if analysis.details or analysis.search.all_complete_orfs
            else "分析完成 · 没有完整 ORF"
        )
        self.status_var.set(status)
        self.overview_tab.render(analysis)
        self.candidates_tab.render(analysis)
        self.splicing_tab.render(analysis)
        self.uniprot_tab.render(analysis)
        self.raw_tab.render(analysis)
        self.notebook.select(self.overview_tab)
        self.root.lift()

    def _fail_analysis(self, message: str, details: str) -> None:
        self._busy = False
        self._stop_busy()
        self.status_var.set("分析失败")
        self.overview_tab.show_error(message, details)
        self.notebook.select(self.overview_tab)
        messagebox.showerror(
            "分析失败",
            f"{message}\n\n可查看“概览”页中的错误详情，或检查输入/网络后重试。",
            parent=self.root,
        )

    # ---------- 交互 ----------
    def show_first_detail(self) -> None:
        if not self.current_analysis or not self.current_analysis.details:
            return
        self.show_candidate_detail(self.current_analysis.details[0].candidate.orf.orf_id)

    def show_candidate_detail(self, orf_id: str) -> None:
        if not self.current_analysis:
            return
        self.splicing_tab.select_orf(orf_id)
        self.notebook.select(self.splicing_tab)

    def select_uniprot_candidate(self, orf_id: str) -> None:
        self.uniprot_tab.select_orf(orf_id)
        self.notebook.select(self.uniprot_tab)

    def copy_to_clipboard(self, content: str, label: str) -> None:
        if not content:
            return
        self.root.clipboard_clear()
        self.root.clipboard_append(content)
        self.status_var.set(f"已复制{label}（{len(content)} 个字符）")

    def export_report(self) -> None:
        if not self.current_analysis:
            messagebox.showinfo("导出报告", "请先完成一次分析。", parent=self.root)
            return
        path = filedialog.asksaveasfilename(
            parent=self.root,
            title="导出分析报告",
            defaultextension=".txt",
            filetypes=[("文本报告", "*.txt"), ("所有文件", "*.*")],
            initialfile="DNA序列分析报告.txt",
        )
        if not path:
            return
        try:
            report = build_text_report(self.current_analysis)
            with open(path, "w", encoding="utf-8-sig", newline="\n") as handle:
                handle.write(report)
            self.status_var.set(f"报告已导出: {path}")
            messagebox.showinfo("导出成功", f"分析报告已保存：\n{path}", parent=self.root)
        except Exception as exc:
            messagebox.showerror("导出失败", str(exc), parent=self.root)

    def _on_close(self) -> None:
        if self._busy:
            if not messagebox.askyesno(
                "退出",
                "分析仍在进行，确定要退出吗？\n（后台线程会随窗口关闭而结束。）",
                parent=self.root,
            ):
                return
        if self._job_after_id is not None:
            try:
                self.root.after_cancel(self._job_after_id)
            except tk.TclError:
                pass
        self.root.destroy()

    def _report_tk_exception(self, exc_type, exc_value, exc_tb) -> None:
        details = "".join(traceback.format_exception(exc_type, exc_value, exc_tb))
        print(details, file=sys.stderr)
        try:
            messagebox.showerror("界面错误", f"{exc_type.__name__}: {exc_value}", parent=self.root)
        except tk.TclError:
            pass


# ====================================================================
# 六、启动
# ====================================================================

def enable_high_dpi() -> None:
    if sys.platform != "win32":
        return
    try:
        ctypes.windll.shcore.SetProcessDpiAwareness(1)  # type: ignore[attr-defined]
        return
    except Exception:
        pass
    try:
        ctypes.windll.user32.SetProcessDPIAware()  # type: ignore[attr-defined]
    except Exception:
        pass


def launch_gui() -> int:
    enable_high_dpi()
    root = tk.Tk()
    root.withdraw()
    app = DNAApp(root)

    def open_first_dialog() -> None:
        if root.winfo_exists():
            app.open_input_dialog(initial=True)

    root.after(120, open_first_dialog)
    root.mainloop()
    return 0


def main() -> int:
    try:
        return launch_gui()
    except tk.TclError as exc:
        print(
            "无法启动图形界面（可能是当前环境没有可用的桌面显示）："
            f"{exc}\n"
            "在有桌面环境的 Windows 上运行：py -3 dna_gui.py\n"
            "命令行版本仍可使用：py -3 dna_sequence_analysis.py --demo --skip-uniprot",
            file=sys.stderr,
        )
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
