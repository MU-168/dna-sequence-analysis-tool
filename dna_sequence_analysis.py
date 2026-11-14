#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
DNA 序列分析教学模拟器（分析内核 / 命令行模式）
=====================================================================

图形前端见同目录 dna_gui.py；无参数运行本文件时也会自动打开该前端。

流水线
------
原始 DNA（正向链或 FASTA）
  → 自动生成反向互补链，六读码框检索“ATG…终止密码子”完整 ORF
  → 用 TATAAA 与 AATAAA 作为简易调控信号辅助标记（只作教学）
  → 提取候选区域并在基因方向上转录为 pre-mRNA（T → U）
  → 在 pre-mRNA 上按 RNA 的 GU-AG（即 DNA 的 GT-AG）边界做简化内含子剔除
  → 拼接外显子得到成熟 mRNA
  → 用 NCBI/标准遗传密码表（table 1）从成熟 mRNA 的第一个 AUG 翻译，
    遇到第一个同框终止密码子停止
  → 将蛋白序列提交 UniProt 公开 REST API，显示蛋白名称、物种和 FUNCTION 注释

依赖安装
--------
    Python 3.9+
    pip install biopython
    图形界面使用 Python 标准库 tkinter，通常无需额外安装。

运行
----
图形界面（推荐；需要桌面环境，tkinter 为 Python 标准库）:
    python dna_sequence_analysis.py          # 自动打开图形界面和输入窗口
    python dna_gui.py                        # 直接启动图形界面

命令行模式:
    python dna_sequence_analysis.py --demo                        # 内置测试序列
    python dna_sequence_analysis.py --demo --skip-uniprot         # 不发网络请求
    python dna_sequence_analysis.py --demo --no-similarity        # 只做精确检索
    python dna_sequence_analysis.py --min-aa 1 --dna "ATGAAATAA"
    python dna_sequence_analysis.py --fasta example.fa
    python dna_sequence_analysis.py --fasta example.fa --record-id seq1
    python dna_sequence_analysis.py --dna "YOUR_DNA" --min-aa 50 --require-tata
    python dna_sequence_analysis.py --dna "YOUR_DNA" --remove-repeats

重点声明（务必阅读）
--------------------
1. 这不是基因预测软件。六读码框扫描只寻找 ATG 起始、同框终止的完整 ORF；
   没有启动子模型、Kozak 序列、剪接位点统计模型、密码子使用偏好、转录本
   证据或物种参数。AUGUSTUS、BRAKER、GENSCAN、GeneMark 等专业基因预测器
   使用隐马尔可夫模型/机器学习、剪接位点模型、长度分布和训练集，可以同时
   预测多外显子结构；本脚本与它们不是同一类工具，结果不能相提并论。
2. TATAAA 只是 TATA 盒的六碱基字面匹配，不检查与转录起始位点的距离、方向、
   染色质状态或其他启动子元件；AATAAA 也只是 poly(A) 信号的简化字面标记，
   不能证明转录终止或 poly(A) 加尾真实发生。真实调控元件远比这两个六联体复杂。
3. 在未剪接的 DNA 上扫描 ORF、然后再剪接，是一种教学用简化流程。真实真核基因
   的内含子会打断外显子阅读框，因此“先找 ORF 再剪接”不能可靠复原真实基因结构。
   本脚本仅用于演示模块之间的数据流。
4. 剪接模块只寻找 GU（DNA 中为 GT）和下游 AG，并施加最小/最大内含子长度和
   最小外显子长度。它忽略分支点、多嘧啶区、剪接增强子/沉默子、可变剪接、
   GC-AG 和 AT-AC 等次要剪接类型。凡满足简单规则的 GU…AG 都可能被误认成
   内含子，因此会漏剪、错剪或过度剪接；绝不能用于正式基因结构注释。
5. “转座子筛查”实际上只是寻找短串联重复（如 ACAA×6），仅标记“疑似重复/
   疑似转座子样区域”。串联重复不等于转座子；没有重复序列数据库、结构域、
   末端重复或转座酶信息。默认只标记、不剔除，只有配置开关打开时才强行删除。
6. 翻译从成熟 mRNA 的第一个 AUG 开始，使用标准遗传密码，未考虑 Kozak 上下文、
   上游 ORF、IRES、硒代半胱氨酸、稀有氨基酸、线粒体遗传密码、RNA 编辑、
   信号肽切割、二硫键、糖基化、磷酸化或其他翻译后修饰。遇到同框终止密码子
   即停止；若没有终止密码子，会翻译到末端并给出警告。
7. UniProt 查询只能返回数据库中已有条目的注释。全新蛋白、突变体、短肽或
   非天然序列通常没有精确命中；即使命中，数据库注释也不等于本序列经过实验验证。
   注意：UniProt 当前 `/uniprotkb/search` 端点不接受 `sequence:` 作为检索字段
   （服务器返回 400：“'sequence' is not a valid search field”）。因此程序默认
   直接使用同一端点支持且已经实测可用的 `checksum:<md5> + length` 精确检索。
   若将来 UniProt 恢复 sequence 字段，可把 UNIPROT_USE_SEQUENCE_FIELD_FIRST
   设为 True 做兼容探测。
8. 无精确命中且开启相似性搜索时，程序会把翻译后的蛋白序列提交给 GenomeNet
   的公开 BLASTP 服务（Swiss-Prot 数据库，https://www.genome.jp/tools/blast/），
   再把命中的 accession 送回 UniProtKB 获取名称、物种和 FUNCTION 注释。
   相似性搜索是外部 Web 服务，可能被限流、排队或改变页面格式；结果按本地
   Smith-Waterman 计算的“一致度 × 查询覆盖度”综合匹配度排序，只是可能相关的
   功能线索，不是全长精确匹配，也不能证明该 DNA 真实表达该蛋白。原始 DNA
   不会上传，只有蛋白序列会发送给 UniProt 和（可关闭的）GenomeNet。
   该公开网页服务仅适合低频教学查询，请勿批量、循环提交；正式用途应改用
   本地 BLAST+、UniProt BLAST 或获得许可的专业服务。
9. 为了模块演示简单，本程序把候选区域（包含启动子侧翼）整段转录为
   pre-mRNA。真实 pre-mRNA 从转录起始位点（TSS）开始，启动子 DNA 本身
   不会被转录；因此输出 pre-mRNA 中出现的 UAUAAA/AAUAAA 只是 DNA 字面
   标记的 RNA 投影，不应被理解为 RNA 层面的调控元件。
10. 绝对不要把本程序的任何结果用于论文、临床诊断、药物设计、育种决策或任何
    正式科研/生产用途。专业分析应使用经过验证的基因预测、比对、结构域和
    进化分析工具，并由具备生物学背景的人员审查。

内置演示序列
------------
合成一小段“启动子-5'UTR-内含子-编码区-3'UTR”教学构建体。ORF 编码人
metallothionein-1H（UniProt P80294，61 aa），以便演示 UniProt 精确匹配；
其 UTR、内含子和重复区域是人工设计的，不代表任何真实基因组位点。
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import socket
import sys
import textwrap
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable, Optional
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode
from urllib.request import Request, urlopen

from Bio import SeqIO
from Bio.Align import PairwiseAligner, substitution_matrices
from Bio.Data import CodonTable
from Bio.Seq import Seq


# ====================================================================
# 一、集中配置区（可直接修改；命令行可覆盖部分参数）
# ====================================================================

# ORF 与候选筛选 ------------------------------------------------------
MIN_ORF_AA = 30                 # 最小 ORF 氨基酸数（不含终止密码子）
REQUIRE_TATA_BOX = False        # True: 只保留 ORF 上游带 TATAAA 的候选
MAX_CANDIDATES_TO_ANALYZE = 5   # 最多对多少个候选执行完整下游流程/API 查询
MAX_ORF_TABLE_ROWS = 20         # 候选汇总表最多显示多少行

# 简易调控信号 --------------------------------------------------------
TATA_BOX_MOTIF = "TATAAA"       # 教学用 TATA 盒字面标记（DNA 形式）
POLY_A_SIGNAL_MOTIF = "AATAAA"  # 教学用 poly(A) 信号字面标记（DNA 形式）
TATA_SEARCH_UPSTREAM_BP = 220   # 只在 ATG 上游该距离内统计 TATAAA
POLY_A_SEARCH_DOWNSTREAM_BP = 140  # 只在终止密码子下游该距离内统计 AATAAA
MAX_REGULATORY_SITES_TO_PRINT = 5  # 每类调控信号最多列出的条数

# 候选区域与序列显示 --------------------------------------------------
# 这两个侧翼长度决定转录出的 pre-mRNA 包含多少 5'/3' 侧序列。
CANDIDATE_UPSTREAM_FLANK = 220
CANDIDATE_DOWNSTREAM_FLANK = 140
MAX_SEQUENCE_PRINT_BP = 3000    # 超过此长度时仅打印序列首尾片段
SEQUENCE_LINE_WIDTH = 80

# 标准遗传密码 --------------------------------------------------------
GENETIC_CODE_TABLE_ID = 1       # NCBI 标准遗传密码（table 1）
START_CODON_DNA = "ATG"

# pre-mRNA 剪接（教学简化） ------------------------------------------
# 在 RNA 层面识别 GU…AG；对应 DNA 记法即 GT…AG。
SPLICE_DONOR_RNA = "GU"
SPLICE_ACCEPTOR_RNA = "AG"
SPLICE_MIN_INTRON_LEN = 40      # 最短内含子长度（nt）
SPLICE_MAX_INTRON_LEN = 5000    # 最长内含子长度（nt）
SPLICE_MIN_EXON_LEN = 10        # 每个外显子至少保留的 nt 数
SPLICE_MAX_INTRONS = 8          # 每次最多剔除多少个内含子

# 串联重复/疑似转座子粗筛 --------------------------------------------
SCREEN_TANDEM_REPEATS = True
REPEAT_MIN_UNIT = 2             # 重复单元最短长度（nt）
REPEAT_MAX_UNIT = 12            # 重复单元最长长度（nt）
REPEAT_MIN_COPIES = 4           # 至少连续重复多少次
REPEAT_MIN_TOTAL_BP = 20        # 重复区间总长度下限
REMOVE_REPEAT_REGIONS = False   # 默认只标记，不剔除；命令行可打开
MAX_REPEAT_REGIONS_TO_PRINT = 10

# UniProt REST API ----------------------------------------------------
UNIPROT_BASE_URL = "https://rest.uniprot.org"
UNIPROT_SEARCH_PATH = "/uniprotkb/search"
UNIPROT_TIMEOUT_SECONDS = 20.0
UNIPROT_MIN_INTERVAL_SECONDS = 1.0   # 相邻 HTTP 请求的最小间隔（限流）
UNIPROT_MAX_RETRIES = 2              # 429/5xx/网络瞬时错误的重试次数
UNIPROT_PAGE_SIZE = 5                # 每次最多取回 5 条
UNIPROT_MAX_HITS_TO_PRINT = 8
UNIPROT_FUNCTION_MAX_CHARS = 700
# 实测：当前 UniProt REST 不接受 sequence 检索字段，会返回 HTTP 400。
# 默认直接使用同一 /uniprotkb/search 端点支持且可用的 MD5 checksum 精确检索。
# 如果以后官方恢复 sequence 字段，可改为 True 让客户端先尝试该字段。
UNIPROT_USE_SEQUENCE_FIELD_FIRST = False
# 避免把很长的蛋白序列塞进 GET URL；开启 sequence 兼容探测时，超长蛋白仍走 checksum。
UNIPROT_SEQUENCE_FIELD_MAX_AA = 500

# 无精确命中时的相似蛋白检索 -----------------------------------------
# 方案：GenomeNet BLASTP（Swiss-Prot 数据库，公开 Web 服务，无需邮箱）
#       → 取得相似条目的 UniProt accession
#       → 从 UniProtKB 获取名称、物种和 FUNCTION 注释
#       → 本地对查询与命中序列做 Smith-Waterman 比对，计算一致度/覆盖度。
# 注意：这是相似性推断，不是精确命中，也不能证明真实功能。
UNIPROT_SIMILARITY_ENABLED = True
UNIPROT_BLAST_ENDPOINT = "https://www.genome.jp/tools-bin/blastplus"
UNIPROT_BLAST_DATABASE = "swissprot"
UNIPROT_BLAST_TIMEOUT_SECONDS = 90.0
UNIPROT_BLAST_MIN_INTERVAL_SECONDS = 2.0
UNIPROT_SIMILARITY_MAX_HITS = 8
UNIPROT_SIMILARITY_FETCH_LIMIT = 12
UNIPROT_SIMILARITY_MIN_IDENTITY = 20.0    # 本地比对一致度下限（%）
UNIPROT_SIMILARITY_MIN_COVERAGE = 20.0    # 查询序列覆盖度下限（%）
UNIPROT_BLAST_FILTER = "F"                # GenomeNet BLAST 低复杂度过滤：F=关闭，T=开启

# 内置教学用 DNA（合成构建体；ORF 编码人 metallothionein-1H，P80294）---
# 结构：启动子(含 TATAAA) + 5'UTR 外显子1 + GU…AG 内含子 + 5'UTR 外显子2
#       + 编码区(ATG…TAA) + 3'UTR(含 AATAAA 与 ACAA 串联重复)。
DEMO_DNA_SEQUENCE = (
    # 5' 上游：A/T-rich 合成启动子（仅一个 TATAAA）、外显子 1、GU…AG 内含子、外显子 2 的 5'UTR
    "AATTTAATAATTTTATATTTTATTAAAAATTATATTTATTTTATTTATAATTAATTATAAATATTTAAATTAAATATAT"
    "TAATAATTTAAATATTGTATTTTAATTATAATTTAATTAATAATAATTAATATATAATTTATTATTTAAATTTTATATTT"
    "AATAAAATAAGATATTTTTT"
    # 主体 ORF：ATG…TAA；编码人 metallothionein-1H（P80294）61 aa
    "ATGGATCCTAATTGCTCTTGCGAAGCTGGCGGCTCTTGCGCTTGCGCTGGCTCTTGCAAATGCAAAAAATGCAAATGCA"
    "CTTCTTGCAAAAAATCTTGCTGCTCTTGCTGCCCTTTAGGCTGCGCTAAATGCGCTCAAGGCTGCATTTGCAAAGGCGC"
    "TTCTGAAAAATGCTCTTGCTGCGCTTAA"
    # 3'UTR：含一个 AATAAA 与 ACAA×6 串联重复
    "TATAAATAAAAAAATATTTTTTATACAAACAAACAAACAAACAAACAATTTTTTATTAT"
)
DEMO_SOURCE_LABEL = (
    "内置演示序列（合成启动子/UTR/内含子；ORF 编码人 metallothionein-1H）"
)


# ====================================================================
# 二、数据结构
# ====================================================================

@dataclass
class RuntimeOptions:
    """本次运行可被命令行覆盖的选项。"""
    min_orf_aa: int = MIN_ORF_AA
    require_tata_box: bool = REQUIRE_TATA_BOX
    remove_repeat_regions: bool = REMOVE_REPEAT_REGIONS
    max_candidates: int = MAX_CANDIDATES_TO_ANALYZE


@dataclass
class ORFRecord:
    """一个完整 ORF（ATG 起始、含终止密码子）。坐标均为 0-based 半开区间。"""
    orf_id: str
    strand: str                    # "+" 或 "-"
    frame: int                     # 1/2/3，相对于该链 5'→3' 方向的读码框编号
    genomic_start: int             # 正向基因组链上的 0-based inclusive 起点
    genomic_end: int               # 正向基因组链上的 0-based exclusive 终点（含终止子）
    oriented_start: int            # 在该基因方向序列中的 0-based 起点
    oriented_end: int              # 在该基因方向序列中的 0-based 终点（含终止子）
    protein: str                   # 不含终止密码子的翻译产物
    start_codon: str
    stop_codon: str

    @property
    def aa_length(self) -> int:
        return len(self.protein)

    @property
    def nt_length(self) -> int:
        return self.genomic_end - self.genomic_start


@dataclass
class MotifHit:
    """一个在候选 pre-mRNA/候选 DNA 方向上的字面调控基序。"""
    kind: str                      # "TATA" 或 "polyA"
    motif_dna: str
    start: int                     # pre-mRNA 0-based inclusive
    end: int                       # pre-mRNA 0-based exclusive
    genomic_start: int             # 正向基因组链 0-based inclusive
    genomic_end: int               # 正向基因组链 0-based exclusive
    distance_bp: int               # 与 ORF 起点的距离（TATA）或终点距离（polyA）


@dataclass
class CandidateGene:
    """一个待分析的候选基因区域。"""
    orf: ORFRecord
    region_start: int              # 正向基因组链上的候选区域起点（0-based）
    region_end: int                # 正向基因组链上的候选区域终点（exclusive）
    oriented_dna: str              # 按基因方向排列的候选区域 DNA
    pre_mrna: str                  # 转录后的 pre-mRNA（U 代替 T）
    orf_pre_start: int             # ORF 在 pre-mRNA 中的 0-based 起点
    orf_pre_end: int               # ORF 在 pre-mRNA 中的 0-based 终点（含终止子）
    tata_hits: list[MotifHit] = field(default_factory=list)
    polya_hits: list[MotifHit] = field(default_factory=list)
    candidate_id: str = ""

    @property
    def has_tata(self) -> bool:
        return bool(self.tata_hits)

    @property
    def has_polya_signal(self) -> bool:
        return bool(self.polya_hits)


@dataclass
class GeneSearchResult:
    """ORF 检索与候选筛选的汇总结果。"""
    all_complete_orfs: list[ORFRecord]
    length_passing_candidates: list[CandidateGene]
    selected_candidates: list[CandidateGene]
    frame_counts: dict[str, int]
    tata_pass_count: int
    tata_filter_enabled: bool


@dataclass
class TandemRepeatRegion:
    """粗筛得到的疑似串联重复区域。"""
    start: int                     # pre-mRNA/输入 RNA 的 0-based inclusive
    end: int                       # 0-based exclusive
    unit: str                      # 代表性重复单元
    copies: int                    # 代表性单元的连续拷贝数估计

    @property
    def length(self) -> int:
        return self.end - self.start


@dataclass
class SpliceIntron:
    """在 pre-mRNA 坐标中记录一个被剔除的内含子。"""
    start: int                     # 0-based inclusive（GU 的 G）
    end: int                       # 0-based exclusive（AG 之后的下一碱基）
    sequence: str
    donor: str
    acceptor: str

    @property
    def length(self) -> int:
        return self.end - self.start


@dataclass
class SplicingResult:
    input_mrna: str
    introns: list[SpliceIntron]
    exons: list[tuple[int, int]]
    mature_mrna: str
    notes: list[str] = field(default_factory=list)


@dataclass
class TranslationResult:
    protein: str
    start_index: Optional[int]
    stop_codon_start: Optional[int]
    stop_codon: Optional[str]
    cds_rna: str
    complete: bool
    warnings: list[str] = field(default_factory=list)


@dataclass
class UniProtHit:
    accession: str
    entry_id: str
    protein_name: str
    organism: str
    function: str
    reviewed: bool
    annotation_score: Optional[float]
    length: Optional[int]
    # 匹配类型与相似性指标。精确命中时一致度/覆盖度按 100% 处理。
    match_type: str = "exact"       # "exact" / "similar"
    identity_percent: Optional[float] = None
    coverage_percent: Optional[float] = None
    combined_match_percent: Optional[float] = None
    evalue: Optional[float] = None
    bitscore: Optional[float] = None
    alignment_length: Optional[int] = None


@dataclass
class UniProtResult:
    status: str                    # "matched" / "no_match" / "error"
    method: str = ""
    query: str = ""
    hits: list[UniProtHit] = field(default_factory=list)
    message: str = ""


# ====================================================================
# 三、通用工具与输入处理
# ====================================================================

def setup_console() -> None:
    """尽量让 Windows 控制台的 stdout/stderr 都以 UTF-8 输出中文。"""
    for stream in (sys.stdout, sys.stderr):
        try:
            if hasattr(stream, "reconfigure"):
                stream.reconfigure(encoding="utf-8", errors="replace")  # type: ignore[attr-defined]
        except Exception:
            # 控制台编码无法改变时继续，最坏情况只会出现替换字符。
            pass


def sanitize_dna(raw: str) -> tuple[str, list[str]]:
    """
    清理并校验 DNA 字符串。

    允许空白字符（会被移除）和 A/C/G/T/N。RNA 的 U 不在接受范围内，因为
    本程序要求输入 DNA，并由程序自己完成 T→U 转录。遇到其他 IUPAC 简并碱基
    会直接报错，避免后续统计和翻译出现含糊结果。
    """
    if raw is None:
        raise ValueError("DNA 序列为空。")
    if ">" in raw:
        raise ValueError("检测到 FASTA 标题符号 '>'。若为 FASTA 文件，请使用 --fasta。")
    seq = "".join(str(raw).split()).upper()
    if not seq:
        raise ValueError("DNA 序列为空。")
    if "U" in seq:
        raise ValueError("输入含 U，看起来像 RNA。请提供 DNA，程序会自动转录。")
    invalid = sorted(set(seq) - set("ACGTN"))
    if invalid:
        raise ValueError(
            "序列包含不支持的字符：" + ", ".join(repr(c) for c in invalid)
            + "。当前仅接受 A/C/G/T/N。"
        )
    if len(seq) < 3:
        raise ValueError("DNA 至少需要 3 nt 才能进行密码子级分析。")
    warnings: list[str] = []
    if "N" in seq:
        warnings.append("序列含 N；含 N 的密码子会翻译为 X，通常不会被当作起始/终止密码子。")
    return seq, warnings


def load_fasta_sequence(path: str, record_id: Optional[str] = None) -> tuple[str, str, list[str]]:
    """
    读取 FASTA。默认使用第一条记录；给出 record_id 时选择对应记录。
    返回 (DNA序列, 来源描述, 警告列表)。
    """
    fasta_path = Path(path)
    if not fasta_path.exists():
        raise ValueError(f"FASTA 文件不存在：{path}")
    try:
        records = list(SeqIO.parse(str(fasta_path), "fasta"))
    except Exception as exc:
        raise ValueError(f"FASTA 解析失败：{exc}") from exc
    if not records:
        raise ValueError("FASTA 中没有解析到任何记录。")

    warnings: list[str] = []
    if record_id is not None:
        chosen = None
        for record in records:
            if record.id == record_id or record.description.split()[0] == record_id:
                chosen = record
                break
        if chosen is None:
            available = ", ".join(record.id for record in records[:10])
            raise ValueError(f"未找到 FASTA 记录 {record_id!r}；可用记录：{available}")
    else:
        chosen = records[0]
        if len(records) > 1:
            warnings.append(
                f"FASTA 含 {len(records)} 条记录，默认只分析第一条（{chosen.id}）；"
                "可用 --record-id 指定其他记录。"
            )

    seq, seq_warnings = sanitize_dna(str(chosen.seq))
    warnings.extend(seq_warnings)
    header = chosen.description if chosen.description else chosen.id
    return seq, f"FASTA: {header}", warnings


def reverse_complement_dna(dna: str) -> str:
    """用 Biopython 生成反向互补链。"""
    return str(Seq(dna).reverse_complement())


def transcribe_to_pre_mrna(dna: str) -> str:
    """用 Biopython 完成 DNA→RNA 的 T→U 转录。"""
    rna = str(Seq(dna).transcribe())
    return rna.upper()


def wrap_text(text: str, width: int = 100, initial_indent: str = "", subsequent_indent: str = "") -> str:
    """按显示宽度折行，保留中文。"""
    return textwrap.fill(
        text,
        width=width,
        initial_indent=initial_indent,
        subsequent_indent=subsequent_indent,
        break_long_words=True,
        break_on_hyphens=False,
    )


def print_section(title: str) -> None:
    print()
    print("=" * 78)
    print(title)
    print("=" * 78)


def print_subsection(title: str) -> None:
    print()
    print("-" * 78)
    print(title)
    print("-" * 78)


def print_sequence_block(
    label: str,
    sequence: str,
    max_chars: int = MAX_SEQUENCE_PRINT_BP,
    unit: str = "nt",
) -> None:
    """
    打印带坐标的序列。坐标从 1 开始，便于人工核对。
    长序列只显示首尾，避免控制台被淹没。
    """
    print(f"  {label}（{len(sequence)} {unit}）:")
    if len(sequence) <= max_chars:
        shown_positions = [(0, len(sequence))]
        truncated = False
    else:
        head_len = max_chars // 2
        tail_len = max_chars - head_len
        shown_positions = [(0, head_len), (len(sequence) - tail_len, len(sequence))]
        truncated = True

    for block_start, block_end in shown_positions:
        if block_start > 0:
            print("    ...（中间序列省略）...")
        for offset in range(block_start, block_end, SEQUENCE_LINE_WIDTH):
            chunk = sequence[offset:min(offset + SEQUENCE_LINE_WIDTH, block_end)]
            print(f"    {offset + 1:>6}  {chunk}")
    if truncated:
        print(
            f"    （仅显示前/后约 {max_chars} {unit}；"
            "修改 MAX_SEQUENCE_PRINT_BP 可显示全文。）"
        )


def strand_label(strand: str) -> str:
    return "正链(+)" if strand == "+" else "负链(-)"


def interval_1based(start: int, end: int) -> str:
    """把 0-based 半开区间 [start,end) 显示为 1-based 闭区间。"""
    return f"{start + 1}-{end}"


def short_query_for_display(query: str, limit: int = 130) -> str:
    if len(query) <= limit:
        return query
    return query[:limit - 20] + "…" + query[-15:]


# ====================================================================
# 四、六读码框 ORF 检索
# ====================================================================

def find_complete_orfs_in_six_frames(
    dna: str,
    genetic_code: int = GENETIC_CODE_TABLE_ID,
) -> list[ORFRecord]:
    """
    在正向链和反向互补链的各 3 个读码框中寻找完整 ORF。

    定义（本教学程序）：
    - 必须以 ATG 开始；
    - 必须以同框的 TAA/TAG/TGA 结束；
    - 终止密码子计入 ORF 核苷酸区间，但不翻译进蛋白；
    - 同一个“起始/终止间隔”中若存在多个 ATG，只取最靠 5' 的那个，
      以得到该区间最长的完整 ORF（嵌套的短 ORF 不重复输出）；
    - 没有终止密码子的 3' 端不完整 ORF 不输出。

    这一步不使用 TATA 盒、不需要 TATA 盒，也不考虑内含子。
    """
    table = CodonTable.unambiguous_dna_by_id[genetic_code]
    stop_codons = set(table.stop_codons)
    dna_plus = dna.upper()
    dna_minus = reverse_complement_dna(dna_plus)
    total_len = len(dna_plus)
    records: list[ORFRecord] = []

    for strand, oriented in (("+", dna_plus), ("-", dna_minus)):
        for frame_index in range(3):
            current_start: Optional[int] = None
            for pos in range(frame_index, len(oriented) - 2, 3):
                codon = oriented[pos:pos + 3]
                if current_start is None:
                    if codon == START_CODON_DNA:
                        current_start = pos
                    continue

                if codon in stop_codons:
                    oriented_start = current_start
                    oriented_end = pos + 3
                    coding_dna = oriented[oriented_start:pos]
                    protein = translate_dna_codons(coding_dna, table)
                    if strand == "+":
                        genomic_start, genomic_end = oriented_start, oriented_end
                    else:
                        # 反向链坐标映射回正向基因组坐标。
                        genomic_start = total_len - oriented_end
                        genomic_end = total_len - oriented_start

                    records.append(
                        ORFRecord(
                            orf_id=f"ORF-{len(records) + 1:03d}",
                            strand=strand,
                            frame=frame_index + 1,
                            genomic_start=genomic_start,
                            genomic_end=genomic_end,
                            oriented_start=oriented_start,
                            oriented_end=oriented_end,
                            protein=protein,
                            start_codon=START_CODON_DNA,
                            stop_codon=codon,
                        )
                    )
                    current_start = None

    return records


def translate_dna_codons(coding_dna: str, table: CodonTable.CodonTable) -> str:
    """翻译一段不含终止密码子的 DNA；模糊密码子记为 X。"""
    amino_acids: list[str] = []
    for pos in range(0, len(coding_dna) - 2, 3):
        codon = coding_dna[pos:pos + 3].replace("U", "T")
        amino_acids.append(table.forward_table.get(codon, "X"))
    return "".join(amino_acids)


# ====================================================================
# 五、候选区域提取与调控信号标记
# ====================================================================

def find_motif_positions(sequence: str, motif: str) -> list[tuple[int, int]]:
    """返回所有重叠匹配的 (start,end)；motif 按字面匹配，不具概率含义。"""
    positions: list[tuple[int, int]] = []
    motif = motif.upper()
    start = 0
    while True:
        hit = sequence.find(motif, start)
        if hit == -1:
            break
        positions.append((hit, hit + len(motif)))
        start = hit + 1
    return positions


def transcript_to_genomic_interval(
    region_start: int,
    region_end: int,
    strand: str,
    start: int,
    end: int,
) -> tuple[int, int]:
    """把候选区域方向（pre-mRNA 坐标）的区间映射回正向基因组坐标。"""
    if strand == "+":
        return region_start + start, region_start + end
    # 候选区域是正向基因组片段的 RC：pre-mRNA 的 0 对应该基因组片段的末端，
    # 因此 [start, end) 映射为 [region_end - end, region_end - start)。
    return region_end - end, region_end - start


def build_candidate_gene(dna: str, orf: ORFRecord) -> CandidateGene:
    """
    围绕一个 ORF 截取上游/下游侧翼，并在基因方向上转录得到 pre-mRNA。
    同时统计 ATG 上游窗口内的 TATAAA 和终止子下游窗口内的 AATAAA。

    注意：这些命中只是位置上的字面匹配，不能说明启动子或 poly(A) 位点真实存在。
    """
    seq_len = len(dna)
    if orf.strand == "+":
        region_start = max(0, orf.genomic_start - CANDIDATE_UPSTREAM_FLANK)
        region_end = min(seq_len, orf.genomic_end + CANDIDATE_DOWNSTREAM_FLANK)
        oriented_dna = dna[region_start:region_end]
        orf_offset = orf.genomic_start - region_start
        orf_end_offset = orf.genomic_end - region_start
    else:
        # 对负链基因而言，转录方向的上游在正向基因组的较大坐标一侧。
        region_start = max(0, orf.genomic_start - CANDIDATE_DOWNSTREAM_FLANK)
        region_end = min(seq_len, orf.genomic_end + CANDIDATE_UPSTREAM_FLANK)
        oriented_dna = reverse_complement_dna(dna[region_start:region_end])
        orf_offset = region_end - orf.genomic_end
        orf_end_offset = region_end - orf.genomic_start

    pre_mrna = transcribe_to_pre_mrna(oriented_dna)
    tata_hits: list[MotifHit] = []
    polya_hits: list[MotifHit] = []

    for start, end in find_motif_positions(oriented_dna, TATA_BOX_MOTIF):
        if end <= orf_offset:
            distance = orf_offset - end
            if distance <= TATA_SEARCH_UPSTREAM_BP:
                gs, ge = transcript_to_genomic_interval(
                    region_start, region_end, orf.strand, start, end
                )
                tata_hits.append(
                    MotifHit(
                        kind="TATA",
                        motif_dna=TATA_BOX_MOTIF,
                        start=start,
                        end=end,
                        genomic_start=gs,
                        genomic_end=ge,
                        distance_bp=distance,
                    )
                )

    for start, end in find_motif_positions(oriented_dna, POLY_A_SIGNAL_MOTIF):
        if start >= orf_end_offset:
            distance = start - orf_end_offset
            if distance <= POLY_A_SEARCH_DOWNSTREAM_BP:
                gs, ge = transcript_to_genomic_interval(
                    region_start, region_end, orf.strand, start, end
                )
                polya_hits.append(
                    MotifHit(
                        kind="polyA",
                        motif_dna=POLY_A_SIGNAL_MOTIF,
                        start=start,
                        end=end,
                        genomic_start=gs,
                        genomic_end=ge,
                        distance_bp=distance,
                    )
                )

    return CandidateGene(
        orf=orf,
        region_start=region_start,
        region_end=region_end,
        oriented_dna=oriented_dna,
        pre_mrna=pre_mrna,
        orf_pre_start=orf_offset,
        orf_pre_end=orf_end_offset,
        tata_hits=tata_hits,
        polya_hits=polya_hits,
    )


def search_gene_candidates(dna: str, options: RuntimeOptions) -> GeneSearchResult:
    """
    执行六读码框 ORF 检索、最小长度过滤、可选 TATA 过滤、排序和候选编号。
    """
    all_orfs = find_complete_orfs_in_six_frames(dna)
    frame_counts = {
        f"+{frame}": 0 for frame in range(1, 4)
    }
    frame_counts.update({f"-{frame}": 0 for frame in range(1, 4)})
    for orf in all_orfs:
        frame_counts[f"{orf.strand}{orf.frame}"] += 1

    length_passing_candidates = [
        build_candidate_gene(dna, orf)
        for orf in all_orfs
        if orf.aa_length >= options.min_orf_aa
    ]

    # 长度过滤后再按 TATA 开关筛选；TATA 绝不参与 ORF 的初始定义。
    if options.require_tata_box:
        retained = [candidate for candidate in length_passing_candidates if candidate.has_tata]
    else:
        retained = list(length_passing_candidates)

    retained.sort(
        key=lambda candidate: (
            -candidate.orf.aa_length,
            candidate.orf.genomic_start,
            candidate.orf.strand,
            candidate.orf.frame,
        )
    )
    tata_pass_count = sum(1 for candidate in length_passing_candidates if candidate.has_tata)
    selected = retained[:options.max_candidates]
    for index, candidate in enumerate(selected, start=1):
        candidate.candidate_id = f"C{index}"

    return GeneSearchResult(
        all_complete_orfs=all_orfs,
        length_passing_candidates=sorted(
            length_passing_candidates,
            key=lambda candidate: (
                -candidate.orf.aa_length,
                candidate.orf.genomic_start,
                candidate.orf.strand,
            ),
        ),
        selected_candidates=selected,
        frame_counts=frame_counts,
        tata_pass_count=tata_pass_count,
        tata_filter_enabled=options.require_tata_box,
    )


# ====================================================================
# 六、串联重复/疑似转座子粗筛
# ====================================================================

def _merge_intervals(intervals: Iterable[tuple[int, int]]) -> list[tuple[int, int]]:
    """合并重叠或首尾相接的 0-based 半开区间。"""
    ordered = sorted((start, end) for start, end in intervals if end > start)
    if not ordered:
        return []
    merged: list[tuple[int, int]] = [ordered[0]]
    for start, end in ordered[1:]:
        last_start, last_end = merged[-1]
        if start <= last_end:
            merged[-1] = (last_start, max(last_end, end))
        else:
            merged.append((start, end))
    return merged


def find_tandem_repeats(
    sequence: str,
    min_unit: int = REPEAT_MIN_UNIT,
    max_unit: int = REPEAT_MAX_UNIT,
    min_copies: int = REPEAT_MIN_COPIES,
    min_total_bp: int = REPEAT_MIN_TOTAL_BP,
) -> list[TandemRepeatRegion]:
    """
    极简串联重复扫描：对每个单元长度，寻找至少 min_copies 次连续重复。
    这不是转座子检测，也不构建重复序列数据库。

    实现细节：从每个位置尝试 unit，若连续拷贝数足够则记录并向尾部跳转，
    减少同一重复区的重复报告；最后合并重叠区间。
    """
    sequence = sequence.upper()
    n = len(sequence)
    raw_hits: list[TandemRepeatRegion] = []
    for unit_len in range(min_unit, max_unit + 1):
        i = 0
        min_span = max(min_total_bp, unit_len * min_copies)
        while i <= n - min_span:
            unit = sequence[i:i + unit_len]
            if "N" in unit:
                i += 1
                continue
            copies = 1
            next_pos = i + unit_len
            while next_pos + unit_len <= n and sequence[next_pos:next_pos + unit_len] == unit:
                copies += 1
                next_pos += unit_len
            if copies >= min_copies and copies * unit_len >= min_total_bp:
                start = i
                # 向前回溯，确保报告的重复区间完整。
                while start - unit_len >= 0 and sequence[start - unit_len:start] == unit:
                    start -= unit_len
                    copies += 1
                end = start + copies * unit_len
                raw_hits.append(
                    TandemRepeatRegion(start=start, end=end, unit=unit, copies=copies)
                )
                i = max(i + 1, end)
            else:
                i += 1

    # 去重并按区间合并；代表单元选择跨度最大的原始命中。
    unique: dict[tuple[int, int, str], TandemRepeatRegion] = {}
    for hit in raw_hits:
        unique[(hit.start, hit.end, hit.unit)] = hit
    ordered = sorted(unique.values(), key=lambda h: (h.start, -(h.end - h.start)))
    if not ordered:
        return []

    merged: list[TandemRepeatRegion] = []
    for hit in ordered:
        if not merged:
            merged.append(hit)
            continue
        previous = merged[-1]
        if hit.start <= previous.end:
            merged_start = previous.start
            merged_end = max(previous.end, hit.end)
            # 保留跨度更大者作为代表性单元，但区间始终使用合并后的边界。
            representative = hit if hit.length > previous.length else previous
            merged[-1] = TandemRepeatRegion(
                start=merged_start,
                end=merged_end,
                unit=representative.unit,
                copies=representative.copies,
            )
        else:
            merged.append(hit)
    return merged


def remove_marked_repeats(
    sequence: str, repeats: list[TandemRepeatRegion]
) -> tuple[str, list[tuple[int, int]]]:
    """按配置剔除重复标记区域；返回 (新序列, 被剔除的合并区间)。"""
    merged = _merge_intervals((repeat.start, repeat.end) for repeat in repeats)
    if not merged:
        return sequence, []
    pieces: list[str] = []
    cursor = 0
    for start, end in merged:
        pieces.append(sequence[cursor:start])
        cursor = end
    pieces.append(sequence[cursor:])
    return "".join(pieces), merged


# ====================================================================
# 七、pre-mRNA 层面 GU-AG 简化剪接
# ====================================================================

def splice_pre_mrna(
    pre_mrna: str,
    min_intron_len: int = SPLICE_MIN_INTRON_LEN,
    max_intron_len: int = SPLICE_MAX_INTRON_LEN,
    min_exon_len: int = SPLICE_MIN_EXON_LEN,
    max_introns: int = SPLICE_MAX_INTRONS,
) -> SplicingResult:
    """
    在 RNA 上做极简的 GU…AG 内含子剔除。

    规则：
    - 供体为 GU（对应 DNA 的 GT），受体为下游 AG；
    - 内含子长度在 [min_intron_len, max_intron_len]；
    - 内含子两侧都保留至少 min_exon_len 个外显子碱基；
    - 从左到右扫描；对每个供体选择第一个满足长度约束的受体；
    - 一个内含子被选中后，从受体之后继续寻找下一个供体；
    - 最多剔除 max_introns 个。

    这不是剪接位点预测模型。没有分支点、多嘧啶区、剪接增强子或可变剪接。
    真实 RNA 中大量 GU/AG 可能被错误配对，因此结果仅用于教学演示。
    """
    rna = pre_mrna.upper()
    if "T" in rna:
        # pre-mRNA 应是 RNA；程序内部从转录步骤得到。若外部调用传入 DNA，
        # 给出明确错误而不是悄悄混淆 DNA/RNA 层。
        raise ValueError("splice_pre_mrna 需要 RNA（U），检测到 T。请先转录。")

    introns: list[SpliceIntron] = []
    notes: list[str] = []
    search_from = 0
    exon_start = 0

    while len(introns) < max_introns:
        donor = rna.find(SPLICE_DONOR_RNA, search_from)
        if donor == -1:
            break
        if donor - exon_start < min_exon_len:
            search_from = donor + 1
            continue

        # 内含子长度 = (acceptor + 2) - donor，因此有效受体最早可出现在这里。
        acceptor_search_from = donor + max(2, min_intron_len - 2)
        acceptor: Optional[int] = None
        while True:
            candidate_acceptor = rna.find(SPLICE_ACCEPTOR_RNA, acceptor_search_from)
            if candidate_acceptor == -1:
                break
            intron_len = candidate_acceptor + 2 - donor
            if intron_len > max_intron_len:
                break
            if intron_len >= min_intron_len:
                if len(rna) - (candidate_acceptor + 2) < min_exon_len:
                    # 更靠后的受体只会让最后的外显子更短，没必要继续。
                    break
                acceptor = candidate_acceptor
                break
            acceptor_search_from = candidate_acceptor + 1

        if acceptor is None:
            search_from = donor + 1
            continue

        introns.append(
            SpliceIntron(
                start=donor,
                end=acceptor + 2,
                sequence=rna[donor:acceptor + 2],
                donor=rna[donor:donor + 2],
                acceptor=rna[acceptor:acceptor + 2],
            )
        )
        exon_start = acceptor + 2
        search_from = acceptor + 2

    if len(introns) >= max_introns:
        remaining_donor = rna.find(SPLICE_DONOR_RNA, search_from)
        if remaining_donor != -1:
            notes.append(f"已达到 SPLICE_MAX_INTRONS={max_introns}，后续位点未处理。")

    exons: list[tuple[int, int]] = []
    cursor = 0
    for intron in introns:
        if cursor < intron.start:
            exons.append((cursor, intron.start))
        cursor = intron.end
    if cursor < len(rna):
        exons.append((cursor, len(rna)))

    if introns:
        mature_mrna = "".join(rna[start:end] for start, end in exons)
    else:
        mature_mrna = rna
        notes.append("未找到满足当前长度约束的 GU…AG 内含子，pre-mRNA 原样保留。")

    return SplicingResult(
        input_mrna=rna,
        introns=introns,
        exons=exons,
        mature_mrna=mature_mrna,
        notes=notes,
    )


# ====================================================================
# 八、成熟 mRNA 翻译（标准遗传密码，终止密码子停止）
# ====================================================================

def translate_mature_mrna(
    rna: str,
    genetic_code: int = GENETIC_CODE_TABLE_ID,
) -> TranslationResult:
    """
    从成熟 mRNA 的第一个 AUG 开始，使用标准遗传密码逐密码子翻译。
    遇到第一个同框终止密码子即停止，终止密码子不进入蛋白序列。

    未实现 Kozak 上下文、uORF 评估、IRES、翻译移码等真实翻译调控。
    """
    sequence = rna.upper()
    table = CodonTable.unambiguous_rna_by_id[genetic_code]
    stop_codons = set(table.stop_codons)
    start = sequence.find("AUG")
    if start == -1:
        return TranslationResult(
            protein="",
            start_index=None,
            stop_codon_start=None,
            stop_codon=None,
            cds_rna="",
            complete=False,
            warnings=["成熟 mRNA 中未找到 AUG 起始密码子，无法翻译。"],
        )

    amino_acids: list[str] = []
    warnings: list[str] = []
    ambiguous_count = 0
    stop_start: Optional[int] = None
    stop_codon: Optional[str] = None
    pos = start
    while pos + 3 <= len(sequence):
        codon = sequence[pos:pos + 3]
        if codon in stop_codons:
            stop_start = pos
            stop_codon = codon
            break
        amino_acid = table.forward_table.get(codon)
        if amino_acid is None:
            amino_acid = "X"
            ambiguous_count += 1
        amino_acids.append(amino_acid)
        pos += 3

    if stop_start is None:
        warnings.append(
            "成熟 mRNA 中未发现同框终止密码子；已翻译到序列末端，结果可能不完整。"
        )
        cds_end = (len(sequence) // 3) * 3
        if cds_end > start:
            # 只保留从起点开始完整的密码子，尾端不足 3 nt 的部分忽略。
            last_complete = start + ((len(sequence) - start) // 3) * 3
            cds_rna = sequence[start:last_complete]
        else:
            cds_rna = ""
        complete = False
    else:
        cds_rna = sequence[start:stop_start + 3]
        complete = True

    if ambiguous_count:
        warnings.append(f"翻译过程中有 {ambiguous_count} 个模糊密码子被记为 X。")

    return TranslationResult(
        protein="".join(amino_acids),
        start_index=start,
        stop_codon_start=stop_start,
        stop_codon=stop_codon,
        cds_rna=cds_rna,
        complete=complete,
        warnings=warnings,
    )


# ====================================================================
# 九、UniProt REST API 对接
# ====================================================================

class UniProtError(RuntimeError):
    """UniProt 查询类错误的基类。"""


class UniProtAPIError(UniProtError):
    """服务器返回非 2xx 状态。"""

    def __init__(self, status: int, detail: str):
        self.status = status
        self.detail = detail
        super().__init__(f"HTTP {status}: {detail}")


class UniProtNetworkError(UniProtError):
    """网络连接、DNS、读超时等异常。"""


class SimilaritySearchError(UniProtError):
    """相似性搜索服务失败、超时或结果无法解析。"""


class UniProtClient:
    """
    UniProtKB /uniprotkb/search 查询客户端。

    限流：
    - 所有 HTTP 请求至少间隔 UNIPROT_MIN_INTERVAL_SECONDS；
    - 每个蛋白通常只需一次 /uniprotkb/search 请求，候选数也有限；
    - 结果在本次运行内按蛋白序列缓存。

    检索策略：
    1. 实测当前 UniProt REST 的 /uniprotkb/search 不支持 sequence 检索字段，
       会返回 HTTP 400；因此默认直接使用该端点支持的
       `checksum:<md5>` + length 精确检索，不发出必败请求。
       checksum 字段按完整序列的 MD5 精确匹配，避免对大型数据库盲目翻页。
    2. 保留可选的 sequence 字段兼容探测：将
       UNIPROT_USE_SEQUENCE_FIELD_FIRST 设为 True 后，程序仍会先尝试
       `sequence:"<完整氨基酸序列>"`，失败再回退到 checksum。
    3. 返回条目仍会在客户端逐个核对 sequence.value == 输入蛋白，避免
       校验和碰撞或服务端模糊匹配造成误报。
    4. 精确未命中且 enable_similarity=True 时，使用 GenomeNet BLASTP
       （Swiss-Prot）获取相似蛋白候选，再通过 UniProtKB 获取注释。
       本地用 BLOSUM62 Smith-Waterman 比对计算一致度与覆盖度，并按
       “一致度 × 查询覆盖度”排序。该功能依赖外部公开 Web 服务，可能失败；
       结果属于相似性推断，不是实验验证。
    5. 捕获超时、网络、HTTP、JSON 解析、BLAST 页面格式变化和无匹配等情况，
       返回结构化状态。
    """

    def __init__(
        self,
        base_url: str = UNIPROT_BASE_URL,
        search_path: str = UNIPROT_SEARCH_PATH,
        timeout: float = UNIPROT_TIMEOUT_SECONDS,
        min_interval: float = UNIPROT_MIN_INTERVAL_SECONDS,
        max_retries: int = UNIPROT_MAX_RETRIES,
        page_size: int = UNIPROT_PAGE_SIZE,
        try_sequence_field_first: bool = UNIPROT_USE_SEQUENCE_FIELD_FIRST,
        enable_similarity: bool = UNIPROT_SIMILARITY_ENABLED,
    ) -> None:
        self.base_url = base_url.rstrip("/")
        self.search_path = search_path
        self.timeout = timeout
        self.min_interval = max(0.0, min_interval)
        self.max_retries = max(0, max_retries)
        self.page_size = max(1, page_size)
        self.try_sequence_field_first = try_sequence_field_first
        self.enable_similarity = enable_similarity
        self._last_request_started = 0.0
        self._last_blast_started = 0.0
        self._sequence_field_supported: Optional[bool] = None
        self._cache: dict[str, UniProtResult] = {}
        self._aligner: Optional[PairwiseAligner] = None

    def _throttle(self) -> None:
        """确保相邻请求之间达到最小时间间隔。"""
        now = time.monotonic()
        wait_seconds = self.min_interval - (now - self._last_request_started)
        if wait_seconds > 0:
            time.sleep(wait_seconds)

    @staticmethod
    def _error_detail_from_body(body: bytes) -> str:
        text = body.decode("utf-8", errors="replace").strip()
        if not text:
            return "服务器未返回错误详情"
        try:
            parsed = json.loads(text)
        except json.JSONDecodeError:
            return text[:400]
        if isinstance(parsed, dict):
            messages = parsed.get("messages")
            if isinstance(messages, list):
                return " | ".join(str(item) for item in messages)
            if "message" in parsed:
                return str(parsed["message"])
        return text[:400]

    def _request_json(self, query: str, page_size: Optional[int] = None) -> dict:
        """
        发起一次限流、带超时/重试的 GET 请求并解析 JSON。
        仅对整个检索请求重试，不会对每条结果额外发请求。
        """
        params = {
            "query": query,
            "format": "json",
            "size": str(page_size or self.page_size),
            "fields": (
                "accession,id,protein_name,organism_name,cc_function,"
                "length,reviewed,sequence"
            ),
        }
        url = f"{self.base_url}{self.search_path}?{urlencode(params)}"
        request = Request(
            url,
            headers={
                "Accept": "application/json",
                "User-Agent": "DNA-Teaching-Simulator/1.0 (educational use)",
            },
            method="GET",
        )

        attempt = 0
        while True:
            self._throttle()
            self._last_request_started = time.monotonic()
            try:
                with urlopen(request, timeout=self.timeout) as response:
                    body = response.read()
                try:
                    parsed = json.loads(body.decode("utf-8"))
                except (UnicodeDecodeError, json.JSONDecodeError) as exc:
                    raise UniProtAPIError(200, f"JSON 解析失败：{exc}") from exc
                if not isinstance(parsed, dict):
                    raise UniProtAPIError(200, "响应 JSON 根节点不是对象。")
                return parsed

            except HTTPError as exc:
                body = exc.read() if hasattr(exc, "read") else b""
                detail = self._error_detail_from_body(body)
                retryable = exc.code == 429 or 500 <= exc.code < 600
                if retryable and attempt < self.max_retries:
                    retry_after = exc.headers.get("Retry-After") if exc.headers else None
                    try:
                        wait = float(retry_after) if retry_after is not None else 2.0 ** attempt
                    except (TypeError, ValueError):
                        wait = 2.0 ** attempt
                    time.sleep(max(wait, self.min_interval))
                    attempt += 1
                    continue
                raise UniProtAPIError(exc.code, detail) from exc

            except (URLError, TimeoutError, socket.timeout, ConnectionError) as exc:
                if attempt < self.max_retries:
                    time.sleep(max(2.0 ** attempt, self.min_interval))
                    attempt += 1
                    continue
                reason = getattr(exc, "reason", None) or exc
                raise UniProtNetworkError(
                    f"网络/超时错误（已尝试 {attempt + 1} 次）：{reason}"
                ) from exc

            except OSError as exc:
                if attempt < self.max_retries:
                    time.sleep(max(2.0 ** attempt, self.min_interval))
                    attempt += 1
                    continue
                raise UniProtNetworkError(f"网络 I/O 错误：{exc}") from exc

    @staticmethod
    def _protein_name(entry: dict) -> str:
        description = entry.get("proteinDescription") or {}
        recommended = description.get("recommendedName") or {}
        full_name = recommended.get("fullName") or {}
        if full_name.get("value"):
            return str(full_name["value"])
        submission_names = description.get("submissionNames") or []
        if submission_names:
            full_name = (submission_names[0] or {}).get("fullName") or {}
            if full_name.get("value"):
                return str(full_name["value"])
        return "（UniProt 未提供推荐蛋白名）"

    @staticmethod
    def _function_annotation(entry: dict) -> str:
        blocks: list[str] = []
        for comment in entry.get("comments") or []:
            if comment.get("commentType") != "FUNCTION":
                continue
            for text_obj in comment.get("texts") or []:
                value = text_obj.get("value")
                if value:
                    blocks.append(str(value))
        return " ".join(blocks).strip()

    @staticmethod
    def _organism(entry: dict) -> str:
        organism = entry.get("organism") or {}
        scientific = organism.get("scientificName") or ""
        common = organism.get("commonName") or ""
        if scientific and common:
            return f"{scientific} ({common})"
        return scientific or common or "（未知物种）"

    def _hit_from_entry(
        self,
        entry: dict,
        *,
        match_type: str = "exact",
        identity_percent: Optional[float] = 100.0,
        coverage_percent: Optional[float] = 100.0,
        combined_match_percent: Optional[float] = 100.0,
        evalue: Optional[float] = None,
        bitscore: Optional[float] = None,
        alignment_length: Optional[int] = None,
    ) -> UniProtHit:
        seq_block = entry.get("sequence") or {}
        entry_type = str(entry.get("entryType") or "")
        score_value = entry.get("annotationScore")
        try:
            score = float(score_value) if score_value is not None else None
        except (TypeError, ValueError):
            score = None
        return UniProtHit(
            accession=str(entry.get("primaryAccession") or "?"),
            entry_id=str(entry.get("uniProtkbId") or ""),
            protein_name=self._protein_name(entry),
            organism=self._organism(entry),
            function=self._function_annotation(entry),
            # 不能简单判断 "reviewed" in entryType：
            # "UniProtKB unreviewed (TrEMBL)" 也包含子串 reviewed。
            reviewed=entry_type.startswith("UniProtKB reviewed"),
            annotation_score=score,
            length=seq_block.get("length"),
            match_type=match_type,
            identity_percent=identity_percent,
            coverage_percent=coverage_percent,
            combined_match_percent=combined_match_percent,
            evalue=evalue,
            bitscore=bitscore,
            alignment_length=alignment_length,
        )

    def _extract_exact_hits(self, payload: dict, protein: str) -> list[UniProtHit]:
        """从 JSON 响应中解析条目，并仅保留完整序列完全一致的命中。"""
        hits: list[UniProtHit] = []
        for entry in payload.get("results") or []:
            seq_block = entry.get("sequence") or {}
            entry_seq = str(seq_block.get("value") or "").upper()
            if entry_seq != protein:
                continue
            hits.append(self._hit_from_entry(entry, match_type="exact"))
        # 优先 Swiss-Prot 已审阅条目，再按注释评分。
        hits.sort(
            key=lambda hit: (
                0 if hit.reviewed else 1,
                -(hit.annotation_score if hit.annotation_score is not None else -1),
            )
        )
        return hits

    # ------------------------------------------------------------------
    # 相似蛋白回退：GenomeNet BLASTP（Swiss-Prot）+ UniProtKB 注释
    # ------------------------------------------------------------------

    def _throttle_blast(self) -> None:
        now = time.monotonic()
        wait_seconds = UNIPROT_BLAST_MIN_INTERVAL_SECONDS - (
            now - self._last_blast_started
        )
        if wait_seconds > 0:
            time.sleep(wait_seconds)

    @staticmethod
    def _parse_genomenet_hits(page: str) -> list[dict]:
        """
        从 GenomeNet BLASTP 结果 HTML 中提取 hits。
        只依赖其公开结果表格中的 checkbox + accession + bit score + E-value。
        """
        if re.search(r"\bno\s+hits\b", page, re.I):
            return []
        title_ok = "BLASTP Search Result" in page or "BLAST Search Result" in page
        if not title_ok:
            raise SimilaritySearchError(
                "GenomeNet BLAST 未返回可识别的结果页，可能被限流或服务暂时不可用。"
            )

        accession_pattern = (
            r"([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})"
            r"(?:-\d+)?"
        )
        parts = re.split(r'<input[^>]*name="ckid"', page, flags=re.I)
        hits: list[dict] = []
        for part in parts[1:]:
            # 每条命中通常位于该 marker 后的一行；限制窗口避免串到下一条。
            window = part[:3000]
            accession_match = re.search(r"\[\s*" + accession_pattern + r"\s*\]", window)
            if not accession_match:
                continue
            accession = accession_match.group(1)
            numeric_anchors = re.findall(
                r"<a\b[^>]*>\s*([0-9]+(?:\.[0-9]+)?)\s*</a>",
                window,
                re.I,
            )
            bitscore = None
            if numeric_anchors:
                try:
                    bitscore = float(numeric_anchors[-1])
                except ValueError:
                    bitscore = None
            # E-value 位于 bits 链接之后，如 “103</a>   4e-29”。
            evalue = None
            evalue_match = re.search(
                r"</a>\s+([0-9]+(?:\.[0-9]+)?(?:[eE][+-]?\d+)?)",
                window,
            )
            if evalue_match:
                try:
                    evalue = float(evalue_match.group(1))
                except ValueError:
                    evalue = None
            hits.append(
                {
                    "accession": accession,
                    "bitscore": bitscore,
                    "evalue": evalue,
                }
            )
        return hits

    def _run_genomenet_blast(self, sequence: str) -> list[dict]:
        """
        提交一次同步 BLASTP 到 GenomeNet。限制在 Swiss-Prot，避免大型库和长耗时。

        这是公开 Web 服务，不要求邮箱；程序设置较长超时和最小请求间隔，且每次
        分析最多对少量候选发起。服务端若限流或改变页面格式，会以错误状态返回。
        """
        self._throttle_blast()
        fields = {
            "sequence": ">query\n" + sequence + "\n",
            "prog": "blastp",
            "dbname": UNIPROT_BLAST_DATABASE,
            "matrix": "blosum62",
            "filter": UNIPROT_BLAST_FILTER,
            "V_value": str(max(UNIPROT_SIMILARITY_FETCH_LIMIT, 20)),
            "B_value": "20",
            "Evalue": "10.0",
            "alignment_view": "0",
            "other_param": "",
        }
        data = urlencode(fields).encode("utf-8")
        request = Request(
            UNIPROT_BLAST_ENDPOINT,
            data=data,
            headers={
                "Content-Type": "application/x-www-form-urlencoded",
                "Accept": "text/html,application/xhtml+xml",
                "User-Agent": (
                    "DNA-Sequence-Analysis-Teaching-Tool/1.0 "
                    "(educational; Python urllib)"
                ),
            },
            method="POST",
        )
        self._last_blast_started = time.monotonic()
        try:
            with urlopen(request, timeout=UNIPROT_BLAST_TIMEOUT_SECONDS) as response:
                raw = response.read()
        except HTTPError as exc:
            detail = ""
            try:
                detail = exc.read().decode("utf-8", errors="replace")[:400]
            except Exception:
                pass
            raise SimilaritySearchError(
                f"GenomeNet BLAST HTTP {exc.code}: {detail or '无错误详情'}"
            ) from exc
        except (URLError, TimeoutError, socket.timeout, ConnectionError, OSError) as exc:
            raise SimilaritySearchError(
                f"GenomeNet BLAST 网络/超时错误：{exc}"
            ) from exc

        page = raw.decode("utf-8", errors="replace")
        return self._parse_genomenet_hits(page)

    def _fetch_entries_by_accessions(self, accessions: list[str]) -> dict[str, dict]:
        """一次检索最多若干 accession 的 UniProtKB JSON 条目。"""
        unique_accessions = list(dict.fromkeys(accessions))
        if not unique_accessions:
            return {}
        query = " OR ".join(
            f"accession:{accession}" for accession in unique_accessions
        )
        page_size = min(max(len(unique_accessions), 1), 25)
        payload = self._request_json(query, page_size=page_size)
        entries: dict[str, dict] = {}
        for entry in payload.get("results") or []:
            accession = str(entry.get("primaryAccession") or "")
            if accession:
                entries[accession] = entry
        return entries

    def _get_pairwise_aligner(self) -> PairwiseAligner:
        if self._aligner is None:
            aligner = PairwiseAligner()
            aligner.mode = "local"
            aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
            aligner.open_gap_score = -11.0
            aligner.extend_gap_score = -1.0
            self._aligner = aligner
        return self._aligner

    @staticmethod
    def _clean_protein(sequence: str) -> str:
        cleaned = (sequence or "").strip().upper().replace("*", "")
        allowed = set("ACDEFGHIKLMNPQRSTVWYBZX")
        return "".join(char if char in allowed else "X" for char in cleaned)

    def _calculate_similarity(
        self,
        query: str,
        subject: str,
    ) -> Optional[dict]:
        """
        对查询和候选序列做局部 Smith-Waterman 比对，计算：
        - identity_percent：一致残基数 / 比对列数（包括 gap 列）；
        - coverage_percent：比对上的查询残基数 / 查询总长度；
        - combined_match_percent：identity × coverage / 100。
        """
        clean_query = self._clean_protein(query)
        clean_subject = self._clean_protein(subject)
        if not clean_query or not clean_subject:
            return None
        try:
            alignments = self._get_pairwise_aligner().align(clean_query, clean_subject)
            if not alignments:
                return None
            alignment = alignments[0]
            aligned_query = str(alignment[0])
            aligned_subject = str(alignment[1])
        except (ValueError, IndexError, MemoryError):
            return None

        columns = min(len(aligned_query), len(aligned_subject))
        if columns <= 0:
            return None
        aligned_query = aligned_query[:columns]
        aligned_subject = aligned_subject[:columns]
        identities = sum(
            1
            for query_char, subject_char in zip(aligned_query, aligned_subject)
            if query_char != "-" and query_char == subject_char
        )
        query_aligned = sum(1 for char in aligned_query if char != "-")
        if query_aligned == 0:
            return None
        identity = identities / columns * 100.0
        coverage = query_aligned / len(clean_query) * 100.0
        combined = identity * coverage / 100.0
        return {
            "identity_percent": round(identity, 1),
            "coverage_percent": round(coverage, 1),
            "combined_match_percent": round(combined, 1),
            "alignment_length": columns,
        }

    def _lookup_similar_proteins(
        self,
        sequence: str,
        preceding_note: str,
    ) -> UniProtResult:
        """
        精确匹配失败后的相似性回退。
        返回的 hits 已按综合匹配度从高到低排序。
        """
        query_description = (
            f"GenomeNet BLASTP / {UNIPROT_BLAST_DATABASE}; "
            f"query length {len(sequence)} aa"
        )
        prefix = "UniProt 全长精确检索未命中"
        if preceding_note:
            prefix += f"（{preceding_note}）"
        if not self.enable_similarity:
            return self._result_no_match(
                method="md5_checksum",
                query="checksum 精确检索；相似性搜索已关闭",
                note=f"{prefix}；相似性搜索已关闭。",
            )

        try:
            blast_hits = self._run_genomenet_blast(sequence)
        except SimilaritySearchError as exc:
            return UniProtResult(
                status="error",
                method="genomenet_blast",
                query=query_description,
                message=f"{prefix}；相似性搜索失败：{exc}",
            )

        if not blast_hits:
            return self._result_no_match(
                method="genomenet_blast",
                query=query_description,
                note=(
                    f"{prefix}；GenomeNet BLASTP"
                    f"（{UNIPROT_BLAST_DATABASE}）也没有找到显著相似蛋白。"
                ),
            )

        # GenomeNet 的 bit score 排序取得候选 accession；限制一次注释数量。
        accessions = [
            hit["accession"]
            for hit in blast_hits
            if hit.get("accession")
        ][:UNIPROT_SIMILARITY_FETCH_LIMIT]
        blast_metrics = {
            hit["accession"]: hit for hit in blast_hits if hit.get("accession")
        }
        try:
            entries = self._fetch_entries_by_accessions(accessions)
        except UniProtAPIError as exc:
            return UniProtResult(
                status="error",
                method="genomenet_blast",
                query=query_description,
                message=(
                    f"{prefix}；BLAST 找到候选，但获取 UniProtKB 注释失败：{exc}"
                ),
            )
        except UniProtNetworkError as exc:
            return UniProtResult(
                status="error",
                method="genomenet_blast",
                query=query_description,
                message=(
                    f"{prefix}；BLAST 找到候选，但获取 UniProtKB 注释时网络异常：{exc}"
                ),
            )

        similar_hits: list[UniProtHit] = []
        for accession in accessions:
            entry = entries.get(accession)
            if not entry:
                continue
            entry_seq = str((entry.get("sequence") or {}).get("value") or "")
            metrics = self._calculate_similarity(sequence, entry_seq)
            if not metrics:
                continue
            identity = float(metrics["identity_percent"])
            coverage = float(metrics["coverage_percent"])
            if (
                identity < UNIPROT_SIMILARITY_MIN_IDENTITY
                or coverage < UNIPROT_SIMILARITY_MIN_COVERAGE
            ):
                continue
            blast_info = blast_metrics.get(accession, {})
            similar_hits.append(
                self._hit_from_entry(
                    entry,
                    match_type="similar",
                    identity_percent=identity,
                    coverage_percent=coverage,
                    combined_match_percent=float(metrics["combined_match_percent"]),
                    evalue=blast_info.get("evalue"),
                    bitscore=blast_info.get("bitscore"),
                    alignment_length=int(metrics["alignment_length"]),
                )
            )

        # 按“综合匹配度”排序，次排序一致度、覆盖度，再优先 reviewed。
        similar_hits.sort(
            key=lambda hit: (
                -(hit.combined_match_percent or 0.0),
                -(hit.identity_percent or 0.0),
                -(hit.coverage_percent or 0.0),
                0 if hit.reviewed else 1,
                hit.accession,
            )
        )
        similar_hits = similar_hits[:UNIPROT_SIMILARITY_MAX_HITS]

        if not similar_hits:
            return self._result_no_match(
                method="genomenet_blast",
                query=query_description,
                note=(
                    f"{prefix}；GenomeNet BLASTP 找到的候选未达到本地一致性过滤阈值"
                    f"（一致度 ≥ {UNIPROT_SIMILARITY_MIN_IDENTITY:.0f}%，"
                    f"覆盖度 ≥ {UNIPROT_SIMILARITY_MIN_COVERAGE:.0f}%）。"
                ),
            )

        notes = [
            f"{prefix}；已通过 GenomeNet BLASTP"
            f"（数据库：{UNIPROT_BLAST_DATABASE}）检索相似蛋白，并从 UniProtKB"
            " 获取名称、物种与 FUNCTION 注释。",
            "以下结果是基于序列相似性的推断，不是全长精确匹配；"
            "按“一致度 × 查询覆盖度”的综合匹配度从高到低排列。",
            "注意：相似性命中不代表该蛋白在真实细胞中表达，也不能替代实验验证。",
        ]
        return self._result_matched(
            similar_hits,
            method="genomenet_blast",
            query=query_description,
            note=" ".join(note.strip() for note in notes),
        )

    def _result_matched(self, hits: list[UniProtHit], method: str, query: str, note: str) -> UniProtResult:
        return UniProtResult(
            status="matched",
            method=method,
            query=query,
            hits=hits[:UNIPROT_MAX_HITS_TO_PRINT],
            message=note,
        )

    def _result_no_match(self, method: str, query: str, note: str) -> UniProtResult:
        return UniProtResult(
            status="no_match",
            method=method,
            query=query,
            hits=[],
            message=note,
        )

    def lookup_protein(self, protein: str) -> UniProtResult:
        """查询一个蛋白序列的 UniProtKB 注释；不抛出预期网络异常。"""
        sequence = (protein or "").strip().upper().rstrip("*")
        sequence = "".join(character for character in sequence if not character.isspace())
        if not sequence:
            return UniProtResult(
                status="error",
                message="蛋白序列为空，未发起 UniProt 查询。",
            )
        if sequence in self._cache:
            return self._cache[sequence]

        fallback_note = ""
        primary_error: Optional[UniProtError] = None

        if self.try_sequence_field_first and self._sequence_field_supported is not False:
            if len(sequence) <= UNIPROT_SEQUENCE_FIELD_MAX_AA:
                primary_query = f'sequence:"{sequence}"'
                try:
                    payload = self._request_json(primary_query)
                    self._sequence_field_supported = True
                    primary_hits = self._extract_exact_hits(payload, sequence)
                    if primary_hits:
                        result = self._result_matched(
                            primary_hits,
                            method="sequence_field",
                            query=primary_query,
                            note="已通过 UniProt sequence 字段检索到全长精确匹配。",
                        )
                        self._cache[sequence] = result
                        return result
                    fallback_note = (
                        "sequence 字段查询已成功，但未返回全长完全一致的条目；"
                        "继续用 MD5 校验和做精确检索。"
                    )
                except UniProtAPIError as exc:
                    primary_error = exc
                    if exc.status == 400:
                        # 当前公开 API 会拒绝 sequence 字段；缓存该能力结果，
                        # 避免后续每个候选都重复发送同一个无效请求。
                        self._sequence_field_supported = False
                        fallback_note = (
                            "UniProt 当前 REST API 拒绝 sequence 检索字段"
                            f"（HTTP 400：{exc.detail}）。已自动改用 MD5 checksum "
                            "精确检索，仍使用 /uniprotkb/search 端点。"
                        )
                    else:
                        result = UniProtResult(
                            status="error",
                            method="sequence_field",
                            query=primary_query,
                            message=f"UniProt sequence 字段查询失败：{exc}",
                        )
                        self._cache[sequence] = result
                        return result
                except UniProtNetworkError as exc:
                    result = UniProtResult(
                        status="error",
                        method="sequence_field",
                        query=primary_query,
                        message=f"UniProt 查询网络异常：{exc}",
                    )
                    self._cache[sequence] = result
                    return result
            else:
                fallback_note = (
                    "蛋白长度超过 "
                    f"{UNIPROT_SEQUENCE_FIELD_MAX_AA} aa，为避免 GET URL 过长，"
                    "直接使用 MD5 校验和精确检索。"
                )

        # ---------- 回退方案：MD5 checksum 精确检索 ----------
        md5_hex = hashlib.md5(sequence.encode("ascii", errors="strict")).hexdigest().upper()
        checksum_query = f"checksum:{md5_hex} AND length:[{len(sequence)} TO {len(sequence)}]"
        try:
            payload = self._request_json(checksum_query)
        except UniProtAPIError as exc:
            detail = f"UniProt checksum 查询失败：{exc}"
            if primary_error is not None:
                detail += f"；最初的 sequence 字段错误：{primary_error}"
            result = UniProtResult(
                status="error",
                method="md5_checksum",
                query=checksum_query,
                message=detail,
            )
            self._cache[sequence] = result
            return result
        except UniProtNetworkError as exc:
            result = UniProtResult(
                status="error",
                method="md5_checksum",
                query=checksum_query,
                message=f"UniProt 查询网络异常：{exc}",
            )
            self._cache[sequence] = result
            return result

        hits = self._extract_exact_hits(payload, sequence)
        if hits:
            notes = [fallback_note] if fallback_note else ["MD5 校验和精确匹配成功。"]
            notes.append(
                "注意：命中仅表示该蛋白序列存在于 UniProtKB；数据库注释不能证明"
                "这段 DNA 在真实细胞中表达出该蛋白。"
            )
            result = self._result_matched(
                hits,
                method="md5_checksum",
                query=checksum_query,
                note=" ".join(notes),
            )
            self._cache[sequence] = result
            return result

        result = self._lookup_similar_proteins(
            sequence,
            preceding_note=fallback_note,
        )
        self._cache[sequence] = result
        return result


# ====================================================================
# 十、命令行与结果打印
# ====================================================================

def print_original_sequence_info(
    dna: str,
    source: str,
    warnings: list[str],
) -> None:
    """打印原始序列组成、正链和反向互补链。"""
    print_section("[1] 原始 DNA 输入与基础序列处理")
    print(f"  来源: {source}")
    print(f"  长度: {len(dna)} bp")
    counts = {base: dna.count(base) for base in ("A", "C", "G", "T", "N")}
    canonical = counts["A"] + counts["C"] + counts["G"] + counts["T"]
    gc = counts["G"] + counts["C"]
    gc_percent = (gc / canonical * 100.0) if canonical else 0.0
    print(
        "  碱基组成: "
        + "  ".join(f"{base}={count}" for base, count in counts.items())
        + f"  GC={gc_percent:.2f}%"
    )
    for warning in warnings:
        print(f"  [提示] {warning}")
    print(
        "  说明: 自动生成反向互补链，并分别在正链与负链的 3 个读码框中扫描。"
    )
    print(
        "  坐标约定: 后续所有“输入正向坐标”均以本次输入序列的第 1 个碱基为 1，"
        "不是 hg38/GRCh38 等参考基因组坐标。"
    )

    print_subsection("正向 DNA（5'→3'）")
    print_sequence_block("正向链", dna)
    print_subsection("反向互补链（5'→3'）")
    print_sequence_block("反向互补链", reverse_complement_dna(dna))


def print_candidate_summary(
    search_result: GeneSearchResult,
    min_aa: int,
    max_candidates: int,
) -> None:
    """打印六读码框统计和候选表。"""
    print_section("[2] 六读码框 ORF 检索与候选筛选")
    print(f"  完整 ORF（ATG 起始且有同框终止密码子）总数: {len(search_result.all_complete_orfs)}")
    frame_text = "  ".join(
        f"{frame}:{count}" for frame, count in search_result.frame_counts.items()
    )
    print(f"  各读码框完整 ORF 数: {frame_text}")
    print(
        "  说明: 同一读码框的同一终止区间只输出最长的 ATG…终止 ORF；"
        "无终止子的不完整 ORF 不输出。"
    )
    print(f"  最小氨基酸长度阈值: ≥{min_aa} aa")
    print(
        f"  通过最小长度阈值的候选: {len(search_result.length_passing_candidates)} 个；"
        f"其中检出 ATG 上游 TATAAA 的: {search_result.tata_pass_count} 个"
    )
    if search_result.tata_filter_enabled:
        filtered_out = (
            len(search_result.length_passing_candidates) - search_result.tata_pass_count
        )
        print(f"  TATA 过滤开关: 开启（过滤掉 {filtered_out} 个无上游 TATAAA 的候选）")
    else:
        print("  TATA 过滤开关: 关闭（TATAAA 仅作为辅助标记，不参与 ORF 定义）")
    print(
        "  注意: TATAAA/AATAAA 是字面基序标记，不是严谨的启动子或 poly(A) 预测；"
        "ORF 本身不要求 TATA 盒。"
    )

    if not search_result.length_passing_candidates:
        print("\n  [结果] 没有任何 ORF 通过最小氨基酸长度阈值。")
        print("         可尝试降低 MIN_ORF_AA 或使用 --min-aa；当前输入也可能没有完整 ORF。")
        return

    shown = search_result.length_passing_candidates[:MAX_ORF_TABLE_ROWS]

    print()
    print(
        "  候选表（坐标相对于本次输入序列，1-based 闭区间，含终止密码子；"
        "“通过筛选”指满足当前 TATA 开关）:"
    )
    header = (
        f"  {'序号':<5} {'ORF':<9} {'链':<5} {'框':<4} "
        f"{'输入正向坐标':<18} {'aa':>4}  {'TATA':<6} {'AATAAA':<6} {'通过筛选'}"
    )
    print(header)
    print("  " + "-" * (len(header) - 2))
    for index, candidate in enumerate(shown, start=1):
        orf = candidate.orf
        coords = f"{orf.genomic_start + 1}-{orf.genomic_end}"
        tata_text = f"是({len(candidate.tata_hits)})" if candidate.has_tata else "否"
        polya_text = f"是({len(candidate.polya_hits)})" if candidate.has_polya_signal else "否"
        passes_filter = candidate.has_tata if search_result.tata_filter_enabled else True
        print(
            f"  {index:<5} {orf.orf_id:<9} {orf.strand:<5} {orf.frame:<4} "
            f"{coords:<18} {orf.aa_length:>4}  {tata_text:<6} {polya_text:<6} "
            f"{'是' if passes_filter else '否'}"
        )
    if len(search_result.length_passing_candidates) > len(shown):
        print(
            f"  ... 共 {len(search_result.length_passing_candidates)} 个，"
            f"表格仅显示前 {len(shown)} 个。"
        )

    selected_count = len(search_result.selected_candidates)
    if selected_count:
        print(
            f"\n  本次将对 {selected_count} 个候选执行 pre-mRNA 剪接、翻译和 UniProt 查询"
            f"（上限 --max-candidates={max_candidates}）。"
        )
        print(
            "  候选编号: "
            + ", ".join(candidate.candidate_id for candidate in search_result.selected_candidates)
        )
        print(
            "  编号对应: "
            + ", ".join(
                f"{candidate.candidate_id}={candidate.orf.orf_id}"
                for candidate in search_result.selected_candidates
            )
        )
        if search_result.tata_filter_enabled and search_result.tata_pass_count > selected_count:
            print(
                f"  提示: 有 {search_result.tata_pass_count} 个候选通过 TATA 筛选，"
                f"但只详细分析排序最靠前的 {selected_count} 个。"
            )
        elif not search_result.tata_filter_enabled and len(search_result.length_passing_candidates) > selected_count:
            print(
                f"  提示: 有 {len(search_result.length_passing_candidates)} 个候选通过长度阈值，"
                f"但只详细分析排序最靠前的 {selected_count} 个。"
            )
    else:
        if search_result.tata_filter_enabled:
            print(
                "\n  [结果] TATA 过滤后没有候选。可关闭 REQUIRE_TATA_BOX 或使用"
                " --no-require-tata 重新运行。"
            )
        else:
            print("\n  [结果] 没有可进入下游流程的候选。")


def print_candidate_header(candidate: CandidateGene) -> None:
    orf = candidate.orf
    print_section(
        f"[3] 候选 {candidate.candidate_id} · {strand_label(orf.strand)} · "
        f"读码框 {orf.strand}{orf.frame} · {orf.aa_length} aa"
    )
    print(f"  ORF 编号: {orf.orf_id}")
    print(
        f"  ORF 输入正向坐标: {interval_1based(orf.genomic_start, orf.genomic_end)}"
        f"（{orf.nt_length} bp，含终止密码子 {orf.stop_codon}）"
    )
    print(
        f"  候选区域输入正向坐标: {interval_1based(candidate.region_start, candidate.region_end)}"
        f"（{candidate.region_end - candidate.region_start} bp，含上下游侧翼）"
    )
    print(
        f"  ORF 在 pre-mRNA 中的坐标: "
        f"{interval_1based(candidate.orf_pre_start, candidate.orf_pre_end)}"
    )
    print(f"  原始起始/终止密码子: {orf.start_codon} / {orf.stop_codon}")
    print(f"  ORF 原始翻译产物长度: {orf.aa_length} aa")


def print_regulatory_markers(candidate: CandidateGene) -> None:
    print_subsection("[3.1] 简易调控信号（TATAAA / AATAAA，仅作教学标记）")
    print(
        f"  TATA 盒字面基序 {TATA_BOX_MOTIF}（pre-mRNA 中显示为 UAUAAA）: "
        f"{len(candidate.tata_hits)} 个上游命中"
    )
    if candidate.tata_hits:
        for index, hit in enumerate(candidate.tata_hits[:MAX_REGULATORY_SITES_TO_PRINT], start=1):
            print(
                f"    T{index}: pre-mRNA {interval_1based(hit.start, hit.end)} | "
                f"输入正向 {interval_1based(hit.genomic_start, hit.genomic_end)} | "
                f"距 ATG 上游 {hit.distance_bp} bp"
            )
        if len(candidate.tata_hits) > MAX_REGULATORY_SITES_TO_PRINT:
            print(f"    ... 另有 {len(candidate.tata_hits) - MAX_REGULATORY_SITES_TO_PRINT} 个命中未列出。")
    else:
        print("    未在 ORF 上游指定窗口内发现 TATAAA。")
    print(
        f"  poly(A) 信号字面基序 {POLY_A_SIGNAL_MOTIF}（pre-mRNA 中显示为 AAUAAA）: "
        f"{len(candidate.polya_hits)} 个下游命中"
    )
    if candidate.polya_hits:
        for index, hit in enumerate(candidate.polya_hits[:MAX_REGULATORY_SITES_TO_PRINT], start=1):
            print(
                f"    P{index}: pre-mRNA {interval_1based(hit.start, hit.end)} | "
                f"输入正向 {interval_1based(hit.genomic_start, hit.genomic_end)} | "
                f"距终止密码子下游 {hit.distance_bp} bp"
            )
        if len(candidate.polya_hits) > MAX_REGULATORY_SITES_TO_PRINT:
            print(f"    ... 另有 {len(candidate.polya_hits) - MAX_REGULATORY_SITES_TO_PRINT} 个命中未列出。")
    else:
        print("    未在终止密码子下游指定窗口内发现 AATAAA。")
    print(
        "  注意: 单个六碱基字面匹配不能证明启动子活性或 poly(A) 加尾；"
        "真实调控元件依赖序列上下文、距离、染色质和反式因子。"
    )


def print_pre_mrna(candidate: CandidateGene) -> None:
    print_subsection("[3.2] 候选区域转录为 pre-mRNA")
    print(
        f"  处理: 按基因方向取候选区域 DNA，使用 Biopython 将 T 转录为 U；"
        f"pre-mRNA 长度 {len(candidate.pre_mrna)} nt。"
    )
    print(
        "  教学简化: 候选区域包含启动子上游侧翼，本程序将整段序列一起转录；"
        "真实 pre-mRNA 从转录起始位点开始，启动子 DNA 不会出现在 RNA 中。"
    )
    print_sequence_block("pre-mRNA（5'→3'）", candidate.pre_mrna)


def print_repeat_screen(
    repeats: list[TandemRepeatRegion],
    removed_intervals: list[tuple[int, int]],
    removal_enabled: bool,
    splice_input_length: int,
) -> None:
    print_subsection("[3.3] 串联重复 / 疑似转座子粗筛")
    if not SCREEN_TANDEM_REPEATS:
        print("  筛查开关 SCREEN_TANDEM_REPEATS=False，已跳过。")
        return
    print(
        "  说明: 这里只做短串联重复（低复杂度重复）粗筛，不等于转座子鉴定；"
        "不构建重复库，不检测末端重复或转座酶。"
    )
    if not repeats:
        print("  未发现满足当前参数的串联重复区域。")
        return
    print(f"  发现 {len(repeats)} 个疑似串联重复区域（pre-mRNA 1-based 坐标）:")
    for index, repeat in enumerate(repeats[:MAX_REPEAT_REGIONS_TO_PRINT], start=1):
        print(
            f"    R{index}: pre-mRNA {interval_1based(repeat.start, repeat.end)}"
            f"（{repeat.length} nt），代表性单元 {repeat.unit} × 约 {repeat.copies} 次"
        )
    if len(repeats) > MAX_REPEAT_REGIONS_TO_PRINT:
        print(f"    ... 另有 {len(repeats) - MAX_REPEAT_REGIONS_TO_PRINT} 个区域未列出。")

    if removal_enabled:
        print(
            f"  已按 REMOVE_REPEAT_REGIONS=True 剔除 {len(removed_intervals)} 个合并区域；"
            f"用于剪接的 pre-mRNA 长度由原始值变为 {splice_input_length} nt。"
        )
        print(
            "  注意: 强行删除重复区域可能破坏外显子/剪接位点，属于非生物学的人工操作。"
        )
    else:
        print("  当前配置仅标记、不剔除（REMOVE_REPEAT_REGIONS=False）。")


def print_splicing_result(
    splicing: SplicingResult,
    repeat_removal_note: str,
) -> None:
    print_subsection("[3.4] pre-mRNA 剪接（简化的 GU-AG / DNA 记法 GT-AG）")
    print(
        "  说明: 在 RNA 层面识别 GU…AG；不是 DNA 层操作，也不使用真实剪接位点模型。"
    )
    print(f"  用于剪接的序列: {repeat_removal_note}，长度 {len(splicing.input_mrna)} nt。")
    print(f"  识别到并剔除的内含子数: {len(splicing.introns)}")

    for index, intron in enumerate(splicing.introns, start=1):
        donor_context = splicing.input_mrna[max(0, intron.start - 10):intron.start]
        acceptor_following = splicing.input_mrna[
            intron.end:min(len(splicing.input_mrna), intron.end + 10)
        ]
        print(
            f"    内含子 {index}: pre-mRNA {interval_1based(intron.start, intron.end)}"
            f"（{intron.length} nt） | 供体 {intron.donor} | 受体 {intron.acceptor}"
        )
        print(
            f"      供体位点: ...{donor_context}[{intron.donor}]...；"
            f"受体位点: ...[{intron.acceptor}]{acceptor_following}..."
        )

    if splicing.exons:
        exon_parts = [
            f"E{index} {interval_1based(start, end)}（{end - start} nt）"
            for index, (start, end) in enumerate(splicing.exons, start=1)
        ]
        print("  拼接后的外显子: " + " | ".join(exon_parts))
        if len(splicing.exons) > 1:
            junction_parts = []
            for left, right in zip(splicing.exons, splicing.exons[1:]):
                left_context = splicing.input_mrna[max(0, left[1] - 8):left[1]]
                right_context = splicing.input_mrna[right[0]:right[0] + 8]
                junction_parts.append(f"{left_context}|{right_context}")
            print("  外显子连接点（各取连接点两侧 8 nt）: " + "   ".join(junction_parts))

    for note in splicing.notes:
        print(f"  [提示] {note}")
    print(
        "  重要: 任意符合 GU…AG + 长度约束的位点都可能被误判；不含分支点、"
        "多嘧啶区、剪接增强子或可变剪接分析。"
    )
    print_sequence_block("成熟 mRNA（5'→3'，剪接后）", splicing.mature_mrna)


def print_translation(
    translation: TranslationResult,
    orf_protein: str,
) -> None:
    print_subsection("[3.5] 成熟 mRNA 翻译（标准遗传密码 table 1）")
    if translation.start_index is None:
        print("  [结果] 成熟 mRNA 中无 AUG，无蛋白产物。")
        for warning in translation.warnings:
            print(f"  [警告] {warning}")
        return

    print(f"  第一个 AUG 起点（成熟 mRNA 坐标）: {translation.start_index + 1}（0-based: {translation.start_index}）")
    print(f"  使用的 CDS 长度: {len(translation.cds_rna)} nt（含终止密码子）")
    if translation.stop_codon is not None and translation.stop_codon_start is not None:
        print(
            f"  同框终止密码子: {translation.stop_codon}，"
            f"位置 {translation.stop_codon_start + 1}（成熟 mRNA 1-based 起始）"
        )
    else:
        print("  同框终止密码子: 未找到（翻译到末端或最后一个完整密码子）")
    print(f"  蛋白长度: {len(translation.protein)} aa")
    print_sequence_block(
        "蛋白质氨基酸序列",
        translation.protein,
        max_chars=len(translation.protein) + 1,
        unit="aa",
    )

    if translation.protein and translation.protein != orf_protein:
        print(
            "  [警告] 成熟 mRNA 的翻译结果与原始基因组 ORF 的翻译结果不同；"
            "可能是剪接/重复剔除改变了阅读框，也可能是原始 ORF 本身就跨越内含子。"
        )
    for warning in translation.warnings:
        print(f"  [警告] {warning}")


def print_uniprot_result(result: UniProtResult) -> None:
    print_subsection("[3.6] 蛋白质功能注释查询")
    print(f"  精确检索接口: {UNIPROT_BASE_URL}{UNIPROT_SEARCH_PATH}")
    if result.method == "genomenet_blast":
        print(
            "  相似性检索: GenomeNet BLASTP / Swiss-Prot "
            f"({UNIPROT_BLAST_ENDPOINT}) + UniProtKB 注释"
        )
    print(
        f"  限流/超时: UniProt 请求间隔 ≥ {UNIPROT_MIN_INTERVAL_SECONDS:.1f} s；"
        f"超时 {UNIPROT_TIMEOUT_SECONDS:.0f} s；最多重试 {UNIPROT_MAX_RETRIES} 次"
    )
    if result.method:
        method_text = {
            "sequence_field": "UniProt sequence 字段检索",
            "md5_checksum": "MD5 checksum 精确检索（/uniprotkb/search）",
            "genomenet_blast": "GenomeNet BLASTP 相似性检索 + UniProtKB 注释",
        }.get(result.method, result.method)
        print(f"  实际检索方式: {method_text}")
    if result.query:
        print(f"  查询式: {short_query_for_display(result.query)}")
    if result.message:
        print(wrap_text(result.message, width=96, initial_indent="  ", subsequent_indent="  "))

    if result.status == "error":
        print("  [结果] 查询失败；主流程其他部分仍已完成。")
        return
    if result.status == "no_match":
        print("  [结果] 无匹配条目。")
        return
    if result.status != "matched":
        print(f"  [结果] 未预期状态: {result.status}")
        return

    is_similarity = result.method == "genomenet_blast"
    if is_similarity:
        top_hit = result.hits[0]
        top_evalue = f"{top_hit.evalue:.2g}" if top_hit.evalue is not None else "—"
        print("  查询状态: 找到相似蛋白（非精确全长命中）")
        print(
            "  排名方式: 综合匹配度（一致度 × 查询覆盖度）从高到低；"
            f"命中 {len(result.hits)} 条"
        )
        print(
            f"  最相似条目（推测参考）: {top_hit.accession} · {top_hit.protein_name}"
            f" | {top_hit.organism} | 一致度 {top_hit.identity_percent:.1f}%"
            f" | 覆盖度 {top_hit.coverage_percent:.1f}%"
            f" | 综合匹配度 {top_hit.combined_match_percent:.1f}%"
            f" | E-value {top_evalue}"
        )
    else:
        print("  查询状态: 成功（全长序列精确匹配）")
        print(f"  命中 {len(result.hits)} 条（最多显示 {UNIPROT_MAX_HITS_TO_PRINT} 条）:")
    for index, hit in enumerate(result.hits, start=1):
        reviewed_text = "reviewed / Swiss-Prot" if hit.reviewed else "unreviewed / TrEMBL"
        print()
        print(f"    {index}. {hit.accession} | {hit.protein_name}")
        print(
            f"       物种: {hit.organism} | {reviewed_text}"
            + (f" | 长度: {hit.length} aa" if hit.length is not None else "")
        )
        if is_similarity:
            evalue_text = f"{hit.evalue:.2g}" if hit.evalue is not None else "—"
            bitscore_text = f"{hit.bitscore:.1f}" if hit.bitscore is not None else "—"
            print(
                f"       一致度: {hit.identity_percent:.1f}% | "
                f"查询覆盖度: {hit.coverage_percent:.1f}% | "
                f"综合匹配度: {hit.combined_match_percent:.1f}%"
            )
            print(
                f"       相似性证据: E-value {evalue_text} | bit score {bitscore_text}"
                + (f" | 比对长度 {hit.alignment_length} aa"
                   if hit.alignment_length is not None else "")
            )
        print(f"       UniProt: https://www.uniprot.org/uniprotkb/{hit.accession}/entry")
        function = hit.function.strip() if hit.function else "（数据库未提供 FUNCTION 注释）"
        if len(function) > UNIPROT_FUNCTION_MAX_CHARS:
            function = function[:UNIPROT_FUNCTION_MAX_CHARS] + "…"
        print(wrap_text(function, width=100, initial_indent="       功能: ", subsequent_indent="             "))
    if is_similarity:
        print(
            "\n  [提示] 以上为 GenomeNet BLASTP 相似性推断，并非全长精确命中。"
            "相似蛋白的 FUNCTION 注释只能作为功能线索，不能据此断言本序列的功能。"
        )
    else:
        print(
            "\n  [提示] UniProt 命中只说明数据库中存在相同序列条目；不能证明本输入 DNA "
            "在细胞内真实表达该蛋白，也不能替代实验验证。"
        )


def run_pipeline(
    dna: str,
    source: str,
    input_warnings: list[str],
    options: RuntimeOptions,
    skip_uniprot: bool,
    uniprot_client: Optional[UniProtClient],
) -> None:
    """执行并打印完整教学流程。"""
    print("=" * 78)
    print("DNA 序列分析教学模拟器")
    print("ORF → pre-mRNA → GU-AG 简化剪接 → 翻译 → UniProt / 相似蛋白注释")
    print("仅用于教学演示，不可用于科研、临床或任何正式决策")
    print("=" * 78)
    print(
        f"运行参数: 最小 ORF={options.min_orf_aa} aa | "
        f"TATA 过滤={'开' if options.require_tata_box else '关'} | "
        f"剔除重复={'开' if options.remove_repeat_regions else '关'} | "
        f"候选上限={options.max_candidates} | "
        f"UniProt={'关闭' if skip_uniprot else '开启'} | "
        f"相似性 BLAST="
        f"{'关闭' if skip_uniprot or not (uniprot_client and uniprot_client.enable_similarity) else '开启'}"
    )

    print_original_sequence_info(dna, source, input_warnings)

    search_result = search_gene_candidates(dna, options)
    print_candidate_summary(search_result, options.min_orf_aa, options.max_candidates)

    if not search_result.selected_candidates:
        print_section("流程结束")
        print("  没有候选进入下游分析。请检查最小 ORF 长度、TATA 过滤开关或输入序列。")
        print("  本程序只为教学演示，不构成任何生物预测结论。")
        return

    for candidate in search_result.selected_candidates:
        print_candidate_header(candidate)
        print_regulatory_markers(candidate)
        print_pre_mrna(candidate)

        repeats: list[TandemRepeatRegion] = []
        if SCREEN_TANDEM_REPEATS:
            repeats = find_tandem_repeats(candidate.pre_mrna)
        elif options.remove_repeat_regions:
            print(
                "  [警告] REMOVE_REPEAT_REGIONS 开启，但 SCREEN_TANDEM_REPEATS 关闭；"
                "没有可剔除的标记区域。"
            )

        removed_intervals: list[tuple[int, int]] = []
        if options.remove_repeat_regions and repeats:
            splice_input, removed_intervals = remove_marked_repeats(candidate.pre_mrna, repeats)
            repeat_removal_note = "已剔除重复标记区域的 pre-mRNA"
        else:
            splice_input = candidate.pre_mrna
            repeat_removal_note = "未剔除重复区域的 pre-mRNA"

        print_repeat_screen(
            repeats=repeats,
            removed_intervals=removed_intervals,
            removal_enabled=options.remove_repeat_regions,
            splice_input_length=len(splice_input),
        )

        splicing = splice_pre_mrna(splice_input)
        print_splicing_result(splicing, repeat_removal_note)

        translation = translate_mature_mrna(splicing.mature_mrna)
        print_translation(translation, candidate.orf.protein)

        # 没有蛋白就没有必要发起 API 请求。
        if translation.start_index is None or not translation.protein:
            print_subsection("[3.6] UniProtKB 功能注释查询")
            print("  [结果] 无蛋白产物，跳过 UniProt 查询。")
            continue

        if skip_uniprot:
            print_subsection("[3.6] UniProtKB 功能注释查询")
            print("  已通过 --skip-uniprot 关闭此步骤。")
            continue

        if uniprot_client is None:
            print_subsection("[3.6] UniProtKB 功能注释查询")
            print("  [结果] UniProt 客户端未初始化，跳过查询。")
            continue

        uniprot_result = uniprot_client.lookup_protein(translation.protein)
        print_uniprot_result(uniprot_result)

    print_section("流程结束与生物局限性声明")
    print(
        "  本程序把 ORF 扫描、六碱基调控基序、简单 GU-AG 剪接和标准翻译串联起来，"
        "目的是展示数据流和代码模块，不是预测真实基因结构。"
    )
    print(
        "  真实基因预测需要 AUGUSTUS/BRAKER/GENSCAN/GeneMark 等模型、RNA-seq/"
        "蛋白同源证据、剪接位点模型和人工审查；功能注释需要 BLAST/InterProScan/"
        "HMMER/结构域数据库等交叉验证。"
    )
    print(
        "  相似蛋白 BLAST 结果只是序列相似性推断，可能受数据库、阈值和外部 Web "
        "服务状态影响，不等于功能实验。"
    )
    print("  任何结果都不能用于正式科研、临床诊断、药物设计或生产决策。")
    for warning in input_warnings:
        print(f"  [输入提示] {warning}")


def parse_args(argv: Optional[list[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "DNA→六读码框ORF→pre-mRNA→GU-AG简化剪接→成熟mRNA→蛋白→UniProt注释"
            "（教学演示程序，结果不可用于正式科研）"
        ),
        epilog=(
            "无参数运行会打开 Tkinter 图形界面。命令行/脚本化分析请显式提供 "
            "--demo、--dna 或 --fasta。"
        ),
    )
    input_group = parser.add_mutually_exclusive_group()
    input_group.add_argument("--dna", type=str, help="原始 DNA 字符串；不含 FASTA 头")
    input_group.add_argument("--fasta", type=str, help="FASTA 文件路径")
    input_group.add_argument(
        "--demo",
        action="store_true",
        help="使用内置合成教学 DNA（命令行模式；无参数运行会打开 GUI）",
    )
    parser.add_argument("--record-id", type=str, default=None, help="FASTA 中的记录 ID（默认第一条）")

    parser.add_argument(
        "--min-aa",
        type=int,
        default=MIN_ORF_AA,
        help=f"最小 ORF 氨基酸长度（默认 {MIN_ORF_AA}）",
    )
    tata_group = parser.add_mutually_exclusive_group()
    tata_group.add_argument(
        "--require-tata",
        dest="require_tata",
        action="store_true",
        help="只保留 ORF 上游带 TATAAA 的候选",
    )
    tata_group.add_argument(
        "--no-require-tata",
        dest="require_tata",
        action="store_false",
        help="不要求 TATA 盒（默认，TATA 只作辅助标记）",
    )
    parser.set_defaults(require_tata=REQUIRE_TATA_BOX)

    repeat_group = parser.add_mutually_exclusive_group()
    repeat_group.add_argument(
        "--remove-repeats",
        dest="remove_repeats",
        action="store_true",
        help="剔除被标记的串联重复区域（默认不剔除）",
    )
    repeat_group.add_argument(
        "--keep-repeats",
        dest="remove_repeats",
        action="store_false",
        help="只标记、不剔除重复区域（默认）",
    )
    parser.set_defaults(remove_repeats=REMOVE_REPEAT_REGIONS)

    parser.add_argument(
        "--max-candidates",
        type=int,
        default=MAX_CANDIDATES_TO_ANALYZE,
        help=f"最多完整分析的候选数（默认 {MAX_CANDIDATES_TO_ANALYZE}）",
    )
    parser.add_argument("--skip-uniprot", action="store_true", help="跳过 UniProt 网络查询")
    parser.add_argument(
        "--no-similarity",
        action="store_true",
        help="精确匹配失败时不做 GenomeNet BLASTP 相似蛋白检索（较快）",
    )
    return parser.parse_args(argv)


def main(argv: Optional[list[str]] = None) -> int:
    setup_console()
    args = parse_args(argv)

    if args.min_aa < 0:
        print("[错误] --min-aa 不能为负数。", file=sys.stderr)
        return 2
    if args.max_candidates < 1:
        print("[错误] --max-candidates 必须 ≥ 1。", file=sys.stderr)
        return 2
    if args.record_id and not args.fasta:
        print("[错误] --record-id 只能与 --fasta 一起使用。", file=sys.stderr)
        return 2

    try:
        if args.dna is not None:
            dna, input_warnings = sanitize_dna(args.dna)
            source = "命令行 --dna 原始 DNA 字符串"
        elif args.fasta:
            dna, source, input_warnings = load_fasta_sequence(args.fasta, args.record_id)
        else:
            dna, demo_warnings = sanitize_dna(DEMO_DNA_SEQUENCE)
            input_warnings = [
                "命令行使用内置合成教学序列。" if args.demo
                else "未提供 --dna 或 --fasta，使用内置合成教学序列。",
                "演示构建体的人工内含子位于 5'UTR，编码区被设计为无额外 GT 供体位点，"
                "以让简化剪接示例可重复；这不代表真实基因的剪接复杂性。",
            ] + demo_warnings
            source = DEMO_SOURCE_LABEL
    except ValueError as exc:
        print(f"[输入错误] {exc}", file=sys.stderr)
        return 2

    options = RuntimeOptions(
        min_orf_aa=args.min_aa,
        require_tata_box=args.require_tata,
        remove_repeat_regions=args.remove_repeats,
        max_candidates=args.max_candidates,
    )

    uniprot_client: Optional[UniProtClient] = None
    if not args.skip_uniprot:
        uniprot_client = UniProtClient(enable_similarity=not args.no_similarity)

    try:
        run_pipeline(
            dna=dna,
            source=source,
            input_warnings=input_warnings,
            options=options,
            skip_uniprot=args.skip_uniprot,
            uniprot_client=uniprot_client,
        )
    except KeyboardInterrupt:
        print("\n[中断] 用户取消了运行。", file=sys.stderr)
        return 130
    except Exception as exc:
        # 保留完整流程的健壮性：顶层异常给出可读错误，而不是只抛出长 traceback。
        print(f"\n[运行错误] {type(exc).__name__}: {exc}", file=sys.stderr)
        print("如需定位代码问题，可临时在 main 中取消异常捕获或使用 traceback。", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    # 双击 / 无参数运行：优先启动图形界面。
    # 命令行分析仍通过显式参数进入（例如 --demo、--dna、--fasta），从而不破坏
    # 原有可重复的批处理/脚本化用法。若当前环境无法创建 Tk 窗口，则自动回退
    # 到命令行内置演示，保证程序仍可运行。
    if len(sys.argv) == 1:
        try:
            from dna_gui import launch_gui
        except Exception as gui_import_error:
            print(
                "[提示] 图形界面模块导入失败，将回退到命令行演示："
                f"{type(gui_import_error).__name__}: {gui_import_error}",
                file=sys.stderr,
            )
        else:
            try:
                raise SystemExit(launch_gui())
            except SystemExit:
                raise
            except Exception as gui_runtime_error:
                print(
                    "[提示] 图形界面启动失败，将回退到命令行演示："
                    f"{type(gui_runtime_error).__name__}: {gui_runtime_error}",
                    file=sys.stderr,
                )
        raise SystemExit(main(["--demo"]))
    raise SystemExit(main())
