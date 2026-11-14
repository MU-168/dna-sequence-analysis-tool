# DNA 序列分析教学工具

一个基于 **Biopython + Tkinter** 的 DNA 序列分析教学程序，用图形界面演示从原始 DNA 到候选基因、pre-mRNA、剪接、蛋白翻译和 UniProt 注释的完整数据流。

> **重要：这是教学简化模拟器，不是科研或临床基因预测软件。结果不能用于论文结论、临床诊断、药物设计或生产决策。详见 [生物局限性](#生物局限性与免责声明)。**

## 功能

- 原始 DNA 字符串或 FASTA 输入，自动生成反向互补链。
- 正链、负链各 3 个读码框，共 6 框扫描 `ATG…终止密码子` 完整 ORF。
- 最小 ORF 氨基酸长度过滤。
- `TATAAA`（TATA 盒）和 `AATAAA`（poly(A) 信号）简易字面标记；可选“只保留带 TATA 盒的候选”。
- 候选区域按基因方向转录为 pre-mRNA（T → U）。
- 基于 RNA `GU…AG`（DNA 的 `GT-AG`）的简化内含子剔除与外显子拼接。
- 串联重复粗筛；默认只标记，可选剔除。
- 成熟 mRNA 从第一个 `AUG` 开始按标准遗传密码翻译，遇同框终止密码子停止。
- 蛋白注释：
  - **精确匹配**：UniProtKB `/uniprotkb/search` 的 MD5 checksum 精确检索。
  - **相似蛋白回退**：没有精确命中时，调用 GenomeNet 的公开 BLASTP 服务（Swiss-Prot 数据库），再从 UniProtKB 拉取命中条目的名称、物种和 FUNCTION 注释。
  - 基于 BLOSUM62 的本地 Smith-Waterman 比对计算一致度和查询覆盖度，结果按“综合匹配度 = 一致度 × 覆盖度”从高到低排列。
- Tkinter 桌面前端：输入弹窗、参数设置、概览、候选表、剪接与翻译、UniProt/相似蛋白、原始序列。
- 导出文本报告、复制成熟 mRNA/蛋白序列、可点击的 UniProt accession。
- 保留原始命令行分析模式。

## 安装

需要 **Python 3.10+**（开发测试环境为 Python 3.13 + Biopython 1.88）。

```bash
python -m pip install -r requirements.txt
```

`tkinter` 是 Python 标准库的一部分。Windows 官方 Python 通常默认包含；Debian/Ubuntu 如果缺失，可安装：

```bash
sudo apt-get install python3-tk
```

## 运行

### 图形界面

Windows 可以双击：

```text
启动DNA分析工具.bat
```

或在任意平台运行：

```bash
py -3 dna_sequence_analysis.py     # Windows
python3 dna_sequence_analysis.py   # macOS / Linux
```

无参数运行会先弹出 DNA 输入窗口。也可以在 Windows/Linux/macOS 的终端中直接启动：

```bash
python3 dna_gui.py
```

输入窗口支持：

- 粘贴原始 DNA（`A/C/G/T/N`）
- 填入内置教学示例
- 从 FASTA 文件载入（多记录时可选择 record ID）
- 设置最小 ORF 长度
- 开启/关闭 TATA 筛选、重复区域剔除
- 开启/关闭 UniProt 查询和相似性 BLAST（默认开启）

### 命令行示例

```bash
# 内置演示序列（只做本地分析，不联网）
python3 dna_sequence_analysis.py --demo --skip-uniprot

# 直接输入 DNA
python3 dna_sequence_analysis.py --dna "ATG..." --min-aa 30

# 从 FASTA 读取
python3 dna_sequence_analysis.py --fasta examples/demo_sequence.fa --skip-uniprot

# 不开启相似性 BLAST，精确未命中时直接返回无匹配
python3 dna_sequence_analysis.py --fasta examples/demo_sequence.fa --no-similarity

# 查看全部参数
python3 dna_sequence_analysis.py --help
```

## 输入与输出

### 输入

- 原始 DNA：仅 `A/C/G/T/N`，不能带 FASTA 头；RNA 的 `U` 会被拒绝。
- FASTA：默认分析第一条记录，可用 `--record-id` 指定。
- 内置示例是一个合成教学构建体，ORF 编码人 metallothionein-1H（UniProt P80294），便于验证精确匹配。

### 图形界面标签页

| 标签页 | 内容 |
| --- | --- |
| 概览 | 序列长度、GC、ORF 数、候选数、内含子、蛋白长度、UniProt/BLAST 摘要 |
| ORF 候选 | 通过长度阈值的候选表，含链、读码框、输入坐标、aa 长度、TATA/AATAAA |
| 剪接与翻译 | pre-mRNA、重复标记、内含子/外显子、成熟 mRNA、蛋白序列 |
| UniProt / 相似蛋白 | 精确命中或 BLAST 相似蛋白排名、名称、物种、功能、匹配度 |
| 原始序列 | 正向 DNA、反向互补链和输入提示 |

### 命令行输出

CLI 会按输入信息、ORF 与候选、候选 pre-mRNA、重复筛查、剪接、翻译、UniProt/相似蛋白分节打印完整结果。

## 网络与数据说明

程序涉及两个外部公开服务：

1. **UniProtKB REST API**
   - 精确匹配：`https://rest.uniprot.org/uniprotkb/search`
   - 当前接口不支持 `sequence:` 检索字段；程序使用 `checksum:<MD5> AND length:[N TO N]` 做全长精确检索。
   - 仅提交翻译后的蛋白序列，不提交原始 DNA。

2. **GenomeNet BLASTP**
   - 网页服务：<https://www.genome.jp/tools/blast/>
   - 仅当 UniProt 精确未命中且相似性搜索开启时触发。
   - 数据库：Swiss-Prot。
   - 程序把蛋白序列提交给 GenomeNet，解析命中 accession，再从 UniProtKB 获取注释。
   - 这是外部公开网页服务，可能排队、限流或改变页面结构；程序有超时、间隔控制和错误捕获，但不能保证服务持续可用。

相似性排序：

```text
综合匹配度 = 一致度 × 查询覆盖度
```

排序优先按综合匹配度，其次按一致度、覆盖度，并优先 reviewed 条目。最多展示 8 条相似蛋白。这样可以避免“只匹配到很短的 100% 一致片段”被错误排到全长高相似蛋白前面。

## 配置

主要参数集中在 [dna_sequence_analysis.py](./dna_sequence_analysis.py) 顶部的“集中配置区”，包括：

- `MIN_ORF_AA`
- `REQUIRE_TATA_BOX`
- `CANDIDATE_UPSTREAM_FLANK` / `CANDIDATE_DOWNSTREAM_FLANK`
- `SPLICE_MIN_INTRON_LEN` / `SPLICE_MAX_INTRON_LEN`
- `REMOVE_REPEAT_REGIONS`
- `UNIPROT_SIMILARITY_ENABLED`
- `UNIPROT_BLAST_*`
- `UNIPROT_SIMILARITY_MIN_IDENTITY` / `UNIPROT_SIMILARITY_MIN_COVERAGE`
- `UNIPROT_MAX_HITS_TO_PRINT`

## 文件结构

```text
dna-analysis-tool/
├── dna_sequence_analysis.py   # 分析内核、命令行入口
├── dna_gui.py                 # Tkinter 图形界面
├── 启动DNA分析工具.bat        # Windows 双击启动器
├── examples/
│   └── demo_sequence.fa       # 内置演示序列的 FASTA
├── requirements.txt
├── run_gui.sh                 # macOS/Linux 启动脚本
├── LICENSE                    # MIT License
├── .gitignore
├── .gitattributes
└── README.md
```

GUI 与命令行共用同一套 ORF、剪接、翻译和注释内核。

## 快速验证

内置演示序列预期关键结果：

- 输入长度：424 bp
- 通过长度阈值的候选：1 个
- ORF：61 aa
- 识别 1 个 75 nt 内含子
- 蛋白与 P80294（Metallothionein-1H）一致
- UniProt 精确检索命中 P80294

将蛋白序列做单氨基酸突变后，精确 checksum 不再命中；相似性搜索会返回多条 metallothionein 同家族蛋白，并按综合匹配度排序（实际结果受 GenomeNet 状态影响）。

## 生物局限性与免责声明

本程序是教学演示，不是专业基因预测工具。请务必了解：

1. **ORF 扫描**只依据 ATG 起始、同框终止，不考虑启动子模型、Kozak 序列、剪接位点统计模型、RNA-seq 证据、物种参数或内含子对 ORF 的打断。与 AUGUSTUS、BRAKER、GENSCAN、GeneMark 等专业基因预测器不是同一类工具。
2. **TATAAA/AATAAA** 只是六碱基字面标记，不验证真实启动子活性或 poly(A) 位点。
3. **GT-AG 剪接**只做简单的 GU…AG 搜索，忽略分支点、多嘧啶区、剪接增强子/沉默子、可变剪接和次要剪接类型。
4. **串联重复筛查**不等于转座子鉴定，只作疑似区域标记。
5. **翻译**从第一个 AUG 开始，使用标准遗传密码，不考虑 Kozak、uORF、IRES、线粒体遗传密码、RNA 编辑或翻译后修饰。
6. **UniProt 精确命中**并不证明输入的 DNA 在真实细胞中表达了该蛋白。
7. **相似蛋白 BLAST** 只是基于序列相似性的功能线索；命中条目的 FUNCTION 注释不能自动归因到候选 ORF。
8. 不要把本程序的任何结果用于正式科研、临床诊断、药物设计、育种或生产决策。

真实研究建议使用经过验证的基因预测、序列比对、结构域注释和 RNA-seq/蛋白证据流程，并由具备生物学背景的人员审查。

## License

本项目采用 [MIT License](./LICENSE)。

你可以自由使用、修改和再分发本项目的代码，但需保留版权声明和许可证文本。程序输出仅用于教学演示，不构成科研结论或临床建议。
