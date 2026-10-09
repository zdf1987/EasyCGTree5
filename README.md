# EasyCGTree 5.0

**An Easy Tool for Constructing Core-Gene Tree** — a pipeline for prokaryotic phylogenomic analysis based on core gene sets,
with the graphical interface **EasyCGTree_GUI 1.0** for Windows, Linux and macOS.

[中文说明见下文](#中文说明)

- Repository on GitHub: https://github.com/zdf1987/EasyCGTree5
- Repository on Gitee (mirror, recommended in mainland China): https://gitee.com/zdf1987/EasyCGTree5

## What it does

EasyCGTree builds genome-based maximum-likelihood trees from microbial genomes (protein and/or DNA sequences in FASTA
format): CDS prediction (Prodigal), search of core genes with profile HMMs (HMMER), screening, alignment (MUSCLE) and
trimming (trimAl), and tree inference by the supermatrix, supertree (ASTRAL) or consensus approach (FastTree, IQ-TREE).
Version 5.0 adds the graphical interface, macOS support, NCBI annotation files as input, trees from nucleotide
alignments, and SNP analyses (core-genome SNPs, SNP trees, strain- and clade-specific SNPs).

## Download

1. **Program package** — from the **Releases** page ([GitHub](https://github.com/zdf1987/EasyCGTree5/releases) /
   [Gitee](https://gitee.com/zdf1987/EasyCGTree5/releases)):

   | System | Package |
   |---|---|
   | Windows 10/11 (64-bit) | `EasyCGTree5-Windows-x64.zip` |
   | Linux (64-bit) | `EasyCGTree5-Linux-x64.tar.gz` |
   | macOS, Apple silicon (M1, M2 ...) | `EasyCGTree5-macOS-arm64.zip` |
   | macOS, Intel | `EasyCGTree5-macOS-x64.zip` |

   Decompress it into a folder without spaces or non-English characters in its path.
   Windows users also need Perl, e.g. [Strawberry Perl](https://strawberryperl.com/).

2. **Profile HMM sets** — the folder `HMM` of the package contains only a README. Download the sets you need from the
   folder [`HMM`](HMM) of this repository (at least `bac120.hmm`, the default) and put them into the folder `HMM` of
   the program. See [HMM/README.txt](HMM/README.txt) for the list of sets.

3. Start **EasyCGTree_GUI**, choose *Help → Check installation*, choose the folder with your genomes and click *Run*.

## Documentation

- Manual: [English](doc/EasyCGTree5_Manual_EN.pdf) · [中文](doc/EasyCGTree5_Manual_CN.pdf)
- Command line: `perl EasyCGTree.pl -help` (and the same for the other scripts)
- Running or building the graphical interface from the source code: see the manual (Section 2.5) and
  [build/README_BUILD.txt](build/README_BUILD.txt)

## Citation

Zhang Dao-Feng, He Wei, Shao Zongze, Ahmed Iftikhar, Zhang Yuqin, Li Wen-Jun & Zhao Zhe. EasyCGTree: a pipeline for
prokaryotic phylogenomic analysis based on core gene sets. *BMC Bioinformatics* 2023 24:390.
https://doi.org/10.1186/s12859-023-05527-2

Contact: zdf1987@163.com

## License

EasyCGTree is free software, released under the GNU General Public License v3.0 ([LICENSE](LICENSE)).
The third-party programs in the folder `bin` are distributed under their own licenses.

---

## 中文说明

**EasyCGTree**（An Easy Tool for Constructing Core-Gene Tree）是基于核心基因集的原核生物系统基因组学分析软件，
5.0 版包含可在Windows、Linux 和 macOS 系统中运行的图形界面 **EasyCGTree_GUI 1.0**。

- GitHub 仓库：https://github.com/zdf1987/EasyCGTree5
- Gitee 仓库（镜像，中国大陆用户推荐）：https://gitee.com/zdf1987/EasyCGTree5

### 功能

以 FASTA 格式的微生物基因组（蛋白和/或 DNA 序列）为输入，依次完成 CDS 预测（Prodigal）、用 profile HMM 检索核心基因
（HMMER）、筛选、比对（MUSCLE）与修剪（trimAl），并用超矩阵法、超树法（ASTRAL）或一致树法（FastTree、IQ-TREE）构建
最大似然树。5.0 版新增图形界面、macOS 支持、NCBI 注释文件输入、基于核酸比对建树，以及 SNP 分析（核心基因组 SNP、
SNP 树、菌株与分支特异 SNP统计）。

### 下载

1. **程序包**：从 **Releases** 页面下载（[GitHub](https://github.com/zdf1987/EasyCGTree5/releases) /
   [Gitee](https://gitee.com/zdf1987/EasyCGTree5/releases)），文件名见上文英文板块。解压到路径中不含空格和中文字符的文件夹。
   Windows 用户还需安装 Perl，如 [Strawberry Perl](https://strawberryperl.com/)。
2. **Profile HMM 集**：程序包的 `HMM` 文件夹中只有说明文件。请从本仓库的 [`HMM`](HMM) 文件夹下载所需的 HMM 集
   （至少下载默认的 `bac120.hmm`），放入程序的 `HMM` 文件夹。HMM 集列表见 [HMM/README.txt](HMM/README.txt)。
3. 启动 **EasyCGTree_GUI**，选择“帮助 → 检查安装”，再选择基因组所在的文件夹，点击“运行”。

### 文档

- 使用手册：[English](doc/EasyCGTree5_Manual_EN.pdf) · [中文](doc/EasyCGTree5_Manual_CN.pdf)
- 命令行帮助：`perl EasyCGTree.pl -help`（其他脚本也是-help）
- 从源代码运行或编译图形界面：见手册第 2.5 节及 [build/README_BUILD.txt](build/README_BUILD.txt)

### 引用

见上文 Citation。联系方式：zdf1987@163.com

### 许可证

EasyCGTree 是自由软件，以 GNU 通用公共许可证第 3 版（GPL-3.0，见 [LICENSE](LICENSE)）发布。`bin` 文件夹中的第三方程序遵循其各自的许可证。
