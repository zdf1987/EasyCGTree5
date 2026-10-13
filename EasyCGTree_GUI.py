#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
EasyCGTree_GUI - graphical interface for EasyCGTree 5.0
Version 1.0 by Dao-Feng Zhang

The GUI only collects the parameters, shows the command lines and runs the Perl scripts of EasyCGTree;
all the analysis is done by the scripts. Requirements: Python 3.9+ and PySide6 ('pip install PySide6'),
and Perl (Linux/macOS: included; Windows: e.g. Strawberry Perl).

Run:  python EasyCGTree_GUI.py
"""

import sys
import os
import re
import json
import shlex
import shutil
import signal
import codecs
import threading
import subprocess
from pathlib import Path

try:
    from PySide6.QtCore import Qt, QObject, Signal, QSettings, QUrl, QTimer, QRect, QSize, QPoint, QByteArray
    from PySide6.QtGui import (QAction, QActionGroup, QDesktopServices, QFontDatabase, QTextCursor,
                               QGuiApplication, QColor, QBrush, QIcon, QPixmap, QPalette)
    from PySide6.QtWidgets import (
        QApplication, QMainWindow, QWidget, QTabWidget, QVBoxLayout, QHBoxLayout, QGridLayout, QGroupBox,
        QLabel, QLineEdit, QPushButton, QToolButton, QComboBox, QSpinBox, QDoubleSpinBox, QCheckBox,
        QPlainTextEdit, QFileDialog, QMessageBox, QSplitter, QTreeView, QTableWidget, QTableWidgetItem,
        QDialog, QDialogButtonBox, QScrollArea, QRadioButton, QButtonGroup, QListWidget, QListWidgetItem,
        QHeaderView, QMenu, QFrame, QAbstractItemView, QLayout, QStackedWidget, QStyleFactory)
    try:
        from PySide6.QtWidgets import QFileSystemModel
    except ImportError:  # Qt 6: QFileSystemModel is in QtGui
        from PySide6.QtGui import QFileSystemModel
except ImportError:
    sys.stderr.write("EasyCGTree_GUI needs PySide6. Please install it with:\n    pip install PySide6\n")
    sys.exit(1)


GUI_VERSION = "1.0"
IS_WIN = sys.platform.startswith("win")
IS_MAC = sys.platform == "darwin"
EXE = ".exe" if IS_WIN else ""

# ======================================================================================================
# About page. This block can be edited freely (e.g. when the GUI is published).
# ======================================================================================================
ABOUT = {
    "title": "EasyCGTree_GUI",
    "authors": "Dao-Feng Zhang, Wei He, Zongze Shao, Iftikhar Ahmed, Yuqin Zhang, Wen-Jun Li, Zhe Zhao",
    "citation": ("Zhang DF, He W, Shao Z, Ahmed I, Zhang Y, Li WJ, Zhao Z. EasyCGTree: a pipeline for prokaryotic "
                 "phylogenomic analysis based on core gene sets. BMC Bioinformatics. 2023;24:390."),
    "doi": "https://doi.org/10.1186/s12859-023-05527-2",
    "url": "https://github.com/zdf1987/EasyCGTree5",
    "url2": "https://gitee.com/zdf1987/EasyCGTree5",
    "contact": "zdf1987@163.com",
}
# The HMM sets are downloaded separately from the folder 'HMM' of the repositories.
REPO_HMM = {"GitHub": "https://github.com/zdf1987/EasyCGTree5/tree/main/HMM",
            "Gitee": "https://gitee.com/zdf1987/EasyCGTree5/tree/main/HMM"}

# ======================================================================================================
# Translations: key -> (English, Chinese)
# ======================================================================================================
TR = {
    # window, menus
    "app_title": ("EasyCGTree_GUI", "EasyCGTree_GUI"),
    "menu_file": ("&File", "文件(&F)"),
    "menu_load": ("Load parameters...", "载入参数..."),
    "menu_save": ("Save parameters...", "保存参数..."),
    "menu_exit": ("Exit", "退出"),
    "menu_settings": ("&Settings", "设置(&S)"),
    "menu_paths": ("Paths and programs...", "路径与程序..."),
    "menu_treeopts": ("Tree program command lines...", "建树程序命令行..."),
    "menu_language": ("Language", "语言"),
    "menu_help": ("&Help", "帮助(&H)"),
    "menu_check": ("Check installation...", "检查安装..."),
    "menu_about": ("About...", "关于..."),
    "tab_main": ("Main pipeline", "主流程"),
    "tab_snp": ("SNP analysis", "SNP分析"),
    "tab_hmm": ("HMM building", "HMM构建"),
    # common
    "global": ("Global settings", "全局参数"),
    "status": ("Status", "状态"),
    "files": ("Files", "文件浏览"),
    "cmdline": ("Command line", "命令行"),
    "concise": ("Concise", "简洁版"),
    "full": ("Full", "完整版"),
    "concise_tip": ("Only the options that differ from the defaults", "只显示与默认值不同的参数"),
    "full_tip": ("All options, including the defaults", "显示全部参数，包括默认值"),
    "copy": ("Copy", "复制"),
    "copied": ("Command line copied to the clipboard.", "命令行已复制到剪贴板。"),
    "run": ("Run", "运行"),
    "stop": ("Stop", "停止"),
    "log": ("Log", "运行日志"),
    "clear_log": ("Clear", "清空"),
    "browse": ("Browse...", "浏览..."),
    "choose": ("Choose...", "选择..."),
    "up": ("Up", "上级"),
    "home": ("Output directory", "输出目录"),
    "refresh": ("Refresh", "刷新"),
    "open": ("Open", "打开"),
    "open_system": ("Open with the system application", "用系统程序打开"),
    "copy_path": ("Copy path", "复制路径"),
    "show_folder": ("Show in the file manager", "在文件管理器中显示"),
    "help_title": ("Help", "说明"),
    "previous": ("(previous / default)", "（沿用之前/默认）"),
    "placeholder_previous": ("empty = previous setting", "留空 = 沿用之前的设置"),
    "ok": ("OK", "确定"),
    "cancel": ("Cancel", "取消"),
    "close": ("Close", "关闭"),
    "error": ("Error", "错误"),
    "warning": ("Warning", "警告"),
    "info": ("Information", "提示"),
    "yes": ("Yes", "是"),
    "no": ("No", "否"),
    # main tab
    "input": ("Input directory", "输入文件夹"),
    "h_input": ("Directory with the genomes, one file per genome: genome (DNA) sequences, proteomes, or NCBI annotation "
                "files (the protein and CDS files of one assembly are paired by locus_tag).\n\nThe results are written "
                "to the directory that contains the input directory; working files go to '<input>_TEM'.\n\n"
                "Note: EasyCGTree renames the input files (spaces, '_' and '.' are replaced by '-').",
                "存放基因组的文件夹，每个基因组一个文件：基因组（DNA）序列、蛋白组，或NCBI注释文件（同一组装号的蛋白与CDS文件按"
                "locus_tag配对）。\n\n结果写入输入文件夹所在的目录；中间文件写入'<输入文件夹>_TEM'。\n\n"
                "注意：EasyCGTree会重命名输入文件（空格、'_'和'.'替换为'-'）。"),
    "check_input": ("Check input", "检查输入"),
    "h_check_input": ("Shows how the input files will be grouped and prepared (genome names, sequence types, NCBI files), "
                      "without renaming or running anything ('-dry_run').",
                      "显示输入文件将如何分组和处理（基因组名称、序列类型、NCBI文件），不重命名也不运行任何程序（'-dry_run'）。"),
    "thread": ("Threads", "线程数"),
    "h_thread": ("Number of threads used by HMMER, MUSCLE, FastTree, IQ-TREE and ASTRAL. The maximum is the number of "
                 "CPU cores of this computer; the scripts also lower a larger value themselves (with a warning).",
                 "HMMER、MUSCLE、FastTree、IQ-TREE和ASTRAL使用的线程数。上限为本机CPU核心数；若设置更大，脚本也会自动降为核心数并给出提示后继续运行。"),
    "task": ("Task", "任务"),
    "h_task": ("Which part of the pipeline is run. Tick the tasks; they must follow each other, so ticking two tasks also "
               "ticks the ones between them. All five ticked = 'all' (one command); otherwise one EasyCGTree command per task is made, and the commands are run one after the other. The later tasks use the results of the earlier ones in "
               "'<input>_TEM' (see 'Status'). Options that the chosen tasks do not use are disabled.\n\n"
               "predict: CDS prediction; hmmsearch: HMM search (reuses the prediction, so several HMM sets can be tried); "
               "refine: screening of the hits; alignment: alignment and trimming; tree_infer: tree inference.\n\n"
               "If a prediction made with the same input files and settings already exists, a notice is shown and 'predict' "
               "is unticked; close the notice to tick it again.",
               "运行流程的哪一部分。勾选所需任务；任务必须连续，勾选两个任务时会自动勾选中间的任务。五个全选即为'all'（一条命令）；否则每个任务生成一条EasyCGTree命令，依次运行。后面的任务使用"
               "'<输入文件夹>_TEM'中前面任务的结果（见'状态'）。所选任务用不到的参数会变灰。\n\n"
               "predict：CDS预测；hmmsearch：HMM检索（复用已有预测，可尝试多个HMM集）；refine：筛选检索结果；"
               "alignment：比对与修剪；tree_infer：建树。\n\n"
               "如果已存在用相同输入文件和设置得到的预测结果，会显示提示并自动取消勾选'predict'；关闭提示后才能重新勾选。"),
    "task_all": ("all (complete pipeline)", "all（完整流程）"),
    "task_predict": ("predict (CDS prediction)", "predict（CDS预测）"),
    "task_hmmsearch": ("hmmsearch (HMM search)", "hmmsearch（HMM检索）"),
    "task_refine": ("refine (screening)", "refine（筛选）"),
    "task_alignment": ("alignment (alignment and trimming)", "alignment（比对与修剪）"),
    "task_tree_infer": ("tree_infer (tree inference)", "tree_infer（建树）"),
    "runmode": ("Run", "运行内容"),
    "h_runmode": ("Run the EasyCGTree pipeline, or only summarize the gene prevalence (Gene_Prevelence.pl).",
                  "运行EasyCGTree主流程，或只汇总基因分布情况（Gene_Prevelence.pl）。"),
    "mode_pipeline": ("EasyCGTree pipeline", "EasyCGTree主流程"),
    "mode_prevalence": ("Gene prevalence only", "仅基因分布汇总"),
    "g_predict": ("1. CDS prediction", "1. CDS预测"),
    "seq": ("Sequence type", "序列类型"),
    "h_seq": ("Sequences used for the alignments and the tree: protein (prot) or nucleotide (nucl).\n\n'nucl' needs DNA "
              "for all genomes (genome sequences or NCBI CDS files).",
              "用于比对和建树的序列：蛋白（prot）或核酸（nucl）。\n\n'nucl'要求所有基因组都提供DNA序列（基因组序列或NCBI CDS文件）。"),
    "seq_prot": ("prot (protein)", "prot（蛋白）"),
    "seq_nucl": ("nucl (nucleotide)", "nucl（核酸）"),
    "keep_CDS_nucl": ("Keep CDS DNA sequences", "保留CDS核酸序列"),
    "h_keep_CDS_nucl": ("Also keep the DNA sequences of the CDS predicted by Prodigal, e.g. for a later nucleotide tree "
                        "or SNP analysis. Automatic with 'nucl'.",
                        "同时保留Prodigal预测的CDS核酸序列，便于之后做核酸树或SNP分析。选择'nucl'时自动保留。"),
    "g_hmmsearch": ("2. HMM search", "2. HMM检索"),
    "hmm": ("Profile HMM set", "HMM基因集"),
    "h_hmm": ("Core gene set used to find the genes: a set from the 'HMM' directory of EasyCGTree (see '../HMM/README.txt' "
              "for more details). The HMM sets are downloaded separately from the folder 'HMM' of the EasyCGTree5 repository "
              "on GitHub or Gitee (Help - HMM sets online). New sets can be made in the 'HMM building' tab.",
              "用于检索核心基因的HMM集：EasyCGTree的'HMM'目录中的集合（查看'../HMM/README.txt'获取更多信息）。HMM集需从GitHub或Gitee上"
              "EasyCGTree5仓库的'HMM'文件夹单独下载（帮助 - 在线下载HMM集）。新的集合可在'HMM构建'页面中构建。"),
    "evalue": ("E-value", "E值"),
    "h_evalue": ("E-value threshold of the hmmsearch hits. A hit must also score at least 1/4 of the HMM length.",
                 "hmmsearch命中的E值阈值。命中的得分还需不低于HMM长度的1/4。"),
    "g_refine": ("3. Screening of the hits", "3. 检索结果筛选"),
    "genome_cutoff": ("Genome cutoff", "基因组阈值"),
    "h_genome_cutoff": ("Genomes with fewer genes than (largest number of genes found in one genome x cutoff) are excluded. "
                        "Lower it to keep incomplete genomes. Recommended 0.5-1.",
                        "所含基因数少于（单个基因组中检出的最大基因数×阈值）的基因组被排除。降低阈值可保留不完整的基因组。推荐0.5-1。"),
    "gene_cutoff": ("Gene cutoff", "基因阈值"),
    "h_gene_cutoff": ("Genes found in fewer than (number of selected genomes x cutoff) genomes are excluded. "
                      "With a consensus tree (cs) EasyCGTree always uses 1; with a supertree (st) 1 is recommended. Recommended range 0.5-1.",
                      "在少于（入选基因组数×阈值）个基因组中检出的基因被排除。一致树（cs）时EasyCGTree总是使用1，超树（st）时建议使用1。推荐范围0.5-1。"),
    "gene_cutoff_cs": ("(1 with a consensus tree)", "（一致树时为1）"),
    "g_alignment": ("4. Alignment and trimming", "4. 比对与修剪"),
    "trim": ("trimAl method", "trimAl方法"),
    "h_trim": ("Method of trimAl for trimming the alignments: gappyout (removes columns rich in gaps), strict, strictplus "
               "(most stringent), nogaps (removes every column that contains a gap; genes without any column left are "
               "skipped). The sequence type is set in step 1.",
               "trimAl修剪比对的方法：gappyout（去除空位较多的列）、strict、strictplus（最严格）、nogaps（去除所有含空位的列；"
               "修剪后不剩任何列的基因会被跳过）。序列类型在第1步中设置。"),
    "g_tree": ("5. Tree inference", "5. 建树"),
    "tree": ("Tree approach", "建树方法"),
    "h_tree": ("sm: supermatrix (concatenated alignment, one tree); st: supertree (gene trees combined by ASTRAL), a gene "
               "cutoff of 1 is recommended; cs: majority-rule consensus of the gene trees (by IQ-TREE), requires a gene cutoff "
               "of 1 (set automatically).",
               "sm：超矩阵（串联比对，建一棵树）；st：超树（用ASTRAL整合基因树），建议将基因阈值设为1（不是1也能得到结果）；cs：基因树的多数一致树（由IQ-TREE生成），要求将基因阈值设为1（自动设置）。"),
    "tree_sm": ("sm (supermatrix)", "sm（超矩阵）"),
    "tree_st": ("st (supertree)", "st（超树）"),
    "tree_cs": ("cs (consensus tree)", "cs（一致树）"),
    "tree_app": ("Tree program", "建树程序"),
    "h_tree_app": ("FastTree (fast) or IQ-TREE (model selection, bootstrap; slower). The command lines can be changed in "
                   "Settings - Tree program command lines.",
                   "FastTree（快速）或IQ-TREE（模型选择、自举检验；较慢）。命令行可在'设置-建树程序命令行'中修改。"),
    "g_prev_main": ("Gene prevalence (Gene_Prevelence.pl)", "基因分布汇总（Gene_Prevelence.pl）"),
    "prev_after": ("Summarize afterwards", "运行后汇总"),
    "h_prev_after": ("After the pipeline, summarize which genes of the HMM set were found in which genome, their copy "
                     "numbers, and how many genomes/genes pass the cutoffs ('<input>.<hmm>.gene_*.txt').",
                     "主流程结束后，汇总HMM集中的基因在各基因组中的检出情况、拷贝数，以及在阈值下保留的基因组/基因数"
                     "（'<输入文件夹>.<hmm>.gene_*.txt'）。"),
    "prev_force": ("Search again if needed (-force)", "必要时重新检索（-force）"),
    "h_prev_force": ("If the existing HMM search used another HMM set or a smaller E-value, run the search again. This "
                     "removes the later results (screening, alignments, trees) in '<input>_TEM'.",
                     "如果已有的HMM检索使用了其他HMM集或更小的E值，则重新检索。这会删除'<输入文件夹>_TEM'中后续步骤的结果"
                     "（筛选、比对、树）。"),
    # SNP tab
    "h_snp_input": ("Input directory used with EasyCGTree. All genomes must be DNA (genome sequences or NCBI CDS files). "
                    "If the nucleotide alignments are missing, EasyCGTree is run first with '-seq nucl'.",
                    "EasyCGTree使用的输入文件夹。所有基因组都必须是DNA（基因组序列或NCBI CDS文件）。如果缺少核酸比对，会先以'-seq nucl'运行EasyCGTree。"),
    "g_ecg": ("EasyCGTree settings (used only if alignments have to be made)", "EasyCGTree参数（仅在需要生成比对时使用）"),
    "h_g_ecg": ("Used only if the nucleotide alignments are missing or were made with other settings; then the needed "
                "EasyCGTree tasks are run with '-seq nucl'. The fields are filled from the record of the previous EasyCGTree "
                "run in the output directory, or with the defaults if there is none.",
                "仅在缺少核酸比对或比对所用设置不同时使用；此时会以'-seq nucl'运行所需的EasyCGTree步骤。各项自动填入输出目录中"
                "上一次EasyCGTree运行记录里的设置；没有运行记录时填入默认值。"),
    "g_snp1": ("1. SNP extraction (EasyCGTree_SNP.pl)", "1. SNP提取（EasyCGTree_SNP.pl）"),
    "ignore": ("Ignored strains", "忽略的菌株"),
    "h_ignore": ("These strains are not considered when SNP sites are detected (sites that vary only because of them are "
                 "not SNPs), but they stay in the SNP alignment and the tree. At least two strains must stay not ignored.",
                 "检测SNP位点时不考虑这些菌株（仅因它们而变化的位点不算SNP），但它们仍保留在SNP比对和树中。至少要保留2个菌株不被忽略。"),
    "max_missing": ("Max. missing fraction", "最大缺失比例"),
    "h_max_missing": ("Largest fraction of the (not ignored) strains with a gap or ambiguous base at a SNP site. 0: only "
                      "core SNPs present in all strains.",
                      "SNP位点上允许有空位或简并碱基的（未忽略）菌株的最大比例。0：只取所有菌株都有碱基的核心SNP。"),
    "aln": ("Alignments", "比对"),
    "h_aln": ("trimmed: the alignments after trimAl (as used for the EasyCGTree trees); original: before trimming.",
              "trimmed：trimAl修剪后的比对（与EasyCGTree建树所用相同）；original：修剪前的比对。"),
    "aln_trimmed": ("trimmed", "trimmed（修剪后）"),
    "aln_original": ("original", "original（修剪前）"),
    "ref": ("Reference strain", "参考菌株"),
    "h_ref": ("Optional, one strain. The positions of the SNPs in the CDS of this strain and the codon positions are reported too.",
              "可选，只能选1个菌株。同时给出SNP在该菌株CDS中的位置和密码子位置。"),
    "g_snp2": ("2. SNP tree", "2. SNP建树"),
    "h_snp_tree_app": ("IQ-TREE (default) or FastTree. FastTree has no correction for using only variable sites "
                       "(branch lengths are substitutions per SNP site).",
                       "IQ-TREE（默认）或FastTree。FastTree无法校正只使用变异位点带来的偏差（枝长为每个SNP位点的替换数）。"),
    "asc": ("Correction", "校正方式"),
    "h_asc": ("fconst: the invariant A/C/G/T sites of the analysed alignment are given to IQ-TREE, so branch lengths are "
              "substitutions per site of the core genes; lewis: ascertainment bias correction of the model (+ASC).",
              "fconst：把分析区域中不变位点的A/C/G/T计数提供给IQ-TREE，枝长为每个核心基因位点的替换数；lewis：模型的确定偏差校正（+ASC）。"),
    "g_snp3": ("3. Strain- and clade-specific SNPs (EasyCGTree_SpecificSNP.pl)", "3. 菌株与分支特异SNP（EasyCGTree_SpecificSNP.pl）"),
    "mode": ("Definition", "判定方式"),
    "h_mode": ("exclusive: all strains of the group share a base that no other strain has.\nstrict: in addition all other "
               "strains share one base (only sites without missing data); counts do not depend on the rooting.",
               "exclusive：组内所有菌株共有一种其他菌株都没有的碱基。\nstrict：此外，其余菌株也都是同一种碱基（只用无缺失的位点）；计数与树的定根无关。"),
    "outgroup": ("Outgroup", "外群"),
    "h_outgroup": ("Optional. Root the tree with these strains (they must form a clade) before counting the clade-specific "
                   "SNPs. Recommended with 'exclusive'. With another tree, the strains of that tree are offered.",
                   "可选。先用这些菌株给树定根（它们须构成一个分支），再统计分支特异SNP。使用'exclusive'时推荐设置。指定了其他树时，列出的是该树中的菌株。"),
    "other_tree": ("Other tree", "其他树文件"),
    "h_other_tree": ("Optional. Annotate this tree instead of the SNP tree (e.g. the supermatrix tree of EasyCGTree): a "
                     "phylogenetic tree in Newick format whose strain names are exactly those of the SNP alignment. It should "
                     "be rooted as wanted (e.g. with the outgroup) beforehand, or give the outgroup strains below and the tree "
                     "is rerooted with them. The rooting matters only for 'exclusive'; 'strict' does not depend on it.\n\n"
                     "With another tree, step 2 (SNP tree) is not needed and is switched off, and step 1 is run only if needed: "
                     "if there are no SNP results yet, or their settings differ from the ones above.",
                     "可选。标注这棵树而不是SNP树（如EasyCGTree的超矩阵树）：Newick格式的进化树，菌株名须与SNP比对中的完全一致。"
                     "树应事先按需要定好根（如以外群为根），或在下面的'外群'中指定外群，由脚本重新定根。定根只影响'exclusive'的结果，"
                     "'strict'与根无关。\n\n指定其他树后，不再需要第2步（SNP建树），该步被关闭；第1步只在需要时运行："
                     "还没有SNP结果，或上面的设置与已有SNP结果的不同。"),
    "msg_no_tree_names": ("No strain names were found in the given tree file. Please check that it is a Newick tree.",
                          "在所给的树文件中没有找到菌株名。请确认它是Newick格式的树。"),
    "msg_tree_names": ("The strains of the given tree differ from those of the SNP alignment: {n}",
                       "所给树中的菌株与SNP比对中的不一致：{n}"),
    "auto1_run": ("Run automatically: {why}.", "将自动运行：{why}。"),
    "auto1_skip": ("Not run: the existing SNP results were made with these settings and are used.",
                   "不运行：已有的SNP结果正是用这些设置得到的，直接使用。"),
    "why_no_predict": ("no CDS prediction yet", "还没有CDS预测"),
    "why_no_fnn": ("the CDS DNA sequences were not kept", "没有保留CDS核酸序列"),
    "why_no_step": ("'{s}' has not been run yet", "还没有运行'{s}'"),
    "why_diff": ("{k} {v} (before: {o})", "{k}为{v}（之前：{o}）"),
    "why_no_nucl": ("no nucleotide alignments after the last screening", "上次筛选后没有核酸比对"),
    "why_snp_ecg": ("EasyCGTree is run again", "EasyCGTree需要重新运行"),
    "why_snp_none": ("no SNP results yet", "还没有SNP结果"),
    "why_snp_old": ("the alignments are newer than the SNP results", "比对结果比SNP结果新"),
    "plan_ecg": ("Compared with the previous run: {why}. EasyCGTree is run again from '{t}' (with '-seq nucl'), and the "
                 "SNPs are extracted again.",
                 "与上次运行相比：{why}。将从'{t}'开始重新运行EasyCGTree（'-seq nucl'），并重新提取SNP。"),
    "plan_ecg_tem": ("The results in '<input>_TEM' from this step on are replaced; the folder is shared with the Main "
                     "pipeline tab.",
                     "'<输入文件夹>_TEM'中从这一步开始的结果会被替换；该文件夹与主流程页共用。"),
    "plan_snp_diff": ("The SNP settings differ from the last extraction: {why}.", "SNP设置与上次提取不同：{why}。"),
    "plan_no_snp": ("No SNP results found: EasyCGTree_SNP.pl will be run first with its default settings (tick step 1 "
                    "to use the settings above).",
                    "没有找到SNP结果：将先以默认设置运行EasyCGTree_SNP.pl（勾选第1步才会使用上面的设置）。"),
    "plan_tree_names": ("The strains of the given tree are not those of the SNP alignment. Only in the tree: {t}. "
                        "Only in the alignment: {a}.",
                        "所给树中的菌株与SNP比对中的不一致。仅在树中：{t}。仅在比对中：{a}。"),
    "plan_backup": ("Result files with the same name (trees, SNP tables) are kept: the old ones are renamed to "
                    "'<file>.<date>.bak'.",
                    "同名的结果文件（树、SNP表格）会保留：旧文件改名为'<文件名>.<日期>.bak'。"),
    "g_prev_snp": ("Gene prevalence (Gene_Prevelence.pl)", "基因分布汇总（Gene_Prevelence.pl）"),
    "need_step": ("Please select at least one step.", "请至少选择一个步骤。"),
    # HMM tab
    "gc": ("Gene family directory", "基因家族文件夹"),
    "h_gc": ("Directory with one FASTA file of protein sequences per gene family. The file name (without extension) "
             "becomes the gene name. The input files are not changed.",
             "每个基因家族一个蛋白序列FASTA文件的文件夹。文件名（不含扩展名）即基因名。不会修改输入文件。"),
    "h_thread_hmm": ("Number of threads used by MUSCLE.", "MUSCLE使用的线程数。"),
    "g_hmm_aln": ("Alignment", "比对"),
    "hmm_aln": ("Already aligned", "已比对"),
    "h_hmm_aln": ("The files are aligned FASTA files; they are not aligned again.", "输入文件已是比对好的FASTA文件，不再重新比对。"),
    "super5": ("'super5' from (sequences)", "使用'super5'的序列数"),
    "h_super5": ("Families with at least this number of sequences are aligned with 'muscle -super5' (faster, for large "
                 "families); smaller ones with 'muscle -align'.",
                 "序列数不少于该值的家族用'muscle -super5'比对（较快，适合大家族）；其余用'muscle -align'。"),
    "g_hmm_out": ("Output", "输出"),
    "hmm_name": ("Name of the HMM set", "HMM集名称"),
    "h_hmm_name": ("The set is written to 'HMM/<name>.hmm' in the EasyCGTree directory and used as '-hmm <name>'. "
                   "Filled in with the name of the gene family directory (characters other than letters, digits, '.', '_' "
                   "and '-' become '_'); it can be changed.",
                   "写入EasyCGTree目录下的'HMM/<名称>.hmm'，使用时为'-hmm <名称>'。自动填入基因家族文件夹的名称（字母、数字、'.'、'_'、'-'以外的字符改为'_'），可以修改。"),
    "hmm_force": ("Overwrite an existing set", "覆盖同名集合"),
    "h_hmm_force": ("Overwrite 'HMM/<name>.hmm' if it exists. Otherwise an existing set is never overwritten.",
                    "如果'HMM/<名称>.hmm'已存在则覆盖。否则不会覆盖已有的集合。"),
    "hmm_sets": ("HMM sets in the EasyCGTree directory", "EasyCGTree目录中的HMM集"),
    "hmm_dir_button": ("HMM directory", "HMM目录"),
    # status texts
    "st_no_input": ("No input directory selected.", "尚未选择输入文件夹。"),
    "st_no_tem": ("No results yet ('{tem}' does not exist).", "尚无结果（'{tem}'不存在）。"),
    "st_tem": ("Working directory: {tem}", "工作目录：{tem}"),
    "st_old": ("No record file (results of EasyCGTree < 5.0?).", "无记录文件（EasyCGTree 5.0以前的结果？）。"),
    "st_predict": ("1. CDS prediction: {v}", "1. CDS预测：{v}"),
    "st_hmmsearch": ("2. HMM search: {v}", "2. HMM检索：{v}"),
    "st_refine": ("3. Screening: {v}", "3. 筛选：{v}"),
    "st_alignment": ("4. Alignments: {v}", "4. 比对：{v}"),
    "st_trees": ("Trees: {v}", "树文件：{v}"),
    "st_snp": ("SNP results: {v}", "SNP结果：{v}"),
    "st_done": ("done", "已完成"),
    "st_not_done": ("not done", "未完成"),
    "st_none": ("none", "无"),
    "st_genomes": ("{n} genomes", "{n}个基因组"),
    "st_running": ("Running...", "运行中..."),
    "st_finished": ("Finished.", "完成。"),
    "st_failed": ("Stopped with an error (exit code {rc}). See the log.", "出错停止（退出码{rc}），请查看日志。"),
    "st_stopped": ("Stopped by the user.", "已被用户停止。"),
    "st_hmm_incompatible": ("cannot be read by the hmmsearch in bin (HMMER {v})", "bin中的hmmsearch（HMMER {v}）无法读取"),
    # messages
    "msg_required": ("Please set '{label}'.", "请设置'{label}'。"),
    "msg_bad_number": ("'{label}': '{v}' is not a valid number.", "'{label}'：'{v}'不是有效的数字。"),
    "msg_range": ("'{label}' must be between {a} and {b}.", "'{label}'应在{a}和{b}之间。"),
    "msg_no_perl": ("Perl was not found. EasyCGTree needs Perl.\n\nWindows: install Strawberry Perl "
                    "(https://strawberryperl.com), or set the Perl program in Settings - Paths and programs.\n"
                    "Linux/macOS: Perl is normally installed; please check.",
                    "未找到Perl。EasyCGTree需要Perl。\n\nWindows：请安装Strawberry Perl（https://strawberryperl.com），"
                    "或在'设置-路径与程序'中指定Perl程序。\nLinux/macOS：一般已自带Perl，请检查。"),
    "msg_no_home": ("The EasyCGTree directory (with EasyCGTree.pl and 'bin') was not found. Please set it in "
                    "Settings - Paths and programs.",
                    "未找到EasyCGTree目录（含EasyCGTree.pl和'bin'）。请在'设置-路径与程序'中设置。"),
    "msg_no_script": ("The script '{s}' was not found in '{d}'.", "在'{d}'中未找到脚本'{s}'。"),
    "msg_running": ("A job is running. Stop it?", "有任务正在运行，是否停止？"),
    "msg_quit_running": ("A job is running. Stop it and quit?", "有任务正在运行，是否停止并退出？"),
    "msg_saved": ("Parameters saved to '{f}'.", "参数已保存到'{f}'。"),
    "msg_load_failed": ("The parameters could not be loaded: {e}", "无法载入参数：{e}"),
    "msg_no_strains": ("No genome names were found. Please check the input directory (and the output directory).",
                       "未找到基因组名称。请检查输入文件夹（以及输出文件夹）。"),
    "msg_dryrun_failed": ("The input check stopped with an error:", "输入检查出错："),
    "msg_set_input": ("Please choose the input directory first.", "请先选择输入文件夹。"),
    # dialogs
    "dlg_strains": ("Choose strains", "选择菌株"),
    "dlg_dryrun": ("Input check", "输入检查"),
    "dry_genome": ("Genome", "基因组"),
    "dry_type": ("Type", "类型"),
    "dry_files": ("Files", "文件"),
    "dry_ignored": ("Ignored files", "忽略的文件"),
    "dry_summary": ("{n} genomes were found in '{d}'.", "在'{d}'中找到{n}个基因组。"),
    "type_nucl": ("genome (DNA), CDS predicted by Prodigal", "基因组（DNA），由Prodigal预测CDS"),
    "type_prot": ("protein sequences", "蛋白序列"),
    "type_ncbi_pair": ("NCBI protein + CDS, paired by locus_tag", "NCBI蛋白+CDS，按locus_tag配对"),
    "type_ncbi_cds": ("NCBI CDS, translated", "NCBI CDS，翻译得到蛋白"),
    "dlg_settings": ("Paths and programs", "路径与程序"),
    "set_home": ("EasyCGTree directory", "EasyCGTree目录"),
    "h_set_home": ("Directory containing EasyCGTree.pl and the other scripts, 'bin' and 'HMM'. If the GUI is in this "
                   "directory, it is found automatically.",
                   "包含EasyCGTree.pl等脚本、'bin'和'HMM'的目录。如果GUI就放在该目录中，会自动识别。"),
    "set_bin": ("bin directory", "bin目录"),
    "h_set_bin": ("Directory with the programs (hmmsearch, muscle5, ...). Empty: '<EasyCGTree directory>/bin'.",
                  "存放程序（hmmsearch、muscle5等）的目录。留空：'<EasyCGTree目录>/bin'。"),
    "set_perl": ("Perl program", "Perl程序"),
    "h_set_perl": ("Empty: found automatically ({p}).", "留空：自动查找（{p}）。"),
    "perl_not_found": ("not found", "未找到"),
    "set_programs": ("Programs (an empty path means: the program in the bin directory)",
                     "程序（路径留空表示使用bin目录中的程序）"),
    "col_program": ("Program", "程序"),
    "col_path": ("Custom path", "自定义路径"),
    "col_used": ("Used", "实际使用"),
    "reset": ("Reset", "恢复默认"),
    "dlg_check": ("Installation check", "安装检查"),
    "col_item": ("Item", "项目"),
    "col_status": ("Status", "状态"),
    "col_details": ("Details", "详情"),
    "checking": ("Checking...", "检查中..."),
    "copy_report": ("Copy report", "复制报告"),
    "dlg_treeopts": ("Tree program command lines", "建树程序命令行"),
    "treeopts_info": ("These command lines are stored in '{f}'. The first word is the program (in the bin directory, "
                      "or set in Paths and programs); the input file and the number of threads are added by EasyCGTree.",
                      "这些命令行保存在'{f}'中。第一个词为程序名（位于bin目录，或在'路径与程序'中设置）；输入文件和线程数由EasyCGTree自动添加。"),
    "treeopts_missing": ("'{f}' was not found.", "未找到'{f}'。"),
    "dlg_table": ("Table", "表格"),
    "dlg_about": ("About", "关于"),
    "about_text": ("Graphical interface for EasyCGTree {v}.", "EasyCGTree {v}的图形界面。"),
    # ---- version 1.1 ----
    "tk_predict": ("1 predict", "1 predict（预测）"),
    "tk_hmmsearch": ("2 hmmsearch", "2 hmmsearch（检索）"),
    "tk_refine": ("3 refine", "3 refine（筛选）"),
    "tk_alignment": ("4 alignment", "4 alignment（比对）"),
    "tk_tree_infer": ("5 tree_infer", "5 tree_infer（建树）"),
    "tasks_all": ("= all", "= all（全部）"),
    "outdir": ("Output directory", "输出文件夹"),
    "h_outdir": ("Folder for all results: the working directory '<input>_TEM', the trees, tables and log files (the layout "
                 "inside it is the same as before). Empty: the folder that contains the input directory ('-outdir').",
                 "所有结果的存放位置：工作目录'<输入文件夹>_TEM'、树文件、表格和日志文件（其中的相对层次不变）。留空：输入文件夹所在的上级文件夹（'-outdir'）。"),
    "h_hmm_outdir": ("Folder for the summary table, the log file and the working directory of BuildHMM. Default: the HMM "
                     "directory of EasyCGTree. The HMM set itself is always written to the HMM directory ('-outdir').",
                     "BuildHMM的汇总表、日志文件和工作目录的存放位置。默认：EasyCGTree的HMM文件夹。HMM集合本身总是写入HMM文件夹（'-outdir'）。"),
    "default_outdir": ("default: {d}", "默认：{d}"),
    "place_output": ("Output", "输出"),
    "place_input": ("Input", "输入"),
    "place_families": ("Gene families", "基因家族"),
    "notice_close": ("Close notice", "关闭提示"),
    "notice_pred": ("A CDS prediction for these input files with the same settings already exists ({info}). 'predict' was "
                    "unticked, so the existing prediction is reused. Close this notice to tick 'predict' again.",
                    "已存在用相同输入文件和设置得到的CDS预测结果（{info}）。已自动取消勾选'predict'，将直接使用已有预测。关闭此提示后才能重新勾选'predict'。"),
    "st_out": ("Output directory: {d}", "输出文件夹：{d}"),
    "st_no_dir": ("The folder does not exist: {d}", "文件夹不存在：{d}"),
    "in_dir": ("Input directory: {d}", "输入文件夹：{d}"),
    "in_files": ("Sequence files: {n}", "序列文件：{n}个"),
    "in_genomes": ("Genomes: {n} (from {f} files)", "基因组：{n}个（来自{f}个文件）"),
    "col_file": ("File", "文件"),
    "col_seqs": ("Sequences", "序列数"),
    "col_length": ("Total length", "总长度"),
    "col_gene": ("Gene family", "基因家族"),
    "col_mean_length": ("Mean length", "平均长度"),
    "type_short_nucl": ("DNA", "DNA"),
    "type_short_prot": ("protein", "蛋白"),
    "type_short_ncbi_pair": ("NCBI pair", "NCBI成对"),
    "type_short_ncbi_cds": ("NCBI CDS", "NCBI CDS"),
    "ignored": ("ignored", "忽略"),
    "fam_files": ("Gene family files: {n}", "基因家族文件：{n}个"),
    "ecg_from_record": ("Filled in from the record of the previous EasyCGTree run.", "已根据上一次EasyCGTree运行记录填入。"),
    "ecg_from_defaults": ("No previous EasyCGTree run found: default values.", "未找到之前的EasyCGTree运行记录：使用默认值。"),
    "quick_start": ("Quick start: just choose the input directory and click 'Run' – the default settings work for most analyses.",
                    "快速开始：选好输入文件夹后直接点击'运行'即可，默认设置适用于大多数分析。"),
    "strain_removed": ("removed by the screening", "已在筛选中被剔除"),
    "strains_removed_note": ("{r} genome(s) in grey were removed in the last screening (too few genes for the genome "
                             "cutoff) and are not in the alignments; they can't be chosen. A run with another genome "
                             "cutoff may change this.",
                             "灰色的{r}个基因组在上一次筛选中被剔除（基因数不满足基因组阈值），不在比对中，不能选择。若用其他基因组阈值重新运行，结果可能会变。"),
    "strains_none": ("Clear", "全部取消"),
    "strains_multi": ("{n} genomes; {k} ticked.", "共{n}个基因组；已勾选{k}个。"),
    "strains_keep": ("At least {m} genomes must stay unticked.", "至少要保留{m}个不勾选。"),
    "strains_single": ("{n} genomes; select one.", "共{n}个基因组；只能选择1个。"),
    "msg_ignore_too_many": ("Too many strains are ignored: at least two strains must be left for the SNP detection.",
                            "忽略的菌株太多：至少要保留2个菌株用于SNP检测。"),
    "msg_hmm_name_bad": ("The database name '{n}' may only contain letters, digits, '.', '_' and '-'.",
                         "数据库名称'{n}'只能包含字母、数字、'.'、'_'和'-'。"),
    "msg_hmm_exists": ("The HMM database '{n}' already exists. Please choose another name, or tick '-force' to overwrite it.",
                       "HMM数据库'{n}'已存在。请换一个名称，或勾选'-force'覆盖。"),
    "menu_hmm_online": ("HMM sets online", "在线下载HMM集"),
    "menu_hmm_folder": ("Open the HMM folder", "打开HMM文件夹"),
    "msg_hmm_missing": ("The HMM set '{n}' was not found in '{d}'.\n\nThe HMM sets are not included in the program package. "
                        "Please download '{n}.hmm' from the folder 'HMM' of the EasyCGTree5 repository and put it into this "
                        "folder (see HMM/README.txt; Help - HMM sets online):\n{gh}\n{ge}",
                        "在'{d}'中未找到HMM集'{n}'。\n\n程序包中不包含HMM集。请从EasyCGTree5仓库的'HMM'文件夹下载'{n}.hmm'，"
                        "放入该文件夹（见HMM/README.txt；帮助 - 在线下载HMM集）：\n{gh}\n{ge}"),
    "hmm_no_summary": ("No BuildHMM summary in this folder yet.", "该文件夹中还没有BuildHMM的汇总文件。"),
}

LANG = {"code": "en"}


def tr(key, **kw):
    """Translated text of a key (English if the key or the translation is missing)."""
    pair = TR.get(key)
    if pair is None:
        text = key
    else:
        text = pair[1] if LANG["code"] == "zh" and pair[1] else pair[0]
    if kw:
        try:
            text = text.format(**kw)
        except (KeyError, IndexError, ValueError):
            pass
    return text


# ======================================================================================================
# Helpers without Qt (paths, records, command lines)
# ======================================================================================================

SCRIPTS = {
    "main": "EasyCGTree.pl",
    "snp": "EasyCGTree_SNP.pl",
    "specific": "EasyCGTree_SpecificSNP.pl",
    "hmm": "BuildHMM.pl",
    "prevalence": "Gene_Prevelence.pl",
}
BASE_PROGRAMS = ["prodigal", "hmmsearch", "hmmbuild", "muscle5", "trimal"]
TASKS = ["predict", "hmmsearch", "refine", "alignment", "tree_infer"]
TREE_KEYS = [("iq", "IQ-TREE Command-line="), ("fast", "FastTree Command-line="), ("astral", "astral Command-line=")]


def env_key(program):
    """Name of the environment variable used by the scripts for a program (e.g. iqtree3 -> ECG_IQTREE3)."""
    name = re.sub(r"\.exe$", "", program, flags=re.I)
    return "ECG_" + re.sub(r"[^A-Z0-9_]", "_", name.upper())


def quote_arg(a):
    a = str(a)
    if IS_WIN:
        if a == "" or re.search(r'[\s"&|<>^()%!;,\'`$@{}]', a):
            return '"' + a.replace('"', '`"') + '"'
        return a
    return shlex.quote(a)


def read_record(tem):
    rec = {}
    p = Path(tem) / "EasyCGTree_record.txt"
    try:
        with open(p, encoding="utf-8", errors="replace") as fh:
            for line in fh:
                line = line.rstrip("\r\n")
                if line.startswith("#") or "=" not in line:
                    continue
                k, v = line.split("=", 1)
                rec[k] = v
    except OSError:
        pass
    return rec


def norm_path(d):
    p = Path(d).expanduser()
    try:
        return p.resolve()
    except OSError:
        return p.absolute()


def io_dirs(input_dir, outdir=""):
    """(input directory, output directory, working directory, input name). The output directory is
    '-outdir' if given, otherwise the directory that contains the input directory."""
    p = norm_path(input_dir)
    out = norm_path(outdir) if outdir else p.parent
    return p, out, out / (p.name + "_TEM"), p.name


def cpu_count():
    try:
        return len(os.sched_getaffinity(0))
    except (AttributeError, OSError):
        return os.cpu_count() or 1


def input_files(d):
    """Files of an input directory as EasyCGTree sees them (no hidden files, no sub-directories)."""
    try:
        return sorted(f.name for f in Path(d).iterdir() if f.is_file() and not f.name.startswith(".") and re.search(r"\w", f.name))
    except OSError:
        return []


def read_seq_types(tem):
    """Rows of Genome_SeqType.txt as dicts."""
    rows = []
    st = Path(tem) / "Genome_SeqType.txt"
    try:
        with open(st, encoding="utf-8", errors="replace") as fh:
            head = None
            for line in fh:
                f = line.rstrip("\r\n").split("\t")
                if head is None:
                    head = f
                    continue
                if f and f[0]:
                    rows.append(dict(zip(head, f)))
    except OSError:
        pass
    return rows


def prediction_match(input_dir, outdir, seq, keep):
    """(True, details) if the CDS prediction in the working directory fits the input files and the settings."""
    if not input_dir or not Path(input_dir).is_dir():
        return False, ""
    inp, out, tem, name = io_dirs(input_dir, outdir)
    rec = read_record(tem)
    rows = read_seq_types(tem)
    if not rec.get("predict_done") or not rows or not (tem / "TEM0_CDS").is_dir():
        return False, ""
    used = set()
    for r in rows:
        for col in ("Source_files", "Ignored_files"):
            used |= {x for x in r.get(col, "").split(",") if x and x != "-"}
        if "Source_files" not in r:
            used.add(r.get("Genome", "") + ".fas")
    files = set(input_files(inp))
    if files != used:
        return False, ""
    t0 = (tem / "Genome_SeqType.txt").stat().st_mtime
    if any((inp / f).stat().st_mtime > t0 for f in files):
        return False, ""
    if seq == "nucl" or keep:
        for r in rows:
            if not (tem / "TEM0_CDS" / (r.get("Genome", "") + ".fnn")).is_file():
                return False, ""
    return True, "%s; %d genomes" % (rec.get("predict_done", "")[:16], len(rows))


def fasta_stats(path):
    """(number of sequences, total length) of a FASTA file."""
    n, total = 0, 0
    try:
        with open(path, "rb") as fh:
            for line in fh:
                if line.startswith(b">"):
                    n += 1
                else:
                    total += len(line.strip())
    except OSError:
        pass
    return n, total


def list_files(directory, prefix, suffix):
    try:
        return sorted(f.name for f in Path(directory).iterdir() if f.is_file() and f.name.startswith(prefix) and f.name.endswith(suffix))
    except OSError:
        return []


def strains_of(input_dir, outdir=""):
    """Genome names known for an input directory (from the CDS prediction of EasyCGTree)."""
    if not input_dir:
        return []
    _, _, tem, _ = io_dirs(input_dir, outdir)
    names = []
    st = tem / "Genome_SeqType.txt"
    if st.is_file():
        with open(st, encoding="utf-8", errors="replace") as fh:
            for line in fh:
                f = line.rstrip("\r\n").split("\t")
                if f and f[0] and f[0] != "Genome":
                    names.append(f[0])
    elif (tem / "TEM0_CDS").is_dir():
        names = sorted(p.stem for p in (tem / "TEM0_CDS").glob("*.faa"))
    return names


def screened_out(input_dir, outdir=""):
    """Genomes removed by the last screening (EasyCGTree 'refine', '-genome_cutoff'), from GenomeGeneScreened.txt."""
    if not input_dir:
        return set()
    tem = io_dirs(input_dir, outdir)[2]
    f = tem / "GenomeGeneScreened.txt"
    if not (f.is_file() and read_record(tem).get("refine_done")):
        return set()
    kept = set()
    with open(f, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("Genome_List="):
                for x in line[len("Genome_List="):].split():
                    kept.add(re.sub(r"\.fas$", "", x.rsplit(":", 1)[0]))
    return {n for n in strains_of(input_dir, outdir) if n not in kept} if kept else set()


def newick_leaves(path):
    """Leaf names of a Newick tree (branch lengths, support values, comments and quotes are handled)."""
    try:
        with open(path, encoding="utf-8", errors="replace") as fh:
            text = fh.read()
    except OSError:
        return []
    text = re.sub(r"\[[^\]]*\]", "", text)
    names, i, n = [], 0, len(text)
    while i < n:
        if text[i] not in "(,":
            i += 1
            continue
        i += 1
        while i < n and text[i].isspace():
            i += 1
        if i < n and text[i] == "'":            # quoted label ('' is a quote)
            j, lab = i + 1, ""
            while j < n:
                if text[j] == "'":
                    if j + 1 < n and text[j + 1] == "'":
                        lab += "'"
                        j += 2
                        continue
                    break
                lab += text[j]
                j += 1
            names.append(lab)
            i = j + 1
        elif i < n and text[i] not in "(),:;":
            j = i
            while j < n and text[j] not in "(),:;":
                j += 1
            lab = text[i:j].strip()
            if lab:
                names.append(lab)
            i = j
    return names


def fasta_ids(path):
    ids = []
    try:
        with open(path, encoding="utf-8", errors="replace") as fh:
            for line in fh:
                if line.startswith(">"):
                    ids.append(line[1:].strip().split()[0] if line[1:].strip() else "")
    except OSError:
        pass
    return ids


def _num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def snp_plan(input_dir, outdir, ecg, snp):
    """What EasyCGTree_SNP.pl will have to do (mirrors its 'planEasyCGTree'), and whether the SNP extraction
    must be repeated because its settings differ from the last one.
    ecg: hmm, evalue, genome_cutoff, gene_cutoff, trim; snp: ignore, max_missing, aln, ref.
    Returns {'ecg_start': task or None, 'ecg_why': [...], 'snp_why': [...], 'snp_fas': path or None}."""
    _, out, tem, name = io_dirs(input_dir, outdir)
    rec = read_record(tem)
    start, why = None, []
    rows = read_seq_types(tem)
    if not tem.is_dir() or not rows or not rec.get("predict_done"):
        start = "predict"
        why.append(tr("why_no_predict"))
    elif any(not (tem / "TEM0_CDS" / (r.get("Genome", "") + ".fnn")).is_file() for r in rows if r.get("Input_type") != "prot"):
        start = "predict"
        why.append(tr("why_no_fnn"))
    hmm = re.sub(r"\.hmm$", "", os.path.basename(str(ecg.get("hmm") or "")), flags=re.I)
    if start is None:
        if not rec.get("hmmsearch_done"):
            start = "hmmsearch"; why.append(tr("why_no_step", s="hmmsearch"))
        elif hmm and hmm != rec.get("hmm", ""):
            start = "hmmsearch"; why.append(tr("why_diff", k="hmm", v=hmm, o=rec.get("hmm", "?")))
        elif _num(ecg.get("evalue")) is not None and _num(ecg.get("evalue")) > (_num(rec.get("evalue")) or 0):
            start = "hmmsearch"; why.append(tr("why_diff", k="evalue", v=ecg.get("evalue"), o=rec.get("evalue", "?")))
    if start is None:
        if not rec.get("refine_done"):
            start = "refine"; why.append(tr("why_no_step", s="refine"))
        else:
            for k in ("genome_cutoff", "gene_cutoff"):
                if _num(ecg.get(k)) is not None and _num(ecg.get(k)) != _num(rec.get(k)):
                    start = "refine"; why.append(tr("why_diff", k=k, v=ecg.get(k), o=rec.get(k, "NA")))
            ev = _num(ecg.get("evalue"))
            if ev is not None and ev != _num(rec.get("refine_evalue", rec.get("evalue"))):
                start = "refine"; why.append(tr("why_diff", k="evalue", v=ecg.get("evalue"), o=rec.get("refine_evalue", rec.get("evalue", "?"))))
    if start is None:
        if not rec.get("alignment_done") or rec.get("alignment_seq", "") != "nucl":
            start = "alignment"; why.append(tr("why_no_nucl"))
        elif ecg.get("trim") and ecg.get("trim") != rec.get("trim", ""):
            start = "alignment"; why.append(tr("why_diff", k="trim", v=ecg.get("trim"), o=rec.get("trim", "NA")))
        elif not (tem / "TEM5_Alignment").is_dir() or not (tem / "TEM6_AlnTrimmed").is_dir():
            start = "alignment"; why.append(tr("why_no_step", s="alignment"))
    # the SNP extraction
    swhy, fas = [], None
    if rec.get("snp_done") and rec.get("snp_hmm") and rec.get("snp_genes"):
        f = out / ("%s.%s.%s.snp.fas" % (name, rec["snp_hmm"], rec["snp_genes"]))
        if f.is_file() and f.with_name(f.name[:-len("snp.fas")] + "snp_positions.txt").is_file():
            fas = f
    if start is not None:
        swhy.append(tr("why_snp_ecg"))
    elif fas is None:
        swhy.append(tr("why_snp_none"))
    else:
        if rec.get("alignment_done", "") > rec.get("snp_done", ""):
            swhy.append(tr("why_snp_old"))
        norm = lambda v: ",".join(sorted((x.strip() for x in str(v or "").split(",") if x.strip()), key=str.lower)).lower()
        if norm(snp.get("ignore")) != norm(rec.get("snp_ignore")):
            swhy.append(tr("why_diff", k="ignore", v=snp.get("ignore") or "-", o=rec.get("snp_ignore") or "-"))
        if _num(snp.get("max_missing")) != _num(rec.get("snp_max_missing")):
            swhy.append(tr("why_diff", k="max_missing", v=snp.get("max_missing"), o=rec.get("snp_max_missing", "NA")))
        if (snp.get("aln") or "trimmed") != rec.get("snp_aln", ""):
            swhy.append(tr("why_diff", k="aln", v=snp.get("aln"), o=rec.get("snp_aln", "NA")))
        if (snp.get("ref") or "").lower() != (rec.get("snp_ref") or "").lower():
            swhy.append(tr("why_diff", k="ref", v=snp.get("ref") or "-", o=rec.get("snp_ref") or "-"))
        if hmm and hmm != rec.get("snp_hmm"):
            swhy.append(tr("why_diff", k="hmm", v=hmm, o=rec.get("snp_hmm")))
    return {"ecg_start": start, "ecg_why": why, "snp_why": swhy, "snp_fas": fas if not swhy else None}


def dry_run_names(input_dir, perl, script, env):
    """Genome names from a dry run of EasyCGTree (before any CDS prediction)."""
    if not (perl and input_dir and Path(script).is_file()):
        return []
    rc, out = capture([perl, str(script), "-input", str(input_dir), "-dry_run"], env=env, timeout=300)
    m = re.search(r"^EasyCGTree_DRY_RUN_JSON=(.*)$", out, re.M)
    if rc != 0 or not m:
        return []
    try:
        return [g["name"] for g in json.loads(m.group(1)).get("genomes", [])]
    except (ValueError, KeyError, TypeError):
        return []


def safe_hmm_name(name):
    """A database name that BuildHMM accepts (letters, digits, '.', '_', '-')."""
    name = re.sub(r"\.hmm$", "", name or "", flags=re.I)
    return re.sub(r"[^\w.\-]+", "_", name, flags=re.A).strip("_") or "myHMM"


def hmm_sets(hmm_dir):
    """[(name, number of HMMs, format line)] of the .hmm files in a directory."""
    out = []
    try:
        files = sorted(Path(hmm_dir).glob("*.hmm"))
    except OSError:
        return out
    for f in files:
        n, fmt = 0, ""
        try:
            with open(f, encoding="utf-8", errors="replace") as fh:
                for i, line in enumerate(fh):
                    if i == 0:
                        fmt = line.strip()
                    if line.startswith("NAME "):
                        n += 1
        except OSError:
            continue
        out.append((f.stem, n, fmt))
    return out


def hmm_format_version(fmt_line):
    """'HMMER3/f [3.4 | Aug 2023]' -> 'f'."""
    m = re.match(r"HMMER3/(\w)", fmt_line or "")
    return m.group(1) if m else ""


def read_tree_options(path):
    opts = {}
    try:
        with open(path, encoding="utf-8", errors="replace") as fh:
            for line in fh:
                line = line.rstrip("\r\n")
                for key, prefix in TREE_KEYS:
                    if line.strip().startswith(prefix):
                        opts[key] = line.strip()[len(prefix):]
    except OSError:
        pass
    return opts


def write_tree_options(path, values):
    """Replace the three command lines in tree_app-options.txt, keeping everything else."""
    with open(path, encoding="utf-8", errors="replace") as fh:
        lines = fh.read().splitlines()
    bak = str(path) + ".bak"
    if not os.path.exists(bak):
        shutil.copyfile(path, bak)
    done = set()
    for i, line in enumerate(lines):
        for key, prefix in TREE_KEYS:
            if line.strip().startswith(prefix) and key in values:
                lines[i] = prefix + values[key].strip()
                done.add(key)
    for key, prefix in TREE_KEYS:
        if key in values and key not in done:
            lines.append(prefix + values[key].strip())
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("\n".join(lines) + "\n")


def find_perl(configured=""):
    """Perl interpreter: the configured one, else from PATH or common places."""
    if configured:
        return configured if os.path.isfile(configured) else ""
    p = shutil.which("perl")
    if p:
        return p
    cands = []
    if IS_WIN:
        cands = [r"C:\Strawberry\perl\bin\perl.exe", r"C:\Perl64\bin\perl.exe", r"C:\Perl\bin\perl.exe"]
    else:
        cands = ["/usr/bin/perl", "/usr/local/bin/perl", "/opt/homebrew/bin/perl"]
    for c in cands:
        if os.path.isfile(c):
            return c
    return ""


def no_window_flags():
    return 0x08000000 if IS_WIN else 0  # CREATE_NO_WINDOW


def capture(argv, env=None, timeout=60, cwd=None):
    """Run a short command, return (exit code, output text). Exit code -1: could not be started; -2: timeout."""
    try:
        r = subprocess.run(argv, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, stdin=subprocess.DEVNULL,
                           env=env, cwd=cwd, timeout=timeout, creationflags=no_window_flags())
        return r.returncode, r.stdout.decode("utf-8", errors="replace")
    except subprocess.TimeoutExpired as e:
        out = e.output.decode("utf-8", errors="replace") if e.output else ""
        return -2, out
    except OSError as e:
        return -1, str(e)


def app_dir():
    """Directory of the GUI: of the script, or of the program made with PyInstaller (for a macOS app bundle,
    the directory that contains 'EasyCGTree_GUI.app')."""
    if getattr(sys, "frozen", False):
        exe = Path(sys.executable).resolve()
        for p in exe.parents:
            if p.suffix == ".app":
                return p.parent
        return exe.parent
    return Path(__file__).resolve().parent


class Config:
    """Paths and program locations (stored with QSettings as one JSON text)."""

    def __init__(self):
        self.qs = QSettings("EasyCGTree", "EasyCGTree_GUI")
        try:
            data = json.loads(self.qs.value("config", "{}") or "{}")
        except (ValueError, TypeError):
            data = {}
        self.home = data.get("home", "")
        self.bin = data.get("bin", "")
        self.perl = data.get("perl", "")
        self.overrides = dict(data.get("overrides", {}))
        self.language = data.get("language", "en")
        self.full_cmd = bool(data.get("full_cmd", False))
        self.params = data.get("params", {})
        self.geometry = data.get("geometry", "")
        self.splitters = data.get("splitters", {}) if isinstance(data.get("splitters"), dict) else {}
        self.right_page = data.get("right_page", 0)
        here = app_dir()
        if (here / SCRIPTS["main"]).is_file() and (not self.home or not (Path(self.home) / SCRIPTS["main"]).is_file()):
            self.home = str(here)

    def save(self):
        data = {"home": self.home, "bin": self.bin, "perl": self.perl, "overrides": self.overrides,
                "language": self.language, "full_cmd": self.full_cmd, "params": self.params, "geometry": self.geometry,
                "splitters": self.splitters, "right_page": self.right_page}
        self.qs.setValue("config", json.dumps(data))
        self.qs.sync()

    def home_ok(self):
        return bool(self.home) and (Path(self.home) / SCRIPTS["main"]).is_file()

    def bin_dir(self):
        return Path(self.bin) if self.bin else Path(self.home or ".") / "bin"

    def hmm_dir(self):
        return Path(self.home or ".") / "HMM"

    def script(self, key):
        return Path(self.home or ".") / SCRIPTS[key]

    def tree_options_file(self):
        f = self.bin_dir() / "tree_app-options.txt"
        return f if f.is_file() else Path(self.home or ".") / "bin" / "tree_app-options.txt"

    def tree_programs(self):
        opts = read_tree_options(self.tree_options_file())
        progs = []
        for key, _ in TREE_KEYS:
            words = opts.get(key, "").split()
            if words and words[0] not in progs:
                progs.append(words[0])
        return progs or ["iqtree3", "FastTreeMP", "astral-weighted"]

    def programs(self):
        return BASE_PROGRAMS + [p for p in self.tree_programs() if p not in BASE_PROGRAMS]

    def program_path(self, name):
        o = self.overrides.get(name, "")
        return Path(o) if o else self.bin_dir() / (re.sub(r"\.exe$", "", name, flags=re.I) + EXE)

    def perl_path(self):
        return find_perl(self.perl)

    def env_overrides(self):
        """Environment variables for the scripts: ECG_BIN and ECG_<PROGRAM>."""
        env = {}
        if self.bin and Path(self.bin).resolve() != (Path(self.home or ".") / "bin").resolve():
            env["ECG_BIN"] = str(Path(self.bin))
        for name, path in sorted(self.overrides.items()):
            if path:
                env[env_key(name)] = path
        return env

    def run_env(self):
        env = dict(os.environ)
        env.update(self.env_overrides())
        env["PYTHONIOENCODING"] = "utf-8"
        return env


def format_commands(cmds, config):
    """Command lines for the terminal (POSIX shells, or PowerShell on Windows)."""
    lines = []
    env = config.env_overrides()
    perl = config.perl_path()
    perl_disp = "perl" if (not config.perl and perl and shutil.which("perl") == perl) or not perl else perl
    if IS_WIN and env:
        for k, v in env.items():
            lines.append('$env:%s = "%s"' % (k, v))
    prefix = "" if IS_WIN else "".join("%s=%s " % (k, shlex.quote(v)) for k, v in env.items())
    for argv in cmds:
        words = [quote_arg(perl_disp)] + [quote_arg(a) for a in argv]
        line = " ".join(words)
        if IS_WIN and words[0].startswith('"'):
            line = "& " + line
        lines.append(prefix + line)
    return "\n".join(lines)


# ======================================================================================================
# Running the scripts
# ======================================================================================================

class Runner(QObject):
    """Runs a list of commands one after the other in a background thread."""
    output = Signal(str)
    started_cmd = Signal(str)
    finished = Signal(int, bool)  # exit code, stopped by the user

    def __init__(self):
        super().__init__()
        self.proc = None
        self.thread = None
        self._stop = False

    def running(self):
        return self.thread is not None and self.thread.is_alive()

    def start(self, cmds, env, cwd):
        self._stop = False
        self.thread = threading.Thread(target=self._work, args=(cmds, env, cwd), daemon=True)
        self.thread.start()

    def _work(self, cmds, env, cwd):
        rc = 0
        for argv in cmds:
            if self._stop:
                break
            self.started_cmd.emit(" ".join(quote_arg(a) for a in argv))
            try:
                kw = {}
                if IS_WIN:
                    kw["creationflags"] = 0x00000200 | no_window_flags()  # CREATE_NEW_PROCESS_GROUP
                else:
                    kw["start_new_session"] = True
                self.proc = subprocess.Popen(argv, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                             stdin=subprocess.DEVNULL, env=env, cwd=cwd, **kw)
            except OSError as e:
                self.output.emit("\nERROR: can't start '%s': %s\n" % (argv[0], e))
                rc = 255
                break
            dec = codecs.getincrementaldecoder("utf-8")(errors="replace")
            while True:
                chunk = self.proc.stdout.read1(4096) if hasattr(self.proc.stdout, "read1") else self.proc.stdout.read(1)
                if not chunk:
                    break
                self.output.emit(dec.decode(chunk))
            rest = dec.decode(b"", final=True)
            if rest:
                self.output.emit(rest)
            rc = self.proc.wait()
            self.proc = None
            if rc != 0:
                break
        self.finished.emit(rc, self._stop)

    def stop(self):
        self._stop = True
        p = self.proc
        if p is None or p.poll() is not None:
            return
        try:
            if IS_WIN:
                subprocess.run(["taskkill", "/F", "/T", "/PID", str(p.pid)], stdout=subprocess.DEVNULL,
                               stderr=subprocess.DEVNULL, creationflags=no_window_flags())
            else:
                os.killpg(p.pid, signal.SIGTERM)
                threading.Timer(3.0, self._kill_hard, args=(p,)).start()
        except (OSError, ProcessLookupError):
            pass

    @staticmethod
    def _kill_hard(p):
        if p.poll() is None:
            try:
                os.killpg(p.pid, signal.SIGKILL)
            except (OSError, ProcessLookupError):
                pass


class Capture(QObject):
    """Runs one short command in the background and calls back with (exit code, output)."""
    done = Signal(int, str)

    def start(self, argv, env, cwd=None, timeout=600):
        threading.Thread(target=lambda: self.done.emit(*capture(argv, env, timeout, cwd)), daemon=True).start()


# ======================================================================================================
# Small widgets
# ======================================================================================================

def fixed_font():
    return QFontDatabase.systemFont(QFontDatabase.SystemFont.FixedFont)


class I18n:
    """Remembers (setter, key) pairs so that all texts can be changed when the language changes."""

    def __init__(self):
        self.items = []

    def add(self, setter, key, **kw):
        self.items.append((setter, key, kw))
        setter(tr(key, **kw))

    def apply(self):
        for setter, key, kw in self.items:
            try:
                setter(tr(key, **kw))
            except RuntimeError:  # widget already deleted
                pass


class _WheelGuard:
    """The mouse wheel changes the value only after the field has been clicked (or got the focus), and
    only while the pointer is on it; otherwise the wheel scrolls the parameter area."""

    def _guard_init(self):
        self._wheel_ok = False
        self.setFocusPolicy(Qt.FocusPolicy.StrongFocus)

    def mousePressEvent(self, e):
        self._wheel_ok = True
        super().mousePressEvent(e)

    def focusInEvent(self, e):
        self._wheel_ok = True
        super().focusInEvent(e)

    def focusOutEvent(self, e):
        self._wheel_ok = False
        super().focusOutEvent(e)

    def leaveEvent(self, e):
        self._wheel_ok = False
        super().leaveEvent(e)

    def wheelEvent(self, e):
        if self._wheel_ok:
            super().wheelEvent(e)
        else:
            e.ignore()


class GComboBox(_WheelGuard, QComboBox):
    def __init__(self, *a):
        super().__init__(*a)
        self._guard_init()


class GSpinBox(_WheelGuard, QSpinBox):
    def __init__(self, *a):
        super().__init__(*a)
        self._guard_init()


class GDoubleSpinBox(_WheelGuard, QDoubleSpinBox):
    def __init__(self, *a):
        super().__init__(*a)
        self._guard_init()


class FlowLayout(QLayout):
    """Places the widgets in a row and wraps them into further rows when the space is too narrow."""

    def __init__(self, parent=None, hspacing=10, vspacing=4):
        super().__init__(parent)
        self._items = []
        self._h = hspacing
        self._v = vspacing
        self.setContentsMargins(0, 0, 0, 0)

    def addItem(self, item):
        self._items.append(item)

    def count(self):
        return len(self._items)

    def itemAt(self, i):
        return self._items[i] if 0 <= i < len(self._items) else None

    def takeAt(self, i):
        return self._items.pop(i) if 0 <= i < len(self._items) else None

    def expandingDirections(self):
        return Qt.Orientation(0)

    def hasHeightForWidth(self):
        return True

    def heightForWidth(self, width):
        return self._arrange(QRect(0, 0, width, 0), True)

    def setGeometry(self, rect):
        super().setGeometry(rect)
        self._arrange(rect, False)

    def sizeHint(self):
        return self.minimumSize()

    def minimumSize(self):
        size = QSize()
        for it in self._items:
            size = size.expandedTo(it.minimumSize())
        m = self.contentsMargins()
        return size + QSize(m.left() + m.right(), m.top() + m.bottom())

    def _arrange(self, rect, test_only):
        m = self.contentsMargins()
        r = rect.adjusted(m.left(), m.top(), -m.right(), -m.bottom())
        x, y, line = r.x(), r.y(), 0
        for it in self._items:
            if it.isEmpty():
                continue
            hint = it.sizeHint()
            if x > r.x() and x + hint.width() > r.right() + 1:
                x = r.x()
                y += line + self._v
                line = 0
            if not test_only:
                it.setGeometry(QRect(QPoint(x, y), hint))
            x += hint.width() + self._h
            line = max(line, hint.height())
        return y + line - rect.y() + m.bottom()


def narrow_combo(combo, chars=10):
    """A combo box that may become narrow; its list still shows the whole texts."""
    combo.setSizeAdjustPolicy(QComboBox.SizeAdjustPolicy.AdjustToMinimumContentsLengthWithIcon)
    combo.setMinimumContentsLength(chars)
    fit_popup(combo)


def fit_popup(combo):
    try:
        combo.view().setMinimumWidth(combo.view().sizeHintForColumn(0) + 30)
    except (TypeError, AttributeError, RuntimeError):
        pass


class HelpButton(QToolButton):
    def __init__(self, key, parent=None, **kw):
        super().__init__(parent)
        self.key = key
        self.kw = kw
        self.setText("?")
        self.setObjectName("helpButton")
        self.setFixedSize(18, 18)
        self.setCursor(Qt.CursorShape.PointingHandCursor)
        self.clicked.connect(self.show_help)

    def show_help(self):
        QMessageBox.information(self, tr("help_title"), tr(self.key, **self.kw))


class LogView(QPlainTextEdit):
    """Read-only log that handles carriage returns (progress lines) like a terminal."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setReadOnly(True)
        self.setFont(fixed_font())
        self.setLineWrapMode(QPlainTextEdit.LineWrapMode.NoWrap)
        self.setMaximumBlockCount(50000)
        self._line = ""
        self._cr = False
        self._blank = True

    def reset(self):
        self.clear()
        self._line = ""
        self._cr = False
        self._blank = True

    def _set_last_line(self, text):
        c = self.textCursor()
        c.movePosition(QTextCursor.MoveOperation.End)
        c.movePosition(QTextCursor.MoveOperation.StartOfBlock, QTextCursor.MoveMode.KeepAnchor)
        c.removeSelectedText()
        c.insertText(text)

    def feed(self, text):
        for piece in re.split(r"(\r\n|\r|\n)", text):
            if piece in ("\n", "\r\n"):
                # At most one empty line in a row (the scripts print many blank lines).
                if self._line.strip() == "" and self._blank:
                    self._line = ""
                    self._cr = False
                    continue
                self._blank = self._line.strip() == ""
                c = self.textCursor()
                c.movePosition(QTextCursor.MoveOperation.End)
                c.insertText("\n")
                self._line = ""
                self._cr = False
            elif piece == "\r":
                self._cr = True
            elif piece:
                if self._cr:
                    self._line = piece
                    self._cr = False
                else:
                    self._line += piece
                self._set_last_line(self._line)
        sb = self.verticalScrollBar()
        sb.setValue(sb.maximum())


class PathEdit(QWidget):
    """Line edit with a 'Browse...' button for a directory or a file."""

    def __init__(self, i18n, kind="dir", parent=None, file_filter=""):
        super().__init__(parent)
        self.kind = kind
        self.file_filter = file_filter
        lay = QHBoxLayout(self)
        lay.setContentsMargins(0, 0, 0, 0)
        self.edit = QLineEdit()
        self.edit.setMinimumWidth(60)
        self.button = QPushButton()
        i18n.add(self.button.setText, "browse")
        lay.addWidget(self.edit, 1)
        lay.addWidget(self.button)
        self.button.clicked.connect(self.browse)

    def browse(self):
        start = self.edit.text() or str(Path.home())
        if self.kind == "dir":
            p = QFileDialog.getExistingDirectory(self, tr("browse"), start)
        else:
            p, _ = QFileDialog.getOpenFileName(self, tr("browse"), start, self.file_filter)
        if p:
            self.edit.setText(os.path.normpath(p))

    def text(self):
        return self.edit.text().strip()

    def setText(self, t):
        self.edit.setText(t or "")


# ======================================================================================================
# Parameters
# ======================================================================================================

class Param:
    """One option of a script: widget, value, and its part of the command line."""

    def __init__(self, tab, key, flag, kind, default=None, label=None, help=None, choices=None, minimum=None,
                 maximum=None, step=None, decimals=2, optional=False, required=False, multi=True, file_filter=""):
        self.tab = tab
        self.key = key
        self.flag = flag
        self.kind = kind
        self.default = default
        self.label_key = label or key
        self.help_key = help or ("h_" + key)
        self.choices = choices or []
        self.minimum = minimum
        self.maximum = maximum
        self.step = step
        self.decimals = decimals
        self.optional = optional
        self.required = required
        self.multi = multi
        self.keep_free = 0          # strains: how many must stay unselected
        self.file_filter = file_filter
        self.active = True
        self.label = None
        self.widget = None
        self.field = None
        self.helpbtn = None

    # ---- widgets ----
    def build(self):
        i18n = self.tab.i18n
        self.label = QLabel()
        i18n.add(self.label.setText, self.label_key)
        k = self.kind
        if k in ("dir", "file"):
            self.field = PathEdit(i18n, k, file_filter=self.file_filter)
            self.field.edit.textChanged.connect(self.tab.changed)
            self.widget = self.field
        elif k == "int":
            self.field = GSpinBox()
            self.field.setMinimumWidth(100)
            self.field.setRange(self.minimum if self.minimum is not None else 0, self.maximum if self.maximum is not None else 1000000)
            self.field.setValue(int(self.default))
            self.field.valueChanged.connect(self.tab.changed)
            self.widget = self.field
        elif k == "float":
            self.field = GDoubleSpinBox()
            self.field.setMinimumWidth(100)
            self.field.setDecimals(self.decimals)
            self.field.setRange(self.minimum if self.minimum is not None else 0.0, self.maximum if self.maximum is not None else 1e9)
            self.field.setSingleStep(self.step or 0.05)
            self.field.setValue(float(self.default))
            self.field.valueChanged.connect(self.tab.changed)
            self.widget = self.field
        elif k in ("efloat", "text"):
            self.field = QLineEdit()
            self.field.setMinimumWidth(60)
            if k == "efloat":
                self.field.setFixedWidth(120)
            if self.default is not None:
                self.field.setText(str(self.default))
            if self.optional:
                i18n.add(self.field.setPlaceholderText, "placeholder_previous")
            self.field.textChanged.connect(self.tab.changed)
            self.widget = self.field
        elif k == "choice":
            self.field = GComboBox()
            if self.optional:
                self.field.addItem(tr("previous"), None)
                i18n.add(lambda t, f=self.field: f.setItemText(0, t), "previous")
            for value, lkey in self.choices:
                idx = self.field.count()
                self.field.addItem(tr(lkey) if lkey else value, value)
                if lkey:
                    i18n.add(lambda t, f=self.field, i=idx: f.setItemText(i, t), lkey)
            # as wide as its longest choice (recomputed when the language changes)
            self.field.setSizeAdjustPolicy(QComboBox.SizeAdjustPolicy.AdjustToContents)
            fit_popup(self.field)
            self.set_value(self.default)
            self.field.currentIndexChanged.connect(self.tab.changed)
            self.widget = self.field
        elif k == "bool":
            # a check box carries its own text and gets a row of its own (not in the label column)
            self.field = QCheckBox()
            i18n.add(self.field.setText, self.label_key)
            self.field.setChecked(bool(self.default))
            self.field.toggled.connect(self.tab.changed)
            self.widget = self.field
        elif k == "hmm":
            self.widget = QWidget()
            lay = QHBoxLayout(self.widget)
            lay.setContentsMargins(0, 0, 0, 0)
            self.field = GComboBox()
            self.field.setEditable(True)
            narrow_combo(self.field)
            if self.optional:
                i18n.add(self.field.lineEdit().setPlaceholderText, "placeholder_previous")
            b = QPushButton()
            i18n.add(b.setText, "browse")
            b.clicked.connect(self.browse_hmm)
            lay.addWidget(self.field, 1)
            lay.addWidget(b)
            self.reload_hmm()
            self.field.editTextChanged.connect(self.tab.changed)
        elif k == "tasks":
            self.widget = QWidget()
            lay = FlowLayout(self.widget)
            self.boxes = []
            self._tsync = False
            for i, t in enumerate(TASKS):
                b = QCheckBox()
                i18n.add(b.setText, "tk_" + t)
                i18n.add(b.setToolTip, "task_" + t)
                b.setChecked(True)
                b.toggled.connect(lambda on, i=i: self.task_toggled(i, on))
                self.boxes.append(b)
                lay.addWidget(b)
            self.all_label = QLabel()
            self.all_label.setStyleSheet("color: #2a6f2a; font-weight: bold")
            sp = self.all_label.sizePolicy()
            sp.setRetainSizeWhenHidden(True)        # the '?' after it does not move
            self.all_label.setSizePolicy(sp)
            i18n.add(self.all_label.setText, "tasks_all")
            lay.addWidget(self.all_label)
            self.field = self.widget
            self.update_all_label()
        elif k == "strains":
            self.widget = QWidget()
            lay = QHBoxLayout(self.widget)
            lay.setContentsMargins(0, 0, 0, 0)
            self.field = QLineEdit()
            b = QPushButton()
            i18n.add(b.setText, "choose")
            b.clicked.connect(self.choose_strains)
            lay.addWidget(self.field, 1)
            lay.addWidget(b)
            self.field.textChanged.connect(self.tab.changed)
        self.helpbtn = HelpButton(self.help_key)
        return self.label, self.widget, self.helpbtn

    def reload_hmm(self):
        cur = self.field.currentText() if self.field.count() or self.field.currentText() else (self.default or "")
        self.field.blockSignals(True)
        self.field.clear()
        if self.optional:
            self.field.addItem("")
        for name, n, fmt in hmm_sets(self.tab.win.config.hmm_dir()):
            self.field.addItem(name)
        self.field.setEditText(cur)
        self.field.blockSignals(False)

    def browse_hmm(self):
        p, _ = QFileDialog.getOpenFileName(self.widget, tr("browse"), str(self.tab.win.config.hmm_dir()), "HMM (*.hmm);;* (*)")
        if p:
            self.field.setEditText(os.path.normpath(p))

    def choose_strains(self):
        tab = self.tab
        outdir = tab.v("outdir") if "outdir" in tab.params else ""
        own = tab.strain_names(self.key)
        if own is not None:            # e.g. the leaves of a given tree
            if not own:
                QMessageBox.information(self.widget, tr("info"), tr("msg_no_tree_names"))
                return
            current = [s.strip() for s in self.field.text().split(",") if s.strip()]
            dlg = StrainDialog(own, current, self.multi, self.widget, keep=self.keep_free)
            if dlg.exec():
                self.field.setText(",".join(dlg.selected()))
            return
        names = strains_of(tab.input_dir(), outdir)
        if not names and tab.input_dir():
            cfg = tab.win.config
            names = dry_run_names(tab.input_dir(), cfg.perl_path(), cfg.script("main"), cfg.run_env())
        if not names:
            QMessageBox.information(self.widget, tr("info"), tr("msg_no_strains"))
            return
        current = [s.strip() for s in self.field.text().split(",") if s.strip()]
        removed = screened_out(tab.input_dir(), outdir)
        dlg = StrainDialog(names, current, self.multi, self.widget, keep=self.keep_free, removed=removed)
        if dlg.exec():
            self.field.setText(",".join(dlg.selected()))

    def task_toggled(self, i, on):
        """Keep the chosen tasks consecutive: ticking two tasks also ticks the ones between them;
        unticking a task in the middle also unticks the later ones; at least one task stays ticked."""
        if self._tsync:
            return
        self._tsync = True
        sel = [k for k, b in enumerate(self.boxes) if b.isChecked()]
        if on:
            for k in range(min(sel), max(sel) + 1):
                if self.boxes[k].isEnabled():
                    self.boxes[k].setChecked(True)
        elif not sel:
            self.boxes[i].setChecked(True)
        elif any(k < i for k in sel) and any(k > i for k in sel):
            for k in sel:
                if k > i:
                    self.boxes[k].setChecked(False)
        self._tsync = False
        self.update_all_label()
        self.tab.changed()

    def update_all_label(self):
        if self.kind == "tasks":
            on = all(b.isChecked() for b in self.boxes)
            self.all_label.setVisible(on)

    def set_active(self, on):
        self.active = bool(on)
        for w in (self.label, self.widget):
            if w is not None:
                w.setEnabled(self.active)

    # ---- values ----
    def value(self):
        k = self.kind
        if k in ("dir", "file"):
            t = self.field.text()
            return os.path.abspath(os.path.expanduser(t)) if t else ""
        if k == "int":
            return int(self.field.value())
        if k == "float":
            return round(float(self.field.value()), self.decimals)
        if k in ("efloat", "text", "strains"):
            return self.field.text().strip()
        if k == "choice":
            return self.field.currentData()
        if k == "bool":
            return bool(self.field.isChecked())
        if k == "hmm":
            return self.field.currentText().strip()
        if k == "tasks":
            return [t for t, b in zip(TASKS, self.boxes) if b.isChecked()]
        return None

    def set_value(self, v):
        k = self.kind
        if k in ("dir", "file", "efloat", "text", "strains"):
            self.field.setText("" if v is None else str(v))
        elif k == "int":
            try:
                self.field.setValue(int(v))
            except (TypeError, ValueError):
                pass
        elif k == "float":
            try:
                self.field.setValue(float(v))
            except (TypeError, ValueError):
                pass
        elif k == "choice":
            for i in range(self.field.count()):
                if self.field.itemData(i) == v:
                    self.field.setCurrentIndex(i)
                    return
        elif k == "bool":
            self.field.setChecked(bool(v))
        elif k == "hmm":
            self.field.setEditText("" if v is None else str(v))
        elif k == "tasks":
            if isinstance(v, str):
                v = TASKS if v == "all" else [x for x in v.split(",") if x]
            v = [x for x in (v or TASKS) if x in TASKS] or TASKS
            self._tsync = True
            for t, b in zip(TASKS, self.boxes):
                b.setChecked(t in v)
            self._tsync = False
            self.update_all_label()
            self.tab.changed()   # like a click on the boxes (ignored while the tab is being built or restored)

    def is_default(self, v):
        if self.kind == "tasks":
            return list(v) == TASKS
        d = self.default
        if d is None:
            return False
        if self.kind in ("float", "efloat"):
            try:
                return abs(float(v) - float(d)) <= 1e-12 * max(1.0, abs(float(d)))
            except (TypeError, ValueError):
                return False
        return v == d

    def fmt(self, v):
        if self.kind == "tasks":
            return "all" if list(v) == TASKS else ",".join(v)
        if self.kind == "float":
            return ("%.*f" % (self.decimals, v)).rstrip("0").rstrip(".") or "0"
        return str(v)

    def args(self, full):
        """Command-line words of this option."""
        if not self.active or self.flag is None:
            return []
        v = self.value()
        if v is None or v == "" or v == []:
            return []
        if self.kind == "bool":
            return ["-" + self.flag] if v else []
        if not full and self.is_default(v):
            return []
        return ["-" + self.flag, self.fmt(v)]

    def validate(self):
        """Error message, or None."""
        if not self.active:
            return None
        v = self.value()
        label = tr(self.label_key)
        if self.required and not v:
            return tr("msg_required", label=label)
        if self.kind == "hmm" and v and not os.path.isfile(v):
            name = re.sub(r"\.hmm$", "", os.path.basename(v), flags=re.I)
            hdir = self.tab.win.config.hmm_dir()
            if not (hdir / (name + ".hmm")).is_file():
                return tr("msg_hmm_missing", n=name, d=str(hdir), gh=REPO_HMM["GitHub"], ge=REPO_HMM["Gitee"])
        if self.kind == "efloat" and v:
            try:
                x = float(v)
            except ValueError:
                return tr("msg_bad_number", label=label, v=v)
            if self.minimum is not None and self.maximum is not None and not (self.minimum <= x <= self.maximum):
                return tr("msg_range", label=label, a=self.minimum, b=self.maximum)
            if x < 0:
                return tr("msg_bad_number", label=label, v=v)
        return None


# ======================================================================================================
# Dialogs
# ======================================================================================================

class StrainDialog(QDialog):
    """List of the genomes; several can be ticked (multi), or one can be selected.
    'keep': at least this many genomes must stay unticked."""

    def __init__(self, names, current, multi, parent=None, keep=0, removed=()):
        super().__init__(parent)
        self.removed = set(removed)
        self.setWindowTitle(tr("dlg_strains"))
        self.multi = multi
        self.keep = keep
        self._sync = False
        lay = QVBoxLayout(self)
        self.info = QLabel()
        self.info.setWordWrap(True)
        lay.addWidget(self.info)
        self.list = QListWidget()
        cur = set(current)
        self.names = []
        for n in names:
            it = QListWidgetItem(n)
            if n in self.removed:
                # removed by the screening: not in the alignments, so it can't be chosen
                it.setText(n + "   (" + tr("strain_removed") + ")")
                it.setFlags(Qt.ItemFlag.NoItemFlags)
                it.setForeground(QBrush(QColor("#999999")))
            elif multi:
                it.setFlags(it.flags() | Qt.ItemFlag.ItemIsUserCheckable)
                it.setCheckState(Qt.CheckState.Checked if n in cur else Qt.CheckState.Unchecked)
            self.list.addItem(it)
            self.names.append(n)
        if multi:
            self.list.setSelectionMode(QAbstractItemView.SelectionMode.NoSelection)
            self.list.itemChanged.connect(self.item_changed)
        else:
            self.list.setSelectionMode(QAbstractItemView.SelectionMode.SingleSelection)
            for i in range(self.list.count()):
                if self.list.item(i).text() in cur:
                    self.list.setCurrentRow(i)
                    break
        lay.addWidget(self.list, 1)
        bb = QDialogButtonBox(QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel)
        clear = bb.addButton(tr("strains_none"), QDialogButtonBox.ButtonRole.ResetRole)
        clear.clicked.connect(self.clear_all)
        bb.accepted.connect(self.accept)
        bb.rejected.connect(self.reject)
        lay.addWidget(bb)
        self.update_info()
        self.resize(460, 520)

    def checked_count(self):
        return sum(1 for i in range(self.list.count())
                   if self.names[i] not in self.removed and self.list.item(i).checkState() == Qt.CheckState.Checked)

    def usable_count(self):
        return sum(1 for n in self.names if n not in self.removed)

    def item_changed(self, it):
        if self._sync:
            return
        if it.checkState() == Qt.CheckState.Checked and self.usable_count() - self.checked_count() < self.keep:
            self._sync = True
            it.setCheckState(Qt.CheckState.Unchecked)   # at least 'keep' strains must stay unticked
            self._sync = False
        self.update_info()

    def clear_all(self):
        self._sync = True
        for i in range(self.list.count()):
            it = self.list.item(i)
            if self.names[i] in self.removed:
                continue
            if self.multi:
                it.setCheckState(Qt.CheckState.Unchecked)
            else:
                it.setSelected(False)
        self.list.clearSelection()
        self._sync = False
        self.update_info()

    def update_info(self):
        n = self.usable_count()
        if self.multi:
            text = tr("strains_multi", n=n, k=self.checked_count())
            if self.keep:
                text += " " + tr("strains_keep", m=self.keep)
        else:
            text = tr("strains_single", n=n)
        if self.removed:
            text += "\n" + tr("strains_removed_note", r=len(self.removed))
        self.info.setText(text)

    def selected(self):
        out = []
        for i in range(self.list.count()):
            it = self.list.item(i)
            if self.names[i] in self.removed:
                continue
            if self.multi:
                if it.checkState() == Qt.CheckState.Checked:
                    out.append(self.names[i])
            elif it.isSelected():
                out.append(self.names[i])
        return out[:1] if not self.multi else out


def table_widget(headers, rows):
    t = QTableWidget(len(rows), len(headers))
    t.setHorizontalHeaderLabels(headers)
    for r, row in enumerate(rows):
        for c, v in enumerate(row):
            it = QTableWidgetItem(str(v))
            it.setFlags(it.flags() & ~Qt.ItemFlag.ItemIsEditable)
            t.setItem(r, c, it)
    t.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Interactive)
    t.resizeColumnsToContents()
    t.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
    return t


class DryRunDialog(QDialog):
    def __init__(self, data, parent=None):
        super().__init__(parent)
        self.setWindowTitle(tr("dlg_dryrun"))
        lay = QVBoxLayout(self)
        g = data.get("genomes", [])
        lay.addWidget(QLabel(tr("dry_summary", n=len(g), d=data.get("input", ""))))
        rows = [[x["name"], tr("type_" + x["type"]), ", ".join(x.get("files", [])), ", ".join(x.get("ignored", [])) or "-"] for x in g]
        lay.addWidget(table_widget([tr("dry_genome"), tr("dry_type"), tr("dry_files"), tr("dry_ignored")], rows), 1)
        msgs = ["NOTE: " + n for n in data.get("notes", [])] + ["WARNING: " + w for w in data.get("warnings", [])]
        if msgs:
            te = QPlainTextEdit("\n".join(msgs))
            te.setReadOnly(True)
            te.setMaximumHeight(140)
            lay.addWidget(te)
        bb = QDialogButtonBox(QDialogButtonBox.StandardButton.Close)
        bb.rejected.connect(self.reject)
        bb.accepted.connect(self.accept)
        lay.addWidget(bb)
        self.resize(900, 500)


class TableDialog(QDialog):
    """Shows a tab-separated (or comma-separated) text file as a table; '#' lines are shown above it."""

    def __init__(self, path, parent=None):
        super().__init__(parent)
        self.path = path
        self.setWindowTitle("%s - %s" % (tr("dlg_table"), os.path.basename(path)))
        lay = QVBoxLayout(self)
        comments, rows = [], []
        sep = "," if path.lower().endswith(".csv") else "\t"
        with open(path, encoding="utf-8", errors="replace") as fh:
            for i, line in enumerate(fh):
                line = line.rstrip("\r\n")
                if line.startswith("#"):
                    comments.append(line)
                elif line.strip():
                    rows.append(line.split(sep))
                if len(rows) > 20000:
                    break
        if comments:
            lab = QLabel("\n".join(comments[:8]))
            lab.setWordWrap(True)
            lay.addWidget(lab)
        if rows:
            n = max(len(r) for r in rows)
            head = rows[0] + [""] * (n - len(rows[0]))
            lay.addWidget(table_widget(head, rows[1:]), 1)
        b = QPushButton(tr("open_system"))
        b.clicked.connect(lambda: QDesktopServices.openUrl(QUrl.fromLocalFile(self.path)))
        lay.addWidget(b)
        self.resize(1000, 600)


class CheckDialog(QDialog):
    """Installation check: scripts, Perl, programs, HMM sets, tree options, script options."""
    result_ready = Signal(list)

    def __init__(self, win, parent=None):
        super().__init__(parent)
        self.win = win
        self.setWindowTitle(tr("dlg_check"))
        lay = QVBoxLayout(self)
        self.info = QLabel(tr("checking"))
        lay.addWidget(self.info)
        self.table = QTableWidget(0, 3)
        self.table.setHorizontalHeaderLabels([tr("col_item"), tr("col_status"), tr("col_details")])
        self.table.horizontalHeader().setSectionResizeMode(2, QHeaderView.ResizeMode.Stretch)
        lay.addWidget(self.table, 1)
        row = QHBoxLayout()
        cb = QPushButton(tr("copy_report"))
        cb.clicked.connect(self.copy_report)
        row.addWidget(cb)
        row.addStretch(1)
        close = QPushButton(tr("close"))
        close.clicked.connect(self.accept)
        row.addWidget(close)
        lay.addLayout(row)
        self.rows = []
        self.result_ready.connect(self.show_rows)
        self.resize(1000, 560)
        threading.Thread(target=lambda: self.result_ready.emit(run_checks(self.win.config, self.win.tabs_list)), daemon=True).start()

    def show_rows(self, rows):
        self.rows = rows
        self.table.setRowCount(len(rows))
        colors = {"OK": QColor(0, 130, 0), "WARNING": QColor(200, 110, 0), "ERROR": QColor(200, 0, 0), "INFO": QColor(60, 60, 60)}
        for r, (item, st, det) in enumerate(rows):
            for c, v in enumerate((item, st, det)):
                it = QTableWidgetItem(v)
                it.setFlags(it.flags() & ~Qt.ItemFlag.ItemIsEditable)
                if c == 1:
                    it.setForeground(QBrush(colors.get(st, QColor(0, 0, 0))))
                self.table.setItem(r, c, it)
        self.table.resizeColumnToContents(0)
        self.table.resizeColumnToContents(1)
        n_err = sum(1 for x in rows if x[1] == "ERROR")
        n_warn = sum(1 for x in rows if x[1] == "WARNING")
        self.info.setText("ERROR: %d   WARNING: %d" % (n_err, n_warn))

    def copy_report(self):
        QGuiApplication.clipboard().setText("\n".join("\t".join(r) for r in self.rows))


def program_probe_args(name):
    n = name.lower()
    if "fasttree" in n:
        return ["-expert"]
    if "iqtree" in n:
        return ["--version"]
    if "astral" in n:
        return ["-h"]
    return {"prodigal": ["-v"], "hmmsearch": ["-h"], "hmmbuild": ["-h"], "muscle5": ["-version"], "trimal": ["--version"]}.get(n, ["-h"])


def version_line(name, text):
    for line in text.splitlines():
        if "HMMER" in line and re.search(r"\d\.\d", line):
            return line.strip("# ").strip()
    for line in text.splitlines():
        if re.search(r"\d+\.\d+", line):
            return line.strip()
    return text.strip().splitlines()[0] if text.strip() else ""


def run_checks(config, tabs):
    """List of (item, status, details)."""
    rows = []
    add = lambda item, st, det="": rows.append((item, st, det))
    home = Path(config.home) if config.home else None
    if not config.home_ok():
        add("EasyCGTree directory", "ERROR", str(home or "-") + ": EasyCGTree.pl not found")
    else:
        add("EasyCGTree directory", "OK", str(home))
    perl = config.perl_path()
    perl_ok = False
    if not perl:
        add("Perl", "ERROR", tr("msg_no_perl").splitlines()[0])
    else:
        rc, out = capture([perl, "-e", "print $]"], timeout=30)
        try:
            ver = float(out.strip())
        except ValueError:
            ver = 0
        if rc == 0 and ver >= 5.014:
            perl_ok = True
            add("Perl", "OK", "%s (%s)" % (perl, out.strip()))
        else:
            add("Perl", "ERROR", "%s: Perl >= 5.14 is needed (%s)" % (perl, out.strip()[:200]))
    # scripts and their options
    for key, script in SCRIPTS.items():
        p = config.script(key)
        if not p.is_file():
            add(script, "ERROR", "not found in %s" % p.parent)
            continue
        if not perl_ok:
            add(script, "WARNING", "found; not checked (no Perl)")
            continue
        rc, out = capture([perl, str(p), "-options_json"], timeout=60)
        m = re.search(r"\{.*\}", out, re.S)
        if rc != 0 or not m:
            add(script, "WARNING", "found, but '-options_json' is not supported (an older version?). "
                "The GUI is made for the version 5.0 scripts.")
            continue
        try:
            spec = json.loads(m.group(0))
        except ValueError:
            add(script, "WARNING", "'-options_json' gave no valid JSON")
            continue
        sopts = {o["name"]: o for o in spec.get("options", [])}
        gui_flags = set()
        for t in tabs:
            gui_flags |= t.flags_for_script(key)
        unknown = sorted(f for f in gui_flags if f not in sopts)
        unused = sorted(f for f in sopts if f not in gui_flags and f not in ("dry_run", "no_tree"))
        diffs = []
        for t in tabs:
            for prm in t.params_for_script(key):
                o = sopts.get(prm.flag)
                if o is None or prm.default is None or o.get("default") is None or prm.kind == "bool" \
                        or prm.flag == "thread" or getattr(prm, "autofilled", False):
                    # the thread default depends on the CPU; autofilled defaults come from the previous run
                    continue
                try:
                    same = float(o["default"]) == float(prm.default)
                except (TypeError, ValueError):
                    same = str(o["default"]) == str(prm.default)
                if not same:
                    diffs.append("-%s (script %s, GUI %s)" % (prm.flag, o["default"], prm.default))
        if unknown:
            add(script, "ERROR", "version %s: options used by the GUI are unknown to the script: %s" % (spec.get("version"), ", ".join("-" + u for u in unknown)))
        elif diffs:
            add(script, "WARNING", "version %s: different defaults: %s" % (spec.get("version"), "; ".join(sorted(set(diffs)))))
        else:
            add(script, "OK", "version %s" % spec.get("version") + ("; options not in the GUI: %s" % ", ".join("-" + u for u in unused) if unused else ""))
    # programs
    bin_dir = config.bin_dir()
    add("bin directory", "OK" if bin_dir.is_dir() else "ERROR", str(bin_dir))
    versions = {}
    for name in config.programs():
        p = config.program_path(name)
        label = "%s%s" % (name, " (custom)" if config.overrides.get(name) else "")
        if not p.is_file():
            add(label, "ERROR", "not found: %s" % p)
            continue
        if not IS_WIN and not os.access(p, os.X_OK):
            add(label, "ERROR", "%s is not executable (chmod +x \"%s\")" % (p, p))
            continue
        rc, out = capture([str(p)] + program_probe_args(name), timeout=30)
        v = version_line(name, out)
        if rc == -1:
            add(label, "ERROR", "%s can't be started: %s" % (p, out.strip()[:200]))
        elif rc == -2:
            add(label, "WARNING", "%s: no answer within 30 s" % p)
        elif not v:
            add(label, "WARNING", "%s: no version information (exit code %d)" % (p, rc))
        else:
            versions[name] = v
            st = "OK"
            det = "%s  [%s]" % (v, p)
            if "iqtree" in name.lower() and re.search(r"version 2\.2\.0", v):
                st, det = "WARNING", det + "  IQ-TREE 2.2.0 can fail with several threads ('Tree taxa and alignment sequence do not match'); IQ-TREE 3 is recommended."
            add(label, st, det)
    hv = {k: re.search(r"HMMER (\d+\.\d+(\.\d+)?)", versions.get(k, "")) for k in ("hmmsearch", "hmmbuild")}
    hv = {k: (m.group(1) if m else "") for k, m in hv.items()}
    if hv["hmmsearch"] and hv["hmmbuild"] and hv["hmmsearch"] != hv["hmmbuild"]:
        add("HMMER", "WARNING", "hmmbuild (%s) and hmmsearch (%s) differ; HMMs built by HMMER >= 3.1 cannot be read by HMMER 3.0." % (hv["hmmbuild"], hv["hmmsearch"]))
    # HMM sets
    sets = hmm_sets(config.hmm_dir())
    if not sets:
        add("HMM sets", "ERROR", "no .hmm file in %s. The HMM sets are downloaded separately from the folder 'HMM' of "
            "the EasyCGTree5 repository (Help - HMM sets online; see HMM/README.txt)." % config.hmm_dir())
    else:
        if not any(n == "bac120" for n, k, fmt in sets):
            add("HMM sets", "WARNING", "the default set 'bac120' is missing; download 'bac120.hmm' from the folder 'HMM' "
                "of the EasyCGTree5 repository, or choose another set.")
        bad = [n for n, k, fmt in sets if hv["hmmsearch"].startswith("3.0") and hmm_format_version(fmt) not in ("", "b", "a")]
        empty = [n for n, k, fmt in sets if k == 0]
        add("HMM sets", "WARNING" if bad or empty else "OK", "%d sets in %s" % (len(sets), config.hmm_dir()) +
            ("; not readable by HMMER 3.0: %s" % ", ".join(bad) if bad else "") + ("; without profiles: %s" % ", ".join(empty) if empty else ""))
    # tree options
    tf = config.tree_options_file()
    opts = read_tree_options(tf)
    missing = [k for k, _ in TREE_KEYS if not opts.get(k)]
    if not tf.is_file():
        add("tree_app-options.txt", "ERROR", "not found: %s" % tf)
    elif missing:
        add("tree_app-options.txt", "WARNING", "%s: missing command line(s): %s" % (tf, ", ".join(missing)))
    else:
        add("tree_app-options.txt", "OK", "%s: %s" % (tf, " | ".join(opts[k] for k, _ in TREE_KEYS)))
    env = config.env_overrides()
    if env:
        add("Environment", "INFO", "; ".join("%s=%s" % kv for kv in env.items()))
    return rows


class SettingsDialog(QDialog):
    def __init__(self, win, parent=None):
        super().__init__(parent)
        self.win = win
        cfg = win.config
        self.i18n = I18n()
        self.setWindowTitle(tr("dlg_settings"))
        lay = QVBoxLayout(self)
        grid = QGridLayout()
        self.home = PathEdit(self.i18n, "dir")
        self.home.setText(cfg.home)
        self.bin = PathEdit(self.i18n, "dir")
        self.bin.setText(cfg.bin)
        self.bin.edit.setPlaceholderText("<EasyCGTree>/bin")
        self.perl = PathEdit(self.i18n, "file")
        self.perl.setText(cfg.perl)
        self.perl.edit.setPlaceholderText(find_perl("") or tr("perl_not_found"))
        for r, (key, w, hk) in enumerate([("set_home", self.home, "h_set_home"), ("set_bin", self.bin, "h_set_bin"),
                                          ("set_perl", self.perl, "h_set_perl")]):
            grid.addWidget(QLabel(tr(key)), r, 0)
            grid.addWidget(w, r, 1)
            grid.addWidget(HelpButton(hk, p=find_perl("") or tr("perl_not_found")), r, 2)
        lay.addLayout(grid)
        lay.addWidget(QLabel(tr("set_programs")))
        self.table = QTableWidget(0, 3)
        self.table.setHorizontalHeaderLabels([tr("col_program"), tr("col_path"), tr("col_used")])
        self.table.horizontalHeader().setSectionResizeMode(1, QHeaderView.ResizeMode.Stretch)
        self.table.horizontalHeader().setSectionResizeMode(2, QHeaderView.ResizeMode.Stretch)
        lay.addWidget(self.table, 1)
        row = QHBoxLayout()
        b1 = QPushButton(tr("browse"))
        b1.clicked.connect(self.browse_program)
        b2 = QPushButton(tr("reset"))
        b2.clicked.connect(self.reset_program)
        b3 = QPushButton(tr("menu_check").rstrip("."))
        b3.clicked.connect(self.check)
        row.addWidget(b1)
        row.addWidget(b2)
        row.addStretch(1)
        row.addWidget(b3)
        lay.addLayout(row)
        bb = QDialogButtonBox(QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel)
        bb.accepted.connect(self.accept_settings)
        bb.rejected.connect(self.reject)
        lay.addWidget(bb)
        self.overrides = dict(cfg.overrides)
        self.fill_table()
        self.home.edit.textChanged.connect(lambda *_: self.fill_table())
        self.bin.edit.textChanged.connect(lambda *_: self.fill_table())
        self.resize(900, 560)

    def temp_config(self):
        c = Config.__new__(Config)
        c.__dict__.update(self.win.config.__dict__)
        c.home = self.home.text()
        c.bin = self.bin.text()
        c.perl = self.perl.text()
        c.overrides = {k: v for k, v in self.overrides.items() if v}
        return c

    def fill_table(self):
        c = self.temp_config()
        progs = c.programs()
        for k in self.overrides:
            if k not in progs:
                progs.append(k)
        self.table.setRowCount(len(progs))
        for r, name in enumerate(progs):
            items = [name, self.overrides.get(name, ""), str(c.program_path(name)) + ("" if c.program_path(name).is_file() else "  (!)")]
            for col, v in enumerate(items):
                it = QTableWidgetItem(v)
                if col != 1:
                    it.setFlags(it.flags() & ~Qt.ItemFlag.ItemIsEditable)
                self.table.setItem(r, col, it)

    def read_table(self):
        for r in range(self.table.rowCount()):
            name = self.table.item(r, 0).text()
            it = self.table.item(r, 1)
            self.overrides[name] = it.text().strip() if it else ""

    def browse_program(self):
        r = self.table.currentRow()
        if r < 0:
            return
        self.read_table()
        p, _ = QFileDialog.getOpenFileName(self, tr("browse"), str(self.temp_config().bin_dir()))
        if p:
            self.overrides[self.table.item(r, 0).text()] = os.path.normpath(p)
            self.fill_table()

    def reset_program(self):
        r = self.table.currentRow()
        if r < 0:
            return
        self.read_table()
        self.overrides[self.table.item(r, 0).text()] = ""
        self.fill_table()

    def apply_to(self, cfg):
        self.read_table()
        cfg.home = self.home.text()
        cfg.bin = self.bin.text()
        cfg.perl = self.perl.text()
        cfg.overrides = {k: v for k, v in self.overrides.items() if v}

    def check(self):
        self.apply_to(self.win.config)
        self.win.config_changed()
        CheckDialog(self.win, self).exec()

    def accept_settings(self):
        self.apply_to(self.win.config)
        self.win.config_changed()
        self.accept()


class TreeOptionsDialog(QDialog):
    def __init__(self, win, parent=None):
        super().__init__(parent)
        self.win = win
        self.file = win.config.tree_options_file()
        self.setWindowTitle(tr("dlg_treeopts"))
        lay = QVBoxLayout(self)
        info = QLabel(tr("treeopts_info", f=self.file) if self.file.is_file() else tr("treeopts_missing", f=self.file))
        info.setWordWrap(True)
        lay.addWidget(info)
        opts = read_tree_options(self.file)
        grid = QGridLayout()
        self.edits = {}
        for r, (key, prefix) in enumerate(TREE_KEYS):
            grid.addWidget(QLabel(prefix.rstrip("=")), r, 0)
            e = QLineEdit(opts.get(key, ""))
            e.setFont(fixed_font())
            grid.addWidget(e, r, 1)
            self.edits[key] = e
        lay.addLayout(grid)
        ex = QLabel("IQ-TREE: iqtree3 -B 1000  |  iqtree3 --alrt 1000 -B 1000  |  iqtree3 -m LG+G4 -B 1000\n"
                    "FastTree: FastTreeMP  |  FastTreeMP -gamma\nASTRAL: astral-weighted -S -x 100")
        ex.setFont(fixed_font())
        lay.addWidget(ex)
        bb = QDialogButtonBox(QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel)
        bb.accepted.connect(self.save)
        bb.rejected.connect(self.reject)
        lay.addWidget(bb)
        bb.button(QDialogButtonBox.StandardButton.Ok).setEnabled(self.file.is_file())
        self.resize(800, 260)

    def save(self):
        try:
            write_tree_options(self.file, {k: e.text() for k, e in self.edits.items() if e.text().strip()})
        except OSError as e:
            QMessageBox.critical(self, tr("error"), str(e))
            return
        self.win.config_changed()
        self.accept()


# ======================================================================================================
# File browser
# ======================================================================================================

class FileBrowser(QWidget):
    """File view with buttons to switch between places (input directory, output directory, ...)."""

    def __init__(self, tab, places):
        super().__init__()
        self.tab = tab
        self.places = places          # [(key, label key, function returning the directory)]
        self.place = places[0][0]
        self.applied = {}
        i18n = tab.i18n
        lay = QVBoxLayout(self)
        lay.setContentsMargins(0, 0, 0, 0)
        # the place buttons are shown by the tab above both the status and the file view
        self.places_bar = QWidget()
        pb = QHBoxLayout(self.places_bar)
        pb.setContentsMargins(0, 0, 0, 0)
        pb.setSpacing(0)
        self.place_buttons = {}
        grp = QButtonGroup(self)
        grp.setExclusive(True)
        for i, (key, lkey, fn) in enumerate(places):
            b = QToolButton()
            b.setObjectName("place")
            b.setProperty("pos", "first" if i == 0 else ("last" if i == len(places) - 1 else "mid"))
            b.setCheckable(True)
            i18n.add(b.setText, lkey)
            b.clicked.connect(lambda checked=False, k=key: self.select_place(k))
            grp.addButton(b)
            pb.addWidget(b)
            self.place_buttons[key] = b
        self.place_buttons[self.place].setChecked(True)
        row = QHBoxLayout()
        up = QToolButton()
        i18n.add(up.setText, "up")
        up.clicked.connect(self.go_up)
        row.addWidget(up)
        self.path = QLineEdit()
        self.path.returnPressed.connect(lambda: self.set_dir(self.path.text()))
        row.addWidget(self.path, 1)
        lay.addLayout(row)
        self.model = QFileSystemModel()
        self.model.setRootPath("")
        self.view = QTreeView()
        self.view.setModel(self.model)
        self.view.setSortingEnabled(True)
        self.view.sortByColumn(0, Qt.SortOrder.AscendingOrder)
        self.view.setColumnHidden(2, True)
        self.view.setColumnWidth(0, 320)
        self.view.doubleClicked.connect(self.activated)
        self.view.setContextMenuPolicy(Qt.ContextMenuPolicy.CustomContextMenu)
        self.view.customContextMenuRequested.connect(self.menu)
        lay.addWidget(self.view, 1)
        self.current = ""

    def place_dir(self, key):
        for k, lkey, fn in self.places:
            if k == key:
                return fn() or ""
        return ""

    def select_place(self, key):
        self.place = key
        self.place_buttons[key].setChecked(True)
        self.applied[key] = None
        self.refresh_place()
        self.tab.refresh_status()

    def refresh_place(self):
        """Show the directory of the selected place (again, if it has changed)."""
        d = self.place_dir(self.place)
        if d != self.applied.get(self.place):
            self.applied[self.place] = d
            if d and os.path.isdir(d):
                self.set_dir(d)
            else:
                self.current = ""
                self.path.setText(d)

    def set_dir(self, d):
        if not d or not os.path.isdir(d):
            return
        d = os.path.normpath(d)
        self.current = d
        self.path.setText(d)
        self.model.setRootPath(d)
        self.view.setRootIndex(self.model.index(d))

    def go_up(self):
        if self.current:
            self.set_dir(os.path.dirname(self.current))

    def activated(self, index):
        p = self.model.filePath(index)
        if self.model.isDir(index):
            self.set_dir(p)
        else:
            open_file(p, self)

    def menu(self, pos):
        idx = self.view.indexAt(pos)
        if not idx.isValid():
            return
        p = self.model.filePath(idx)
        m = QMenu(self)
        a1 = m.addAction(tr("open"))
        a2 = m.addAction(tr("open_system"))
        a3 = m.addAction(tr("show_folder"))
        a4 = m.addAction(tr("copy_path"))
        act = m.exec(self.view.viewport().mapToGlobal(pos))
        if act == a1:
            self.activated(idx)
        elif act == a2:
            QDesktopServices.openUrl(QUrl.fromLocalFile(p))
        elif act == a3:
            QDesktopServices.openUrl(QUrl.fromLocalFile(p if os.path.isdir(p) else os.path.dirname(p)))
        elif act == a4:
            QGuiApplication.clipboard().setText(p)


def open_file(p, parent):
    """Tables are shown in the GUI; everything else (trees, logs, FASTA) with the system application."""
    low = p.lower()
    if low.endswith((".txt", ".tsv", ".csv")):
        try:
            with open(p, encoding="utf-8", errors="replace") as fh:
                head = [fh.readline() for _ in range(5)]
            if any("\t" in h or (low.endswith(".csv") and "," in h) for h in head):
                TableDialog(p, parent).exec()
                return
        except OSError:
            pass
    QDesktopServices.openUrl(QUrl.fromLocalFile(p))


def dir_key(d):
    """Identifies the state of a directory (file names, sizes, times) for caching."""
    out = []
    for f in input_files(d):
        try:
            st = (Path(d) / f).stat()
            out.append((f, st.st_size, int(st.st_mtime)))
        except OSError:
            pass
    return (str(d), tuple(out))


def input_summary_text(input_dir, perl, script, env):
    """Input directory: how the files are grouped into genomes (dry run of EasyCGTree) and their statistics."""
    files = input_files(input_dir)
    stats = {f: fasta_stats(Path(input_dir) / f) for f in files}
    lines = [tr("in_dir", d=input_dir)]
    data = None
    err = ""
    if perl and Path(script).is_file():
        rc, out = capture([perl, str(script), "-input", str(input_dir), "-dry_run"], env=env, timeout=600)
        m = re.search(r"^EasyCGTree_DRY_RUN_JSON=(.*)$", out, re.M)
        if rc == 0 and m:
            try:
                data = json.loads(m.group(1))
            except ValueError:
                data = None
        if data is None:
            err = "\n".join(l for l in out.splitlines() if l.startswith("ERROR") or (l and not l.startswith(("#", " ")) and "Reading command line" not in l))[-1500:]
    if data is None:
        lines.append(tr("in_files", n=len(files)))
        lines.append("")
        lines.append("%-40s %10s %14s" % (tr("col_file"), tr("col_seqs"), tr("col_length")))
        for f in files:
            n, total = stats[f]
            lines.append("%-40s %10s %14s" % (f, format(n, ","), format(total, ",")))
        if err:
            lines += ["", err]
        return "\n".join(lines)
    g = data.get("genomes", [])
    counts = {}
    for x in g:
        counts[x["type"]] = counts.get(x["type"], 0) + 1
    lines.append(tr("in_genomes", f=len(files), n=len(g)) + "  (" + ", ".join("%s: %d" % (tr("type_short_" + t), c) for t, c in sorted(counts.items())) + ")")
    lines.append("")
    lines.append("%-28s %-10s %10s %14s  %s" % (tr("dry_genome"), tr("dry_type"), tr("col_seqs"), tr("col_length"), tr("col_file")))
    for x in g:
        fl = x.get("files", [])
        for i, f in enumerate(fl):
            n, total = stats.get(f, (0, 0))
            lines.append("%-28s %-10s %10s %14s  %s" % (x["name"] if i == 0 else "", x["type"] if i == 0 else "", format(n, ","), format(total, ","), f))
        for f in x.get("ignored", []):
            lines.append("%-28s %-10s %10s %14s  %s" % ("", "", "", "", f + "  (" + tr("ignored") + ")"))
    for n in data.get("notes", []):
        lines.append("NOTE: " + n)
    for w in data.get("warnings", []):
        lines.append("WARNING: " + w)
    return "\n".join(lines)


def family_summary_text(gc_dir):
    """Gene family directory (HMM building): sequences and mean length of each family."""
    files = input_files(gc_dir)
    lines = [tr("in_dir", d=gc_dir), tr("fam_files", n=len(files)), ""]
    lines.append("%-30s %10s %12s  %s" % (tr("col_gene"), tr("col_seqs"), tr("col_mean_length"), tr("col_file")))
    for f in files:
        n, total = fasta_stats(Path(gc_dir) / f)
        gene = re.sub(r"\.(fas|fasta|fa|faa|afa|aln|fst|txt|seq|pep)$", "", f, flags=re.I)
        gene = re.sub(r"[\s/\\:*?\"<>|()\[\],;]+", "-", gene).strip("-")
        lines.append("%-30s %10s %12s  %s" % (gene, format(n, ","), format(int(round(total / n)) if n else 0, ","), f))
    return "\n".join(lines)


# ======================================================================================================
# Tabs
# ======================================================================================================

class BaseTab(QWidget):
    """Parameter groups (left), status and files (right), command line and Run/Stop (bottom), log."""
    name = ""
    summary_ready = Signal(object, str)

    def __init__(self, win):
        super().__init__()
        self.win = win
        self.i18n = win.i18n
        self.params = {}
        self.groups = {}
        self.grids = []
        self._building = True
        self._sum_cache = {}
        self._sum_pending = None
        self.summary_ready.connect(self.on_summary)
        outer = QVBoxLayout(self)
        outer.setContentsMargins(8, 8, 8, 8)
        # Four panels: parameters | status/files (top), command line | log (bottom). The top and bottom
        # halves can be resized against each other, and the two panels of each half side by side.
        vsplit = QSplitter(Qt.Orientation.Vertical)
        top = QSplitter(Qt.Orientation.Horizontal)
        for sp in (vsplit, top):
            sp.setChildrenCollapsible(False)
            sp.setHandleWidth(8)
        self.form = QWidget()
        self.form.setObjectName("paramForm")
        self.form_layout = QVBoxLayout(self.form)
        self.form_layout.setSpacing(10)
        scroll = QScrollArea()
        scroll.setObjectName("paramScroll")
        scroll.setWidgetResizable(True)
        scroll.setFrameShape(QFrame.Shape.NoFrame)
        scroll.setWidget(self.form)
        scroll.setMinimumWidth(320)
        top.addWidget(scroll)
        # right: 'Status' and 'Files' as two pages of one panel; the place buttons (Output/Input ...) above
        # them choose the folder for both pages
        right = QFrame()
        right.setObjectName("panel")
        rlay = QVBoxLayout(right)
        rlay.setContentsMargins(10, 8, 10, 10)
        rlay.setSpacing(6)
        head = QHBoxLayout()
        self.view_buttons = {}
        vgrp = QButtonGroup(self)
        vgrp.setExclusive(True)
        self.right_stack = QStackedWidget()
        for i, key in enumerate(("status", "files")):
            b = QToolButton()
            b.setObjectName("segment")
            b.setCheckable(True)
            self.i18n.add(b.setText, key)
            b.clicked.connect(lambda checked=False, k=i: self.show_right_page(k))
            vgrp.addButton(b)
            head.addWidget(b)
            self.view_buttons[i] = b
        head.addSpacing(12)
        head.addStretch(1)
        self.browser = FileBrowser(self, self.browser_places())
        head.addWidget(self.browser.places_bar)
        rlay.addLayout(head)
        self.status = QPlainTextEdit()
        self.status.setObjectName("statusView")
        self.status.setReadOnly(True)
        self.status.setLineWrapMode(QPlainTextEdit.LineWrapMode.NoWrap)
        self.status.setFont(fixed_font())
        self.right_stack.addWidget(self.status)
        self.right_stack.addWidget(self.browser)
        rlay.addWidget(self.right_stack, 1)
        right.setMinimumWidth(300)
        top.addWidget(right)
        self.show_right_page(win.config.right_page if win.config.right_page in (0, 1) else 0)
        top.setStretchFactor(0, 1)
        top.setStretchFactor(1, 1)
        top.setSizes([1000, 1000])          # left and right equally wide
        self.hsplit = top
        vsplit.addWidget(top)
        # bottom: command line (with Run/Stop) | log
        bottom = QSplitter(Qt.Orientation.Horizontal)
        bottom.setChildrenCollapsible(False)
        bottom.setHandleWidth(8)
        cg, chead, cl = self.panel("cmdline")
        # first row: title and Concise|Full
        self.rb_concise = QRadioButton()
        self.rb_full = QRadioButton()
        cbg = QButtonGroup(self)
        cbg.setExclusive(True)
        for b, key in ((self.rb_concise, "concise"), (self.rb_full, "full")):
            self.i18n.add(b.setText, key)
            self.i18n.add(b.setToolTip, key + "_tip")
            cbg.addButton(b)
        (self.rb_full if win.config.full_cmd else self.rb_concise).setChecked(True)
        self.rb_full.toggled.connect(self.changed)
        chead.insertSpacing(1, 14)
        chead.insertWidget(2, self.rb_concise)
        chead.insertSpacing(3, 6)
        chead.insertWidget(4, self.rb_full)
        # second row: Run, Stop and the state of the run on the left, Copy on the right
        brow = QHBoxLayout()
        self.run_btn = QPushButton()
        self.run_btn.setObjectName("runButton")
        self.i18n.add(self.run_btn.setText, "run")
        self.run_btn.clicked.connect(lambda: self.run())
        self.stop_btn = QPushButton()
        self.stop_btn.setObjectName("stopButton")
        self.i18n.add(self.stop_btn.setText, "stop")
        self.stop_btn.clicked.connect(self.win.stop)
        self.stop_btn.setEnabled(False)
        self.state_label = QLabel()
        self.state_label.setObjectName("stateLabel")
        self._state = None
        brow.addWidget(self.run_btn)
        brow.addWidget(self.stop_btn)
        brow.addSpacing(8)
        brow.addWidget(self.state_label, 1)
        self.copy_btn = QPushButton()
        self.i18n.add(self.copy_btn.setText, "copy")
        self.copy_btn.clicked.connect(self.copy_cmd)
        brow.addWidget(self.copy_btn)
        cl.addLayout(brow)
        self.cmd_view = QPlainTextEdit()
        self.cmd_view.setObjectName("cmdView")
        self.cmd_view.setReadOnly(True)
        self.cmd_view.setFont(fixed_font())
        cl.addWidget(self.cmd_view, 1)
        cg.setMinimumWidth(320)
        bottom.addWidget(cg)
        lg, lhead, ll = self.panel("log")
        clear = QPushButton()
        self.i18n.add(clear.setText, "clear_log")
        lhead.addWidget(clear)
        self.log = LogView()
        self.log.setObjectName("logView")
        clear.clicked.connect(self.log.reset)
        ll.addWidget(self.log, 1)
        lg.setMinimumWidth(280)
        bottom.addWidget(lg)
        bottom.setStretchFactor(0, 1)
        bottom.setStretchFactor(1, 1)
        bottom.setSizes([1000, 1000])
        self.bsplit = bottom
        vsplit.addWidget(bottom)
        vsplit.setStretchFactor(0, 3)
        vsplit.setStretchFactor(1, 2)
        vsplit.setSizes([600, 400])
        self.vsplit = vsplit
        outer.addWidget(vsplit)
        QTimer.singleShot(0, self.restore_splitters)   # once the window is shown
        self.status_timer = QTimer(self)
        self.status_timer.setSingleShot(True)
        self.status_timer.setInterval(400)
        self.status_timer.timeout.connect(self.refresh_status)
        self.build()
        self.form_layout.addStretch(1)
        self.align_labels()
        self._building = False
        self.changed()

    # ---- building helpers ----
    def strain_names(self, key):
        """Names offered for a strain option, or None for the genomes of the input directory."""
        return None

    def panel(self, title_key):
        """A white panel with a title row (like the Status/Files panel): (frame, title row, body layout)."""
        frame = QFrame()
        frame.setObjectName("panel")
        lay = QVBoxLayout(frame)
        lay.setContentsMargins(10, 8, 10, 10)
        lay.setSpacing(6)
        head = QHBoxLayout()
        title = QLabel()
        title.setObjectName("panelTitle")
        self.i18n.add(title.setText, title_key)
        head.addWidget(title)
        head.addStretch(1)
        lay.addLayout(head)
        return frame, head, lay

    def equal_halves(self):
        for sp in (self.hsplit, self.bsplit):
            w = sum(sp.sizes()) or 2000
            sp.setSizes([w // 2, w - w // 2])

    def splitters(self):
        return {"v": self.vsplit, "top": self.hsplit, "bottom": self.bsplit}

    def restore_splitters(self):
        """Sizes of the four panels as the user left them (shared by all tabs), else equal halves."""
        saved = self.win.config.splitters or {}
        self.equal_halves()
        for k, sp in self.splitters().items():
            sizes = saved.get(k)
            if isinstance(sizes, list) and len(sizes) == 2 and all(isinstance(x, int) and x > 0 for x in sizes):
                sp.setSizes(sizes)

    def save_splitters(self):
        self.win.config.splitters = {k: list(sp.sizes()) for k, sp in self.splitters().items()}

    def show_right_page(self, i):
        self.right_stack.setCurrentIndex(i)
        self.view_buttons[i].setChecked(True)
        self.win.config.right_page = i

    def set_run_state(self, key, kind, **kw):
        """Text and colour of the state label next to Run/Stop (kept when the language changes)."""
        self._state = (key, kind, kw)
        self.state_label.setText(tr(key, **kw))
        self.state_label.setProperty("state", kind)
        st = self.state_label.style()
        st.unpolish(self.state_label)
        st.polish(self.state_label)

    def retranslate_state(self):
        if self._state:
            key, kind, kw = self._state
            self.state_label.setText(tr(key, **kw))

    def browser_places(self):
        return [("output", "place_output", self.output_dir), ("input", "place_input", self.input_dir)]

    def add_group(self, gkey, rows, help_key=None, checkable=False, checked=True):
        g = QGroupBox()
        self.i18n.add(g.setTitle, gkey)
        if checkable:
            g.setCheckable(True)
            g.setChecked(checked)
            g.toggled.connect(self.changed)
        grid = QGridLayout(g)
        grid.setColumnStretch(1, 1)
        self.grids.append(grid)
        r = 0
        if help_key:
            lab = QLabel()
            lab.setWordWrap(True)
            lab.setStyleSheet("color: gray")
            self.i18n.add(lab.setText, help_key)
            grid.addWidget(lab, r, 0, 1, 3)
            r += 1
        for item in rows:
            if isinstance(item, Param):
                lab, w, hb = item.build()
                self.params[item.key] = item
                # the '?' follows the field; numbers and choices keep their own width
                if item.kind == "tasks":
                    w.layout().addWidget(hb)          # after the task boxes and '= all'
                    grid.addWidget(lab, r, 0)
                    grid.addWidget(w, r, 1, 1, 2)
                else:
                    cell = QWidget()
                    hl = QHBoxLayout(cell)
                    hl.setContentsMargins(0, 0, 0, 0)
                    hl.setSpacing(6)
                    short = item.kind in ("bool", "int", "float", "efloat", "choice")
                    hl.addWidget(w, 0 if short else 1)
                    hl.addWidget(hb, 0, Qt.AlignmentFlag.AlignVCenter)
                    if short:
                        hl.addStretch(1)
                    if item.kind == "bool":
                        grid.addWidget(cell, r, 0, 1, 3)
                    else:
                        grid.addWidget(lab, r, 0)
                        grid.addWidget(cell, r, 1, 1, 2)
            else:  # a widget spanning the row
                grid.addWidget(item, r, 0, 1, 3)
            r += 1
        self.form_layout.addWidget(g)
        self.groups[gkey] = g
        return g

    LABEL_MAX = 210   # longer option names are wrapped

    def align_labels(self):
        """The fields of all groups start at the same position: one width for the label column."""
        labels = [p.label for p in self.params.values() if p.label is not None and p.kind != "bool"]
        try:
            for lab in labels:
                lab.setWordWrap(False)
            w = max([lab.sizeHint().width() for lab in labels] + [0])
            w = min(int(w), self.LABEL_MAX)
            for lab in labels:
                if lab.sizeHint().width() > w:
                    lab.setWordWrap(True)
            for grid in self.grids:
                grid.setColumnMinimumWidth(0, w)
            for p in self.params.values():
                if isinstance(p.field, QComboBox):
                    fit_popup(p.field)
                    p.field.updateGeometry()
        except (TypeError, ValueError, AttributeError):
            pass

    def thread_param(self, *scripts, help=None):
        n = cpu_count()
        return mark(Param(self, "thread", "thread", "int", min(4, n), minimum=1, maximum=n, help=help or "h_thread"), *scripts)

    def outdir_param(self, *scripts, help="h_outdir"):
        return mark(Param(self, "outdir", "outdir", "dir", None, label="outdir", help=help), *scripts)

    def p(self, key):
        return self.params[key]

    def v(self, key):
        return self.params[key].value()

    # ---- directories ----
    def input_dir(self):
        return ""

    def output_dir(self):
        """The '-outdir' directory, or the directory containing the input directory."""
        d = self.input_dir()
        if "outdir" in self.params and self.v("outdir"):
            return self.v("outdir")
        return str(io_dirs(d)[1]) if d else ""

    def outdir_args(self, full):
        if self.v("outdir"):
            return ["-outdir", self.v("outdir")]
        if full and self.input_dir():
            return ["-outdir", self.output_dir()]
        return []

    def update_outdir_hint(self):
        if "outdir" in self.params:
            d = self.input_dir()
            self.p("outdir").field.edit.setPlaceholderText(tr("default_outdir", d=str(io_dirs(d)[1])) if d else "")

    # ---- state ----
    def changed(self, *args):
        if self._building:
            return
        self.update_outdir_hint()
        self.update_state()
        self.win.config.full_cmd = self.rb_full.isChecked()
        try:
            text = format_commands(self.commands(self.rb_full.isChecked()), self.win.config)
        except Exception as e:  # never break the GUI because of the preview
            text = "(%s)" % e
        self.cmd_view.setPlainText(text)
        self.status_timer.start()

    def update_state(self):
        pass

    def after_refresh(self):
        """Tab-specific updates that depend on the files (prediction notice, values from the record)."""
        pass

    def refresh_status(self):
        self.after_refresh()
        self.browser.refresh_place()
        mode = self.browser.place
        if mode == "input":
            self.show_input_summary()
        else:
            self.status.setPlainText(self.status_text(mode))

    def status_text(self, mode):
        return ""

    # ---- summary of the input directory (computed in the background) ----
    def compute_input_summary(self, d):
        return ""

    def show_input_summary(self):
        d = self.input_dir()
        if not d:
            self.status.setPlainText(tr("st_no_input"))
            return
        if not os.path.isdir(d):
            self.status.setPlainText(tr("st_no_dir", d=d))
            return
        key = dir_key(d)
        if key in self._sum_cache:
            self.status.setPlainText(self._sum_cache[key])
            return
        if self._sum_pending == key:
            return
        self._sum_pending = key
        self.status.setPlainText(tr("checking"))
        threading.Thread(target=lambda: self.summary_ready.emit(key, self.compute_input_summary(d)), daemon=True).start()

    def on_summary(self, key, text):
        self._sum_cache[key] = text
        if self._sum_pending == key:
            self._sum_pending = None
        d = self.input_dir()
        if self.browser.place == "input" and d and key[0] == str(d):
            self.status.setPlainText(text)

    def copy_cmd(self):
        QGuiApplication.clipboard().setText(self.cmd_view.toPlainText())
        self.win.statusBar().showMessage(tr("copied"), 4000)

    def validate(self):
        for prm in self.params.values():
            e = prm.validate()
            if e:
                return e
        return None

    def commands(self, full):
        return []

    def run(self, cmds=None):
        err = self.validate()
        if err:
            QMessageBox.warning(self, tr("warning"), err)
            return
        self.win.start_job(self, cmds if cmds is not None else self.commands(self.rb_full.isChecked()))

    def job_cwd(self):
        d = self.output_dir()
        if d and os.path.isdir(d):
            return d
        d = self.input_dir()
        return str(io_dirs(d)[1]) if d else str(Path.home())

    # ---- parameters for presets and the installation check ----
    def get_state(self):
        st = {k: p.value() for k, p in self.params.items()}
        for k, g in self.groups.items():
            if g.isCheckable():
                st["group:" + k] = g.isChecked()
        st["place"] = self.browser.place
        return st

    def set_state(self, st):
        self._building = True
        for k, v in (st or {}).items():
            if k.startswith("group:"):
                g = self.groups.get(k[6:])
                if g is not None and g.isCheckable():
                    g.setChecked(bool(v))
            elif k in self.params:
                self.params[k].set_value(v)
        self.set_extra_state(st or {})
        self._building = False
        if (st or {}).get("place") in self.browser.place_buttons:
            self.browser.place = st["place"]
            self.browser.place_buttons[st["place"]].setChecked(True)
        self.changed()

    def set_extra_state(self, st):
        pass

    def params_for_script(self, key):
        return [p for p in self.params.values() if p.flag and key in getattr(p, "scripts", ())]

    def flags_for_script(self, key):
        return {p.flag for p in self.params_for_script(key)}

    def reload_hmm_lists(self):
        for p in self.params.values():
            if p.kind == "hmm":
                p.reload_hmm()


def mark(param, *scripts):
    param.scripts = scripts
    return param


def pipeline_status(input_dir, outdir="", snp=False):
    if not input_dir:
        return tr("st_no_input")
    inp, out, tem, name = io_dirs(input_dir, outdir)
    lines = [tr("st_out", d=out)]
    if not tem.is_dir():
        lines.append(tr("st_no_tem", tem=tem))
    else:
        lines.append(tr("st_tem", tem=tem))
        rec = read_record(tem)
        if not rec:
            lines.append(tr("st_old"))
        else:
            dn = lambda k: rec.get(k, "")[:16]
            pred = (tr("st_done") + " (%s; %s)" % (dn("predict_done"), tr("st_genomes", n=rec.get("genomes", "?")))) if rec.get("predict_done") else tr("st_not_done")
            lines.append(tr("st_predict", v=pred))
            hs = ("%s, E-value %s (%s)" % (rec.get("hmm"), rec.get("evalue"), dn("hmmsearch_done"))) if rec.get("hmmsearch_done") else tr("st_not_done")
            lines.append(tr("st_hmmsearch", v=hs))
            rf = ("genome_cutoff %s, gene_cutoff %s (%s)" % (rec.get("genome_cutoff"), rec.get("gene_cutoff"), dn("refine_done"))) if rec.get("refine_done") else tr("st_not_done")
            lines.append(tr("st_refine", v=rf))
            al = ("%s, trimAl -%s (%s)" % (rec.get("alignment_seq"), rec.get("trim"), dn("alignment_done"))) if rec.get("alignment_done") else tr("st_not_done")
            lines.append(tr("st_alignment", v=al))
    trees = list_files(out, name + ".", ".tree")
    lines.append(tr("st_trees", v=", ".join(trees) if trees else tr("st_none")))
    if snp:
        sn = list_files(out, name + ".", ".snp.fas")
        lines.append(tr("st_snp", v=", ".join(sn) if sn else tr("st_none")))
    return "\n".join(lines)


class MainTab(BaseTab):
    name = "main"
    TASK_USES = {
        "predict": {"seq", "keep_CDS_nucl"},
        "hmmsearch": {"hmm", "evalue", "thread"},
        "refine": {"genome_cutoff", "gene_cutoff", "tree"},
        "alignment": {"trim", "seq", "thread"},
        "tree_infer": {"tree", "tree_app", "thread"},
    }

    def build(self):
        self._notice_key = None        # settings for which the user closed the notice
        self._auto_unchecked = False   # 'predict' was unticked because of the notice
        g = QWidget()
        hl = FlowLayout(g)
        self.rb_pipe = QRadioButton()
        self.rb_prev = QRadioButton()
        self.i18n.add(self.rb_pipe.setText, "mode_pipeline")
        self.i18n.add(self.rb_prev.setText, "mode_prevalence")
        grp = QButtonGroup(self)
        grp.addButton(self.rb_pipe)
        grp.addButton(self.rb_prev)
        self.rb_pipe.setChecked(True)
        self.rb_pipe.toggled.connect(self.changed)
        lab = QLabel()
        self.i18n.add(lab.setText, "runmode")
        hl.addWidget(lab)
        hl.addWidget(self.rb_pipe)
        hl.addWidget(self.rb_prev)
        hl.addWidget(HelpButton("h_runmode"))
        self.check_btn = QPushButton()
        self.i18n.add(self.check_btn.setText, "check_input")
        self.check_btn.clicked.connect(self.check_input)
        cw = QWidget()
        cl = QHBoxLayout(cw)
        cl.setContentsMargins(0, 0, 0, 0)
        cl.addWidget(self.check_btn)
        cl.addWidget(HelpButton("h_check_input"))
        cl.addStretch(1)
        # notice: an existing prediction fits the input and the settings
        self.notice = QFrame()
        self.notice.setStyleSheet("QFrame { background: #fff3cd; border: 1px solid #e0a800; border-radius: 4px; }"
                                  "QLabel { border: none; color: #6b4e00; font-weight: bold; }")
        nl = QHBoxLayout(self.notice)
        self.notice_label = QLabel()
        self.notice_label.setWordWrap(True)
        nl.addWidget(self.notice_label, 1)
        nb = QPushButton()
        self.i18n.add(nb.setText, "notice_close")
        nb.clicked.connect(self.close_notice)
        nl.addWidget(nb)
        self.notice.setVisible(False)
        # hint for a quick start, below the input directory
        self.quick = QLabel()
        self.quick.setWordWrap(True)
        self.quick.setStyleSheet("QLabel { background: #e3f1e3; color: #1d5e20; font-weight: bold; "
                                 "border: 1px solid #8cc28c; border-radius: 4px; padding: 4px 6px; }")
        self.i18n.add(self.quick.setText, "quick_start")
        self.add_group("global", [
            g,
            mark(Param(self, "input", "input", "dir", required=True), "main", "prevalence"),
            self.quick,
            self.outdir_param("main", "prevalence"),
            cw,
            mark(Param(self, "task", "task", "tasks", "all"), "main"),
            self.notice,
            self.thread_param("main", "prevalence"),
        ])
        self.add_group("g_predict", [
            mark(Param(self, "seq", "seq", "choice", "prot", choices=[("prot", "seq_prot"), ("nucl", "seq_nucl")]), "main"),
            mark(Param(self, "keep_CDS_nucl", "keep_CDS_nucl", "bool", False), "main"),
        ])
        self.add_group("g_hmmsearch", [
            mark(Param(self, "hmm", "hmm", "hmm", "bac120"), "main", "prevalence"),
            mark(Param(self, "evalue", "evalue", "efloat", "1e-10"), "main", "prevalence"),
        ])
        self.add_group("g_refine", [
            mark(Param(self, "genome_cutoff", "genome_cutoff", "float", 0.8, minimum=0, maximum=1), "main", "prevalence"),
            mark(Param(self, "gene_cutoff", "gene_cutoff", "float", 0.8, minimum=0, maximum=1), "main", "prevalence"),
        ])
        self.add_group("g_alignment", [
            mark(Param(self, "trim", "trim", "choice", "strict", choices=[("gappyout", None), ("strict", None), ("strictplus", None), ("nogaps", None)]), "main"),
        ])
        self.add_group("g_tree", [
            mark(Param(self, "tree", "tree", "choice", "sm", choices=[("sm", "tree_sm"), ("st", "tree_st"), ("cs", "tree_cs")]), "main"),
            mark(Param(self, "tree_app", "tree_app", "choice", "fasttree", choices=[("fasttree", None), ("iqtree", None)]), "main"),
        ])
        self.add_group("g_prev_main", [
            Param(self, "prev_after", None, "bool", False, label="prev_after"),
            mark(Param(self, "prev_force", "force", "bool", False, label="prev_force"), "prevalence"),
        ])

    def input_dir(self):
        return self.v("input")

    def tasks(self):
        return self.v("task")

    def update_state(self):
        prev_only = self.rb_prev.isChecked()
        self.quick.setVisible(not prev_only)
        tasks = self.tasks()
        uses = set()
        for t in tasks:
            uses |= self.TASK_USES[t]
        for key in ("seq", "keep_CDS_nucl", "hmm", "evalue", "genome_cutoff", "gene_cutoff", "trim", "tree", "tree_app", "thread"):
            on = key in uses
            if prev_only:
                on = key in ("hmm", "evalue", "genome_cutoff", "gene_cutoff", "thread")
            self.p(key).set_active(on)
        self.p("task").set_active(not prev_only)
        if not prev_only and self.v("seq") == "nucl":
            self.p("keep_CDS_nucl").set_active(False)
        if not prev_only and self.v("tree") == "cs" and self.p("gene_cutoff").active:
            self.p("gene_cutoff").set_active(False)
        self.p("prev_after").set_active(not prev_only and tasks != ["predict"])
        self.p("prev_force").set_active(prev_only or self.v("prev_after"))

    # ---- notice about an existing prediction ----
    def after_refresh(self):
        boxes = self.p("task").boxes
        if self.rb_prev.isChecked():
            self.notice.setVisible(False)
            return
        match, info = prediction_match(self.v("input"), self.v("outdir"), self.v("seq"), self.v("keep_CDS_nucl"))
        key = (self.v("input"), self.v("outdir"), self.v("seq"), self.v("keep_CDS_nucl"), info)
        if match and key != self._notice_key:
            self.notice_label.setText(tr("notice_pred", info=info))
            self.notice.setVisible(True)
            if boxes[0].isChecked():
                prm = self.p("task")
                prm._tsync = True
                boxes[0].setChecked(False)
                if not any(b.isChecked() for b in boxes):
                    boxes[1].setChecked(True)
                prm._tsync = False
                prm.update_all_label()
                self._auto_unchecked = True
                boxes[0].setEnabled(False)
                self.changed()
            boxes[0].setEnabled(False)
        else:
            self.notice.setVisible(False)
            boxes[0].setEnabled(True)
            if not match and self._auto_unchecked:
                self._auto_unchecked = False
                if not boxes[0].isChecked():
                    boxes[0].setChecked(True)   # the task box fills the gap to the other ticked tasks

    def prediction_exists(self):
        d = self.v("input")
        return bool(d) and prediction_match(d, self.v("outdir"), self.v("seq"), self.v("keep_CDS_nucl"))[0]

    def close_notice(self):
        self._notice_key = (self.v("input"), self.v("outdir"), self.v("seq"), self.v("keep_CDS_nucl"),
                            prediction_match(self.v("input"), self.v("outdir"), self.v("seq"), self.v("keep_CDS_nucl"))[1])
        self.notice.setVisible(False)
        self.p("task").boxes[0].setEnabled(True)
        self._auto_unchecked = False

    # ---- commands ----
    def prevalence_cmd(self, full):
        """Gene_Prevelence.pl uses the HMM and E-value of the existing HMM search by default, so '-hmm' and
        '-evalue' are given only if they differ from the EasyCGTree defaults (and are used by the run)."""
        cfg = self.win.config
        argv = [str(cfg.script("prevalence"))]
        if self.v("input"):
            argv += ["-input", self.v("input")]
        argv += self.outdir_args(full)
        if self.rb_prev.isChecked() or "hmmsearch" in self.tasks():
            for key in ("hmm", "evalue"):
                prm = self.p(key)
                if prm.value() and not prm.is_default(prm.value()):
                    argv += ["-" + key, prm.value()]
        for key in ("genome_cutoff", "gene_cutoff"):
            prm = self.p(key)
            if full or not prm.is_default(prm.value()):
                argv += ["-" + key, prm.fmt(prm.value())]
        th = self.p("thread")
        if full or not th.is_default(th.value()):
            argv += ["-thread", str(th.value())]
        if self.v("prev_force"):
            argv.append("-force")
        return argv

    def commands(self, full):
        if self.rb_prev.isChecked():
            return [self.prevalence_cmd(full)]
        keys = ("seq", "keep_CDS_nucl", "hmm", "evalue", "genome_cutoff", "gene_cutoff", "trim", "tree", "tree_app", "thread")
        base = [str(self.win.config.script("main"))] + self.p("input").args(full) + self.outdir_args(full)
        tasks = self.tasks()
        if tasks == TASKS or not tasks:
            argv = base + self.p("task").args(full)
            for key in keys:
                argv += self.p(key).args(full)
            cmds = [argv]
        else:
            # Several tasks: one EasyCGTree command per task, run one after the other
            # (each with the options used by that task).
            cmds = []
            for t in tasks:
                uses = set(self.TASK_USES[t])
                if t == "hmmsearch" and "predict" not in tasks and not self.prediction_exists():
                    uses |= self.TASK_USES["predict"]   # EasyCGTree predicts the CDS first
                argv = base + ["-task", t]
                for key in keys:
                    if key in uses:
                        argv += self.p(key).args(full)
                cmds.append(argv)
        if self.v("prev_after") and self.tasks() != ["predict"]:
            cmds.append(self.prevalence_cmd(full))
        return cmds

    def check_input(self):
        d = self.v("input")
        if not d:
            QMessageBox.information(self, tr("info"), tr("msg_set_input"))
            return
        perl = self.win.require_perl()
        if not perl:
            return
        argv = [perl, str(self.win.config.script("main")), "-input", d, "-dry_run"]
        if self.v("seq") == "nucl":
            argv += ["-seq", "nucl"]
        self.check_btn.setEnabled(False)
        self._cap = Capture()
        self._cap.done.connect(self.show_dry_run)
        self._cap.start(argv, self.win.config.run_env())

    def show_dry_run(self, rc, out):
        self.check_btn.setEnabled(True)
        m = re.search(r"^EasyCGTree_DRY_RUN_JSON=(.*)$", out, re.M)
        if rc == 0 and m:
            try:
                DryRunDialog(json.loads(m.group(1)), self).exec()
                return
            except ValueError:
                pass
        err = "\n".join(l for l in out.splitlines() if l.strip() and not l.startswith("#####") and "Reading command line" not in l)
        QMessageBox.warning(self, tr("warning"), tr("msg_dryrun_failed") + "\n\n" + err[-3000:])

    def set_extra_state(self, st):
        if st.get("mode") == "prevalence":
            self.rb_prev.setChecked(True)
        else:
            self.rb_pipe.setChecked(True)

    def get_state(self):
        st = super().get_state()
        st["mode"] = "prevalence" if self.rb_prev.isChecked() else "pipeline"
        return st

    def compute_input_summary(self, d):
        cfg = self.win.config
        return input_summary_text(d, cfg.perl_path(), cfg.script("main"), cfg.run_env())

    def status_text(self, mode):
        return pipeline_status(self.v("input"), self.v("outdir"))


class SnpTab(BaseTab):
    name = "snp"
    ECG_DEFAULTS = {"hmm": "bac120", "evalue": "1e-10", "genome_cutoff": 0.8, "gene_cutoff": 0.8, "trim": "strict"}

    def build(self):
        self._filled_key = None
        self.add_group("global", [
            mark(Param(self, "input", "input", "dir", required=True, help="h_snp_input"), "snp", "specific", "prevalence"),
            self.outdir_param("snp", "specific", "prevalence"),
            self.thread_param("snp", "specific", "prevalence"),
        ])
        self.ecg_source = QLabel()
        self.ecg_source.setStyleSheet("color: #2a6f2a")
        self.add_group("g_ecg", [
            self.ecg_source,
            mark(Param(self, "hmm", "hmm", "hmm", "bac120"), "snp", "specific", "prevalence"),
            mark(Param(self, "evalue", "evalue", "efloat", "1e-10"), "snp"),
            mark(Param(self, "genome_cutoff", "genome_cutoff", "float", 0.8, minimum=0, maximum=1), "snp"),
            mark(Param(self, "gene_cutoff", "gene_cutoff", "float", 0.8, minimum=0, maximum=1), "snp"),
            mark(Param(self, "trim", "trim", "choice", "strict", choices=[("gappyout", None), ("strict", None), ("strictplus", None), ("nogaps", None)]), "snp"),
        ], help_key="h_g_ecg")
        self.auto1 = QLabel()          # step 1 decided automatically (another tree is annotated)
        self.auto1.setWordWrap(True)
        self.auto1.setStyleSheet("color: #1d5e20; font-weight: bold")
        self.auto1.setVisible(False)
        self.add_group("g_snp1", [
            self.auto1,
            mark(Param(self, "ignore", "ignore", "strains"), "snp"),
            mark(Param(self, "max_missing", "max_missing", "float", 0.0, minimum=0, maximum=1), "snp"),
            mark(Param(self, "aln", "aln", "choice", "trimmed", choices=[("trimmed", "aln_trimmed"), ("original", "aln_original")]), "snp"),
            mark(Param(self, "ref", "ref", "strains", multi=False), "snp"),
        ], checkable=True)
        self.add_group("g_snp2", [
            mark(Param(self, "tree_app", "tree_app", "choice", "iqtree", choices=[("iqtree", None), ("fasttree", None)], help="h_snp_tree_app"), "snp", "specific"),
            mark(Param(self, "asc", "asc", "choice", "fconst", choices=[("fconst", None), ("lewis", None)]), "snp"),
        ], checkable=True)
        self.add_group("g_snp3", [
            mark(Param(self, "mode", "mode", "choice", "exclusive", choices=[("exclusive", None), ("strict", None)]), "specific"),
            mark(Param(self, "outgroup", "outgroup", "strains"), "specific"),
            mark(Param(self, "other_tree", "tree", "file", file_filter="Trees (*.tree *.treefile *.contree *.nwk *.newick *.tre);;* (*)"), "specific"),
        ], checkable=True)
        self.add_group("g_prev_snp", [
            Param(self, "prev_after", None, "bool", False, label="prev_after"),
        ])
        self.p("ignore").keep_free = 2    # SNPs need at least two strains
        # notice: what will be run again, and why
        self.plan_box = QFrame()
        self.plan_box.setStyleSheet("QFrame { background: #fff3cd; border: 1px solid #e0a800; border-radius: 4px; }"
                                    "QLabel { border: none; color: #6b4e00; }")
        pl = QVBoxLayout(self.plan_box)
        self.plan_label = QLabel()
        self.plan_label.setWordWrap(True)
        pl.addWidget(self.plan_label)
        self.plan_box.setVisible(False)
        self.form_layout.addWidget(self.plan_box)
        self._auto = False            # step 1 automatic
        self._saved = {}              # check states of steps 1 and 2 before the automatic mode

    def input_dir(self):
        return self.v("input")

    def after_refresh(self):
        """Fill the EasyCGTree settings from the record of the previous run (or the defaults)."""
        d = self.v("input")
        rec = read_record(io_dirs(d, self.v("outdir"))[2]) if d else {}
        vals = dict(self.ECG_DEFAULTS)
        if rec.get("hmmsearch_done"):
            vals["hmm"] = rec.get("hmm") or vals["hmm"]
            vals["evalue"] = rec.get("evalue") or vals["evalue"]
        if rec.get("refine_done"):
            for k in ("genome_cutoff", "gene_cutoff"):
                try:
                    vals[k] = float(rec.get(k))
                except (TypeError, ValueError):
                    pass
        if rec.get("alignment_done") and rec.get("trim"):
            vals["trim"] = rec["trim"]
        src = "ecg_from_record" if rec.get("hmmsearch_done") else "ecg_from_defaults"
        key = (d, self.v("outdir"), tuple(sorted((k, str(v)) for k, v in vals.items())))
        self.ecg_source.setText(tr(src))
        if key != self._filled_key:
            self._filled_key = key
            self._building = True
            for k, v in vals.items():
                self.p(k).autofilled = True
                self.p(k).default = v
                self.p(k).set_value(v)
            self._building = False
            self.changed()

    # ---- steps ----
    def other_mode(self):
        """Another tree is annotated: step 2 is not needed and step 1 is run only if needed."""
        return self.groups["g_snp3"].isChecked() and bool(self.v("other_tree"))

    def ecg_values(self):
        return {k: self.v(k) for k in ("hmm", "evalue", "genome_cutoff", "gene_cutoff", "trim")}

    def snp_values(self):
        return {k: self.v(k) for k in ("ignore", "max_missing", "aln", "ref")}

    def plan(self):
        d = self.v("input")
        if not d:
            return None
        try:
            return snp_plan(d, self.v("outdir"), self.ecg_values(), self.snp_values())
        except OSError:
            return None

    def step_on(self, n):
        g = self.groups["g_snp%d" % n]
        if n == 1 and self._auto:
            pl = self.plan()
            return pl is None or bool(pl["snp_why"])
        return g.isChecked()

    def set_auto(self, on):
        g1, g2 = self.groups["g_snp1"], self.groups["g_snp2"]
        if on == self._auto:
            return
        self._auto = on
        for g in (g1, g2):
            g.blockSignals(True)
        if on:
            self._saved = {"g1": g1.isChecked(), "g2": g2.isChecked()}
            g1.setCheckable(False)          # the settings stay editable; the step is run if needed
            g2.setChecked(False)
            g2.setEnabled(False)
        else:
            g1.setCheckable(True)
            g1.setChecked(self._saved.get("g1", True))
            g2.setEnabled(True)
            g2.setChecked(self._saved.get("g2", True) and g1.isChecked())
        for g in (g1, g2):
            g.blockSignals(False)
        self.auto1.setVisible(on)

    def get_state(self):
        st = super().get_state()
        if self._auto:     # the check states the user chose before the automatic mode
            st["group:g_snp1"] = self._saved.get("g1", True)
            st["group:g_snp2"] = self._saved.get("g2", True)
        return st

    def set_state(self, st):
        self.set_auto(False)
        super().set_state(st)

    def strain_names(self, key):
        if key == "outgroup" and self.other_mode():
            return newick_leaves(self.v("other_tree"))
        return None

    def update_state(self):
        self.set_auto(self.other_mode())
        s1 = self._auto or self.groups["g_snp1"].isChecked()
        g2 = self.groups["g_snp2"]
        if not self._auto:
            if not s1 and g2.isChecked():
                g2.blockSignals(True)
                g2.setChecked(False)
                g2.blockSignals(False)
            g2.setEnabled(s1)
        s2 = g2.isChecked()
        s3 = self.groups["g_snp3"].isChecked()
        for k in ("ignore", "max_missing", "aln", "ref"):
            self.p(k).set_active(s1)
        self.p("asc").set_active(s2 and self.v("tree_app") == "iqtree")
        self.p("tree_app").set_active(s2 or (s3 and not self._auto))
        for k in ("mode", "outgroup", "other_tree"):
            self.p(k).set_active(s3)
        for k in ("evalue", "genome_cutoff", "gene_cutoff", "trim"):
            self.p(k).set_active(s1 or s3)
        self.update_plan()

    def update_plan(self):
        """Notice: which steps will be run again and why; names of the given tree."""
        pl = self.plan()
        lines = []
        s1 = self.step_on(1)
        if pl is not None:
            if self._auto:
                if pl["snp_why"]:
                    self.auto1.setText(tr("auto1_run", why="; ".join(pl["snp_why"])))
                else:
                    self.auto1.setText(tr("auto1_skip"))
            if s1 or (self.groups["g_snp3"].isChecked() and pl["snp_fas"] is None):
                if pl["ecg_start"]:
                    lines.append(tr("plan_ecg", t=pl["ecg_start"], why="; ".join(pl["ecg_why"])))
                    if pl["ecg_start"] != "predict":
                        lines.append(tr("plan_ecg_tem"))
                elif s1 and not self._auto and pl["snp_why"] and pl["snp_why"] != [tr("why_snp_none")]:
                    lines.append(tr("plan_snp_diff", why="; ".join(pl["snp_why"])))
            if not s1 and self.groups["g_snp3"].isChecked() and pl["snp_fas"] is None and not self._auto:
                lines.append(tr("plan_no_snp"))
            # the strains of the given tree must be those of the SNP alignment
            if self.other_mode():
                leaves = newick_leaves(self.v("other_tree"))
                if pl["snp_fas"] is not None:
                    taxa = fasta_ids(pl["snp_fas"])
                else:
                    removed = screened_out(self.v("input"), self.v("outdir"))
                    taxa = [n for n in strains_of(self.v("input"), self.v("outdir")) if n not in removed]
                if leaves and taxa:
                    only_t = sorted(set(leaves) - set(taxa))
                    only_a = sorted(set(taxa) - set(leaves))
                    if only_t or only_a:
                        lines.append(tr("plan_tree_names", t=", ".join(only_t) or "-", a=", ".join(only_a) or "-"))
                elif not leaves:
                    lines.append(tr("msg_no_tree_names"))
            rerun = (pl["ecg_start"] and (s1 or self.groups["g_snp3"].isChecked())) or \
                (s1 and pl["snp_fas"] is None and pl["snp_why"] != [tr("why_snp_none")])
            if rerun:
                lines.append(tr("plan_backup"))
        self.plan_label.setText("\n".join(lines))
        self.plan_box.setVisible(bool(lines))

    def validate(self):
        if not (self.groups["g_snp1"].isChecked() or self.groups["g_snp3"].isChecked() or self.v("prev_after")):
            return tr("need_step")
        if self.other_mode():
            ot = self.v("other_tree")
            if not os.path.isfile(ot) or not newick_leaves(ot):
                return tr("msg_no_tree_names")
            pl = self.plan()
            if pl and pl["snp_fas"] is not None:
                diff = set(newick_leaves(ot)) ^ set(fasta_ids(pl["snp_fas"]))
                if diff:
                    return tr("msg_tree_names", n=", ".join(sorted(diff)[:10]))
        if self.step_on(1) and self.v("ignore"):
            removed = screened_out(self.v("input"), self.v("outdir"))
            names = [n for n in strains_of(self.v("input"), self.v("outdir")) if n not in removed]
            ign = {x.strip() for x in self.v("ignore").split(",") if x.strip()}
            if names and len([n for n in names if n not in ign]) < 2:
                return tr("msg_ignore_too_many")
        return super().validate()

    def commands(self, full):
        cfg = self.win.config
        cmds = []
        s1 = self.step_on(1)
        s2 = self.groups["g_snp2"].isChecked() and not self._auto
        s3 = self.groups["g_snp3"].isChecked()
        if s1:
            argv = [str(cfg.script("snp"))] + self.p("input").args(full) + self.outdir_args(full)
            for k in ("ignore", "max_missing", "aln", "ref"):
                argv += self.p(k).args(full)
            if s2:
                argv += self.p("tree_app").args(full) + self.p("asc").args(full)
            else:
                argv.append("-no_tree")
            for k in ("hmm", "evalue", "genome_cutoff", "gene_cutoff", "trim", "thread"):
                argv += self.p(k).args(full)
            cmds.append(argv)
        if s3:
            argv = [str(cfg.script("specific"))] + self.p("input").args(full) + self.outdir_args(full)
            # with another tree, '-tree_app' (which SNP tree) does not matter
            for k in (("mode", "outgroup", "other_tree", "hmm", "thread") if self._auto else
                      ("tree_app", "mode", "outgroup", "other_tree", "hmm", "thread")):
                argv += self.p(k).args(full)
            cmds.append(argv)
        if self.v("prev_after"):
            argv = [str(cfg.script("prevalence"))] + self.p("input").args(full) + self.outdir_args(full)
            for k in ("hmm", "thread"):
                argv += self.p(k).args(full)
            cmds.append(argv)
        return cmds

    def compute_input_summary(self, d):
        cfg = self.win.config
        return input_summary_text(d, cfg.perl_path(), cfg.script("main"), cfg.run_env())

    def status_text(self, mode):
        return pipeline_status(self.v("input"), self.v("outdir"), snp=True)


class HmmTab(BaseTab):
    name = "hmm"

    def build(self):
        self._auto_name = ""
        self.add_group("global", [
            mark(Param(self, "gc", "gc", "dir", required=True), "hmm"),
            self.outdir_param("hmm", help="h_hmm_outdir"),
            self.thread_param("hmm", help="h_thread_hmm"),
        ])
        self.add_group("g_hmm_aln", [
            mark(Param(self, "aln", "aln", "bool", False, label="hmm_aln", help="h_hmm_aln"), "hmm"),
            mark(Param(self, "super5", "super5", "int", 1000, minimum=2, maximum=10000000), "hmm"),
        ])
        self.add_group("g_hmm_out", [
            mark(Param(self, "name", "name", "text", "", label="hmm_name", help="h_hmm_name"), "hmm"),
            mark(Param(self, "force", "force", "bool", False, label="hmm_force", help="h_hmm_force"), "hmm"),
        ])

    def browser_places(self):
        return [("input", "place_families", self.input_dir), ("output", "place_output", self.output_dir),
                ("hmmdir", "hmm_dir_button", lambda: str(self.win.config.hmm_dir()))]

    def input_dir(self):
        return self.v("gc")

    def output_dir(self):
        """Default output directory: the HMM directory of EasyCGTree."""
        return self.v("outdir") or str(self.win.config.hmm_dir())

    def outdir_args(self, full):
        # BuildHMM.pl itself would write next to the gene family directory, so '-outdir' is always given
        return ["-outdir", self.output_dir()] if self.v("gc") else []

    def update_outdir_hint(self):
        self.p("outdir").field.edit.setPlaceholderText(tr("default_outdir", d=str(self.win.config.hmm_dir())))

    def update_state(self):
        self.p("super5").set_active(not self.v("aln"))
        d = self.v("gc")
        auto = safe_hmm_name(os.path.basename(d.rstrip("/\\"))) if d else ""
        name = self.p("name")
        # the name follows the gene family directory until the user types another one
        if not self._building and (name.value() in ("", self._auto_name)) and auto != self._auto_name:
            self._building = True
            name.set_value(auto)
            self._building = False
        self._auto_name = auto
        name.field.setPlaceholderText(auto)

    def db_name(self):
        return self.v("name") or self._auto_name

    def validate(self):
        e = super().validate()
        if e:
            return e
        n = self.db_name()
        if n and not re.fullmatch(r"[A-Za-z0-9_.\-]+", re.sub(r"\.hmm$", "", n, flags=re.I)):
            return tr("msg_hmm_name_bad", n=n)
        if n and not self.v("force") and (self.win.config.hmm_dir() / (re.sub(r"\.hmm$", "", n, flags=re.I) + ".hmm")).exists():
            return tr("msg_hmm_exists", n=n)
        return None

    def commands(self, full):
        argv = [str(self.win.config.script("hmm"))] + self.p("gc").args(full) + self.outdir_args(full)
        if self.v("gc") and self.db_name():
            argv += ["-name", self.db_name()]      # always given, so the database name is visible
        for k in ("aln", "super5", "force", "thread"):
            argv += self.p(k).args(full)
        return [argv]

    def compute_input_summary(self, d):
        return family_summary_text(d)

    def status_text(self, mode):
        cfg = self.win.config
        if mode == "output":
            d = self.v("gc")
            if not d:
                return tr("st_no_input")
            f = Path(self.output_dir()) / (Path(d.rstrip("/\\")).name + ".BuildHMM_summary.txt")
            if not f.is_file():
                return tr("st_out", d=self.output_dir()) + "\n" + tr("hmm_no_summary")
            with open(f, encoding="utf-8", errors="replace") as fh:
                return tr("st_out", d=self.output_dir()) + "\n" + fh.read()
        sets = hmm_sets(cfg.hmm_dir())
        hs = self.win.hmmsearch_version()
        lines = [tr("hmm_sets") + ": " + str(cfg.hmm_dir())]
        for name, n, fmt in sets:
            note = ""
            if hs.startswith("3.0") and hmm_format_version(fmt) not in ("", "a", "b"):
                note = "  <- " + tr("st_hmm_incompatible", v=hs)
            lines.append("  %-14s %4d   %s%s" % (name, n, fmt, note))
        return "\n".join(lines)



class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.config = Config()
        LANG["code"] = self.config.language if self.config.language in ("en", "zh") else "en"
        self.i18n = I18n()
        self.runner = Runner()
        self.runner.output.connect(self.on_output)
        self.runner.started_cmd.connect(self.on_started)
        self.runner.finished.connect(self.on_finished)
        self.job_tab = None
        self._hs_version = None
        self.i18n.add(self.setWindowTitle, "app_title")
        self.tabs = QTabWidget()
        self.main_tab = MainTab(self)
        self.snp_tab = SnpTab(self)
        self.hmm_tab = HmmTab(self)
        self.tabs_list = [self.main_tab, self.snp_tab, self.hmm_tab]
        for t, key in zip(self.tabs_list, ("tab_main", "tab_snp", "tab_hmm")):
            idx = self.tabs.addTab(t, tr(key))
            self.i18n.add(lambda text, i=idx: self.tabs.setTabText(i, text), key)
        self.setCentralWidget(self.tabs)
        self._tab_index = 0
        self.tabs.currentChanged.connect(self.tab_switched)
        self.build_menus()
        # The analysis parameters start with their defaults every time (only the environment settings are
        # kept); parameter sets can be saved and loaded with File - Save/Load parameters.
        self.resize(1300, 900)
        if self.config.geometry:
            try:
                self.restoreGeometry(QByteArray.fromBase64(self.config.geometry.encode("ascii")))
            except (TypeError, ValueError, AttributeError):
                pass
        self.statusBar()
        self.update_statusbar()
        if not self.config.home_ok():
            QTimer.singleShot(300, lambda: QMessageBox.warning(self, tr("warning"), tr("msg_no_home")))

    def build_menus(self):
        mb = self.menuBar()
        m_file = mb.addMenu("")
        self.i18n.add(m_file.setTitle, "menu_file")
        a = m_file.addAction("")
        self.i18n.add(a.setText, "menu_load")
        a.triggered.connect(self.load_params)
        a = m_file.addAction("")
        self.i18n.add(a.setText, "menu_save")
        a.triggered.connect(self.save_params)
        m_file.addSeparator()
        a = m_file.addAction("")
        self.i18n.add(a.setText, "menu_exit")
        a.triggered.connect(self.close)
        m_set = mb.addMenu("")
        self.i18n.add(m_set.setTitle, "menu_settings")
        a = m_set.addAction("")
        self.i18n.add(a.setText, "menu_paths")
        a.triggered.connect(lambda: SettingsDialog(self, self).exec())
        a = m_set.addAction("")
        self.i18n.add(a.setText, "menu_treeopts")
        a.triggered.connect(lambda: TreeOptionsDialog(self, self).exec())
        m_lang = m_set.addMenu("")
        self.i18n.add(m_lang.setTitle, "menu_language")
        ag = QActionGroup(self)
        for code, label in (("en", "English"), ("zh", "中文")):
            act = m_lang.addAction(label)
            act.setCheckable(True)
            act.setChecked(LANG["code"] == code)
            ag.addAction(act)
            act.triggered.connect(lambda checked=False, c=code: self.set_language(c))
        m_help = mb.addMenu("")
        self.i18n.add(m_help.setTitle, "menu_help")
        a = m_help.addAction("")
        self.i18n.add(a.setText, "menu_check")
        a.triggered.connect(lambda: CheckDialog(self, self).exec())
        m_hmm = m_help.addMenu("")
        self.i18n.add(m_hmm.setTitle, "menu_hmm_online")
        for site, url in REPO_HMM.items():
            act = m_hmm.addAction(site)
            act.triggered.connect(lambda checked=False, u=url: QDesktopServices.openUrl(QUrl(u)))
        a = m_hmm.addAction("")
        self.i18n.add(a.setText, "menu_hmm_folder")
        a.triggered.connect(lambda: open_file(str(self.config.hmm_dir()), self))
        m_help.addSeparator()
        a = m_help.addAction("")
        self.i18n.add(a.setText, "menu_about")
        a.triggered.connect(self.about)

    # ---- settings ----
    def tab_switched(self, i):
        """All tabs share the panel sizes and the page shown on the right (Status or Files)."""
        old = self.tabs.widget(self._tab_index)
        new = self.tabs.widget(i)
        self._tab_index = i
        if old is not None and old is not new:
            old.save_splitters()
            new.restore_splitters()
            new.show_right_page(self.config.right_page if self.config.right_page in (0, 1) else 0)

    def set_language(self, code):
        LANG["code"] = code
        self.config.language = code
        self.i18n.apply()
        for t in self.tabs_list:
            t.retranslate_state()
            t.align_labels()
            t.changed()
        self.update_statusbar()

    def config_changed(self):
        self._hs_version = None
        for t in self.tabs_list:
            t.reload_hmm_lists()
            t.changed()
        self.update_statusbar()
        self.config.save()

    def hmmsearch_version(self):
        if self._hs_version is None:
            p = self.config.program_path("hmmsearch")
            m = None
            if p.is_file():
                rc, out = capture([str(p), "-h"], timeout=20)
                m = re.search(r"HMMER (\d+\.\d+(\.\d+)?)", out)
            self._hs_version = m.group(1) if m else ""
        return self._hs_version

    def update_statusbar(self):
        perl = self.config.perl_path()
        home = self.config.home if self.config.home_ok() else "?"
        self.statusBar().showMessage("EasyCGTree: %s    Perl: %s" % (home, perl or tr("perl_not_found")))

    def require_perl(self):
        if not self.config.home_ok():
            QMessageBox.warning(self, tr("warning"), tr("msg_no_home"))
            return ""
        perl = self.config.perl_path()
        if not perl:
            QMessageBox.warning(self, tr("warning"), tr("msg_no_perl"))
        return perl

    # ---- jobs ----
    def start_job(self, tab, cmds):
        if self.runner.running():
            QMessageBox.information(self, tr("info"), tr("st_running"))
            return
        perl = self.require_perl()
        if not perl or not cmds:
            return
        for argv in cmds:
            if not Path(argv[0]).is_file():
                QMessageBox.warning(self, tr("warning"), tr("msg_no_script", s=os.path.basename(argv[0]), d=os.path.dirname(argv[0])))
                return
        self.job_tab = tab
        for t in self.tabs_list:
            t.run_btn.setEnabled(False)
        tab.stop_btn.setEnabled(True)
        tab.set_run_state("st_running", "running")
        tab.log.feed("\n" if tab.log.toPlainText() else "")
        self.runner.start([[perl] + argv for argv in cmds], self.config.run_env(), tab.job_cwd())

    def stop(self):
        if self.runner.running():
            self.runner.stop()

    def on_started(self, cmd):
        if self.job_tab:
            self.job_tab.log.feed("\n$ %s\n" % cmd)

    def on_output(self, text):
        if self.job_tab:
            self.job_tab.log.feed(text)

    def on_finished(self, rc, stopped):
        tab = self.job_tab
        for t in self.tabs_list:
            t.run_btn.setEnabled(True)
            t.stop_btn.setEnabled(False)
        if tab:
            key, kind, kw = (("st_stopped", "stopped", {}) if stopped else
                             (("st_finished", "ok", {}) if rc == 0 else ("st_failed", "failed", {"rc": rc})))
            tab.set_run_state(key, kind, **kw)
            tab.log.feed("\n[%s]\n" % tr(key, **kw))
        for t in self.tabs_list:
            t.reload_hmm_lists()
            t.refresh_status()
        self.job_tab = None

    # ---- presets ----
    def all_states(self):
        return {t.name: t.get_state() for t in self.tabs_list}

    def save_params(self):
        p, _ = QFileDialog.getSaveFileName(self, tr("menu_save"), str(Path.home() / "EasyCGTree_parameters.json"), "JSON (*.json)")
        if not p:
            return
        with open(p, "w", encoding="utf-8") as fh:
            json.dump({"EasyCGTree_GUI": GUI_VERSION, "tabs": self.all_states()}, fh, indent=1, ensure_ascii=False)
        self.statusBar().showMessage(tr("msg_saved", f=p), 5000)

    def load_params(self):
        p, _ = QFileDialog.getOpenFileName(self, tr("menu_load"), str(Path.home()), "JSON (*.json);;* (*)")
        if not p:
            return
        try:
            with open(p, encoding="utf-8") as fh:
                data = json.load(fh)
            for t in self.tabs_list:
                t.set_state(data.get("tabs", {}).get(t.name, {}))
        except (OSError, ValueError, AttributeError) as e:
            QMessageBox.warning(self, tr("warning"), tr("msg_load_failed", e=e))

    def about(self):
        a = ABOUT
        text = ("<h3>%s %s</h3><p>%s</p><p>%s</p><p><b>%s</b></p><p>%s</p><p><a href='%s'>%s</a><br>"
                "GitHub: <a href='%s'>%s</a><br>Gitee: <a href='%s'>%s</a></p>%s"
                % (a["title"], GUI_VERSION, tr("about_text", v="5.0"), a["authors"], "Citation", a["citation"],
                   a["doi"], a["doi"], a["url"], a["url"], a["url2"], a["url2"],
                   ("<p>%s</p>" % a["contact"]) if a["contact"] else ""))
        QMessageBox.about(self, tr("dlg_about"), text)

    def closeEvent(self, event):
        if self.runner.running():
            if QMessageBox.question(self, tr("warning"), tr("msg_quit_running")) != QMessageBox.StandardButton.Yes:
                event.ignore()
                return
            self.runner.stop()
        self.config.params = {}       # analysis parameters are not kept
        self.tabs.currentWidget().save_splitters()
        try:
            self.config.geometry = bytes(self.saveGeometry().toBase64()).decode("ascii")
        except (TypeError, ValueError, AttributeError):
            pass
        self.config.save()
        event.accept()


# ======================================================================================================
# Appearance: a light, restrained theme (navy accent) on Qt's Fusion style, the same on all systems
# ======================================================================================================

C = {"bg": "#eceff3", "panel": "#ffffff", "border": "#d3d9e1", "field": "#c3ccd8", "text": "#1e2a36",
     "muted": "#5f6b7a", "accent": "#1f4e79", "accent_hover": "#2a6296", "accent_soft": "#e3ecf5",
     "code": "#f6f8fa", "ok": "#2f6b3a", "warn": "#8a5a00", "bad": "#9b2c2c"}

STYLE = """
QMainWindow, QDialog { background: %(bg)s; }
QWidget#paramForm { background: %(bg)s; }
QScrollArea#paramScroll { background: transparent; }
QTabWidget::pane { border: none; border-top: 1px solid %(border)s; }
QTabBar::tab { background: transparent; color: %(muted)s; padding: 7px 18px; margin-right: 2px;
               border: none; border-bottom: 3px solid transparent; font-weight: 600; }
QTabBar::tab:hover { color: %(text)s; }
QTabBar::tab:selected { color: %(accent)s; border-bottom: 3px solid %(accent)s; }
QGroupBox { background: %(panel)s; border: 1px solid %(border)s; border-radius: 6px;
            margin-top: 12px; padding: 10px 8px 8px 8px; }
QGroupBox { color: %(accent)s; font-weight: 600; }
QGroupBox QWidget { font-weight: normal; }
QGroupBox::title { subcontrol-origin: margin; subcontrol-position: top left; left: 10px; padding: 0 4px; }
QLabel#panelTitle { color: %(accent)s; font-weight: 600; padding: 4px 2px; }
QFrame#panel { background: %(panel)s; border: 1px solid %(border)s; border-radius: 6px; }
QToolButton#segment { border: none; border-bottom: 2px solid transparent; padding: 4px 12px;
                      color: %(muted)s; font-weight: 600; background: transparent; }
QToolButton#segment:hover { color: %(text)s; }
QToolButton#segment:checked { color: %(accent)s; border-bottom: 2px solid %(accent)s; }
QToolButton#place { border: 1px solid %(field)s; border-left: none; background: #f5f7fa; padding: 3px 12px;
                    color: %(text)s; }
QToolButton#place[pos="first"] { border-left: 1px solid %(field)s; border-top-left-radius: 4px;
                                 border-bottom-left-radius: 4px; }
QToolButton#place[pos="last"] { border-top-right-radius: 4px; border-bottom-right-radius: 4px; }
QToolButton#place:hover { background: %(accent_soft)s; }
QToolButton#place:checked { background: %(accent)s; color: white; border-color: %(accent)s; }
QPushButton { background: #f7f9fb; border: 1px solid #b9c3cf; border-radius: 4px; padding: 4px 14px;
              color: %(text)s; }
QPushButton:hover { background: %(accent_soft)s; border-color: %(accent)s; }
QPushButton:pressed { background: #d4e2f0; }
QPushButton:disabled { color: #9aa5b1; border-color: #d6dce3; background: #f3f5f7; }
QPushButton#runButton { background: %(accent)s; color: white; border: 1px solid %(accent)s;
                        font-weight: 600; padding: 5px 24px; }
QPushButton#runButton:hover { background: %(accent_hover)s; }
QPushButton#runButton:disabled { background: #9fb1c4; border-color: #9fb1c4; color: #eef2f6; }
QPushButton#stopButton { color: %(bad)s; border-color: %(bad)s; background: white; padding: 5px 18px; }
QPushButton#stopButton:hover { background: #f8e9e9; }
QPushButton#stopButton:disabled { color: #c9b3b3; border-color: #e3d3d3; background: #f7f5f5; }
QToolButton#helpButton { background: %(accent_soft)s; color: %(accent)s; border: 1px solid #c5d5e6;
                         border-radius: 9px; font-weight: 700; padding: 0px; }
QToolButton#helpButton:hover { background: #cddcee; border-color: %(accent)s; }
QLabel#stateLabel { color: %(muted)s; }
QLabel#stateLabel[state="running"] { color: %(accent)s; font-weight: 600; }
QLabel#stateLabel[state="ok"] { color: %(ok)s; font-weight: 600; }
QLabel#stateLabel[state="failed"] { color: %(bad)s; font-weight: 600; }
QLabel#stateLabel[state="stopped"] { color: %(warn)s; font-weight: 600; }
QLineEdit, QPlainTextEdit, QTreeView, QListWidget, QTableWidget {
    background: white; border: 1px solid %(field)s; border-radius: 4px;
    selection-background-color: #cfe0f1; selection-color: %(text)s; }
QLineEdit { padding: 2px 4px; }
QLineEdit:focus, QPlainTextEdit:focus { border: 1px solid %(accent)s; }
QPlainTextEdit#cmdView { background: %(code)s; }
QPlainTextEdit#logView { background: #fbfcfd; }
QHeaderView::section { background: #f1f4f7; color: %(muted)s; border: none;
                       border-bottom: 1px solid %(border)s; padding: 4px 6px; font-weight: 600; }
QSplitter::handle { background: transparent; }
QSplitter::handle:hover { background: #c9d6e3; }
QMenuBar { background: %(panel)s; border-bottom: 1px solid %(border)s; }
QMenuBar::item:selected { background: %(accent_soft)s; }
QMenu { background: %(panel)s; border: 1px solid %(border)s; }
QMenu::item:selected { background: %(accent_soft)s; color: %(text)s; }
QStatusBar { background: %(panel)s; border-top: 1px solid %(border)s; color: %(muted)s; }
QToolTip { background: #24313f; color: white; border: none; padding: 4px 6px; }
""" % C

ICON_PNG_B64 = (
    "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"
    "AElFTkSuQmCC")


def app_icon():
    pm = QPixmap()
    pm.loadFromData(QByteArray.fromBase64(ICON_PNG_B64.encode("ascii")), "PNG")
    return QIcon(pm)


def apply_theme(app):
    """Fusion style with a light palette: the same look on Windows, macOS and Linux, also when the
    system uses a dark theme (the parameter forms are designed for a light background)."""
    try:
        app.setStyle(QStyleFactory.create("Fusion"))
        pal = QPalette()
        role = QPalette.ColorRole
        for r, c in ((role.Window, C["bg"]), (role.WindowText, C["text"]), (role.Base, "#ffffff"),
                     (role.AlternateBase, "#f6f8fa"), (role.Text, C["text"]), (role.Button, "#f7f9fb"),
                     (role.ButtonText, C["text"]), (role.Highlight, C["accent"]),
                     (role.HighlightedText, "#ffffff"), (role.ToolTipBase, "#24313f"),
                     (role.ToolTipText, "#ffffff"), (role.Link, C["accent"]), (role.PlaceholderText, "#8a96a3")):
            pal.setColor(r, QColor(c))
        for r in (role.WindowText, role.Text, role.ButtonText):
            pal.setColor(QPalette.ColorGroup.Disabled, r, QColor("#9aa5b1"))
        app.setPalette(pal)
        app.setStyleSheet(STYLE)
        app.setWindowIcon(app_icon())
    except Exception:   # never fail to start because of the appearance
        pass


def main():
    app = QApplication(sys.argv)
    app.setApplicationName("EasyCGTree_GUI")
    apply_theme(app)
    w = MainWindow()
    w.show()
    return app.exec()


if __name__ == "__main__":
    sys.exit(main())
