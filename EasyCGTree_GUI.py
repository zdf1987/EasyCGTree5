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
        self.setAutoRaise(True)
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
            self.field.setRange(self.minimum if self.minimum is not None else 0, self.maximum if self.maximum is not None else 1000000)
            self.field.setValue(int(self.default))
            self.field.valueChanged.connect(self.tab.changed)
            self.widget = self.field
        elif k == "float":
            self.field = GDoubleSpinBox()
            self.field.setDecimals(self.decimals)
            self.field.setRange(self.minimum if self.minimum is not None else 0.0, self.maximum if self.maximum is not None else 1e9)
            self.field.setSingleStep(self.step or 0.05)
            self.field.setValue(float(self.default))
            self.field.valueChanged.connect(self.tab.changed)
            self.widget = self.field
        elif k in ("efloat", "text"):
            self.field = QLineEdit()
            self.field.setMinimumWidth(60)
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
            narrow_combo(self.field)
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
            self.all_label.setText(tr("tasks_all") if on else "")
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
        # one header row: title, Concise|Full, (space), Copy, Run, Stop; the command text below
        self.rb_concise = QToolButton()
        self.rb_full = QToolButton()
        cbg = QButtonGroup(self)
        cbg.setExclusive(True)
        for b, key, pos in ((self.rb_concise, "concise", "first"), (self.rb_full, "full", "last")):
            b.setObjectName("place")
            b.setProperty("pos", pos)
            b.setCheckable(True)
            self.i18n.add(b.setText, key)
            self.i18n.add(b.setToolTip, key + "_tip")
            cbg.addButton(b)
        (self.rb_full if win.config.full_cmd else self.rb_concise).setChecked(True)
        self.rb_full.toggled.connect(self.changed)
        chead.insertSpacing(1, 8)
        chead.insertWidget(2, self.rb_concise)
        chead.insertWidget(3, self.rb_full)
        self.copy_btn = QPushButton()
        self.i18n.add(self.copy_btn.setText, "copy")
        self.copy_btn.clicked.connect(self.copy_cmd)
        chead.addWidget(self.copy_btn)
        chead.addSpacing(10)
        self.run_btn = QPushButton()
        self.run_btn.setObjectName("runButton")
        self.i18n.add(self.run_btn.setText, "run")
        self.run_btn.clicked.connect(lambda: self.run())
        self.stop_btn = QPushButton()
        self.stop_btn.setObjectName("stopButton")
        self.i18n.add(self.stop_btn.setText, "stop")
        self.stop_btn.clicked.connect(self.win.stop)
        self.stop_btn.setEnabled(False)
        chead.addWidget(self.run_btn)
        chead.addWidget(self.stop_btn)
        self.cmd_view = QPlainTextEdit()
        self.cmd_view.setObjectName("cmdView")
        self.cmd_view.setReadOnly(True)
        self.cmd_view.setFont(fixed_font())
        cl.addWidget(self.cmd_view, 1)
        cg.setMinimumWidth(420)
        bottom.addWidget(cg)
        lg, lhead, ll = self.panel("log")
        self.state_label = QLabel()          # state of the run, next to the title of the log
        self.state_label.setObjectName("stateLabel")
        self._state = None
        lhead.insertSpacing(1, 10)
        lhead.insertWidget(2, self.state_label)
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
                if item.kind == "bool":
                    grid.addWidget(w, r, 0, 1, 2)
                else:
                    grid.addWidget(lab, r, 0)
                    grid.addWidget(w, r, 1)
                grid.addWidget(hb, r, 2)
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
    "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")


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
