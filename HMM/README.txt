EasyCGTree 5.0 - profile HMM database (PHD)
============================================

Put the profile HMM sets (.hmm files) that you need into this folder.

The HMM sets are NOT included in the program package, because they are large (about 220 MB together) and most
users need only a few of them. Download each set you need separately from the folder 'HMM' of the EasyCGTree5
repository, and save it in this folder without changing its name (e.g. 'bac120.hmm'):

    GitHub:  https://github.com/zdf1987/EasyCGTree5/tree/main/HMM
    Gitee:   https://gitee.com/zdf1987/EasyCGTree5/tree/main/HMM    (recommended in mainland China)

How to download one file: open the file in the folder 'HMM' of the repository, then click 'Download raw file'
(GitHub) or 'Download' (Gitee).

The default set of EasyCGTree is 'bac120': download 'bac120.hmm' to run EasyCGTree with its default settings.
If a set is missing, EasyCGTree and EasyCGTree_GUI stop with a message and point to this file.
A set can also be used from another folder by giving the path of the .hmm file ('-hmm /path/to/set.hmm').
New sets can be built from your own gene families with BuildHMM.pl (or the 'HMM building' tab of EasyCGTree_GUI);
they are written into this folder.


Available sets
--------------
Set (file)          Size     Genes / description

**Prokaryotes**
rp1.hmm             1.0 MB   16 ubiquitous ribosomal protein genes (18 domains) in Prokaryotes.
                             (http://dx.doi.org/10.1038/nature14486)
rp2.hmm             1.4 MB   23 ubiquitous ribosomal protein genes (27 domains) in Prokaryotes.
                             (http://dx.doi.org/10.1038/nature12352)

**Bacteria**
bac120.hmm         18.5 MB   120 ubiquitous genes in the domain Bacteria (default of EasyCGTree).
                             (http://dx.doi.org/10.1038/s41564-017-0012-7)
essential.hmm      15.2 MB   107 essential single-copy core genes in Bacteria.
                             (https://cdnsciencepub.com/doi/10.1139/gen-2015-0175)
ubcg1.hmm          11.6 MB   92 core genes named as up-to-date bacterial core genes.
                             (http://dx.doi.org/10.1007/s12275-018-8014-6)
ubcg2.hmm           9.4 MB   81 core genes named as up-to-date bacterial core genes.
                             (https://doi.org/10.1007/s12275-021-1231-4)
bac52.hmm           6.5 MB   52 ubiquitous ribosomal protein genes in Bacteria.
                             (https://www.kegg.jp/kegg/annotation/br01610.html)

**Archaea**
ar122.hmm          14.7 MB   122 ubiquitous genes in the domain Archaea.
                             (http://dx.doi.org/10.1038/s41564-017-0012-7)
uacg.hmm           13.2 MB   128 ubiquitous genes in the domain Archaea.
                             (https://link.springer.com/article/10.1007/s12275-023-00064-2)
ar53.hmm            6.6 MB   53 ubiquitous ribosomal protein genes in Archaea.
                             (https://www.kegg.jp/kegg/annotation/br01610.html)

**Others**
AXT.hmm             1.3 MB   The five key genes for astaxanthin synthesis in Prokaryotes.
                             (https://doi.org/10.1016/j.syapm.2025.126624)
CKC20.hmm           4.0 MB   20 core genes distinguishable among the Corynebacterium kroppenstedtii complex (CKC).
                             (https://doi.org/10.1093/jambio/lxad314)
ery288.hmm         44.5 MB   288 core genes of the family Erythrobacteraceae.
                             (http://dx.doi.org/10.1099/ijsem.0.004293)
rhodo268.hmm       40.4 MB   268 core genes of the family Rhodobacteraceae.
                             (https://doi.org/10.1099/ijsem.0.006156)
spi269.hmm         33.9 MB   269 core genes of the family Spirosomaceae.
                             (https://doi.org/10.1007/s12275-022-2102-3)


中文说明
--------
请把需要用到的 profile HMM 集（.hmm 文件）放在本文件夹中。

由于 HMM 文件较大（合计约 220 MB），且大多数用户只会用到其中少数几个，程序包中不包含 HMM 集。
请从 EasyCGTree5 仓库的 'HMM' 文件夹中单独下载所需的 HMM 集，保存到本文件夹，不要修改文件名（如 'bac120.hmm'）：

    GitHub：https://github.com/zdf1987/EasyCGTree5/tree/main/HMM
    Gitee： https://gitee.com/zdf1987/EasyCGTree5/tree/main/HMM    （中国大陆用户推荐）

下载单个文件的方法：在仓库的 'HMM' 文件夹中打开该文件，点击 'Download raw file'（GitHub）或“下载”（Gitee）。

EasyCGTree 的默认 HMM 集为 'bac120'：使用默认设置运行前，请先下载 'bac120.hmm'。
如果缺少所需的 HMM 集，EasyCGTree 和 EasyCGTree_GUI 会停止运行并提示查看本文件。
也可以通过给出 .hmm 文件的路径使用其他位置的 HMM 集（'-hmm /path/to/set.hmm'）。
用自己的基因家族构建的新 HMM 集（BuildHMM.pl 或 EasyCGTree_GUI 的“HMM构建”页面）会写入本文件夹。
各 HMM 集的说明和大小见上表。
