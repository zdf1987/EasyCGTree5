#!/usr/bin/perl
#
# EasyCGTree_SNP - SNP extraction and SNP tree from the core-gene nucleotide alignments of EasyCGTree
# Part of EasyCGTree 5.0, by Dao-Feng Zhang
#
# The script uses the results of EasyCGTree in '<input>_TEM'. When the nucleotide alignments needed
# are missing (or were made with other settings), the required EasyCGTree tasks are run first.
# Only core Perl modules are used. Put this script in the EasyCGTree directory (next to EasyCGTree.pl).
#
use warnings;
use strict;
use Getopt::Long;
use File::Copy qw(copy);
use File::Path qw(mkpath rmtree);
use File::Basename qw(basename dirname);
use File::Spec;
use Cwd qw(abs_path);
use FindBin qw($RealBin);

$| = 1;

my $VERSION = "5.0";    # version of EasyCGTree
my $UPDATE  = "2026-09-25";
my $START_TIME = time;
my ($sec,$min,$hour,$mday,$mon,$year) = localtime($START_TIME);
my $STAMP     = sprintf("%04d-%02d-%02d-%02d-%02d", $year+1900, $mon+1, $mday, $hour, $min);
my $START_STR = sprintf("%02d:%02d:%02d, %04d-%02d-%02d", $hour, $min, $sec, $year+1900, $mon+1, $mday);

my $IS_WIN   = ($^O eq 'MSWin32');
my $EXE      = $IS_WIN ? ".exe" : "";
my $HOME_DIR = $RealBin;
my $BIN_DIR  = (defined $ENV{"ECG_BIN"} && $ENV{"ECG_BIN"} ne "") ? $ENV{"ECG_BIN"} : "$HOME_DIR/bin";   # ECG_BIN: set by EasyCGTree_GUI
my $OPT_FILE = (-f "$BIN_DIR/tree_app-options.txt") ? "$BIN_DIR/tree_app-options.txt" : "$HOME_DIR/bin/tree_app-options.txt";
my $ECG      = "$HOME_DIR/EasyCGTree.pl";

my @usage = qq(
====== EasyCGTree_SNP ======
     EasyCGTree $VERSION, by Dao-Feng Zhang
     Update $UPDATE

SNPs are extracted from the concatenated core-gene nucleotide alignments of EasyCGTree, and a tree is
inferred from them. If the nucleotide alignments are missing in '<input>_TEM', or were made with other
settings than the ones given here, the required EasyCGTree tasks (predict, hmmsearch, refine, alignment)
are run automatically with '-seq nucl'.

Usage: perl EasyCGTree_SNP.pl -input <dir> [Options]

Essential Options:
-input <String>
	The input directory used (or to be used) with EasyCGTree. All genomes must be DNA sequences
	(or NCBI CDS files). Results are written to the directory that contains the input directory.

Optional Options:
-outdir <String>
	Output directory for all results (working directory '<input>_TEM', trees, tables, log files).
	[default: the directory that contains the input directory]
-ignore <String>
	Strains ignored when SNP sites are detected: comma-separated genome names, or a file with one name
	per line. These strains are still included in the SNP alignment and the tree, but sites that are
	variable only because of them are not taken as SNPs.
-max_missing <Decimal, 0..1>
	Largest fraction of the (not ignored) strains allowed to have a gap or an ambiguous base (N, R, Y ...)
	at a SNP site. [default: 0, i.e. core SNPs present in all strains]
-aln <String, 'trimmed' or 'original'>
	Alignments used: trimmed by trimAl (as for the EasyCGTree trees), or the original alignments
	before trimming. [default: trimmed]
-ref <String>
	Reference genome: the positions of the SNPs in its CDS (and the codon positions) are reported too.
-tree_app <String, 'iqtree' or 'fasttree'>
	Application used for tree inference. [default: iqtree]
	The IQ-TREE and FastTree command lines are taken from 'bin/tree_app-options.txt'.
-asc <String, 'fconst' or 'lewis'>
	Correction for using only variable sites (IQ-TREE only). [default: fconst]
	fconst: the numbers of invariant A, C, G and T sites of the analysed alignment are given to IQ-TREE
	        ('-fconst'), so branch lengths are substitutions per site of the core-gene alignment.
	lewis:  ascertainment bias correction of the model ('+ASC'); IQ-TREE also uses it by itself when a model
	        is selected for an alignment without invariant sites.
-no_tree <no value required>
	Only extract the SNPs, do not infer a tree.
-thread <Int>
	Number of threads. [default: 4]
-help <no value required>
	Display this message.
-options_json <no value required>
	Print the options of this script in JSON format (used by EasyCGTree_GUI) and exit.

Options passed to EasyCGTree when its tasks have to be run (also compared with the existing results):
-hmm <String>            Profile HMM. [default: the one recorded in '<input>_TEM', otherwise bac120]
-evalue <Real>           E-value for hmmsearch. [default: EasyCGTree default]
-gene_cutoff <Decimal>   Cutoff for omitting low-prevalence genes. [default: EasyCGTree default]
-genome_cutoff <Decimal> Cutoff for omitting low-quality genomes. [default: EasyCGTree default]
-trim <String>           trimAl standard ('gappyout', 'strict', 'strictplus', 'nogaps'). [default: EasyCGTree default]

Output (in the directory containing the input directory):
<input>.<hmm>.<genes>.snp.<tree_app>.tree   SNP tree
<input>.<hmm>.<genes>.snp.fas               SNP alignment (all strains)
<input>.<hmm>.<genes>.snp_positions.txt     position of each SNP (gene, alignment positions, alleles)
<input>.<hmm>.<genes>.snp_genes.txt         SNPs covered by each gene, in the order of the SNP alignment
<input>.<hmm>.<genes>.snp_distance.txt      pairwise SNP distances

);

########## Options ##########
my %opt;
GetOptions(\%opt, "input:s", "outdir:s", "ignore:s", "max_missing:f", "aln:s", "ref:s", "tree_app:s", "asc:s", "no_tree!",
	"thread:i", "hmm:s", "evalue:f", "gene_cutoff:f", "genome_cutoff:f", "trim:s", "help!", "options_json!")
	or die "@usage\n\nERROR: Unrecognized option (see above).\n\n";
if ($opt{"options_json"}) {
	&printOptionsJson();
	exit 0;
}
if (!%opt || $opt{"help"}) {
	print @usage;
	exit;
}
&usageError("a directory must be specified for option '-input'.") unless defined $opt{"input"} && $opt{"input"} =~ /\S/;

my $thread     = 4;
my $maxMissing = 0;
my $alnMode    = "trimmed";
my $tree_app   = "iqtree";
my $asc        = "fconst";
if (exists $opt{"thread"}) {
	&usageError("'-thread' should be an integer >= 1.") unless defined $opt{"thread"} && $opt{"thread"} >= 1;
	$thread = $opt{"thread"};
}
my $ncpuWarn = "";
{
	my $ncpu = &cpuCount();
	if ($ncpu && $thread > $ncpu) {
		$ncpuWarn = "WARNING: '-thread $thread' is larger than the number of CPU cores of this computer ($ncpu); $ncpu threads are used.\n\n" if exists $opt{"thread"};
		$thread = $ncpu;
	}
}
if (exists $opt{"max_missing"}) {
	&usageError("'-max_missing' should be a decimal between 0 and 1.") unless defined $opt{"max_missing"} && $opt{"max_missing"} >= 0 && $opt{"max_missing"} <= 1;
	$maxMissing = $opt{"max_missing"};
}
if (exists $opt{"aln"}) {
	$alnMode = $opt{"aln"};
	&usageError("'-aln' should be 'trimmed' or 'original'.") unless $alnMode eq "trimmed" || $alnMode eq "original";
}
if (exists $opt{"tree_app"}) {
	$tree_app = $opt{"tree_app"};
	&usageError("'-tree_app' should be 'iqtree' or 'fasttree'.") unless $tree_app eq "iqtree" || $tree_app eq "fasttree";
}
if (exists $opt{"asc"}) {
	$asc = $opt{"asc"};
	&usageError("'-asc' should be 'fconst' or 'lewis'.") unless $asc eq "fconst" || $asc eq "lewis";
}
if (exists $opt{"trim"}) {
	$opt{"trim"} = "nogaps" if $opt{"trim"} eq "nogap";
	&usageError("'-trim' should be 'gappyout', 'strict', 'strictplus' or 'nogaps'.") unless $opt{"trim"} =~ /^(gappyout|strict|strictplus|nogaps)$/;
}
foreach my $k ("gene_cutoff", "genome_cutoff") {
	next unless exists $opt{$k};
	&usageError("'-$k' should be a decimal between 0 and 1.") unless defined $opt{$k} && $opt{$k} >= 0 && $opt{$k} <= 1;
}

########## Paths (as in EasyCGTree) ##########
(my $inputArg = $opt{"input"}) =~ s/[\/\\]+$//;
my $inputDir = File::Spec->rel2abs($inputArg);
$inputDir = abs_path($inputDir) if -d $inputDir;
my $inputName = basename($inputDir);
my $outDir    = dirname($inputDir);
if (defined $opt{"outdir"} && $opt{"outdir"} ne "") {   # '-outdir': all results go there (same layout)
	(my $od = $opt{"outdir"}) =~ s/[\/\\]+$// unless $opt{"outdir"} =~ /^[\/\\]+$/;
	$od = $opt{"outdir"} unless defined $od;
	$outDir = File::Spec->rel2abs($od);
	&usageError("'-outdir $opt{'outdir'}' is a file, not a directory.") if -e $outDir && !-d $outDir;
	$outDir = abs_path($outDir) if -d $outDir;
}
unless (-d $outDir) {
	mkpath($outDir);
	&usageError("The output directory '$outDir' can't be created.") unless -d $outDir;
	$outDir = abs_path($outDir);
}
my $TEMdir    = "$outDir/${inputName}_TEM";
my $RECORD    = "$TEMdir/EasyCGTree_record.txt";
my $TYPEFILE  = "$TEMdir/Genome_SeqType.txt";
my $SNPdir    = "$TEMdir/TEM8_SNP";
&usageError("'$ECG' not found. Please put EasyCGTree_SNP.pl in the EasyCGTree directory.") unless -f $ECG;
&usageError("Neither the input directory '$inputDir' nor the working directory '$TEMdir' exists.") unless -d $inputDir || -d $TEMdir;

my $hmmWanted;
if (defined $opt{"hmm"} && $opt{"hmm"} ne "") {
	($hmmWanted = basename($opt{"hmm"})) =~ s/\.hmm$//i;
}

my %rec0 = &readRecord();
my $logName = "$outDir/$inputName.".($hmmWanted || $rec0{"hmm"} || "bac120").".snp.$tree_app\_$STAMP.log";
open(LOG, ">", $logName) or die "Can't open '$logName': $!\n";
$SIG{__DIE__} = sub { print LOG "\nERROR: $_[0]" if defined fileno(LOG); };

&msg("\n====== EasyCGTree_SNP ======\n	EasyCGTree $VERSION\nby Dao-Feng Zhang\n\n");
&msg("Job Started at: $START_STR\n\nInput directory: $inputDir\nOutput directory: $outDir\nWorking directory: $TEMdir\n\n$ncpuWarn");


############### Step 1: make sure the nucleotide alignments exist ###############
&msg("#============= Step 1: Check the EasyCGTree results ==============#\n\n");
my @tasksToRun = &planEasyCGTree();
if (@tasksToRun) {
	&msg("The following EasyCGTree tasks will be run (with '-seq nucl'): ".join(", ", @tasksToRun).".\n\n");
	# Settings not given by the user are taken from the previous EasyCGTree run (the record is reset
	# when the prediction is run again).
	my %prev = &readRecord();
	my %use = %opt;
	$use{"hmm"} = $prev{"hmm_file"} if !(defined $use{"hmm"} && $use{"hmm"} ne "") && $prev{"hmm_file"} && -f $prev{"hmm_file"};
	foreach my $k ("evalue", "gene_cutoff", "genome_cutoff", "trim") {
		$use{$k} = $prev{$k} if !exists $use{$k} && defined $prev{$k};
	}
	foreach my $t (@tasksToRun) {
		my @cmd = ($^X, $ECG, "-input", $inputDir, "-outdir", $outDir, "-task", $t);
		push @cmd, ("-thread", $thread) if $t eq "hmmsearch" || $t eq "alignment";
		push @cmd, ("-seq", "nucl") if $t eq "predict" || $t eq "alignment";
		if ($t eq "hmmsearch") {
			push @cmd, ("-hmm", $use{"hmm"}) if defined $use{"hmm"} && $use{"hmm"} ne "";
			push @cmd, ("-evalue", $use{"evalue"}) if exists $use{"evalue"};
		}
		if ($t eq "refine") {
			foreach my $k ("gene_cutoff", "genome_cutoff") {
				push @cmd, ("-$k", $use{$k}) if exists $use{$k};
			}
			push @cmd, ("-evalue", $opt{"evalue"}) if exists $opt{"evalue"};
		}
		push @cmd, ("-trim", $use{"trim"}) if $t eq "alignment" && exists $use{"trim"};
		&msg("Running: perl EasyCGTree.pl ".join(" ", @cmd[2..$#cmd])."\n\n");
		my $rc = system(@cmd);
		$rc = ($rc == -1) ? 255 : ($rc & 127) ? 128 + ($rc & 127) : ($rc >> 8);
		&fail("EasyCGTree '-task $t' stopped with an error (exit code $rc). Please see the messages above and the EasyCGTree log file in '$outDir'.") if $rc;
	}
	&msg("\nAll required EasyCGTree tasks have been completed.\n\n");
} else {
	&msg("The nucleotide alignments of EasyCGTree in '$TEMdir' are used.\n\n");
}
my %rec  = &readRecord();
my $hmm  = $rec{"hmm"} || $hmmWanted || "bac120";
my $trim = $rec{"trim"} || "strict";
&msg("Profile HMM: $hmm; trimAl: -$trim; gene_cutoff: ".($rec{'gene_cutoff'} // "NA")."; genome_cutoff: ".($rec{'genome_cutoff'} // "NA")."\n\n");


############### Step 2: read the gene alignments ###############
&msg("#============= Step 2: Read the gene alignments ==============#\n\n");
&removePath($SNPdir) if -e $SNPdir;
&makeDir($SNPdir);

opendir(my $dh, "$TEMdir/TEM6_AlnTrimmed") or &fail("Can't open '$TEMdir/TEM6_AlnTrimmed': $!");
my @genes = sort map { (my $g = $_) =~ s/\.fasta$//; $g } grep { /\.fasta$/ } readdir($dh);
closedir($dh);
&fail("No alignment was found in '$TEMdir/TEM6_AlnTrimmed'.") unless @genes;
my $nGenes = scalar @genes;
my $outPrefix = "$outDir/$inputName.$hmm.$nGenes";

my $trimal = &tool("trimal");
my (%taxa, @geneData);
foreach my $gi (0..$#genes) {
	my $g = $genes[$gi];
	my ($tIds, $tSeq) = &readFasta("$TEMdir/TEM6_AlnTrimmed/$g.fasta");
	my ($oIds, $oSeq) = &readFasta("$TEMdir/TEM5_Alignment/$g.fas.fasta");
	&fail("The original alignment '$TEMdir/TEM5_Alignment/$g.fas.fasta' is missing or empty.") unless @$oIds;
	my $tLen = @$tIds ? length($tSeq->{$tIds->[0]}) : 0;
	my $oLen = length($oSeq->{$oIds->[0]});

	# Columns of the original alignment kept by trimAl (the same trimming is repeated with '-colnumbering').
	my @map;
	my $mapLog = "$SNPdir/$g.colnumbering.txt";
	my $rc = &runCapture($mapLog, $trimal, "-in", "$TEMdir/TEM5_Alignment/$g.fas.fasta", "-out", "$SNPdir/tmp.fas", "-$trim", "-colnumbering");
	if ($rc == 0 && open(my $mh, "<", $mapLog)) {
		while (my $l = <$mh>) {
			if ($l =~ s/^#ColumnsMap\s*//) {
				@map = map { $_ + 1 } ($l =~ /(\d+)/g);    # 1-based
			}
		}
		close $mh;
	}
	unless (@map == $tLen) {
		&msg("WARNING: the columns kept by trimAl could not be determined for gene '$g'; its positions in the original alignment are reported as 'NA'.\n");
		@map = ();
	}
	unlink "$SNPdir/tmp.fas";

	my (%seq, @colTrim, @colOrig);
	if ($alnMode eq "trimmed") {
		%seq = %$tSeq;
		@colTrim = (1..$tLen);
		@colOrig = @map ? @map : ("NA") x $tLen;
	} else {
		%seq = %$oSeq;
		my %t2 = @map ? (map { $map[$_] => $_ + 1 } 0..$#map) : ();
		@colOrig = (1..$oLen);
		@colTrim = map { $t2{$_} || "-" } @colOrig;
	}
	my $len = ($alnMode eq "trimmed") ? $tLen : $oLen;
	$taxa{$_} = 1 foreach keys %seq;
	push @geneData, { "name" => $g, "seq" => \%seq, "len" => $len, "tLen" => $tLen, "oLen" => $oLen,
		"colTrim" => \@colTrim, "colOrig" => \@colOrig, "orig" => $oSeq };
}
my @taxa = sort keys %taxa;
my $nTaxa = scalar @taxa;
&msg("$nGenes gene alignments ('".($alnMode eq "trimmed" ? "trimmed by trimAl" : "original, not trimmed")."') of $nTaxa strains were read.\n");
&fail("At least 4 strains are needed to infer a tree ($nTaxa found).") if $nTaxa < 4 && !$opt{"no_tree"};

# The alignments must be nucleotide sequences.
{
	my $s = join("", map { my $gd = $_; join("", map { $gd->{"seq"}{$_} } sort keys %{ $gd->{"seq"} }) } @geneData[0..($nGenes > 5 ? 4 : $nGenes-1)]);
	$s =~ s/[-.?]//g;
	my $acgt = ($s =~ tr/ACGTUN//);
	if (!length($s) || $acgt / length($s) < 0.9) {
		&fail("The alignments in '$TEMdir' are not nucleotide alignments. Please run EasyCGTree '-task alignment -seq nucl' (this script does so when the record file is present).");
	}
}

########## Ignored strains and reference ##########
my %ignore;
if (defined $opt{"ignore"} && $opt{"ignore"} ne "") {
	my @names;
	if (-f $opt{"ignore"}) {
		open(my $ih, "<", $opt{"ignore"}) or &fail("Can't open '$opt{'ignore'}': $!");
		while (my $l = <$ih>) {
			$l =~ s/\r?\n$//;
			push @names, grep { /\S/ } split /[,\s]+/, $l;
		}
		close $ih;
	} else {
		@names = grep { /\S/ } split /\s*,\s*/, $opt{"ignore"};
	}
	my (@unknown, @screened);
	my @all = &allGenomes();
	foreach my $n (@names) {
		my $t = &matchTaxon($n, \@taxa);
		if (defined $t) { $ignore{$t} = 1; }
		elsif (defined &matchTaxon($n, \@all)) { push @screened, &matchTaxon($n, \@all); }
		else { push @unknown, $n; }
	}
	&msg("NOTE: the following strains given with '-ignore' are not in the alignments, because they were removed when the genomes were screened ('-genome_cutoff'); they need not be ignored: ".join(", ", @screened).".\n") if @screened;
	&fail("The following strains given with '-ignore' were not found in the alignments:\n	".join("\n	", @unknown)."\nAvailable strains:\n	".join("\n	", @taxa)) if @unknown;
}
my @cons = grep { !$ignore{$_} } @taxa;    # strains considered for SNP detection
my $nCons = scalar @cons;
&fail("At least 2 strains must remain when the ignored strains are left out ($nCons remain).") if $nCons < 2;
if (%ignore) {
	&msg("Strains ignored when detecting SNP sites (still included in the SNP alignment and tree): ".join(", ", sort keys %ignore).".\n");
}
my $ref;
if (defined $opt{"ref"} && $opt{"ref"} ne "") {
	$ref = &matchTaxon($opt{"ref"}, \@taxa);
	if (!defined $ref && defined(my $s = &matchTaxon($opt{"ref"}, [&allGenomes()]))) {
		&fail("The reference genome '$s' is not in the alignments, because it was removed when the genomes were screened (too few genes for '-genome_cutoff').\nPlease choose another reference genome, or run with a lower '-genome_cutoff'.\nAvailable strains:\n	".join("\n	", @taxa));
	}
	&fail("The reference genome '$opt{'ref'}' was not found in the alignments.\nAvailable strains:\n	".join("\n	", @taxa)) unless defined $ref;
	&msg("Reference genome for CDS positions: $ref\n");
}
my $maxMiss = int($maxMissing * $nCons + 1e-9);
&msg("A site is taken as a SNP when the $nCons considered strains show at least 2 different bases (A, C, G, T),\nand at most $maxMiss of them have a gap or an ambiguous base (-max_missing $maxMissing).\n\n");


############### Step 3: SNP detection ###############
&msg("#============= Step 3: SNP detection ==============#\n\n");
my %refCds = defined $ref ? &refCdsIds($ref) : ();
my (@snpCols, %snpSeq, %constCount);
my ($concatPos, $nInformative, $nMulti) = (0, 0, 0);
my @geneSummary;
$snpSeq{$_} = "" foreach @taxa;
foreach my $gd (@geneData) {
	my $len = $gd->{"len"};
	my %s;
	foreach my $t (@taxa) {
		$s{$t} = exists $gd->{"seq"}{$t} ? $gd->{"seq"}{$t} : "-" x $len;
		$s{$t} .= "-" x ($len - length($s{$t})) if length($s{$t}) < $len;
	}
	# Ungapped positions of the reference in the original alignment of this gene.
	my @refPos;
	if (defined $ref && exists $gd->{"orig"}{$ref}) {
		my $n = 0;
		foreach my $c (split //, $gd->{"orig"}{$ref}) {
			$n++ if $c =~ /[A-Z]/;
			push @refPos, ($c =~ /[A-Z]/ ? $n : "-");
		}
	}
	my $first = @snpCols + 1;
	my @origList;
	for (my $i = 0; $i < $len; $i++) {
		$concatPos++;
		my (%cnt, $miss);
		$miss = 0;
		foreach my $t (@cons) {
			my $b = substr($s{$t}, $i, 1);
			if ($b =~ /[ACGT]/) { $cnt{$b}++; } else { $miss++; }
		}
		my $nAll = keys %cnt;
		if ($nAll < 2 || $miss > $maxMiss) {
			# Invariant column: counted for '-fconst' only if all strains (including ignored ones) have the same base.
			if ($nAll == 1 && $miss == 0) {
				my ($b) = keys %cnt;
				my $same = 1;
				foreach my $t (keys %ignore) {
					if (substr($s{$t}, $i, 1) ne $b) { $same = 0; last; }
				}
				$constCount{$b}++ if $same;
			}
			next;
		}
		foreach my $t (@taxa) {
			my $b = substr($s{$t}, $i, 1);
			$b = "-" unless $b =~ /[ACGT]/;
			$snpSeq{$t} .= $b;
		}
		my @al = sort { $cnt{$b} <=> $cnt{$a} || $a cmp $b } keys %cnt;
		my $inf = (grep { $cnt{$_} >= 2 } @al) >= 2 ? "yes" : "no";
		$nInformative++ if $inf eq "yes";
		$nMulti++ if $nAll > 2;
		my $oPos = $gd->{"colOrig"}[$i];
		my %row = ("gene" => $gd->{"name"}, "concat" => $concatPos, "trim" => $gd->{"colTrim"}[$i], "orig" => $oPos,
			"alleles" => join("/", @al), "counts" => join(",", map { "$_:$cnt{$_}" } @al), "missing" => $miss, "inf" => $inf);
		if (defined $ref) {
			my $rp = ($oPos ne "NA" && @refPos) ? $refPos[$oPos - 1] : "-";
			$row{"refId"}  = $refCds{$gd->{"name"}} || "-";
			$row{"refPos"} = $rp;
			$row{"codon"}  = ($rp ne "-") ? (($rp - 1) % 3) + 1 : "-";
			my $rb = substr($s{$ref}, $i, 1);
			$row{"refBase"} = ($rb =~ /[ACGT]/) ? $rb : "-";
		}
		push @snpCols, \%row;
		push @origList, $oPos;
	}
	my $nSnp = @snpCols - $first + 1;
	push @geneSummary, [$gd->{"name"}, $gd->{"tLen"}, $gd->{"oLen"}, $concatPos - $len + 1, $concatPos,
		scalar(grep { exists $gd->{"seq"}{$_} } @taxa), $nSnp, ($nSnp ? "$first-".scalar(@snpCols) : "-"), ($nSnp ? join(",", @origList) : "-")];
}
my $nSnp = scalar @snpCols;
&msg("Alignment analysed: $concatPos sites of $nGenes genes.\n");
&msg("SNP sites: $nSnp (parsimony-informative: $nInformative; with more than 2 bases: $nMulti).\n");
&msg("Invariant sites (same base in all strains): ".join(", ", map { "$_ ".($constCount{$_} || 0) } qw(A C G T))."\n\n");
&fail("No SNP was found. A SNP tree cannot be inferred.") unless $nSnp;

########## Output files ##########
# Results of a previous run with the same name are kept as backups. The SNP trees of the previous run
# belong to the old SNP alignment, so they are backed up, too (also with '-no_tree').
&backupFile($_) foreach map { "$outPrefix.$_" } qw(snp.fas snp_positions.txt snp_genes.txt snp_distance.txt snp.iqtree.tree snp.fasttree.tree);
my $snpFas = "$outPrefix.snp.fas";
open(my $sf, ">", $snpFas) or &fail("Can't open '$snpFas': $!");
foreach my $t (@taxa) {
	print $sf ">$t\n";
	for (my $p = 0; $p < length($snpSeq{$t}); $p += 60) {
		print $sf substr($snpSeq{$t}, $p, 60), "\n";
	}
}
close $sf;

my $posFile = "$outPrefix.snp_positions.txt";
open(my $pf, ">", $posFile) or &fail("Can't open '$posFile': $!");
print $pf "# SNP positions. SNP_No: column in '".basename($snpFas)."'; Concat_pos: position in the concatenated ".($alnMode eq "trimmed" ? "trimmed " : "original ")."alignments of all genes;\n";
print $pf "# Trimmed_aln_pos / Original_aln_pos: column in the gene alignment after / before trimming by trimAl ('$TEMdir/TEM6_AlnTrimmed' / 'TEM5_Alignment');\n";
print $pf "# Alleles, Allele_counts and Missing refer to the ".$nCons." strains considered".(%ignore ? " (ignored strains are left out)" : "").".\n";
print $pf "# Ref_CDS_pos: position in the CDS of the reference genome; Codon_pos: position in the codon (1-3).\n" if defined $ref;
print $pf join("\t", qw(SNP_No Gene Gene_No Concat_pos Trimmed_aln_pos Original_aln_pos Alleles Allele_counts Missing Parsimony_informative),
	(defined $ref ? qw(Ref_CDS_ID Ref_CDS_pos Codon_pos Ref_base) : ())), "\n";
my %geneNo = map { $genes[$_] => $_ + 1 } 0..$#genes;
foreach my $i (0..$#snpCols) {
	my $r = $snpCols[$i];
	print $pf join("\t", $i + 1, $r->{"gene"}, $geneNo{$r->{"gene"}}, $r->{"concat"}, $r->{"trim"}, $r->{"orig"}, $r->{"alleles"}, $r->{"counts"}, $r->{"missing"}, $r->{"inf"},
		(defined $ref ? ($r->{"refId"}, $r->{"refPos"}, $r->{"codon"}, $r->{"refBase"}) : ())), "\n";
}
close $pf;

my $geneFile = "$outPrefix.snp_genes.txt";
open(my $gf, ">", $geneFile) or &fail("Can't open '$geneFile': $!");
print $gf "# SNPs covered by each gene, in the order of the concatenation and of the SNP alignment.\n";
print $gf "# SNP_range: columns in '".basename($snpFas)."'; Original_aln_positions: the SNP columns in the original alignment of the gene (before trimming), in the same order.\n";
print $gf join("\t", qw(Gene_No Gene Trimmed_length Original_length Concat_start Concat_end Strains SNPs SNP_range Original_aln_positions)), "\n";
foreach my $i (0..$#geneSummary) {
	print $gf join("\t", $i + 1, @{ $geneSummary[$i] }), "\n";
}
close $gf;

my $distFile = "$outPrefix.snp_distance.txt";
open(my $df, ">", $distFile) or &fail("Can't open '$distFile': $!");
print $df "# Pairwise SNP distances: number of SNP sites at which two strains have different bases (sites with a gap or an ambiguous base in either strain are not counted).\n";
print $df join("\t", "Strain", @taxa), "\n";
foreach my $a (@taxa) {
	my @row;
	foreach my $b (@taxa) {
		my $d = 0;
		if ($a ne $b) {
			my ($x, $y) = ($snpSeq{$a}, $snpSeq{$b});
			for (my $p = 0; $p < $nSnp; $p++) {
				my ($u, $v) = (substr($x, $p, 1), substr($y, $p, 1));
				$d++ if $u ne "-" && $v ne "-" && $u ne $v;
			}
		}
		push @row, $d;
	}
	print $df join("\t", $a, @row), "\n";
}
close $df;
&msg("SNP alignment:        $snpFas\nSNP positions:        $posFile\nSNPs of each gene:    $geneFile\nPairwise distances:   $distFile\n\n");


# Settings of this SNP extraction (used by EasyCGTree_GUI to decide whether it must be repeated).
{
	my $given = (defined $opt{"ignore"} && $opt{"ignore"} ne "") ? join(",", sort { lc($a) cmp lc($b) } grep { /\S/ } split /\s*,\s*/, $opt{"ignore"}) : "";
	&writeRecord("snp_done" => &now(), "snp_hmm" => $hmm, "snp_genes" => $nGenes, "snp_ignore" => $given,
		"snp_max_missing" => $maxMissing, "snp_aln" => $alnMode, "snp_ref" => (defined $ref ? $ref : ""));
}

############### Step 4: SNP tree ###############
unless ($opt{"no_tree"}) {
	&msg("#============= Step 4: SNP tree ==============#\n\n");
	my ($iq, $fast) = &treeOptions();
	my $outTree = "$outPrefix.snp.$tree_app.tree";
	my $pre = "$SNPdir/$inputName.$hmm.$nGenes.snp";
	unlink $outTree;
	if ($tree_app eq "iqtree") {
		&fail("No command line for IQ-TREE was found in '$OPT_FILE'.") unless $iq;
		my ($prog, @args) = split /\s+/, $iq;
		$prog = &tool($prog);
		&fail("The IQ-TREE command line in '$OPT_FILE' contains '-m MF' (model selection only, no tree).") if $iq =~ /(^|\s)-m\s+MF(\s|$)/;
		if ($asc eq "fconst") {
			my $fc = join(",", map { $constCount{$_} || 0 } qw(A C G T));
			push @args, ("-fconst", $fc);
			&msg("Invariant sites added to IQ-TREE: -fconst $fc (A,C,G,T)\n");
		} elsif ($asc eq "lewis") {
			my $mi = -1;
			foreach my $j (0..$#args-1) { $mi = $j + 1 if $args[$j] eq "-m"; }
			if ($mi >= 0) {
				$args[$mi] .= "+ASC" unless $args[$mi] =~ /ASC/;
			} else {
				push @args, ("-m", "MFP+ASC");
			}
			&msg("Ascertainment bias correction: ".($mi >= 0 ? "-m $args[$mi]" : "-m MFP+ASC")."\n");
		}
		push @args, "-redo" unless grep { $_ eq "-redo" } @args;
		push @args, ("-st", "DNA", "-s", $snpFas, "-pre", $pre, "-T", $thread);
		my $rc = &runCmd($prog, @args);
		if ($rc) {
			my $m = "IQ-TREE exited with an error (exit code $rc).".($rc >= 128 ? " The process was killed, probably because it ran out of memory." : "");
			if (open(my $lh, "<", "$pre.log")) {
				my @e = grep { /ERROR/ } <$lh>;
				close $lh;
				$m .= "\nMessages from '$pre.log':\n".join("", @e) if @e;
				$m .= "\nThis is a known problem of IQ-TREE 2.2.0 with several threads. Please update IQ-TREE or use '-thread 1'." if grep { /Tree taxa and alignment sequence do not match/ } @e;
			}
			&fail($m);
		}
		my $iqTree = (-s "$pre.contree") ? "$pre.contree" : "$pre.treefile";
		&checkTree($iqTree, "IQ-TREE", "$pre.log");
		&copyTree($iqTree, $outTree);
		my ($model, $how, $lnl, $bic) = &iqtreeModel("$pre.iqtree");
		&msg("\nUseful Information from IQ-TREE:\nModel	Model_selection	LogL	BIC\n$model	$how	$lnl	$bic\n\n");
	} else {
		&fail("No command line for FastTree was found in '$OPT_FILE'.") unless $fast;
		my ($prog, @args) = split /\s+/, $fast;
		$prog = &tool($prog);
		push @args, ("-nt", "-gtr") unless grep { $_ eq "-nt" } @args;
		$ENV{"OMP_NUM_THREADS"} = $thread;
		my $rc = &runCapture("$pre.fasttree.log", $prog, @args, "-out", $outTree, $snpFas);
		&fail("FastTree exited with an error (exit code $rc).".($rc >= 128 ? " The process was killed, probably because it ran out of memory." : "")." See '$pre.fasttree.log'.") if $rc;
		&msg("NOTE: FastTree has no correction for using only variable sites; branch lengths are substitutions per SNP site.\n");
	}
	&checkTree($outTree, ($tree_app eq "iqtree" ? "IQ-TREE" : "FastTree"), ($tree_app eq "fasttree" ? "$pre.fasttree.log" : "$pre.log"));
	&msg("The SNP tree has been written successfully in '$outTree'.\n\n");
}

my $el = time - $START_TIME;
&msg(sprintf("Running time: %d d %d h %d min %d sec.\n\nEasyCGTree_SNP (EasyCGTree %s)\n\n", int($el/86400), int(($el%86400)/3600), int(($el%3600)/60), $el%60, $VERSION));
close LOG;
exit 0;


############################## Subroutines ##############################

sub msg {
	my $m = shift;
	print $m;
	print LOG $m if defined fileno(LOG);
}

sub fail {
	my $m = shift;
	print "\n\nERROR: $m\n\n";
	print LOG "\n\nERROR: $m\n\n" if defined fileno(LOG);
	close LOG if defined fileno(LOG);
	exit 1;
}

sub usageError {
	print @usage, "\nERROR: $_[0]\n\n";
	exit 1;
}

# Full path of a program: the environment variable ECG_<NAME> (e.g. ECG_IQTREE2; set by EasyCGTree_GUI)
# if it is set, otherwise the 'bin' directory ('.exe' is added on Windows).
sub tool {
	my $name = shift;
	$name =~ s/\.exe$//i;
	(my $key = "ECG_".uc($name)) =~ s/[^A-Z0-9_]/_/g;
	my $env = (defined $ENV{$key} && $ENV{$key} ne "") ? 1 : 0;
	my $p = $env ? $ENV{$key} : "$BIN_DIR/$name$EXE";
	unless (-e $p) {
		&fail($env ? "The program '$p' given by the environment variable $key does not exist." : "The program '$name$EXE' cannot be found in '$BIN_DIR'.");
	}
	if (!$IS_WIN && !-x $p) {
		&fail("The program '$p' is not executable. Please run: chmod +x \"$p\"");
	}
	return $p;
}

sub runCmd {
	my $rc = system(@_);
	&fail("Can't execute '$_[0]': $!") if $rc == -1;
	return ($rc & 127) ? 128 + ($rc & 127) : ($rc >> 8);
}

sub runCapture {
	my ($logfile, @cmd) = @_;
	open(my $oldOut, ">&", \*STDOUT) or die "Can't duplicate STDOUT: $!\n";
	open(my $oldErr, ">&", \*STDERR) or die "Can't duplicate STDERR: $!\n";
	open(STDOUT, ">", $logfile) or die "Can't open '$logfile': $!\n";
	open(STDERR, ">&", \*STDOUT) or die "Can't redirect STDERR: $!\n";
	my $rc = system(@cmd);
	my $err = $!;
	open(STDOUT, ">&", $oldOut) or die "Can't restore STDOUT: $!\n";
	open(STDERR, ">&", $oldErr) or die "Can't restore STDERR: $!\n";
	&fail("Can't execute '$cmd[0]': $err") if $rc == -1;
	return ($rc & 127) ? 128 + ($rc & 127) : ($rc >> 8);
}

sub makeDir {
	my $d = shift;
	mkpath($d) unless -d $d;
	&fail("Permission denied to create directory '$d'.") unless -d $d;
}

sub removePath {
	my $d = shift;
	rmtree($d);
	&fail("'$d' already exists, but permission denied when deleting it.") if -e $d;
}

# An existing result file with the same name is kept: it is renamed to '<file>.<date of the file>.bak'.
sub backupFile {
	my $f = shift;
	return unless -e $f;
	my @t = localtime((stat($f))[9]);
	my $s = sprintf("%04d-%02d-%02d-%02d-%02d", $t[5] + 1900, $t[4] + 1, $t[3], $t[2], $t[1]);
	my $b = "$f.$s.bak";
	my $i = 1;
	$b = "$f.$s." . $i++ . ".bak" while -e $b;
	if (rename($f, $b)) {
		&msg("NOTE: the existing file '".basename($f)."' was renamed to '".basename($b)."'.\n");
	} else {
		&msg("WARNING: the existing file '$f' could not be renamed ($!); it is overwritten.\n");
	}
}

sub now {
	my ($s, $mi, $h, $d, $mo, $y) = localtime();
	return sprintf("%04d-%02d-%02d %02d:%02d:%02d", $y + 1900, $mo + 1, $d, $h, $mi, $s);
}

sub writeRecord {
	my %new = @_;
	my %r = &readRecord();
	@r{keys %new} = values %new;
	open(my $fh, ">", $RECORD) or return;
	print $fh "# EasyCGTree record of the tasks performed in this directory. Please do not edit.\n";
	print $fh "$_=$r{$_}\n" foreach sort keys %r;
	close $fh;
}

sub readRecord {
	my %r;
	open(my $fh, "<", $RECORD) or return %r;
	while (my $l = <$fh>) {
		$l =~ s/\r?\n$//;
		next if $l =~ /^\s*#/ || $l !~ /=/;
		my ($k, $v) = split /=/, $l, 2;
		$r{$k} = $v;
	}
	close $fh;
	return %r;
}

# Decide which EasyCGTree tasks are needed to get nucleotide alignments made with the requested settings.
sub planEasyCGTree {
	my @order = qw(predict hmmsearch refine alignment);
	my %r = &readRecord();
	my $start;    # index of the first task to run
	my @why;

	# 1. CDS prediction with DNA sequences (.fnn) for every genome
	my %types;
	if (open(my $th, "<", $TYPEFILE)) {
		while (my $l = <$th>) {
			$l =~ s/\r?\n$//;
			my @f = split /\t/, $l;
			next if !@f || $f[0] eq "Genome";
			$types{$f[0]} = $f[1];
		}
		close $th;
	}
	if (!-d $TEMdir || !%types || !$r{"predict_done"}) {
		$start = 0;
		push @why, "no CDS prediction of EasyCGTree 5.0 was found";
	} else {
		my @prot = sort grep { $types{$_} eq "prot" } keys %types;
		if (@prot) {
			&fail("The following genomes were provided as protein sequences, so no DNA sequences are available for SNPs:\n	".join("\n	", @prot)."\nPlease replace them with genome (DNA) sequences or NCBI CDS files in '$inputDir' (or remove them).");
		}
		my @noFnn = sort grep { !-s "$TEMdir/TEM0_CDS/$_.fnn" } keys %types;
		if (@noFnn) {
			$start = 0;
			push @why, "the DNA sequences of the CDS were not kept for ".scalar(@noFnn)." genomes";
		}
	}
	&fail("The input directory '$inputDir' does not exist, so the CDS prediction cannot be run.") if defined $start && $start == 0 && !-d $inputDir;

	# 2. HMM search with the requested profile HMM and E-value
	if (!defined $start) {
		if (!$r{"hmmsearch_done"}) {
			$start = 1; push @why, "no HMM search was found";
		} elsif (defined $hmmWanted && $hmmWanted ne ($r{"hmm"} || "")) {
			$start = 1; push @why, "the HMM search was done with '".($r{'hmm'} || "?")."', not '$hmmWanted'";
		} elsif (exists $opt{"evalue"} && $opt{"evalue"} > ($r{"evalue"} || 0)) {
			$start = 1; push @why, "the HMM search was done with a smaller E-value ($r{'evalue'})";
		}
	}
	# 3. Screening with the requested cutoffs
	if (!defined $start) {
		if (!$r{"refine_done"}) {
			$start = 2; push @why, "no screening of the HMM search results was found";
		} else {
			foreach my $k ("gene_cutoff", "genome_cutoff") {
				if (exists $opt{$k} && (!defined $r{$k} || $opt{$k} != $r{$k})) {
					$start = 2; push @why, "'-$k $opt{$k}' differs from the value used before (".($r{$k} // "NA").")";
				}
			}
			if (exists $opt{"evalue"} && $opt{"evalue"} != ($r{"refine_evalue"} // $r{"evalue"})) {
				$start = 2; push @why, "'-evalue $opt{'evalue'}' differs from the value used before";
			}
		}
	}
	# 4. Nucleotide alignments with the requested trimming
	if (!defined $start) {
		if (!$r{"alignment_done"} || ($r{"alignment_seq"} || "") ne "nucl") {
			$start = 3; push @why, "no nucleotide alignments made after the last screening were found";
		} elsif (exists $opt{"trim"} && $opt{"trim"} ne ($r{"trim"} || "")) {
			$start = 3; push @why, "'-trim $opt{'trim'}' differs from the trimming used before (".($r{'trim'} || "NA").")";
		} elsif (!-d "$TEMdir/TEM5_Alignment" || !-d "$TEMdir/TEM6_AlnTrimmed") {
			$start = 3; push @why, "the alignment directories are missing";
		}
	}
	return () unless defined $start;
	&msg("Reason: ".join("; ", @why).".\n");
	return @order[$start..$#order];
}

# Find a strain name given by the user: exact, or after the renaming rules of EasyCGTree.
# All genomes of the input directory (from the CDS prediction), including those removed by the screening.
sub allGenomes {
	my @g;
	if (open(my $th, "<", $TYPEFILE)) {
		while (my $l = <$th>) {
			$l =~ s/\r?\n$//;
			my ($n) = split /\t/, $l;
			push @g, $n if defined $n && $n ne "" && $n ne "Genome";
		}
		close $th;
	}
	return @g;
}

sub matchTaxon {
	my ($n, $taxa) = @_;
	my %lc = map { lc($_) => $_ } @$taxa;
	return $lc{lc $n} if exists $lc{lc $n};
	my $m = $n;
	if ($m =~ /(GC[AF])[_-](\d{9})[._-](\d+)/ && exists $lc{lc "$1-$2-$3"}) {
		return $lc{lc "$1-$2-$3"};
	}
	$m =~ s/\.\w+$//;
	$m =~ s/[\s_.()]+/-/g;
	return $lc{lc $m} if exists $lc{lc $m};
	return undef;
}

# CDS ID (protein ID) of the reference genome for each gene, from the screened HMM search results.
sub refCdsIds {
	my $g = shift;
	my %id;
	opendir(my $d, "$TEMdir/TEM2_HMMsearch_outS") or return %id;
	my ($f) = grep { /__\Q$g\E\.fas$/ } readdir($d);
	closedir($d);
	return %id unless $f;
	open(my $fh, "<", "$TEMdir/TEM2_HMMsearch_outS/$f") or return %id;
	while (my $l = <$fh>) {
		next if $l =~ /^#/ || $l !~ /\S/;
		my @x = split /\s+/, $l;
		my $gene = ($x[3] =~ /-/) ? $x[2] : $x[3];
		$id{$gene} = $x[0];
	}
	close $fh;
	return %id;
}

sub readFasta {
	my $f = shift;
	my (@ids, %seq, $id);
	open(my $fh, "<", $f) or return (\@ids, \%seq);
	while (my $l = <$fh>) {
		$l =~ s/\r?\n$//;
		next unless length $l;
		if ($l =~ s/^\>//) {
			($id) = split /\s+/, $l;
			push @ids, $id unless exists $seq{$id};
			$seq{$id} = "";
		} elsif (defined $id) {
			$l =~ s/\s+//g;
			$seq{$id} .= uc($l);
		}
	}
	close $fh;
	return (\@ids, \%seq);
}

sub treeOptions {
	my ($iq, $fast);
	open(my $oh, "<", $OPT_FILE) or &fail("Can't open '$OPT_FILE': $!");
	while (my $l = <$oh>) {
		$l =~ s/\r?\n$//;
		if ($l =~ s/^\s*IQ-TREE Command-line=//) { $iq = $l; }
		elsif ($l =~ s/^\s*FastTree Command-line=//) { $fast = $l; }
	}
	close $oh;
	return ($iq, $fast);
}

sub checkTree {
	my ($file, $app, $appLog) = @_;
	my $why;
	if (!-e $file) { $why = "was not created"; }
	elsif (-z $file) { $why = "is empty"; }
	else {
		open(my $fh, "<", $file) or &fail("Can't open '$file': $!");
		local $/;
		my $c = <$fh>;
		close $fh;
		$why = "does not contain a complete tree" unless $c =~ /\(.*\)[^;]*;/s;
	}
	return unless $why;
	&fail("The SNP tree '$file' $why.\nThis occurs probably due to memory limitation (out of memory) or an interruption when running $app.".(($appLog && -e $appLog) ? "\nPlease check the messages of $app in '$appLog'." : ""));
}

sub copyTree {
	my ($from, $to) = @_;
	open(my $fi, "<", $from) or &fail("Can't open '$from': $!");
	local $/;
	my $t = <$fi>;
	close $fi;
	$t =~ s/\s+$//;
	open(my $fo, ">", $to) or &fail("Can't open '$to': $!");
	print $fo "$t\n";
	close $fo;
}

sub iqtreeModel {
	my $f = shift;
	my ($model, $how, $lnl, $bic) = ("NA", "NA", "NA", "NA");
	open(my $fh, "<", $f) or return ($model, $how, $lnl, $bic);
	while (my $l = <$fh>) {
		$l =~ s/\r?\n$//;
		if ($l =~ /^Best-fit model according to (\w+):\s*(\S+)/) { ($how, $model) = ("best-fit ($1)", $2); }
		elsif ($l =~ /^Model of substitution:\s*(\S+)/ && $model eq "NA") { ($how, $model) = ("user-specified", $1); }
		elsif ($l =~ /^Log-likelihood of the tree:\s*(\S+)/) { $lnl = $1; }
		elsif ($l =~ /^Bayesian information criterion \(BIC\) score:\s*(\S+)/) { $bic = $1; }
	}
	close $fh;
	return ($model, $how, $lnl, $bic);
}

# Options of this script in JSON format, for EasyCGTree_GUI (JSON::PP is a core module since Perl 5.14).
sub printOptionsJson {
	require JSON::PP;
	my @o = (
		{ "name" => "input", "type" => "dir", "default" => undef, "required" => 1 },
		{ "name" => "outdir", "type" => "dir", "default" => undef },
		{ "name" => "ignore", "type" => "text", "default" => undef },
		{ "name" => "max_missing", "type" => "float", "default" => 0, "min" => 0, "max" => 1 },
		{ "name" => "aln", "type" => "choice", "default" => "trimmed", "values" => ["trimmed", "original"] },
		{ "name" => "ref", "type" => "text", "default" => undef },
		{ "name" => "tree_app", "type" => "choice", "default" => "iqtree", "values" => ["iqtree", "fasttree"] },
		{ "name" => "asc", "type" => "choice", "default" => "fconst", "values" => ["fconst", "lewis"] },
		{ "name" => "no_tree", "type" => "bool", "default" => 0 },
		{ "name" => "thread", "type" => "int", "default" => 4, "min" => 1 },
		{ "name" => "hmm", "type" => "text", "default" => undef },
		{ "name" => "evalue", "type" => "float", "default" => undef },
		{ "name" => "gene_cutoff", "type" => "float", "default" => undef, "min" => 0, "max" => 1 },
		{ "name" => "genome_cutoff", "type" => "float", "default" => undef, "min" => 0, "max" => 1 },
		{ "name" => "trim", "type" => "choice", "default" => undef, "values" => ["gappyout", "strict", "strictplus", "nogaps"] },
	);
	print JSON::PP->new->canonical(1)->encode({ "script" => "EasyCGTree_SNP.pl", "version" => $VERSION, "options" => \@o }), "\n";
}

# Number of CPU cores (0 if unknown). Uses only core Perl: nproc / /proc/cpuinfo (Linux),
# sysctl (macOS), NUMBER_OF_PROCESSORS (Windows).
sub cpuCount {
	my $n = "";
	if ($^O eq 'MSWin32') {
		$n = $ENV{"NUMBER_OF_PROCESSORS"} || "";
	} elsif ($^O eq 'darwin') {
		$n = `sysctl -n hw.ncpu 2>/dev/null` || "";
	} else {
		$n = `nproc 2>/dev/null` || "";
		if ($n !~ /\d/ && open(my $fh, "<", "/proc/cpuinfo")) {
			$n = scalar grep { /^processor\s*:/ } <$fh>;
			close $fh;
		}
	}
	$n =~ s/\s+//g;
	return ($n =~ /^\d+$/ && $n > 0) ? $n : 0;
}
