#!/usr/bin/perl
#
# EasyCGTree - a pipeline for prokaryotic phylogenomic analysis based on core gene sets
# Version 5.0 by Dao-Feng Zhang
#
# Only core Perl modules are used, so no extra installation is needed.
# The same script runs on Linux, macOS and Windows: the programs in 'bin' are called
# with '.exe' on Windows, and all file operations are done in Perl.
#
use warnings;
use strict;
use Getopt::Long;
use File::Copy qw(copy move);
use File::Path qw(mkpath rmtree);
use File::Basename qw(basename dirname);
use File::Spec;
use Cwd qw(abs_path);
use FindBin qw($RealBin);
use Digest::MD5;

$| = 1;    # flush output immediately (progress lines use "\r")

my $VERSION = "5.0";
my $UPDATE  = "2026-10-07";
my $DRY_RUN = 0;

my $START_TIME = time;
my ($sec,$min,$hour,$mday,$mon,$year) = localtime($START_TIME);
my $STAMP     = sprintf("%04d-%02d-%02d-%02d-%02d", $year+1900, $mon+1, $mday, $hour, $min);
my $START_STR = sprintf("%02d:%02d:%02d, %04d-%02d-%02d", $hour, $min, $sec, $year+1900, $mon+1, $mday);

########## Installation paths (independent of the current working directory) ##########
my $IS_WIN   = ($^O eq 'MSWin32');
my $EXE      = $IS_WIN ? ".exe" : "";
my $HOME_DIR = $RealBin;                 # directory holding EasyCGTree.pl
my $BIN_DIR  = (defined $ENV{"ECG_BIN"} && $ENV{"ECG_BIN"} ne "") ? $ENV{"ECG_BIN"} : "$HOME_DIR/bin";   # ECG_BIN: set by EasyCGTree_GUI
my $HMM_DIR  = "$HOME_DIR/HMM";
my $OPT_FILE = (-f "$BIN_DIR/tree_app-options.txt") ? "$BIN_DIR/tree_app-options.txt" : "$HOME_DIR/bin/tree_app-options.txt";

my @usage=qq(
====== EasyCGTree ======
     Version $VERSION by Dao-Feng Zhang
     Update $UPDATE

Usage: perl EasyCGTree.pl [Options]
       (EasyCGTree.pl can be called from any directory, e.g. perl /path/to/EasyCGTree/EasyCGTree.pl -input myGenomes)

Essential Options:
-input <String>
	Input data (DNA and/or protein sequences in FASTA format) directory. It can be located anywhere.
	All results (log file, trees and the working directory '<input>_TEM') are written to the directory that contains the input directory,
	or to the directory given with '-outdir'.
	If you want a tree constructed using DNA/nucleotide sequences (-seq nucl), all the genomes must be provided as DNA/nucleotide sequences.
	One file per genome, except for NCBI annotation files (see below).

	NCBI annotation files (files whose names contain the assembly accession, e.g. 'GCF_000009925.1_ASM992v1_...'):
	- protein file with locus_tag ('translated_cds.faa') + CDS file ('cds_from_genomic.fna') of the same assembly:
	  the sequences are paired by locus_tag (>=90% of the proteins must be paired), and no CDS prediction is needed.
	- CDS file ('cds_from_genomic.fna') alone: the CDS are translated (genetic code 11).
	- A protein file without locus_tag ('protein.faa') cannot be used together with a CDS file of the same assembly.
	- If the genome ('genomic.fna') is also provided, the annotation files are used and the genome is ignored.
	These genomes are named by their accession (e.g. 'GCF-000009925-1'); protein and CDS DNA sequences can both be used.

Optional Options:
-keep_CDS_nucl 	<no value required>
	Also output the CDS as nucleotide/DNA sequences (besides the protein sequences) for the genomes provided as DNA sequences.
	With '-seq nucl' the DNA sequences of CDS are kept automatically, so this option is not needed in that case.
-evalue <Real, 10..0>
	Expect value for screening hmmsearch hits. [default: 1e-10]
-gene_cutoff <Decimal, 0..1>
	Cutoff for omitting low-prevalence genes. [default: 0.8]
-genome_cutoff <Decimal, 0..1>
	Cutoff for omitting low-quality genomes. [default: 0.8]
-outdir <String>
	Output directory for all results (working directory '<input>_TEM', trees, tables, log files).
	[default: the directory that contains the input directory]
-dry_run <no value required>
	Only check the input directory: show how the input files would be grouped and prepared (genome
	names, sequence types, NCBI annotation files), without renaming or running anything.
-help 	<no value required>
	Display this message.
-options_json <no value required>
	Print the options of this script in JSON format (used by EasyCGTree_GUI) and exit.
-hmm <String, 'bac120', 'rp1', et al., or the path of a .hmm file>
	Profile HMM used for gene calling. [default: bac120]
	Available HMMs can be found in '$HMM_DIR'. The HMM sets are downloaded separately from the folder
	'HMM' of https://github.com/zdf1987/EasyCGTree5 or https://gitee.com/zdf1987/EasyCGTree5 (see HMM/README.txt).
-seq <String, 'prot' or 'nucl'>
	Sequence type (amino acid or nucleotide) used for tree inference. [default: prot]
	'nucl' requires the whole input genome dataset (-input) to be DNA sequences (genomes or NCBI CDS files).
-task <String, 'all', 'predict', 'hmmsearch', 'refine', 'alignment', 'tree_infer'>
	Set run mode. [default: all]
	predict:    CDS prediction (Prodigal) of the genomes provided as DNA sequences.
	hmmsearch:  HMM search. Reuses an existing prediction, so several HMM sets can be tried without predicting again
	            (the prediction is run automatically if none exists yet).
	refine:     screening of the HMM search results.
	alignment:  retrieval, alignment and trimming of the gene clusters.
	tree_infer: tree inference.
-tree <String, 'sm', 'st', 'cs'>
	Approach used for tree inference. [default: sm]
	sm, supermatrix; st, supertree; cs, consensus tree
-tree_app <String, 'fasttree', 'iqtree'>
	Application used for tree inference. [default: fasttree]
-trim <String, 'gappyout', 'strict', 'strictplus', 'nogaps'>
	Standard for trimming alignments used by trimAl. [default: strict]
	'nogaps' removes every column that contains a gap ('nogap' is accepted, too).
-thread <Int>
	Number of threads to be used by HMMER, MUSCLE, FastTree, IQ-TREE and ASTRAL. [default: 4]


);
my ($inputDir);


########## Default parameters ################
my $task="all";
my $evalue=1e-10;
my $geneCutoff=0.8;
my $hmm = "bac120";
my $trim="strict";
my $genomeCutoff=0.8;
my $thread=4;
my $tree_type ="sm";
my $tree_app = "fasttree";
my $seq_type = "prot";
my $cds_nucl = "no";
##########################################

my %opt=qw();
GetOptions(\%opt,"input:s","task:s","tree:s","tree_app:s","hmm:s","trim:s","thread:i","evalue:f","gene_cutoff:f", "seq:s", "genome_cutoff:f","help!", "options_json!","keep_CDS_nucl!","dry_run!","outdir:s") or die "@usage\n\n\nERROR: Unrecognized option (see above).\n\n";

if ($opt{"options_json"}) {
	&printOptionsJson();
	exit 0;
}
print "\n##############################\nReading command line... checking options and input data...\n\n";
if (scalar(keys %opt )==0 || exists($opt{"help"})) {
	print join("\n",@usage)."\n\n";
	exit;
}

$cds_nucl = "yes" if exists($opt{"keep_CDS_nucl"});

unless (exists($opt{"input"}) && defined $opt{"input"} && $opt{"input"} =~ /\S/) {
	print join("\n",@usage)."\n\nERROR: a directory must be specified for option '-input' to start.\n\n";
	exit 1;
}

my @tasks= qw(all predict hmmsearch refine alignment tree_infer);
if (exists($opt{"task"})) {
	$task = $opt{"task"};
	unless (grep { $_ eq $task } @tasks) {
		print join("\n",@usage)."\nERROR: Argument '-task $task': the value '$task' should be one of the following: ".join(", ", @tasks).".\n\n";
		exit 1;
	}
}

########## Input directory, output directory and working (TEM) directory ##########
# The input directory may be anywhere; all outputs go to the directory that contains it.
my $inputArg = $opt{"input"};
$inputArg =~ s/[\/\\]+$// unless $inputArg =~ /^[\/\\]+$/;
$inputDir = File::Spec->rel2abs($inputArg);
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
my $TEMdir    = "$outDir/${inputName}_TEM";
my $outPrefix = "$outDir/$inputName";
my $RECORD    = "$TEMdir/EasyCGTree_record.txt";
my $TYPEFILE  = "$TEMdir/Genome_SeqType.txt";

if (-e $inputDir && !-d $inputDir) {
	print join("\n",@usage)."\n\nERROR: Argument '-input $inputArg': The value should be a directory, not a file.\n\n";
	exit 1;
}

########## Dry run: show how the input files would be prepared, then stop ##########
if ($opt{"dry_run"}) {
	unless (-d $inputDir) {
		print "\nERROR: Argument '-input $inputArg': directory '$inputDir' does not exist.\n\n";
		exit 1;
	}
	$DRY_RUN = 1;
	my ($units, $notes) = &planInputs($inputDir);
	my $sq = (defined $opt{"seq"} && $opt{"seq"} eq "nucl") ? "nucl" : "prot";
	my @warn;
	push @warn, "Fewer than 5 genomes (".scalar(@$units)."); EasyCGTree needs at least 5." if @$units < 5;
	my @prot = map { $_->{"name"} } grep { $_->{"type"} eq "prot" } @$units;
	push @warn, "'-seq nucl' needs DNA sequences, but ".scalar(@prot)." genome(s) are protein sequences: ".join(", ", @prot) if $sq eq "nucl" && @prot;
	my %tn = ("nucl" => "genome (DNA): CDS predicted by Prodigal", "prot" => "protein sequences", "ncbi_pair" => "NCBI protein + CDS, paired by locus_tag", "ncbi_cds" => "NCBI CDS, translated");
	print "\nDry run: '$inputDir' contains ".scalar(@$units)." genomes.\n\nGenome\tType\tFiles\tIgnored_files\n";
	my @ju;
	foreach my $u (@$units) {
		my @files = @{ $u->{"orig"} || $u->{"source"} };
		print join("\t", $u->{"name"}, $u->{"type"}, join(",", @files), (@{ $u->{"ignored"} } ? join(",", @{ $u->{"ignored"} }) : "-")), "\n";
		push @ju, { "name" => $u->{"name"}, "type" => $u->{"type"}, "type_text" => $tn{$u->{"type"}}, "files" => \@files, "ignored" => $u->{"ignored"} };
	}
	print "\nTypes: ".join("; ", map { "$_ = $tn{$_}" } sort keys %tn)."\n";
	print "NOTE: $_\n" foreach @$notes;
	print "WARNING: $_\n" foreach @warn;
	require JSON::PP;
	print "\nEasyCGTree_DRY_RUN_JSON=".JSON::PP->new->canonical(1)->encode({ "input" => $inputDir, "genomes" => \@ju, "notes" => $notes, "warnings" => \@warn })."\n";
	exit 0;
}

# Is there already a usable CDS prediction in the working directory?
my $havePrediction = (-d "$TEMdir/TEM0_CDS" && -f $TYPEFILE) ? 1 : 0;
my $doPredict = 0;
if ($task eq "all" || $task eq "predict") {
	$doPredict = 1;
} elsif ($task eq "hmmsearch" && !$havePrediction) {
	$doPredict = 1;
}

my $fn1 = 0;
if ($doPredict) {
	unless (-d $inputDir) {
		print join("\n",@usage)."\n\nERROR: Argument '-input $inputArg': directory '$inputDir' does not exist.\n\n";
		exit 1;
	}
	# The input files are classified, grouped and counted in Step 0 (see sub planInputs).
} elsif (!-d $TEMdir) {
	print join("\n",@usage)."\n\nERROR: '-task $task' needs the results of the previous task(s) in '$TEMdir', but this directory does not exist.\nPlease run the previous task(s) first, or use '-task all'".((defined $opt{"outdir"} && $opt{"outdir"} ne "") ? "" : " (or give the output directory of the previous run with '-outdir')").".\n\n";
	exit 1;
}
unless (-d $outDir) {
	mkpath($outDir);
	unless (-d $outDir) {
		print "\nERROR: the output directory '$outDir' can't be created.\n\n";
		exit 1;
	}
}
$outDir = abs_path($outDir);
$TEMdir = "$outDir/${inputName}_TEM";
$outPrefix = "$outDir/$inputName";
$RECORD = "$TEMdir/EasyCGTree_record.txt";
$TYPEFILE = "$TEMdir/Genome_SeqType.txt";

$evalue=$opt{"evalue"} if exists($opt{"evalue"});

my $treesig;
if (exists($opt{"tree"})) {
	$tree_type = $opt{"tree"};
	if ($tree_type eq "cs") {
		$treesig=1;
	} elsif ($tree_type eq "sm" || $tree_type eq "st") {
	} else {
		print join("\n",@usage)."\nERROR: Argument '-tree $tree_type': the value '$tree_type' cannot be recognized (requested to be 'st', 'cs', or 'sm').\n\n";
		exit 1;
	}
}

if (exists($opt{"trim"})) {
	$trim = $opt{"trim"};
	$trim = "nogaps" if $trim eq "nogap";
	unless ($trim eq "gappyout" || $trim eq "strict" || $trim eq "strictplus" || $trim eq "nogaps") {
		print join("\n",@usage)."\nERROR: Argument '-trim $trim': the value '$trim' cannot be recognized (requested to be 'gappyout', 'strict', 'strictplus', or 'nogaps').\n\n";
		exit 1;
	}
}

if (exists($opt{"tree_app"})) {
	$tree_app = $opt{"tree_app"};
	unless ($tree_app eq "fasttree" || $tree_app eq "iqtree") {
		print join("\n",@usage)."\nERROR: Argument '-tree_app $tree_app': the value '$tree_app' cannot be recognized (requested to be 'fasttree' or 'iqtree').\n\n";
		exit 1;
	}
}

########## Profile HMM: a name in the HMM directory, or the path of a .hmm file ##########
my $hmmGiven = exists($opt{"hmm"}) ? 1 : 0;
my $hmmFile;
if ($hmmGiven) {
	my $hmm0 = $opt{"hmm"};
	if (-f $hmm0) {
		$hmmFile = abs_path($hmm0);
		$hmm = basename($hmm0);
	} else {
		$hmm = $hmm0;
	}
	$hmm =~ s/\.hmm$//i;
	$hmmFile = "$HMM_DIR/$hmm.hmm" unless $hmmFile;
} else {
	$hmmFile = "$HMM_DIR/$hmm.hmm";
}
if (($task eq "all" || $task eq "hmmsearch") && !-f $hmmFile) {
	print "\nERROR: ".&hmmMissingMsg($hmm, $HMM_DIR)."\n\n";
	exit 1;
}

if (exists($opt{"genome_cutoff"})) {
	my $qwe=$opt{"genome_cutoff"};
	unless ($qwe =~ /^(0?\.\d+|0|1|1\.0*)$/ && $qwe >= 0 && $qwe <= 1) {
		print join("\n",@usage)."\nERROR: Argument '-genome_cutoff $qwe': the value '$qwe' should be a decimal between 0 and 1 (recommended 0.5-1).\n\n";
		exit 1;
	}
	$genomeCutoff = $qwe;
}

if (exists($opt{"thread"})) {
	my $BNM=$opt{"thread"};
	unless (defined $BNM && $BNM =~ /^\d+$/ && $BNM >= 1) {
		print join("\n",@usage)."\nERROR: Argument '-thread': The value should be an integer >= 1.\n\n";
		exit 1;
	}
	$thread=$BNM;
}
my @notes;    # remarks printed with the starting information
{
	my $ncpu = &cpuCount();
	if ($ncpu && $thread > $ncpu) {
		push @notes, "WARNING: '-thread $thread' is larger than the number of CPU cores of this computer ($ncpu); $ncpu threads are used." if exists $opt{"thread"};
		$thread = $ncpu;
	}
}

if (exists($opt{"seq"})) {
	my $seqt=$opt{"seq"};
	unless (defined $seqt && ($seqt eq "prot" || $seqt eq "nucl")) {
		print join("\n",@usage)."\nERROR: Argument '-seq': The value should be 'prot' or 'nucl'.\n\n";
		exit 1;
	}
	$seq_type = $seqt;
	$cds_nucl = "yes" if $seqt eq "nucl";
}

if (exists($opt{"gene_cutoff"}))	 {
	my $rewq=$opt{"gene_cutoff"};
	unless ($rewq =~ /^(0?\.\d+|0|1|1\.0*)$/ && $rewq >= 0 && $rewq <= 1) {
		print join("\n",@usage)."\nERROR: Argument '-gene_cutoff $rewq': the value '$rewq' should be a decimal between 0 and 1 (recommended 0.5-1).\n\n";
		exit 1;
	}
	$geneCutoff= $rewq;
}

########## Record of previous tasks (stored in the TEM directory) ##########
my %rec;
%rec = &readRecord() if -f $RECORD && !($task eq "all" || $task eq "predict");

# Later tasks must use the same profile HMM as the existing HMM search results.
if ($task eq "refine" || $task eq "alignment" || $task eq "tree_infer") {
	if (defined $rec{"hmm"}) {
		if ($hmmGiven) {
			if ($hmm ne $rec{"hmm"}) {
				print "\nERROR: Argument '-hmm $opt{'hmm'}': the existing HMM search results in '$TEMdir' were obtained with the profile HMM '$rec{'hmm'}', not '$hmm'.\n",
					"Please use '-hmm $rec{'hmm'}' (or leave out '-hmm'), or run '-task hmmsearch -hmm $hmm' first.\n\n";
				exit 1;
			}
			if (-f $hmmFile && $rec{"hmm_md5"} && &md5File($hmmFile) ne $rec{"hmm_md5"}) {
				print "\nERROR: Argument '-hmm $opt{'hmm'}': the file '$hmmFile' has changed since the HMM search was performed with it.\n",
					"Please run '-task hmmsearch -hmm $hmm' again.\n\n";
				exit 1;
			}
		} else {
			$hmm = $rec{"hmm"};
			push @notes, "The profile HMM '$hmm' recorded in '$RECORD' is used.";
		}
	} else {
		push @notes, "WARNING: no record of the HMM search was found in '$TEMdir' (it may have been created by EasyCGTree < 5.0). The profile HMM cannot be checked, '$hmm' is assumed.";
	}
	# Hits worse than the E-value used by hmmsearch were never reported, so a larger value has no effect.
	if ($task eq "refine" && defined $rec{"evalue"}) {
		if (!exists($opt{"evalue"})) {
			$evalue = $rec{"evalue"};
		} elsif ($evalue > $rec{"evalue"}) {
			push @notes, "WARNING: '-evalue $evalue' is larger than the E-value used by hmmsearch ($rec{'evalue'}); hits between the two values are not available.";
		}
	}
}

my $logFile = ($task eq "predict") ? "$outPrefix.predict_$STAMP.log" : "$outPrefix.$hmm.$tree_type.$tree_app\_$STAMP.log";
open (LOG, ">", $logFile) or die "Can't open '$logFile': $!\n";
# Every error message also goes to the log file.
$SIG{__DIE__} = sub { print LOG "\nERROR: $_[0]" if defined fileno(LOG); };

my $order;
my @nouseop;
my %taskOptions = (
	"predict"    => [ [qw(seq keep_CDS_nucl)],                  [qw(hmm thread evalue tree tree_app gene_cutoff genome_cutoff trim)] ],
	"hmmsearch"  => [ [qw(hmm thread evalue seq keep_CDS_nucl)], [qw(tree tree_app gene_cutoff genome_cutoff trim)] ],
	"refine"     => [ [qw(gene_cutoff genome_cutoff evalue tree)], [qw(thread tree_app trim seq keep_CDS_nucl)] ],
	"alignment"  => [ [qw(thread trim seq)],                     [qw(evalue tree tree_app gene_cutoff genome_cutoff keep_CDS_nucl)] ],
	"tree_infer" => [ [qw(thread tree tree_app)],                [qw(evalue gene_cutoff genome_cutoff trim seq keep_CDS_nucl)] ],
);
if ($task eq "all") {
	$order= "The user ordered a complete analysis (-task all) of the pipeline!\n\n";
} elsif ($task eq "predict") {
	$order= "The user ordered a CDS prediction task (-task predict) of the pipeline!\n\n";
} elsif ($task eq "hmmsearch") {
	$order= "The user ordered a HMM search task (-task hmmsearch) of the pipeline!\n\n";
} elsif ($task eq "refine") {
	$order= "The user ordered a task (-task refine) of screening significant homologs from the pre-existing HMM searching result and generating gene clusters!\n\n";
} elsif ($task eq "alignment") {
	$order= "The user ordered a task (-task alignment) of aligning sequences of each gene cluster!\n\n";
} elsif ($task eq "tree_infer") {
	$order= "The user ordered a task (-task tree_infer) of inferring phylogenetic tree from pre-existing alignment(s)!\n\n";
}
if ($taskOptions{$task}) {
	my %noneed = map { $_ => 1 } @{ $taskOptions{$task}[1] };
	# An existing prediction is reused by 'hmmsearch', so the prediction options have no effect.
	if ($task eq "hmmsearch" && !$doPredict) {
		$noneed{"seq"}=1;
		$noneed{"keep_CDS_nucl"}=1;
	}
	foreach my $op (sort keys %opt) {
		push @nouseop, $op if $noneed{$op};
	}
}

if (@nouseop) {
	&msg("\nWARNING: '-task $task' does not require the following settings, which will be ignored in current analysis:");
	foreach my $zz (@nouseop) {
		my $val = (defined $opt{$zz} && $opt{$zz} ne "1") ? " $opt{$zz}" : "";
		$val = "" if $zz eq "keep_CDS_nucl";
		&msg(" -$zz$val");
	}
	&msg(".\n");
}

if ($treesig) {
	&msg("\nWARNING: '-tree $tree_type' requires '-gene_cutoff' to be set as '1'.\nThe value '1' will be used instead.\n") if $geneCutoff <1 && ($task eq "all" || $task eq "refine");
	$geneCutoff=1;
}


&msg("\n\n############### Starting Information ###############\n\n");
&msg("====== EasyCGTree ======\n	Version $VERSION\nby Dao-Feng Zhang\n\n");
&msg("#########################\n\nOptions: \n");
&msg("-input $inputDir\n-task $task\n-tree $tree_type\n-tree_app $tree_app\n-trim $trim\n-hmm $hmm\n-thread $thread\n-evalue $evalue\n-gene_cutoff $geneCutoff\n-genome_cutoff $genomeCutoff\n-seq $seq_type\n");
&msg("-keep_CDS_nucl\n") if $cds_nucl eq "yes" && $seq_type ne "nucl";
&msg("\nOutput directory: $outDir\nWorking directory: $TEMdir\n");
foreach my $nn (@notes) {
	&msg("\n$nn\n");
}
&msg("\n#########################\n\nJob Started at: $START_STR\n\n*****$order\n\n");


############### Task 0: CDS prediction #################
if ($doPredict) {

	if ($task eq "hmmsearch") {
		&msg("NOTE: no CDS prediction was found in '$TEMdir'; it will be performed first.\n\n");
	}

#============= Step 0: CDS prediction ==============#

	&msg("############### Task 0: CDS prediction ###############\n#============= Step 0: CDS prediction ==============# \n\nStep 0: Preparation of proteome from the genome data in directory '$inputDir': \n\n");

	# Classify all input files and decide how each genome is prepared, before doing anything else
	# (problems in the input are reported before the old working directory is removed).
	my ($units, $planNotes) = &planInputs($inputDir);
	$fn1 = scalar @$units;
	unless ($fn1 >=5) {
		&fail("Number of genomes (taxa on a tree) in input directory '$inputDir' ($fn1) must be >=5 to start a run.");
	}
	my @protIn = map { $_->{"source"}[0] } grep { $_->{"type"} eq "prot" } @$units;
	if ($seq_type eq "nucl" && @protIn) {
		&fail("Argument '-seq nucl' requires that all the input genome data are provided as nucleotide/DNA sequences.\nHOWEVER, the following ".scalar(@protIn)." file(s) in '$inputDir' contain protein sequences:\n	".join("\n	", @protIn)."\nPlease replace them with their genome (DNA) sequences or NCBI CDS sequences (cds_from_genomic.fna), remove them, or use '-seq prot'.");
	}
	&removePath($TEMdir) if -e $TEMdir;
	&makeDir($TEMdir);
	&makeDir("$TEMdir/TEM0_CDS");
	&writeRecord("version" => $VERSION);
	foreach my $nn (@$planNotes) {
		&msg("NOTE: $nn\n");
	}
	&msg("\n") if @$planNotes;

	my $prodigal;
	$prodigal = &tool("prodigal") if grep { $_->{"type"} eq "nucl" } @$units;
	my $genomeNum0=0;
	my %typeCount;
	open (TYP, ">", $TYPEFILE) or die "Step 0: Can't open '$TYPEFILE': $!\n";
	print TYP "Genome\tInput_type\tProtein_seq\tCDS_DNA_seq\tSource_files\tIgnored_files\n";
	foreach my $unit (@$units) {
		$genomeNum0++;
		my $rat= &givePercentage($genomeNum0,$fn1);
		my $onam = $unit->{"name"};
		my $type = $unit->{"type"};
		my @src  = map { "$inputDir/$_" } @{ $unit->{"source"} };
		my $faa = "$TEMdir/TEM0_CDS/$onam.faa";
		my $fnn = "$TEMdir/TEM0_CDS/$onam.fnn";
		if ($type eq "nucl") {
			print ("\r","Predicting coding sequences (CDS) from '$src[0]': $genomeNum0/$fn1 ($rat).  ");
			# A normalized copy (upper case, no '\r') is given to Prodigal.
			my $tem = "$TEMdir/tem_genome.fna";
			&copyFastaUpper($src[0], $tem);
			my @cmd0 = ($prodigal, "-i", $tem, "-o", "$TEMdir/tem.geneinfo", "-m", "-q", "-a", $faa);
			push @cmd0, ("-d", $fnn) if $cds_nucl eq "yes";
			my $rc = &runCmd(@cmd0);
			unless ($rc == 0 && -s $faa) {
				&fail("Step 0: Prodigal failed to predict CDS from '$src[0]' (exit code $rc).\nPlease make sure that 'prodigal$EXE' is present in '$BIN_DIR', and that the file is a complete genome\n(Prodigal needs at least 20,000 bp of sequence in its normal mode).");
			}
		} elsif ($type eq "prot") {
			print ("\r","Preparing from '$src[0]': $genomeNum0/$fn1 ($rat).  ");
			&copyFastaUpper($src[0], $faa);
		} elsif ($type eq "ncbi_pair") {
			print ("\r","Pairing NCBI protein and CDS sequences by locus_tag: '$onam' $genomeNum0/$fn1 ($rat).  ");
			my ($n, $total) = &ncbiPair($src[0], $src[1], $faa, $fnn);
			$unit->{"info"} = "$n of $total proteins paired";
		} elsif ($type eq "ncbi_cds") {
			print ("\r","Translating NCBI CDS sequences: '$onam' $genomeNum0/$fn1 ($rat).  ");
			my ($n, $total) = &ncbiTranslate($src[0], $faa, $fnn);
			$unit->{"info"} = "$n of $total CDS translated";
		}
		$typeCount{$type}++;
		print TYP join("\t", $onam, $type, "$onam.faa", (-s $fnn ? "$onam.fnn" : "-"),
			join(",", @{ $unit->{"source"} }), (@{ $unit->{"ignored"} } ? join(",", @{ $unit->{"ignored"} }) : "-")), "\n";
	}
	close TYP;
	unlink "$TEMdir/tem_genome.fna", "$TEMdir/tem.geneinfo";
	&msg("\nStep 0: $fn1 genomes were prepared in '$TEMdir/TEM0_CDS':\n");
	&msg("        $typeCount{'nucl'} genomes were provided as DNA sequences, and their CDS have been predicted by Prodigal.\n") if $typeCount{"nucl"};
	&msg("        DNA sequences of the CDS predicted by Prodigal were kept (.fnn files).\n") if $cds_nucl eq "yes" && $typeCount{"nucl"};
	&msg("        $typeCount{'prot'} genomes were provided as protein sequences and have been used directly.\n") if $typeCount{"prot"};
	&msg("        $typeCount{'ncbi_pair'} genomes were provided as NCBI protein + CDS sequences, paired by locus_tag.\n") if $typeCount{"ncbi_pair"};
	&msg("        $typeCount{'ncbi_cds'} genomes were provided as NCBI CDS sequences and have been translated.\n") if $typeCount{"ncbi_cds"};
	foreach my $unit (@$units) {
		print LOG "        $unit->{'name'}: $unit->{'info'}\n" if $unit->{"info"};
	}
	&msg("        The sequence type and source files of each genome are listed in '$TYPEFILE'.\n\n");
	&writeRecord("predict_done" => &now(), "predict_cds_nucl" => $cds_nucl, "genomes" => $fn1);
}


############### Task 1: HMM search #################
if ($task eq "all" || $task eq "hmmsearch") {

	unless (-d "$TEMdir/TEM0_CDS" && -f $TYPEFILE) {
		&fail("Step 1: The CDS prediction results in '$TEMdir/TEM0_CDS' are missing. Please run '-task predict' first.");
	}
	# Results of an earlier HMM search (maybe with another HMM) and everything derived from them are removed.
	&cleanAfterPrediction();

#============= Step 1: HMM search ==============#
	&msg("\n############### Task 1: HMM search ###############\n#============= Step 1: HMM search ==============# \n");
	&makeDir("$TEMdir/TEM1_HMMsearch_out");
	my $genomeNum = 0;

	open(HMM, "<", $hmmFile) or die "Step 1: Can't open '$hmmFile': $!\n";
	my (%hmm,$name);
	foreach my $hmmin (<HMM>) {
		$hmmin=~ s/\r?\n$//;
		next unless $hmmin;
		if ($hmmin =~ s/^NAME\s+//) {
			$hmm{$hmmin}="$hmmin	Null	MO0W2E1XX";
			$name=$hmmin;
		}  elsif ($hmmin =~ s/^ACC\s+//) {
			$hmm{$name}=~ s/Null/$hmmin/;
		}  elsif ($hmmin =~ s/^LENG\s+//) {
			$hmm{$name}=~ s/MO0W2E1XX/$hmmin/;
		}
	}
	close HMM;
	&fail("Step 1: No profile HMM was found in '$hmmFile'.") unless %hmm;
	open(MMIN, ">", "$TEMdir/HMMinfo.txt") or die "Step 1: Can't open '$TEMdir/HMMinfo.txt': $!\n";
	print MMIN "NAME	ACC	Length\n";
	foreach my $kku (sort keys %hmm) {
		print MMIN "$hmm{$kku}\n";
	}
	close MMIN;

	opendir(DIR1, "$TEMdir/TEM0_CDS") || die "Step 1: Can't open input directory '$TEMdir/TEM0_CDS': $!\n";
	my @query1 = sort grep { /\.faa$/ } readdir(DIR1);
	closedir(DIR1);
	&fail("Step 1: Directory '$TEMdir/TEM0_CDS' does not contain any protein (.faa) files. Please run '-task predict' again.") unless @query1;
	my $nfaa = scalar @query1;

	# Warn when the input directory has been changed since the prediction.
	if (!$doPredict && -d $inputDir) {
		my %used;
		if (open(my $th, "<", $TYPEFILE)) {
			while (my $l = <$th>) {
				$l =~ s/\r?\n$//;
				my @f = split /\t/, $l;
				next if !@f || $f[0] eq "Genome";
				if (defined $f[4]) {
					$used{$_} = 1 foreach (split /,/, $f[4]), (split /,/, ($f[5] || ""));
				} else {
					$used{"$f[0].fas"} = 1;    # made by an earlier 5.0 build without source columns
				}
			}
			close $th;
		}
		my @newIn = grep { !$used{$_} } &listInputFiles($inputDir);
		&msg("\nWARNING: ".scalar(@newIn)." file(s) in '$inputDir' were not included in the existing CDS prediction (e.g. '$newIn[0]').\nRun '-task predict' (or '-task all') if they should be analysed.\n\n") if @newIn;
	}

	&msg("Step 1: HMM searching against $nfaa proteomes in directory '$TEMdir/TEM0_CDS' with '$hmmFile': \n\n");
	my $hmmsearch = &tool("hmmsearch");
	foreach my $file (@query1) {
		$genomeNum++;
		my $onam=$file;
		$onam=~ s/\.faa$/\.fas/;
		my $rat= &givePercentage($genomeNum,$nfaa);
		print ("\r","	Searching against '$file': $genomeNum/$nfaa ($rat).           ");
		my $rc = &runCmd($hmmsearch, "--tblout", "$TEMdir/TEM1_HMMsearch_out/$onam", "-E", $evalue, "--cpu", $thread, "-o", File::Spec->devnull(), $hmmFile, "$TEMdir/TEM0_CDS/$file");
		unless ($rc == 0 && -e "$TEMdir/TEM1_HMMsearch_out/$onam") {
			&fail("Step 1: hmmsearch failed on '$TEMdir/TEM0_CDS/$file' (exit code $rc).\nPlease make sure that 'hmmsearch$EXE' is present in '$BIN_DIR', and that the HMM file '$hmmFile' is correct.");
		}
	}
	print "\n\n";
	&msg("=====Results of $genomeNum HMM searches were written in directory '$TEMdir/TEM1_HMMsearch_out/'!=====\n");
	&writeRecord("hmm" => $hmm, "hmm_file" => $hmmFile, "hmm_md5" => &md5File($hmmFile), "evalue" => $evalue, "hmmsearch_done" => &now());
}


############### Task 2: Filtrate HMM Search Results #################
#============= Step 2: Filtrate HMM Search Results ==============#
if ($task eq "all" || $task eq "refine") {

	&msg("\n\n############### Task 2: Filtration HMM Search Results ###############\n#============= Step 2: Filtration HMM Search Results ==============# \n");

	unless (-d "$TEMdir/TEM1_HMMsearch_out" && -f "$TEMdir/HMMinfo.txt") {
		&fail("Step 2: The HMM search results ('$TEMdir/TEM1_HMMsearch_out' and '$TEMdir/HMMinfo.txt') are missing. Please run '-task hmmsearch' first.");
	}
	&removePath("$TEMdir/TEM2_HMMsearch_outS") if -e "$TEMdir/TEM2_HMMsearch_outS";
	&makeDir("$TEMdir/TEM2_HMMsearch_outS");

	opendir(DIR3, "$TEMdir/TEM1_HMMsearch_out") || die "Step 2: Can't open input directory '$TEMdir/TEM1_HMMsearch_out': $!\n";
	my @query2 = sort grep { /\.fas$/ } readdir(DIR3);
	closedir(DIR3);
	&fail("Step 2: Input directory '$TEMdir/TEM1_HMMsearch_out' does not contain any files.") unless @query2;

	open(MMIN, "<", "$TEMdir/HMMinfo.txt") or die "Step 2: Can't open '$TEMdir/HMMinfo.txt': $!\n";
	my %score;
	foreach my $hmmin (<MMIN>) {
		$hmmin=~ s/\r?\n$//;
		next unless $hmmin;
		my @hmmi =split /\t/, $hmmin;
		$score{$hmmi[0]} = $hmmi[2];
	}
	close MMIN;

	foreach my $file3 (@query2) {
	  if ($geneCutoff >=0.5) {################################################
		my (%out3);
		open (FIL3, "<", "$TEMdir/TEM1_HMMsearch_out/$file3") or die "Can't open '$TEMdir/TEM1_HMMsearch_out/$file3': $!\n";
		my $geneC=0;
		my $head="";
		foreach my $in3 (<FIL3>) {
			$in3=~ s/\r?\n$//;
			next unless $in3;
			if ($in3 =~ /^\#/) {
				$head .= "$in3\n";
				next;
			}
			my @in3 =split /\s+/, $in3;
			next if $in3[4] > $evalue;
			my $sco = ($score{$in3[2]} || 0)/4;
			next if $in3[5] < $sco;
			my $keyr = ($in3[3] =~ /-/) ? $in3[2] : $in3[3];
			unless ($out3{$keyr}) {
				$out3{$keyr}=$in3;
				$geneC++;
			}
		}
		close FIL3;
		my $outfile3="$geneC"."__$file3";
		open(OU3, ">", "$TEMdir/TEM2_HMMsearch_outS/$outfile3") or die "Step 2: Can't open '$TEMdir/TEM2_HMMsearch_outS/$outfile3': $!\n";
		print OU3 "$head";
		foreach my $in (sort keys %out3) {
			print OU3 "$out3{$in}\n";
		}
		close OU3;
	  } else {####################################################################
		open (FIL3, "<", "$TEMdir/TEM1_HMMsearch_out/$file3") or die "Can't open '$TEMdir/TEM1_HMMsearch_out/$file3': $!\n";
		my (%hmmSin,%hmmE,%hmmSinO);
		my $geneC=0;
		my $head="";
		foreach my $in11 (<FIL3>) {
			$in11=~ s/\r?\n$//;
			next unless $in11;
			if ($in11 =~ /^\#/) {
				$head .= "$in11\n";
				next;
			}
			my @in3 =split /\s+/, $in11;
			next if $in3[4] > $evalue;
			my $sco = ($score{$in3[2]} || 0)/4;
			next if $in3[5] < $sco;
			my $keyr=$in3[0];
			if (!$hmmSin{$keyr} || $hmmE{$keyr}>$in3[4]) {
				$hmmSin{$keyr}=$in11;
				$hmmE{$keyr}=$in3[4];
			}
		}
		close FIL3;
		my %ee;
		foreach my $ch (sort keys %hmmE) {
			my @in3 =split /\s+/, $hmmSin{$ch};
			my $keyr=$in3[2];
			unless ($hmmSinO{$keyr}) {
				$geneC++;
				$hmmSinO{$keyr}=$hmmSin{$ch};
				$ee{$keyr}=$in3[4];
			} elsif ($ee{$keyr} > $in3[4]) {
				$hmmSinO{$keyr} =$hmmSin{$ch};
				$ee{$keyr} = $in3[4];
			}
		}
		my $outfile3="$geneC"."__$file3";
		open(OU3, ">", "$TEMdir/TEM2_HMMsearch_outS/$outfile3") or die "Step 2: Can't open '$TEMdir/TEM2_HMMsearch_outS/$outfile3': $!\n";
		print OU3 "$head";
		foreach my $in (sort keys %hmmSinO) {
			print OU3 "$hmmSinO{$in}\n";
		}
		close OU3;
	  }########################################################################################
	}
	&msg("=====Screened HMM search results have been written to '$TEMdir/TEM2_HMMsearch_outS'!=====\n");

###################################################
	opendir(DIR, "$TEMdir/TEM2_HMMsearch_outS") || die "Step 2: Can't open input directory '$TEMdir/TEM2_HMMsearch_outS': $!\n";
	my @query3 = sort grep { /.\.fas$/ } readdir(DIR);
	closedir(DIR);
	&fail("Task 2: Directory '$TEMdir/TEM2_HMMsearch_outS' does not contain any files.") unless @query3;

	my (@querylist0,@querylistN,%geneIndex);
	my (@finalgenelist,@outgenelist,@finalgenome,@outgenome);
	my (@finalfile);
	my $geneNum=0;
	my $genomeNum=0;

	foreach my $kk (@query3) {
		$genomeNum++;
		my @yy=split /__/,$kk;
		$geneNum = $yy[0] if $geneNum < $yy[0];
	}
	my $cutoff2=int($geneNum*$genomeCutoff);
	&msg("\n These genomes harboring more than $cutoff2 ($geneNum * $genomeCutoff) genes, which will be used in following analysis. This is the list (genome: gene_number):\n");

	foreach my $kk (@query3) {
		my @yy=split /__/,$kk;
		if ($yy[0]>=  $cutoff2) {
			push @finalfile, $kk ;
			push @finalgenome, "$yy[1]:$yy[0]";
			&msg("$yy[1]	$yy[0]\n");
		} else {
			push @outgenome, "$yy[1]	$yy[0]\n";
		}
	}

	my $dis=scalar @outgenome;
	my $ingenomeN=scalar @finalgenome;

	&msg(".\n=====In total, $ingenomeN of $genomeNum genomes were selected!=====\n\n");
	if (@outgenome) {
		&msg("The following $dis genomes were excluded (genome/gene_number): @outgenome.\n\n");
	}
	&fail("Step 2: Fewer than 4 genomes passed the genome cutoff ($ingenomeN), so no tree can be inferred. Please check the input genomes or lower '-genome_cutoff'.") if $ingenomeN < 4;

	foreach my $file3 (@finalfile) {
		open (FIL3, "<", "$TEMdir/TEM2_HMMsearch_outS/$file3") or die "Step 2: Can't open '$TEMdir/TEM2_HMMsearch_outS/$file3': $!\n";
		foreach my $in3 (<FIL3>) {
			$in3=~ s/\r?\n$//;
			next unless $in3;
			next if $in3 =~ /^\#/;
			my @in3 =split /\s+/, $in3;
			my $keyr = ($in3[3] =~ /-/) ? $in3[2] : $in3[3];
			unless (exists $geneIndex{$keyr}) {
				$geneIndex{$keyr} = scalar @querylist0;
				push @querylist0,$keyr;
				push @querylistN, 0;
			}
			$querylistN[$geneIndex{$keyr}]++;
		}
		close FIL3;
	}
	my $cutoff=int($ingenomeN*$geneCutoff);

	&msg("\nThese genes were present in >= $cutoff ($ingenomeN * $geneCutoff) of the $ingenomeN selected genomes, which will be used for tree inference (gene/prevalence):\n");

	foreach my $kk (0..$#querylistN) {
		if ($querylistN[$kk]<  $cutoff) {
			push @outgenelist, "$querylist0[$kk]	$querylistN[$kk]\n";
		} else {
			&msg("$querylist0[$kk]	$querylistN[$kk]\n");
			push @finalgenelist, $querylist0[$kk];
		}
	}
	my $ingene =scalar @finalgenelist;
	my $outgene =scalar @outgenelist;

	&msg(".\n=====In total, $ingene of $geneNum genes were selected!=====\n");
	if (@outgenelist) {
		&msg("\nThe following $outgene/$geneNum genes were excluded (gene/prevalence):\n@outgenelist\n");
	}
	&fail("Step 2: No gene passed the gene cutoff. Please lower '-gene_cutoff' or use another HMM set.") unless @finalgenelist;
	open (TXT3, ">", "$TEMdir/GenomeGeneScreened.txt") or die "Step 2: Can't open '$TEMdir/GenomeGeneScreened.txt': $!\n";
	print TXT3 "Genome_List=@finalgenome\n\nGene_List=@finalgenelist\n";
	close TXT3;
	&msg("The genomes and genes used in following analysis were also listed in '$TEMdir/GenomeGeneScreened.txt'\n");
	&writeRecord("refine_done" => &now(), "gene_cutoff" => $geneCutoff, "genome_cutoff" => $genomeCutoff, "refine_evalue" => $evalue);
	&deleteRecord(qw(alignment_done alignment_seq trim));
}

############### Task 3: Multiple Sequences Alignment #################
#============= Step 3: Retrieve Sequences from Each proteome ==============#
if ($task eq "all" || $task eq "alignment") {

	&msg("\n\n############### Task 3: Multiple Sequences Alignment #################\n#============= Step 3: Retrieve Sequences from Each proteome ==============#\n\n");

	open (TXT4, "<", "$TEMdir/GenomeGeneScreened.txt") or &fail("Step 3: Can't open '$TEMdir/GenomeGeneScreened.txt': $!\nIf it is absent, you should run '-task refine' again.");
	my (@fgenome,@fgene);
	foreach my $k3 (<TXT4>) {
		$k3=~ s/\r?\n$//;
		next unless $k3;
		if ( $k3 =~ s/^Genome_List=//) {
			@fgenome = split /\s+/, $k3;
		} elsif( $k3 =~ s/^Gene_List=//) {
			@fgene = split /\s+/, $k3;
		}
	}
	close TXT4;
	my %fgene = map { $_ => 1 } @fgene;
	my $fg=scalar @fgenome;
	my $fgn=scalar @fgene;

	# '-seq nucl' needs the DNA sequences (.fnn) of the CDS of every selected genome.
	if ($seq_type eq "nucl") {
		my %types = &readTypes();
		my (@protGenome, @noFnn);
		foreach my $gg (@fgenome) {
			my ($g) = split /\:/, $gg;
			$g =~ s/\.fas$//;
			next if -s "$TEMdir/TEM0_CDS/$g.fnn";
			if (defined $types{$g} && $types{$g} eq "prot") {
				push @protGenome, $g;
			} else {
				push @noFnn, $g;
			}
		}
		if (@protGenome || @noFnn) {
			my $mm = "Step 3: '-seq nucl' requires the DNA sequences of the CDS ('.fnn' files in '$TEMdir/TEM0_CDS') of all the selected genomes, but:\n";
			$mm .= "- the following genome(s) were provided as protein sequences, so no DNA sequences are available:\n	".join("\n	", @protGenome)."\n  Please replace them with genome (DNA) sequences in '$inputDir' or remove them, then run from '-task predict' again, or use '-seq prot'.\n" if @protGenome;
			$mm .= "- the DNA sequences of the CDS of the following genome(s) were not kept during the prediction:\n	".join("\n	", @noFnn)."\n  Please run '-task predict -seq nucl' (or '-task predict -keep_CDS_nucl') and the following tasks again.\n" if @noFnn;
			&fail($mm);
		}
	}

	&msg("According to the record in '$TEMdir/GenomeGeneScreened.txt', $fgn common genes will be extracted from $fg genomes.\n");
	&msg("Sequence type: ".($seq_type eq "nucl" ? "nucleotide (CDS DNA sequences)" : "protein")."\n");

	foreach my $dd (qw(TEM3_GeneSeqs TEM4_GeneCluster TEM5_Alignment TEM6_AlnTrimmed)) {
		&removePath("$TEMdir/$dd") if -e "$TEMdir/$dd";
	}
	&makeDir("$TEMdir/TEM3_GeneSeqs");

	my $numg=0;
	foreach my $i1 (0..$#fgenome) {
		$numg++;
		my @hh=split /\:/, $fgenome[$i1];
		my $fgenome =$hh[0];
		my $onam=$fgenome;
		$onam=~ s/\.fas$//;
		$onam .= ($seq_type eq "nucl") ? ".fnn" : ".faa";
		my ($ids, $geno) = &readFasta("$TEMdir/TEM0_CDS/$onam");

		my $seqnum=0;
		open(OUT, ">", "$TEMdir/TEM3_GeneSeqs/$fgenome") or die "Step 3: Can't open '$TEMdir/TEM3_GeneSeqs/$fgenome': $!\n";
		my $hh = "$hh[1]__$hh[0]";
		open (FIL, "<", "$TEMdir/TEM2_HMMsearch_outS/$hh") or die "Step 3: Can't open '$TEMdir/TEM2_HMMsearch_outS/$hh': $!\n";
		my $seqTag="$hh[0]";
		$seqTag=~ s/\.fas$/_/;
		foreach my $i2 (<FIL>) {
			$i2=~ s/\r?\n$//;
			next unless $i2;
			next if $i2=~ /^\#/;
			my @in0 = split /\s+/, $i2;
			my $keyr = ($in0[3] =~ /-/) ? $in0[2] : $in0[3];
			my $geneSeq = $geno->{$in0[0]};
			unless (defined $geneSeq && length $geneSeq) {
				&msg("	WARNING: sequence '$in0[0]' ($keyr) was not found in '$TEMdir/TEM0_CDS/$onam'.\n");
				next;
			}
			$geneSeq =~ s/\*/-/g;
			print OUT ">$seqTag"."$keyr\n$geneSeq\n";
			$seqnum++;
		}
		close FIL;
		close OUT;
		print "	$numg: Got $seqnum/$fgn sequences from genome: $fgenome\n";
	}
	&msg("=====\nTotally, $numg sequence files have been written to '$TEMdir/TEM3_GeneSeqs'!=====\n");


#============= Step 4: Gather the Orthologs in Single File ==============#
	&msg("\n\n\n#============= Step 4: Gather the Orthologs in Single File ==============# \n\n");

	&makeDir("$TEMdir/TEM4_GeneCluster");
	my %outClu;
	my %cluCount;
	my $num5=0;
	foreach my $in5 (@fgenome) {
		my $tt =$in5;
		$tt =~ s/\:\d+$//;
		open (FILE, "<", "$TEMdir/TEM3_GeneSeqs/$tt") or die "Step 4: Can't open '$TEMdir/TEM3_GeneSeqs/$tt': $!\n";
		$tt =~ s/\.fas$//;
		my ($id,$sig);
		foreach my $i1 (<FILE>) {
			$i1=~ s/\r?\n$//;
			next unless $i1;
			if ( $i1 =~ /^\>/ ) {
				# Header: '>genome_gene'. Genome names never contain '_' (see subs planInputs and renameFiles),
				# so everything after the first '_' is the gene name, even if it contains '_'.
				(my $hd = $i1) =~ s/^\>//;
				(undef, $id) = split /_/, $hd, 2;
				$sig = (defined $id && $fgene{$id}) ? 1 : 0;
			} elsif ($sig) {
				$outClu{$id} .= ">$tt\n$i1\n";
				$cluCount{$id}++;
				$sig=0;
			}
		}
		close FILE;
		$num5++;
		my $rat =&givePercentage($num5, scalar @fgenome);
		print ("\r", "	Step 4: Reading the gene set ($rat) '$TEMdir/TEM3_GeneSeqs/$tt'. 		");
	}
	print "\n\n";
	my $n=0;
	my $cc=scalar keys %outClu;
	foreach my $in (sort keys %outClu) {
		$n++;
		open (OUT, ">", "$TEMdir/TEM4_GeneCluster/$in.fas") or die "Step 4: Can't open '$TEMdir/TEM4_GeneCluster/$in.fas': $!\n";
		print OUT "$outClu{$in}";
		close OUT;
		my $rat= &givePercentage($n, $cc);
		print ("\r", "	Step 4: Writing each gene cluster in single file. ($rat)		");
		print LOG "	Step 4: Writing $cluCount{$in} sequences of gene cluster '$in' in file '$TEMdir/TEM4_GeneCluster/$in.fas'\n";
	}
	print ".\n";
	print LOG "In total, $n gene clusters were written to '$TEMdir/TEM4_GeneCluster'.\n";

#============= Step 5: Multiple Sequence Alignment ==============#
	&msg("\n\n\n#============= Step 5: Multiple Sequence Alignment ==============# \n\n");

	&makeDir("$TEMdir/TEM5_Alignment");
	my $num7 = 0;
	my $muscle = &tool("muscle5");
	my @alnGenes;
	foreach (@fgene) {
		unless (-s "$TEMdir/TEM4_GeneCluster/$_.fas") {
			&msg("WARNING: Step 5: no sequence was retrieved for gene '$_'; it is skipped.\n");
			next;
		}
		my $outname = $_ . ".fas.fasta";
		my $an = ($fg < 200) ? "-align" : "-super5";
		my $rc = &runCmd($muscle, $an, "$TEMdir/TEM4_GeneCluster/$_.fas", "-output", "$TEMdir/TEM5_Alignment/$outname", "-threads", $thread);
		unless ($rc == 0 && -s "$TEMdir/TEM5_Alignment/$outname") {
			&fail("Step 5: MUSCLE failed to align '$TEMdir/TEM4_GeneCluster/$_.fas' (exit code $rc).".($rc >= 128 ? " The process was killed, probably because it ran out of memory." : "")."\nPlease make sure that 'muscle5$EXE' is present in '$BIN_DIR'.");
		}
		$num7++;
		push @alnGenes, $_;
		print "\nStep 5: $_ has been completed! It is the $num7/$fgn file.\n\n";
	}
	&msg("\nTotally, $num7 alignments have been created in '$TEMdir/TEM5_Alignment'!\n\n");

#============= Step 6: Trim Alignments ==============#
	&msg("\n\n\n#============= Step 6: Trim Alignments ==============#\n\n");

	&makeDir("$TEMdir/TEM6_AlnTrimmed");
	my $num8 = 0;
	my @emptyTrim;
	my $uuu=scalar @alnGenes;
	my $trimal = &tool("trimal");
	open (TRI, ">", "$TEMdir/trimAl_Details.txt") or die "Step 6: Can't open '$TEMdir/trimAl_Details.txt': $!\n";
	print TRI "Gene	Original_Length	Length_trimmed	Percentage\n";
	foreach my $filex (@alnGenes) {
		my $file = "$filex.fas.fasta";
		(my $outname = $file) =~ s/\.fas//;
		my $LenO = &alignmentLength("$TEMdir/TEM5_Alignment/$file");

		my $rc = &runCmd($trimal, "-in", "$TEMdir/TEM5_Alignment/$file", "-out", "$TEMdir/TEM6_AlnTrimmed/$outname", "-$trim");
		if ($rc != 0) {
			&fail("Step 6: trimAl failed to trim '$TEMdir/TEM5_Alignment/$file' (exit code $rc).\nPlease make sure that 'trimal$EXE' is present in '$BIN_DIR'.");
		}
		unless (-s "$TEMdir/TEM6_AlnTrimmed/$outname") {
			# trimAl writes no file if no column is left (e.g. '-trim nogaps' and a gap in every column).
			unlink "$TEMdir/TEM6_AlnTrimmed/$outname";
			print TRI "$filex	$LenO	0	0\n";
			push @emptyTrim, $filex;
			&msg("\nWARNING: Step 6: no column of '$file' is left after trimming with '-$trim'; gene '$filex' is skipped.\n");
			next;
		}
		my $LenA = &alignmentLength("$TEMdir/TEM6_AlnTrimmed/$outname");
		my $perc = $LenO ? sprintf("%.4f", $LenA/$LenO) : 0;
		print TRI "$filex	$LenO	$LenA	$perc\n";

		$num8++;
		my $rat =&givePercentage($num8, $uuu);
		print ("\r", "	Step 6: Trimming gene cluster ($rat): $filex.");
	}
	close TRI;
	&fail("Step 6: no alignment is left after trimming with '-trim $trim'. Please use another '-trim' value.") unless $num8;
	&msg("\nWARNING: ".scalar(@emptyTrim)." gene(s) without any column left after trimming were skipped: ".join(", ", @emptyTrim).".\n") if @emptyTrim;
	&msg("\nTotally, $num8 alignments have been trimmed successfully and written to '$TEMdir/TEM6_AlnTrimmed'!\n\n");
	&writeRecord("alignment_done" => &now(), "alignment_seq" => $seq_type, "trim" => $trim);
}


############### Task 4: Tree Inference #################
my ($fast,$iq,$astral);
my $iqModelOnly = 0;
my ($concatFile, @geneAlnFiles);
if ($task eq "all" || $task eq "tree_infer") {
	open (OPT, "<", $OPT_FILE) or &fail("Can't open '$OPT_FILE': $!");
	foreach my $i3 (<OPT>) {
		$i3=~ s/\r?\n$//;
		next unless $i3;
		if ( $i3 =~ s/^\s*IQ-TREE Command-line=// ) {
			$iq= $i3;
		} elsif ($i3 =~ s/^\s*FastTree Command-line=//) {
			$fast= $i3;
		} elsif ($i3 =~ s/^\s*astral Command-line=//) {
			$astral=$i3;
		}
	}
	close OPT;
	&fail("No command line for FastTree was found in '$OPT_FILE'.") if $tree_app eq "fasttree" && !$fast;
	&fail("No command line for IQ-TREE was found in '$OPT_FILE' (IQ-TREE is needed by '-tree_app iqtree' and '-tree cs').") if !$iq && ($tree_app eq "iqtree" || $tree_type eq "cs");
	&fail("No command line for ASTRAL was found in '$OPT_FILE'.") if $tree_type eq "st" && !$astral;
	# '-m MF' only selects the substitution model; IQ-TREE then infers no tree.
	$iqModelOnly = 1 if $tree_app eq "iqtree" && $iq =~ /(^|\s)-m\s+MF(\s|$)/;
	if ($iqModelOnly && $tree_type ne "sm") {
		&fail("The IQ-TREE command line in '$OPT_FILE' contains '-m MF' (model selection only, no tree), so no gene tree can be inferred for '-tree $tree_type'.\nPlease use '-m MFP' or remove the '-m' option.");
	}

	unless (-d "$TEMdir/TEM6_AlnTrimmed") {
		&fail("Step 7: Directory '$TEMdir/TEM6_AlnTrimmed' does not exist. Please run '-task alignment' first.");
	}
	if ($task eq "tree_infer" && %rec && !$rec{"alignment_done"}) {
		&msg("\nWARNING: according to '$RECORD', the alignments were not produced after the last '-task hmmsearch' or '-task refine'.\nThey may not match the current results. Run '-task alignment' first if needed.\n\n");
	}

	opendir(DIR, "$TEMdir/TEM6_AlnTrimmed") || die "Step 7: Can't open input directory '$TEMdir/TEM6_AlnTrimmed': $!\n";
	@geneAlnFiles = sort grep { /\.fasta$/ } readdir(DIR);
	closedir(DIR);
	&fail("Step 7: Directory '$TEMdir/TEM6_AlnTrimmed' does not contain any alignment files. Maybe you need to run the former task(s).") unless @geneAlnFiles;

	# Sequence type of the alignments (FastTree needs '-nt' for nucleotides).
	my $alnType = $rec{"alignment_seq"};
	$alnType = $seq_type if $task eq "all";
	$alnType = &checkSeqType("$TEMdir/TEM6_AlnTrimmed/$geneAlnFiles[0]") unless defined $alnType && ($alnType eq "nucl" || $alnType eq "prot");
	&msg("Sequence type of the alignments: ".($alnType eq "nucl" ? "nucleotide" : "protein")."\n");

	my ($fastProg, @fastArgs);
	if ($tree_app eq "fasttree") {
		($fastProg, @fastArgs) = split /\s+/, $fast;
		$fastProg = &tool($fastProg);
		if ($alnType eq "nucl" && !grep { $_ eq "-nt" } @fastArgs) {
			push @fastArgs, ("-nt", "-gtr");
			&msg("NOTE: '-nt -gtr' was added to the FastTree command line for the nucleotide alignments.\n");
		}
		$ENV{"OMP_NUM_THREADS"} = $thread;    # used by FastTreeMP
	}
	my ($iqProg, @iqArgs);
	if ($iq) {
		($iqProg, @iqArgs) = split /\s+/, $iq;
		$iqProg = &tool($iqProg) if $tree_app eq "iqtree" || $tree_type eq "cs";
		# Without '-redo', IQ-TREE refuses to overwrite the results of a previous run.
		push @iqArgs, "-redo" unless grep { $_ eq "-redo" } @iqArgs;
	}

	if ($tree_type eq "sm") {
#============= Step 7: Concatenate Sequences (If Necessary) ==============#

		&msg("\n\n############### Task 4: Tree Inference #################\n#============= Step 7: Concatenate Sequences ==============# \n\n");

		my $iii=scalar @geneAlnFiles;
		my $num10 = 0;
		my (%out,$totalLen);
		$totalLen = 0;
		foreach my $file (@geneAlnFiles) {
			my ($ids, $seqs) = &readFasta("$TEMdir/TEM6_AlnTrimmed/$file");
			my $len = 0;
			foreach my $id (@$ids) {
				$seqs->{$id} =~ s/\*/-/g;
				$len = length($seqs->{$id}) if length($seqs->{$id}) > $len;
			}
			# Taxa missing from this gene are filled with gaps, and new taxa get gaps for the previous genes.
			foreach my $id (@$ids) {
				$out{$id} = "-" x $totalLen unless exists $out{$id};
			}
			foreach my $id (keys %out) {
				my $s = exists $seqs->{$id} ? $seqs->{$id} : "";
				$s .= "-" x ($len - length($s)) if length($s) < $len;
				$out{$id} .= $s;
			}
			$totalLen += $len;
			$num10++;
			my $rat = &givePercentage($num10, $iii);
			print ("\r", "	Step 7: Reading the files: $num10/$iii ($rat)  .");
		}
		print ".\n\n";

		$concatFile = "$TEMdir/$inputName.$hmm.$num10.concatenation.fas";
		open (OUT, ">", $concatFile) or die "Can't open '$concatFile': $!\n";
		my $nn=0;
		foreach my $id (sort keys %out) {
			$nn++;
			print OUT ">$id\n";
			for (my $p = 0; $p < length($out{$id}); $p += 60) {
				print OUT substr($out{$id}, $p, 60), "\n";
			}
		}
		close OUT;
		&msg("Got $nn concatenated sequences ($totalLen sites) written to '$concatFile'!\n");

#============= Step 8: Phylogeny Inference ==============#
		&msg("\n\n\n#============= Step 8: Phylogeny Inference ==============#\n\n");
		print "Inferring phylogeny of the supermatrix approach......\n\n";

		my $outTree = "$outPrefix.$hmm.$num10.supermatrix.$tree_app.tree";
		&backupFile($outTree);
		if ($tree_app eq "fasttree") {
			unlink $outTree;
			my $rc = &runCapture("$concatFile.fasttree.log", $fastProg, @fastArgs, "-out", $outTree, $concatFile);
			&FastreeResponse("$concatFile.fasttree.log", $rc);
			&checkTree($outTree, "The supermatrix tree", "FastTree", "$concatFile.fasttree.log");
		} elsif ($tree_app eq "iqtree") {
			unlink "$concatFile.contree", "$concatFile.treefile";
			my $rc = &runCmd($iqProg, @iqArgs, "-s", $concatFile, "-pre", $concatFile, "-T", $thread);
			&IQtreeResponse("$concatFile.log", $rc);
			if ($iqModelOnly) {
				&msg("NOTE: IQ-TREE was run with '-m MF' (model selection only), so no tree was inferred.\n\n");
			} else {
				my $iqTree = (-s "$concatFile.contree") ? "$concatFile.contree" : "$concatFile.treefile";
				&checkTree($iqTree, "The supermatrix tree", "IQ-TREE", "$concatFile.log");
				&copyTree($iqTree, $outTree);
			}
		}
		unless ($iqModelOnly) {
			&checkTree($outTree, "The supermatrix tree", ($tree_app eq "iqtree" ? "IQ-TREE" : "FastTree"));
			&msg("The supermatrix tree has been written successfully in '$outTree'.\n\n");
		}

	} elsif ($tree_type eq "st" || $tree_type eq "cs") {
	#============= Step 7: Construct gene Trees ==============#
		&msg("\n\n############### Task 4: Tree Inference #################\n#============= Step 7: Construct gene Trees ==============#\n\n");

		if ($tree_type eq "cs" ) {
			my $seqnum=0;
			my $report="**Sequence numbers in some files:\n";
			foreach my $file11 (@geneAlnFiles) {
				my $nnum=&countSeq("$TEMdir/TEM6_AlnTrimmed/$file11");
				$report .="$TEMdir/TEM6_AlnTrimmed/$file11	$nnum\n";
				if ($seqnum==0) {
					$seqnum=$nnum;
				} elsif ($seqnum != $nnum) {
					&fail("Task 4 (-task tree_infer): Option '-tree cs' requires the alignment files in directory '$TEMdir/TEM6_AlnTrimmed/' to contain the same number of sequences.\n$report\nPlease run EasyCGTree from 'Task 2 (-task refine)' with '-tree cs', or run the whole pipeline.\n\n****Tips: if this ERROR still occurs while you have run 'Task 2' or the whole pipeline, the reason must be that some proteomes have no sequence retained during trimming of some gene clusters by trimAl. Please remove related cluster(s) from '$TEMdir/TEM6_AlnTrimmed/' manually and run 'Task 4 (-task tree_infer)' again.");
				}
			}
		}

		&removePath("$TEMdir/TEM7_GeneTrees") if -e "$TEMdir/TEM7_GeneTrees";
		&makeDir("$TEMdir/TEM7_GeneTrees");

		my $num10=scalar @geneAlnFiles;
		my $num11 = 0;
		my @geneTrees;
		foreach my $file11 (@geneAlnFiles) {
			$num11++;
			print "Step 7: Inferring phylogeny from '$TEMdir/TEM6_AlnTrimmed/$file11': $num11/$num10.\n";
			my $pre = "$TEMdir/TEM7_GeneTrees/$file11";
			if ($tree_app eq "fasttree") {
				my $rc = &runCapture("$pre.tree.log", $fastProg, @fastArgs, "-out", "$pre.tree", "$TEMdir/TEM6_AlnTrimmed/$file11");
				&FastreeResponse("$pre.tree.log", $rc);
				&checkTree("$pre.tree", "The gene tree", "FastTree", "$pre.tree.log");
				push @geneTrees, "$pre.tree";
			} elsif ($tree_app eq "iqtree") {
				# All IQ-TREE output files go directly to TEM7_GeneTrees.
				my $rc = &runCmd($iqProg, @iqArgs, "-s", "$TEMdir/TEM6_AlnTrimmed/$file11", "-pre", $pre, "-T", $thread);
				&IQtreeResponse("$pre.log", $rc);
				my $iqTree = (-s "$pre.contree") ? "$pre.contree" : "$pre.treefile";
				&checkTree($iqTree, "The gene tree", "IQ-TREE", "$pre.log");
				push @geneTrees, $iqTree;
			}
		}
		&msg("Totally, $num11 phylogenies have been written to '$TEMdir/TEM7_GeneTrees'.\n\n");

	#============= Step 8: Generate consensus Tree / supertree ==============#
		&msg("\n\n\n#============= Step 8: Generate ".($tree_type eq "cs" ? "consensus tree" : "supertree")." ==============#\n\n");

		# All gene trees in one file. Genome names are used directly (no renaming is needed
		# for IQ-TREE or ASTRAL, which both accept long names).
		my $intree = "$TEMdir/intree";
		open (OU11, ">", $intree) or die "Step 8: Can't open '$intree': $!\n";
		foreach my $gt (@geneTrees) {
			open (FIL11, "<", $gt) or die "Step 8: Can't open '$gt': $!\n";
			foreach my $in11 (<FIL11>) {
				$in11=~ s/\r?\n$//;
				next unless $in11=~ /\S/;
				print OU11 "$in11\n";
			}
			close FIL11;
		}
		close OU11;

		if ($tree_type eq "cs") {
			# Extended majority-rule consensus by IQ-TREE; support values are the percentage of gene trees.
			my $conPre = "$TEMdir/consensus";
			unlink "$conPre.contree";
			my $rc = &runCapture("$conPre.run.log", $iqProg, "-con", "-t", $intree, "-pre", $conPre, "-redo");
			&IQtreeResponse("$conPre.run.log", $rc);
			&checkTree("$conPre.contree", "The consensus tree", "IQ-TREE", "$conPre.run.log");
			my $outTree = "$outPrefix.$hmm.$num11.consensus.$tree_app.tree";
			&backupFile($outTree);
			&copyTree("$conPre.contree", $outTree);
			&msg("The consensus tree (support values: percentage of the $num11 gene trees) has been written successfully in '$outTree'.\n\n");
		} elsif ($tree_type eq "st" ) {
			my ($astralProg, @astralArgs) = split /\s+/, $astral;
			$astralProg = &tool($astralProg);
			my $outTree = "$outPrefix.$hmm.$num11.supertree.$tree_app.tree";
			&backupFile($outTree);
			unlink $outTree;
			my $rc = &runCmd($astralProg, @astralArgs, "-t", $thread, "-i", $intree, "-o", $outTree);
			if ($rc != 0) {
				&fail("Step 8: ASTRAL exited with an error (exit code $rc).".($rc >= 128 ? " The process was killed, probably because it ran out of memory." : ""));
			}
			&checkTree($outTree, "The supertree", "ASTRAL");
			&msg("The supertree has been written successfully in '$outTree'.\n\n");
		}
	}
}

############### Information from IQ-TREE #################
if ($tree_app eq "iqtree" && ($task eq "all" || $task eq "tree_infer")) {
	# The model is read from the IQ-TREE report (.iqtree), whose format is the same in IQ-TREE 2 and 3.
	# (IQ-TREE 3 no longer writes the table of tested models to the .log file.)
	if ($tree_type eq "sm") {
		&msg("\n\n*******************************\n\nUseful Information from IQ-TREE:\n\nModel	Model_selection	LogL	BIC\n\n");
		my ($model, $how, $lnl, $bic) = &iqtreeModel("$concatFile.iqtree");
		&msg("$model	$how	$lnl	$bic\n\n");
	} elsif ($tree_type eq "st" || $tree_type eq "cs") {
		&msg("\n\n*******************************\n\nUseful Information from IQ-TREE:\n\nNo.	Gene	Model	Model_selection	LogL	BIC\n\n");
		my $num22=0;
		foreach my $file11 (@geneAlnFiles) {
			$num22++;
			my $genet =$file11;
			$genet =~ s/\.fasta$//;
			my ($model, $how, $lnl, $bic) = &iqtreeModel("$TEMdir/TEM7_GeneTrees/$file11.iqtree");
			&msg("$num22	$genet	$model	$how	$lnl	$bic\n");
		}
		&msg("\n");
	}
	&msg("*******************************\n\n");
}

############### Ending Information ###############
&msg("\n\n\n############### Ending Information ###############\n\n");
my ($sec1,$min1,$hour1,$mday1,$mon1,$year1) = localtime(time);
&msg(sprintf("Job was finished at: %02d:%02d:%02d, %04d-%02d-%02d.\n", $hour1, $min1, $sec1, $year1+1900, $mon1+1, $mday1));
my $el = time - $START_TIME;
&msg(sprintf("Running time: %d d %d h %d min %d sec.\n\n\nEasyCGTree Version %s by Dao-Feng Zhang\n\n", int($el/86400), int(($el%86400)/3600), int(($el%3600)/60), $el%60, $VERSION));
close LOG;
exit 0;


############################## Subroutines ##############################

# Print a message on the screen and in the log file.
sub msg {
	my $m = shift;
	print $m;
	print LOG $m if defined fileno(LOG);
}

# Report an error (screen and log file) and stop.
sub fail {
	my $m = shift;
	print "\n\nERROR: $m\n\n";
	print LOG "\n\nERROR: $m\n\n" if defined fileno(LOG);
	close LOG if defined fileno(LOG);
	exit 1;
}

sub now {
	my @t = localtime(time);
	return sprintf("%04d-%02d-%02d %02d:%02d:%02d", $t[5]+1900, $t[4]+1, $t[3], $t[2], $t[1], $t[0]);
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

# Run a program without a shell (safe for paths with spaces). Returns the exit code;
# 128 + signal number if the program was killed (e.g. 137 = killed, often out of memory).
sub runCmd {
	my @cmd = @_;
	my $rc = system(@cmd);
	if ($rc == -1) {
		&fail("Can't execute '$cmd[0]': $!");
	}
	return ($rc & 127) ? 128 + ($rc & 127) : ($rc >> 8);
}

# Like runCmd, but the screen output (stdout and stderr) of the program is written to a file.
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

# Remove everything in the working directory except the CDS prediction.
sub cleanAfterPrediction {
	opendir(my $dh, $TEMdir) or die "Can't open '$TEMdir': $!\n";
	my @items = grep { $_ ne "." && $_ ne ".." } readdir($dh);
	closedir($dh);
	my %keep = map { $_ => 1 } ("TEM0_CDS", basename($TYPEFILE), basename($RECORD));
	foreach my $it (@items) {
		next if $keep{$it};
		&removePath("$TEMdir/$it");
	}
	&deleteRecord(qw(hmm hmm_file hmm_md5 evalue hmmsearch_done refine_done gene_cutoff genome_cutoff refine_evalue alignment_done alignment_seq trim));
}

########## Record file: one 'key=value' per line ##########
sub readRecord {
	my %r;
	return %r unless -f $RECORD;
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

sub saveRecord {
	my %r = @_;
	open(my $fh, ">", $RECORD) or die "Can't write '$RECORD': $!\n";
	print $fh "# EasyCGTree record of the tasks performed in this directory. Please do not edit.\n";
	foreach my $k (sort keys %r) {
		print $fh "$k=$r{$k}\n";
	}
	close $fh;
}

sub writeRecord {
	my %new = @_;
	my %r = &readRecord();
	@r{keys %new} = values %new;
	&saveRecord(%r);
	%rec = %r;
}

sub deleteRecord {
	my @keys = @_;
	my %r = &readRecord();
	delete @r{@keys};
	&saveRecord(%r);
	%rec = %r;
}

# Input sequence type of each genome, from the CDS prediction ('nucl' or 'prot').
sub readTypes {
	my %t;
	open(my $fh, "<", $TYPEFILE) or return %t;
	while (my $l = <$fh>) {
		$l =~ s/\r?\n$//;
		my @f = split /\t/, $l;
		next if !@f || $f[0] eq "Genome";
		$t{$f[0]} = $f[1];
	}
	close $fh;
	return %t;
}

sub md5File {
	my $f = shift;
	open(my $fh, "<", $f) or return "";
	binmode($fh);
	my $d = Digest::MD5->new->addfile($fh)->hexdigest;
	close $fh;
	return $d;
}

sub givePercentage {
	my $g1 =shift;
	my $g2 =shift;
	return "0.00%" unless $g2;
	return sprintf("%.2f%%", $g1*100/$g2);
}

# Files of the input directory (hidden files such as '.DS_Store' and sub-directories are skipped).
sub listInputFiles {
	my $d = shift;
	opendir(my $dh, $d) or die "Can't open input directory '$d': $!\n";
	my @f = sort grep { /\w/ && !/^\./ && -f "$d/$_" } readdir($dh);
	closedir($dh);
	return @f;
}

# Read a FASTA file. Returns (\@ids, \%seq); IDs are the first word of the header,
# sequences are in upper case without white space.
sub readFasta {
	my $f = shift;
	my (@ids, %seq, $id);
	open(my $fh, "<", $f) or die "Can't open '$f': $!\n";
	while (my $l = <$fh>) {
		$l =~ s/\r?\n$//;
		next unless length $l;
		if ($l =~ s/^\>//) {
			($id) = split /\s+/, $l;
			$id = "" unless defined $id;
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

# Copy a FASTA file with sequence lines in upper case and without '\r'.
sub copyFastaUpper {
	my ($in, $out) = @_;
	open(my $fi, "<", $in) or die "Can't open '$in': $!\n";
	open(my $fo, ">", $out) or die "Can't open '$out': $!\n";
	while (my $l = <$fi>) {
		$l =~ s/\r?\n$//;
		if ($l =~ /^\>/) {
			print $fo "$l\n";
		} else {
			$l =~ s/\s+//g;
			print $fo uc($l), "\n" if length $l;
		}
	}
	close $fi;
	close $fo;
}

sub alignmentLength {
	my $f = shift;
	my ($ids, $seqs) = &readFasta($f);
	return 0 unless @$ids;
	return length($seqs->{$ids->[0]});
}

# Stop with an explanation if a tree file is missing, empty or not a Newick tree.
sub checkTree {
	my ($file, $desc, $app, $appLog) = @_;
	my $why;
	if (!-e $file) {
		$why = "was not created";
	} elsif (-z $file) {
		$why = "is empty";
	} else {
		open(my $fh, "<", $file) or &fail("Can't open '$file': $!");
		local $/;
		my $c = <$fh>;
		close $fh;
		$why = "does not contain a complete tree" unless $c =~ /\(.*\)[^;]*;/s;
	}
	return unless $why;
	my $m = "$desc '$file' $why.\nThis occurs probably due to memory limitation (out of memory) or an interruption when running $app.\n";
	$m .= "Please check the messages of $app in '$appLog'.\n" if $appLog && -e $appLog;
	$m .= "Please try Task 4 (-task tree_infer) again after ending some processes of other applications currently running on your computer.\nOtherwise, we strongly recommend to try again on a more powerful machine";
	$m .= " or use FastTree instead with the '-tree_app fasttree' option" if $app eq "IQ-TREE";
	&fail("$m.");
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

# Check the messages and exit code of FastTree.
sub FastreeResponse {
	my ($logfile, $rc) = @_;
	my $response11 = "";
	if (open(my $fh, "<", $logfile)) {
		local $/;
		$response11 = <$fh>;
		close $fh;
	}
	$response11 = "" unless defined $response11;
	if ($response11 =~ /truncated/i && $response11 =~ /be too long/i) {
		&fail("Step 7: EasyCGTree stopped when executing 'FastTree':\n$response11\n\nPlease check whether the genome/proteome files were formatted correctly.\nPlease find more information in section 3.1 of the Manual.");
	} elsif ($response11 =~ /out of memory/i || $rc >= 128) {
		&fail("Step 7: EasyCGTree stopped when executing 'FastTree' (exit code $rc):\n$response11\n\nThis ERROR occurred probably because 'FastTree' needs larger memory space than that of your computer. Please find a more powerful PC or server to complete the run.");
	} elsif ($rc != 0) {
		&fail("Step 7: 'FastTree' exited with an error (exit code $rc):\n$response11");
	}
}

# Check the exit code of IQ-TREE and show the error messages from its log file.
sub IQtreeResponse {
	my ($logfile, $rc) = @_;
	return if $rc == 0;
	my @err;
	my $txt = "";
	if (open(my $fh, "<", $logfile)) {
		while (my $l = <$fh>) {
			$txt .= $l;
			push @err, $l if $l =~ /ERROR/;
		}
		close $fh;
	}
	my $m = "IQ-TREE exited with an error (exit code $rc).";
	$m .= " The process was killed, probably because it ran out of memory." if $rc >= 128;
	$m .= "\nMessages from '$logfile':\n".join("", @err) if @err;
	if ($txt =~ /Tree taxa and alignment sequence do not match/) {
		$m .= "\nThis is a known problem of IQ-TREE 2.2.0 when it runs with several threads. Please update IQ-TREE (e.g. IQ-TREE 3, https://iqtree.github.io), or use '-thread 1'.";
	}
	&fail($m);
}

# Model, how it was chosen, log-likelihood and BIC of the final tree from an IQ-TREE report (.iqtree).
sub iqtreeModel {
	my $f = shift;
	my ($model, $how, $lnl, $bic) = ("NA", "NA", "NA", "NA");
	open(my $fh, "<", $f) or return ($model, $how, $lnl, $bic);
	while (my $l = <$fh>) {
		$l =~ s/\r?\n$//;
		if ($l =~ /^Best-fit model according to (\w+):\s*(\S+)/) {
			($how, $model) = ("best-fit ($1)", $2);
		} elsif ($l =~ /^Model of substitution:\s*(\S+)/ && $model eq "NA") {
			($how, $model) = ("user-specified", $1);
		} elsif ($l =~ /^Log-likelihood of the tree:\s*(\S+)/) {
			$lnl = $1;
		} elsif ($l =~ /^Bayesian information criterion \(BIC\) score:\s*(\S+)/) {
			$bic = $1;
		}
	}
	close $fh;
	return ($model, $how, $lnl, $bic);
}

sub countSeq {
	my $dna = shift;
	open (my $fh, "<", $dna) || die "Can't open '$dna':$!.\n";
	my $count=0;
	while (my $in = <$fh>) {
		$count++ if $in =~ /^\>/;
	}
	close $fh;
	return $count;
}

########## Input files: classification, grouping and NCBI annotation files ##########

# Classify the input files and decide how each genome is prepared. Returns (\@units, \@notes);
# each unit is {name, type, source => [files], ignored => [files]}, with type:
#   nucl      - genome (DNA): CDS predicted by Prodigal
#   prot      - protein sequences, used directly
#   ncbi_pair - NCBI protein + CDS files with locus_tag (e.g. translated_cds.faa + cds_from_genomic.fna)
#   ncbi_cds  - NCBI CDS file with locus_tag (cds_from_genomic.fna), translated by EasyCGTree
# Files whose names contain the same assembly accession (GCF_/GCA_) belong to the same genome.
# Such genomes are named by the accession; ordinary files are renamed as before (sub renameFiles).
sub planInputs {
	my $dir = shift;
	my @files = &listInputFiles($dir);
	&fail("Input directory '$dir' does not contain any files.") unless @files;
	my (%kind, @bad, %groups, @singles, @notes);
	foreach my $f (@files) {
		my $st = &checkSeqType("$dir/$f");
		unless ($st eq "nucl" || $st eq "prot") {
			push @bad, $st;
			next;
		}
		my ($nh, $nlt) = &locusTagInfo("$dir/$f");
		my $lt = ($nh && $nlt >= 0.9 * $nh) ? 1 : 0;
		$kind{$f} = ($st eq "nucl") ? ($lt ? "cds_lt" : "genome") : ($lt ? "prot_lt" : "prot");
		my $acc = &accessionOf($f);
		if (defined $acc) {
			push @{ $groups{$acc} }, $f;
		} else {
			push @singles, $f;
		}
	}
	if (@bad) {
		&fail("Argument '-input $dir': the input directory should only contain files of FASTA-formatted sequences, while:\n".join("", @bad));
	}

	my %kindName = ("cds_lt" => "CDS DNA sequences with locus_tag", "prot_lt" => "protein sequences with locus_tag",
		"prot" => "protein sequences without locus_tag", "genome" => "genome (DNA) sequences");
	my @units;
	foreach my $acc (sort keys %groups) {
		my @g = @{ $groups{$acc} };
		my %by;
		push @{ $by{$kind{$_}} }, $_ foreach @g;
		foreach my $k (sort keys %by) {
			if (@{ $by{$k} } > 1) {
				&fail("The following files in '$dir' all contain $kindName{$k} of the assembly '$acc':\n	".join("\n	", @{ $by{$k} })."\nPlease keep only one of them.");
			}
		}
		# A single ordinary file is treated as before (named after the file).
		if (@g == 1 && ($kind{$g[0]} eq "genome" || $kind{$g[0]} eq "prot")) {
			push @singles, $g[0];
			next;
		}
		my ($cds) = @{ $by{"cds_lt"}  || [] };
		my ($plt) = @{ $by{"prot_lt"} || [] };
		my ($p)   = @{ $by{"prot"}    || [] };
		my ($gen) = @{ $by{"genome"}  || [] };
		my %u = ("name" => $acc, "ignored" => []);
		if ($cds) {
			if ($p) {
				&fail("For the assembly '$acc', the CDS file '$cds' (with locus_tag) is provided together with the protein file '$p', whose sequence headers contain no locus_tag (e.g. NCBI 'protein.faa').\nThe two files cannot be paired by locus_tag. Please delete '$p' from '$dir':\nthe CDS will then be translated by EasyCGTree (or provide NCBI 'translated_cds.faa' instead).");
			}
			if ($plt) {
				$u{"type"} = "ncbi_pair";
				$u{"source"} = [$plt, $cds];
			} else {
				$u{"type"} = "ncbi_cds";
				$u{"source"} = [$cds];
			}
			if ($gen) {
				push @{ $u{"ignored"} }, $gen;
				push @notes, "'$gen' (genome of '$acc') is ignored, because the annotated CDS file '$cds' of the same assembly is used.";
			}
		} elsif ($gen) {
			$u{"type"} = "nucl";
			$u{"source"} = [$gen];
			foreach my $x (grep { defined } ($plt, $p)) {
				push @{ $u{"ignored"} }, $x;
				push @notes, "'$x' (proteins of '$acc') is ignored, because the genome '$gen' of the same assembly is used for the CDS prediction.";
			}
		} else {
			$u{"type"} = "prot";
			$u{"source"} = [$plt];
			if ($p) {
				push @{ $u{"ignored"} }, $p;
				push @notes, "'$p' (proteins of '$acc') is ignored, because '$plt' of the same assembly is used.";
			}
		}
		push @units, \%u;
	}

	# Ordinary files: renamed as before; the genome name comes from the file name.
	my @renamed = &renameFiles($dir, @singles);
	foreach my $i (0..$#singles) {
		(my $name = $renamed[$i]) =~ s/\.fas$//;
		my $k = $kind{$singles[$i]};
		my $type = ($k eq "genome") ? "nucl" : ($k eq "cds_lt") ? "ncbi_cds" : "prot";
		push @units, { "name" => $name, "type" => $type, "source" => [$renamed[$i]], "orig" => [$singles[$i]], "ignored" => [] };
	}

	my %seen;
	foreach my $u (@units) {
		my $k = lc $u->{"name"};
		if ($seen{$k}) {
			&fail("Two genomes in '$dir' would get the same name '$u->{'name'}' ('".join("', '", @{ $seen{$k} })."' and '".join("', '", @{ $u->{'source'} })."'). Please rename or remove one of them.");
		}
		$seen{$k} = $u->{"source"};
	}
	@units = sort { $a->{"name"} cmp $b->{"name"} } @units;
	return (\@units, \@notes);
}

# Assembly accession in a file name, e.g. 'GCF_000009925.1_ASM992v1_protein.faa' -> 'GCF-000009925-1'.
sub accessionOf {
	my $f = shift;
	return "$1-$2-$3" if $f =~ /(GC[AF])[_-](\d{9})[._-](\d+)/;
	return undef;
}

# Number of sequences (among the first 200) and how many of their headers contain '[locus_tag=...]'.
# ('[old_locus_tag=...]' does not match, because '[' must come directly before 'locus_tag'.)
sub locusTagInfo {
	my $f = shift;
	my ($n, $lt) = (0, 0);
	open(my $fh, "<", $f) or die "Can't open '$f': $!\n";
	while (my $l = <$fh>) {
		next unless $l =~ /^\>/;
		$n++;
		$lt++ if $l =~ /\[locus_tag=[^\]]+\]/;
		last if $n >= 200;
	}
	close $fh;
	return ($n, $lt);
}

# Read an NCBI FASTA file: a list of {key, pseudo, partial5, seq}. The key is the locus_tag
# ('_2', '_3', ... is added if a locus_tag occurs more than once); undef without locus_tag.
sub readNcbi {
	my $f = shift;
	my (@r, %count, $cur);
	open(my $fh, "<", $f) or die "Can't open '$f': $!\n";
	while (my $l = <$fh>) {
		$l =~ s/\r?\n$//;
		if ($l =~ /^\>/) {
			my $key;
			if ($l =~ /\[locus_tag=([^\]]+)\]/) {
				$key = $1;
				$key =~ s/\s+/_/g;
				$count{$key}++;
				$key .= "_$count{$key}" if $count{$key} > 1;
			}
			$cur = { "key" => $key, "pseudo" => ($l =~ /\[pseudo=true\]/ ? 1 : 0),
				"partial5" => ($l =~ /\[partial=[^\]]*5'/ ? 1 : 0), "seq" => "" };
			push @r, $cur;
		} elsif ($cur) {
			$l =~ s/\s+//g;
			$cur->{"seq"} .= uc($l);
		}
	}
	close $fh;
	return \@r;
}

# NCBI protein + CDS files: pair the sequences by locus_tag and write them with the locus_tag as ID.
# At least 90% of the proteins must have a CDS with the same locus_tag.
sub ncbiPair {
	my ($protFile, $cdsFile, $faa, $fnn) = @_;
	my $P = &readNcbi($protFile);
	my $C = &readNcbi($cdsFile);
	my %cds = map { $_->{"key"} => $_ } grep { defined $_->{"key"} && !$_->{"pseudo"} } @$C;
	my $total = scalar @$P;
	my @pairs;
	foreach my $p (@$P) {
		next unless defined $p->{"key"} && $cds{$p->{"key"}} && length $p->{"seq"};
		push @pairs, [$p, $cds{$p->{"key"}}];
	}
	my $n = scalar @pairs;
	if ($total == 0 || $n < 0.9 * $total) {
		&fail("Step 0: Only $n of the $total proteins in '$protFile' have a CDS with the same locus_tag in '$cdsFile' (at least 90% are required).\nThe two files probably do not belong to the same genome. Please check them.");
	}
	open(my $fa, ">", $faa) or die "Can't open '$faa': $!\n";
	open(my $fn, ">", $fnn) or die "Can't open '$fnn': $!\n";
	foreach my $pc (@pairs) {
		(my $ps = $pc->[0]{"seq"}) =~ s/\*$//;
		print $fa ">$pc->[0]{'key'}\n$ps\n";
		print $fn ">$pc->[0]{'key'}\n$pc->[1]{'seq'}\n";
	}
	close $fa;
	close $fn;
	return ($n, $total);
}

# NCBI CDS file only: translate each CDS (genetic code 11) and write protein and DNA sequences
# with the locus_tag as ID. Pseudogenes, CDS whose length is not a multiple of 3 and CDS with
# internal stop codons are skipped.
sub ncbiTranslate {
	my ($cdsFile, $faa, $fnn) = @_;
	my $C = &readNcbi($cdsFile);
	my ($total, $n, $pseudo, $noTag, $badLen, $stop, $tgaOnly) = (0, 0, 0, 0, 0, 0, 0);
	open(my $fa, ">", $faa) or die "Can't open '$faa': $!\n";
	open(my $fn, ">", $fnn) or die "Can't open '$fnn': $!\n";
	foreach my $c (@$C) {
		$total++;
		if ($c->{"pseudo"}) { $pseudo++; next; }
		unless (defined $c->{"key"}) { $noTag++; next; }
		my $s = $c->{"seq"};
		if (length($s) == 0 || length($s) % 3) { $badLen++; next; }
		my ($prot, $nStop, $nonTGA) = &translateCDS($s, $c->{"partial5"});
		if ($nStop) {
			$stop++;
			$tgaOnly++ unless $nonTGA;
			next;
		}
		print $fa ">$c->{'key'}\n$prot\n";
		print $fn ">$c->{'key'}\n$s\n";
		$n++;
	}
	close $fa;
	close $fn;
	my $coding = $total - $pseudo;
	if ($coding && $tgaOnly > 0.1 * $coding) {
		&fail("Step 0: In '$cdsFile', $tgaOnly of $coding CDS contain internal TGA codons, which are stop codons in the genetic code 11.\nThis genome probably uses the genetic code 4 (TGA = Trp, e.g. Mycoplasma/Spiroplasma), which EasyCGTree does not support for NCBI CDS files.\nPlease provide the NCBI protein sequences with locus_tag ('translated_cds.faa') together with this CDS file instead.");
	}
	&fail("Step 0: No CDS could be translated from '$cdsFile'. Please check the file.") unless $n;
	return ($n, $total);
}

# Translate a CDS with the genetic code 11: the first codon is read as Met if it is a start codon
# (unless the CDS is 5'-partial, i.e. has no real start), a final stop codon is removed. Returns (protein, number of internal stops, any internal stop other than TGA).
my %CODON;
sub translateCDS {
	my ($s, $partial5) = @_;
	unless (%CODON) {
		my @b = split //, "TCAG";
		my @aa = split //, "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
		my $i = 0;
		foreach my $x (@b) { foreach my $y (@b) { foreach my $z (@b) { $CODON{"$x$y$z"} = $aa[$i++]; } } }
	}
	my %start = map { $_ => 1 } qw(ATG GTG TTG CTG ATT ATC ATA);
	my @cod = unpack("(A3)*", $s);
	my @p = map { defined $CODON{$_} ? $CODON{$_} : "X" } @cod;
	$p[0] = "M" if $start{$cod[0]} && !$partial5;
	if (@p && $p[-1] eq "*") {
		pop @p;
		pop @cod;
	}
	my ($nStop, $nonTGA) = (0, 0);
	foreach my $i (0..$#p) {
		next unless $p[$i] eq "*";
		$nStop++;
		$nonTGA = 1 if $cod[$i] ne "TGA";
	}
	return (join("", @p), $nStop, $nonTGA);
}

# Rename ordinary input files so that the names are safe for all the programs:
# white space, '_', '.' and brackets are replaced, and the extension becomes '.fas'.
# Returns the new names (in the same order).
sub renameFiles {
	my ($in, @files) = @_;
	my (%target, @new);
	foreach my $file (@files) {
		my $file1=$file;
		$file1 =~ s/\s/-/g;
		$file1 =~ s/_/-/g;
		$file1=~ s/\.\w+\z/\.fas/;
		$file1=~ s/\./-/g;
		$file1=~ s/-fas\z/\.fas/;
		$file1 .=".fas" unless $file1=~ /\./;
		$file1=~ s/[()]//g;
		if ($target{lc $file1}) {
			&fail("Both '$target{lc $file1}' and '$file' in '$in' would be renamed to '$file1'. Please give them clearly different names.");
		}
		if ($file1 ne $file && -e "$in/$file1" && !$DRY_RUN) {
			&fail("'$file' in '$in' should be renamed to '$file1', but a file with this name already exists. Please rename one of them.");
		}
		$target{lc $file1} = $file;
		push @new, $file1;
	}
	my $number1 = 0;
	foreach my $i (0..$#files) {
		next if $new[$i] eq $files[$i];
		$number1++;
		next if $DRY_RUN;
		&fail("Can't rename '$in/$files[$i]': $!") unless CORE::rename("$in/$files[$i]", "$in/$new[$i]");
	}
	&msg("$number1 input files ".($DRY_RUN ? "would be" : "have been")." renamed.\n") if $number1;
	return @new;
}

# 'nucl' or 'prot' for a FASTA file (upper and lower case are both accepted), otherwise an error message.
sub checkSeqType {
	my $dna = shift;
	open (my $fh, "<", $dna)|| die "Can't open '$dna':$!.\n";
	my $count=0;
	my $sig=0;
	my $seq="";
	while (my $in = <$fh>) {
		$in=~ s/\r?\n$//;
		if ($in =~ /^\>/) {
			$count++;
			$sig=1;
			last if $count > 5;
		} elsif ($sig) {
			$seq .= $in;
		}
		last if length($seq) > 1000000;
	}
	close $fh;
	return "$dna seems not like a FASTA-formatted file. Please check it and other files.\n" unless $sig;
	$seq = uc($seq);
	$seq =~ s/[\s\-\.\*]//g;
	my $len = length($seq);
	return "$dna does not contain any sequence. Please check it and other files.\n" unless $len;
	my $acgt = ($seq =~ tr/ACGTUN//);
	return ($acgt/$len > 0.9) ? "nucl" : "prot";
}

# Options of this script in JSON format, for EasyCGTree_GUI (JSON::PP is a core module since Perl 5.14).
sub printOptionsJson {
	require JSON::PP;
	my @o = (
		{ "name" => "input", "type" => "dir", "default" => undef, "required" => 1 },
		{ "name" => "task", "type" => "choice", "default" => "all", "values" => ["all", "predict", "hmmsearch", "refine", "alignment", "tree_infer"] },
		{ "name" => "outdir", "type" => "dir", "default" => undef },
		{ "name" => "seq", "type" => "choice", "default" => "prot", "values" => ["prot", "nucl"], "tasks" => ["all", "predict", "hmmsearch", "alignment"] },
		{ "name" => "keep_CDS_nucl", "type" => "bool", "default" => 0, "tasks" => ["all", "predict", "hmmsearch"] },
		{ "name" => "hmm", "type" => "text", "default" => "bac120", "tasks" => ["all", "hmmsearch", "refine", "alignment", "tree_infer"] },
		{ "name" => "evalue", "type" => "float", "default" => 1e-10, "tasks" => ["all", "hmmsearch", "refine"] },
		{ "name" => "genome_cutoff", "type" => "float", "default" => 0.8, "min" => 0, "max" => 1, "tasks" => ["all", "refine"] },
		{ "name" => "gene_cutoff", "type" => "float", "default" => 0.8, "min" => 0, "max" => 1, "tasks" => ["all", "refine"] },
		{ "name" => "trim", "type" => "choice", "default" => "strict", "values" => ["gappyout", "strict", "strictplus", "nogaps"], "tasks" => ["all", "alignment"] },
		{ "name" => "tree", "type" => "choice", "default" => "sm", "values" => ["sm", "st", "cs"], "tasks" => ["all", "refine", "tree_infer"] },
		{ "name" => "tree_app", "type" => "choice", "default" => "fasttree", "values" => ["fasttree", "iqtree"], "tasks" => ["all", "tree_infer"] },
		{ "name" => "thread", "type" => "int", "default" => 4, "min" => 1, "tasks" => ["all", "hmmsearch", "alignment", "tree_infer"] },
		{ "name" => "dry_run", "type" => "bool", "default" => 0 },
	);
	print JSON::PP->new->canonical(1)->encode({ "script" => "EasyCGTree.pl", "version" => $VERSION, "options" => \@o }), "\n";
}

sub usageError {
	print join("\n",@usage), "\nERROR: $_[0]\n\n";
	exit 1;
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

# Message for a missing profile HMM: the HMM sets are downloaded separately (folder 'HMM' of the repository).
sub hmmMissingMsg {
	my ($name, $dir) = @_;
	return "The profile HMM '$name' was not found in '$dir'.\n"
		. "The HMM sets are not included in the program package. Please download the ones you need (e.g. '$name.hmm')\n"
		. "from the folder 'HMM' of the EasyCGTree5 repository, and put them into '$dir':\n"
		. "	GitHub: https://github.com/zdf1987/EasyCGTree5 (folder 'HMM')\n"
		. "	Gitee:  https://gitee.com/zdf1987/EasyCGTree5 (folder 'HMM')\n"
		. "See also '$dir/README.txt'.";
}
