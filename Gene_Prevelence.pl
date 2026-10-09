#!/usr/bin/perl
#
# Gene_Prevelence - distribution of the genes of a profile HMM set in the genomes of an EasyCGTree run
# Part of EasyCGTree 5.0, by Dao-Feng Zhang
#
# The HMM search results of EasyCGTree ('<input>_TEM/TEM1_HMMsearch_out') are summarized: which genes
# were found in which genome, in how many copies, and how many genomes/genes would pass the cutoffs of
# EasyCGTree. If no matching HMM search results exist, 'EasyCGTree.pl -task hmmsearch' is run first.
# Only core Perl modules are used; the script can be called from any directory.
#
use warnings;
use strict;
use Getopt::Long;
use File::Path qw(mkpath);
use File::Basename qw(basename dirname);
use File::Spec;
use Cwd qw(abs_path);
use FindBin qw($RealBin);
use Digest::MD5;

$| = 1;

my $VERSION = "5.0";    # version of EasyCGTree
my $UPDATE  = "2026-10-09";
my $START_TIME = time;
my ($sec,$min,$hour,$mday,$mon,$year) = localtime($START_TIME);
my $STAMP     = sprintf("%04d-%02d-%02d-%02d-%02d", $year+1900, $mon+1, $mday, $hour, $min);
my $START_STR = sprintf("%02d:%02d:%02d, %04d-%02d-%02d", $hour, $min, $sec, $year+1900, $mon+1, $mday);
my $ECG     = "$RealBin/EasyCGTree.pl";
my $HMM_DIR = "$RealBin/HMM";

my @usage = qq(
====== Gene_Prevelence (EasyCGTree) ======
     EasyCGTree $VERSION, by Dao-Feng Zhang
     Update $UPDATE

Summarizes the presence and copy number of the genes of a profile HMM set in each genome, and shows how
many genomes and genes would pass the cutoffs of EasyCGTree ('-genome_cutoff', '-gene_cutoff').

Usage: perl Gene_Prevelence.pl -input <dir> [Options]

Essential Options:
-input <String>
	The input directory used (or to be used) with EasyCGTree.

Optional Options:
-outdir <String>
	Output directory for all results (working directory '<input>_TEM', trees, tables, log files).
	[default: the directory that contains the input directory]
-hmm <String>
	Profile HMM set ('bac120', 'rp1' ... or the path of a .hmm file). [default: the one of the existing
	HMM search results in '<input>_TEM'; bac120 if there are none]
-evalue <Real>
	E-value for counting hits. [default: the E-value of the existing HMM search, or 1e-10]
-force <no value required>
	If the existing HMM search results in '<input>_TEM' were obtained with another HMM (or a smaller
	E-value), run the HMM search again with EasyCGTree. This replaces the results of the later EasyCGTree
	tasks (refine, alignment, trees) in '<input>_TEM'; without '-force' the script stops instead.
-genome_cutoff <Decimal, 0..1>
	Genome cutoff for the preview of the EasyCGTree screening. [default: 0.8]
-gene_cutoff <Decimal, 0..1>
	Gene cutoff for the preview of the EasyCGTree screening. [default: 0.8]
-thread <Int>
	Number of threads (used only when the HMM search has to be run). [default: 4]
-help <no value required>
	Display this message.
-options_json <no value required>
	Print the options of this script in JSON format (used by EasyCGTree_GUI) and exit.

A hit is counted as in EasyCGTree: E-value <= '-evalue' and score >= 1/4 of the HMM length.

Output (in the directory containing the input directory; <out> = '<input>.<hmm>'):
<out>.gene_prevalence.txt   genome x gene: IDs of the proteins found ('/' between copies, '-' if absent;
                            '*': a protein used twice, see below)
<out>.gene_copies.txt       genome x gene: number of copies ('*': see below)
<out>.gene_summary.txt      each gene: genomes with the gene, prevalence, single-/multi-copy genomes
<out>.genome_summary.txt    each genome: genes found, completeness, single-/multi-copy and missing genes
<out>.shared_hits.txt       pairs of genes (HMMs) whose best hits are the same protein in some genomes

Shared best hits: an HMM of a protein domain, or of a paralogous gene, can hit the same protein as another
HMM of the set. EasyCGTree uses the best hit of each gene; if the best hits of two genes are one protein, it
is used for both genes, so its sequence is counted twice in the tree. Such proteins are marked with '*' in
the gene_prevalence and gene_copies tables, and '<out>.shared_hits.txt' lists the pairs of genes: in how
many genomes the protein is used twice, and in how many genomes the two genes are different proteins (e.g.
two domains fused in one gene in some genomes, but separate genes in others: pattern 'mixed'). This helps
to find HMMs of a set that overlap.

);

########## Options ##########
my %opt;
GetOptions(\%opt, "input:s", "outdir:s", "hmm:s", "evalue:f", "force!", "genome_cutoff:f", "gene_cutoff:f", "thread:i", "help!", "options_json!")
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
my $thread = 4;
if (exists $opt{"thread"}) {
	&usageError("'-thread' should be an integer >= 1.") unless defined $opt{"thread"} && $opt{"thread"} >= 1;
	$thread = $opt{"thread"};
}
my ($genomeCutoff, $geneCutoff) = (0.8, 0.8);
foreach my $k ("genome_cutoff", "gene_cutoff") {
	next unless exists $opt{$k};
	&usageError("'-$k' should be a decimal between 0 and 1.") unless defined $opt{$k} && $opt{$k} >= 0 && $opt{$k} <= 1;
}
$genomeCutoff = $opt{"genome_cutoff"} if exists $opt{"genome_cutoff"};
$geneCutoff   = $opt{"gene_cutoff"}   if exists $opt{"gene_cutoff"};
if (exists $opt{"evalue"}) {
	&usageError("'-evalue' should be a positive number.") unless defined $opt{"evalue"} && $opt{"evalue"} > 0;
}

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
&usageError("Neither the input directory '$inputDir' nor the working directory '$TEMdir' exists.") unless -d $inputDir || -d $TEMdir;

# Requested HMM: a name in the HMM directory, or a file.
my ($hmmWanted, $hmmWantedFile);
if (defined $opt{"hmm"} && $opt{"hmm"} ne "") {
	if (-f $opt{"hmm"}) {
		$hmmWantedFile = abs_path($opt{"hmm"});
		($hmmWanted = basename($opt{"hmm"})) =~ s/\.hmm$//i;
	} else {
		($hmmWanted = $opt{"hmm"}) =~ s/\.hmm$//i;
		$hmmWantedFile = "$HMM_DIR/$hmmWanted.hmm";
		unless (-f $hmmWantedFile) {
			print "\nERROR: ".&hmmMissingMsg($hmmWanted, $HMM_DIR)."\n\n";
			exit 1;
		}
	}
}

my %rec = &readRecord();
open(LOG, ">", "$outDir/$inputName.".($hmmWanted || $rec{"hmm"} || "hmm").".gene_prevalence_$STAMP.log") or die "Can't open the log file in '$outDir': $!\n";
$SIG{__DIE__} = sub { print LOG "\nERROR: $_[0]" if defined fileno(LOG); };

&msg("\n====== Gene_Prevelence (EasyCGTree) ======\n	EasyCGTree $VERSION\nby Dao-Feng Zhang\n\n");
&msg("Job Started at: $START_STR\n\nInput directory: $inputDir\nWorking directory: $TEMdir\n\n");


############### Step 1: HMM search results ###############
&msg("#============= Step 1: HMM search results ==============#\n\n");
my $haveSearch = (-d "$TEMdir/TEM1_HMMsearch_out" && -f "$TEMdir/HMMinfo.txt") ? 1 : 0;
my ($runSearch, $why) = (0, "");
if (!$haveSearch) {
	($runSearch, $why) = (1, "no HMM search results were found");
} elsif (!defined $rec{"hmm"}) {
	# Results of EasyCGTree < 5.0: the HMM cannot be checked.
	&msg("WARNING: no record of the HMM search was found in '$TEMdir' (made by EasyCGTree < 5.0?).\nThe existing results are used; the profile HMM cannot be checked.\n\n");
} else {
	if (defined $hmmWanted && $hmmWanted ne $rec{"hmm"}) {
		$why = "the existing HMM search results were obtained with '$rec{'hmm'}', not '$hmmWanted'";
	} elsif (defined $hmmWantedFile && $rec{"hmm_md5"} && &md5File($hmmWantedFile) ne $rec{"hmm_md5"}) {
		$why = "the file '$hmmWantedFile' has changed since the HMM search was performed";
	} elsif (exists $opt{"evalue"} && defined $rec{"evalue"} && $opt{"evalue"} > $rec{"evalue"}) {
		$why = "'-evalue $opt{'evalue'}' is larger than the E-value used by the existing HMM search ($rec{'evalue'}), so hits are missing";
	}
	if ($why) {
		unless ($opt{"force"}) {
			&fail("$why.\nRunning the HMM search again would replace the results of the later EasyCGTree tasks (refine, alignment, trees) in '$TEMdir'.\nUse '-force' to do so, or run without '-hmm'/'-evalue' to summarize the existing results, or copy the input directory to analyse it separately.");
		}
		$runSearch = 1;
	}
}
if ($runSearch) {
	&fail("'$ECG' was not found. Please put Gene_Prevelence.pl in the EasyCGTree directory.") unless -f $ECG;
	&fail("The input directory '$inputDir' does not exist, so the HMM search cannot be run.") unless -d $inputDir || -d "$TEMdir/TEM0_CDS";
	&msg("The HMM search will be run with EasyCGTree ($why).\n\n");
	my @cmd = ($^X, $ECG, "-input", $inputDir, "-outdir", $outDir, "-task", "hmmsearch", "-thread", $thread);
	push @cmd, ("-hmm", $hmmWantedFile) if defined $hmmWantedFile;
	push @cmd, ("-evalue", $opt{"evalue"}) if exists $opt{"evalue"};
	&msg("Running: perl EasyCGTree.pl ".join(" ", @cmd[2..$#cmd])."\n\n");
	my $rc = system(@cmd);
	$rc = ($rc == -1) ? 255 : ($rc & 127) ? 128 + ($rc & 127) : ($rc >> 8);
	&fail("EasyCGTree '-task hmmsearch' stopped with an error (exit code $rc). Please see the messages above.") if $rc;
	%rec = &readRecord();
	&msg("\n");
}
my $hmm    = $rec{"hmm"} || $hmmWanted || "unknown";
my $evalue = exists $opt{"evalue"} ? $opt{"evalue"} : ($rec{"evalue"} || 1e-10);
&msg("Profile HMM: $hmm".($rec{"hmm_file"} ? " ($rec{'hmm_file'})" : "")."\nE-value: $evalue\n\n");


############### Step 2: read the results ###############
&msg("#============= Step 2: Gene distribution ==============#\n\n");
# HMM information: NAME, ACC, length
my (@genes, %acc, %leng);
open(my $mi, "<", "$TEMdir/HMMinfo.txt") or &fail("Can't open '$TEMdir/HMMinfo.txt': $!");
while (my $l = <$mi>) {
	$l =~ s/\r?\n$//;
	next unless $l =~ /\S/;
	my @x = split /\t/, $l;
	next if $x[0] eq "NAME";
	push @genes, $x[0];
	$acc{$x[0]}  = (defined $x[1] && $x[1] ne "Null") ? $x[1] : "-";
	$leng{$x[0]} = $x[2] || 0;
}
close $mi;
&fail("No profile HMM is listed in '$TEMdir/HMMinfo.txt'.") unless @genes;

opendir(my $dh, "$TEMdir/TEM1_HMMsearch_out") or &fail("Can't open '$TEMdir/TEM1_HMMsearch_out': $!");
my @files = sort grep { /\.fas$/ } readdir($dh);
closedir($dh);
&fail("'$TEMdir/TEM1_HMMsearch_out' does not contain any HMM search results.") unless @files;

my (%hits, %keyCount, @genomes);
foreach my $f (@files) {
	(my $g = $f) =~ s/\.fas$//;
	push @genomes, $g;
	my %keys;
	open(my $fh, "<", "$TEMdir/TEM1_HMMsearch_out/$f") or &fail("Can't open '$TEMdir/TEM1_HMMsearch_out/$f': $!");
	while (my $l = <$fh>) {
		$l =~ s/\r?\n$//;
		next if $l =~ /^#/ || $l !~ /\S/;
		my @x = split /\s+/, $l;
		next if @x < 6 || $x[4] > $evalue;
		next if $x[5] < ($leng{$x[2]} || 0) / 4;
		push @{ $hits{$g}{$x[2]} }, $x[0] unless grep { $_ eq $x[0] } @{ $hits{$g}{$x[2]} || [] };
		# Gene key as in EasyCGTree (the accession if the HMM has one, otherwise the name).
		$keys{ ($x[3] =~ /-/) ? $x[2] : $x[3] } = 1;
	}
	close $fh;
	$keyCount{$g} = scalar keys %keys;
}
my $nG = scalar @genomes;
my $nGenes = scalar @genes;
&msg("$nG genomes and $nGenes genes ($hmm).\n\n");

########## Genes whose best hits are the same protein ##########
# EasyCGTree uses the best hit (first in the HMMER table) of each gene. If the best hits of two genes are one
# protein (e.g. two domains of one protein, or paralogs one of which is missing), this protein is used twice.
my (%pairTwice, %pairWhere, %pairDiff, %pattern, %usedTwice);
foreach my $g (@genomes) {
	my %byBest;    # protein -> genes whose best hit it is
	foreach my $gene (@genes) {
		push @{ $byBest{ $hits{$g}{$gene}[0] } }, $gene if $hits{$g}{$gene};
	}
	foreach my $p (sort keys %byBest) {
		my @gs = @{ $byBest{$p} };
		next if @gs < 2;
		$usedTwice{$g}{$_} = $p foreach @gs;
		for my $i (0..$#gs - 1) {
			for my $j ($i + 1..$#gs) {
				my $pair = join("\t", $gs[$i], $gs[$j]);
				$pairTwice{$pair}++;
				push @{ $pairWhere{$pair} }, "$g:$p";
			}
		}
	}
}
foreach my $pair (keys %pairTwice) {
	my ($g1, $g2) = split /\t/, $pair;
	# genomes with both genes whose best hits are different proteins
	$pairDiff{$pair} = grep { $hits{$_}{$g1} && $hits{$_}{$g2} && $hits{$_}{$g1}[0] ne $hits{$_}{$g2}[0] } @genomes;
	$pattern{$pair} = $pairDiff{$pair} ? "mixed" : "always one protein";
}
my @pairs = sort { $pairTwice{$b} <=> $pairTwice{$a} || $a cmp $b } keys %pairTwice;

########## Output tables ##########
my $prefix = "$outDir/$inputName.$hmm";
my (%present, %single, %multi, %maxc, %gPresent, %gSingle, %gMulti, @gMissing);
foreach my $g (@genomes) {
	foreach my $gene (@genes) {
		my $c = $hits{$g}{$gene} ? scalar @{ $hits{$g}{$gene} } : 0;
		next unless $c;
		$present{$gene}++;
		$gPresent{$g}++;
		if ($c == 1) { $single{$gene}++; $gSingle{$g}++; } else { $multi{$gene}++; $gMulti{$g}++; }
		$maxc{$gene} = $c if !$maxc{$gene} || $c > $maxc{$gene};
	}
}

my ($t1, $t2) = ("", "");
foreach my $t (\$t1, \$t2) {
	$$t .= join("\t", "Genome", "Genes_present", @genes)."\n";
	$$t .= join("\t", "ACC", "-", map { $acc{$_} } @genes)."\n";
}
foreach my $g (@genomes) {
	# '*': the best hit is also the best hit of another gene (used twice by EasyCGTree)
	$t1 .= join("\t", $g, $gPresent{$g} || 0, map { my $gn = $_; $hits{$g}{$gn}
		? join("/", map { ($usedTwice{$g} && ($usedTwice{$g}{$gn} // "") eq $_) ? "$_*" : $_ } @{ $hits{$g}{$gn} }) : "-" } @genes)."\n";
	$t2 .= join("\t", $g, $gPresent{$g} || 0, map { $hits{$g}{$_} ? scalar(@{ $hits{$g}{$_} }).(($usedTwice{$g} && $usedTwice{$g}{$_}) ? "*" : "") : 0 } @genes)."\n";
}
if (@pairs) {
	my $note = "# *: the best hit of this gene is also the best hit of another gene, so EasyCGTree uses the protein twice (see '".basename("$prefix.shared_hits.txt")."').\n";
	$t1 .= $note;
	$t2 .= $note;
}
&writeFile("$prefix.gene_prevalence.txt", $t1);
&writeFile("$prefix.gene_copies.txt", $t2);

my $s1 = "# Genes of '$hmm' in $nG genomes (E-value <= $evalue, score >= HMM length / 4).\n";
$s1 .= join("\t", qw(Gene ACC HMM_length Genomes_present Prevalence(%) Single_copy_genomes Multi_copy_genomes Max_copies Absent_in))."\n";
foreach my $gene (@genes) {
	my @abs = grep { !$hits{$_}{$gene} } @genomes;
	$s1 .= join("\t", $gene, $acc{$gene}, $leng{$gene}, $present{$gene} || 0, sprintf("%.1f", 100 * ($present{$gene} || 0) / $nG),
		$single{$gene} || 0, $multi{$gene} || 0, $maxc{$gene} || 0, (@abs ? join(",", @abs) : "-"))."\n";
}
&writeFile("$prefix.gene_summary.txt", $s1);

my $s2 = "# Genomes: genes of '$hmm' found (E-value <= $evalue, score >= HMM length / 4).\n";
$s2 .= join("\t", qw(Genome Genes_present Completeness(%) Single_copy_genes Multi_copy_genes Missing_genes Missing Multi_copy))."\n";
foreach my $g (@genomes) {
	my @miss = grep { !$hits{$g}{$_} } @genes;
	my @mc   = grep { $hits{$g}{$_} && @{ $hits{$g}{$_} } > 1 } @genes;
	$s2 .= join("\t", $g, $gPresent{$g} || 0, sprintf("%.1f", 100 * ($gPresent{$g} || 0) / $nGenes), $gSingle{$g} || 0, $gMulti{$g} || 0,
		scalar(@miss), (@miss ? join(",", @miss) : "-"), (@mc ? join(",", map { "$_(".scalar(@{ $hits{$g}{$_} }).")" } @mc) : "-"))."\n";
}
&writeFile("$prefix.genome_summary.txt", $s2);

my $s3 = "# Pairs of genes (HMMs) of '$hmm' whose best hits are the same protein in some genomes (E-value <= $evalue, score >= HMM length / 4), in $nG genomes.\n";
$s3 .= "# EasyCGTree uses the best hit of each gene, so such a protein is used for both genes (marked with '*' in the gene_prevalence and gene_copies tables).\n";
$s3 .= "# Genomes_used_twice: genomes in which the best hits of both genes are the same protein.\n";
$s3 .= "# Genomes_best_hits_differ: genomes in which both genes are found and their best hits are different proteins.\n";
$s3 .= "# Pattern: 'always one protein' (e.g. HMMs of two domains of one protein); 'mixed': one protein in some genomes, two in others\n";
$s3 .= "#          (e.g. domains fused in one gene in some genomes and separate genes in others, or paralogs one of which is missing).\n";
$s3 .= join("\t", qw(Gene_1 Gene_2 Genomes_used_twice Genomes_best_hits_differ Pattern Used_twice_in))."\n";
$s3 .= join("\t", $_, $pairTwice{$_}, $pairDiff{$_}, $pattern{$_}, join(",", @{ $pairWhere{$_} }))."\n" foreach @pairs;
$s3 .= "# None: no protein is the best hit of two genes.\n" unless @pairs;
&writeFile("$prefix.shared_hits.txt", $s3);


########## Overview on screen ##########
my @pr = map { 100 * ($present{$_} || 0) / $nG } @genes;
my @bins = ([100, 100, "in all genomes"], [90, 99.99, "in 90-<100% of the genomes"], [80, 89.99, "in 80-<90%"], [50, 79.99, "in 50-<80%"], [0, 49.99, "in <50%"]);
&msg("Gene prevalence:\n");
foreach my $bin (@bins) {
	my $n = grep { $_ >= $bin->[0] && $_ <= $bin->[1] } @pr;
	&msg(sprintf("	%4d genes %s\n", $n, $bin->[2]));
}
&msg(sprintf("	%4d genes are multi-copy in at least one genome\n", scalar grep { $multi{$_} } @genes));
my @comp = sort { $a <=> $b } map { 100 * ($gPresent{$_} || 0) / $nGenes } @genomes;
&msg(sprintf("Genome completeness (genes found / %d): min %.1f%%, median %.1f%%, max %.1f%%\n\n", $nGenes, $comp[0], $comp[int($#comp/2)], $comp[-1]));
if (@pairs) {
	&msg("Shared best hits: in ".scalar(@pairs)." pair(s) of genes, the best hits of both genes are the same protein in some genomes,\n");
	&msg("so EasyCGTree would use this protein twice (marked with '*' in the tables; genomes: used twice / best hits differ):\n");
	foreach my $pair (@pairs[0..($#pairs < 9 ? $#pairs : 9)]) {
		(my $p = $pair) =~ s/\t/ + /;
		&msg("	$p: $pairTwice{$pair} / $pairDiff{$pair}  ($pattern{$pair})\n");
	}
	&msg("	... (".(@pairs - 10)." more)\n") if @pairs > 10;
	&msg("See '$prefix.shared_hits.txt'.\n\n");
} else {
	&msg("Shared best hits: none (no protein is the best hit of two genes).\n\n");
}

########## Preview of the EasyCGTree screening (Task 2 'refine') ##########
my $maxKey = 0;
foreach my $g (@genomes) { $maxKey = $keyCount{$g} if $keyCount{$g} > $maxKey; }
my $cut2 = int($maxKey * $genomeCutoff);
my @sel = grep { $keyCount{$_} >= $cut2 } @genomes;
my @out = grep { $keyCount{$_} < $cut2 } @genomes;
my %keyPrev;
foreach my $g (@sel) {
	my %seen;
	open(my $fh, "<", "$TEMdir/TEM1_HMMsearch_out/$g.fas") or next;
	while (my $l = <$fh>) {
		next if $l =~ /^#/ || $l !~ /\S/;
		my @x = split /\s+/, $l;
		next if @x < 6 || $x[4] > $evalue || $x[5] < ($leng{$x[2]} || 0) / 4;
		my $k = ($x[3] =~ /-/) ? $x[2] : $x[3];
		$keyPrev{$k}++ unless $seen{$k}++;
	}
	close $fh;
}
my $cut = int(scalar(@sel) * $geneCutoff);
my $nKeys = scalar keys %keyPrev;
my $nPass = grep { $keyPrev{$_} >= $cut } keys %keyPrev;
&msg("Preview of the EasyCGTree screening (-genome_cutoff $genomeCutoff -gene_cutoff $geneCutoff):\n");
&msg("	Genomes with >= $cut2 genes ($maxKey * $genomeCutoff): ".scalar(@sel)." of $nG selected".(@out ? "; excluded: ".join(", ", map { "$_ ($keyCount{$_})" } @out) : "")."\n");
&msg("	Genes present in >= $cut (".scalar(@sel)." * $geneCutoff) selected genomes: $nPass (of the $nKeys genes found in the selected genomes)\n");
&msg("	(with '-tree cs', EasyCGTree uses '-gene_cutoff 1': ".scalar(grep { $keyPrev{$_} >= scalar(@sel) } keys %keyPrev)." genes)\n") if $geneCutoff < 1;
&msg("\n");

&msg("Output:\n	$prefix.gene_prevalence.txt\n	$prefix.gene_copies.txt\n	$prefix.gene_summary.txt\n	$prefix.genome_summary.txt\n	$prefix.shared_hits.txt\n\n");
my $el = time - $START_TIME;
&msg(sprintf("Running time: %d d %d h %d min %d sec.\n\nGene_Prevelence (EasyCGTree %s)\n\n", int($el/86400), int(($el%86400)/3600), int(($el%3600)/60), $el%60, $VERSION));
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

sub writeFile {
	my ($f, $c) = @_;
	open(my $fh, ">", $f) or &fail("Can't open '$f': $!");
	print $fh $c;
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

sub md5File {
	my $f = shift;
	open(my $fh, "<", $f) or return "";
	binmode($fh);
	my $d = Digest::MD5->new->addfile($fh)->hexdigest;
	close $fh;
	return $d;
}

# Options of this script in JSON format, for EasyCGTree_GUI (JSON::PP is a core module since Perl 5.14).
sub printOptionsJson {
	require JSON::PP;
	my @o = (
		{ "name" => "input", "type" => "dir", "default" => undef, "required" => 1 },
		{ "name" => "outdir", "type" => "dir", "default" => undef },
		{ "name" => "hmm", "type" => "text", "default" => undef },
		{ "name" => "evalue", "type" => "float", "default" => undef },
		{ "name" => "force", "type" => "bool", "default" => 0 },
		{ "name" => "genome_cutoff", "type" => "float", "default" => 0.8, "min" => 0, "max" => 1 },
		{ "name" => "gene_cutoff", "type" => "float", "default" => 0.8, "min" => 0, "max" => 1 },
		{ "name" => "thread", "type" => "int", "default" => 4, "min" => 1 },
	);
	print JSON::PP->new->canonical(1)->encode({ "script" => "Gene_Prevelence.pl", "version" => $VERSION, "options" => \@o }), "\n";
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
