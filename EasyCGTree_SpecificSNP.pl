#!/usr/bin/perl
#
# EasyCGTree_SpecificSNP - strain- and clade-specific SNPs on a tree, based on EasyCGTree_SNP
# Part of EasyCGTree 5.0, by Dao-Feng Zhang
#
# The SNP alignment and SNP positions made by EasyCGTree_SNP.pl are used (it is run first if they are
# missing). Each strain name in the tree gets '_<number of strain-specific SNPs>', and each clade gets
# '<number of clade-specific SNPs>/<support value>' as its label.
# Only core Perl modules are used. Put this script in the EasyCGTree directory.
#
use warnings;
use strict;
no warnings 'recursion';
use Getopt::Long;
use File::Path qw(mkpath);
use File::Basename qw(basename dirname);
use File::Spec;
use Cwd qw(abs_path);
use FindBin qw($RealBin);

$| = 1;

my $VERSION = "5.0";    # version of EasyCGTree
my $UPDATE  = "2026-09-28";
my $START_TIME = time;
my ($sec,$min,$hour,$mday,$mon,$year) = localtime($START_TIME);
my $STAMP     = sprintf("%04d-%02d-%02d-%02d-%02d", $year+1900, $mon+1, $mday, $hour, $min);
my $START_STR = sprintf("%02d:%02d:%02d, %04d-%02d-%02d", $hour, $min, $sec, $year+1900, $mon+1, $mday);
my $SNP_SCRIPT = "$RealBin/EasyCGTree_SNP.pl";

my @usage = qq(
====== EasyCGTree_SpecificSNP ======
     EasyCGTree $VERSION, by Dao-Feng Zhang
     Update $UPDATE

Strain- and clade-specific SNPs are counted on a tree and written into it:
  - each strain name gets '_N'          (N: SNPs specific to the strain), e.g. 'GCF-000009925-1_12'
  - each clade gets the label 'N/S'     (N: SNPs specific to the clade; S: support value of the branch)
The SNP alignment and positions of EasyCGTree_SNP.pl are used; if they are missing, EasyCGTree_SNP.pl
(and, through it, EasyCGTree) is run first.

Usage: perl EasyCGTree_SpecificSNP.pl -input <dir> [Options]

Essential Options:
-input <String>
	The input directory used with EasyCGTree / EasyCGTree_SNP. Results are written to the directory
	that contains the input directory.

Optional Options:
-outdir <String>
	Output directory for all results (working directory '<input>_TEM', trees, tables, log files).
	[default: the directory that contains the input directory]
-hmm <String>
	Profile HMM of the SNP results to use (needed only when SNP results of several HMMs exist).
-tree_app <String, 'iqtree' or 'fasttree'>
	SNP tree to annotate: '<input>.<hmm>.<genes>.snp.<tree_app>.tree'. [default: iqtree]
-tree <String>
	Annotate another tree file instead (e.g. the supermatrix tree of EasyCGTree). Its strains must be
	the same as those of the SNP alignment.
-outgroup <String>
	Root the tree with these strains (comma-separated; they must form one clade in the unrooted tree).
	Without '-outgroup', the tree is used as it is rooted in the file (IQ-TREE and FastTree trees are
	rooted arbitrarily).
-mode <String, 'exclusive' or 'strict'>
	Definition of a specific SNP. [default: exclusive]
	exclusive: all strains of the group (a strain or a clade) have the same base, and no other strain has
	           this base (other strains may have different bases).
	           Strains with a gap or an ambiguous base at a site are not counted as having any base there.
	strict:    as 'exclusive', and all other strains have one and the same base, i.e. the SNP splits all
	           strains exactly into the group and the rest. Only sites without gaps or ambiguous bases are
	           used. The count of each branch does not depend on the rooting of the tree.
-thread <Int>
	Number of threads (used only when EasyCGTree_SNP.pl has to be run). [default: 4]
-help <no value required>
	Display this message.
-options_json <no value required>
	Print the options of this script in JSON format (used by EasyCGTree_GUI) and exit.

Output (in the directory containing the input directory; <out> = '<input>.<hmm>.<genes>.specific_snp.<tree_app>',
or '<tree file name>.specific_snp' with '-tree'):
<out>.tree            tree with the numbers of specific SNPs
<out>.clade_ID.tree   the same tree with the clade IDs (C1, C2, ...) as labels, to find the clades in the tables
<out>.positions.txt   position of each specific SNP (group, gene, alignment positions, bases)
<out>.groups.txt      each strain and clade: support, number of specific SNPs, genes and positions, members

);

########## Options ##########
my %opt;
GetOptions(\%opt, "input:s", "outdir:s", "hmm:s", "tree_app:s", "tree:s", "outgroup:s", "mode:s", "thread:i", "help!", "options_json!")
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
my $tree_app = $opt{"tree_app"} || "iqtree";
&usageError("'-tree_app' should be 'iqtree' or 'fasttree'.") unless $tree_app eq "iqtree" || $tree_app eq "fasttree";
my $mode = $opt{"mode"} || "exclusive";
&usageError("'-mode' should be 'exclusive' or 'strict'.") unless $mode eq "exclusive" || $mode eq "strict";
my $thread = $opt{"thread"} || 4;
&usageError("'-thread' should be an integer >= 1.") unless $thread =~ /^\d+$/ && $thread >= 1;
if (defined $opt{"tree"} && $opt{"tree"} ne "") {
	&usageError("The tree file '$opt{'tree'}' does not exist.") unless -f $opt{"tree"};
}
my $hmmWanted;
if (defined $opt{"hmm"} && $opt{"hmm"} ne "") {
	($hmmWanted = basename($opt{"hmm"})) =~ s/\.hmm$//i;
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

open(LOG, ">", "$outDir/$inputName.".($hmmWanted || "snp").".specific_snp_$STAMP.log") or die "Can't open the log file in '$outDir': $!\n";
$SIG{__DIE__} = sub { print LOG "\nERROR: $_[0]" if defined fileno(LOG); };

&msg("\n====== EasyCGTree_SpecificSNP ======\n	EasyCGTree $VERSION\nby Dao-Feng Zhang\n\n");
&msg("Job Started at: $START_STR\n\nInput directory: $inputDir\n\n");


############### Step 1: SNP results ###############
&msg("#============= Step 1: SNP results ==============#\n\n");
my $set = &chooseSnpSet();
unless ($set) {
	&fail("No SNP results were found, and '$SNP_SCRIPT' is not available to create them.") unless -f $SNP_SCRIPT;
	&msg("No SNP results were found in '$outDir'; EasyCGTree_SNP.pl will be run first.\n\n");
	my @cmd = ($^X, $SNP_SCRIPT, "-input", $inputDir, "-outdir", $outDir, "-thread", $thread, "-tree_app", $tree_app);
	push @cmd, ("-hmm", $opt{"hmm"}) if defined $hmmWanted;
	push @cmd, "-no_tree" if defined $opt{"tree"} && $opt{"tree"} ne "";
	&msg("Running: perl EasyCGTree_SNP.pl ".join(" ", @cmd[2..$#cmd])."\n\n");
	my $rc = system(@cmd);
	$rc = ($rc == -1) ? 255 : ($rc & 127) ? 128 + ($rc & 127) : ($rc >> 8);
	&fail("EasyCGTree_SNP.pl stopped with an error (exit code $rc). Please see the messages above.") if $rc;
	$set = &chooseSnpSet();
	&fail("EasyCGTree_SNP.pl finished, but its results were not found in '$outDir'.") unless $set;
	&msg("\n");
}
my $prefix = $set->{"prefix"};
my $snpFas = "$prefix.snp.fas";
my $posFile = "$prefix.snp_positions.txt";
&msg("SNP alignment: $snpFas\nSNP positions: $posFile\n");

# Warn if the alignments of EasyCGTree are newer than the SNP results.
if (-d "$TEMdir/TEM6_AlnTrimmed" && opendir(my $dh, "$TEMdir/TEM6_AlnTrimmed")) {
	my $t0 = (stat($snpFas))[9];
	my @newer = grep { /\.fasta$/ && (stat("$TEMdir/TEM6_AlnTrimmed/$_"))[9] > $t0 } readdir($dh);
	closedir($dh);
	&msg("WARNING: the alignments in '$TEMdir/TEM6_AlnTrimmed' are newer than the SNP results.\nPlease run EasyCGTree_SNP.pl again if EasyCGTree has been run again since.\n") if @newer;
}

my ($taxa, $seq) = &readFasta($snpFas);
my @taxa = sort @$taxa;
my $nSnp = @taxa ? length($seq->{$taxa[0]}) : 0;
&fail("The SNP alignment '$snpFas' is empty.") unless $nSnp;
my ($posCols, $posRows) = &readTable($posFile);
&fail("'$posFile' lists ".scalar(@$posRows)." SNPs, but '$snpFas' has $nSnp columns. Please run EasyCGTree_SNP.pl again.") unless @$posRows == $nSnp;
my @refCols = grep { my $c = $_; grep { $_ eq $c } @$posCols } qw(Ref_CDS_ID Ref_CDS_pos Codon_pos Ref_base);
&msg(scalar(@taxa)." strains and $nSnp SNP sites were read.\n\n");


############### Step 2: tree ###############
&msg("#============= Step 2: Tree ==============#\n\n");
my ($treeFile, $outBase);
if (defined $opt{"tree"} && $opt{"tree"} ne "") {
	$treeFile = abs_path($opt{"tree"});
	(my $b = basename($treeFile)) =~ s/\.(tree|treefile|contree|nwk|newick|tre)$//i;
	$outBase = "$outDir/$b.specific_snp";
} else {
	$treeFile = "$prefix.snp.$tree_app.tree";
	unless (-f $treeFile) {
		&fail("The SNP tree '$treeFile' was not found.\nPlease run: perl EasyCGTree_SNP.pl -input $inputArg -tree_app $tree_app (with the options used before),\nor give a tree with '-tree'.");
	}
	$outBase = "$prefix.specific_snp.$tree_app";
}
&msg("Tree: $treeFile\n");
my $root = &readTree($treeFile);

# The strains of the tree and of the SNP alignment must be the same.
my @leaves = map { $_->{"name"} } &leafNodes($root);
{
	my %inTree = map { $_ => 1 } @leaves;
	my %inAln  = map { $_ => 1 } @taxa;
	my @onlyTree = sort grep { !$inAln{$_} } keys %inTree;
	my @onlyAln  = sort grep { !$inTree{$_} } keys %inAln;
	my %cnt;
	my @dup = grep { $cnt{$_}++ == 1 } @leaves;
	&fail("The tree contains the strain name(s) ".join(", ", @dup)." more than once.") if @dup;
	if (@onlyTree || @onlyAln) {
		my $m = "The strains of the tree and of the SNP alignment differ.\n";
		$m .= "Only in the tree:\n	".join("\n	", @onlyTree)."\n" if @onlyTree;
		$m .= "Only in the SNP alignment:\n	".join("\n	", @onlyAln)."\n" if @onlyAln;
		&fail($m);
	}
}

# Support values: 'SH-aLRT/UFBoot' labels (IQ-TREE .treefile) are reduced to the last value.
my ($nSlash, $nSup) = (0, 0);
foreach my $n (&allNodes($root)) {
	next if $n->{"leaf"} || !defined $n->{"support"} || $n->{"support"} eq "";
	$nSup++;
	if ($n->{"support"} =~ /\//) {
		$nSlash++;
		$n->{"support"} = (split /\//, $n->{"support"})[-1];
	}
}
&msg("NOTE: the tree has support labels like 'SH-aLRT/UFBoot'; only the last value (UFBoot) is used.\n") if $nSlash;
&msg("NOTE: the tree has no support values; clades are labelled with the number of specific SNPs only.\n") unless $nSup;

if (defined $opt{"outgroup"} && $opt{"outgroup"} ne "") {
	my (@og, @unknown, @screened);
	my @all = &allGenomes();
	foreach my $x (grep { /\S/ } split /\s*,\s*/, $opt{"outgroup"}) {
		my $t = &matchTaxon($x, \@taxa);
		if (defined $t) { push @og, $t; }
		elsif (defined &matchTaxon($x, \@all)) { push @screened, &matchTaxon($x, \@all); }
		else { push @unknown, $x; }
	}
	&fail("The following outgroup strains were not found in the tree:\n	".join("\n	", @unknown)) if @unknown;
	if (@screened) {
		&fail("No outgroup strain is left: ".join(", ", @screened)." ".(@screened > 1 ? "were" : "was")." removed when the genomes were screened ('-genome_cutoff') and ".(@screened > 1 ? "are" : "is")." not in the tree.\nPlease choose other outgroup strains, or run with a lower '-genome_cutoff'.") unless @og;
		&msg("WARNING: the outgroup strain(s) ".join(", ", @screened)." are not in the tree, because they were removed when the genomes were screened ('-genome_cutoff'); the tree is rooted with the others.\n");
	}
	&fail("The outgroup cannot contain all the strains.") if @og >= @taxa;
	$root = &rerootOutgroup($root, \@og);
	&msg("The tree has been rooted with the outgroup: ".join(", ", @og).".\n");
} else {
	&msg("NOTE: no outgroup was given ('-outgroup'), so the tree is used as it is rooted in the file.\n") if $mode eq "exclusive";
}
&msg("Definition of specific SNPs: '$mode'.\n\n");

# Leaf sets, clade IDs and the lookup table from a set of strains to its node.
my %byKey;
my $cladeNo = 0;
foreach my $n (&allNodes($root)) {
	next if $n == $root;
	my @l = sort map { $_->{"name"} } &leafNodes($n);
	$n->{"members"} = \@l;
	$n->{"id"} = $n->{"leaf"} ? $n->{"name"} : "C".(++$cladeNo);
	$n->{"snps"} = [];
	$byKey{join("\t", @l)} = $n;
}


############### Step 3: specific SNPs ###############
&msg("#============= Step 3: Strain- and clade-specific SNPs ==============#\n\n");
my ($nAssigned, $nStrainSnp, $nCladeSnp) = (0, 0, 0);
for (my $i = 0; $i < $nSnp; $i++) {
	my (%al, $miss);
	$miss = 0;
	foreach my $t (@taxa) {
		my $b = substr($seq->{$t}, $i, 1);
		if ($b =~ /[ACGT]/) { push @{ $al{$b} }, $t; } else { $miss++; }
	}
	my $nAl = keys %al;
	next if $nAl < 2 || ($mode eq "strict" && ($nAl != 2 || $miss));
	my $any = 0;
	foreach my $b (sort keys %al) {
		my $node = $byKey{join("\t", @{ $al{$b} })};
		next unless $node;
		my $others = join(",", map { "$_:".scalar(@{ $al{$_} }) } sort { @{ $al{$b} } <=> @{ $al{$a} } || $a cmp $b } grep { $_ ne $b } keys %al);
		push @{ $node->{"snps"} }, [$i, $b, $others, $miss];
		if ($node->{"leaf"}) { $nStrainSnp++; } else { $nCladeSnp++; }
		$any = 1;
	}
	$nAssigned++ if $any;
}
my @groups = sort { ($b->{"leaf"} <=> $a->{"leaf"}) || ($a->{"leaf"} ? $a->{"id"} cmp $b->{"id"} : substr($a->{"id"},1) <=> substr($b->{"id"},1)) } grep { $_ != $root } &allNodes($root);

&msg("Strain-specific SNPs (strain	SNPs):\n");
foreach my $g (grep { $_->{"leaf"} } @groups) {
	&msg("$g->{'id'}	".scalar(@{ $g->{'snps'} })."\n");
}
&msg("\nClade-specific SNPs (clade	SNPs/support	strains):\n");
foreach my $g (grep { !$_->{"leaf"} } @groups) {
	&msg("$g->{'id'}	".scalar(@{ $g->{'snps'} })."/".(defined $g->{'support'} && $g->{'support'} ne "" ? $g->{'support'} : "-")."	".scalar(@{ $g->{'members'} })."\n");
}
&msg("\nIn total: $nStrainSnp strain-specific and $nCladeSnp clade-specific SNPs.\n");
&msg("$nAssigned of the $nSnp SNP sites are specific to at least one strain or clade; the other ".($nSnp - $nAssigned)." sites\nare shared by strains that do not form a clade in this tree (homoplasy or conflict with the tree)".($mode eq "strict" ? ", have more than 2 bases or missing data," : "")." or have no group-specific base.\n");
&msg("NOTE: the two branches at the root are one branch of the unrooted tree; the SNPs separating the two sides are given on both.\n") if @{ $root->{"children"} } == 2;
&msg("\n");


############### Output ###############
my $outTree = "$outBase.tree";
&writeFile($outTree, &toNewick($root, sub {
	my $n = shift;
	return $n->{"name"}."_".scalar(@{ $n->{"snps"} }) if $n->{"leaf"};
	return "" if $n == $root;
	my $s = scalar(@{ $n->{"snps"} });
	return (defined $n->{"support"} && $n->{"support"} ne "") ? "$s/$n->{'support'}" : $s;
}).";\n");
&writeFile("$outBase.clade_ID.tree", &toNewick($root, sub {
	my $n = shift;
	return $n->{"name"} if $n->{"leaf"};
	return ($n == $root) ? "" : $n->{"id"};
}).";\n");

my @posHead = qw(Group_ID Group_type SNP_No Gene Gene_No Concat_pos Trimmed_aln_pos Original_aln_pos Group_base Other_bases Missing);
my $pos = "# Strain- and clade-specific SNPs ('-mode $mode'), in the order of the groups and of the SNP alignment.\n";
$pos .= "# Group_ID: strain name, or clade ID (see '".basename("$outBase.clade_ID.tree")."' and '".basename("$outBase.groups.txt")."').\n";
$pos .= "# SNP_No: column in '".basename($snpFas)."'; the other positions as in '".basename($posFile)."'.\n";
$pos .= "# Group_base: base of the group; Other_bases: bases of the other strains (number of strains); Missing: strains with a gap or ambiguous base.\n";
$pos .= join("\t", @posHead, @refCols)."\n";
foreach my $g (@groups) {
	foreach my $s (@{ $g->{"snps"} }) {
		my ($i, $b, $others, $miss) = @$s;
		my $r = $posRows->[$i];
		$pos .= join("\t", $g->{"id"}, ($g->{"leaf"} ? "strain" : "clade"), $i + 1,
			(map { defined $r->{$_} ? $r->{$_} : "NA" } qw(Gene Gene_No Concat_pos Trimmed_aln_pos Original_aln_pos)),
			$b, ($others eq "" ? "-" : $others), $miss, (map { defined $r->{$_} ? $r->{$_} : "NA" } @refCols))."\n";
	}
}
&writeFile("$outBase.positions.txt", $pos);

my $grp = "# Strains and clades of '".basename($outTree)."' with their specific SNPs ('-mode $mode').\n";
$grp .= "# Genes: gene(number of specific SNPs), in the order of the concatenation; Original_aln_positions: positions of the specific SNPs\n";
$grp .= "# in the original alignment of each gene (before trimming), in the same order ('gene:pos,pos;gene:pos').\n";
$grp .= join("\t", qw(Group_ID Group_type Support Strains Specific_SNPs SNP_No Genes Original_aln_positions Members))."\n";
foreach my $g (@groups) {
	my (@genes, %byGene, @snpNo);
	foreach my $s (@{ $g->{"snps"} }) {
		my $r = $posRows->[$s->[0]];
		my $gene = defined $r->{"Gene"} ? $r->{"Gene"} : "NA";
		push @genes, $gene unless exists $byGene{$gene};
		push @{ $byGene{$gene} }, (defined $r->{"Original_aln_pos"} ? $r->{"Original_aln_pos"} : "NA");
		push @snpNo, $s->[0] + 1;
	}
	$grp .= join("\t", $g->{"id"}, ($g->{"leaf"} ? "strain" : "clade"),
		((!$g->{"leaf"} && defined $g->{"support"} && $g->{"support"} ne "") ? $g->{"support"} : "-"),
		scalar(@{ $g->{"members"} }), scalar(@{ $g->{"snps"} }),
		(@snpNo ? &ranges(@snpNo) : "-"),
		(@genes ? join(",", map { "$_(".scalar(@{ $byGene{$_} }).")" } @genes) : "-"),
		(@genes ? join(";", map { "$_:".join(",", @{ $byGene{$_} }) } @genes) : "-"),
		join(",", @{ $g->{"members"} }))."\n";
}
&writeFile("$outBase.groups.txt", $grp);

&msg("Annotated tree:        $outTree\nTree with clade IDs:   $outBase.clade_ID.tree\nSpecific SNP positions: $outBase.positions.txt\nStrains and clades:    $outBase.groups.txt\n\n");
my $el = time - $START_TIME;
&msg(sprintf("Running time: %d d %d h %d min %d sec.\n\nEasyCGTree_SpecificSNP (EasyCGTree %s)\n\n", int($el/86400), int(($el%86400)/3600), int(($el%3600)/60), $el%60, $VERSION));
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

sub writeFile {
	my ($f, $c) = @_;
	&backupFile($f);
	open(my $fh, ">", $f) or &fail("Can't open '$f': $!");
	print $fh $c;
	close $fh;
}

# Compact list of numbers: 1,2,3,7,9,10 -> 1-3,7,9-10
sub ranges {
	my @n = @_;
	my @r;
	my ($s, $e) = ($n[0], $n[0]);
	foreach my $x (@n[1..$#n], undef) {
		if (defined $x && $x == $e + 1) { $e = $x; next; }
		push @r, ($s == $e ? $s : "$s-$e");
		($s, $e) = ($x, $x) if defined $x;
	}
	return join(",", @r);
}

sub readRecord {
	my %r;
	open(my $fh, "<", "$TEMdir/EasyCGTree_record.txt") or return %r;
	while (my $l = <$fh>) {
		$l =~ s/\r?\n$//;
		next if $l =~ /^\s*#/ || $l !~ /=/;
		my ($k, $v) = split /=/, $l, 2;
		$r{$k} = $v;
	}
	close $fh;
	return %r;
}

# SNP results '<input>.<hmm>.<genes>.snp.fas' (+ '.snp_positions.txt') in the output directory.
sub chooseSnpSet {
	opendir(my $dh, $outDir) or &fail("Can't open '$outDir': $!");
	my @sets;
	foreach my $f (sort readdir($dh)) {
		next unless index($f, "$inputName.") == 0 && $f =~ /\.snp\.fas$/;
		my $mid = substr($f, length($inputName) + 1);
		$mid =~ s/\.snp\.fas$//;
		next unless $mid =~ /^(.+)\.(\d+)$/;
		my $p = "$outDir/$inputName.$mid";
		push @sets, { "hmm" => $1, "genes" => $2, "prefix" => $p } if -f "$p.snp_positions.txt";
	}
	closedir($dh);
	@sets = grep { $_->{"hmm"} eq $hmmWanted } @sets if defined $hmmWanted;
	return undef unless @sets;
	return $sets[0] if @sets == 1;
	# Several SNP results: take the one that matches the current EasyCGTree results.
	my %r = &readRecord();
	my $ng = 0;
	if (opendir(my $d6, "$TEMdir/TEM6_AlnTrimmed")) {
		$ng = grep { /\.fasta$/ } readdir($d6);
		closedir($d6);
	}
	my @m = grep { defined $r{"hmm"} && $_->{"hmm"} eq $r{"hmm"} && $_->{"genes"} == $ng } @sets;
	return $m[0] if @m == 1;
	&fail("Several SNP results were found in '$outDir':\n	".join("\n	", map { basename($_->{'prefix'}).".snp.fas" } @sets)."\nPlease choose one with '-hmm'.");
}

sub readFasta {
	my $f = shift;
	my (@ids, %seq, $id);
	open(my $fh, "<", $f) or &fail("Can't open '$f': $!");
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

# Tab-separated table with '#' comment lines and a header line. Returns (\@columns, \@rows of hashes).
sub readTable {
	my $f = shift;
	open(my $fh, "<", $f) or &fail("Can't open '$f': $!");
	my (@cols, @rows);
	while (my $l = <$fh>) {
		$l =~ s/\r?\n$//;
		next if $l =~ /^#/ || $l !~ /\S/;
		my @x = split /\t/, $l;
		if (!@cols) { @cols = @x; next; }
		my %h;
		@h{@cols} = @x;
		push @rows, \%h;
	}
	close $fh;
	return (\@cols, \@rows);
}

# All genomes of the input directory (from the CDS prediction), including those removed by the screening.
sub allGenomes {
	my @g;
	if (open(my $th, "<", "$TEMdir/Genome_SeqType.txt")) {
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
	if ($n =~ /(GC[AF])[_-](\d{9})[._-](\d+)/ && exists $lc{lc "$1-$2-$3"}) {
		return $lc{lc "$1-$2-$3"};
	}
	my $m = $n;
	$m =~ s/\.\w+$//;
	$m =~ s/[\s_.()]+/-/g;
	return $lc{lc $m} if exists $lc{lc $m};
	return undef;
}

########## Newick trees ##########
# Node: {leaf, name (leaves), support (label of internal nodes), len, children}

sub readTree {
	my $f = shift;
	open(my $fh, "<", $f) or &fail("Can't open '$f': $!");
	local $/;
	my $s = <$fh>;
	close $fh;
	$s = "" unless defined $s;
	$s =~ s/\[[^\]]*\]//g;    # comments
	$s =~ s/\s+//g;
	$s =~ s/;.*$//s;
	&fail("'$f' does not contain a tree.") unless $s =~ /^\(/;
	my $p = 0;
	my $root = eval { &parseNode(\$s, \$p) };
	&fail("Can't read the tree in '$f': $@") if $@;
	&fail("Can't read the tree in '$f': unexpected text at position $p.") if $p != length($s);
	return $root;
}

sub parseNode {
	my ($sr, $pr) = @_;
	my $n = { "children" => [] };
	if (substr($$sr, $$pr, 1) eq "(") {
		$$pr++;
		while (1) {
			push @{ $n->{"children"} }, &parseNode($sr, $pr);
			my $c = substr($$sr, $$pr, 1);
			if ($c eq ",") { $$pr++; next; }
			if ($c eq ")") { $$pr++; last; }
			die "unexpected '".($c eq "" ? "end of tree" : $c)."' at position $$pr\n";
		}
	}
	my $lab = "";
	if (substr($$sr, $$pr, 1) eq "'") {
		$$pr++;
		while ($$pr < length($$sr) && substr($$sr, $$pr, 1) ne "'") { $lab .= substr($$sr, $$pr++, 1); }
		$$pr++;
	} else {
		while ($$pr < length($$sr) && substr($$sr, $$pr, 1) !~ /[:,();]/) { $lab .= substr($$sr, $$pr++, 1); }
	}
	if (@{ $n->{"children"} }) {
		$n->{"leaf"} = 0;
		$n->{"support"} = $lab;
	} else {
		die "a strain without name at position $$pr\n" if $lab eq "";
		$n->{"leaf"} = 1;
		$n->{"name"} = $lab;
	}
	if (substr($$sr, $$pr, 1) eq ":") {
		$$pr++;
		my $l = "";
		while ($$pr < length($$sr) && substr($$sr, $$pr, 1) !~ /[,();]/) { $l .= substr($$sr, $$pr++, 1); }
		$n->{"len"} = $l;
	}
	return $n;
}

sub toNewick {
	my ($n, $label) = @_;
	my $s = $n->{"leaf"} ? $label->($n) : "(".join(",", map { &toNewick($_, $label) } @{ $n->{"children"} }).")".$label->($n);
	$s .= ":$n->{'len'}" if defined $n->{"len"} && $n->{"len"} ne "";
	return $s;
}

sub allNodes {
	my $n = shift;
	return ($n, map { &allNodes($_) } @{ $n->{"children"} });
}

sub leafNodes {
	my $n = shift;
	return $n if $n->{"leaf"};
	return map { &leafNodes($_) } @{ $n->{"children"} };
}

# Root the tree on the branch that separates the outgroup from the other strains.
# Support values belong to branches, so they stay with their branches.
sub rerootOutgroup {
	my ($root, $og) = @_;
	my %og = map { $_ => 1 } @$og;
	# Undirected graph: node index -> [ [neighbour, edge], ... ]
	my @nodes = &allNodes($root);
	my %idx;
	@idx{map { "$_" } @nodes} = (0..$#nodes);
	my @adj = map { [] } @nodes;
	my $link = sub {
		my ($i, $j, $e) = @_;
		push @{ $adj[$i] }, [$j, $e];
		push @{ $adj[$j] }, [$i, $e];
	};
	foreach my $p (@nodes) {
		foreach my $c (@{ $p->{"children"} }) {
			&$link($idx{"$c"}, $idx{"$p"}, { "len" => $c->{"len"}, "sup" => ($c->{"leaf"} ? undef : $c->{"support"}) });
		}
	}
	# A bifurcating root is not a real node of the unrooted tree: join its two branches.
	if (@{ $root->{"children"} } == 2) {
		my $r = $idx{"$root"};
		my ($x, $y) = @{ $adj[$r] };
		my %e;
		my @l = grep { defined && $_ ne "" } ($x->[1]{"len"}, $y->[1]{"len"});
		$e{"len"} = @l ? &num(($x->[1]{"len"} || 0) + ($y->[1]{"len"} || 0)) : undef;
		($e{"sup"}) = grep { defined && $_ ne "" } ($x->[1]{"sup"}, $y->[1]{"sup"});
		$adj[$r] = [];
		foreach my $z ($x, $y) {
			my $other = ($z == $x) ? $y : $x;
			$adj[$z->[0]] = [ grep { $_->[0] != $r } @{ $adj[$z->[0]] } ];
			push @{ $adj[$z->[0]] }, [$other->[0], \%e];
		}
	}
	# Strains on the side of 'u' when the branch u-v is cut.
	my %memo;
	my $side;
	$side = sub {
		my ($u, $v) = @_;
		return $memo{"$u,$v"} if $memo{"$u,$v"};
		my @l = $nodes[$u]{"leaf"} ? ($nodes[$u]{"name"}) : map { @{ $side->($_->[0], $u) } } grep { $_->[0] != $v } @{ $adj[$u] };
		return $memo{"$u,$v"} = \@l;
	};
	my $nOg = scalar @$og;
	my ($inU, $outV, $edge);
	EDGE: foreach my $u (0..$#nodes) {
		foreach my $nb (@{ $adj[$u] }) {
			my $l = $side->($nb->[0], $u);
			if (@$l == $nOg && !grep { !$og{$_} } @$l) {
				($inU, $outV, $edge) = ($u, $nb->[0], $nb->[1]);
				last EDGE;
			}
		}
	}
	&fail("The outgroup (".join(", ", @$og).") does not form a clade in the unrooted tree, so the tree cannot be rooted with it.") unless defined $edge;
	my $copy;
	$copy = sub {
		my ($x, $from, $e) = @_;
		my $o = $nodes[$x];
		my %n = ("leaf" => $o->{"leaf"}, "name" => $o->{"name"}, "len" => $e->{"len"}, "children" => []);
		$n{"support"} = $e->{"sup"} unless $o->{"leaf"};
		foreach my $nb (@{ $adj[$x] }) {
			push @{ $n{"children"} }, $copy->($nb->[0], $x, $nb->[1]) unless $nb->[0] == $from;
		}
		return \%n;
	};
	my $half = (defined $edge->{"len"} && $edge->{"len"} ne "") ? &num($edge->{"len"} / 2) : undef;
	my %e2 = ("len" => $half, "sup" => $edge->{"sup"});
	my $in  = $copy->($inU,  $outV, \%e2);
	my $out = $copy->($outV, $inU,  \%e2);
	return { "leaf" => 0, "support" => "", "children" => [$in, $out] };
}

sub num {
	my $x = sprintf("%.10f", shift);
	$x =~ s/0+$//;
	$x =~ s/\.$//;
	return $x;
}

# Options of this script in JSON format, for EasyCGTree_GUI (JSON::PP is a core module since Perl 5.14).
sub printOptionsJson {
	require JSON::PP;
	my @o = (
		{ "name" => "input", "type" => "dir", "default" => undef, "required" => 1 },
		{ "name" => "outdir", "type" => "dir", "default" => undef },
		{ "name" => "hmm", "type" => "text", "default" => undef },
		{ "name" => "tree_app", "type" => "choice", "default" => "iqtree", "values" => ["iqtree", "fasttree"] },
		{ "name" => "tree", "type" => "file", "default" => undef },
		{ "name" => "outgroup", "type" => "text", "default" => undef },
		{ "name" => "mode", "type" => "choice", "default" => "exclusive", "values" => ["exclusive", "strict"] },
		{ "name" => "thread", "type" => "int", "default" => 4, "min" => 1 },
	);
	print JSON::PP->new->canonical(1)->encode({ "script" => "EasyCGTree_SpecificSNP.pl", "version" => $VERSION, "options" => \@o }), "\n";
}
