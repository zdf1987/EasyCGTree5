#!/usr/bin/perl
#
# BuildHMM - build a profile HMM database for EasyCGTree from gene families
# Part of EasyCGTree 5.0, by Dao-Feng Zhang
#
# Each FASTA file in the input directory holds the protein sequences of one gene family. The families are
# aligned with MUSCLE (unless they are already aligned), a profile HMM is built for each with hmmbuild, and
# all HMMs are written into one file in the 'HMM' directory of EasyCGTree, ready for 'EasyCGTree.pl -hmm'.
# Only core Perl modules are used; the script can be called from any directory.
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
my $UPDATE  = "2026-10-07";
my $START_TIME = time;
my ($sec,$min,$hour,$mday,$mon,$year) = localtime($START_TIME);
my $STAMP     = sprintf("%04d-%02d-%02d-%02d-%02d", $year+1900, $mon+1, $mday, $hour, $min);
my $START_STR = sprintf("%02d:%02d:%02d, %04d-%02d-%02d", $hour, $min, $sec, $year+1900, $mon+1, $mday);

my $IS_WIN  = ($^O eq 'MSWin32');
my $EXE     = $IS_WIN ? ".exe" : "";
my $BIN_DIR = (defined $ENV{"ECG_BIN"} && $ENV{"ECG_BIN"} ne "") ? $ENV{"ECG_BIN"} : "$RealBin/bin";   # ECG_BIN: set by EasyCGTree_GUI
my $HMM_DIR = "$RealBin/HMM";

my @usage = qq(
====== BuildHMM (EasyCGTree) ======
     EasyCGTree $VERSION, by Dao-Feng Zhang
     Update $UPDATE

Builds a profile HMM database from gene families for EasyCGTree ('-hmm').

Usage: perl BuildHMM.pl -gc <dir> [Options]

Essential Options:
-gc <String>
	Input directory of gene families: one FASTA file of protein sequences per gene family.
	The file name (without extension) becomes the name of the gene (e.g. 'rpoB.fas' -> 'rpoB').
	White space and the characters / \\ : * ? " < > | ( ) [ ] , ; are replaced by '-' in the names;
	the input files are not changed.

Optional Options:
-outdir <String>
	Directory for the summary, log file and working directory. [default: the directory containing '-gc']
	(The HMM database is always written to the 'HMM' directory of EasyCGTree.)
-name <String>
	Name of the HMM database: 'HMM/<name>.hmm', used as 'EasyCGTree.pl -hmm <name>'.
	[default: the name of the input directory]
-force <no value required>
	Overwrite 'HMM/<name>.hmm' if it already exists.
-aln <no value required>
	The files are already aligned (aligned FASTA); they are not aligned again.
-super5 <Int>
	Families with at least this number of sequences are aligned with 'muscle -super5' (faster, for large
	families); smaller ones with 'muscle -align'. [default: 1000]
-thread <Int>
	Number of threads used by MUSCLE. [default: 4]
-help <no value required>
	Display this message.
-options_json <no value required>
	Print the options of this script in JSON format (used by EasyCGTree_GUI) and exit.

Output:
'$HMM_DIR/<name>.hmm'               the HMM database
'<gc>.BuildHMM_summary.txt'          sequences, alignment length and HMM length of each family
'<gc>_BuildHMM_TEM/'                 cleaned sequences, alignments and the HMM of each family
(the summary, log and working directory are written next to the input directory)

);

########## Options ##########
my %opt;
GetOptions(\%opt, "gc:s", "outdir:s", "name:s", "force!", "aln!", "super5:i", "thread:i", "help!", "options_json!")
	or die "@usage\n\nERROR: Unrecognized option (see above).\n\n";
if ($opt{"options_json"}) {
	&printOptionsJson();
	exit 0;
}
if (!%opt || $opt{"help"}) {
	print @usage;
	exit;
}
&usageError("a directory must be specified for option '-gc'.") unless defined $opt{"gc"} && $opt{"gc"} =~ /\S/;
my $thread = 4;
if (exists $opt{"thread"}) {
	&usageError("'-thread' should be an integer >= 1.") unless defined $opt{"thread"} && $opt{"thread"} >= 1;
	$thread = $opt{"thread"};
}
{
	my $ncpu = &cpuCount();
	if ($ncpu && $thread > $ncpu) {
		print "\nWARNING: '-thread $thread' is larger than the number of CPU cores of this computer ($ncpu); $ncpu threads are used.\n" if exists $opt{"thread"};
		$thread = $ncpu;
	}
}
my $super5 = 1000;
if (exists $opt{"super5"}) {
	&usageError("'-super5' should be an integer >= 2.") unless defined $opt{"super5"} && $opt{"super5"} >= 2;
	$super5 = $opt{"super5"};
}

(my $inputArg = $opt{"gc"}) =~ s/[\/\\]+$//;
my $inputDir = File::Spec->rel2abs($inputArg);
&usageError("The input directory '$inputDir' does not exist.") unless -d $inputDir;
$inputDir = abs_path($inputDir);
my $inputName = basename($inputDir);
my $outDir    = dirname($inputDir);
if (defined $opt{"outdir"} && $opt{"outdir"} ne "") {   # '-outdir': all results go there (same layout)
	(my $od = $opt{"outdir"}) =~ s/[\/\\]+$// unless $opt{"outdir"} =~ /^[\/\\]+$/;
	$od = $opt{"outdir"} unless defined $od;
	$outDir = File::Spec->rel2abs($od);
	&usageError("'-outdir $opt{'outdir'}' is a file, not a directory.") if -e $outDir && !-d $outDir;
	$outDir = abs_path($outDir) if -d $outDir;
}
if (!-d $outDir) {
	mkpath($outDir);
	&usageError("The output directory '$outDir' can't be created.") unless -d $outDir;
}
my $TEMdir    = "$outDir/${inputName}_BuildHMM_TEM";

my $dbName;
if (defined $opt{"name"} && $opt{"name"} ne "") {
	$dbName = $opt{"name"};
} else {
	# default: the name of the input directory, with characters other than letters, digits, '.', '_'
	# and '-' replaced by '_' (e.g. spaces)
	($dbName = $inputName) =~ s/[^\w.\-]+/_/g;
	$dbName =~ s/^_+|_+$//g;
	$dbName = "myHMM" if $dbName eq "";
}
$dbName =~ s/\.hmm$//i;
&usageError("'-name $dbName': the name may only contain letters, digits, '.', '_' and '-'.") unless $dbName =~ /^[\w.\-]+$/;
my $dbFile = "$HMM_DIR/$dbName.hmm";
if (-e $dbFile && !$opt{"force"}) {
	&usageError("'$dbFile' already exists. Please choose another name with '-name', or use '-force' to overwrite it.");
}

open(LOG, ">", "$outDir/$inputName.BuildHMM_$STAMP.log") or die "Can't open the log file in '$outDir': $!\n";
$SIG{__DIE__} = sub { print LOG "\nERROR: $_[0]" if defined fileno(LOG); };

&msg("\n====== BuildHMM (EasyCGTree) ======\n	EasyCGTree $VERSION\nby Dao-Feng Zhang\n\n");
&msg("Job Started at: $START_STR\n\nInput directory: $inputDir\nHMM database: $dbFile\n".($opt{"aln"} ? "The input files are treated as alignments (-aln).\n" : "")."\n");


############### Step 1: read and check the gene families ###############
&msg("#============= Step 1: Check the gene families ==============#\n\n");
opendir(my $dh, $inputDir) or &fail("Can't open '$inputDir': $!");
my @files = sort grep { !/^\./ && -f "$inputDir/$_" } readdir($dh);
closedir($dh);
&fail("The input directory '$inputDir' does not contain any files.") unless @files;

my (@fam, %nameOf, @problems);
foreach my $f (@files) {
	my $name = $f;
	$name =~ s/\.(fas|fasta|fa|faa|afa|aln|fst|txt|seq|pep)$//i;
	$name =~ s/[\s\/\\:*?"<>|()\[\],;]+/-/g;
	$name =~ s/^-+|-+$//g;
	if ($name eq "") {
		push @problems, "'$f': no usable name can be derived from the file name.";
		next;
	}
	if ($nameOf{lc $name}) {
		push @problems, "'$nameOf{lc $name}' and '$f' would both get the gene name '$name'.";
		next;
	}
	$nameOf{lc $name} = $f;

	my ($ids, $seqs, $bad, $renamed) = &readFamily("$inputDir/$f");
	if ($bad) {
		push @problems, "'$f': $bad";
		next;
	}
	my $all = join("", map { $seqs->{$_} } @$ids);
	(my $res = $all) =~ s/[-.]//g;
	if (!length $res) {
		push @problems, "'$f': no sequence was found.";
		next;
	}
	my $acgt = ($res =~ tr/ACGTUN//);
	if ($acgt / length($res) > 0.9) {
		push @problems, "'$f' contains nucleotide sequences; protein sequences are needed (EasyCGTree searches proteomes).";
		next;
	}
	if ($opt{"aln"}) {
		my %len = map { length($seqs->{$_}) => 1 } @$ids;
		if (keys %len > 1) {
			push @problems, "'$f' is not an alignment (the sequences have different lengths), but '-aln' was given.";
			next;
		}
	}
	push @fam, { "file" => $f, "name" => $name, "ids" => $ids, "seqs" => $seqs, "n" => scalar @$ids, "renamedIds" => $renamed };
}
&fail("Problems were found in the input files:\n	".join("\n	", @problems)) if @problems;
my $nFam = scalar @fam;
&msg("$nFam gene families were found:\n");
foreach my $fm (@fam) {
	&msg("	$fm->{'name'}	$fm->{'n'} sequences	(file '$fm->{'file'}')".($fm->{'renamedIds'} ? "; duplicate sequence IDs were renamed" : "")."\n");
}
my @few = grep { $_->{"n"} < 3 } @fam;
&msg("\nWARNING: ".scalar(@few)." families have fewer than 3 sequences (".join(", ", map { $_->{'name'} } @few)."); their HMMs will be less sensitive.\n") if @few;
&msg("\n");

&removePath($TEMdir) if -e $TEMdir;
&makeDir("$TEMdir/TEM0_Sequences");
&makeDir("$TEMdir/TEM1_Alignment");
&makeDir("$TEMdir/TEM2_HMMs");

# Cleaned copies: upper case, no white space; '*' at the end removed, other '*' -> 'X'; duplicate IDs made unique.
foreach my $fm (@fam) {
	my $out = "$TEMdir/TEM0_Sequences/$fm->{'name'}.fas";
	open(my $fo, ">", $out) or &fail("Can't open '$out': $!");
	foreach my $id (@{ $fm->{"ids"} }) {
		my $s = $fm->{"seqs"}{$id};
		$s =~ s/\./-/g;
		$s =~ s/-//g unless $opt{"aln"};
		print $fo ">$id\n$s\n";
	}
	close $fo;
}


############### Step 2: alignment ###############
&msg("#============= Step 2: Align the gene families ==============#\n\n");
my $muscle;
$muscle = &tool("muscle5") unless $opt{"aln"};
my $num = 0;
foreach my $fm (@fam) {
	$num++;
	my $in  = "$TEMdir/TEM0_Sequences/$fm->{'name'}.fas";
	my $aln = "$TEMdir/TEM1_Alignment/$fm->{'name'}.afa";
	if ($opt{"aln"}) {
		copy($in, $aln) or &fail("Can't copy '$in': $!");
		$fm->{"method"} = "given";
	} elsif ($fm->{"n"} == 1) {
		copy($in, $aln) or &fail("Can't copy '$in': $!");
		$fm->{"method"} = "single sequence";
	} else {
		my $an = ($fm->{"n"} >= $super5) ? "-super5" : "-align";
		print "\r	Aligning ($num/$nFam): $fm->{'name'} ($fm->{'n'} sequences, muscle $an)          ";
		my $rc = &runCapture("$TEMdir/TEM1_Alignment/$fm->{'name'}.muscle.log", $muscle, $an, $in, "-output", $aln, "-threads", $thread);
		unless ($rc == 0 && -s $aln) {
			&fail("MUSCLE failed to align '$in' (exit code $rc).".($rc >= 128 ? " The process was killed, probably because it ran out of memory." : "")."\nPlease see '$TEMdir/TEM1_Alignment/$fm->{'name'}.muscle.log'.");
		}
		$fm->{"method"} = "muscle $an";
	}
	my ($ids, $seqs) = &readFamily($aln);
	$fm->{"alnLen"} = @$ids ? length($seqs->{$ids->[0]}) : 0;
}
print "\n" unless $opt{"aln"};
&msg(($opt{"aln"} ? "The given alignments are used." : "$nFam families have been aligned.")." Alignments: '$TEMdir/TEM1_Alignment'\n\n");


############### Step 3: profile HMMs ###############
&msg("#============= Step 3: Build the profile HMMs ==============#\n\n");
my $hmmbuild = &tool("hmmbuild");
$num = 0;
foreach my $fm (@fam) {
	$num++;
	my $aln = "$TEMdir/TEM1_Alignment/$fm->{'name'}.afa";
	my $hmm = "$TEMdir/TEM2_HMMs/$fm->{'name'}.hmm";
	print "\r	Building HMM ($num/$nFam): $fm->{'name'}                    ";
	my $rc = &runCapture("$TEMdir/TEM2_HMMs/$fm->{'name'}.hmmbuild.log", $hmmbuild, "--amino", "--informat", "afa", "-n", $fm->{"name"}, $hmm, $aln);
	unless ($rc == 0 && -s $hmm) {
		&fail("hmmbuild failed for '$aln' (exit code $rc).\nPlease see '$TEMdir/TEM2_HMMs/$fm->{'name'}.hmmbuild.log'.");
	}
	open(my $fh, "<", $hmm) or &fail("Can't open '$hmm': $!");
	while (my $l = <$fh>) {
		if ($l =~ /^LENG\s+(\d+)/) { $fm->{"hmmLen"} = $1; last; }
	}
	close $fh;
}
print "\n";

# All HMMs in one file.
&makeDir($HMM_DIR);
my $tmpDb = "$TEMdir/$dbName.hmm";
open(my $db, ">", $tmpDb) or &fail("Can't open '$tmpDb': $!");
foreach my $fm (@fam) {
	open(my $fh, "<", "$TEMdir/TEM2_HMMs/$fm->{'name'}.hmm") or &fail("Can't open the HMM of '$fm->{'name'}': $!");
	print $db $_ while <$fh>;
	close $fh;
}
close $db;

# The database must be readable by the hmmsearch of EasyCGTree (HMMs of HMMER 3.1+ cannot be read by HMMER 3.0).
my $hmmsearch = &tool("hmmsearch");
my $testOut = "$TEMdir/test_hmmsearch.tbl";
my $testSeq = "$TEMdir/test_sequences.fas";
open(my $ts, ">", $testSeq) or &fail("Can't open '$testSeq': $!");
foreach my $fm (@fam) {
	(my $s = $fm->{"seqs"}{ $fm->{"ids"}[0] }) =~ s/[-.]//g;
	print $ts ">$fm->{'name'}\n$s\n";
}
close $ts;
my $rc = &runCapture("$TEMdir/test_hmmsearch.log", $hmmsearch, "--tblout", $testOut, "--cpu", 1, $tmpDb, $testSeq);
if ($rc != 0) {
	my $m = "";
	if (open(my $lh, "<", "$TEMdir/test_hmmsearch.log")) { local $/; $m = <$lh>; close $lh; }
	my $hint = ($m =~ /Unrecognized format/i)
		? "The HMMs were built by a newer HMMER than the 'hmmsearch' in '$BIN_DIR' can read.\nPlease use 'hmmbuild' and 'hmmsearch' of the same HMMER version (e.g. replace both with HMMER 3.4)."
		: "Please see '$TEMdir/test_hmmsearch.log'.";
	&fail("The new HMM database could not be read by '$hmmsearch' (exit code $rc).\n$hint");
}
copy($tmpDb, $dbFile) or &fail("Can't write '$dbFile': $!");
my $nNames = 0;
if (open(my $fh, "<", $dbFile)) { while (<$fh>) { $nNames++ if /^NAME\s/; } close $fh; }
&fail("'$dbFile' contains $nNames HMMs instead of $nFam.") unless $nNames == $nFam;

# Summary table
my $summary = "$outDir/$inputName.BuildHMM_summary.txt";
open(my $so, ">", $summary) or &fail("Can't open '$summary': $!");
print $so "# HMM database '$dbFile' built from '$inputDir'.\n";
print $so join("\t", qw(Gene Input_file Sequences Alignment Alignment_length HMM_length Note)), "\n";
foreach my $fm (@fam) {
	my @note;
	push @note, "fewer than 3 sequences" if $fm->{"n"} < 3;
	push @note, "duplicate sequence IDs renamed" if $fm->{"renamedIds"};
	print $so join("\t", $fm->{"name"}, $fm->{"file"}, $fm->{"n"}, $fm->{"method"}, $fm->{"alnLen"}, $fm->{"hmmLen"} // "NA", (@note ? join("; ", @note) : "-")), "\n";
}
close $so;

&msg("$nFam profile HMMs have been built and written to '$dbFile'.\nSummary of the families: '$summary'\n\n");
&msg("Use it with:  perl EasyCGTree.pl -input <your genomes> -hmm $dbName\n\n");
my $el = time - $START_TIME;
&msg(sprintf("Running time: %d d %d h %d min %d sec.\n\nBuildHMM (EasyCGTree %s)\n\n", int($el/86400), int(($el%86400)/3600), int(($el%3600)/60), $el%60, $VERSION));
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

# Read one FASTA file. Returns (\@ids, \%seq, $note): IDs are the first word of the header (duplicates get
# '_2', '_3' ...); sequences in upper case without white space, '*' at the end removed and other '*' -> 'X'.
# Also returns an error message (undef if the file is fine) and 1 if duplicate IDs were renamed.
sub readFamily {
	my $f = shift;
	my (@ids, %seq, %count, $id, $renamed);
	open(my $fh, "<", $f) or return ([], {}, "can't be opened: $!");
	while (my $l = <$fh>) {
		$l =~ s/\r?\n$//;
		next unless $l =~ /\S/;
		if ($l =~ s/^\>//) {
			($id) = split /\s+/, $l;
			$id = "seq" unless defined $id && $id ne "";
			$count{$id}++;
			if ($count{$id} > 1) {
				$id .= "_$count{$id}";
				$renamed = 1;
			}
			push @ids, $id;
			$seq{$id} = "";
		} elsif (defined $id) {
			$l =~ s/\s+//g;
			$seq{$id} .= uc($l);
		} else {
			close $fh;
			return ([], {}, "does not look like a FASTA file (no '>' before the first sequence line).");
		}
	}
	close $fh;
	return ([], {}, "does not contain any sequence.") unless @ids;
	foreach my $i (@ids) {
		$seq{$i} =~ s/\*+$//;
		$seq{$i} =~ s/\*/X/g;
		return ([], {}, "the sequence '$i' is empty.") if $seq{$i} !~ /[A-Z]/;
	}
	return (\@ids, \%seq, undef, $renamed ? 1 : 0);
}

# Options of this script in JSON format, for EasyCGTree_GUI (JSON::PP is a core module since Perl 5.14).
sub printOptionsJson {
	require JSON::PP;
	my @o = (
		{ "name" => "gc", "type" => "dir", "default" => undef, "required" => 1 },
		{ "name" => "outdir", "type" => "dir", "default" => undef },
		{ "name" => "name", "type" => "text", "default" => undef },
		{ "name" => "force", "type" => "bool", "default" => 0 },
		{ "name" => "aln", "type" => "bool", "default" => 0 },
		{ "name" => "super5", "type" => "int", "default" => 1000, "min" => 2 },
		{ "name" => "thread", "type" => "int", "default" => 4, "min" => 1 },
	);
	print JSON::PP->new->canonical(1)->encode({ "script" => "BuildHMM.pl", "version" => $VERSION, "options" => \@o }), "\n";
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
