EasyCGTree_GUI 1.0 - graphical interface for EasyCGTree 5.0
===========================================================

The program packages (Releases page) contain EasyCGTree_GUI as a ready-to-use program. This file describes
how to run it from its source code (EasyCGTree_GUI.py) and gives an overview; see the manual (doc/) for details.

Installation
------------
1. Python 3.9 or newer, and PySide6:       pip install PySide6
2. Perl:  Linux/macOS: already installed.  Windows: e.g. Strawberry Perl (https://strawberryperl.com).
3. Put EasyCGTree_GUI.py into the EasyCGTree directory (next to EasyCGTree.pl, 'bin' and 'HMM'),
   together with the five version-5.0 scripts:
   EasyCGTree.pl, EasyCGTree_SNP.pl, EasyCGTree_SpecificSNP.pl, BuildHMM.pl, Gene_Prevelence.pl
4. Start:   python EasyCGTree_GUI.py      (Windows: py EasyCGTree_GUI.py, or double-click if .py opens Python)

Program without Python: the folder 'build' contains scripts that make EasyCGTree_GUI a program for Windows,
macOS or Linux (with PyInstaller), and a GitHub workflow that builds all of them; see build/README_BUILD.txt.

The GUI can also be kept elsewhere: set the EasyCGTree directory in Settings - Paths and programs.

Overview
--------
- Three tabs: Main pipeline, SNP analysis, HMM building. In each tab: global settings first, then the
  steps in the order they are run. Options not used by the chosen task/steps are greyed out.
  '?' next to an option explains it.
- Gene prevalence (Gene_Prevelence.pl) is part of the Main pipeline tab (after the run, or alone) and of
  the SNP analysis tab.
- Command line: the command(s) that 'Run' will execute; 'Concise' leaves out default values, 'Full'
  shows all options. 'Copy' copies them for a terminal (Windows: PowerShell syntax).
- Output directory: all results go there ('-outdir'); empty = the folder that contains the input folder.
- Files: the buttons 'Output' / 'Input' (HMM tab: 'Gene families' / 'Output' / 'HMM directory') switch
  the file view; the path is shown next to them. Double-click opens folders, tables (shown in the GUI) and
  other files (system application).
- Status follows the chosen place: Output = what has been done in '<input>_TEM' (EasyCGTree_record.txt);
  Input = genomes, sequence types, number and total length of the sequences per file.
- Task: tick the tasks; ticking two tasks also ticks the ones between them; all five = 'all'. Several tasks
  (not all) are run as one EasyCGTree command per task, one after the other (shown in Command line).
- If a CDS prediction with the same input files and settings already exists, a yellow notice is shown and
  'predict' is unticked; close the notice to tick it again.
- Threads: limited to the number of CPU cores of the computer.
- SNP analysis: the EasyCGTree settings are filled in from the record of the previous run (or the defaults).
- Mouse wheel: changes a value only after the field was clicked and while the pointer is on it; otherwise
  it scrolls the parameter list.
- 'Check input' (Main pipeline): shows how the input files will be grouped and prepared, without changing
  anything.
- Settings - Paths and programs: EasyCGTree directory, bin directory, Perl, and the path of each program
  (passed to the scripts as environment variables ECG_BIN / ECG_<PROGRAM>). 'Check installation' tests
  the scripts, Perl, every program (version), the HMM sets and tree_app-options.txt.
- Settings - Tree program command lines: edits bin/tree_app-options.txt (a backup .bak is kept).
- File - Save/Load parameters: parameter sets as JSON. At each start the analysis parameters have their default
  values; the environment (language, paths and programs, concise/full command, window size) is kept.
- Settings - Language: English / Chinese.
- Log: repeated empty lines are shown as one; 'Clear log' is in the log box.
- Help - About: the text is in the ABOUT block at the top of EasyCGTree_GUI.py and can be edited.
