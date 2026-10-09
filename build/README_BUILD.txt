Building EasyCGTree_GUI as a program (no Python needed by the users)
==================================================================

A program can only be built on the operating system it is made for: the Windows program on Windows,
the macOS app on a Mac, the Linux program on Linux. There are two ways:

A. On GitHub, all systems at once (recommended)
   The repository EasyCGTree5 contains two workflows (folder '.github/workflows'):
   - 'Build EasyCGTree_GUI' (build-gui.yml): only the GUI programs.
     'Actions' -> 'Build EasyCGTree_GUI' -> 'Run workflow'; after about 5-10 minutes, download the programs
     under 'Artifacts' of the finished run: EasyCGTree_GUI-linux-x86_64, -windows-x64, -macos-arm64
     (Apple silicon), -macos-x86_64 (Intel Mac).
   - 'Build release packages' (build-release.yml): the four complete program packages
     (EasyCGTree5-Linux-x64.tar.gz, -Windows-x64.zip, -macOS-arm64.zip, -macOS-x64.zip) with the scripts,
     the GUI, the third-party programs, HMM/README.txt, the manuals, README.md and LICENSE.
     The macOS programs of 'bin' are built from their source code during the run (build_macos_bin.sh).
     The Linux and Windows programs are taken from a release with the tag 'bin-programs' of the repository,
     with the assets 'bin-linux-x64.zip' and 'bin-windows-x64.zip' (programs at the top level of each zip).
     Make that release once ('Releases' -> 'Draft a new release', tag 'bin-programs', upload the two zips,
     tick 'Set as a pre-release'), then run the workflow and attach the packages from 'Artifacts' to the
     release of EasyCGTree.

B. On each computer
   Windows: double-click build_windows.bat   -> dist\EasyCGTree_GUI.exe
   macOS:   bash build_macos.sh               -> dist/EasyCGTree_GUI.app
   Linux:   bash build_linux.sh               -> dist/EasyCGTree_GUI
   (Python 3.9 or newer and internet access are needed only for building.)

The third-party programs for macOS ('bin')
   On a Mac (Xcode command line tools and Homebrew needed): bash build_macos_bin.sh  -> bin-macos-<arch>/
   It builds HMMER, Prodigal, MUSCLE, trimAl, FastTree and wASTRAL from their source code (linked only against
   macOS libraries) and downloads IQ-TREE; copy the files into the folder 'bin' of EasyCGTree.

Using the program
   Copy the program into the EasyCGTree directory (next to EasyCGTree.pl, 'bin' and 'HMM'); it then finds
   EasyCGTree by itself. Elsewhere, set the EasyCGTree directory in Settings - Paths and programs.
   The program contains Python and Qt, but not Perl: Perl is still needed to run the EasyCGTree scripts
   (Linux/macOS: already installed; Windows: e.g. Strawberry Perl, or a portable Perl set in the settings).

Notes
   - Size: about 50-80 MB per program (Python and the Qt libraries are included).
   - Windows/Linux programs are single files; they start a little slower (a few seconds) because they unpack
     themselves first.
   - macOS: the app is not signed by Apple. At the first start, right-click it and choose 'Open' (or allow it
     in System Settings - Privacy & Security).
   - Windows: SmartScreen may warn about an unknown program ('More info' -> 'Run anyway').
   - Linux: the program built on Ubuntu 22.04 runs on distributions with glibc 2.35 or newer (Ubuntu 22.04+,
     Debian 12+, Fedora 36+ ...). For older systems, build it there with build_linux.sh.
