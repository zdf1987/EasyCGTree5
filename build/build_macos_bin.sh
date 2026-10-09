#!/bin/bash
# Build the third-party programs of EasyCGTree ('bin' folder) for macOS.
# Run on a Mac with the Xcode command line tools (xcode-select --install) and Homebrew (https://brew.sh).
# An Apple silicon Mac (M1, M2 ...) builds programs for Apple silicon, an Intel Mac for Intel.
#   bash build_macos_bin.sh [output folder]        (default: ./bin-macos-<arch>)
# The programs are linked only against macOS system libraries, so they run on Macs without Homebrew.
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
ARCH="$(uname -m)"                                   # arm64 or x86_64
OUT="${1:-$HERE/bin-macos-$ARCH}"
mkdir -p "$OUT"; OUT="$(cd "$OUT" && pwd)"
W="$(mktemp -d)"; trap 'rm -rf "$W"' EXIT
export MACOSX_DEPLOYMENT_TARGET=11.0
J="$(sysctl -n hw.ncpu)"

brew list libomp >/dev/null 2>&1 || brew install libomp
OMP="$(brew --prefix libomp)"
if [ ! -f "$OMP/lib/libomp.a" ]; then
  echo "ERROR: $OMP/lib/libomp.a not found (a static OpenMP library is needed for MUSCLE and FastTree)." >&2; exit 1
fi
OMPFLAGS=(-Xpreprocessor -fopenmp -I"$OMP/include")
cd "$W"

echo "== HMMER 3.4"
curl -fsSL -o hmmer.tar.gz http://eddylab.org/software/hmmer/hmmer-3.4.tar.gz
tar xzf hmmer.tar.gz
(cd hmmer-3.4 && ./configure --quiet && make -j"$J" >/dev/null)
cp hmmer-3.4/src/hmmsearch hmmer-3.4/src/hmmbuild "$OUT/"

echo "== Prodigal 2.6.3"
git clone -q --depth 1 -b v2.6.3 https://github.com/hyattpd/Prodigal
(cd Prodigal && make CC=clang >/dev/null)
cp Prodigal/prodigal "$OUT/"

echo "== MUSCLE 5.3"
git clone -q --depth 1 -b v5.3 https://github.com/rcedgar/muscle
(
  cd muscle/src
  echo '"v5.3"' > gitver.txt
  sed -i '' 's/MemBytesToStr(n)/MemBytesToStr((double) n)/' test_malloc.cpp
  FILES=$(grep -o 'ClCompile Include="[^"]*"' muscle.vcxproj | sed 's/.*="//; s/"//')
  clang++ -O3 -std=c++17 -DNDEBUG "${OMPFLAGS[@]}" $FILES "$OMP/lib/libomp.a" -o "$OUT/muscle5"
)

echo "== trimAl 1.5.1"
git clone -q --depth 1 -b v1.5.1 https://github.com/inab/trimal
(cd trimal/source && clang++ -O2 -o "$OUT/trimal" main.cpp alignment.cpp statisticsGaps.cpp utils.cpp \
   similarityMatrix.cpp statisticsConservation.cpp sequencesMatrix.cpp compareFiles.cpp)

echo "== FastTree 2.2 (FastTreeMP)"
git clone -q --depth 1 https://github.com/morgannprice/fasttree
clang -O3 -DOPENMP "${OMPFLAGS[@]}" -funsafe-math-optimizations fasttree/FastTree.c "$OMP/lib/libomp.a" -lm \
  -o "$OUT/FastTreeMP"

echo "== ASTER 1.25 (wASTRAL -> astral-weighted)"
git clone -q --depth 1 -b v1.25 https://github.com/chaoszhang/ASTER
clang++ -std=gnu++11 -O3 -pthread ASTER/src/astral-hybrid.cpp -o "$OUT/astral-weighted"

echo "== IQ-TREE 3.1.4 (official build)"
if [ "$ARCH" = "arm64" ]; then IQ=iqtree-3.1.4-macOS-arm; else IQ=iqtree-3.1.4-macOS-intel; fi
curl -fsSL -o iq.zip "https://github.com/iqtree/iqtree3/releases/download/v3.1.4/$IQ.zip"
unzip -q iq.zip
cp "$(find . -path "*$IQ*/bin/iqtree3" | head -1)" "$OUT/"

cp "$HERE/../bin/tree_app-options.txt" "$HERE/../bin/THIRD_PARTY.txt" "$OUT/" 2>/dev/null || true
chmod +x "$OUT"/hmmsearch "$OUT"/hmmbuild "$OUT"/prodigal "$OUT"/muscle5 "$OUT"/trimal "$OUT"/FastTreeMP \
         "$OUT"/astral-weighted "$OUT"/iqtree3

echo "== Check: architecture and libraries"
BAD=0
for f in hmmsearch hmmbuild prodigal muscle5 trimal FastTreeMP astral-weighted iqtree3; do
  echo "$f: $(lipo -archs "$OUT/$f")"
  if otool -L "$OUT/$f" | tail -n +2 | grep -v -E '^\s*/usr/lib/|^\s*/System/' ; then
    echo "  ^ ERROR: $f needs a library that is not part of macOS" >&2; BAD=1
  fi
done
"$OUT/hmmsearch" -h | sed -n 2p
"$OUT/muscle5" -version | head -1
"$OUT/trimal" --version
"$OUT/iqtree3" --version | head -1
[ $BAD = 0 ] || exit 1
echo
echo "Done: $OUT"
