#!/usr/bin/env bash
# Builds the altered FluidX3D of this repository. It never starts the solver and never touches another build.
# usage: ./build.sh [VARIANT] [extra compiler flags...]
#   ./build.sh                                                     -> bin/FluidX3D-wind                   (ground as upstream: bounce-back)
#   ./build.sh wall-model -DGROUND_WALL_MODEL                      -> bin/FluidX3D-wind-wall-model        (wall-stress ground, see src/defines.hpp)
#   ./build.sh wall-model-fp32 -DGROUND_WALL_MODEL -DFP32_STORAGE  -> bin/FluidX3D-wind-wall-model-fp32   (the same with 32-bit storage)
# The variant only names the binary and its object folder; the flags decide what is built. Objects go to temp/<variant>/,
# so variants do not mix. Portable: macOS links the OpenCL framework, Linux links libOpenCL
# (the bundled x86-64 copy in src/OpenCL/lib is on the linker path; on other machines install ocl-icd-opencl-dev). Needs g++ (or clang++ as g++) with C++17.
set -euo pipefail
cd "$(dirname "$0")"
VARIANT="${1:-}"; [ $# -gt 0 ] && shift
NAME="FluidX3D-wind${VARIANT:+-$VARIANT}"
OBJ="temp/${VARIANT:-default}"
CXX="${CXX:-g++}"
CFLAGS="-std=c++17 -pthread -O -Wno-comment -Wno-deprecated-declarations -I./src/OpenCL/include"
case "$(uname -s)" in
	Darwin) LIBS="-framework OpenCL" ;;
	*) LIBS="-L./src/OpenCL/lib -lOpenCL" ;;
esac
mkdir -p "$OBJ" bin
pids=()
for f in src/*.cpp; do
	o="$OBJ/$(basename "${f%.cpp}").o"
	# rebuild when the source, any header or this script is newer than the object, or the flags changed
	if [ ! -f "$o" ] || [ -n "$(find src build.sh -newer "$o" \( -name '*.cpp' -o -name '*.hpp' -o -name build.sh \) -print -quit)" ] || [ "$(cat "$OBJ/flags.txt" 2>/dev/null)" != "$*" ]; then
		$CXX -c "$f" -o "$o" $CFLAGS "$@" & pids+=($!)
	fi
done
for p in "${pids[@]:-}"; do [ -z "$p" ] || wait "$p"; done
printf '%s' "$*" > "$OBJ/flags.txt"
$CXX "$OBJ"/*.o -o "bin/$NAME" -pthread $LIBS
echo "built bin/$NAME (flags: ${*:-none})"
