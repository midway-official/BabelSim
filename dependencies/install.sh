#!/usr/bin/env bash
# Bootstrap the system toolchain, then build the bundled dependencies and BabelSim.
set -euo pipefail

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
PROJECT_DIR=$(cd -- "$SCRIPT_DIR/.." && pwd -P)
JOBS=4
CHECK_ONLY=0
DEPS_ONLY=0
CLEAN_BUILD=0

usage() {
    cat <<'EOF'
用法：./dependencies/install.sh [选项]

默认：检查系统依赖，在 Ubuntu/Debian 上安装缺少的包，再用本地源码编译到 build/。
请用普通用户运行；安装系统包时才会调用 sudo。

  --check       仅检查依赖和随附源码，不安装系统包、不构建项目
  --deps-only   检查并安装系统依赖，不构建项目
  --clean       构建前执行 make clean，从头编译 PETSc/Hypre 和 BabelSim
  --jobs N      依赖和项目的并行编译数，默认 4
  -h, --help    显示帮助

系统依赖齐全时无需 sudo 或联网。PETSc、Eigen、Hypre 使用 dependencies/ 内的源码。
EOF
}

fail() { printf '错误：%s\n' "$*" >&2; exit 1; }

while (($#)); do
    case "$1" in
        --check) CHECK_ONLY=1 ;;
        --deps-only) DEPS_ONLY=1 ;;
        --clean) CLEAN_BUILD=1 ;;
        --jobs)
            (($# >= 2)) || fail '--jobs 需要正整数参数'
            JOBS=$2
            [[ "$JOBS" =~ ^[1-9][0-9]*$ ]] || fail '--jobs 必须为正整数'
            shift ;;
        -h|--help) usage; exit 0 ;;
        *) fail "未知选项：$1（使用 --help 查看用法）" ;;
    esac
    shift
done
if ((CHECK_ONLY + DEPS_ONLY + CLEAN_BUILD > 1)); then
    fail '--check、--deps-only、--clean 不能同时使用'
fi

# Check the delivery before changing the system. Bash is sufficient to bootstrap Python.
for path in Makefile dependencies/build_petsc.py dependencies/sources.json \
            dependencies/petsc/configure dependencies/petsc/include/petscversion.h \
            dependencies/eigen/Eigen/Core dependencies/eigen/Eigen/Cholesky; do
    [[ -f "$PROJECT_DIR/$path" ]] || fail "随附源码或构建文件缺失：$path"
done

PROBE_DIR=$(mktemp -d "${TMPDIR:-/tmp}/babelsim-dependencies.XXXXXXXX")
trap 'rm -rf -- "$PROBE_DIR"' EXIT
MISSING=()
PACKAGES=()

missing() {
    MISSING+=("$1")
    printf '[缺少或不可用] %s\n' "$1"
    shift
    local package existing found
    for package in "$@"; do
        found=0
        for existing in "${PACKAGES[@]}"; do
            [[ "$package" != "$existing" ]] || found=1
        done
        ((found)) || PACKAGES+=("$package")
    done
}

check_tool() {
    local tool=$1
    shift
    if command -v "$tool" >/dev/null 2>&1; then
        printf '[通过] %s\n' "$tool"
    else
        missing "$tool" "$@"
    fi
}

# Compile and run small programs: having wrapper commands alone is insufficient.
cat > "$PROBE_DIR/mpi.c" <<'EOF'
#include <mpi.h>
#include <stdio.h>
int main(int argc, char **argv) {
    char version[MPI_MAX_LIBRARY_VERSION_STRING];
    int length = 0;
    if (MPI_Init(&argc, &argv) != MPI_SUCCESS) return 1;
    if (MPI_Get_library_version(version, &length) != MPI_SUCCESS) return 1;
    printf("%.*s\n", length, version);
    return MPI_Finalize() != MPI_SUCCESS;
}
EOF
cat > "$PROBE_DIR/mpi.cpp" <<'EOF'
#include <mpi.h>
#include <cstdio>
#include <filesystem>
#include <optional>
int main(int argc, char **argv) {
    char version[MPI_MAX_LIBRARY_VERSION_STRING];
    int length = 0;
    std::optional<std::filesystem::path> path{std::filesystem::current_path()};
    if (!path || MPI_Init(&argc, &argv) != MPI_SUCCESS) return 1;
    if (MPI_Get_library_version(version, &length) != MPI_SUCCESS) return 1;
    std::printf("%.*s\n", length, version);
    return MPI_Finalize() != MPI_SUCCESS;
}
EOF
cat > "$PROBE_DIR/mpi.f90" <<'EOF'
program check_mpi
  use mpi
  implicit none
  integer :: ierr, length
  character(len=MPI_MAX_LIBRARY_VERSION_STRING) :: version
  call MPI_Init(ierr)
  if (ierr /= MPI_SUCCESS) stop 1
  call MPI_Get_library_version(version, length, ierr)
  if (ierr /= MPI_SUCCESS) stop 1
  write(*, '(a)') version(:length)
  call MPI_Finalize(ierr)
  if (ierr /= MPI_SUCCESS) stop 1
end program check_mpi
EOF
cat > "$PROBE_DIR/blas.c" <<'EOF'
extern void dgemm_(const char *, const char *, const int *, const int *, const int *,
                   const double *, const double *, const int *, const double *,
                   const int *, const double *, double *, const int *);
extern void dgesv_(const int *, const int *, double *, const int *, int *,
                   double *, const int *, int *);
int main(void) {
    int n = 1, pivot = 0, info = -1;
    double a = 2, b = 3, c = 0, alpha = 1, beta = 0, rhs = 6;
    dgemm_("N", "N", &n, &n, &n, &alpha, &a, &n, &b, &n, &beta, &c, &n);
    dgesv_(&n, &n, &a, &n, &pivot, &rhs, &n, &info);
    return info != 0 || c != 6 || rhs != 3;
}
EOF

check_mpi() {
    local compiler=$1 source=$2 binary=$3
    shift 3
    if command -v "$compiler" >/dev/null 2>&1 && \
       "$compiler" "$@" "$PROBE_DIR/$source" -o "$PROBE_DIR/$binary" \
           > "$PROBE_DIR/$binary.log" 2>&1 && \
       "$PROBE_DIR/$binary" > "$PROBE_DIR/$binary.version" 2>> "$PROBE_DIR/$binary.log"; then
        printf '[通过] %s 编译、链接和 MPI 初始化\n' "$compiler"
    else
        missing "$compiler 编译、链接或 MPI 初始化" openmpi-bin libopenmpi-dev
        [[ ! -f "$PROBE_DIR/$binary.log" ]] || tail -n 8 "$PROBE_DIR/$binary.log" >&2
    fi
}

mpi_version() {
    local version
    version=$(< "$1")
    # MPI's Fortran bindings may pad the returned string with whitespace.
    while [[ "$version" == *[[:space:]] ]]; do version=${version%?}; done
    printf '%s' "$version"
}

check_dependencies() {
    MISSING=()
    PACKAGES=()
    rm -f -- "$PROBE_DIR/mpi-c.version" "$PROBE_DIR/mpi-cxx.version" "$PROBE_DIR/mpi-fortran.version"
    printf '\n检查系统构建依赖：\n'
    for tool in gcc g++ gcc-ar; do check_tool "$tool" build-essential; done
    check_tool gfortran gfortran
    check_tool pkg-config pkg-config
    check_tool tar tar
    check_tool gzip gzip
    if command -v make >/dev/null 2>&1 && [[ $(make --version 2>/dev/null) == *'GNU Make'* ]]; then
        printf '[通过] GNU Make\n'
    else
        missing 'GNU Make' build-essential
    fi
    if command -v python3 >/dev/null 2>&1 && \
       python3 -c 'import sys; sys.exit(sys.version_info < (3, 9))'; then
        printf '[通过] Python >= 3.9\n'
    else
        missing 'Python >= 3.9' python3
    fi
    check_tool mpirun openmpi-bin
    check_mpi mpicc mpi.c mpi-c
    check_mpi mpic++ mpi.cpp mpi-cxx -std=c++17 -DOMPI_SKIP_MPICXX=1 -DMPICH_SKIP_MPICXX=1
    check_mpi mpifort mpi.f90 mpi-fortran
    if [[ -s "$PROBE_DIR/mpi-c.version" && -s "$PROBE_DIR/mpi-cxx.version" && \
          -s "$PROBE_DIR/mpi-fortran.version" ]]; then
        if [[ $(mpi_version "$PROBE_DIR/mpi-c.version") == "$(mpi_version "$PROBE_DIR/mpi-cxx.version")" && \
              $(mpi_version "$PROBE_DIR/mpi-c.version") == "$(mpi_version "$PROBE_DIR/mpi-fortran.version")" ]]; then
            printf '[通过] C/C++/Fortran 的 MPI 实现和版本一致\n'
        else
            missing 'C/C++/Fortran 使用了不同的 MPI 实现或版本' openmpi-bin libopenmpi-dev
        fi
    fi
    if command -v gcc >/dev/null 2>&1 && \
       gcc "$PROBE_DIR/blas.c" -llapack -lblas -o "$PROBE_DIR/blas" \
           > "$PROBE_DIR/blas.log" 2>&1 && \
       "$PROBE_DIR/blas" >> "$PROBE_DIR/blas.log" 2>&1; then
        printf '[通过] BLAS/LAPACK 链接和数值调用\n'
    else
        missing 'BLAS/LAPACK 开发库或链接环境' libopenblas-dev
        [[ ! -f "$PROBE_DIR/blas.log" ]] || tail -n 8 "$PROBE_DIR/blas.log" >&2
    fi
}

check_dependencies
if ((${#MISSING[@]})); then
    printf '\n需要的系统包：%s\n' "${PACKAGES[*]}"
    ((CHECK_ONLY == 0)) || exit 1
    command -v apt-get >/dev/null 2>&1 || \
        fail '自动安装支持 Ubuntu/Debian 的 apt-get；请安装上述依赖后重新运行'
    PRIVILEGE=()
    if ((EUID != 0)); then
        command -v sudo >/dev/null 2>&1 || fail '安装系统包需要 root 或 sudo'
        PRIVILEGE=(sudo)
    fi
    printf '\n安装缺少的系统依赖：\n'
    "${PRIVILEGE[@]}" apt-get update
    "${PRIVILEGE[@]}" env DEBIAN_FRONTEND=noninteractive \
        apt-get install -y --no-install-recommends "${PACKAGES[@]}"
    check_dependencies
    ((${#MISSING[@]} == 0)) || fail '安装后仍有不可用的依赖；请检查上面的编译/运行诊断'
fi

python3 - "$SCRIPT_DIR" <<'PY'
import hashlib
import json
from pathlib import Path
import sys

dependencies = Path(sys.argv[1])
metadata = json.loads((dependencies / "sources.json").read_text())
archive = (dependencies.parent / metadata["hypre"]["archive"]).resolve()
if not archive.is_relative_to(dependencies) or not archive.is_file():
    sys.exit("错误：随附 Hypre 源码包缺失或路径无效")
if hashlib.sha256(archive.read_bytes()).hexdigest() != metadata["hypre"]["sha256"]:
    sys.exit("错误：随附 Hypre 源码包 SHA-256 校验失败")
print("[通过] 本地 PETSc/Eigen 源码入口及 Hypre 源码包校验")
PY

if ((CHECK_ONLY || DEPS_ONLY)); then
    printf '\n依赖检查完成。\n'
    exit 0
fi

# Pin this entry point to the delivery sources even if PETSC_DIR is exported.
MAKE_ARGS=(-C "$PROJECT_DIR" "BUILD=build" "PETSC_DIR=$PROJECT_DIR/build/petsc"
           "PETSC_SOURCE=$SCRIPT_DIR/petsc" "EIGEN_DIR=$SCRIPT_DIR/eigen"
           "PETSC_ARCH=arch-babelsim-build" "PETSC_JOBS=$JOBS"
           "CXX=mpic++" "PETSC_CC=mpicc" "PETSC_FC=mpifort")
if ((CLEAN_BUILD)); then
    make "${MAKE_ARGS[@]}" clean
fi
printf '\n从本地源码编译（并行数：%s）：\n' "$JOBS"
make "${MAKE_ARGS[@]}" "-j$JOBS" all
printf '\n构建完成：\n  %s\n  %s\n  %s\n  %s\n' \
    "$PROJECT_DIR/build/petsc/lib" "$PROJECT_DIR/build/libbabelsim.a" \
    "$PROJECT_DIR/build/babelsim-solve" "$PROJECT_DIR/build/babelsim-post"
