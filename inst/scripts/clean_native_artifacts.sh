#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  clean_native_artifacts.sh [--delete] [PACKAGE_DIR]

Clears stale compiled native artifacts that can break SIRUP builds when switching
between macOS and Linux/Discovery builds.

By default this is a dry run. Pass --delete to remove matching files.

Examples:
  SIRUP/inst/scripts/clean_native_artifacts.sh
  SIRUP/inst/scripts/clean_native_artifacts.sh --delete
  SIRUP/inst/scripts/clean_native_artifacts.sh --delete SIRUP
EOF
}

delete=false
package_dir=""

while [ "$#" -gt 0 ]; do
  case "$1" in
    --delete)
      delete=true
      shift
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    -*)
      echo "Unknown option: $1" >&2
      usage >&2
      exit 2
      ;;
    *)
      package_dir="$1"
      shift
      ;;
  esac
done

if [ -z "$package_dir" ]; then
  if [ -d "SIRUP/src" ]; then
    package_dir="SIRUP"
  elif [ -d "src" ] && [ -f "DESCRIPTION" ]; then
    package_dir="."
  else
    package_dir="SIRUP"
  fi
fi

if [ ! -d "$package_dir" ]; then
  echo "Package directory not found: $package_dir" >&2
  exit 1
fi

src_dir="${package_dir%/}/src"
if [ ! -d "$src_dir" ]; then
  echo "No src directory found at: $src_dir"
  exit 0
fi

find_args=(
  "$src_dir"
  -maxdepth 1
  -type f
  "("
    -name "*.o"
    -o -name "*.so"
    -o -name "*.dll"
    -o -name "*.dylib"
    -o -name "*.o.tmp"
  ")"
)

matches=$(find "${find_args[@]}" -print)
if [ -z "$matches" ]; then
  echo "No stale native artifacts found."
  exit 0
fi

printf '%s\n' "$matches"
count=$(printf '%s\n' "$matches" | wc -l | tr -d ' ')

if [ "$delete" = true ]; then
  find "${find_args[@]}" -delete
  echo "Removed $count stale native artifact(s)."
else
  echo
  echo "Dry run: found $count stale native artifact(s)."
  echo "Run again with --delete to remove them."
fi
