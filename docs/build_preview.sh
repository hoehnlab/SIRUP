#!/bin/bash
# Build a static copy of the docs site for the lab, until the package is public.
#
#   bash docs/build_preview.sh
#
# The site goes to ../docs_preview (next to the package folder, not inside it):
#   ../docs_preview/sirup-docs/      the site as plain HTML files
#   ../docs_preview/sirup-docs.zip   the same site in one zip file
#   ../docs_preview/README.txt       how to open it, and which commit it came from
# Open sirup-docs/index.html in a browser. The pages work without a web server.
#
# build_docs.sbatch runs this script after each full docs build. To use another
# folder, set PREVIEW_DIR. To use another Python with mkdocs, set MKDOCS_PYTHON.
set -euo pipefail

cd "$(dirname "$0")/.."
PREVIEW_DIR="${PREVIEW_DIR:-$(cd .. && pwd)/docs_preview}"
MKDOCS_PYTHON="${MKDOCS_PYTHON:-python3}"

mkdir -p "$PREVIEW_DIR"
work="$(mktemp -d)"
trap 'rm -rf "$work"' EXIT

# --no-directory-urls writes links to page.html files, so they work from a folder.
"$MKDOCS_PYTHON" -m mkdocs build --no-directory-urls -q -d "$work/sirup-docs"

rsync -rlt --delete "$work/sirup-docs/" "$PREVIEW_DIR/sirup-docs/"
(cd "$work" && zip -qr sirup-docs.zip sirup-docs)
mv "$work/sirup-docs.zip" "$PREVIEW_DIR/sirup-docs.zip"

if [ -n "$(git status --porcelain 2>/dev/null)" ]; then
  state="with changes that are not committed"
else
  state="from a clean working tree"
fi
cat > "$PREVIEW_DIR/README.txt" <<README
SIRUP documentation preview

Open sirup-docs/index.html in a web browser. Or unzip sirup-docs.zip and open
index.html in the unzipped folder.

Built:  $(date '+%Y-%m-%d %H:%M')
Commit: $(git rev-parse --short HEAD 2>/dev/null || echo unknown), $state

This is a snapshot for the lab. It changes each time the docs are rebuilt.
README
echo "Preview written to $PREVIEW_DIR"
