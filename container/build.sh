#!/bin/bash
set -euo pipefail

# ── Build the Psomagen Bulk RNA-Seq container ──
# Reads version from VERSION file, produces psomagen_v<VERSION>.sif

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
VERSION="$(cat "${SCRIPT_DIR}/VERSION" | tr -d '[:space:]')"
SIF_NAME="psomagen_v${VERSION}.sif"
DEF_FILE="${SCRIPT_DIR}/psomagen.def"

# Safety: refuse to overwrite unless --force
if [[ -f "${SIF_NAME}" && "${1:-}" != "--force" ]]; then
    echo "ERROR: ${SIF_NAME} already exists. Use --force to overwrite."
    exit 1
fi

echo "Building ${SIF_NAME} from ${DEF_FILE} ..."
apptainer build --fakeroot "${SIF_NAME}" "${DEF_FILE}"
echo "Done: ${SIF_NAME}"
