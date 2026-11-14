#!/usr/bin/env bash
# Cross-platform GUI launcher for macOS/Linux.
set -euo pipefail
cd -- "$(dirname -- "$0")"
exec python3 dna_sequence_analysis.py
