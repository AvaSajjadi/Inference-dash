#!/bin/bash
# Find Rscript in the Nix store and save its path
RSCRIPT=$(which Rscript)
if [ -z "$RSCRIPT" ]; then
    # Try to find it in common Nix locations
    RSCRIPT=$(find /nix -name Rscript -type f 2>/dev/null | head -1)
fi

if [ -n "$RSCRIPT" ]; then
    echo "$RSCRIPT" > /app/.rscript_path
    chmod 644 /app/.rscript_path
    echo "Rscript found at: $RSCRIPT"
else
    echo "ERROR: Rscript not found" >&2
    exit 1
fi
