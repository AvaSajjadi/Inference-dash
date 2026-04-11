#!/bin/bash
# Startup wrapper that sources Nix environment and starts the app

# Source Nix environment to get R in PATH
if [ -f /nix/var/nix/profiles/default/etc/profile.d/nix.sh ]; then
    source /nix/var/nix/profiles/default/etc/profile.d/nix.sh
fi

# Add common Nix paths to PATH
export PATH="/nix/var/nix/profiles/default/bin:$PATH"

# Find Rscript and export its full path
RSCRIPT=$(which Rscript)
if [ -z "$RSCRIPT" ]; then
    RSCRIPT=$(find /nix -name Rscript -type f 2>/dev/null | head -1)
fi

if [ -n "$RSCRIPT" ]; then
    echo "Found Rscript at: $RSCRIPT"
    export RSCRIPT_PATH="$RSCRIPT"
    export PATH="$(dirname $RSCRIPT):$PATH"
else
    echo "ERROR: Rscript not found anywhere"
    exit 1
fi

echo "Exported RSCRIPT_PATH=$RSCRIPT_PATH"
echo "PATH=$PATH"

# Start the app with Rscript path in environment
python app2.py
