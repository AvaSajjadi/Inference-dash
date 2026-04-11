#!/bin/bash
# Startup wrapper that sources Nix environment and starts the app

echo "=== Starting Nix environment setup ==="

# Check if R is in PATH already
if command -v R &> /dev/null; then
    echo "✓ R found in PATH: $(which R)"
fi

if command -v Rscript &> /dev/null; then
    echo "✓ Rscript found in PATH: $(which Rscript)"
    export RSCRIPT_PATH="$(which Rscript)"
fi

# Try sourcing Nix profile
if [ -f /nix/var/nix/profiles/default/etc/profile.d/nix.sh ]; then
    echo "Sourcing Nix profile..."
    source /nix/var/nix/profiles/default/etc/profile.d/nix.sh
fi

# Add common Nix paths to PATH
export PATH="/nix/var/nix/profiles/default/bin:/nix/store/*/bin:$PATH"

# If still not found, search Nix store
if [ -z "$RSCRIPT_PATH" ]; then
    echo "Searching /nix/store for Rscript..."
    RSCRIPT=$(find /nix/store -maxdepth 3 -name Rscript -type f 2>/dev/null | head -1)
    if [ -n "$RSCRIPT" ]; then
        echo "Found: $RSCRIPT"
        export RSCRIPT_PATH="$RSCRIPT"
    fi
fi

# Last resort: check which after all updates
if [ -z "$RSCRIPT_PATH" ]; then
    echo "Checking which Rscript..."
    RSCRIPT=$(which Rscript 2>/dev/null)
    if [ -n "$RSCRIPT" ]; then
        echo "Found: $RSCRIPT"
        export RSCRIPT_PATH="$RSCRIPT"
    fi
fi

echo "=== Environment Setup Complete ==="
echo "RSCRIPT_PATH=${RSCRIPT_PATH:-NOT FOUND}"
echo "R location: $(which R 2>/dev/null || echo 'NOT FOUND')"
echo "Rscript location: $(which Rscript 2>/dev/null || echo 'NOT FOUND')"

if [ -z "$RSCRIPT_PATH" ]; then
    echo "ERROR: Rscript still not found"
    ls -la /nix/var/nix/profiles/default/bin/ 2>/dev/null || echo "Cannot list nix bin dir"
    exit 1
fi

# Start the app
python app2.py
