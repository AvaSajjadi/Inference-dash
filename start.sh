#!/bin/bash
# Startup wrapper that sources Nix environment and starts the app

# Source Nix environment to get R in PATH
if [ -f /nix/var/nix/profiles/default/etc/profile.d/nix.sh ]; then
    source /nix/var/nix/profiles/default/etc/profile.d/nix.sh
fi

# Add common Nix paths to PATH
export PATH="/nix/var/nix/profiles/default/bin:$PATH"
export PATH="/nix/store/*/bin:$PATH"

# Verify R is available
if ! command -v Rscript &> /dev/null; then
    echo "WARNING: Rscript still not in PATH, searching for it..."
    RSCRIPT=$(find /nix -name Rscript -type f 2>/dev/null | head -1)
    if [ -n "$RSCRIPT" ]; then
        echo "Found Rscript at: $RSCRIPT"
        export PATH="$(dirname $RSCRIPT):$PATH"
    fi
fi

echo "PATH=$PATH"
echo "Rscript location: $(which Rscript)"

# Start the app
python app2.py
