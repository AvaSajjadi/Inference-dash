#!/usr/bin/env python3
"""Find Rscript location, checking multiple methods."""

import os
import subprocess
import shutil
from pathlib import Path

def find_rscript():
    """Find Rscript by trying multiple methods."""

    # Method 1: Check environment variable
    rscript = os.environ.get("RSCRIPT_PATH")
    if rscript and Path(rscript).exists():
        return rscript

    # Method 2: Check PATH with shutil.which
    rscript = shutil.which("Rscript")
    if rscript:
        return rscript

    # Method 3: Try executing R/Rscript and get the path
    try:
        result = subprocess.run(["bash", "-c", "which Rscript"],
                              capture_output=True, text=True, timeout=5)
        if result.returncode == 0:
            rscript = result.stdout.strip()
            if rscript and Path(rscript).exists():
                return rscript
    except:
        pass

    # Method 4: Common hardcoded paths
    common_paths = [
        "/root/.nix-profile/bin/Rscript",
        "/usr/bin/Rscript",
        "/usr/local/bin/Rscript",
        "/nix/var/nix/profiles/default/bin/Rscript",
        Path("/nix/store").glob("*/bin/Rscript"),  # Search Nix store
    ]

    for path in common_paths:
        if isinstance(path, Path):
            # Handle glob results
            for match in path:
                if match.exists():
                    return str(match)
        else:
            if Path(path).exists():
                return path

    # Method 5: Last resort - search filesystem
    try:
        result = subprocess.run(
            ["find", "/nix/store", "-maxdepth", "3", "-name", "Rscript", "-type", "f"],
            capture_output=True, text=True, timeout=10
        )
        if result.returncode == 0 and result.stdout:
            rscript = result.stdout.split('\n')[0]
            if rscript and Path(rscript).exists():
                return rscript
    except:
        pass

    return None

if __name__ == "__main__":
    rscript = find_rscript()
    if rscript:
        print(f"Found: {rscript}")
        exit(0)
    else:
        print("NOT FOUND")
        exit(1)
