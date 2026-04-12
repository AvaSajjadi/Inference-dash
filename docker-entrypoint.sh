#!/bin/bash
set -e

echo "Starting Inference Dash..."

# Verify CIE is available
echo "Verifying CIE package..."
Rscript -e "
if (require('CIE', character.only = TRUE, quietly = TRUE)) {
    cat('CIE is available\n')
} else {
    cat('Warning: CIE package not found. CIE analysis will not be available.\n')
    cat('ORNOR analysis will still work.\n')
}
" 2>&1

# Start the application
echo "Starting Flask application..."
echo "Current directory: $(pwd)"
echo "Python version: $(python3 --version)"
echo "Testing imports..."
python3 << 'PYEOF' || {
    echo "❌ Failed to import app2"
    exit 1
}
try:
    print("Importing app2...")
    from app2 import app
    print("✅ app2 imported successfully")
except Exception as e:
    print(f"❌ Error: {e}")
    import traceback
    traceback.print_exc()
    raise
PYEOF

echo "Launching Flask..."
exec python3 app2.py
