#!/bin/bash

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
" 2>&1 || true

# Start the application
echo "Starting Flask application..."
echo "Current directory: $(pwd)"
echo "Python version: $(python3 --version)"

# Launch Flask - let it handle any startup errors
exec python3 app2.py
