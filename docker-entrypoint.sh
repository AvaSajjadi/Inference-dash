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
exec python3 app2.py
