FROM rocker/r-base:latest

# Accept optional GitHub token for CIE installation
ARG GITHUB_TOKEN=""

# Install Python and system dependencies
# Note: curl is already provided by rocker/r-base; installing it separately
# causes libcurl4t64 version conflicts on newer Debian releases.
RUN apt-get update && apt-get install -y --no-install-recommends \
    git \
    libxml2-dev \
    libcurl4-gnutls-dev \
    libuv1-dev \
    build-essential \
    libgsl-dev \
    cython3 \
    && rm -rf /var/lib/apt/lists/*

# Install Python packages separately; python3 and pip come from the rocker base.
# python3-dev is needed to compile Cython extensions (nlbayes).
RUN apt-get update && apt-get install -y --no-install-recommends \
    python3 python3-pip python3-dev \
    && rm -rf /var/lib/apt/lists/* \
    || (apt-get install -y --no-install-recommends python3 python3-pip && rm -rf /var/lib/apt/lists/*)

# Verify R installation
RUN which R && R --version && which Rscript && Rscript --version

WORKDIR /app

# ── R SLOW ZONE ────────────────────────────────────────────────────────────────
# These layers are cached as long as R-packages/ doesn't change.
# Placing them before COPY . /app means Python-only edits never bust this cache.

# Copy all locally installed R packages (includes CIE and all dependencies)
COPY --chown=root:root R-packages/ /usr/local/lib/R/site-library/

# Force reinstall any packages built with an older R ABI.
# R-packages/ may contain binaries compiled for a different R version; checkBuilt=TRUE
# detects and reinstalls all such packages from CRAN automatically.
RUN Rscript -e "update.packages(ask=FALSE, checkBuilt=TRUE, repos='https://cloud.r-project.org')"

# ── PYTHON / APP FAST ZONE ─────────────────────────────────────────────────────
# BUILD_ID busts only the layers below — R compilation above stays cached.
ARG BUILD_ID="2026-05-03-09"
RUN echo "Build: $BUILD_ID"

# Copy app files
COPY . /app

# Install Python requirements
RUN pip install --no-cache-dir --break-system-packages -r requirements.txt 2>&1 | tail -20

# Setup nlbayes: copy source and install with proper flags
COPY nlbayes-python-src /app/nlbayes-src
RUN cd /app/nlbayes-src && \
    pip install --no-cache-dir --break-system-packages cython>=3.0 numpy scipy scikit-learn && \
    python3 setup.py build_ext --inplace && \
    pip install --no-cache-dir --break-system-packages . && \
    python3 -c "import sys; sys.path.insert(0, '/app'); from nlbayes import ModelORNOR; print('✅ nlbayes ready')" || echo "⚠️ nlbayes import issue (will retry at runtime)"

# Install R packages (pass GITHUB_TOKEN if provided)
RUN if [ -n "$GITHUB_TOKEN" ]; then \
        GITHUB_PAT=$GITHUB_TOKEN Rscript install.R; \
    else \
        Rscript install.R; \
    fi

# Verify packages are installed
RUN Rscript -e "library(dplyr); library(magrittr); library(data.table); cat('All packages loaded successfully\n')"

# Expose port
EXPOSE 8080

# Set environment
ENV PORT=8080
ENV PYTHONPATH="/app:${PYTHONPATH}"

# Create uploads directory
RUN mkdir -p /app/uploads && chmod 777 /app/uploads

# Test that Rscript works
RUN echo 'cat("Rscript test successful\n")' | Rscript -

# Copy entrypoint script
COPY docker-entrypoint.sh /app/docker-entrypoint.sh
RUN chmod +x /app/docker-entrypoint.sh

# Start the app via entrypoint (attempts CIE install at runtime)
ENTRYPOINT ["/app/docker-entrypoint.sh"]
