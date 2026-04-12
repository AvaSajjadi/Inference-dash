FROM rocker/r-base:latest

# Accept optional GitHub token for CIE installation
ARG GITHUB_TOKEN=""

# Install Python and system dependencies
RUN apt-get update && apt-get install -y \
    python3 python3-pip \
    curl git \
    libxml2-dev \
    build-essential \
    && rm -rf /var/lib/apt/lists/*

# Verify R installation
RUN which R && R --version && which Rscript && Rscript --version

# Set working directory
WORKDIR /app

# Copy app files
COPY . /app

# Copy all locally installed R packages (includes CIE and all dependencies)
COPY --chown=root:root R-packages/ /usr/local/lib/R/site-library/

# Install Python requirements (use --break-system-packages for Python 3.13+)
RUN pip install --no-cache-dir --break-system-packages -r requirements.txt

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

# Create uploads directory
RUN mkdir -p /app/uploads && chmod 777 /app/uploads

# Test that Rscript works
RUN echo 'cat("Rscript test successful\n")' | Rscript -

# Copy entrypoint script
COPY docker-entrypoint.sh /app/docker-entrypoint.sh
RUN chmod +x /app/docker-entrypoint.sh

# Start the app via entrypoint (attempts CIE install at runtime)
ENTRYPOINT ["/app/docker-entrypoint.sh"]
