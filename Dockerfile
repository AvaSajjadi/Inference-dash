FROM rocker/r-base:latest

# Install Python and system dependencies
RUN apt-get update && apt-get install -y \
    python3 python3-pip \
    curl git \
    libxml2-dev \
    && rm -rf /var/lib/apt/lists/*

# Verify R installation
RUN which R && R --version && which Rscript && Rscript --version

# Set working directory
WORKDIR /app

# Copy app files
COPY . /app

# Install Python requirements
RUN pip install --no-cache-dir -r requirements.txt

# Install R packages
RUN Rscript install.R

# Verify packages are installed
RUN Rscript -e "library(dplyr); library(magrittr); library(data.table); cat('All packages loaded successfully\n')"

# Expose port
EXPOSE 8080

# Set environment
ENV PORT=8080

# Test that Rscript works
RUN echo 'cat("Rscript test successful\n")' | Rscript -

# Start the app
CMD ["python3", "app2.py"]
