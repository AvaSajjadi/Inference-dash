FROM rocker/r-base:latest

# Install Python and system dependencies
RUN apt-get update && apt-get install -y \
    python3 python3-pip \
    curl git \
    libxml2-dev \
    && rm -rf /var/lib/apt/lists/*

# Set working directory
WORKDIR /app

# Copy app files
COPY . /app

# Install Python requirements
RUN pip install --no-cache-dir -r requirements.txt

# Install R packages
RUN Rscript install.R

# Expose port
EXPOSE 8080

# Set environment
ENV PORT=8080

# Start the app
CMD ["python3", "app2.py"]
