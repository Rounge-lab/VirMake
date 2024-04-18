FROM debian:latest

# Dockerfile for VirMake

# Update and install packages
RUN apt update && apt upgrade && apt -y install time wget

# Install Miniconda
ARG TARGETARCH
ARG CONDA_VER="latest"
ARG CONDA_URL="https://repo.anaconda.com/miniconda/Miniconda3-$CONDA_VER"
WORKDIR /opt/conda
RUN if [ "$TARGETARCH" = "amd64" ]; then \
      URL="$CONDA_URL-Linux-x86_64.sh"; \
    elif [ "$TARGETARCH" = "arm64" ]; then \
      URL="$CONDA_URL-Linux-aarch64.sh"; \
    else \
      echo "Unsupported target architecture"; \
      exit 1; \
    fi && \
    wget "$URL" -O miniconda.sh && \
    bash /opt/conda/miniconda.sh -b -u -p /opt/conda && \
    rm /opt/conda/miniconda.sh

# Copy VirMake repo
WORKDIR /VirMake
COPY . .

# Enter here
ENTRYPOINT ["/bin/bash"]
