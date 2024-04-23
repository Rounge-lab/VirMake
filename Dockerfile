FROM debian:latest

# Dockerfile for VirMake

# Update and install packages
RUN apt update && apt -y upgrade && apt -y install time wget

# Install Miniconda
WORKDIR /opt/conda
RUN wget "https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-$(uname -m).sh" -O miniconda.sh && \
    bash /opt/conda/miniconda.sh -b -u -p /opt/conda && \
    rm /opt/conda/miniconda.sh

# Copy VirMake repo
WORKDIR /VirMake
COPY . .

# Enter here
ENTRYPOINT ["/bin/bash"]
