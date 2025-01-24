FROM debian:latest

# Dockerfile for VirMake

# Set environment variables
ENV DEBIAN_FRONTEND="noninteractive"
ENV MACHINE=x86_64

# Update Debian and install required packages
RUN apt-get update && apt-get -y upgrade && apt-get -y install wget time

# Install Miniconda
WORKDIR /conda
RUN wget -nv "https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-$MACHINE.sh" -O conda.sh && \
    bash /conda/conda.sh -b -u -p /conda && \
    rm /conda/conda.sh

# Copy VirMake repo files
WORKDIR /VirMake
ADD LICENSE README.md entrypoint.sh setup.py .
ADD resources resources
ADD utils utils
ADD workflow workflow

# Setup
RUN bash -c "source /conda/bin/activate; python setup.py"

# Enter here
ENTRYPOINT ["/VirMake/entrypoint.sh"]
