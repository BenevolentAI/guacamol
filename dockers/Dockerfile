FROM ubuntu:18.04

# metainformation
LABEL org.opencontainers.image.version = "0.5.5"
LABEL org.opencontainers.image.authors = "BenevolentAI"
LABEL org.opencontainers.image.source = "https://github.com/BenevolentAI/guacamol"
LABEL org.opencontainers.image.licenses = "MIT"
LABEL org.opencontainers.image.base.name="docker.io/library/ubuntu:18.04"

RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential \
        cmake ca-certificates \
        libglib2.0-0 libxext6 libsm6 libxrender1 \
        wget curl bash bzip2 && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# MiniConda
RUN curl -LO --silent https://repo.continuum.io/miniconda/Miniconda3-4.5.11-Linux-x86_64.sh && \
    bash Miniconda3-4.5.11-Linux-x86_64.sh -p /miniconda -b && \
    rm Miniconda3-4.5.11-Linux-x86_64.sh

ENV PATH=/miniconda/bin:${PATH}

# Add the source code
RUN mkdir -p /app
ADD . /app

# python deps for running tests
RUN pip install --upgrade pip && pip install --no-cache-dir -r /app/dockers/requirements.txt

# install guacamol
RUN pip install --upgrade pip && pip install --no-cache-dir /app/

# Launch inside the folder
WORKDIR /app/
