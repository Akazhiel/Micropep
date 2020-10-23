# We will use Ubuntu for our image

FROM ubuntu:18.04

WORKDIR /Micropeptide

# Updating Ubuntu packages

RUN apt-get update && apt-get install -y \
    wget \
    bzip2

# Anaconda installing

RUN wget https://repo.anaconda.com/archive/Anaconda3-2020.07-Linux-x86_64.sh \
    && bash Anaconda3-2020.07-Linux-x86_64.sh -b \
    && rm Anaconda3-2020.07-Linux-x86_64.sh

# Set path to conda

ENV PATH /root/anaconda3/bin:$PATH

# Copy conda environment yml

COPY * .

# Updating Anaconda packages

RUN conda config --add channels conda-forge \
    && conda config --add channels bioconda
RUN conda env create -f environment.yml

# Pull the environment name out of the environment.yml
RUN echo "source activate $(head -1 environment.yml | cut -d' ' -f2)" > ~/.bashrc
ENV PATH /opt/conda/envs/$(head -1 environment.yml | cut -d' ' -f2)/bin:$PATH
RUN rm environment.yml

# Copying pipeline scripts
RUN python setup.py install