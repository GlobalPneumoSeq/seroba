# This container will install SeroBA from master

# base image
FROM ubuntu:bionic as app

# set workdir to default for building; set to /data at the end
WORKDIR /

# Install ubuntu dependencies
RUN apt-get update && apt-get -y upgrade && apt-get -y install git \
  wget \
  unzip \
  zlib1g-dev \
  libbz2-dev \
  libjpeg-dev \
  liblzma-dev \
  python3 \
  python3-pip \
  libpython3-dev \
  python3-setuptools \
  python-minimal && \
  apt-get clean && apt-get autoclean && rm -rf /var/lib/apt/lists/*

# Copy repository into the image and install dependencies
COPY . /seroba/
RUN cd seroba && \
  pip3 install pysam==0.15.0 &&\
  /seroba/install_dependencies.sh

# set path
ENV PATH="/seroba:/seroba/build:/seroba/build/MUMmer3.23:/seroba/build/bowtie2-2.3.1-legacy:/seroba/build/cdhit-4.6.8:${PATH}"

# install seroba and create database (/seroba/database/)
RUN cd /seroba && \
  python3 setup.py install && \
  seroba createDBs database/ 71

RUN mkdir /data
WORKDIR /data

# new base for testing
FROM app as test

# print out various help options and version
RUN seroba version && \
  seroba --help

# run built-in test
RUN cd /seroba && \
  python3 setup.py test