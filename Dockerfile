FROM staphb/ariba:2.14.6 as ariba_stage

FROM ubuntu:bionic as app

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

# Copy ARIBA installation from ariba_stage
COPY --from=ariba_stage /usr/local/ /usr/local/
COPY --from=ariba_stage /usr/bin/ /usr/bin/

# Copy repository into the image and install dependencies
COPY . /seroba/
RUN cd /seroba && \
  /seroba/install_dependencies.sh

# set path
ENV PATH="/seroba:/seroba/build:/seroba/build/bin:/seroba/build/MUMmer3.23:/seroba/build/bowtie2-2.3.1-legacy:/seroba/build/cdhit-4.6.8:${PATH}"

# install seroba and create database
RUN cd /seroba && \
  python3 setup.py install && \
  seroba createDBs database/ 71

RUN mkdir /data
WORKDIR /data

FROM app as test

RUN seroba version && \
  seroba --help

RUN cd /seroba && \
  python3 setup.py test