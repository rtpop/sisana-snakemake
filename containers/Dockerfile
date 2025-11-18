FROM python:3.10.16

## Author: Romana Pop
## Date November 2025


## ------------ ##
## How to build ##
## ------------ ##

# To build the Docker container image, run the following command
# this can and should only be done by the owner (i.e., rtpop)
# sudo docker build --file containers/Dockerfile --no-cache -t rtpop/sisana_container:tag --progress=plain . 2>&1 | tee containers/sisana_container/logs/sisana_container:tag

## --------------------------- ##
## Install system requirements ##
## --------------------------- ##

RUN mkdir -p /opt/software
WORKDIR /opt

RUN apt-get update -qq && \
        apt-get install -y --allow-downgrades \
                apt-utils

RUN apt-cache policy \
        apt-utils \
        build-essential \
        freetype2-demos \
        freetype2-doc \
        git \
        libcurl4-openssl-dev \
        libcairo2-dev \
        libfontconfig1-dev \
        libfreetype-dev \
        libfreetype6 \
        librsvg2-bin \
        librsvg2-dev \
        libudunits2-dev \
        make \
        pkg-config \
        procps \
        libssl-dev \
        libopenmpi-dev \
        libgsl-dev \
        libcurl4-gnutls-dev \
        libmysqlclient-dev \
        libpq-dev

RUN apt-get update -qq && \
        apt-get install -y --allow-downgrades \
                build-essential \
                cmake \
                git \
                libcurl4-openssl-dev \
                libcairo2-dev \
                libssl-dev \
                vim

## ------------ ##
## Pip installs ##
## ------------ ##
# install SiSaNa
RUN pip install sisana