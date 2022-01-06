# syntax=docker/dockerfile:experimental
FROM ubuntu:20.04

RUN apt update

# allow new apt repositories
RUN apt install -y \
	apt-transport-https \
	software-properties-common

RUN apt install -y \
	git \
	libcurl4-openssl-dev \
	libssl-dev \
	libxml2-dev \
    python3=3.8* \
    python3-pip \
    r-base

# install cellforest
RUN pip3 install
		dataforest \
		cellforest

# install cellforestR
RUN	(git clone --branch hotfix/cellforestR-dependencies https://github.com/TheAustinator/cellforest.git && \
	 cd cellforest && \
	 R -e 'install.packages("devtools")' && \
	 R -e 'library("devtools"); devtools::install("cellforestR", dependencies = TRUE)'\
	)