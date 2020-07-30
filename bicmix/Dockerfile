
# Base image
FROM ubuntu:20.04

# Maintainer and author
LABEL maintainer="Bryan Quach <bquach@rti.org>"

#### Basic Ubuntu setup ####

# Set a default user. Available via runtime flag `--user docker` 
# Add user to 'staff' group and make home directory 
RUN useradd docker \
    && mkdir /home/docker \
    && chown docker:docker /home/docker \
    && addgroup docker staff

ENV DEBIAN_FRONTEND noninteractive
    
RUN echo 'deb http://mirror.math.princeton.edu/pub/ubuntu/ focal main' >> /etc/apt/sources.list \
    && apt-get update \ 
    && apt-get install -y --no-install-recommends \
        dirmngr \
        software-properties-common \
        lsb-release \
        ed \
        less \
        locales \
        vim \
        wget \
        ca-certificates \
        fonts-texgyre \
        make \
        libglib2.0-dev \
    && rm -rf /var/lib/apt/lists/*

# Configure default locale, see https://github.com/rocker-org/rocker/issues/19
RUN echo "en_US.UTF-8 UTF-8" >> /etc/locale.gen \
    && locale-gen en_US.utf8 \
    && /usr/sbin/update-locale LANG=en_US.UTF-8

ENV LC_ALL en_US.UTF-8
ENV LANG en_US.UTF-8

#### Install dependencies ####

RUN apt-get update \ 
    && apt-get install -y --no-install-recommends \
        unzip \
        gcc \
        g++ \
        gsl-bin \
        libgsl-dev \
        libboost-dev \
        libeigen3-dev

WORKDIR /

#### Install bicmix ####
RUN wget https://www.cs.princeton.edu/~bee/software/BicMix-Code-for-distribution.zip \
    && unzip BicMix-Code-for-distribution.zip \
    && rm BicMix-Code-for-distribution.zip \
    && mv BicMix-Code-for-distribution /opt/ \
    && cd /opt/BicMix-Code-for-distribution/ \
    && perl -i -pe 's/Your_library/usr\/include/g' Makefile \
    && mv /usr/include/eigen3/Eigen/ /usr/include/ \
    && make \
    && mv /opt/BicMix-Code-for-distribution/BicMix /opt/ \
    && cd / \
    && rm -rf /opt/BicMix-Code-for-distribution
    
WORKDIR /

# Add /opt to executables path list
ENV PATH=pkg-config:$PATH:/opt/

# Default command 
CMD ["/opt/BicMix"]
