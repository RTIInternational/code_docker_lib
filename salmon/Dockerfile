
# Base image
FROM ubuntu:18.04

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
    
RUN echo 'deb http://mirror.math.princeton.edu/pub/ubuntu/ bionic main' >> /etc/apt/sources.list \
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
        libglib2.0 \
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
        gcc \
        g++ \
        autoconf \
        zlib1g-dev \
        libboost1.65-all-dev \
        libc6-dev \
    && rm -rf /var/lib/apt/lists/* \
    && apt-get autoremove

### Install MashMap ###
RUN wget https://github.com/marbl/MashMap/archive/v2.0.tar.gz \
    && tar -xzf v2.0.tar.gz \
    && rm v2.0.tar.gz \
    && cd MashMap-2.0 \
    && mkdir /opt/lib \
    && cp /usr/lib/x86_64-linux-gnu/libboost_math_c99.a /opt/lib/ \
    && ./bootstrap.sh \
    && ./configure --prefix=/opt/ --with-boost=/opt \
    && make \
    && make install \
    && cd / \
    && rm -rf MashMap-2.0
    
#### Install Salmon ####

RUN wget https://github.com/COMBINE-lab/salmon/releases/download/v0.14.1/salmon-0.14.1_linux_x86_64.tar.gz \
    && tar -xzf salmon-0.14.1_linux_x86_64.tar.gz \
    && rm salmon-0.14.1_linux_x86_64.tar.gz \
    && mv salmon-latest_linux_x86_64/bin/* /opt/bin/ \
    && mv salmon-latest_linux_x86_64/lib/* /opt/lib/ \
    && rm -r salmon-latest_linux_x86_64/

WORKDIR /

ENV PATH=pkg-config:$PATH:/opt/bin/
ENV LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/lib/

CMD ["salmon", "--help"]