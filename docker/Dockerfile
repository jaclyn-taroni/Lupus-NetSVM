FROM ubuntu:15.10
MAINTAINER "Jaclyn Taroni" <jtaroni@mail.med.upenn.edu>

# Adapted from Greg Way (@gwaygenomics)
# https://github.com/greenelab/hgsc_subtypes.git Sep 20 2016

# Install Ubuntu dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
		git \
                ed \		
		less \
		locales \
		vim-tiny \
		wget \
                ca-certificates \
		bison \
		flex \
		apt-utils \
		gawk \
		libcurl4-openssl-dev \
		libxml2-dev\
		libssl-dev

# Install R
RUN apt-get install -y --no-install-recommends \
		r-base \
		r-base-dev \
		r-recommended \
		r-cran-xml \
		r-cran-rjava

# Install Sleipnir
RUN apt-get install -y --no-install-recommends \
		mercurial \
		gengetopt \
		libboost-regex-dev \
		libboost-graph-dev \
		liblog4cpp5-dev \
		build-essential \
		libgsl0-dev

RUN wget -nv http://download.joachims.org/svm_perf/current/svm_perf.tar.gz \
--directory-prefix svmperf/
RUN cd svmperf && tar -xf svm_perf.tar.gz && make

RUN hg clone https://bitbucket.org/libsleipnir/sleipnir

RUN cd sleipnir && ./gen_auto && ./gen_tools_am
RUN cd sleipnir && ./configure --with-svm-perf=~/svmperf/ && make && make install

