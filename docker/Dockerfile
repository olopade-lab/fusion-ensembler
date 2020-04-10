FROM ubuntu:16.04
MAINTAINER annawoodard@uchicago.edu
# adapted from https://github.com/trinityrnaseq/trinityrnaseq/blob/master/Docker/Dockerfile

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y gcc g++ perl python automake make \
       wget git curl libdb-dev \
       zlib1g-dev bzip2 libncurses5-dev \
       texlive-latex-base \
       default-jre \
       python-pip python-dev \
       vim \
       gfortran \
       build-essential libghc-zlib-dev libncurses-dev libbz2-dev liblzma-dev libpcre3-dev libxml2-dev \
       libblas-dev gfortran git unzip ftp libzmq3-dev ftp fort77 libreadline-dev \
       libcurl4-openssl-dev libx11-dev libxt-dev \
       x11-common libcairo2-dev libpng-dev libreadline-dev libjpeg-dev pkg-config libtbb-dev \
       cmake rsync libssl-dev tzdata \
       && apt-get clean



## perl stuff
RUN curl -L https://cpanmin.us | perl - App::cpanminus
RUN cpanm install DB_File
RUN cpanm install URI::Escape
RUN cpanm install PerlIO::gzip


## set up tool config and deployment area
ENV SRC /usr/local/src
ENV BIN /usr/local/bin


## arriba
WORKDIR $SRC
RUN wget -q -O - "https://github.com/suhrig/arriba/releases/download/v1.2.0/arriba_v1.2.0.tar.gz" | tar -xzf -
-O bowtie-1.2.1.1-linux-x86_64.zip && \
        unzip bowtie-1.2.1.1-linux-x86_64.zip && \
	mv bowtie-1.2.1.1/bowtie* $BIN


## bowtie
WORKDIR $SRC
RUN wget https://sourceforge.net/projects/bowtie-bio/files/bowtie/1.2.1.1/bowtie-1.2.1.1-linux-x86_64.zip/download -O bowtie-1.2.1.1-linux-x86_64.zip && \
        unzip bowtie-1.2.1.1-linux-x86_64.zip && \
	mv bowtie-1.2.1.1/bowtie* $BIN


# blast
WORKDIR $SRC
RUN wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.5.0/ncbi-blast-2.5.0+-x64-linux.tar.gz && \
    tar xvf ncbi-blast-2.5.0+-x64-linux.tar.gz && \
    cp ncbi-blast-2.5.0+/bin/* $BIN && \
    rm -r ncbi-blast-2.5.0+


## bowtie2
WORKDIR $SRC
RUN wget https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.3.4.1/bowtie2-2.3.4.1-linux-x86_64.zip/download -O bowtie2-2.3.4.1-linux-x86_64.zip && \
    unzip bowtie2-2.3.4.1-linux-x86_64.zip && \
    mv bowtie2-2.3.4.1-linux-x86_64/bowtie2* $BIN && \
    rm *.zip && \
    rm -r bowtie2-2.3.4.1-linux-x86_64



## samtools
RUN wget https://github.com/samtools/samtools/releases/download/1.7/samtools-1.7.tar.bz2 && \
    tar xvf samtools-1.7.tar.bz2 && \
    cd samtools-1.7/ && \
    ./configure && make && make install


# ## jellyfish
# RUN wget https://github.com/gmarcais/Jellyfish/releases/download/v2.3.0/jellyfish-2.3.0.tar.gz && \
#     tar xvf jellyfish-2.3.0.tar.gz && \
#     cd jellyfish-2.3.0/ && \
#     ./configure && make && make install


## RSEM
RUN mkdir /usr/local/lib/site_perl
WORKDIR $SRC
RUN wget https://github.com/deweylab/RSEM/archive/v1.3.3.tar.gz && \
	 tar xvf v1.3.3.tar.gz && \
     cd RSEM-1.3.3 && \
     make && \
     cp rsem-* $BIN && \
     cp convert-sam-for-rsem $BIN && \
     cp rsem_perl_utils.pm /usr/local/lib/site_perl/ && \
     cd ../ && rm -r RSEM-1.3.3


## salmon
WORKDIR $SRC
ENV SALMON_VERSION=1.0.0
RUN wget https://github.com/COMBINE-lab/salmon/releases/download/v${SALMON_VERSION}/Salmon-${SALMON_VERSION}_linux_x86_64.tar.gz && \
    tar xvf Salmon-${SALMON_VERSION}_linux_x86_64.tar.gz && \
    ln -s $SRC/Salmon-latest_linux_x86_64/bin/salmon $BIN/.

WORKDIR $SRC
RUN ln -sf /usr/local/src/salmon-latest_linux_x86_64/bin/salmon $BIN/.


## STAR
ENV STAR_VERSION=2.7.3a
RUN STAR_URL="https://github.com/alexdobin/STAR/archive/${STAR_VERSION}.tar.gz" &&\
    wget -P $SRC $STAR_URL && \
    tar -xvf $SRC/${STAR_VERSION}.tar.gz -C $SRC && \
    mv $SRC/STAR-${STAR_VERSION}/bin/Linux_x86_64_static/STAR /usr/local/bin


## blat
RUN wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/blat/blat -P $BIN && \
    chmod a+x $BIN/blat


## Trinity
ENV TRINITY_VERSION="2.9.1"
ENV TRINITY_CO="e806a0f341e4d3917c2b181c4ecd4a4c3f4367fc"

WORKDIR $SRC
RUN git clone --recursive https://github.com/trinityrnaseq/trinityrnaseq.git && \
    cd trinityrnaseq && \
    git checkout ${TRINITY_CO} && \
    git submodule init && git submodule update && \
    rm -rf ./trinity_ext_sample_data && \
    make && make plugins && \
    make install && \
    cd ../ && rm -r trinityrnaseq

ENV TRINITY_HOME /usr/local/bin/trinityrnaseq
ENV PATH=${TRINITY_HOME}:${PATH}


# COPY run.sh /usr/local/bin
# RUN chmod a+x /usr/local/bin/run.sh

# cleanup
WORKDIR $SRC
RUN rm -r *.tar.gz *.bz2
RUN apt-get clean
