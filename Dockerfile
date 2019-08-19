FROM ubuntu:18.04 AS base

RUN apt-get clean && \
    apt-get update && \
    apt-get upgrade -y && \
    apt-get install -y wget \
                       build-essential \
                       g++ \
                       gcc \
                       make \
                       libbz2-dev \
                       zlib1g-dev \
                       libncurses5-dev \
                       libncursesw5-dev \
                       liblzma-dev \
                       libssl-dev \
                       zip \
                       bzip2 \
                       python3.6 \
                       python3-dev \
                       python3-pip \
                       git \
                       bsdmainutils \
                       openjdk-11-jre-headless && \
    apt-get clean && \
    apt-get purge && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*


# HTSLIB
RUN cd /usr/local/bin/ && \
    wget https://github.com/samtools/htslib/releases/download/1.9/htslib-1.9.tar.bz2 && \
    tar -xvf htslib-1.9.tar.bz2 && \
    rm htslib-1.9.tar.bz2 && \
    cd htslib-1.9 && \
    make
ENV PATH "$PATH:/usr/local/bin/htslib-1.9"


# SAMtools
RUN cd /usr/local/bin && \
    wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2 && \
    tar -xvf samtools-1.9.tar.bz2 && \
    rm samtools-1.9.tar.bz2 && \
    cd samtools-1.9 && \
    make
ENV PATH "$PATH:/usr/local/bin/samtools-1.9"

# BCFtools
RUN cd /usr/local/bin && \
    wget --continue https://github.com/samtools/bcftools/releases/download/1.9/bcftools-1.9.tar.bz2 && \
    tar -xvf bcftools-1.9.tar.bz2 && \
    rm bcftools-1.9.tar.bz2 && \
    cd bcftools-1.9 && \
    make
ENV PATH "$PATH:/usr/local/bin/bcftools-1.9"

#
## BEDtools
#RUN cd /usr/local/bin && \
#    wget https://github.com/arq5x/bedtools2/releases/download/v2.28.0/bedtools-2.28.0.tar.gz && \
#    tar -zxvf bedtools-2.28.0.tar.gz && \
#    cd bedtools2 && \
#    make

# FastQValidator
RUN cd /usr/local/bin && \
    git clone https://github.com/statgen/libStatGen.git && \
    git clone https://github.com/statgen/fastQValidator.git && \
    cd libStatGen && \
    make all && \
    cd ../fastQValidator && \
    make all
ENV PATH "$PATH:/usr/local/bin/fastQValidator/bin/"

# Picard and GATK 4 tools
RUN mkdir -p /jarfiles && \
    wget https://github.com/broadinstitute/picard/releases/download/2.20.3/picard.jar -P /jarfiles/ && \
    wget https://github.com/broadinstitute/gatk/releases/download/4.1.2.0/gatk-4.1.2.0.zip -P /jarfiles/

# VerifyBamID
RUN cd /usr/local/bin && \
    wget https://github.com/statgen/verifyBamID/releases/download/v1.1.3/verifyBamID.1.1.3.tgz && \
    tar -xvf verifyBamID.1.1.3.tgz && \
    cd verifyBamID_1.1.3 && \
    make
ENV PATH "$PATH:/usr/local/bin/verifyBamID_1.1.3/bin/"

# Fastqc
RUN cd /usr/local/bin && \
    wget http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.8.zip && \
    unzip fastqc_v0.11.8.zip && \
    rm fastqc_v0.11.8.zip && \
    chmod +x FastQC/fastqc
ENV PATH "$PATH:/usr/local/bin/FastQC"

# MultiQC
RUN git clone https://github.com/ewels/MultiQC.git && \
    cd MultiQC && \
    python3 setup.py install
ENV LC_ALL "C.UTF-8"
ENV LANG "C.UTF-8"

# Qualimap
RUN cd /usr/local/bin && \
    wget https://bitbucket.org/kokonech/qualimap/downloads/qualimap_v2.2.1.zip && \
    unzip fastqc_v0.11.8.zip && \
    rm fastqc_v0.11.8.zip && \
    chmod +x FastQC/fastqc
ENV PATH "$PATH:/usr/local/bin/FastQC"



RUN pip3 install numpy pandas scipy matplotlib seaborn biopython
RUN mkdir -p /home/src/ && mkdir -p /home/output/
COPY file_check.py verify_reference.sh fastq_read_pair_id_check.c run_file_check.sh /home/src/
WORKDIR /home/src/
RUN gcc -Wall -pthread  fastq_read_pair_id_check.c -lm -lz -std=c99 -Wextra -o fastq_pair_check
ENV PATH "$PATH:/home/src/"