
# Set the base image to Ubuntu
FROM ubuntu:22.04

# to override the interactive installation of dependencies such as tzdata
ENV DEBIAN_FRONTEND=noninteractive 

# OS updates and utilities
RUN apt-get update && apt-get upgrade -y && apt-get install build-essential unzip zlib1g-dev wget ncbi-entrez-direct fastqc trimmomatic ncbi-blast+ bowtie2 samtools -y

# setup SRA toolkit
RUN cd /opt && \
    wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/3.0.1/sratoolkit.3.0.1-ubuntu64.tar.gz  && \
    tar -xvf sratoolkit.3.0.1-ubuntu64.tar.gz -C .

RUN cd
# setup FASTX-Toolkit
RUN cd /opt && \
    wget http://hannonlab.cshl.edu/fastx_toolkit/fastx_toolkit_0.0.13_binaries_Linux_2.6_amd64.tar.bz2  && \
    tar -xvf fastx_toolkit_0.0.13_binaries_Linux_2.6_amd64.tar.bz2 -C .

RUN cd
# setup kraken2
RUN cd /opt && \
    wget https://github.com/DerrickWood/kraken2/archive/refs/tags/v2.1.2.zip && \
    unzip v2.1.2.zip && \
    cd kraken2-2.1.2 && \
    mkdir kraken2 && \
    ./install_kraken2.sh kraken2/

# adding the executables to PATH
ENV PATH /opt/sratoolkit.3.0.1-ubuntu64/bin/:$PATH
ENV PATH /opt/bin/:$PATH
ENV PATH /opt/kraken2-2.1.2/kraken2/:$PATH

RUN cd
# setting up the Perl modules required for the pipeline and default configuration
RUN cd /home && \
    wget https://github.com/maurya-anand/ARA/archive/refs/heads/main.zip && \
    unzip main.zip && \
    cd ARA-main/ && \
    perl setup.pl
    
# configuration file generated. Installation finished
RUN echo -e "Installation completed. ARA is ready to use."
