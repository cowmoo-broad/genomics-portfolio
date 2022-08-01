#This is a sample Image
FROM ubuntu
MAINTAINER cowmoomoo@gmail.com

RUN apt-get update
RUN apt-get install -y wget
RUN apt-get install -y default-jre
RUN apt-get install -y fastqc
RUN wget https://github.com/alexdobin/STAR/archive/2.7.10a.tar.gz
RUN tar -xzf 2.7.10a.tar.gz -C /opt/
RUN rm -rf 2.7.10a.tar.gz
RUN apt install build-essential -y --no-install-recommends
RUN apt-get install -y zlib1g zlib1g-dev
RUN cd /opt/STAR-2.7.10a/source/ && make
RUN cd / && wget https://github.com/deweylab/RSEM/archive/refs/tags/v1.3.3.tar.gz
RUN tar -xzf v1.3.3.tar.gz
RUN rm -rf v1.3.3.tar.gz
RUN cd /RSEM-1.3.3/ && make install
RUN apt-get install -y software-properties-common
RUN add-apt-repository -y ppa:deadsnakes/ppa
RUN apt update
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get install -y python3.9
RUN wget https://downloads.sourceforge.net/project/bowtie-bio/bowtie/1.3.1/bowtie-1.3.1-linux-x86_64.zip
RUN apt-get install -y unzip
RUN unzip bowtie-1.3.1-linux-x86_64.zip
RUN rm -rf bowtie-1.3.1-linux-x86_64.zip
RUN wget https://github.com/COMBINE-lab/salmon/releases/download/v1.9.0/salmon-1.9.0_linux_x86_64.tar.gz
RUN tar -xzf salmon-1.9.0_linux_x86_64.tar.gz
RUN rm -rf salmon-1.9.0_linux_x86_64.tar.gz
RUN wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | tee -a /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc
RUN add-apt-repository -y "deb https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)-cran40/"
RUN apt-get update
RUN apt-get install -y r-base
RUN R -e "install.packages('BiocManager')"
RUN R -e "install.packages('rmarkdown')"
RUN apt-get install -y pandoc texlive-latex-base texlive-fonts-recommended texlive-latex-extra
RUN apt-get install -y git
RUN R -e "BiocManager::install('tximport')"
RUN R -e "install.packages('readr')"
RUN R -e "BiocManager::install('tximportData')"
RUN apt-get install -y libcurl4-gnutls-dev libxml2-dev libssl-dev
RUN R -e "BiocManager::install('DESeq2')"
#RUN R -e "BiocManager::install('TxDb.Hsapiens.UCSC.hg19.knownGene')"
#RUN R -e "BiocManager::install('tximeta', update=TRUE, ask=FALSE, force=TRUE)"
RUN R -e "BiocManager::install('apeglm')"
RUN R -e "BiocManager::install('ashr')"
RUN R -e "BiocManager::install('pheatmap')"
RUN git clone https://github.com/cowmoo-broad/genomics-portfolio.git && cd genomics-portfolio && git checkout 7ba237f
COPY reports/sample_report.Rmd /genomics-portfolio/reports/sample_report.Rmd