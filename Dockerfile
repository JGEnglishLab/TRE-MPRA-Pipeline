#Fully functional
FROM ubuntu:16.04 AS fastq-join

RUN apt-get update && \
	apt-get install -y git libgsl0-dev zlib1g-dev build-essential && \
	apt-get clean && \
	apt-get autoremove && \
	rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

RUN git clone https://github.com/ExpressionAnalysis/ea-utils.git && \
	cd ea-utils/clipper && make && \
	make install && \
	rm -rf ../../ea-utils/

FROM continuumio/miniconda3 AS conda-stuff

COPY environment.yml .
RUN conda env create -f environment.yml
RUN pip install conda-pack

RUN conda-pack -n sm -o /tmp/env.tar && \
  mkdir /venv && cd /venv && tar xf /tmp/env.tar && \
  rm /tmp/env.tar

RUN /venv/bin/conda-unpack


#FROM ubuntu:20.04 AS runtime
FROM ubuntu:latest AS runtime
#FROM rocker/r-ubuntu:20.04 AS runtime
#FROM r-base:4.3.1 AS rtime

COPY --from=fastq-join /usr/bin/fastq-join /usr/bin
COPY --from=conda-stuff /opt/conda/bin/python /usr/bin
COPY --from=conda-stuff /venv/bin/starcode /usr/bin
#COPY --from=conda-stuff /venv/bin/snakemake /usr/bin
COPY --from=conda-stuff /venv/bin/trim_galore /usr/bin
COPY --from=conda-stuff /venv/bin/pip /usr/bin

ENV PATH="${PATH}:/usr/bin"
ENV DEBIAN_FRONTEND=noninteractive

#RUN apt-get update && apt-get install -y --no-install-recommends build-essential r-base python3.9 python3-pip python3-setuptools python3-dev

#NEW
RUN apt-get update && apt-get install -y --no-install-recommends build-essential python3.9 python3-pip python3-setuptools python3-dev

#Install R
#This should get 4.3.1
RUN apt-get update && apt-get -y install --no-install-recommends --no-install-suggests \
        ca-certificates software-properties-common gnupg2 gnupg1 \
      && apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9 \
      && add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)-cran40/" \
      && apt-get install r-base r-base-dev -y



#END NEW

RUN apt-get update && apt-get install -y --no-install-recommends libcurl4-openssl-dev libssl-dev libfontconfig1-dev libxml2-dev libharfbuzz-dev libfribidi-dev libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev
RUN apt-get update && apt-get install -y pigz
RUN apt-get update && apt-get install bc
RUN apt-get update && apt-get install -y isal

#RUN apt-get install -y wget
WORKDIR /app

#Install Python packages
RUN pip install snakemake
RUN pip install biopython
RUN pip install pandas
RUN pip install numpy
RUN pip install cutadapt
RUN pip install datetime
RUN pip install indexed-gzip
RUN pip install pulp==2.7.0

#Install R Packages
COPY install.R /app
RUN Rscript install.R

#Copy the app
COPY . /app
ENTRYPOINT ["python3"]

