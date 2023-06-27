FROM --platform=linux/amd64 ubuntu:latest

ENV DEBIAN_FRONTEND=noninteractive
ENV TZ 'Europe/Zurich'

RUN mkdir -p  /app
RUN mkdir -p /app/datadir
RUN mkdir -p /app/workdir

# Move liftOver program (linux version) to image 
COPY . /app
RUN mv /app/liftOver /usr/bin/

RUN apt-get update

# install R & Python
RUN apt-get -y install python3 python3-pip python3-setuptools python3-dev r-base libssl-dev libsasl2-dev

# install required Python libraries
RUN pip3 install click pandas numpy

# install required R libraries

## Install LabelSeg package using remotes
RUN R -e "install.packages('remotes', repos = 'https://cloud.r-project.org/')"
RUN R -e "remotes::install_github('baudisgroup/LabelSeg')"

## Other packages
RUN R -q -e "if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager',repos = 'https://cloud.r-project.org/')"
RUN R -q -e "BiocManager::install(c('aroma.core','DNAcopy','affxparser'))"

RUN R -q -e "install.packages(c('aroma.affymetrix','ACNE','pastecs','genlasso','R.utils','matrixStats','mongolite','tibble','plyr'),repos = 'https://cloud.r-project.org/')"
RUN R -q -e "install.packages('sfit',repos = 'https://henrikbengtsson.r-universe.dev')"

WORKDIR /app