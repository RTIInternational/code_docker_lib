# Base Image
FROM python:2.7.15-stretch

# Metadata
LABEL base.image="qcparser:latest"
LABEL version="1"
LABEL software="QCParser"
LABEL software.version="latest"
LABEL description="Bioinformatics utility parsing various NGS software log files to extract QC content and aggregate them into unified format"
LABEL tags="NGS QC"

# Maintainer
MAINTAINER Alex Waldrop <awaldrop@rti.org>

# update the OS related packages
RUN apt-get update

# install required dependencies for QCParser
RUN pip install numpy
RUN pip install scipy
RUN pip install pandas
RUN pip install matplotlib

# make directory to store tools such as QCParser
RUN mkdir tools

# get the QCParser from GitHub
RUN git clone https://github.com/alexwaldrop/QCParser.git /tools/QCParser

ENV PATH /tools/QCParser:$PATH

CMD ["python", "QCParser.py"]