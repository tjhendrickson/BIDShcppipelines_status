# Use Ubuntu 14.04 LTS
FROM ubuntu:trusty-20170119

# Install python and nibabel
RUN apt-get update && \
    apt-get install -y python3 python3-pip && \
    pip3 install nibabel && \
    apt-get remove -y python3-pip && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# Install python Dependencies
RUN apt-get update && apt-get install -y --no-install-recommends python-pip python-six python-nibabel python-setuptools
RUN pip install pybids==0.5.1
RUN pip install --upgrade pybids

## Install the validator
RUN apt-get update && \
    apt-get install -y curl && \
    curl -sL https://deb.nodesource.com/setup_6.x | bash - && \
    apt-get remove -y curl && \
    apt-get install -y nodejs

RUN npm install -g bids-validator@0.25.7

ENV PYTHONPATH=""

COPY run.py /run.py

COPY version /version

ENTRYPOINT ["/run.py"]
