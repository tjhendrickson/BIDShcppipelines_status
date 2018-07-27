# Use Ubuntu 14.04 LTS
Bootstrap: docker
From: ubuntu:trusty-20170119

%post
    # Install the validator 0.26.11, along with pybids 0.6.0
    apt-get update
    apt-get install -y curl
    
    curl -sSL http://neuro.debian.net/lists/trusty.us-ca.full >> /etc/apt/sources.list.d/neurodebian.sources.list
    sed -i -e 's,main *$,main contrib non-free,g' /etc/apt/sources.list.d/neurodebian.sources.list
    grep -q 'deb .* multiverse$' /etc/apt/sources.list || sed -i -e 's,universe *$,universe multiverse,g' /etc/apt/sources.list
    apt-key adv --recv-keys --keyserver hkp://pgp.mit.edu:80 0xA5D32F012649A5A9
    curl -sL https://deb.nodesource.com/setup_6.x | bash -
    
    # Install python and nibabel
    apt-get update
    apt-get install -y python-pip python-nibabel python-setuptools git nodejs
    
    npm install -g bids-validator@0.26.11
    pip install git+https://github.com/INCF/pybids.git@0.6.0

    export PYTHONPATH=""

    #make /bids_dir and /output_dir
    mkdir /bids_dir
    mkdir /output_dir

    chmod +x run.py

%files
    run.py /run.py
    version /version
    LICENSE /LICENSE


%runscript
    exec /run.py "$@"
