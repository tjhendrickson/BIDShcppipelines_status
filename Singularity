# Use Ubuntu 14.04 LTS
Bootstrap: docker
From: ubuntu:trusty-20170119

%post
    # Install the validator 0.26.11, along with pybids 0.6.0
    #mkdir /dev/fuse
    chmod 777 /dev/fuse
    apt-get update -qq
    apt-get install -y curl debian-keyring
    curl -sL https://deb.nodesource.com/setup_6.x | bash -
    
    # Install python and nibabel
    apt-get update -qq
    apt-get install -y python-pip python-nibabel python-setuptools git nodejs
    npm install -g bids-validator@0.26.11
    pip install git+https://github.com/INCF/pybids.git@0.6.0 colorama
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
