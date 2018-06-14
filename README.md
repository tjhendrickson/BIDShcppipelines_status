## BIDS-apps HCPPipeline Processing Status

This repository was designed in order to determine the processing status of the HCPPipeline processing stream steps (i.e. PreFreeSurfer, FreeSurfer, PostFreeSurfer, fMRIVolume, fMRISurface, ICAFIX, PostFix, RestingStateStats, TaskfMRIAnalysis, and Diffusion). 

### Description
The HCP Pipelines product is a set of tools (primarily, but not exclusively,
shell scripts) for processing MRI images for the [Human Connectome Project](https://www.humanconnectome.org/).
Among other things, these tools implement the Minimal Preprocessing Pipeline
(MPP) described in [Glasser et al. 2013](https://www.ncbi.nlm.nih.gov/pubmed/23668970).

### Usage
This App has the following command line arguments:

```
sudo docker run -ti --rm -v /home/timothy/sandbox_DO_NOT_DELETE/BIDS/142_CIFASD_4/BIDS_output:/bids_dir -v /home/timothy/sandbox_DO_NOT_DELETE/BIDS/142_CIFASD_4/HCP_output:/output_dir tjhendrickson/bidshcppipeline_status:v0.1 -h
usage: run.py [-h]
              [--participant_label PARTICIPANT_LABEL [PARTICIPANT_LABEL ...]]
              [--stages {PreFreeSurfer,FreeSurfer,PostFreeSurfer,fMRIVolume,fMRISurface,ICAFIX,PostFix,RestingStateStats,DiffusionPreprocessing} [{PreFreeSurfer,FreeSurfer,PostFreeSurfer,fMRIVolume,fMRISurface,ICAFIX,PostFix,RestingStateStats,DiffusionPreprocessing} ...]]
              [-v]
              bids_dir output_dir {participant,group}

HCP Pipeline status BIDS App (structural, functional MRI, diffusion, resting
state).

positional arguments:
  bids_dir              The directory with the input dataset formatted
                        according to the BIDS standard.
  output_dir            The directory where the HCP output files are stored.
                        The processing status csv file will be outptted here
  {participant,group}   Level of the analysis that will be performed. Multiple
                        participant level analyses can be run independently
                        (in parallel) using the same output_dir.

optional arguments:
  -h, --help            show this help message and exit
  --participant_label PARTICIPANT_LABEL [PARTICIPANT_LABEL ...]
                        The label(s) of the participant(s) that should be
                        analyzed. The label corresponds to
                        sub-<participant_label> from the BIDS spec (so it does
                        not include "sub-"). If this parameter is not provided
                        all subjects should be analyzed. Multiple participants
                        can be specified with a space separated list.
  --stages {PreFreeSurfer,FreeSurfer,PostFreeSurfer,fMRIVolume,fMRISurface,ICAFIX,PostFix,RestingStateStats,DiffusionPreprocessing} [{PreFreeSurfer,FreeSurfer,PostFreeSurfer,fMRIVolume,fMRISurface,ICAFIX,PostFix,RestingStateStats,DiffusionPreprocessing} ...]
                        Which stages to run. Space separated list.
  -v, --version         show program's version number and exit
```

To run processing status across all participants and all stages (i.e. PreFreeSurfer, FreeSurfer, PostFreeSurfer, 
fMRIVolume, fMRISurface, ICAFIX, PostFix, RestingStateStats, and DiffusionPreprocessing):
```
sudo docker run -ti --rm -v /home/timothy/sandbox_DO_NOT_DELETE/BIDS/142_CIFASD_4/BIDS_output:/bids_dir 
-v /home/timothy/sandbox_DO_NOT_DELETE/BIDS/142_CIFASD_4/HCP_output:/output_dir tjhendrickson/bidshcppipeline_status:v0.1 
/bids_dir /output_dir participant
```
Or to run processing status for one participant and one stage:
```
sudo docker run -ti --rm -v /home/timothy/sandbox_DO_NOT_DELETE/BIDS/142_CIFASD_4/BIDS_output:/bids_dir 
-v /home/timothy/sandbox_DO_NOT_DELETE/BIDS/142_CIFASD_4/HCP_output:/output_dir tjhendrickson/bidshcppipeline_status:v0.1 
/bids_dir /output_dir participant --participant_label 01 --stages PreFreeSurfer
```







