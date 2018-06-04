#!/usr/bin/env python3
import argparse
import os
import subprocess
import nibabel
import numpy
from glob import glob
from bids.grabbids import BIDSLayout
from functools import partial
from collections import OrderedDict
import csv

__version__ = open(os.path.join(os.path.dirname(os.path.realpath(__file__)),
                                'version')).read()

#create csv
with open('HCP_processing_status.csv', 'w') as csvfile:
    fieldnames = ['subject_id', 'session_id', 'PreFreeSurfer', 'Finish Date', 'FreeSurfer', 'Finish Date', 'PostFreeSurfer', \
                  'fMRIVolume', 'Finish Date', 'fMRISurface', 'Finish Date', 'ICAFIX', 'Finish Date', 'PostFix', ]


def run(command, env={}):
    merged_env = os.environ
    merged_env.update(env)
    process = subprocess.Popen(command, stdout=subprocess.PIPE,
                               stderr=subprocess.STDOUT, shell=True,
                               env=merged_env)
    while True:
        line = process.stdout.readline()
        line = str(line, 'utf-8')[:-1]
        print(line)
        if line == '' and process.poll() != None:
            break
    if process.returncode != 0:
        raise Exception("Non zero return code: %d"%process.returncode)

def run_pre_freesurfer(**args):
    args.update(os.environ)

    cmd = subprocess.check_output("cat {path}/{subject}/MNINonLinear/xfms/log.txt | grep 'END' ", shell=True)
    cmd = cmd.format(**args)
    output = run(cmd, cwd=args["path"], env={"OMP_NUM_THREADS": str(args["n_cpus"])})
    if len(output) > 0:



def run_freesurfer(**args):
    args.update(os.environ)
    args["subjectDIR"] = os.path.join(args["path"], args["subject"], "T1w")
    cmd = ' '

    cmd = cmd.format(**args)

    if not os.path.exists(os.path.join(args["subjectDIR"], "fsaverage")):
        shutil.copytree(os.path.join(os.environ["SUBJECTS_DIR"], "fsaverage"),
                        os.path.join(args["subjectDIR"], "fsaverage"))
    if not os.path.exists(os.path.join(args["subjectDIR"], "lh.EC_average")):
        shutil.copytree(os.path.join(os.environ["SUBJECTS_DIR"], "lh.EC_average"),
                        os.path.join(args["subjectDIR"], "lh.EC_average"))
    if not os.path.exists(os.path.join(args["subjectDIR"], "rh.EC_average")):
        shutil.copytree(os.path.join(os.environ["SUBJECTS_DIR"], "rh.EC_average"),
                        os.path.join(args["subjectDIR"], "rh.EC_average"))

    run(cmd, cwd=args["path"], env={"NSLOTS": str(args["n_cpus"]),
                                    "OMP_NUM_THREADS": str(args["n_cpus"])})

def run_post_freesurfer(**args):
    args.update(os.environ)
    cmd = ' '

    cmd = cmd.format(**args)
    run(cmd, cwd=args["path"], env={"OMP_NUM_THREADS": str(args["n_cpus"])})

def run_generic_fMRI_volume_processsing(**args):
    args.update(os.environ)
    cmd = ' '
    cmd = cmd.format(**args)
    run(cmd, cwd=args["path"], env={"OMP_NUM_THREADS": str(args["n_cpus"])})

def run_generic_fMRI_surface_processsing(**args):
    args.update(os.environ)
    cmd = ' '
    cmd = cmd.format(**args)
    run(cmd, cwd=args["path"], env={"OMP_NUM_THREADS": str(args["n_cpus"])})

def run_diffusion_processsing(**args):
    args.update(os.environ)
    cmd = ' '
    cmd = cmd.format(**args)
    run(cmd, cwd=args["path"], env={"OMP_NUM_THREADS": str(args["n_cpus"])})

def run_ICAFIX_processing(**args):
    args.update(os.environ)
    cmd = ' '
    run(cmd, cwd=args["path"], env={"OMP_NUM_THREADS": str(args["n_cpus"]),
                                    "XAPPLRESDIR": "/usr/local/R2014a/v83/X11/app-defaults"})

def run_PostFix_processing(**args):
    args.update(os.environ)
    cmd = ' '
    cmd = cmd.format(**args)
    run(cmd, cwd=args["path"], env={"OMP_NUM_THREADS": str(args["n_cpus"]),
                                    "XAPPLRESDIR": "/usr/local/R2013a/v81/X11/app-defaults",
                                    "matlab_compiler_runtime": "/usr/local/R2013a/v81"})

def run_RestingStateStats_processing(**args):
    args.update(os.environ)
    cmd = ' '
    cmd = cmd.format(**args)
    run(cmd, cwd=args["path"], env={"OMP_NUM_THREADS": str(args["n_cpus"]),
                                "XAPPLRESDIR": "/usr/local/R2013a/v81/X11/app-defaults",
                                "matlab_compiler_runtime": "/usr/local/R2013a/v81"})

parser = argparse.ArgumentParser(description='Example BIDS App entrypoint script.')
parser.add_argument('bids_dir', help='The directory with the input dataset '
                    'formatted according to the BIDS standard.')
parser.add_argument('output_dir', help='The directory where the output files '
                    'should be stored. If you are running group level analysis '
                    'this folder should be prepopulated with the results of the'
                    'participant level analysis.')
parser.add_argument('analysis_level', help='Level of the analysis that will be performed. '
                    'Multiple participant level analyses can be run independently '
                    '(in parallel) using the same output_dir.',
                    choices=['participant', 'group'])
parser.add_argument('--participant_label', help='The label(s) of the participant(s) that should be analyzed. The label '
                   'corresponds to sub-<participant_label> from the BIDS spec '
                   '(so it does not include "sub-"). If this parameter is not '
                   'provided all subjects should be analyzed. Multiple '
                   'participants can be specified with a space separated list.',
                   nargs="+")
parser.add_argument('--stages', help='Which stages to run. Space separated list.',
                   nargs="+", choices=['PreFreeSurfer', 'FreeSurfer',
                                       'PostFreeSurfer', 'fMRIVolume',
                                       'fMRISurface', 'ICAFIX', 'PostFix', 'RestingStateStats', 'DiffusionPreprocessing'],
                   default=['PreFreeSurfer', 'FreeSurfer', 'PostFreeSurfer',
                            'fMRIVolume', 'fMRISurface', 'ICAFIX', 'PostFix', 'RestingStateStats',
                            'DiffusionPreprocessing'])
parser.add_argument('-v', '--version', action='version',
                    version='HCPPipeline Status {}'.format(__version__))


args = parser.parse_args()

run("bids-validator " + args.bids_dir)

layout = BIDSLayout(args.bids_dir)

subjects_to_analyze = []
# only for a subset of subjects
if args.participant_label:
    subjects_to_analyze = args.participant_label
# for all subjects
else:
    subject_dirs = glob(os.path.join(args.output_dir, "sub-*"))
    subjects_to_analyze = [subject_dir.split("-")[-1] for subject_dir in subject_dirs]

# running participant level
if args.analysis_level == "participant":
    for subject_label in subjects_to_analyze:
        # if subject label has sessions underneath those need to be outputted into different directories
        if glob(os.path.join(args.bids_dir, "sub-" + subject_label, "ses-*")):
            ses_dirs = glob(os.path.join(args.bids_dir, "sub-" + subject_label, "ses-*"))
            ses_to_analyze = [ses_dir.split("-")[-1] for ses_dir in ses_dirs]
            for ses_label in ses_to_analyze:
                t1ws = [f.filename for f in layout.get(subject=subject_label, session=ses_label,
                                                       type='T1w',
                                                       extensions=["nii.gz", "nii"])]
                t2ws = [f.filename for f in layout.get(subject=subject_label, session=ses_label,
                                                       type='T2w',
                                                       extensions=["nii.gz", "nii"])]
                assert (len(t1ws) > 0), "No T1w files found for subject %s and session %s!" % (subject_label, ses_label)
                assert (len(t2ws) > 0), "No T2w files found for subject %s and session %s!" % (subject_label, ses_label)

                struct_stages_dict = OrderedDict([("PreFreeSurfer", partial(run_pre_freesurfer,
                                                                            path=args.output_dir + "/sub-%s" % (
                                                                                subject_label),
                                                                            subject="ses-%s" % (ses_label),
                                                                            n_cpus=args.n_cpus)),
                                                  ("FreeSurfer", partial(run_freesurfer,
                                                                         path=args.output_dir + "/sub-%s" % (
                                                                             subject_label),
                                                                         subject="ses-%s" % (ses_label),
                                                                         n_cpus=args.n_cpus)),
                                                  ("PostFreeSurfer", partial(run_post_freesurfer,
                                                                             path=args.output_dir + "/sub-%s" % (
                                                                                 subject_label),
                                                                             subject="ses-%s" % (ses_label),
                                                                             n_cpus=args.n_cpus))
                                                  ])
        else:

