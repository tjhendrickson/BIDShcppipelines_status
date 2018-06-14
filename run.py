#!/usr/bin/env python2
from __future__ import print_function
import argparse
import os
import shutil
import nibabel
from glob import glob
from shutil import rmtree
import subprocess
import nibabel as nip
from bids.grabbids import BIDSLayout
from subprocess import Popen, PIPE
from functools import partial
from collections import OrderedDict
import pdb
import csv


__version__ = open(os.path.join(os.path.dirname(os.path.realpath(__file__)),
                                'version')).read()

#month dictionary
monthDict = {'Jan': '01', 'Feb': '02', 'Mar': '03', 'Apr': '04', 'May': '05', 'Jun': '06', 'Jul': '07', 'Aug': '08',
            'Sep': '09', 'Oct': '10', 'Nov': '11', 'Dec': '12'}


def run(command, env={}, cwd=None):
    merged_env = os.environ
    merged_env.update(env)
    merged_env.pop("DEBUG", None)
    print(command)
    process = Popen(command, stdout=PIPE, stderr=subprocess.STDOUT,
                    shell=True, env=merged_env, cwd=cwd)
    while True:
        line = process.stdout.readline()
        print(line)
        line = str(line)[:-1]
        if line == '' and process.poll() != None:
            break
    if process.returncode != 0:
        raise Exception("Non zero return code: %d"%process.returncode)


def run_pre_freesurfer(**args):
    args.update(os.environ)
    cmd = "{path}/{subject}/MNINonLinear/xfms/log.txt"
    cmd = cmd.format(**args)
    if not os.path.exists(cmd):
        pre_FS = 'No'
        pre_FS_finish = 'NA'
    else:
        cmd = subprocess.check_output("stat " + cmd + "| grep 'Modify' ", shell=True)
        output = cmd

        if len(output) > 0:
            finish_year = output.split(" ")[1].split("-")[0]
            finish_month = output.split(" ")[1].split("-")[1]
            finish_day = output.split(" ")[1].split("-")[2]
            pre_FS_finish = finish_year + '/' + finish_month + '/' + finish_day
            pre_FS = 'Yes'
        else:
            pre_FS = 'No'
            pre_FS_finish = 'NA'
    return pre_FS, pre_FS_finish


def run_freesurfer(**args):
    args.update(os.environ)
    args["subjectDIR"] = os.path.join(args["path"], args["subject"], "T1w")
    cmd = "{subjectDIR}/{subject}/scripts/recon-all.log"
    cmd = cmd.format(**args)
    if not os.path.exists(cmd):
        FS = 'No'
        FS_finish = 'NA'
    else:
        cmd = subprocess.check_output("stat " + cmd + "| grep 'Modify' ", shell=True)
        output = cmd
        if output.count("\n") == 4:
            finish_year = output.split(" ")[1].split("-")[0]
            finish_month = output.split(" ")[1].split("-")[1]
            finish_day = output.split(" ")[1].split("-")[2]
            FS_finish = finish_year + '/' + finish_month + '/' + finish_day
            FS = 'Yes'
        else:
            FS = 'No'
            FS_finish = 'NA'
    return FS, FS_finish


def run_post_freesurfer(**args):
    args.update(os.environ)
    cmd = "{path}/{subject}/MNINonLinear/fsaverage_LR32k/{subject}.32k_fs_LR.wb.spec"
    cmd_set = cmd.format(**args)
    if not os.path.exists(cmd_set):
        post_FS_finish = 'NA'
        post_FS = 'No'
    else:
        output = os.path.exists(cmd_set)
        if output == True:
            fd = " {path}/{subject}/MNINonLinear/fsaverage_LR32k/{subject}.32k_fs_LR.wb.spec  "
            fd = cmd.format(**args)
            fd = subprocess.check_output("stat " + fd + "| grep 'Modify' ", shell=True)
            finish_year = fd.split(" ")[1].split("-")[0]
            finish_month = fd.split(" ")[1].split("-")[1]
            finish_day = fd.split(" ")[1].split("-")[2]
            post_FS_finish = finish_year + '/' + finish_month + '/' + finish_day
            post_FS = 'Yes'
        else:
            post_FS_finish = 'NA'
            post_FS = 'No'
    return post_FS, post_FS_finish

def run_generic_fMRI_volume_processsing(**args):
    args.update(os.environ)
    cmd = "{path}/{subject}/MNINonLinear/Results/{fmriname}/{fmriname}.nii.gz"
    cmd_set = cmd.format(**args)
    if not os.path.exists(cmd_set):
        volumefMRI_finish = 'NA'
        volumefMRI = 'No'
    else:
        output = os.path.exists(cmd_set)

        if output == True:
            fd = "{path}/{subject}/MNINonLinear/Results/{fmriname}/{fmriname}.nii.gz"
            fd = cmd.format(**args)
            fd = subprocess.check_output("stat " + fd + "| grep 'Modify' ", shell=True)
            finish_year = fd.split(" ")[1].split("-")[0]
            finish_month = fd.split(" ")[1].split("-")[1]
            finish_day = fd.split(" ")[1].split("-")[2]
            volumefMRI_finish = finish_year + '/' + finish_month + '/' + finish_day
            volumefMRI = 'Yes'
        else:
            volumefMRI_finish = 'NA'
            volumefMRI = 'No'
    return volumefMRI, volumefMRI_finish

def run_generic_fMRI_surface_processsing(**args):
    args.update(os.environ)
    cmd = "{path}/{subject}/MNINonLinear/Results/{fmriname}/{fmriname}_Atlas.dtseries.nii"
    cmd_set = cmd.format(**args)
    if not os.path.exists(cmd_set):
        surfacefMRI_finish = 'NA'
        surfacefMRI = 'No'
    else:
        output = os.path.exists(cmd_set)

        if output == True:
            fd = "{path}/{subject}/MNINonLinear/Results/{fmriname}/{fmriname}_Atlas.dtseries.nii"
            fd = cmd.format(**args)
            fd = subprocess.check_output("stat " + fd + "| grep 'Modify' ", shell=True)
            finish_year = fd.split(" ")[1].split("-")[0]
            finish_month = fd.split(" ")[1].split("-")[1]
            finish_day = fd.split(" ")[1].split("-")[2]
            surfacefMRI_finish = finish_year + '/' + finish_month + '/' + finish_day
            surfacefMRI = 'Yes'
        else:
            surfacefMRI_finish = 'NA'
            surfacefMRI = 'No'
    return surfacefMRI, surfacefMRI_finish

def run_ICAFIX_processing(**args):
    args.update(os.environ)
    cmd = "{path}/{subject}/MNINonLinear/Results/{fmriname}/{fmriname}_hp{high_pass}.ica/fix4melview_{high_pass}_thr10.txt"
    cmd_set = cmd.format(**args)
    if not os.path.exists(cmd_set):
        ICAFIX_finish = 'NA'
        ICAFIX = 'No'
    else:
        output = os.path.exists(cmd_set)

        if output == True:
            fd = "{path}/{subject}/MNINonLinear/Results/{fmriname}/{fmriname}_hp{high_pass}.ica/fix4melview_{high_pass}_thr10.txt"
            fd = cmd.format(**args)
            fd = subprocess.check_output("stat " + fd + "| grep 'Modify' ", shell=True)
            finish_year = fd.split(" ")[1].split("-")[0]
            finish_month = fd.split(" ")[1].split("-")[1]
            finish_day = fd.split(" ")[1].split("-")[2]
            ICAFIX_finish = finish_year + '/' + finish_month + '/' + finish_day
            ICAFIX = 'Yes'
        else:
            ICAFIX_finish = 'NA'
            ICAFIX = 'No'
    return ICAFIX, ICAFIX_finish

def run_PostFix_processing(**args):
    args.update(os.environ)
    cmd = "{path}/{subject}/MNINonLinear/Results/{fmriname}/{subject}_{fmriname}_ICA_Classification_singlescreen.scene"
    cmd_set = cmd.format(**args)
    if not os.path.exists(cmd_set):
        PostFix_finish = 'NA'
        PostFix = 'No'
    else:
        output = os.path.exists(cmd_set)

        if output == True:
            fd = "{path}/{subject}/MNINonLinear/Results/{fmriname}/{subject}_{fmriname}_ICA_Classification_singlescreen.scene"
            fd = cmd.format(**args)
            fd = subprocess.check_output("stat " + fd + "| grep 'Modify' ", shell=True)
            finish_year = fd.split(" ")[1].split("-")[0]
            finish_month = fd.split(" ")[1].split("-")[1]
            finish_day = fd.split(" ")[1].split("-")[2]
            PostFix_finish = finish_year + '/' + finish_month + '/' + finish_day
            PostFix = 'Yes'
        else:
            PostFix_finish = 'NA'
            PostFix = 'No'
    return PostFix, PostFix_finish

def run_RestingStateStats_processing(**args):
    args.update(os.environ)
    cmd = "{path}/{subject}/MNINonLinear/Results/{fmriname}/{fmriname}_Atlas_stats.dscalar.nii"
    cmd_set = cmd.format(**args)
    if not os.path.exists(cmd_set):
        RSS_finish = 'NA'
        RSS = 'No'
    else:
        output = os.path.exists(cmd_set)

        if output == True:
            fd = "{path}/{subject}/MNINonLinear/Results/{fmriname}/{fmriname}_Atlas_stats.dscalar.nii"
            fd = cmd.format(**args)
            fd = subprocess.check_output("stat " + fd + "| grep 'Modify' ", shell=True)
            finish_year = fd.split(" ")[1].split("-")[0]
            finish_month = fd.split(" ")[1].split("-")[1]
            finish_day = fd.split(" ")[1].split("-")[2]
            RSS_finish = finish_year + '/' + finish_month + '/' + finish_day
            RSS = 'Yes'
        else:
            RSS_finish = 'NA'
            RSS = 'No'
    return RSS, RSS_finish


def run_diffusion_processsing(**args):
    args.update(os.environ)
    cmd = "{path}/{subject}/MNINonLinear/Results/Diffusion/eddy/eddy_unwarped_images.eddy_post_eddy_shell_alignment_parameters"
    cmd_set = cmd.format(**args)
    if not os.path.exists(cmd_set):
        Diffusion_finish = 'NA'
        Diffusion = 'No'
    else:
        output = os.path.exists(cmd_set)

        if output == True:
            fd = "{path}/{subject}/MNINonLinear/Results/Diffusion/eddy/eddy_unwarped_images.eddy_post_eddy_shell_alignment_parameters"
            fd = cmd.format(**args)
            fd = subprocess.check_output("stat " + fd + "| grep 'Modify' ", shell=True)
            finish_year = fd.split(" ")[1].split("-")[0]
            finish_month = fd.split(" ")[1].split("-")[1]
            finish_day = fd.split(" ")[1].split("-")[2]
            Diffusion_finish = finish_year + '/' + finish_month + '/' + finish_day
            Diffusion = 'Yes'
        else:
            Diffusion_finish = 'NA'
            Diffusion = 'No'
    return Diffusion, Diffusion_finish

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
parser.add_argument('--n_cpus', help='Number of CPUs/cores available to use.',
                   default=1, type=int)
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


#create csv
csvfile = open(args.output_dir + '/HCP_processing_status.csv', 'w')
fieldnames = ['subject_id', 'session_id', 'PreFreeSurfer', 'Finish Date', 'FreeSurfer', 'Finish Date',
              'PostFreeSurfer', 'Finish Date', 'fMRIVolume', 'Finish Date', 'fMRISurface', 'Finish Date',
              'ICAFIX', 'Finish Date', 'PostFix', 'Finish Date', 'RestingStateStats', 'Finish Date', 'DiffusionProcessing', 'Finish Date']
writer = csv.writer(csvfile, delimiter=',')
writer.writerow(fieldnames)

# running participant level
if args.analysis_level == "participant":
    for subject_label in subjects_to_analyze:
        # if subject label has sessions underneath those need to be outputted into different directories
        if glob(os.path.join(args.bids_dir, "sub-" + subject_label, "ses-*")):
            ses_dirs = glob(os.path.join(args.bids_dir, "sub-" + subject_label, "ses-*"))
            ses_to_analyze = [ses_dir.split("-")[-1] for ses_dir in ses_dirs]
            for ses_label in ses_to_analyze:
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
                                                                             n_cpus=args.n_cpus))])

                for stage, stage_func in struct_stages_dict.iteritems():
                    if stage in args.stages:
                        if stage == "PreFreeSurfer":
                            pre_FS, pre_FS_finish = stage_func()
                        elif stage == "FreeSurfer":
                            FS, FS_finish = stage_func()
                        else:
                            post_FS, post_FS_finish = stage_func()
                bolds = [f.filename for f in layout.get(subject=subject_label, session=ses_label,
                                                        type='bold',
                                                        extensions=["nii.gz", "nii"])]
                highpass = "2000"
                training_data = "HCP_hp2000"

                for fmritcs in bolds:
                    fmriname = fmritcs.split("%s/func/" % ses_label)[-1].split(".")[0]
                    assert fmriname

                    func_stages_dict = OrderedDict([("fMRIVolume", partial(run_generic_fMRI_volume_processsing,
                                                                           path=args.output_dir + "/sub-%s" % (
                                                                               subject_label),
                                                                           subject="ses-%s" % (ses_label),
                                                                           fmriname=fmriname,
                                                                           n_cpus=args.n_cpus)),
                                                    ("fMRISurface", partial(run_generic_fMRI_surface_processsing,
                                                                            path=args.output_dir + "/sub-%s" % (
                                                                                subject_label),
                                                                            subject="ses-%s" % (ses_label),
                                                                            fmriname=fmriname,
                                                                            n_cpus=args.n_cpus))])
                    if 'rest' in fmriname:
                        rest_stages_dict = OrderedDict([("ICAFIX", partial(run_ICAFIX_processing,
                                                                           path=args.output_dir + "/sub-%s" % (
                                                                               subject_label),
                                                                           n_cpus=args.n_cpus,
                                                                           subject="ses-%s" % (ses_label),
                                                                           fmriname=fmriname,
                                                                           high_pass=highpass,
                                                                           training_data=training_data)),
                                                        ("PostFix", partial(run_PostFix_processing,
                                                                            path=args.output_dir + "/sub-%s" % (
                                                                                subject_label),
                                                                            n_cpus=args.n_cpus,
                                                                            subject="ses-%s" % (ses_label),
                                                                            fmriname=fmriname)),
                                                        ("RestingStateStats", partial(run_RestingStateStats_processing,
                                                                                      path=args.output_dir + "/sub-%s" % subject_label,
                                                                                      n_cpus=args.n_cpus,
                                                                                      subject="ses-%s" % ses_label,
                                                                                      fmriname=fmriname))])
                    # TODO: finish task fMRI portion
                    #else:
                        #task_stages_dict

                for stage, stage_func in func_stages_dict.iteritems():
                    if stage in args.stages:
                        if stage == "fMRIVolume":
                            volumefMRI, volumefMRI_finish = stage_func()
                        else:
                            surfacefMRI, surfacefMRI_finish = stage_func()
                for stage, stage_func in rest_stages_dict.iteritems():
                    if stage in args.stages:
                        if stage == "ICAFIX":
                            ICAFIX, ICAFIX_finish = stage_func()
                        elif stage == "PostFix":
                            PostFix, PostFix_finish = stage_func()
                        else:
                            RSS, RSS_finish = stage_func()

                dif_stages_dict = OrderedDict([("DiffusionPreprocessing", partial(run_diffusion_processsing,
                                                                                  path=args.output_dir + "/sub-%s" % (
                                                                                      subject_label),
                                                                                  n_cpus=args.n_cpus,
                                                                                  subject="ses-%s" % (ses_label)))])

                for stage, stage_func in dif_stages_dict.iteritems():
                    if stage in args.stages:
                        Diffusion, Diffusion_finish = stage_func()

                row = [subject_label, ses_label, pre_FS, pre_FS_finish, FS, FS_finish, post_FS, post_FS_finish,
                       volumefMRI, volumefMRI_finish, surfacefMRI, surfacefMRI_finish, ICAFIX, ICAFIX_finish,
                       PostFix, PostFix_finish, RSS, RSS, Diffusion, Diffusion_finish]
                writer.writerow(row)

                #writer.writerow({writer.fieldnames[0]: subject_label, writer.fieldnames[1]: ses_label,
                #                 writer.fieldnames[2]: pre_FS, writer.fieldnames[3]: pre_FS_finish,
                #                 writer.fieldnames[4]: FS, writer.fieldnames[5]: FS_finish,
                #                 writer.fieldnames[6]: post_FS, writer.fieldnames[7]: post_FS_finish,
                #                 writer.fieldnames[8]: volumefMRI, writer.fieldnames[9]: volumefMRI_finish,
                #                 writer.fieldnames[10]: surfacefMRI, writer.fieldnames[11]: surfacefMRI_finish,
                #                 writer.fieldnames[12]: ICAFIX, writer.fieldnames[13]: ICAFIX_finish,
                #                 writer.fieldnames[14]: PostFix, writer.fieldnames[15]: PostFix_finish,
                #                 writer.fieldnames[16]: RSS, writer.fieldnames[17]: RSS_finish,
                #                 #writer.fieldnames[18]: TaskfMRI, writer.fieldnames[19]: TaskfMRI_finish,
                #                 writer.fieldnames[20]: Diffusion, writer.fieldnames[21]: Diffusion_finish})

        #else:
            #ses_label = ' '
            """
            writer.writerow({writer.fieldnames[0]: subject_label, writer.fieldnames[1]: ses_label,
                             writer.fieldnames[2]: pre_FS, writer.fieldnames[3]: pre_FS_finish,
                             writer.fieldnames[4]: FS, writer.fieldnames[5]: FS_finish,
                             writer.fieldnames[6]: post_FS, writer.fieldnames[7]: post_FS_finish,
                             writer.fieldnames[8]: volumefMRI, writer.fieldnames[9]: volumefMRI_finish,
                             writer.fieldnames[10]: surfacefMRI, writer.fieldnames[11]: surfacefMRI_finish,
                             writer.fieldnames[12]: ICAFIX, writer.fieldnames[13]: ICAFIX_finish,
                             writer.fieldnames[14]: PostFix, writer.fieldnames[15]: PostFix_finish,
                             writer.fieldnames[16]: RSS, writer.fieldnames[17]: RSS_finish,
                             writer.fieldnames[18]: TaskfMRI, writer.fieldnames[19]: TaskfMRI_finish,
                             writer.fieldnames[20]: Diffusion, writer.fieldnames[21]: Diffusion_finish})
            """
csvfile.close()