#!/usr/bin/env python
from __future__ import print_function
import argparse
import os
from glob import glob
import subprocess
from bids.grabbids import BIDSLayout
from subprocess import Popen, PIPE
from functools import partial
from collections import OrderedDict
import json
import pprint
import pdb


__version__ = open(os.path.join(os.path.dirname(os.path.realpath(__file__)),
                                'version')).read()

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
    else:
        cmd = subprocess.check_output("stat " + cmd + "| grep 'Modify' ", shell=True)
        output = cmd

        if len(output) > 0:
            pre_FS = 'Yes'
        else:
            pre_FS = 'No'
    return pre_FS


def run_freesurfer(**args):
    args.update(os.environ)
    args["subjectDIR"] = os.path.join(args["path"], args["subject"], "T1w")
    cmd = "{subjectDIR}/{subject}/scripts/recon-all.log"
    cmd = cmd.format(**args)
    if not os.path.exists(cmd):
        FS = 'No'
    else:
        output = subprocess.check_output("stat " + cmd + "| grep 'Modify' ", shell=True)
        if len(output) > 0:
            FS = 'Yes'
        else:
            FS = 'No'
    return FS


def run_post_freesurfer(**args):
    args.update(os.environ)
    cmd = "{path}/{subject}/MNINonLinear/fsaverage_LR32k/{subject}.32k_fs_LR.wb.spec"
    cmd = cmd.format(**args)
    if not os.path.exists(cmd):
        post_FS = 'No'
    else:
        output = subprocess.check_output("stat " + cmd + "| grep 'Modify' ", shell=True)
        if len(output) > 0:
            post_FS = 'Yes'
        else:
            post_FS = 'No'
    return post_FS


def run_generic_fMRI_volume_processsing(**args):
    args.update(os.environ)
    cmd = "{path}/{subject}/MNINonLinear/Results/{fmriname}/{fmriname}.nii.gz"
    cmd = cmd.format(**args)
    if not os.path.exists(cmd):
        volumefMRI = 'No'
    else:
        output = subprocess.check_output("stat " + cmd + "| grep 'Modify' ", shell=True)
        if len(output) > 0:
            volumefMRI = 'Yes'
        else:
            volumefMRI = 'No'
    return volumefMRI


def run_generic_fMRI_surface_processsing(**args):
    args.update(os.environ)
    cmd = "{path}/{subject}/MNINonLinear/Results/{fmriname}/{fmriname}_Atlas.dtseries.nii"
    cmd = cmd.format(**args)
    if not os.path.exists(cmd):
        surfacefMRI = 'No'
    else:
        output = subprocess.check_output("stat " + cmd + "| grep 'Modify' ", shell=True)
        if len(output) > 0:
            surfacefMRI = 'Yes'
        else:
            surfacefMRI = 'No'
    return surfacefMRI


def run_ICAFIX_processing(**args):
    args.update(os.environ)
    cmd = "{path}/{subject}/MNINonLinear/Results/{fmriname}/{fmriname}_Atlas_hp{high_pass}_clean.dtseries.nii "
    cmd = cmd.format(**args)
    if not os.path.exists(cmd):
        ICAFIX = 'No'
    else:
        output = subprocess.check_output("stat " + cmd + "| grep 'Modify' ", shell=True)
        if len(output) > 0:
            ICAFIX = 'Yes'
    return ICAFIX

def run_MSMAll_processing(**args):
    args.update(os.environ)
    cmd = "{path}/{subject}/MNINonLinear/Results/{fmriname}/{fmriname}_{MSMAll_name}_hp{high_pass}_clean.dtseries.nii"
    cmd = cmd.format(**args)
    if not os.path.exists(cmd):
        MSMAll = 'No'
    else:
        output = subprocess.check_output("stat " + cmd + "| grep 'Modify' ", shell=True)
        if len(output) > 0:
            MSMAll = 'Yes'
    return MSMAll

def run_RestingStateStats_processing(**args):
    args.update(os.environ)
    cmd = "{path}/{subject}/MNINonLinear/Results/{fmriname}/RestingStateStats/{fmriname}*.pconn.nii"
    cmd = cmd.format(**args)
    cmd = glob(cmd)
    if not cmd:
        RSS = 'No'
    else:
        output = subprocess.check_output("stat " + cmd[0] + "| grep 'Modify' ", shell=True)
        if len(output) > 0:
            RSS = 'Yes'
        else:
            RSS = 'No'
    return RSS

def run_TaskfMRIAnalysis_processing(**args):
    args.update(os.environ)
    cmd = "{path}/{subject}/MNINonLinear/Results/{fmriname}/{shortfmriname}*.feat"
    cmd = cmd.format(**args)
    cmd = glob(cmd)
    if not cmd:
        tfMRI = 'No'
    else:
        tfMRI = 'Yes'
    return tfMRI


def run_diffusion_processsing(**args):
    args.update(os.environ)
    cmd = "{path}/{subject}/Diffusion/eddy/eddy_unwarped_images.eddy_post_eddy_shell_alignment_parameters"
    cmd = cmd.format(**args)
    if not os.path.exists(cmd):
        Diffusion = 'No'
    else:
        output = subprocess.check_output("stat " + cmd + "| grep 'Modify' ", shell=True)
        if len(output) > 0:
            Diffusion = 'Yes'
        else:
            Diffusion = 'No'
    return Diffusion

def snapshot(json_file):
    with open(json_file, 'r') as f:
        data = json.load(f)
    ii=0
    jj=0
    kk=0
    ll=0
    mm=0
    session_total = len(data["Scanning Sessions"])
    fmri_total = len([ item for sublist in data["fMRINames"] for subsublist in sublist for item in subsublist])
    rsfMRI_total = len([ item for sublist in data["fMRINames"] for subsublist in sublist for item in subsublist if 'rest' == item.split("task-")[1].split("_")[0]])
    tfMRI_total = fmri_total - rsfMRI_total
    
    failed_sMRIs = []
    failed_fMRIs = []
    failed_dMRIs = []
    failed_rsfMRI = []
    failed_tfMRI = []
    for ses_counter, session_id in enumerate(data["Scanning Sessions"]):
        if data["PostFreeSurferFinish"][ses_counter][0] == 'Yes':
            ii = ii + 1
        else:
            failed_sMRIs.append(session_id)
        if data["DiffusionPreProcessingFinish"][ses_counter] == 'Yes':
            kk = kk + 1
        else:
            failed_dMRIs.append(session_id)
        for fMRI_counter, fMRI_status in enumerate(data["fMRISurfaceFinish"][ses_counter][0]):
            if fMRI_status == 'Yes':
                jj = jj + 1
            else:
                failed_fMRIs.append(data["fMRINames"][ses_counter][0][fMRI_counter])
        for RSS_counter, RSS_status in enumerate(data["RestingStateStatsFinish"][ses_counter][0]):
            if RSS_status == 'Yes':
                ll = ll + 1
            else:
                failed_rsfMRI.append(data["fMRINames"][ses_counter][0][RSS_counter])
        for tfMRI_counter, tfMRI_status in enumerate(data["TaskfMRIAnalysisFinish"][ses_counter][0]):
            if tfMRI_status == 'Yes':
                mm = mm + 1
            else:
                failed_tfMRI.append(data["fMRINames"][ses_counter][0][tfMRI_counter])

    failed_fMRIs = [str(item) for item in failed_fMRIs]
    failed_sMRIs = [str(item) for item in failed_sMRIs]
    failed_dMRIs = [str(item) for item in failed_dMRIs]
    failed_rsfMRI = [str(item) for item in failed_rsfMRI]
    failed_tfMRI = [str(item) for item in failed_tfMRI]

    if ii < session_total:
        sMRI_summary = "%s out of %s sMRI sessions completed." % (ii, session_total)
        sMRI_output = sorted(failed_sMRIs)
    else:
        sMRI_summary="All sMRI sessions have completed structural preprocessing"
        sMRI_output = ""
    if jj < fmri_total:
        fMRI_summary= "%s out of %s fMRI scans completed." %( jj, fmri_total)
        fMRI_output= sorted(failed_fMRIs)
    else:
        fMRI_summary = "All fMRI scans completed."
        fMRI_output = ""
    if kk < session_total:
        dMRI_summary = "%s out of %s dMRI sessions completed." % (kk, session_total)
        dMRI_output = sorted(failed_dMRIs)
    else:
        dMRI_summary = "All dMRI scans completed"
        dMRI_output = ""
    if ll < rsfMRI_total:
        rsfMRI_summary = "%s out of %s rsfMRI sessions completed processing." % (ll, rsfMRI_total)
        rsfMRI_output = sorted(failed_rsfMRI)
    else:
        rsfMRI_summary = "All rsfMRI processing completed"
        rsfMRI_output = ""
    if mm < tfMRI_total:
        tfMRI_summary = "%s out of %s tfMRI sessions completed first level processing." % (mm, tfMRI_total)
        tfMRI_output = sorted(failed_tfMRI)
    else:
        tfMRI_summary = "All tfMRI first level processing completed."
        tfMRI_output = ""
    print("Total number of scanning sessions: %s" %(session_total))
    pp = pprint.PrettyPrinter(indent = 4)
    print()
    print(sMRI_summary)
    print("sMRI failures: ")
    pp.pprint(sMRI_output)
    print()
    print(fMRI_summary)
    print("fMRI failures:")
    pp.pprint(fMRI_output)
    print()
    print(rsfMRI_summary)
    print("rsfMRI processing failures:")
    pp.pprint(rsfMRI_output)
    print()
    print(tfMRI_summary)
    print("tfMRI processing failures:")
    pp.pprint(tfMRI_output)
    print()
    print(dMRI_summary)
    print("dMRI failures: ")
    pp.pprint(dMRI_output)

    
    return sMRI_output, fMRI_output, rsfMRI_output, tfMRI_output, dMRI_output
    #needed_processing(failed_sMRIs,failed_fMRIs, failed_rsfMRI, failed_tfMRI, failed_dMRIs)

def main(output_dir, json_file):
    with open(output_dir + '/HCP_processing_status.json', 'w') as json_file:
        data.update({"Scanning Sessions": session_list})
        data.update({"PreFreeSurferFinish": pre_FS_list})
        data.update({"PreFreeSurferFinishTotal": pre_FS_num})
        data.update({"FreeSurferFinish": FS_list})
        data.update({"FreeSurferFinishTotal": FS_num})
        data.update({"PostFreeSurferFinish": post_FS_list})
        data.update({"PostFreeSurferFinishTotal": post_FS_num})
        data.update({"fMRINames": fmriname_list})
        data.update({"fMRIVolumeFinish": volumefMRI_list})
        data.update({"fMRIVolumeFinishTotal": fMRIVolume_num})
        data.update({"fMRISurfaceFinish": surfacefMRI_list})
        data.update({"fMRISurfaceFinishTotal": fMRISurface_num})
        data.update({"ICAFIXFinish": ICAFIX_list})
        data.update({"ICAFIXFinishTotal": ICAFIX_num})
        data.update({"MSMAllFinish": MSMAll_list})
        data.update({"MSMAllFinishTotal": MSMAll_num})
        data.update({"RestingStateStatsFinish": RestingStateStats_list})
        data.update({"RestingStateStatsFinishTotal": RestingStateStats_num})
        data.update({"TaskfMRIAnalysisFinish": tfMRI_list})
        data.update({"TaskfMRIAnalysisFinishTotal": tfMRI_num})
        data.update({"DiffusionPreProcessingFinish": Diffusion_list})
        data.update({"DiffusionPreProcessingTotal": Diffusion_num})
        json.dump(data, json_file)
        os.system("chmod a+rw " +args.output_dir + '/HCP_processing_status.json')
        json_file.close()
    sMRI_output, fMRI_output, rsfMRI_output, tfMRI_output, dMRI_output = snapshot(output_dir + '/HCP_processing_status.json')
    with open(output_dir + '/HCP_processing_status.json', 'w') as json_file:
        data.update({"sMRI failures": sMRI_output})
        data.update({"fMRI failures": fMRI_output})
        data.update({"rsfMRI failures": rsfMRI_output})
        data.update({"tfMRI failures": tfMRI_output})
        data.update({"dMRI failures": dMRI_output})
        json.dump(data, json_file)
        os.system("chmod a+rw " +args.output_dir + '/HCP_processing_status.json')
        json_file.close()
    

parser = argparse.ArgumentParser(description='HCP Pipeline status BIDS App (structural, functional MRI, diffusion, resting state).')
parser.add_argument('bids_dir', help='The directory with the input dataset '
                    'formatted according to the BIDS standard.')
parser.add_argument('output_dir', help='The directory where the HCP output files are stored. '
                    'The processing status csv file will be outptted here')
parser.add_argument('--participant_label', 
                    help='The label(s) of the participant(s) that should be analyzed. The label '
                   'corresponds to sub-<participant_label> from the BIDS spec '
                   '(so it does not include "sub-"). If this parameter is not '
                   'provided all subjects should be analyzed. Multiple '
                   'participants can be specified with a space separated list.',
                   nargs="+")
parser.add_argument('--stages', 
                   help='Which stages to run. Space separated list.',
                   nargs="+", 
                   choices=['PreFreeSurfer', 'FreeSurfer', 'PostFreeSurfer', 'fMRIVolume', 
                   'fMRISurface', 'ICAFIX', 'RestingStateStats', 'DiffusionPreprocessing', 'MSMAll', 'TaskfMRIAnalysis'],
                   default=['PreFreeSurfer', 'FreeSurfer', 'PostFreeSurfer','fMRIVolume', 
                   'fMRISurface', 'ICAFIX', 'RestingStateStats', 'DiffusionPreprocessing', 'MSMAll', 'TaskfMRIAnalysis'])
parser.add_argument('-v', '--version', action='version',
                    version='HCPPipeline Status {}'.format(__version__))

args = parser.parse_args()
layout = BIDSLayout(args.bids_dir)
subjects_to_analyze = []

# only for a subset of subjects
if args.participant_label:
    subjects_to_analyze = args.participant_label
# for all subjects
else:
    subject_dirs = glob(os.path.join(args.bids_dir, "sub-*"))
    subjects_to_analyze = [subject_dir.split("-")[-1] for subject_dir in subject_dirs]

data = {}

session_list = []
pre_FS_list = []
FS_list = []
post_FS_list = []
fmriname_list = []
volumefMRI_list = []
surfacefMRI_list = []
RestingStateStats_list = []
ICAFIX_list = []
MSMAll_list = []
tfMRI_list = []
Diffusion_list = []

pre_FS_num = 0
FS_num = 0
post_FS_num = 0
fMRITotal_num = 0
rfMRITotal_num = 0
tfMRITotal_num = 0
fMRIVolume_num = 0
fMRISurface_num = 0
ICAFIX_num = 0
MSMAll_num = 0
tfMRI_num = 0
Diffusion_num = 0
RestingStateStats_num = 0

for subject_label in subjects_to_analyze:
    # if subject label has sessions underneath those need to be outputted into different directories
    if glob(os.path.join(args.bids_dir, "sub-" + subject_label, "ses-*")):
        ses_dirs = glob(os.path.join(args.bids_dir, "sub-" + subject_label, "ses-*"))
        ses_to_analyze = [ses_dir.split("-")[-1] for ses_dir in ses_dirs]
        for ses_label in ses_to_analyze:
            session_list.append(["sub-" + subject_label + "/ses-" + ses_label])
            struct_stages_dict = OrderedDict([("PreFreeSurfer", partial(run_pre_freesurfer,
                                                                        path=args.output_dir + "/sub-%s" % (
                                                                            subject_label),
                                                                        subject="ses-%s" % (ses_label))),
                                                ("FreeSurfer", partial(run_freesurfer,
                                                                        path=args.output_dir + "/sub-%s" % (
                                                                            subject_label),
                                                                        subject="ses-%s" % (ses_label))),
                                                ("PostFreeSurfer", partial(run_post_freesurfer,
                                                                            path=args.output_dir + "/sub-%s" % (
                                                                                subject_label),
                                                                            subject="ses-%s" % (ses_label)))])

            for stage, stage_func in struct_stages_dict.iteritems():
                if stage in args.stages:
                    if stage == "PreFreeSurfer":
                        pre_FS = stage_func()
                        pre_FS_list.append([pre_FS])
                        if pre_FS == 'Yes':
                            pre_FS_num += 1
                    elif stage == "FreeSurfer":
                        FS = stage_func()
                        FS_list.append([FS])
                        if FS == 'Yes':
                            FS_num += 1
                    else:
                        post_FS= stage_func()
                        post_FS_list.append([post_FS])
                        if post_FS == 'Yes':
                            post_FS_num += 1
            
            bolds = [f.filename for f in layout.get(subject=subject_label, session=ses_label,
                                                    type='bold',
                                                    extensions=["nii.gz", "nii"])]
            highpass = "2000"
            training_data = "HCP_hp2000"
            MSMAll_name = "Atlas_MSMAll_2_d40_WRN"

            session_fmriname_list = []
            session_surfacefMRI_list = []
            session_volumefMRI_list = []
            session_RestingStateStats_list = []
            session_ICAFIX_list = []
            session_MSMAll_list = []
            session_tfMRI_list = []

            for fmritcs in bolds:
                fmriname = fmritcs.split("%s/func/" % ses_label)[-1].split(".")[0]
                assert fmriname
                session_fmriname_list.append(fmriname)
                fMRITotal_num += 1
                func_stages_dict = OrderedDict([("fMRIVolume", partial(run_generic_fMRI_volume_processsing,
                                                                        path=args.output_dir + "/sub-%s" % (subject_label),
                                                                        subject="ses-%s" % (ses_label),
                                                                        fmriname=fmriname)),
                                                ("fMRISurface", partial(run_generic_fMRI_surface_processsing,
                                                                        path=args.output_dir + "/sub-%s" % (subject_label),
                                                                        subject="ses-%s" % (ses_label),
                                                                        fmriname=fmriname)),
                                                ("ICAFIX", partial(run_ICAFIX_processing,
                                                                        path=args.output_dir + "/sub-%s" % (subject_label),
                                                                        subject="ses-%s" % (ses_label),
                                                                        fmriname=fmriname,
                                                                        high_pass=highpass,
                                                                        training_data=training_data)),
                                                ("MSMAll", partial(run_MSMAll_processing,
                                                                   path=args.output_dir + "/sub-%s" % (subject_label),
                                                                   subject="ses-%s" %(ses_label),
                                                                   fmriname=fmriname,
                                                                   high_pass=highpass,
                                                                   MSMAll_name=MSMAll_name))])
                if 'rest' in fmriname:
                    rest_stages_dict = OrderedDict([("RestingStateStats", partial(run_RestingStateStats_processing,
                                                                                    path=args.output_dir + "/sub-%s" % subject_label,
                                                                                    subject="ses-%s" % ses_label,
                                                                                    fmriname=fmriname))])
                else:
                    shortfmriname = fmriname.split("_")[2].split("-")[1]
                    task_stages_dict = OrderedDict([("tfMRIAnalysis", partial(run_TaskfMRIAnalysis_processing,
                                                                                    path=args.output_dir + "/sub-%s" % subject_label,
                                                                                    subject="ses-%s" % ses_label,
                                                                                    fmriname=fmriname,
                                                                                    shortfmriname=shortfmriname))])
                
                for stage, stage_func in func_stages_dict.iteritems():
                    if stage in args.stages:
                        if stage == "fMRIVolume":
                            volumefMRI = stage_func()
                            session_volumefMRI_list.append(volumefMRI)
                            if volumefMRI == 'Yes':
                                fMRIVolume_num += 1
                        elif stage == 'fMRISurface':
                            surfacefMRI = stage_func()
                            session_surfacefMRI_list.append(surfacefMRI)
                            if surfacefMRI == 'Yes':
                                fMRISurface_num += 1
                        elif stage == 'ICAFIX':
                            ICAFIX = stage_func()
                            session_ICAFIX_list.append(ICAFIX)
                            if ICAFIX == 'Yes':
                                ICAFIX_num += 1
                        elif stage == 'MSMAll':
                            MSMAll = stage_func()
                            session_MSMAll_list.append(MSMAll)
                            if MSMAll == 'Yes':
                                MSMAll_num += 1

                if 'rest' in fmriname:
                    for stage, stage_func in rest_stages_dict.iteritems():
                        if stage in args.stages:
                            if stage == 'RestingStateStats':
                                RSS = stage_func()
                                session_RestingStateStats_list.append(RSS)
                                if RSS == 'Yes':
                                    RestingStateStats_num += 1
                        rfMRITotal_num += 1
                else:
                    for stage, stage_func in task_stages_dict.iteritems():
                        if stage in args.stages:
                            if stage == 'TaskfMRIAnalysis':
                                tfMRI = stage_func()
                                session_tfMRI_list.append(tfMRI)
                                if tfMRI == 'Yes':
                                    tfMRI_num += 1
                        tfMRITotal_num += 1

            fmriname_list.append([session_fmriname_list])
            volumefMRI_list.append([session_volumefMRI_list])
            surfacefMRI_list.append([session_surfacefMRI_list])
            ICAFIX_list.append([session_ICAFIX_list])
            RestingStateStats_list.append([session_RestingStateStats_list])
            MSMAll_list.append([session_MSMAll_list])
            tfMRI_list.append([session_tfMRI_list])

            dif_stages_dict = OrderedDict([("DiffusionPreprocessing", partial(run_diffusion_processsing,
                                                                                path=args.output_dir + "/sub-%s" % (
                                                                                    subject_label),
                                                                                subject="ses-%s" % (ses_label)))])
            for stage, stage_func in dif_stages_dict.iteritems():
                if stage in args.stages:
                    Diffusion = stage_func()
                    Diffusion_list.append(Diffusion)
                    if Diffusion == 'Yes':
                        Diffusion_num += 1


