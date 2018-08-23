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
from colorama import Fore, init                                                                                  
import pprint  


__version__ = open(os.path.join(os.path.dirname(os.path.realpath(__file__)),
                                'version')).read()

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
        pre_FS_finish = None
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
            pre_FS_finish = None
    return pre_FS, pre_FS_finish


def run_freesurfer(**args):
    
    args.update(os.environ)
    args["subjectDIR"] = os.path.join(args["path"], args["subject"], "T1w")
    cmd = "{subjectDIR}/{subject}/scripts/recon-all.log"
    cmd = cmd.format(**args)
    if not os.path.exists(cmd):
        FS = 'No'
        FS_finish = None
    else:
        output = subprocess.check_output("stat " + cmd + "| grep 'Modify' ", shell=True)
        if len(output) > 0:
            finish_year = output.split(" ")[1].split("-")[0]
            finish_month = output.split(" ")[1].split("-")[1]
            finish_day = output.split(" ")[1].split("-")[2]
            FS = 'Yes'
            FS_finish = finish_year + '/' + finish_month + '/' + finish_day
        else:
            FS = 'No'
            FS_finish = None
    return FS, FS_finish


def run_post_freesurfer(**args):
    args.update(os.environ)
    cmd = "{path}/{subject}/MNINonLinear/fsaverage_LR32k/{subject}.32k_fs_LR.wb.spec"
    cmd = cmd.format(**args)
    if not os.path.exists(cmd):
        post_FS_finish = None
        post_FS = 'No'
    else:
        output = subprocess.check_output("stat " + cmd + "| grep 'Modify' ", shell=True)
        if len(output) > 0:
            finish_year = output.split(" ")[1].split("-")[0]
            finish_month = output.split(" ")[1].split("-")[1]
            finish_day = output.split(" ")[1].split("-")[2]
            post_FS_finish = finish_year + '/' + finish_month + '/' + finish_day
            post_FS = 'Yes'
        else:
            post_FS_finish = None
            post_FS = 'No'
    return post_FS, post_FS_finish


def run_generic_fMRI_volume_processsing(**args):
    args.update(os.environ)
    cmd = "{path}/{subject}/MNINonLinear/Results/{fmriname}/{fmriname}.nii.gz"
    cmd = cmd.format(**args)
    if not os.path.exists(cmd):
        volumefMRI_finish = None
        volumefMRI = 'No'
    else:
        output = subprocess.check_output("stat " + cmd + "| grep 'Modify' ", shell=True)
        if len(output) > 0:
            finish_year = output.split(" ")[1].split("-")[0]
            finish_month = output.split(" ")[1].split("-")[1]
            finish_day = output.split(" ")[1].split("-")[2]
            volumefMRI_finish = finish_year + '/' + finish_month + '/' + finish_day
            volumefMRI = 'Yes'
        else:
            volumefMRI_finish = None
            volumefMRI = 'No'
    return volumefMRI, volumefMRI_finish


def run_generic_fMRI_surface_processsing(**args):
    args.update(os.environ)
    cmd = "{path}/{subject}/MNINonLinear/Results/{fmriname}/{fmriname}_Atlas.dtseries.nii"
    cmd = cmd.format(**args)
    if not os.path.exists(cmd):
        surfacefMRI_finish = None
        surfacefMRI = 'No'
    else:
        output = subprocess.check_output("stat " + cmd + "| grep 'Modify' ", shell=True)
        if len(output) > 0:
            finish_year = output.split(" ")[1].split("-")[0]
            finish_month = output.split(" ")[1].split("-")[1]
            finish_day = output.split(" ")[1].split("-")[2]
            surfacefMRI_finish = finish_year + '/' + finish_month + '/' + finish_day
            surfacefMRI = 'Yes'
        else:
            surfacefMRI_finish = None
            surfacefMRI = 'No'
    return surfacefMRI, surfacefMRI_finish


def run_ICAFIX_processing(**args):
    args.update(os.environ)
    cmd = "{path}/{subject}/MNINonLinear/Results/{fmriname}/{fmriname}_hp{high_pass}.ica/fix4melview_{training_data}_thr10.txt"
    cmd = cmd.format(**args)
    if not os.path.exists(cmd):
        ICAFIX_finish = None
        ICAFIX = 'No'
    else:
        output = subprocess.check_output("stat " + cmd + "| grep 'Modify' ", shell=True)
        if len(output) > 0:
            finish_year = output.split(" ")[1].split("-")[0]
            finish_month = output.split(" ")[1].split("-")[1]
            finish_day = output.split(" ")[1].split("-")[2]
            ICAFIX_finish = finish_year + '/' + finish_month + '/' + finish_day
            ICAFIX = 'Yes'
        else:
            ICAFIX_finish = None
            ICAFIX = 'No'
    return ICAFIX, ICAFIX_finish


def run_PostFix_processing(**args):
    args.update(os.environ)
    cmd = "{path}/{subject}/MNINonLinear/Results/{fmriname}/{subject}_{fmriname}_ICA_Classification_singlescreen.scene"
    cmd = cmd.format(**args)
    if not os.path.exists(cmd):
        PostFix_finish = None
        PostFix = 'No'
    else:
        output = subprocess.check_output("stat " + cmd + "| grep 'Modify' ", shell=True)
        if len(output) > 0:
            finish_year = output.split(" ")[1].split("-")[0]
            finish_month = output.split(" ")[1].split("-")[1]
            finish_day = output.split(" ")[1].split("-")[2]
            PostFix_finish = finish_year + '/' + finish_month + '/' + finish_day
            PostFix = 'Yes'
        else:
            PostFix_finish = None
            PostFix = 'No'
    return PostFix, PostFix_finish


def run_RestingStateStats_processing(**args):
    args.update(os.environ)
    cmd = "{path}/{subject}/MNINonLinear/Results/{fmriname}/{fmriname}_Atlas_stats.dscalar.nii"
    cmd = cmd.format(**args)
    if not os.path.exists(cmd):
        RSS_finish = None
        RSS = 'No'
    else:
        output = subprocess.check_output("stat " + cmd + "| grep 'Modify' ", shell=True)
        if len(output) > 0:
            finish_year = output.split(" ")[1].split("-")[0]
            finish_month = output.split(" ")[1].split("-")[1]
            finish_day = output.split(" ")[1].split("-")[2]
            RSS_finish = finish_year + '/' + finish_month + '/' + finish_day
            RSS = 'Yes'
        else:
            RSS_finish = None
            RSS = 'No'
    return RSS, RSS_finish


def run_diffusion_processsing(**args):
    args.update(os.environ)
    cmd = "{path}/{subject}/Diffusion/eddy/eddy_unwarped_images.eddy_post_eddy_shell_alignment_parameters"
    cmd = cmd.format(**args)
    if not os.path.exists(cmd):
        Diffusion_finish = None
        Diffusion = 'No'
    else:
        output = subprocess.check_output("stat " + cmd + "| grep 'Modify' ", shell=True)
        if len(output) > 0:
            finish_year = output.split(" ")[1].split("-")[0]
            finish_month = output.split(" ")[1].split("-")[1]
            finish_day = output.split(" ")[1].split("-")[2]
            Diffusion_finish = finish_year + '/' + finish_month + '/' + finish_day
            Diffusion = 'Yes'
        else:
            Diffusion_finish = None
            Diffusion = 'No'
    return Diffusion, Diffusion_finish
def snapshot(json_file):                                                                                                                                                                                                                                                                      
                                                                               
    data = json_file                                                                                    
    ii=0                                                                                                       
    jj=0             
    kk=0                                                                                                       
    session_total = len(data["Scanning Sessions"])                                                             
    fmri_total = len([ item for sublist in data["fMRINames"] for subsublist in sublist for item in subsublist])
    failed_structs = []                                                                                        
    failed_fMRIs = []                                                                                          
    failed_dMRIs = []                                                                                          
    for ses_counter, session_id in enumerate(data["Scanning Sessions"]):                                       
        if data["PostFreeSurferFinish"][ses_counter][0] == 'Yes':                                              
            ii = ii + 1                                                                                        
        else:                                                                                                  
            failed_structs.append(session_id)                                                                  
        if data["DiffusionPreProcessingFinish"][ses_counter][0] == 'Yes':                                      
            kk = kk + 1                                                                                        
        else:                                                                                                  
            failed_dMRIs.append(session_id)                                                   
        for fMRI_counter, fMRI_status in enumerate(data["fMRISurfaceFinish"][ses_counter][0]):
            if fMRI_status == 'Yes':                                                          
                jj = jj + 1                                                                   
            else:                                                                             
                failed_fMRIs.append(data["fMRINames"][ses_counter][0][fMRI_counter])          
                                                                                            
    failed_fMRIs = [str(item) for item in failed_fMRIs]                                       
    failed_structs = [str(item) for item in failed_structs]                                   
    failed_dMRIs = [str(item) for item in failed_dMRIs]                                       
    init()                                                                                    
    pp = pprint.PrettyPrinter(indent = 4)                                                     
    print(Fore.WHITE + "Total number of scanning sessions: %s" %(session_total))              
    print                                                                           
    print(Fore.WHITE + "Structural Preprocessing Summary:")                               
    if ii < session_total:                                                                    
        print(Fore.RED + "%s out of %s session completed correctly" % (ii, session_total))    
        print(Fore.RED + "The sessions with failed or yet to be run structural prepocessing:")
        pp.pprint(failed_structs)                                                             
    else:                                                                                     
        print(Fore.LIGHTGREEN_EX + "All sessions have completed structural preprocessing")    
    print                                                                                     
    print(Fore.WHITE + "Basic fMRI Preprocessing Summary:")                                   
    if jj < fmri_total:                                                                       
        print(Fore.RED + "%s out of %s fMRI scans completed correctly" %(jj, fmri_total))     
        print(Fore.RED + "The failed or yet to be run fMRI Scans:")                           
        pp.pprint(failed_fMRIs)                                                               
    else:                                                                                     
        print(Fore.LIGHTGREEN_EX + "All fMRI scans completed correctly")                      
    print(Fore.WHITE + "Diffusion Preprocessing Summary:")                                 
    if kk < session_total:                                                                   
        print(Fore.RED + "%s out of %s sessions completed correctly" % (kk, session_total))  
        print(Fore.RED + "The sessions with failed or yet to be run diffusion prepocessing:")
        pp.pprint(failed_dMRIs)                                                              
    else:                                                                                    
        print(Fore.LIGHTGREEN_EX + "All sessions have completed diffusion preprocessing")
        



        





parser = argparse.ArgumentParser(description='HCP Pipeline status BIDS App (structural, functional MRI, diffusion, resting state).')
parser.add_argument('bids_dir', help='The directory with the input dataset '
                    'formatted according to the BIDS standard.')
parser.add_argument('output_dir', help='The directory where the HCP output files are stored. '
                    'The processing status csv file will be outptted here')
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
layout = BIDSLayout(args.bids_dir)
subjects_to_analyze = []

# only for a subset of subjects
if args.participant_label:
    subjects_to_analyze = args.participant_label
# for all subjects
else:
    subject_dirs = glob(os.path.join(args.output_dir, "sub-*"))
    subjects_to_analyze = [subject_dir.split("-")[-1] for subject_dir in subject_dirs]

data = {}
session_list = []
pre_FS_list = []
pre_FS_finish_list = []
FS_list = []
FS_finish_list = []
post_FS_list = []
post_FS_finish_list = []
fmriname_list = []
volumefMRI_list = []
volumefMRI_finish_list = []
surfacefMRI_list = []
surfacefMRI_finish_list = []
RestingStateStats_list = []
RestingStateStats_finish_list = []
ICAFIX_list = []
ICAFIX_finish_list = []
PostFix_list = []
PostFix_finish_list = []
Diffusion_list = []
Diffusion_finish_list = []
pre_FS_num = 0
FS_num = 0
post_FS_num = 0
fMRITotal_num = 0
rfMRITotal_num = 0
fMRIVolume_num = 0
fMRISurface_num = 0
ICAFIX_num = 0
PostFix_num = 0
Diffusion_num = 0
RestingStateStats_num = 0
# running participant level
if args.analysis_level == "participant":
    for subject_label in subjects_to_analyze:
        # if subject label has sessions underneath those need to be outputted into different directories
        if glob(os.path.join(args.bids_dir, "sub-" + subject_label, "ses-*")):
            ses_dirs = glob(os.path.join(args.bids_dir, "sub-" + subject_label, "ses-*"))
            ses_to_analyze = [ses_dir.split("-")[-1] for ses_dir in ses_dirs]
            for ses_label in ses_to_analyze:
                session_list.append([ses_label])
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
                            pre_FS, pre_FS_finish = stage_func()
                            pre_FS_list.append([pre_FS])
                            pre_FS_finish_list.append([pre_FS_finish])
                            if pre_FS == 'Yes':
                                pre_FS_num += 1
                        elif stage == "FreeSurfer":
                            FS, FS_finish = stage_func()
                            FS_list.append([FS])
                            FS_finish_list.append([FS_finish])
                            if FS == 'Yes':
                                FS_num += 1
                        else:
                            post_FS, post_FS_finish = stage_func()
                            post_FS_list.append([post_FS])
                            post_FS_finish_list.append([post_FS_finish])
                            if post_FS == 'Yes':
                                post_FS_num += 1
                bolds = [f.filename for f in layout.get(subject=subject_label, session=ses_label,
                                                        type='bold',
                                                        extensions=["nii.gz", "nii"])]
                highpass = "2000"
                training_data = "HCP_hp2000"

                session_fmriname_list = []
                session_surfacefMRI_list = []
                session_surfacefMRI_finish_list = []
                session_volumefMRI_list = []
                session_volumefMRI_finish_list = []
                session_RestingStateStats_list = []
                session_RestingStateStats_finish_list = []
                session_ICAFIX_list = []
                session_ICAFIX_finish_list = []
                session_PostFix_list = []
                session_PostFix_finish_list = []

                for fmritcs in bolds:
                    fmriname = fmritcs.split("%s/func/" % ses_label)[-1].split(".")[0]
                    assert fmriname
                    session_fmriname_list.append(fmriname)

                    func_stages_dict = OrderedDict([("fMRIVolume", partial(run_generic_fMRI_volume_processsing,
                                                                           path=args.output_dir + "/sub-%s" % (
                                                                               subject_label),
                                                                           subject="ses-%s" % (ses_label),
                                                                           fmriname=fmriname)),
                                                    ("fMRISurface", partial(run_generic_fMRI_surface_processsing,
                                                                            path=args.output_dir + "/sub-%s" % (
                                                                                subject_label),
                                                                            subject="ses-%s" % (ses_label),
                                                                            fmriname=fmriname))])
                    if 'rest' in fmriname:
                        rest_stages_dict = OrderedDict([("ICAFIX", partial(run_ICAFIX_processing,
                                                                           path=args.output_dir + "/sub-%s" % (
                                                                               subject_label),
                                                                           subject="ses-%s" % (ses_label),
                                                                           fmriname=fmriname,
                                                                           high_pass=highpass,
                                                                           training_data=training_data)),
                                                        ("PostFix", partial(run_PostFix_processing,
                                                                            path=args.output_dir + "/sub-%s" % (
                                                                                subject_label),
                                                                            subject="ses-%s" % (ses_label),
                                                                            fmriname=fmriname)),
                                                        ("RestingStateStats", partial(run_RestingStateStats_processing,
                                                                                      path=args.output_dir + "/sub-%s" % subject_label,
                                                                                      subject="ses-%s" % ses_label,
                                                                                      fmriname=fmriname))])
                    # TODO: finish task fMRI portion
                    for stage, stage_func in func_stages_dict.iteritems():
                        if stage in args.stages:
                            if stage == "fMRIVolume":
                                volumefMRI, volumefMRI_finish = stage_func()
                                session_volumefMRI_list.append(volumefMRI)
                                session_volumefMRI_finish_list.append(volumefMRI_finish)
                                if volumefMRI == 'Yes':
                                    fMRIVolume_num += 1
                            else:
                                surfacefMRI, surfacefMRI_finish = stage_func()
                                session_surfacefMRI_list.append(surfacefMRI)
                                session_surfacefMRI_finish_list.append(surfacefMRI_finish)
                                if surfacefMRI == 'Yes':
                                    fMRISurface_num += 1
                            fMRITotal_num += 1


                    if 'rest' in fmriname:
                        for stage, stage_func in rest_stages_dict.iteritems():
                            if stage in args.stages:
                                if stage == "ICAFIX":
                                    ICAFIX, ICAFIX_finish = stage_func()
                                    session_ICAFIX_list.append(ICAFIX)
                                    session_ICAFIX_finish_list.append(ICAFIX_finish)
                                    if ICAFIX == 'Yes':
                                        ICAFIX_num += 1
                                elif stage == "PostFix":
                                    PostFix, PostFix_finish = stage_func()
                                    session_PostFix_list.append(PostFix)
                                    session_PostFix_finish_list.append(PostFix_finish)
                                    if PostFix == 'Yes':
                                        PostFix_num += 1
                                else:
                                    RSS, RSS_finish = stage_func()
                                    session_RestingStateStats_list.append(RSS)
                                    session_RestingStateStats_finish_list.append(RSS_finish)
                                    if RSS == 'Yes':
                                        RestingStateStats_num += 1
                            rfMRITotal_num += 1

                fmriname_list.append([session_fmriname_list])

                volumefMRI_list.append([session_volumefMRI_list])
                volumefMRI_finish_list.append([session_volumefMRI_finish_list])

                surfacefMRI_list.append([session_surfacefMRI_list])
                surfacefMRI_finish_list.append([session_surfacefMRI_finish_list])

                ICAFIX_list.append([session_ICAFIX_list])
                ICAFIX_finish_list.append([session_ICAFIX_finish_list])

                PostFix_list.append([session_PostFix_list])
                PostFix_finish_list.append([session_PostFix_finish_list])

                RestingStateStats_list.append([session_RestingStateStats_list])
                RestingStateStats_finish_list.append([session_RestingStateStats_finish_list])

                dif_stages_dict = OrderedDict([("DiffusionPreprocessing", partial(run_diffusion_processsing,
                                                                                  path=args.output_dir + "/sub-%s" % (
                                                                                      subject_label),
                                                                                  subject="ses-%s" % (ses_label)))])
                for stage, stage_func in dif_stages_dict.iteritems():
                    if stage in args.stages:
                        Diffusion, Diffusion_finish = stage_func()
                        Diffusion_list.append(Diffusion)
                        Diffusion_finish_list.append(Diffusion_finish)
                        if Diffusion == 'Yes':
                            Diffusion_num += 1

    with open(args.output_dir + '/HCP_processing_status.json', 'w') as json_file:
        data.update({"Scanning Sessions": session_list})
        data.update({"PreFreeSurferFinish": pre_FS_list})
        data.update({"PreFreeSurferFinishDate": pre_FS_finish_list})
        data.update({"PreFreeSurferFinishTotal": pre_FS_num})
        data.update({"FreeSurferFinish": FS_list})
        data.update({"FreeSurferFinishDate": FS_finish_list})
        data.update({"FreeSurferFinishTotal": FS_num})
        data.update({"PostFreeSurferFinish": post_FS_list})
        data.update({"PostFreeSurferFinishDate": post_FS_finish_list})
        data.update({"PostFreeSurferFinishTotal": post_FS_num})
        data.update({"fMRINames": fmriname_list})
        data.update({"fMRIVolumeFinish": volumefMRI_list})
        data.update({"fMRIVolumeFinishDate": volumefMRI_finish_list})
        data.update({"fMRIVolumeFinishTotal": fMRIVolume_num})
        data.update({"fMRISurfaceFinish": surfacefMRI_list})
        data.update({"fMRISurfaceFinishDate": surfacefMRI_finish_list})
        data.update({"fMRISurfaceFinishTotal": fMRISurface_num})
        data.update({"ICAFIXFinish": ICAFIX_list})
        data.update({"ICAFIXFinishDate": ICAFIX_finish_list})
        data.update({"ICAFIXFinishTotal": ICAFIX_num})
        data.update({"PostFixFinish": PostFix_list})
        data.update({"PostFixFinishDate": PostFix_finish_list})
        data.update({"PostFixFinishTotal": PostFix_num})
        data.update({"RestingStateStatsFinish": RestingStateStats_list})
        data.update({"RestingStateStatsFinishDate": RestingStateStats_finish_list})
        data.update({"RestingStateStatsFinishTotal": RestingStateStats_num})
        data.update({"DiffusionPreProcessingFinish": Diffusion_list})
        data.update({"DiffusionPreProcessingFinishDate": Diffusion_finish_list})
        data.update({"DiffusionPreProcessingTotal": Diffusion_num})
        json.dump(data, json_file)
        os.system("chmod 777 " + args.output_dir + '/HCP_processing_status.json')
        snapshot(json_file)
        json_file.close()

