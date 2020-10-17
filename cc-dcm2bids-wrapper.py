#!/usr/bin/env python
# -*- coding: utf-8 -*-


from datetime import datetime
import pydicom
import socket
import glob
import sys
import os
import re
import subprocess
import pandas as pd

from colorama import Fore as fg, Back as bg, Style as st
from tabulate import tabulate
import argparse


from collections import OrderedDict
from copy import deepcopy
from pprint import pprint

from pydicom.datadict import DicomDictionary, keyword_dict

def is_int_geq_zero(value):
    ivalue = int(value)
    if ivalue <= 0:
         raise argparse.ArgumentTypeError("%s is an invalid positive int value" % value)
    return ivalue


description = '''%(prog)s:
 - explore DICOM directories,
 - display summaries for the data discovered,
 - optionally convert to NIFTI (BIDS) data format.

Tested on GE Discovery 750 MR scanner.

%(prog)s is free software licensed under GPL 3.0.  Assembled from
various code snippets found on the web by Jan Nikadon (nikadon <at>
gmail <dot> com).

 '''


epilog = 'Examples:\n\n''' + \
    st.BRIGHT + os.path.basename(__file__) + st.RESET_ALL + ' --help\n' + \
    st.BRIGHT + os.path.basename(__file__) + st.RESET_ALL + ' -dm -pes -lrw -cn\n' + \
    st.BRIGHT + os.path.basename(__file__) + st.RESET_ALL + ' -pes\n' + \
    st.BRIGHT + os.path.basename(__file__) + st.RESET_ALL + ' -lr\n'



par = argparse.ArgumentParser(
    add_help=False,
    prog=st.BRIGHT + os.path.basename(__file__) + st.RESET_ALL,
    description=description,
    formatter_class=argparse.RawDescriptionHelpFormatter,
    # usage='%(prog)s [options]',
    epilog=epilog,
)
par.add_argument( # OK
    "-O",
    "--origin",
    metavar="PATH",
    dest="origin",
    default=os.getcwd(),
    help='Start directory, e.g., ' + fg.GREEN + '/home/username/study/dicom' + st.RESET_ALL +  '. The default is ' + fg.GREEN +  '<current directory>' + st.RESET_ALL + ' e.g., ' + fg.BLUE +   '%(default)s' + st.RESET_ALL,
)
par.add_argument( # OK
    "-D",
    "--dicom-dir",
    metavar="DIRECTORY",
    dest="dicom_dir",
    default='dicom',
    help='DICOM (' + fg.RED + 'source' + st.RESET_ALL + ') directory. The default is ' + fg.GREEN +   '%(default)s' + st.RESET_ALL,
)
par.add_argument( # ??
    "-N",
    "--nifti-dir",
    metavar="DIRECTORY",
    dest="nifti_dir",
    default='nifti',
    help='NIFTI (' + fg.RED + 'target' + st.RESET_ALL + ') directory. The default is ' + fg.GREEN +   '%(default)s' + st.RESET_ALL,
)
par.add_argument( # OK
    "-p",
    "--show-patients",
    action="store_true",
    dest="show_patients",
    default=False,
    help='Display patients (subjects) discovered.',
)
par.add_argument( # OK
    "-P",
    "--patients",
    metavar="PATTERN",
    dest="patients",
    default='sub-*',
    help='Unix-like pattern that matches patient (subject) directories, e.g., ' + fg.GREEN + 'p*' + st.RESET_ALL +  ' or ' + fg.GREEN + 'group-*/sub-*' + st.RESET_ALL +  '. The default is ' + fg.GREEN +  '%(default)s' + st.RESET_ALL,
)
par.add_argument( # OK
    "-e",
    "--show-exams",
    action="store_true",
    dest="show_exams",
    default=False,
    help='Display exams (sessions) discovered.',
)
par.add_argument( # OK
    "-E",
    "--exams",
    metavar="PATTERN",
    dest="exams",
    default='ses-*',
    help='Unix-like pattern that matches exam (session) directories, e.g., ' + fg.GREEN + 'e*' + st.RESET_ALL +  '. The default is ' + fg.GREEN +  '%(default)s' + st.RESET_ALL,
)
par.add_argument( # OK
    "-s",
    "--show-series",
    action="store_true",
    dest="show_series",
    default=False,
    help='Display series (sequences) discovered.',
)
par.add_argument( # ??
    "-S",
    "--series",
    metavar="PATTERN",
    dest="series",
    default='s*',
    help='Unix-like pattern that matches series (sequence) directories, e.g., ' + fg.GREEN +  '%(default)s' + st.RESET_ALL + ' (which also is the default).',
)
par.add_argument( # OK
    "-f",
    "--show-files",
    action="store_true",
    dest="show_files",
    default=False,
    help='Display DICOM files discovered.',
)
par.add_argument( # ??
    "-F",
    "--files",
    metavar="PATTERN",
    nargs='+',
    dest="files",
    default=['i*MRDC*','i*SCPT*'],
    help='Unix-like pattern that matches DICOM files, e.g., ' + fg.GREEN +  '\'*\'' + st.RESET_ALL + '. Multiple patterns are accepted. The default is ' + fg.GREEN +  '\'i*MRDC*\' \'i*SCPT*\'' + st.RESET_ALL + '.',
)
par.add_argument( # OK
    "-l",
    "--log",
    action="store_true",
    dest="log",
    default=False,
    help='Save log to file.',
)
par.add_argument( # OK
    "-L",
    "--log-file",
    metavar="FILE",
    dest="log_file",
    default='z---' + os.path.basename(__file__) + '---' + socket.gethostname() + '---' + str(datetime.utcnow().strftime('%Y%m%d-%H%M%S-%f')[:-3]) + '---log.org',
    help='Log output filename. The default is based on filename executed, hostname, date and time, i.e., ' + fg.GREEN + 'z---<executable-filename>---<hostname>---<YYYYmmdd-HH-MM-SS-sss>---log.org' + st.RESET_ALL + ' e.g., ' + fg.BLUE +   '%(default)s' + st.RESET_ALL + '.',
)
par.add_argument( # OK
    "-t",
    "--table",
    action="store_true",
    dest="table",
    default=False,
    help='Produce detailed table on DICOM files found. The default is ' + fg.GREEN +  '%(default)s' + st.RESET_ALL + '.',
)
par.add_argument( # OK
    "-w",
    "--write-csv",
    action="store_true",
    dest="write_csv",
    default=False,
    help='Write DICOM information to CSV file. Implies ' + fg.RED +  '-r' + st.RESET_ALL + '.',
)
par.add_argument( # OK
    "-W",
    "--write-csv-file",
    metavar="FILE",
    dest="write_csv_file",
    default='z---' + os.path.basename(__file__) + '---' + socket.gethostname() + '---' + str(datetime.utcnow().strftime('%Y%m%d-%H%M%S-%f')[:-3]) + '---log.csv',
    help='File name for DICOM information to CSV. The default is basename oflog-file with \'.csv\' extension.',
)
par.add_argument( # OK
    "-k",
    "--keep-sequence-information",
    action="store_true",
    dest="keep_sequence_information",
    default=False,
    help='Keep series information in the table. The default is ' + fg.GREEN +  '%(default)s' + st.RESET_ALL,
)
par.add_argument( # OK
    "-i",
    "--include-stacked-screen-saves",
    action="store_true",
    dest="include_stacked_screen_saves",
    default=False,
    help='Include series containing stacked screen saves in the table. The default is ' + fg.GREEN +  '%(default)s' + st.RESET_ALL,
)
par.add_argument( # OK
    "-m",
    "--more-info",
    action="store_true",
    dest="more_info",
    default=False,
    help='Display more information derrived from DICOM files.',
)
par.add_argument( # OK
    "-r",
    "--report-files",
    action="store_true",
    dest="report_files",
    default=False,
    help='Display slice timing related information. Implies ' + fg.RED + '-tkm' + st.RESET_ALL + '.',
)
# par.add_argument( # OK
#     "-R",
#     "--report-num-files",
#     metavar="FILE",
#     type=is_int_geq_zero,
#     dest="report_num_files",
#     default=0,
#     help='How many file to include. Implies ' + fg.RED + '-r' + st.RESET_ALL + '.',
# )
par.add_argument( # OK
    "-c",
    "--convert-demo",
    action="store_true",
    dest="convert_demo",
    default=False,
    help='Show commands used to convert discovered series (sequences) from DICOM to NIFTI format. The default is ' + fg.GREEN +  '%(default)s' + st.RESET_ALL + '.',
)
par.add_argument( # OK
    "-C",
    "--convert",
    action="store_true",
    dest="convert",
    default=False,
    help='Convert discovered series (sequences) from DICOM to NIFTI format. The default is ' + fg.GREEN +  '%(default)s' + st.RESET_ALL + '. Implies ' + fg.RED +  '-c' + st.RESET_ALL + '.',
)
par.add_argument( # OK
    "-J",
    "--study-specific-json",
    metavar="FILE",
    dest="study_specific_json",
    default=None,
    help='Study specific JSON file dontaining series description e.g., ' + fg.GREEN + 'study-specific.json' + st.RESET_ALL + '. If none provided session-specific JSON file(s) are expected to be located in <current directory> with names corresponding to session(s) E.g., ' + fg.GREEN +  'ses-01pre.json' + st.RESET_ALL + ', ' + fg.GREEN +  'ses-02post.json' + st.RESET_ALL + ' ETC.',
)
par.add_argument( # FIXME
    "-n",
    "--no-defaceing",
    action="store_true",
    dest="no_defaceing",
    default=False,
    help='Turn off defacing during conversion. The default is ' + fg.GREEN +  '%(default)s' + st.RESET_ALL,
)
par.add_argument( # ??
    "-g",
    "--GE",
    action="store_true",
    dest="ge",
    default=False,
    help='Same as ' + fg.GREEN + '-P "p*" -E "e*" -S "s*" ' + st.RESET_ALL + ' (patterns tailored to raw data obtained from GE MR750 scanner running DV 2.4)',
)
par.add_argument( # ??
    "-b",
    "--BIDS",
    action="store_true",
    dest="bids",
    default=False,
    help='Same as ' + fg.GREEN + '-P "sub-*" -E "ses-*" -S "*" ' + st.RESET_ALL + ' (patterns tailored to BIDS data structure)',
)
par.add_argument( # OK
    "-q",
    "--quiet",
    action="store_true",
    dest="quiet",
    default=False,
    help='Be (quite) quiet.',
)
par.add_argument( # OK
    "-v",
    "--verbose",
    action="count",
    dest="verbose",
    default=0,
    help='Verbose (use multiple times to increase verbbosity up to -vvvv).',
)
par.add_argument( # OK
    "-d",
    "--debug",
    action="count",
    dest="debug",
    default=0,
    help='Debug mode (use multiple times to increase debugging level up to -ddd).',
)
par.add_argument( # OK
    '--version',
    version='0.01',
    action='version',)
par.add_argument( # OK
    '-h',
    '--help',
    action='help',
    help='Show this help message and exit',
)

# -m -d


args = par.parse_args()

if args.ge:
    args.patients = 'p*'
    args.exams    = 'e*'
    args.series   = 's*'
    args.files    = '*'

if args.bids:
    args.patients = 'sub-*'
    args.exams    = 'ses-*'
    args.series   = '*'
    args.files    = '*'

if args.write_csv:
    args.table = True

if args.convert:
    args.convert_demo = True

if args.report_files:
    args.keep_sequence_information = True
    args.table                     = True
    args.more_info                 = True

VERBOSE = args.verbose
DEBUG   = args.debug

if VERBOSE > 0:
    print( "Verosity level: " + str(VERBOSE) )

if DEBUG > 0:
    print( "Debug level: " + str(DEBUG) )




levels = [1]
if args.show_patients:
    levels.append(2)
if args.show_exams:
    levels.append(3)
if args.show_series:
    levels.append(4)
if args.show_files:
    levels.append(5)




def lchop(thestring, beginning):
    if thestring.startswith(beginning):
        return thestring[len(beginning):]
    return thestring

def rchop(thestring, ending):
    if thestring.endswith(ending):
        return thestring[:-len(ending)]
    return thestring

def lrchop(thestring, beginning, ending):
    thestring = lchop(thestring, beginning)
    thestring = rchop(thestring, ending)
    return thestring

ansi_escape = re.compile(r'\x1B\[[0-?]*[ -/]*[@-~]')




def print_this(currLog, currLine):
    if not args.quiet:
        print(currLine)
    if args.log:
        with open(currLog, 'a+') as outFile:
            currLine = currLine + '\n'
            outFile.write( ansi_escape.sub('', currLine) )

currLog = args.log_file

if args.log:
    print_this(currLog, 'CYBERCRAFT: Logging to: ' + currLog)




def is_dicom(file_path):
    '''SRC: https://github.com/pydicom/pydicom/blob/f0e8e14ae1019cacaecf7c5f20fee879640f55e9/pydicom/misc.py'''
    with open(file_path, 'rb') as fp:
        fp.read(0x80)  # preamble
        magic = fp.read(4)
    return magic == b"DICM"


def get_patients(tree, patients_pattern):
    for dicom in tree:
        patients = glob.glob( dicom + os.sep + patients_pattern + os.sep )
        tree[dicom] = OrderedDict()
        for patient in patients:
            # patient = patient.rstrip(os.sep).split(os.sep)[-1]
            tree[dicom][patient] = None
    return tree


def get_exams(tree, exams_pattern):
    for dicom in tree:
        for patient in tree[dicom]:
            exams = glob.glob( patient + os.sep + exams_pattern + os.sep )
            tree[dicom][patient] = OrderedDict()
            for exam in exams:
                tree[dicom][patient][exam] = None
    return tree


def get_series(tree, series_pattern):
    for dicom in tree:
        for patient in tree[dicom]:
            for exam in tree[dicom][patient]:
                series = glob.glob( exam + os.sep + series_pattern + os.sep )
                tree[dicom][patient][exam] = OrderedDict()
                for serie in series:
                    tree[dicom][patient][exam][serie] = None
    return tree


def get_dicoms_many(tree, files_patterns):
    for dicom in tree:
        for patient in tree[dicom]:
            for exam in tree[dicom][patient]:
                for serie in tree[dicom][patient][exam]:
                    dicom_files = []
                    for pattern in files_patterns:
                        if DEBUG > 3:
                            print(serie + os.sep + pattern)
                        dicom_files.extend( glob.glob( serie + os.sep + pattern ) )
                    if DEBUG > 3:
                        print(dicom_files)
                    dicom_files = [ fn for fn in dicom_files if is_dicom(fn) ]
                    tree[dicom][patient][exam][serie] = OrderedDict()
                    for dicom_file in dicom_files:
                        tree[dicom][patient][exam][serie][dicom_file] = None
    return tree


def get_dicoms_first(tree, files_patterns):
    for dicom in tree:
        for patient in tree[dicom]:
            for exam in tree[dicom][patient]:
                for serie in tree[dicom][patient][exam]:
                    dicom_files = []
                    for pattern in files_patterns:
                        if DEBUG > 3:
                            print(serie + pattern)
                        dicom_files.extend( glob.glob( serie + pattern ) )
                    if DEBUG > 3:
                        print(dicom_files)
                    tree[dicom][patient][exam][serie] = OrderedDict()
                    for fn in dicom_files:
                        if is_dicom(fn):
                            tree[dicom][patient][exam][serie][fn] = None
                            break
    return tree




def get_dicom_tree_first(
        origin    = os.getcwd(),
        nifti_dir = 'nifti',
        dicom_dir = 'dicom',
        patients  = 'sub-*',
        exams     = 'ses-*',
        series    = 's*',
        files     = ['*'],
):
    dicom_tree_first = { dicom_dir + os.sep: [] }
    dicom_tree_first = get_patients(     dicom_tree_first, patients )
    dicom_tree_first = get_exams(        dicom_tree_first, exams    )
    dicom_tree_first = get_series(       dicom_tree_first, series   )
    dicom_tree_first = get_dicoms_first( dicom_tree_first, files    )
    return dicom_tree_first


def get_dicom_tree_many(
        origin    = os.getcwd(),
        nifti_dir = 'nifti',
        dicom_dir = 'dicom',
        patients  = 'sub-*',
        exams     = 'ses-*',
        series    = 's*',
        files     = ['*'],
):
    dicom_tree_many = { dicom_dir + os.sep: [] }
    dicom_tree_many = get_patients(    dicom_tree_many, patients )
    dicom_tree_many = get_exams(       dicom_tree_many, exams    )
    dicom_tree_many = get_series(      dicom_tree_many, series   )
    dicom_tree_many = get_dicoms_many( dicom_tree_many, files    )
    return dicom_tree_many





def print_info(
        origin    = os.getcwd(),
        nifti_dir = 'nifti',
        dicom_dir = 'dicom',
        patients  = 'sub-*',
        exams     = 'ses-*',
        series    = 's*',
        files     = ['*'],
):
    print_this(currLog, st.RESET_ALL + fg.GREEN + '* INFO'                                                      + st.RESET_ALL)
    print_this(currLog, st.RESET_ALL + fg.BLUE  + '  - O: '    + st.BRIGHT + fg.RED  + '=' + str(origin)    + '=' + st.RESET_ALL)
    print_this(currLog, st.RESET_ALL + fg.BLUE  + '  - N: '    + st.BRIGHT + fg.RED  + '=' + str(nifti_dir) + '=' + st.RESET_ALL)
    print_this(currLog, st.RESET_ALL + fg.BLUE  + '  - D: '    + st.BRIGHT + fg.RED  + '=' + str(dicom_dir) + '=' + st.RESET_ALL)
    print_this(currLog, st.RESET_ALL + fg.BLUE  + '  - P: '    + st.BRIGHT + fg.RED  + '=' + str(patients)  + '=' + st.RESET_ALL)
    print_this(currLog, st.RESET_ALL + fg.BLUE  + '  - E: '    + st.BRIGHT + fg.RED  + '=' + str(exams)     + '=' + st.RESET_ALL)
    print_this(currLog, st.RESET_ALL + fg.BLUE  + '  - S: '    + st.BRIGHT + fg.RED  + '=' + str(series)    + '=' + st.RESET_ALL)
    print_this(currLog, st.RESET_ALL + fg.BLUE  + '  - F: '    + st.BRIGHT + fg.RED  + '=' + str(files)     + '=' + st.RESET_ALL)


def print_tree_standard(tree,lvlCur=1,lvlPrint=range(1,5)):
    lvlMax=max(lvlPrint)
    if isinstance(tree,dict) and lvlCur <= lvlMax:
        for ii, key in enumerate(tree):
            if lvlCur in lvlPrint:
                currLine = st.RESET_ALL + fg.GREEN + '*' * (lvlCur + 1) + ' ' * (lvlMax-lvlCur+1) + fg.BLUE + key + st.RESET_ALL + fg.RED + ' (' + str(ii+1) + ' of ' + str(len(tree)) + ')' + st.RESET_ALL
                print_this(currLog, currLine)
            print_tree_standard(tree[key],lvlCur = lvlCur+1,lvlPrint = lvlPrint)


def print_tree_detailed(dicom_tree):
    indent =  2
    lvlCur = -1
    if isinstance(dicom_tree,dict):
        for dicom in dicom_tree:
            lvlCur += 1
            dicomID = rchop(dicom,os.sep)
            print_this(currLog, st.BRIGHT    + fg.GREEN + ' ' * (lvlCur*indent+2) + '- ' + dicomID + st.RESET_ALL + fg.CYAN + ' [' + dicom + ']' + st.RESET_ALL)
            if isinstance(dicom_tree[dicom],dict):
                for patient in dicom_tree[dicom]:
                    lvlCur += 1
                    patientID = lrchop(patient,dicom+args.patients.rsplit('*',1)[0],os.sep)
                    print_this(currLog, st.BRIGHT    + fg.GREEN + ' ' * (lvlCur*indent+2) + '- ' + patientID + st.RESET_ALL + fg.CYAN + ' [' + patient + ']' + st.RESET_ALL)
                    if isinstance(dicom_tree[dicom][patient],dict):
                        for exam in dicom_tree[dicom][patient]:
                            lvlCur += 1
                            examID = lrchop(exam,patient+args.exams.rsplit('*',1)[0],os.sep)
                            print_this(currLog, st.BRIGHT    + fg.GREEN + ' ' * (lvlCur*indent+2) + '- ' + examID + st.RESET_ALL + fg.CYAN + ' [' + exam + ']' + st.RESET_ALL)
                            if isinstance(dicom_tree[dicom][patient][exam],dict):
                                for serie in dicom_tree[dicom][patient][exam]:
                                    lvlCur += 1
                                    serieID = lrchop(serie,exam+args.series.rsplit('*',1)[0],os.sep)
                                    print_this(currLog, st.BRIGHT    + fg.GREEN + ' ' * (lvlCur*indent+2) + '- ' + serieID + st.RESET_ALL + fg.CYAN + ' [' + serie + ']' + st.RESET_ALL)
                                    if isinstance(dicom_tree[dicom][patient][exam][serie],dict):
                                        for fn in dicom_tree[dicom][patient][exam][serie]:
                                            lvlCur += 1
                                            fnID = lchop(fn,serie)
                                            print_this(currLog, st.BRIGHT    + fg.GREEN + ' ' * (lvlCur*indent+2) + '- ' + fnID + st.RESET_ALL + fg.CYAN + ' [' + fn + ']' + st.RESET_ALL)
                                            lvlCur -= 1
                                    lvlCur -= 1
                            lvlCur -= 1
                    lvlCur -= 1
            lvlCur -= 1




def convert_dicom_data(dicom_tree,convert_demo=True,convert=False):
    if isinstance(dicom_tree,dict):
        for ii, dicom in enumerate(dicom_tree):
            dicomID = dicom.strip(os.sep)
            currLine = st.RESET_ALL + fg.GREEN + '**' * 1 + ' ' + fg.BLUE + dicomID + st.RESET_ALL + fg.RED + ' (' + str(ii+1) + ' of ' + str(len(dicom_tree)) + ')' + st.RESET_ALL
            print_this(currLog, currLine)

            if isinstance(dicom_tree[dicom],dict):
                for jj, patient in enumerate(dicom_tree[dicom]):
                    patientID = lrchop(patient,dicom+args.patients.rsplit('*',1)[0],os.sep)
                    currLine = st.RESET_ALL + fg.GREEN + '***' * 1 + ' ' + fg.BLUE + patientID + st.RESET_ALL + fg.RED + ' (' + str(jj+1) + ' of ' + str(len(dicom_tree[dicom])) + ')' + st.RESET_ALL
                    print_this(currLog, currLine)

                    if isinstance(dicom_tree[dicom][patient],dict):
                        for kk, exam in enumerate(dicom_tree[dicom][patient]):
                            examID = lrchop(exam,patient+args.exams.rsplit('*',1)[0],os.sep)
                            currLine = st.RESET_ALL + fg.GREEN + '****' * 1 + ' ' + fg.BLUE + examID + st.RESET_ALL + fg.RED + ' (' + str(kk+1) + ' of ' + str(len(dicom_tree[dicom][patient])) + ')' + st.RESET_ALL
                            print_this(currLog, currLine)
                            args_list = []
                            if False:
                                args_list.extend([ 'echo',                         ]) # Executable (for testing)
                            if True:
                                args_list.extend([ 'dcm2bids',                     ]) # Executable
                                args_list.extend([ '-o', args.nifti_dir,           ]) # Output BIDS directory
                                args_list.extend([ '-p', patientID,                ]) # Participant ID
                                args_list.extend([ '-d', exam,                     ]) # DICOM directory (data source)
                                args_list.extend([ '-s', examID,                   ]) # Session ID
                            if args.study_specific_json is None:
                                args_list.extend([ '-c', 'ses-' + examID + '.json' ]) # JSON configuration file (see example/config.json)
                            else:
                                args_list.extend([ '-c', args.study_specific_json, ]) # JSON configuration file (see example/config.json)
                            if not args.no_defaceing:
                                args_list.extend([ '-a', 'cc-deface-args.sh',      ]) # Link to anonymizer script/binary
                            if False:
                                args_list.extend([ '--clobber',                    ]) # Overwrite old temporary dcm2niix output if it exists
                                args_list.extend([ '--forceDcm2niix',              ]) # Overwrite old temporary dcm2niix output if it exists
                            if args.convert_demo:
                                print_this(currLog, '     ' + ' '.join(args_list))
                            if args.convert:
                                subprocess.call(args_list)




def table_dicom_data(dicom_tree):
    if VERBOSE > 1:
        print_this(currLog, st.RESET_ALL + fg.GREEN + '* VISITING' + st.RESET_ALL)

    new_dict_items = {
        0x0021105e: ('DS', '1', "FloatSlopRTIATimer",   '', 'FloatSlopRTIATimer'),
        0x00181060: ('DS', '1', "TriggerTime",          '', 'TriggerTime'),
    }
    DicomDictionary.update(new_dict_items)
    new_names_dict = dict([(val[4], tag) for tag, val in new_dict_items.items()])
    keyword_dict.update(new_names_dict)

    props = [
        'Path',
        'NumFiles',
        'StudyID',
        'StudyDate',
        'StudyTime',
        'PatientID',
        'PatientName',
        'ProtocolName',
        'SeriesDescription',
        'SeriesNumber',
        'SeriesTime',
        'ImagesInAcquisition',
        'InStackPositionNumber',
        'InstanceNumber',
        'SliceLocation',
        'TriggerTime',
        'FloatSlopRTIATimer',
    ]

    db = OrderedDict([(ii, []) for ii in props ])
    if isinstance(dicom_tree,dict):
        for dicom in dicom_tree:
            if isinstance(dicom_tree[dicom],dict):
                for patient in dicom_tree[dicom]:
                    if isinstance(dicom_tree[dicom][patient],dict):
                        for exam in dicom_tree[dicom][patient]:
                            if isinstance(dicom_tree[dicom][patient][exam],dict):
                                for serie in dicom_tree[dicom][patient][exam]:
                                    if isinstance(dicom_tree[dicom][patient][exam][serie],dict):
                                        for fn in dicom_tree[dicom][patient][exam][serie]:
                                            if VERBOSE > 1:
                                                print_this(currLog, fg.CYAN + '  ' + fn + st.RESET_ALL )

                                            ds = pydicom.dcmread(fn, stop_before_pixels = True)
                                            ds.decode()
                                            for prop in props:
                                                if prop == 'Path':
                                                    db[prop] += [ os.path.dirname(fn) ]
                                                elif prop == 'NumFiles':
                                                    db[prop] += [  len( glob.glob( os.path.dirname(fn)  + '/*' ) )  ]
                                                else:
                                                    db[prop] += [''] if not hasattr(ds, prop) else [getattr(ds, prop) ]

                                            if False:
                                                print_this( currLog, str(db) )
                                            if False:
                                                print_this( currLog, str(ds) )

    df = pd.DataFrame(db)
    df = df.sort_values(by=['StudyID','SeriesNumber','InstanceNumber'],axis=0)
    if not args.include_stacked_screen_saves:
        df = df[df.SeriesDescription.str.contains("Stacked Screen Save") == False]
        df = df[df.SeriesDescription.str.contains("Screen Save") == False]
    if not args.keep_sequence_information:
        # Trim path.
        df.Path = df.Path.apply(lambda xx: os.path.dirname(xx) )
        # Remove duplicates (must be done AFTER path trimming).
        df = df.drop_duplicates(subset='Path', keep="last").reset_index()
    if not args.more_info:
        # df = df.drop(labels='index',                 axis='columns')
        df = df.drop(labels='SeriesDescription',     axis='columns')
        df = df.drop(labels='SeriesNumber',          axis='columns')
        df = df.drop(labels='SeriesTime',            axis='columns')
        df = df.drop(labels='NumFiles',              axis='columns')
    if not args.report_files:
        df = df.drop(labels='ImagesInAcquisition',   axis='columns')
        df = df.drop(labels='InStackPositionNumber', axis='columns')
        df = df.drop(labels='InstanceNumber',        axis='columns')
        df = df.drop(labels='SliceLocation',         axis='columns')
        df = df.drop(labels='FloatSlopRTIATimer',    axis='columns')
        df = df.drop(labels='TriggerTime',           axis='columns')

    print_this(currLog, st.RESET_ALL + fg.GREEN + '* TABLE' + st.RESET_ALL)
    print_this(currLog, tabulate(df, headers='keys', tablefmt="orgtbl") )
    if args.write_csv:
        fn = args.log_file
        df.to_csv( path_or_buf = rchop(fn,'.org') + '.csv', sep = ',')
        # else:
        #     print( fg.GREEN + '* No DICOM files have been found!' + st.RESET_ALL )





if __name__ == "__main__":

    if len(sys.argv[1:])==0:
        # parser.print_help()
        par.print_usage()
        print(par.epilog)

    if DEBUG > 0:
        for arg in vars(args):
                print(fg.GREEN + st.BRIGHT + arg + st.RESET_ALL + '\n    ' + str(getattr(args, arg)))


    print_info(
        origin    = args.origin,
        nifti_dir = args.nifti_dir,
        dicom_dir = args.dicom_dir,
        patients  = args.patients,
        exams     = args.exams,
        series    = args.series,
        files     = args.files,
    )
    if args.report_files:
        dicom_tree = get_dicom_tree_many(
            origin    = args.origin,
            nifti_dir = args.nifti_dir,
            dicom_dir = args.dicom_dir,
            patients  = args.patients,
            exams     = args.exams,
            series    = args.series,
            files     = args.files,
        )
    else:
        dicom_tree = get_dicom_tree_first(
            origin    = args.origin,
            nifti_dir = args.nifti_dir,
            dicom_dir = args.dicom_dir,
            patients  = args.patients,
            exams     = args.exams,
            series    = args.series,
            files     = args.files,
        )
    if VERBOSE > 3:
        print(dicom_tree)



    if len(levels) > 1:
        print_this(currLog, st.RESET_ALL + fg.GREEN + '* STANDARD TREE' + st.RESET_ALL)
        print_tree_standard( dicom_tree, lvlCur=1, lvlPrint=levels )
    if VERBOSE > 0:
        print_this(currLog, st.RESET_ALL + fg.GREEN + '* DETAILED TREE' + st.RESET_ALL)
        print_tree_detailed(dicom_tree)
    if args.convert_demo:
        print_this(currLog, st.RESET_ALL + fg.GREEN + '* CONVERT' + st.RESET_ALL)
        convert_dicom_data(dicom_tree,convert_demo=args.convert_demo,convert=args.convert)
    if args.table:
        table_dicom_data(dicom_tree)
