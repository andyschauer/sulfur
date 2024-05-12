#!/usr/bin/env python3
"""
This script opens raw sulfur data files from the mass spectrometer Shrek, organizes based on peaks (e.g.
    working gas vs sample), and writes the reduced data to shrekS_analysis_log.csv. Initial quality control
    is completed (e.g. saturated detectors). If analyses are deemed erroneous, a flag is set to zero, but the
    data are still written to the log file.

    version 1.x - 2019.03.01 - matlab versions used while running shrek in the way Julien Foriel originally set it up (50 ug of Sulfur target).
        The last matlab script under this scenario is saved as shrekS_Archive_190311.m. Note that these versions were a single script (both shrekS
        and shrekScalibrate).
    version 2.x - 2020.08.17 - matlab versions using the new DIY cryofocusing technique. Code is now separated into a mass spec file reader (shrekS)
        and a log file reader / sample calibrator (shrekScalibrate). This eventually migrated to python late in Urusula's project and during Alli's
        project.
    version 3.0 - 2023.12.04 - GasBench versions start here. Stopped using shrekS_standards in favor of lab wide reference_materials.json.
    version 4.0 - 2024.03.11 - changed to the "get path" model for easier use on multiple computers
    version 4.1 - 2024.05.11 - changed lab import to isolab_lib and updated associated functions from older version of lab to current version of isolab_lib
"""

__authors__ = "Andy Schauer, Ursula Jongebloed"
__email__ = "aschauer@uw.edu"
__last_modified__ = "2024.05.11"
__version__ = "4.1"
__copyright__ = "Copyright 2024, Andy Schauer"
__license__ = "Apache 2.0"
__acknowledgements__ = "Alli Moon, Drew Pronovost"



# -------------------- imports --------------------
import csv
import json
import isolab_lib
import os
import re
from shrekS_lib import *
import sys
import time


# -------------------- functions --------------------
def add_to_data_meta():
    for i in data_meta:
        data_meta[i].append(data[i][index]) if i in data else data_meta[i].append(None)


def add_to_data_supp():
    data_supp['file'].append(file)
    data_supp['flag'].append(flag)
    data_supp['notes'].append(note)
    data_supp['rows_per_sample'].append(rows)
    data_supp['pyversions'].append(version)
    data_supp['empty'].append('')

    if data['Comment'][index] == 'amount_from_balance':
        get_Sqty(data['Identifier1'][index], float(data['Amount'][index]))
    else:
        data_supp['Sqty'].append(data['Amount'][index])

    if 'Information' in data:
        if data['Information'][index] is None:
            data_supp['peak_center'].append(None)
        elif re.match('Peak Center found at', data['Information'][index]) is None:
            data_supp['peak_center'].append(None)
        else:
            m = re.findall(r'\d+', data['Information'][index])
            data_supp['peak_center'].append(m[0])
    else:
        data_supp['peak_center'].append(None)


def add_to_data_wg():  # put data for sulfur refence peak into list
    for i in data_wg:
        data_wg[i].append(data[i][index + peak_number_offset]) if i in data else data_wg[i].append(None)


def add_to_data_sam():  # put data for sulfur sample peak into list
    for i in data_sam:
        data_sam[i].append(data[i][index + peak_number_offset]) if i in data else data_sam[i].append(None)


def none_wg():  # set all sulfur reference peak data to None
    for i in data_wg:
        data_wg[i].append(None)


def none_sam():  # set all sulfur sample peak data to None
    for i in data_sam:
        data_sam[i].append(None)


def get_Sqty(material, mass):
    """Calculate the amount of sulfur in a given mass based on an accepted percent sulfur. Percent sulfur precalculated and hardcoded into
        reference_materials.json."""
    if material in refmat_list:
        data_supp['Sqty'].append(mass * eval(material)['fractionS'] * 1000)
    else:
        data_supp['Sqty'].append(None)


# -------------------- setup --------------------
version = os.path.basename(__file__) + ' - ' + time.ctime(os.path.getctime(__file__))

python_directory = isolab_lib.get_path("shrekS", "python")
project_directory = isolab_lib.get_path("shrekS", "project")
new_data_directory = 'rawdata_new'
archive_data_directory = 'rawdata_archive'
junk_data_directory = 'rawdata_junk'
log_file_name = 'shrekS_analysis_log.csv'

if os.path.isdir(project_directory) is False:
    print('directory does not exist...exiting....')
    sys.exit()


# -------------------- reference materials ----------------------------
# load reference material information
with open(isolab_lib.get_path("shrekS", "standards"), 'r') as f:
    refmat = json.load(f)

refmat_list = []
refmat_keys = refmat['sulfates'].keys()
for i in refmat_keys:
    globals()[i] = refmat['sulfates'][i]
    refmat_list.append(i)

refmat_keys = refmat['sulfides'].keys()
for i in refmat_keys:
    globals()[i] = refmat['sulfides'][i]
    refmat_list.append(i)


# -------------------- create list of files to process --------------------
filelist = isolab_lib.make_file_list(os.path.join(project_directory, new_data_directory), 'csv')
if not filelist:
    print(f'\n    No files in {new_data_directory}.')
else:
    filelist = sorted(filelist)
    print(f'\n    Processing {len(filelist)} mass spec files.')


# -------------------- Main loop through each file --------------------
for file in filelist:

    for i in headers_meta:
        data_meta[i] = []

    for i in headers_peak:
        data_wg[i] = []
        data_sam[i] = []

    for i in headers_supp:
        data_supp[i] = []

    headers, data = isolab_lib.read_file(os.path.join(project_directory, new_data_directory, file), ',')

    # test for file problems
    if 'Analysis' not in headers:
        print(f'    Problem with file {file}. Perhaps only headers and no data in file.')
        os.rename(os.path.join(project_directory, new_data_directory, file), os.path.join(project_directory, junk_data_directory, file))
        print('file ' + file + ' was moved to the junk folder')
        continue

    for n in data['Analysis']:
        try:
            data['Analysis'] = [int(index) for index in data['Analysis']]
        except ValueError:
            print(f'File {file} contains strings in the Analysis column')
            os.rename(os.path.join(project_directory, new_data_directory, file), os.path.join(project_directory, junk_data_directory, file))
            print(f'File {file} was moved to the junk folder.')
            continue

    # sample_index_first_row = first row of each sample
    analysis_difference = [j - index for index, j in zip(data['Analysis'][: -1], data['Analysis'][1:])]
    sample_index_first_row = [index for index, e in enumerate(analysis_difference) if e != 0]
    sample_index_first_row = [index + 1 for index in sample_index_first_row]
    sample_index_first_row.append(0)
    sample_index_first_row = sorted(sample_index_first_row)

    # rows_per_sample
    rows_per_sample = [j - index for index, j in zip(sample_index_first_row[:-1], sample_index_first_row[1:])]
    rows_per_sample.append(len(data['Analysis']) - sample_index_first_row[len(sample_index_first_row) - 1])  # adds last sample number of rows

    # sample_index_last_row = last row of each sample
    sample_index_last_row = [j + index for index, j in zip(sample_index_first_row, rows_per_sample)]


    # Loop through each sample
    for index, rows in zip(sample_index_first_row, rows_per_sample):
        flag = 1
        note = ''
        add_to_data_meta()

        if rows == 1:  # One reference peak and no sample peak
            currnote = 'Only one peak detected; '
            print(f"{currnote} - file {file} and Analysis {str(data['Analysis'][index])}")
            note = note + currnote
            peak_number_offset = 0
            add_to_data_wg()
            none_sam()

        elif rows == 2:  # One reference peak and one sample peak
            peak_number_offset = 0
            add_to_data_wg()
            peak_number_offset = 1
            add_to_data_sam()

        elif rows == 3:  # One reference peak, one organic related peak, and one sample peak
            currnote = 'Third peak detected; '
            print(f"{currnote} - file {file} and Analysis {str(data['Analysis'][index])}")
            note = note + currnote
            peak_number_offset = 0
            add_to_data_wg()
            peak_number_offset = 2
            add_to_data_sam()

        elif rows > 3:
            flag = 0
            currnote = f'More than 3 peaks detected (n = {str(rows)}); check trace; '
            note = note + currnote
            print(f"{currnote} - file {file} and Analysis {str(data['Analysis'][index])}")
            none_wg()
            none_sam()

        else:
            flag = 0
            currnote = 'Problem counting peaks; '
            note = note + currnote
            print(f"{currnote} - file {file} and Analysis {str(data['Analysis'][index])}")
            none_wg()
            none_sam()

        if int(data['Ampl64'][index]) > 49900 or int(data['Ampl66'][index]) > 49900:
            flag = 0
            currnote = 'Saturdated detector; '
            note = note + currnote
            print(f"{currnote} - file {file} and Analysis {str(data['Analysis'][index])}")

        add_to_data_supp()



    # export new raw data to analysis log file
    if os.path.isfile(os.path.join(project_directory, log_file_name)) is False:
        with open(os.path.join(project_directory, log_file_name), 'w', newline='') as csvfile:  # if the log file has not been created, create it with column headers and data
            datawriter = csv.writer(csvfile, quoting=csv.QUOTE_MINIMAL)
            datawriter.writerow(shrekS_analysis_log_headers)
            for ii in range(len(data_meta['Analysis'])):
                datawriter.writerow(eval(data_to_write))

    else:
        with open(os.path.join(project_directory, log_file_name), 'a', newline='') as csvfile:
            datawriter = csv.writer(csvfile, quoting=csv.QUOTE_MINIMAL)
            for ii in range(len(data_meta['Analysis'])):
                datawriter.writerow(eval(data_to_write))

    os.rename(os.path.join(project_directory, new_data_directory, file), os.path.join(project_directory, archive_data_directory, file))  # done with datafile, put it in the archive directory


print('\n    Done')
print(f'\n    Sulfur data written to {log_file_name}')
print('\n    Now you may wish to run shrekScalibrate.py')
