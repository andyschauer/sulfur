#!/usr/bin/env python3
"""
Library of functions used by the IsoLab shrekS_* suite of python scripts.

    version 3.0 - 2023.12.04 - Starting off with version 3 to match the other shrekS scripts. This library did not exist
        prior to this version.
    version 4.0 - 2024.03.11 - added get_path function
"""


__authors__ = "Andy Schauer, Ursula Jongebloed"
__email__ = "aschauer@uw.edu"
__last_modified__ = "2024-03-11"
__version__ = "4.0"
__copyright__ = "Copyright 2024, Andy Schauer"
__license__ = "Apache 2.0"
__acknowledgements__ = "Alli Moon, Drew Pronovost"



# -------------------- imports --------------------
import os



# -------------------- functions --------------------
def get_path(desired_path):
    """Make your life easier with this section. These are the paths that seem to change depending on the computer we are working on."""
    shrekS_path_file = os.path.join(os.getcwd(), 'shrekS_path.txt')
    if os.path.isfile(shrekS_path_file):
        # print(' :-) Using existing shrekS path file for a warm and fuzzy experience. (-:')
        with open(shrekS_path_file, 'r') as ppf:
            python_path, project_path, standards_path = ppf.readline().split(',')

    else:
        python_path_check = False
        project_path_check = False
        standards_path_check = False
        print(' )-: shrekS path file does not exist yet. :-(')
        print(" Let's make one... :-| ")
        while python_path_check is False:
            python_path = input(f'Enter the current path to the shrekS python scripts. Perhaps it is {os.getcwd()}. ')
            if os.path.isdir(python_path):
                python_path_check = True
                if python_path[-1] != '/':
                    python_path += '/'
            else:
                print(f'oops, try typing that in again (you typed {python_path}): ')

        while project_path_check is False:
            project_path = input('Enter the current path to your projects: ')
            if os.path.isdir(project_path):
                project_path_check = True
                if project_path[-1] != '/':
                    project_path += '/'
            else:
                print(f'oops, try typing that in again (you typed {project_path}): ')

        while standards_path_check is False:
            standards_path = input('Enter the current path and filename to your reference materials file: ')
            if os.path.isfile(standards_path):
                standards_path_check = True
            else:
                print(f'oops, try typing that in again (you typed {standards_path}): ')

        with open(shrekS_path_file, 'w') as ppf:
            ppf.write(f'{python_path},{project_path},{standards_path}')

    if desired_path == "project":
        return project_path
    elif desired_path == "python":
        return python_path
    elif desired_path == "standards":
        return standards_path
    else:
        unknown_path = input('Enter the path to your project: ')
        return unknown_path



zero = {'names': ['zero'],
        'material': 'reference gas peaks treated as unknowns'}

blank = {'names': ['blank'],
             'material': 'a tin capsule packed with WO3 and Sn powder'}

qtycal = {'names': ['qtycal_Na2SO4', 'qtycal'],
          'material': 'Na2SO4',
          'fractionS': 0.22574,
          'notes': 'sodium sulfate solution based amounts'}

void = {'names': ['void'],
         'material': None,
         'notes': 'no material dropped into EA'}



headers_meta = ['Identifier1', 'Analysis', 'Amount', 'Date', 'Time', 'Comment', 'Identifier2', 'Preparation', 'Method']

headers_peak = ['Ampl64', 'Ampl66', 'AreaAll', 'Area64', 'Area66', 'BGD64', 'BGD66', 'R34S32S', 'R66SO264SO2', 'd34S32S',
             'PeakNr', 'Start', 'Width', 'End']

headers_supp = ['Sqty', 'file', 'flag', 'notes', 'peak_center', 'rows_per_sample', 'pyversions', 'empty']


numlist = ['Analysis', 'Amount', 'Ampl64_wg', 'Ampl66_wg', 'AreaAll_wg', 'Area64_wg', 'Area66_wg',
           'BGD64_wg', 'BGD66_wg', 'R34S32S_wg', 'R66SO264SO2_wg', 'd34S32S_wg', 'PeakNr_wg', 'Start_wg',
           'Width_wg', 'End_wg', 'Ampl64_sam', 'Ampl66_sam', 'AreaAll_sam', 'Area64_sam', 'Area66_sam',
           'BGD64_sam', 'BGD66_sam', 'R34S32S_sam', 'R66SO264SO2_sam', 'd34S32S_sam', 'PeakNr_sam',
           'Start_sam', 'Width_sam', 'End_sam', 'Sqty', 'flag', 'peak_center', 'rows_per_sample']

data_meta = {}
data_wg = {}
data_sam = {}
data_supp = {}

shrekS_analysis_log_headers = []
shrekS_analysis_log_headers.extend(headers_meta)
shrekS_analysis_log_headers.extend([f"{i}_wg" for i in headers_peak])
shrekS_analysis_log_headers.extend([f"{i}_sam" for i in headers_peak])
shrekS_analysis_log_headers.extend(headers_supp)

data_to_write = []
data_to_write.extend([f"data_meta['{i}'][ii]" for i in headers_meta])
data_to_write.extend([f"data_wg['{i}'][ii]" for i in headers_peak])
data_to_write.extend([f"data_sam['{i}'][ii]" for i in headers_peak])
data_to_write.extend([f"data_supp['{i}'][ii]" for i in headers_supp])
data_to_write = str(data_to_write).replace("\"", "")


python_scripts = ['shrekS_lib.py', 'shrekS.py', 'shrekScalibrate.py']
