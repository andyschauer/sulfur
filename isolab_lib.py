#!/usr/bin/env python3
"""
Library of functions generally useful while processing data from IsoLab.
"""

__author__ = "Andy Schauer"
__email__ = "aschauer@uw.edu"
__last_modified__ = "2024-04-16"
__version__ = "2.0"
__copyright__ = "Copyright 2024, Andy Schauer"
__license__ = "Apache 2.0"



# ---------- imports ----------
import json
import numpy as np
import os
import re



# ---------- functions ----------
def get_outliers(data, sigma):
    m = np.nanmean(data)
    s = np.nanstd(data) * sigma
    o = [i for i, e in enumerate(data) if e > m + s or e < m - s]
    return m, s, o



def get_path(instrument, desired_path):
    """Make your life easier with this section. These are the paths that seem to change depending on the computer we are working on."""
    path_file = os.path.join(os.getcwd(), f'{instrument}_path.txt')
    if os.path.isfile(path_file):
        # print(' :-) Using existing shrekCN path file for a warm and fuzzy experience. (-:')
        with open(path_file, 'r') as ppf:
            home, python_path, project_path, standards_path = ppf.readline().split(',')
            python_path = home + python_path
            project_path = home + project_path
            standards_path = home + standards_path

    else:
        python_path_check = False
        project_path_check = False
        standards_path_check = False
        print(' )-: Happy path file does not exist yet. :-(')
        print(" Let's make one... :-| ")
        while python_path_check is False:
            python_path = input(f'Enter the current path to the python scripts. Perhaps it is {os.getcwd()}. ')
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

        with open(path_file, 'w') as ppf:
            ppf.write(f'{python_path},{project_path},{standards_path}')

    if desired_path == "project":
        return project_path
    elif desired_path == "python":
        return python_path
    elif desired_path == "standards":
        return standards_path
    elif desired_path == "refmat_meta":
        with open(standards_path) as file:
            refmat = json.load(file)
            return f"{refmat['file_meta_data']['file']} - {refmat['file_meta_data']['modification_date']}"
    else:
        unknown_path = input('Unknown path, please enter the path to your project: ')
        return unknown_path



def make_file_list(directory, filetype):
    """Create and return a list of files contained within a directory
    of file type."""
    filelist = []
    initial_list = os.listdir(directory)
    for file in initial_list:
        if re.search(filetype, file):
            filelist.append(file)
    return filelist



def read_file(file_to_import, delim=None, header_row=1):
    """Read in a delimited text file containing a single header row
    followed by data and return those headers as a list and the data
    as a dictionary."""
    with open(file_to_import, 'r') as f:
        if header_row > 1:
            f.readline()
        headers = f.readline().split(delim)

        # remove unwanted characters from headers using a regular expression
        p = re.compile(r'[./\s()%]')  # list of characters to match
        for ii in range(len(headers)):
            m = p.findall(headers[ii])
            for em in m:
                headers[ii] = headers[ii].replace(em, '')

        data = {}
        for h in headers:
            data[h] = []

        for line in f:
            row = line.split(delim)
            if delim is None:
                if len(row) < len(headers):
                    row.append(0)
            # populate dictionary with all data in all rows
            for h, v in zip(headers, row):
                v = v.replace('\n', '').replace('1.#IO', '').replace('1.#INF000', '')
                if v == '':
                    v = None
                data[h].append(v)

    return headers, data
