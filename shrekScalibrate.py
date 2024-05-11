#!/usr/bin/env python3
"""
This script opens the shrekS_analysis_log.csv for sample diagnosis and calibration.
    Figures are used to help assess run quality. Ultimately, 34S isotope composition
    delta values are normalized to the VCDT scale and exported along with sulfur
    quantity and other salient bits.

    version 1.x - 2019.03.01 - matlab versions used while running shrek in the way Julien Foriel originally set it up (50 ug of Sulfur target).
        The last matlab script under this scenario is saved as shrekS_Archive_190311.m. Note that these versions were a single script (both shrekS
        and shrekScalibrate).
    version 2.x - 2020.08.17 - matlab versions using the new DIY cryofocusing technique. Code is now separated into a mass spec file reader (shrekS)
        and a log file reader / sample calibrator (shrekScalibrate). This eventually migrated to python late in Urusula's project and during Alli's
        project. Ursula developed a separate python script while andy continued to limp along with matlab. Andy eventually incorporated much of Urusula's
        script during Alli's project.
    version 3.0 - 2023.12.04 - GasBench versions start here. Stopped using shrekS_standards in favor of lab wide reference_materials.json.
    version 4.0 - 2024.03.11 - changed to the "get path" model for easier use on multiple computers 
"""

__authors__ = "Andy Schauer, Ursula Jongebloed"
__email__ = "aschauer@uw.edu"
__last_modified__ = "2024.03.11"
__version__ = "4.0"
__copyright__ = "Copyright 2024, Andy Schauer"
__license__ = "Apache 2.0"
__acknowledgements__ = "Alli Moon, Drew Pronovost"



# -------------------- imports --------------------
import argparse
from bokeh.plotting import figure, show
from bokeh.embed import file_html
from bokeh import palettes
from bokeh.resources import CDN  # , INLINE
import csv
import datetime as dt
import json
import lab
import matplotlib.pyplot as pplt
from natsort import natsorted
import numpy as np
import os
from shrekS_lib import *
import shutil
import sys
import time
import webbrowser


# -------------------- functions --------------------
def rewrite_log_file():
    """Recreate the analysis log file where the only modification is likely to be setting a flag to 0, making a note, and updating the python script version.
       I have chosen this way while I am still using csv files to archive data. JSON or other dictionary based storage would make this much easier in here
       but less immediately transferable to other platforms compared with csv.

       NOTE - The empty column is a hack because I could not get jupyter or python on linux to stop throwing a carriage return within the version string."""

    # Make a copy of the existing log file
    curr_time = dt.datetime.now()
    archive_log_file_name = f"{log_file_name[0:-4]}_Archive_{curr_time.strftime('%Y%m%d%H%M%S')}.csv"
    src = os.path.join(home_dir, project_directory, log_file_name)
    dst = os.path.join(home_dir, project_directory, 'archive/', archive_log_file_name)
    shutil.copy2(src, dst)
    os.remove(os.path.join(home_dir, project_directory, log_file_name))

    # use the existing accepted log file headers (defined in shrekS_standards.py) to construct the export data
    data_to_write = [f"original_data['{i}'][ii]" for i in log_file_headers]
    data_to_write = str(data_to_write).replace('"', '')

    # write data to csv log file
    with open(os.path.join(home_dir, project_directory, log_file_name), 'w', newline='') as csvfile:
        datawriter = csv.writer(csvfile, quoting=csv.QUOTE_MINIMAL)
        datawriter.writerow(log_file_headers)
        for ii in range(len(original_data['Analysis'])):
            datawriter.writerow(eval(data_to_write))


def toggle_flag_to_zero(bad_data_Analysis):
    """So far this function is used manually. Call it after identifying data that should be excluded from further consideration for reasons
       specified in the note. After calling this function n times, once per serial number (Analysis), then call rewrite_log_file() and the flag
       will be changed to zero, the note will be inserted, and the python script version will be updated.

       Pass the bad sample analysis number(s) to this function one at a time, enter the reason when prompted, then run rewrite_log_file().
    """
    bodi = original_data['Analysis'].index(str(bad_data_Analysis))
    note = input('Why are you excluding this sample: ')
    original_data['flag'][bodi] = 0
    original_data['notes'][bodi] = note
    original_data['pyversions'][bodi] = version


def get_outliers(data, sigma):
    m = np.nanmean(data)
    s = np.nanstd(data) * sigma
    o = [i for i, e in enumerate(data) if e > m + s or e < m - s]
    return m, s, o


# -------------------- setup --------------------

# verbose mode
parser = argparse.ArgumentParser()
parser.add_argument("--verbose", help="Include exhaustive diagnostic information and figures in report.", action="store_true")
args = parser.parse_args()
if args.verbose:
    verbose = True
    print('\n    ****    Verbose mode on    ****    ')
else:
    verbose = False

# start_time = time.time()
version = os.path.basename(__file__) + ' - ' + time.ctime(os.path.getctime(__file__))

python_directory = get_path("python")
project_directory = get_path("project")
new_data_directory = 'rawdata_new'
archive_data_directory = 'rawdata_archive'
junk_data_directory = 'rawdata_junk'
log_file_name = 'shrekS_analysis_log.csv'
report_dir = 'report/'
fig_dir = 'figures/'

log_file_name = 'shrekS_analysis_log.csv'
# log_file_name = 'shrekS_analysis_log_fall2023.csv'

if os.path.isdir(project_directory) is False:
    print('\n    The project directory does not exist...exiting....')
    sys.exit()

# copy log file to report directory
shutil.copy2(os.path.join(project_directory, log_file_name), os.path.join(project_directory, report_dir, f"data/{log_file_name}_REPORT_COPY"))

# load reference material information
with open(get_path("standards"), 'r') as f:
    refmat = json.load(f)

refmat_keys = refmat['sulfates'].keys()
for i in refmat_keys:
    globals()[i] = refmat['sulfates'][i]
    globals()[i]['index'] = np.empty(0, dtype="int16")

refmat_keys = refmat['sulfides'].keys()
for i in refmat_keys:
    globals()[i] = refmat['sulfides'][i]
    globals()[i]['index'] = np.empty(0, dtype="int16")


# -------------------- get data --------------------
print('    Importing data...')
headers, data = lab.read_file(os.path.join(project_directory, log_file_name), ',')


# -------------------- organize data --------------------
print('    Organizing data...')

data['empty'] = ['' for i in data['Analysis']]  # currently a carriage return is probably written into the "empty" column, which I don't want

# remove flag 0 analyses from data
original_data = data.copy()
flag0_indices = [i for i, e in enumerate(data['flag']) if int(e) == 0]
flag1_indices = [i for i, e in enumerate(data['flag']) if int(e) == 1]
for header in headers[:-1]:
    data[header] = [data[header][index] for index in flag1_indices]

strlist = set(headers) - set(numlist)

for i in numlist:
    globals()[i] = np.asarray(data[i], dtype=float)

for i in strlist:
    globals()[i] = np.asarray(data[i])

run_index_first_row = [np.where(file == j)[0][0] for j in sorted(set(file))]

# Identifier1 = np.asarray([i.strip() for i in Identifier1])


# -------------------- categorize all Identifier 1 (Identifier1) analyses into appropriate material type --------------------
all_analyses = [i for i in range(len(Analysis))]

void['index'] = [i for i, e in enumerate(Identifier1) if str(e).lower() in (name.lower() for name in void['names'])]
blank['index'] = [i for i, e in enumerate(Identifier1) if str(e).lower() in (name.lower() for name in blank['names'])]
zero['index'] = [i for i, e in enumerate(Identifier1) if str(e).lower() in (name.lower() for name in zero['names'])]

Ag2S['index'] = [i for i, e in enumerate(Identifier1) if str(e).lower() in (name.lower() for name in Ag2S['names'])]
BaSO4['index'] = [i for i, e in enumerate(Identifier1) if str(e).lower() in (name.lower() for name in BaSO4['names'])]
Na2SO4['index'] = [i for i, e in enumerate(Identifier1) if str(e).lower() in (name.lower() for name in Na2SO4['names'])]
ZnS['index'] = [i for i, e in enumerate(Identifier1) if str(e).lower() in (name.lower() for name in ZnS['names'])]
Ag2SO4['index'] = [i for i, e in enumerate(Identifier1) if str(e).lower() in (name.lower() for name in Ag2SO4['names'])]

SodSul_1 = {'index': [i for i, e in enumerate(Identifier1) if e == 'SodSul_1']}
SodSul_1_RotoVap = {'index': [i for i, e in enumerate(Identifier1) if e == 'SodSul_1_RotoVap']}
SodSul_1_Oven60 = {'index': [i for i, e in enumerate(Identifier1) if e == 'SodSul_1_Oven60']}
SodSul_2 = {'index': [i for i, e in enumerate(Identifier1) if e == 'SodSul_2']}
SodSul_3 = {'index': [i for i, e in enumerate(Identifier1) if e == 'SodSul_3']}
SodSul_4 = {'index': [i for i, e in enumerate(Identifier1) if e == 'SodSul_4']}
SodSul_5 = {'index': [i for i, e in enumerate(Identifier1) if e == 'SodSul_5']}
SodSul_6 = {'index': [i for i, e in enumerate(Identifier1) if e == 'SodSul_6']}
SodSul_7 = {'index': [i for i, e in enumerate(Identifier1) if e == 'SodSul_7']}
SodSul_8 = {'index': [i for i, e in enumerate(Identifier1) if e == 'SodSul_8']}


# qtycal['index'] = [i for i, e in enumerate(Identifier1) if str(e).lower() in (name.lower() for name in qtycal['names'])]
# qtycal['index'] = SodSul_1['index']
# qtycal['index'] = SodSul_1_RotoVap['index']
qtycal['index'] = SodSul_1_Oven60['index']

all_stds = Ag2S['index'] + BaSO4['index'] + Na2SO4['index'] + ZnS['index'] + qtycal['index']

samples = {'index': list(set(all_analyses) - set(Ag2S['index']) - set(BaSO4['index']) - set(Na2SO4['index']) - set(ZnS['index']) - set(void['index']) - set(blank['index']) - set(qtycal['index']) - set(zero['index']))}
# samples = {'index': []}
# samples['index'] = np.intersect1d(samples['index'], sample_index_project)

recent_run = np.arange(run_index_first_row[-1], len(Analysis))


# get WO3_Sn amount from Identifier2
WO3_Sn = np.zeros(len(Identifier2))
for i,j in enumerate(Identifier2):
    if j is not None:
        try:
            if j[0:6]=='WO3_Sn':
                WO3_Sn[i] = j[7:]
        except ValueError:
            WO3_Sn[i] = 0



# -------------------- Calculations --------------------
print('    Maths...')
# calculate how long each sample takes
analysis_date_time = [i + ' ' + j for i,j in zip(Date, Time)]
timediff = [(dt.datetime.strptime(i, '%m/%d/%y %H:%M:%S').timestamp() - dt.datetime.strptime(j, '%m/%d/%y %H:%M:%S').timestamp()) for i, j in zip(analysis_date_time[1:], analysis_date_time[:-1]) if len(i)==17 and len(j)==17]
timediff = [i for i in timediff if i>2000 and i<6000]
meantimediff = np.mean(timediff)

# Normalize peak area to reference peak fluctuations
AreaAll_wg_norm_ratio = AreaAll_wg / np.nanmean(AreaAll_wg)
AreaAll_sam_norm = AreaAll_sam / AreaAll_wg_norm_ratio

# blank correction
blank['mean_peak_area'] = np.mean(AreaAll_sam_norm[blank['index']])
blank['mean_d34S'] = np.mean(d34S32S_sam[blank['index']])
d34S_blank_corr = ((d34S32S_sam * AreaAll_sam_norm) - (blank['mean_d34S'] * blank['mean_peak_area'])) / AreaAll_sam_norm
AreaAll_sam_norm_blank_corr = AreaAll_sam_norm - blank['mean_peak_area']


#np.mean(Area66_sam[blank['index']]/Area64_sam[blank['index']])

R66_wg = Area66_wg / Area64_wg
R66_sam = Area66_sam / Area64_sam
d66 = (R66_sam / R66_wg - 1) * 1000
d66_blank_corrected = ((AreaAll_sam * d66) - (np.mean(AreaAll_sam[blank['index']]) * np.mean(d66[blank['index']]))) / (AreaAll_sam - np.mean(AreaAll_sam[blank['index']]))
R66_sam_blank_corrected = (d66_blank_corrected/1000+1)*np.mean(R66_wg)


# -------------------- Sqty vs Peak Area least squares fit --------------------
# Sfit = lab.linreg(Sqty[qtycal['index'] + blank['index']], AreaAll_sam_norm[qtycal['index'] + blank['index']])  # slope, intercept
# Sfit = lab.linreg(Sqty[qtycal['index']], AreaAll_sam_norm_blank_corr[qtycal['index']])  # slope, intercept
Sfit = lab.linreg(Sqty[qtycal['index']], AreaAll_sam_norm[qtycal['index']])  # slope, intercept
Sfit_eq_str = f"y = {round(Sfit[0], 1)} * x + {round(Sfit[1], 1)}"
Sqty_calc = (AreaAll_sam_norm - Sfit[1]) / Sfit[0]
Sqty_calc[blank['index']] = AreaAll_sam_norm[blank['index']] / Sfit[0]
qtycal['area_sort_i'] = np.argsort(AreaAll_sam_norm[qtycal['index']])
# Sfit_residual = AreaAll_sam_norm[all_stds] - (Sfit[0] * Sqty[all_stds] + Sfit[1])

# blank size analysis
blank['Sqty_calc'] = AreaAll_sam_norm[blank['index']] / Sfit[0]

non_blanks = np.setdiff1d(all_analyses, np.union1d(blank['index'],void['index']))

blank_void_index = np.sort(np.union1d(blank['index'],void['index']))
blank_void_diff = np.diff(blank_void_index)
clean_blanks = np.intersect1d(blank_void_index[np.where(blank_void_diff==1)[0]]+1, blank['index'])  # blanks that either had a blank or a void immediately prior
dirty_blanks = np.setdiff1d(blank['index'], clean_blanks)


# microbalance error estimation
microbalance_error = 0.002  # precision in mg
Na2SO4['Sqty_error'] = microbalance_error * 1000 * Na2SO4['fractionS']
Na2SO4['Sqty_lower'] = Sqty[Na2SO4['index']] - Na2SO4['Sqty_error']
Na2SO4['Sqty_upper'] = Sqty[Na2SO4['index']] + Na2SO4['Sqty_error']


# # -------------------- Correct to VCDT --------------------

d34S_vcdt_fit = lab.linreg([np.mean(d34S32S_sam[Ag2S['index']]), np.mean(d34S32S_sam[ZnS['index']])], [Ag2S['d34S'], ZnS['d34S']])

d34S_vcdt = d34S_vcdt_fit[0] * d34S32S_sam + d34S_vcdt_fit[1]


# d34S residual
d34S_residual = np.full(len(d34S32S_sam), np.nan)
d34S_residual[SodSul_1['index']] = d34S32S_sam[SodSul_1['index']] - np.mean(d34S32S_sam[SodSul_1['index']])
d34S_residual[SodSul_2['index']] = d34S32S_sam[SodSul_2['index']] - np.mean(d34S32S_sam[SodSul_2['index']])
d34S_residual[SodSul_3['index']] = d34S32S_sam[SodSul_3['index']] - np.mean(d34S32S_sam[SodSul_3['index']])
d34S_residual[SodSul_4['index']] = d34S32S_sam[SodSul_4['index']] - np.mean(d34S32S_sam[SodSul_4['index']])
d34S_residual[SodSul_5['index']] = d34S32S_sam[SodSul_5['index']] - np.mean(d34S32S_sam[SodSul_5['index']])
d34S_residual[SodSul_6['index']] = d34S32S_sam[SodSul_6['index']] - np.mean(d34S32S_sam[SodSul_6['index']])
d34S_residual[SodSul_7['index']] = d34S32S_sam[SodSul_7['index']] - np.mean(d34S32S_sam[SodSul_7['index']])
d34S_residual[SodSul_8['index']] = d34S32S_sam[SodSul_8['index']] - np.mean(d34S32S_sam[SodSul_8['index']])
d34S_residual[SodSul_1_RotoVap['index']] = d34S32S_sam[SodSul_1_RotoVap['index']] - np.mean(d34S32S_sam[SodSul_1_RotoVap['index']])
d34S_residual[SodSul_1_Oven60['index']] = d34S32S_sam[SodSul_1_Oven60['index']] - np.mean(d34S32S_sam[SodSul_1_Oven60['index']])
d34S_residual[BaSO4['index']] = d34S32S_sam[BaSO4['index']] - np.mean(d34S32S_sam[BaSO4['index']])
d34S_residual[ZnS['index']] = d34S32S_sam[ZnS['index']] - np.mean(d34S32S_sam[ZnS['index']])
d34S_residual[Ag2S['index']] = d34S32S_sam[Ag2S['index']] - np.mean(d34S32S_sam[Ag2S['index']])

# -------------------- figures --------------------
print('    Making figures...')

marker_sizes = [5, 10, 15, 20, 25]
blank['marker_size'] = marker_sizes[0]
blank['marker_color'] = 'black'
qtycal['marker_size'] = marker_sizes[2]
qtycal['marker_color'] = 'yellow'
BaSO4['marker_size'] = marker_sizes[1]
BaSO4['marker_color'] = 'white'
Ag2S['marker_size'] = marker_sizes[1]
Ag2S['marker_color'] = 'black'
ZnS['marker_size'] = marker_sizes[1]
ZnS['marker_color'] = 'tan'
SodSul_1['marker_size'] = marker_sizes[2]
SodSul_1['marker_color'] = 'hotpink'
SodSul_1_RotoVap['marker_size'] = marker_sizes[2]
SodSul_1_RotoVap['marker_color'] = 'pink'
SodSul_1_Oven60['marker_size'] = marker_sizes[2]
SodSul_1_Oven60['marker_color'] = 'deeppink'
SodSul_2['marker_size'] = marker_sizes[3]
SodSul_2['marker_color'] = 'lavender'
SodSul_3['marker_size'] = marker_sizes[3]
SodSul_3['marker_color'] = 'violet'
SodSul_4['marker_size'] = marker_sizes[3]
SodSul_4['marker_color'] = 'magenta'
SodSul_5['marker_size'] = marker_sizes[3]
SodSul_5['marker_color'] = 'blueviolet'
SodSul_6['marker_size'] = marker_sizes[3]
SodSul_6['marker_color'] = 'darkmagenta'
SodSul_7['marker_size'] = marker_sizes[3]
SodSul_7['marker_color'] = 'darkslateblue'
SodSul_8['marker_size'] = marker_sizes[3]
SodSul_8['marker_color'] = 'mediumslateblue'


figures = {}
fig_n = 1

figures[fig_n] = {}
figures[fig_n]['cap'] = f"""Figure {fig_n}. Normalized peak area versus analysis number. Peak Area is normalized by the working gas to account for fluctuations
                        in the mass spec sensitivity."""
figures[fig_n]['fig'] = figure(width=1800, height=700, x_axis_label="Analysis", y_axis_label="Normalized Peak Area (Vs)",
                           tools="pan, box_zoom, reset, save", active_drag="box_zoom")
figures[fig_n]['fig'].circle(Analysis[blank['index']], AreaAll_sam_norm[blank['index']],
                         legend_label="Blank", size=blank['marker_size'], line_color='black', fill_color=blank['marker_color'])
figures[fig_n]['fig'].square(Analysis[SodSul_1['index']], AreaAll_sam_norm[SodSul_1['index']],
                         legend_label="SodSul_1 Microbalance", size=SodSul_1['marker_size'], line_color='black', fill_color=SodSul_1['marker_color'])
figures[fig_n]['fig'].triangle(Analysis[BaSO4['index']], AreaAll_sam_norm[BaSO4['index']],
                           legend_label="BaSO4", size=marker_sizes[1], line_color='black', fill_color=BaSO4['marker_color'])
figures[fig_n]['fig'].triangle(Analysis[Ag2S['index']], AreaAll_sam_norm[Ag2S['index']],
                           legend_label="Ag2S", size=marker_sizes[1], line_color='black', fill_color=Ag2S['marker_color'])
figures[fig_n]['fig'].triangle(Analysis[ZnS['index']], AreaAll_sam_norm[ZnS['index']],
                           legend_label="ZnS", size=marker_sizes[1], line_color='black', fill_color=ZnS['marker_color'])
figures[fig_n]['fig'].diamond(Analysis[SodSul_1_RotoVap['index']], AreaAll_sam_norm[SodSul_1_RotoVap['index']],
                          legend_label="SodSul_1_RotoVap", size=SodSul_1_RotoVap['marker_size'], line_color='black', fill_color=SodSul_1_RotoVap['marker_color'])
figures[fig_n]['fig'].square_pin(Analysis[SodSul_1_Oven60['index']], AreaAll_sam_norm[SodSul_1_Oven60['index']],
                          legend_label="SodSul_1_Oven60", size=SodSul_1_Oven60['marker_size'], line_color='black', fill_color=SodSul_1_Oven60['marker_color'])
figures[fig_n]['fig'].diamond(Analysis[SodSul_2['index']], AreaAll_sam_norm[SodSul_2['index']],
                          legend_label="SodSul_2", size=SodSul_2['marker_size'], line_color='black', fill_color=SodSul_2['marker_color'])
figures[fig_n]['fig'].diamond(Analysis[SodSul_3['index']], AreaAll_sam_norm[SodSul_3['index']],
                          legend_label="SodSul_3", size=SodSul_3['marker_size'], line_color='black', fill_color=SodSul_3['marker_color'])
figures[fig_n]['fig'].diamond(Analysis[SodSul_4['index']], AreaAll_sam_norm[SodSul_4['index']],
                          legend_label="SodSul_4", size=SodSul_4['marker_size'], line_color='black', fill_color=SodSul_4['marker_color'])
figures[fig_n]['fig'].diamond(Analysis[SodSul_5['index']], AreaAll_sam_norm[SodSul_5['index']],
                          legend_label="SodSul_5", size=SodSul_5['marker_size'], line_color='black', fill_color=SodSul_5['marker_color'])
figures[fig_n]['fig'].diamond(Analysis[SodSul_6['index']], AreaAll_sam_norm[SodSul_6['index']],
                          legend_label="SodSul_6", size=SodSul_6['marker_size'], line_color='black', fill_color=SodSul_6['marker_color'])
figures[fig_n]['fig'].diamond(Analysis[SodSul_7['index']], AreaAll_sam_norm[SodSul_7['index']],
                          legend_label="SodSul_7", size=SodSul_7['marker_size'], line_color='black', fill_color=SodSul_7['marker_color'])
figures[fig_n]['fig'].diamond(Analysis[SodSul_8['index']], AreaAll_sam_norm[SodSul_8['index']],
                          legend_label="SodSul_8", size=SodSul_8['marker_size'], line_color='black', fill_color=SodSul_8['marker_color'])
figures[fig_n]['fig'].diamond(Analysis[samples['index']], AreaAll_sam_norm[samples['index']],
                          legend_label="samples", size=12, line_color='black', fill_color='red')
[figures[fig_n]['fig'].line([Analysis[i], Analysis[i]], [np.nanmin(AreaAll_sam_norm), np.nanmax(AreaAll_sam_norm)]) for i in run_index_first_row]


fig_n += 1


figures[fig_n] = {}
figures[fig_n]['cap'] = f"""Figure {fig_n}. Normalized peak area versus target sulfur quantity. Peak Area is normalized by the working gas to account for fluctuations
                        in the mass spec sensitivity. The line of best fit is {Sfit_eq_str}."""
figures[fig_n]['fig'] = figure(width=1800, height=700, x_axis_label="Target sulfur quantity (ug)", y_axis_label="Normalized Peak Area (Vs)", tools="pan, box_zoom, reset, save", active_drag="box_zoom")

figures[fig_n]['fig'].circle(Sqty[blank['index']], AreaAll_sam_norm[blank['index']],
                         legend_label="Blank", size=blank['marker_size'], line_color='black', fill_color=blank['marker_color'])
figures[fig_n]['fig'].square(Sqty[SodSul_1['index']], AreaAll_sam_norm[SodSul_1['index']],
                         legend_label="SodSul_1 Microbalance", size=SodSul_1['marker_size'], line_color='black', fill_color=SodSul_1['marker_color'])
figures[fig_n]['fig'].triangle(Sqty[BaSO4['index']], AreaAll_sam_norm[BaSO4['index']],
                           legend_label="BaSO4", size=marker_sizes[1], line_color='black', fill_color=BaSO4['marker_color'])
figures[fig_n]['fig'].triangle(Sqty[Ag2S['index']], AreaAll_sam_norm[Ag2S['index']],
                           legend_label="Ag2S", size=marker_sizes[1], line_color='black', fill_color=Ag2S['marker_color'])
figures[fig_n]['fig'].triangle(Sqty[ZnS['index']], AreaAll_sam_norm[ZnS['index']],
                           legend_label="ZnS", size=marker_sizes[1], line_color='black', fill_color=ZnS['marker_color'])
figures[fig_n]['fig'].diamond(Sqty[SodSul_1_RotoVap['index']], AreaAll_sam_norm[SodSul_1_RotoVap['index']],
                          legend_label="SodSul_1_RotoVap", size=SodSul_1_RotoVap['marker_size'], line_color='black', fill_color=SodSul_1_RotoVap['marker_color'])
figures[fig_n]['fig'].square_pin(Sqty[SodSul_1_Oven60['index']], AreaAll_sam_norm[SodSul_1_Oven60['index']],
                          legend_label="SodSul_1_Oven60", size=SodSul_1_Oven60['marker_size'], line_color='black', fill_color=SodSul_1_Oven60['marker_color'])
figures[fig_n]['fig'].diamond(Sqty[SodSul_2['index']], AreaAll_sam_norm[SodSul_2['index']],
                          legend_label="SodSul_2", size=SodSul_2['marker_size'], line_color='black', fill_color=SodSul_2['marker_color'])
figures[fig_n]['fig'].diamond(Sqty[SodSul_3['index']], AreaAll_sam_norm[SodSul_3['index']],
                          legend_label="SodSul_3", size=SodSul_3['marker_size'], line_color='black', fill_color=SodSul_3['marker_color'])
figures[fig_n]['fig'].diamond(Sqty[SodSul_4['index']], AreaAll_sam_norm[SodSul_4['index']],
                          legend_label="SodSul_4", size=SodSul_4['marker_size'], line_color='black', fill_color=SodSul_4['marker_color'])
figures[fig_n]['fig'].diamond(Sqty[SodSul_5['index']], AreaAll_sam_norm[SodSul_5['index']],
                          legend_label="SodSul_5", size=SodSul_5['marker_size'], line_color='black', fill_color=SodSul_5['marker_color'])
figures[fig_n]['fig'].diamond(Sqty[SodSul_6['index']], AreaAll_sam_norm[SodSul_6['index']],
                          legend_label="SodSul_6", size=SodSul_6['marker_size'], line_color='black', fill_color=SodSul_6['marker_color'])
figures[fig_n]['fig'].diamond(Sqty[SodSul_7['index']], AreaAll_sam_norm[SodSul_7['index']],
                          legend_label="SodSul_7", size=SodSul_7['marker_size'], line_color='black', fill_color=SodSul_7['marker_color'])
figures[fig_n]['fig'].diamond(Sqty[SodSul_8['index']], AreaAll_sam_norm[SodSul_8['index']],
                          legend_label="SodSul_8", size=SodSul_8['marker_size'], line_color='black', fill_color=SodSul_8['marker_color'])
figures[fig_n]['fig'].diamond(Sqty[samples['index']], AreaAll_sam_norm[samples['index']],
                          legend_label="samples", size=12, line_color='black', fill_color='red')
figures[fig_n]['fig'].line(Sqty_calc[qtycal['index']][qtycal['area_sort_i']], AreaAll_sam_norm[qtycal['index']][qtycal['area_sort_i']])


fig_n += 1


figures[fig_n] = {}
figures[fig_n]['cap'] = f"""Figure {fig_n}. Calculated sulfur quantity versus analysis number."""
figures[fig_n]['fig'] = figure(width=1800, height=700, x_axis_label="Analysis", y_axis_label="Sulfur quantity (ug)",
                           tools="pan, box_zoom, reset, save", active_drag="box_zoom")
figures[fig_n]['fig'].circle(Analysis[blank['index']], Sqty_calc[blank['index']],
                         legend_label="Blank", size=blank['marker_size'], line_color='black', fill_color=blank['marker_color'])
figures[fig_n]['fig'].square(Analysis[SodSul_1['index']], Sqty_calc[SodSul_1['index']],
                         legend_label="SodSul_1 Microbalance", size=SodSul_1['marker_size'], line_color='black', fill_color=SodSul_1['marker_color'])
figures[fig_n]['fig'].triangle(Analysis[BaSO4['index']], Sqty_calc[BaSO4['index']],
                           legend_label="BaSO4", size=marker_sizes[1], line_color='black', fill_color=BaSO4['marker_color'])
figures[fig_n]['fig'].triangle(Analysis[Ag2S['index']], Sqty_calc[Ag2S['index']],
                           legend_label="Ag2S", size=marker_sizes[1], line_color='black', fill_color=Ag2S['marker_color'])
figures[fig_n]['fig'].triangle(Analysis[ZnS['index']], Sqty_calc[ZnS['index']],
                           legend_label="ZnS", size=marker_sizes[1], line_color='black', fill_color=ZnS['marker_color'])
figures[fig_n]['fig'].diamond(Analysis[SodSul_1_RotoVap['index']], Sqty_calc[SodSul_1_RotoVap['index']],
                          legend_label="SodSul_1_RotoVap", size=SodSul_1_RotoVap['marker_size'], line_color='black', fill_color=SodSul_1_RotoVap['marker_color'])
figures[fig_n]['fig'].square_pin(Analysis[SodSul_1_Oven60['index']], Sqty_calc[SodSul_1_Oven60['index']],
                          legend_label="SodSul_1_Oven60", size=SodSul_1_Oven60['marker_size'], line_color='black', fill_color=SodSul_1_Oven60['marker_color'])
figures[fig_n]['fig'].diamond(Analysis[SodSul_2['index']], Sqty_calc[SodSul_2['index']],
                          legend_label="SodSul_2", size=SodSul_2['marker_size'], line_color='black', fill_color=SodSul_2['marker_color'])
figures[fig_n]['fig'].diamond(Analysis[SodSul_3['index']], Sqty_calc[SodSul_3['index']],
                          legend_label="SodSul_3", size=SodSul_3['marker_size'], line_color='black', fill_color=SodSul_3['marker_color'])
figures[fig_n]['fig'].diamond(Analysis[SodSul_4['index']], Sqty_calc[SodSul_4['index']],
                          legend_label="SodSul_4", size=SodSul_4['marker_size'], line_color='black', fill_color=SodSul_4['marker_color'])
figures[fig_n]['fig'].diamond(Analysis[SodSul_5['index']], Sqty_calc[SodSul_5['index']],
                          legend_label="SodSul_5", size=SodSul_5['marker_size'], line_color='black', fill_color=SodSul_5['marker_color'])
figures[fig_n]['fig'].diamond(Analysis[SodSul_6['index']], Sqty_calc[SodSul_6['index']],
                          legend_label="SodSul_6", size=SodSul_6['marker_size'], line_color='black', fill_color=SodSul_6['marker_color'])
figures[fig_n]['fig'].diamond(Analysis[SodSul_7['index']], Sqty_calc[SodSul_7['index']],
                          legend_label="SodSul_7", size=SodSul_7['marker_size'], line_color='black', fill_color=SodSul_7['marker_color'])
figures[fig_n]['fig'].diamond(Analysis[SodSul_8['index']], Sqty_calc[SodSul_8['index']],
                          legend_label="SodSul_8", size=SodSul_8['marker_size'], line_color='black', fill_color=SodSul_8['marker_color'])
figures[fig_n]['fig'].diamond(Analysis[samples['index']], Sqty_calc[samples['index']],
                          legend_label="samples", size=12, line_color='black', fill_color='red')
[figures[fig_n]['fig'].line([Analysis[i], Analysis[i]], [np.nanmin(Sqty_calc), np.nanmax(Sqty_calc)]) for i in run_index_first_row]


fig_n += 1


figures[fig_n] = {}
figures[fig_n]['cap'] = f"""Figure {fig_n}. Yield."""
figures[fig_n]['fig'] = figure(width=1800, height=700, x_axis_label="Analysis", y_axis_label="Yield (Vs / ug)",
                           tools="pan, box_zoom, reset, save", active_drag="box_zoom")
figures[fig_n]['fig'].square(Analysis[SodSul_1['index']], AreaAll_sam_norm[SodSul_1['index']]/Sqty[SodSul_1['index']],
                         legend_label="SodSul_1 Microbalance", size=SodSul_1['marker_size'], line_color='black', fill_color=SodSul_1['marker_color'])
figures[fig_n]['fig'].triangle(Analysis[BaSO4['index']], AreaAll_sam_norm[BaSO4['index']]/Sqty[BaSO4['index']],
                           legend_label="BaSO4", size=marker_sizes[1], line_color='black', fill_color=BaSO4['marker_color'])
figures[fig_n]['fig'].triangle(Analysis[Ag2S['index']], AreaAll_sam_norm[Ag2S['index']]/Sqty[Ag2S['index']],
                           legend_label="Ag2S", size=marker_sizes[1], line_color='black', fill_color=Ag2S['marker_color'])
figures[fig_n]['fig'].triangle(Analysis[ZnS['index']], AreaAll_sam_norm[ZnS['index']]/Sqty[ZnS['index']],
                           legend_label="ZnS", size=marker_sizes[1], line_color='black', fill_color=ZnS['marker_color'])
figures[fig_n]['fig'].diamond(Analysis[SodSul_1_RotoVap['index']], AreaAll_sam_norm[SodSul_1_RotoVap['index']]/Sqty[SodSul_1_RotoVap['index']],
                          legend_label="SodSul_1_RotoVap", size=SodSul_1_RotoVap['marker_size'], line_color='black', fill_color=SodSul_1_RotoVap['marker_color'])
figures[fig_n]['fig'].square_pin(Analysis[SodSul_1_Oven60['index']], AreaAll_sam_norm[SodSul_1_Oven60['index']]/Sqty[SodSul_1_Oven60['index']],
                          legend_label="SodSul_1_Oven60", size=SodSul_1_Oven60['marker_size'], line_color='black', fill_color=SodSul_1_Oven60['marker_color'])
figures[fig_n]['fig'].diamond(Analysis[SodSul_2['index']], AreaAll_sam_norm[SodSul_2['index']]/Sqty[SodSul_2['index']],
                          legend_label="SodSul_2", size=SodSul_2['marker_size'], line_color='black', fill_color=SodSul_2['marker_color'])
figures[fig_n]['fig'].diamond(Analysis[SodSul_3['index']], AreaAll_sam_norm[SodSul_3['index']]/Sqty[SodSul_3['index']],
                          legend_label="SodSul_3", size=SodSul_3['marker_size'], line_color='black', fill_color=SodSul_3['marker_color'])
figures[fig_n]['fig'].diamond(Analysis[SodSul_4['index']], AreaAll_sam_norm[SodSul_4['index']]/Sqty[SodSul_4['index']],
                          legend_label="SodSul_4", size=SodSul_4['marker_size'], line_color='black', fill_color=SodSul_4['marker_color'])
figures[fig_n]['fig'].diamond(Analysis[SodSul_5['index']], AreaAll_sam_norm[SodSul_5['index']]/Sqty[SodSul_5['index']],
                          legend_label="SodSul_5", size=SodSul_5['marker_size'], line_color='black', fill_color=SodSul_5['marker_color'])
figures[fig_n]['fig'].diamond(Analysis[SodSul_6['index']], AreaAll_sam_norm[SodSul_6['index']]/Sqty[SodSul_6['index']],
                          legend_label="SodSul_6", size=SodSul_6['marker_size'], line_color='black', fill_color=SodSul_6['marker_color'])
figures[fig_n]['fig'].diamond(Analysis[SodSul_7['index']], AreaAll_sam_norm[SodSul_7['index']]/Sqty[SodSul_7['index']],
                          legend_label="SodSul_7", size=SodSul_7['marker_size'], line_color='black', fill_color=SodSul_7['marker_color'])
figures[fig_n]['fig'].diamond(Analysis[SodSul_8['index']], AreaAll_sam_norm[SodSul_8['index']]/Sqty[SodSul_8['index']],
                          legend_label="SodSul_8", size=SodSul_8['marker_size'], line_color='black', fill_color=SodSul_8['marker_color'])
[figures[fig_n]['fig'].line([Analysis[i], Analysis[i]], [0, 200]) for i in run_index_first_row]


fig_n += 1


figures[fig_n] = {}
figures[fig_n]['cap'] = f"""Figure {fig_n}. Measured d34S vs working gas."""
figures[fig_n]['fig'] = figure(width=1800, height=700, x_axis_label="Analysis", y_axis_label="d34S measured (permil)",
                           tools="pan, box_zoom, reset, save", active_drag="box_zoom")
figures[fig_n]['fig'].circle(Analysis[blank['index']], d34S32S_sam[blank['index']],
                         legend_label="Blank", size=blank['marker_size'], line_color='black', fill_color=blank['marker_color'])
figures[fig_n]['fig'].square(Analysis[SodSul_1['index']], d34S32S_sam[SodSul_1['index']],
                         legend_label="SodSul_1 Microbalance", size=SodSul_1['marker_size'], line_color='black', fill_color=SodSul_1['marker_color'])
figures[fig_n]['fig'].triangle(Analysis[BaSO4['index']], d34S32S_sam[BaSO4['index']],
                           legend_label="BaSO4", size=marker_sizes[1], line_color='black', fill_color=BaSO4['marker_color'])
figures[fig_n]['fig'].triangle(Analysis[Ag2S['index']], d34S32S_sam[Ag2S['index']],
                           legend_label="Ag2S", size=marker_sizes[1], line_color='black', fill_color=Ag2S['marker_color'])
figures[fig_n]['fig'].triangle(Analysis[ZnS['index']], d34S32S_sam[ZnS['index']],
                           legend_label="ZnS", size=marker_sizes[1], line_color='black', fill_color=ZnS['marker_color'])
figures[fig_n]['fig'].diamond(Analysis[SodSul_1_RotoVap['index']], d34S32S_sam[SodSul_1_RotoVap['index']],
                          legend_label="SodSul_1_RotoVap", size=SodSul_1_RotoVap['marker_size'], line_color='black', fill_color=SodSul_1_RotoVap['marker_color'])
figures[fig_n]['fig'].square_pin(Analysis[SodSul_1_Oven60['index']], d34S32S_sam[SodSul_1_Oven60['index']],
                          legend_label="SodSul_1_Oven60", size=SodSul_1_Oven60['marker_size'], line_color='black', fill_color=SodSul_1_Oven60['marker_color'])
figures[fig_n]['fig'].diamond(Analysis[SodSul_2['index']], d34S32S_sam[SodSul_2['index']],
                          legend_label="SodSul_2", size=SodSul_2['marker_size'], line_color='black', fill_color=SodSul_2['marker_color'])
figures[fig_n]['fig'].diamond(Analysis[SodSul_3['index']], d34S32S_sam[SodSul_3['index']],
                          legend_label="SodSul_3", size=SodSul_3['marker_size'], line_color='black', fill_color=SodSul_3['marker_color'])
figures[fig_n]['fig'].diamond(Analysis[SodSul_4['index']], d34S32S_sam[SodSul_4['index']],
                          legend_label="SodSul_4", size=SodSul_4['marker_size'], line_color='black', fill_color=SodSul_4['marker_color'])
figures[fig_n]['fig'].diamond(Analysis[SodSul_5['index']], d34S32S_sam[SodSul_5['index']],
                          legend_label="SodSul_5", size=SodSul_5['marker_size'], line_color='black', fill_color=SodSul_5['marker_color'])
figures[fig_n]['fig'].diamond(Analysis[SodSul_6['index']], d34S32S_sam[SodSul_6['index']],
                          legend_label="SodSul_6", size=SodSul_6['marker_size'], line_color='black', fill_color=SodSul_6['marker_color'])
figures[fig_n]['fig'].diamond(Analysis[SodSul_7['index']], d34S32S_sam[SodSul_7['index']],
                          legend_label="SodSul_7", size=SodSul_7['marker_size'], line_color='black', fill_color=SodSul_7['marker_color'])
figures[fig_n]['fig'].diamond(Analysis[SodSul_8['index']], d34S32S_sam[SodSul_8['index']],
                          legend_label="SodSul_8", size=SodSul_8['marker_size'], line_color='black', fill_color=SodSul_8['marker_color'])
figures[fig_n]['fig'].diamond(Analysis[samples['index']], d34S32S_sam[samples['index']],
                          legend_label="samples", size=12, line_color='black', fill_color='red')
[figures[fig_n]['fig'].line([Analysis[i], Analysis[i]], [np.nanmin(d34S32S_sam), np.nanmax(d34S32S_sam)]) for i in run_index_first_row]



fig_n += 1



figures[fig_n] = {}
figures[fig_n]['cap'] = f"""Figure {fig_n}. d34S residual vs calculated Sulfur quantity."""
figures[fig_n]['fig'] = figure(width=1800, height=700, x_axis_label="Sqty_calc", y_axis_label="residual d34S (permil)",
                           tools="pan, box_zoom, reset, save", active_drag="box_zoom")
figures[fig_n]['fig'].square(Sqty_calc[SodSul_1['index']], d34S32S_sam[SodSul_1['index']] - np.nanmean(d34S32S_sam[SodSul_1['index']]),
                         legend_label="SodSul_1 Microbalance", size=SodSul_1['marker_size'], line_color='black', fill_color=SodSul_1['marker_color'])
figures[fig_n]['fig'].triangle(Sqty_calc[BaSO4['index']], d34S32S_sam[BaSO4['index']] - np.nanmean(d34S32S_sam[BaSO4['index']]),
                           legend_label="BaSO4", size=marker_sizes[1], line_color='black', fill_color=BaSO4['marker_color'])
figures[fig_n]['fig'].triangle(Sqty_calc[Ag2S['index']], d34S32S_sam[Ag2S['index']] - np.nanmean(d34S32S_sam[Ag2S['index']]),
                           legend_label="Ag2S", size=marker_sizes[1], line_color='black', fill_color=Ag2S['marker_color'])
figures[fig_n]['fig'].triangle(Sqty_calc[ZnS['index']], d34S32S_sam[ZnS['index']] - np.nanmean(d34S32S_sam[ZnS['index']]),
                           legend_label="ZnS", size=marker_sizes[1], line_color='black', fill_color=ZnS['marker_color'])
figures[fig_n]['fig'].diamond(Sqty_calc[SodSul_1_RotoVap['index']], d34S32S_sam[SodSul_1_RotoVap['index']] - np.nanmean(d34S32S_sam[SodSul_1_RotoVap['index']]),
                          legend_label="SodSul_1_RotoVap", size=SodSul_1_RotoVap['marker_size'], line_color='black', fill_color=SodSul_1_RotoVap['marker_color'])
figures[fig_n]['fig'].square_pin(Sqty_calc[SodSul_1_Oven60['index']], d34S32S_sam[SodSul_1_Oven60['index']] - np.nanmean(d34S32S_sam[SodSul_1_Oven60['index']]),
                          legend_label="SodSul_1_Oven60", size=SodSul_1_Oven60['marker_size'], line_color='black', fill_color=SodSul_1_Oven60['marker_color'])
figures[fig_n]['fig'].diamond(Sqty_calc[SodSul_2['index']], d34S32S_sam[SodSul_2['index']] - np.nanmean(d34S32S_sam[SodSul_2['index']]),
                          legend_label="SodSul_2", size=SodSul_2['marker_size'], line_color='black', fill_color=SodSul_2['marker_color'])
figures[fig_n]['fig'].diamond(Sqty_calc[SodSul_3['index']], d34S32S_sam[SodSul_3['index']] - np.nanmean(d34S32S_sam[SodSul_3['index']]),
                          legend_label="SodSul_3", size=SodSul_3['marker_size'], line_color='black', fill_color=SodSul_3['marker_color'])
figures[fig_n]['fig'].diamond(Sqty_calc[SodSul_4['index']], d34S32S_sam[SodSul_4['index']] - np.nanmean(d34S32S_sam[SodSul_4['index']]),
                          legend_label="SodSul_4", size=SodSul_4['marker_size'], line_color='black', fill_color=SodSul_4['marker_color'])
figures[fig_n]['fig'].diamond(Sqty_calc[SodSul_5['index']], d34S32S_sam[SodSul_5['index']] - np.nanmean(d34S32S_sam[SodSul_5['index']]),
                          legend_label="SodSul_5", size=SodSul_5['marker_size'], line_color='black', fill_color=SodSul_5['marker_color'])
figures[fig_n]['fig'].diamond(Sqty_calc[SodSul_6['index']], d34S32S_sam[SodSul_6['index']] - np.nanmean(d34S32S_sam[SodSul_6['index']]),
                          legend_label="SodSul_6", size=SodSul_6['marker_size'], line_color='black', fill_color=SodSul_6['marker_color'])
figures[fig_n]['fig'].diamond(Sqty_calc[SodSul_7['index']], d34S32S_sam[SodSul_7['index']] - np.nanmean(d34S32S_sam[SodSul_7['index']]),
                          legend_label="SodSul_7", size=SodSul_7['marker_size'], line_color='black', fill_color=SodSul_7['marker_color'])
figures[fig_n]['fig'].diamond(Sqty_calc[SodSul_8['index']], d34S32S_sam[SodSul_8['index']] - np.nanmean(d34S32S_sam[SodSul_8['index']]),
                          legend_label="SodSul_8", size=SodSul_8['marker_size'], line_color='black', fill_color=SodSul_8['marker_color'])


fig_n += 1


figures[fig_n] = {}
figures[fig_n]['cap'] = f"""Figure {fig_n}. d34S vcdt vs d34S measured."""
figures[fig_n]['fig'] = figure(width=1800, height=700, x_axis_label="d34S measured", y_axis_label="d34S VCDT (permil)",
                           tools="pan, box_zoom, reset, save", active_drag="box_zoom")
figures[fig_n]['fig'].square(d34S32S_sam[SodSul_1['index']], d34S_vcdt[SodSul_1['index']],
                         legend_label="SodSul_1 Microbalance", size=SodSul_1['marker_size'], line_color='black', fill_color=SodSul_1['marker_color'])
figures[fig_n]['fig'].triangle(d34S32S_sam[BaSO4['index']], d34S_vcdt[BaSO4['index']],
                           legend_label="BaSO4", size=marker_sizes[1], line_color='black', fill_color=BaSO4['marker_color'])
figures[fig_n]['fig'].triangle(d34S32S_sam[Ag2S['index']], d34S_vcdt[Ag2S['index']],
                           legend_label="Ag2S", size=marker_sizes[1], line_color='black', fill_color=Ag2S['marker_color'])
figures[fig_n]['fig'].triangle(d34S32S_sam[ZnS['index']], d34S_vcdt[ZnS['index']],
                           legend_label="ZnS", size=marker_sizes[1], line_color='black', fill_color=ZnS['marker_color'])
figures[fig_n]['fig'].diamond(d34S32S_sam[SodSul_1_RotoVap['index']], d34S_vcdt[SodSul_1_RotoVap['index']],
                          legend_label="SodSul_1_RotoVap", size=SodSul_1_RotoVap['marker_size'], line_color='black', fill_color=SodSul_1_RotoVap['marker_color'])
figures[fig_n]['fig'].square_pin(d34S32S_sam[SodSul_1_Oven60['index']], d34S_vcdt[SodSul_1_Oven60['index']],
                          legend_label="SodSul_1_Oven60", size=SodSul_1_Oven60['marker_size'], line_color='black', fill_color=SodSul_1_Oven60['marker_color'])
figures[fig_n]['fig'].diamond(d34S32S_sam[SodSul_2['index']], d34S_vcdt[SodSul_2['index']],
                          legend_label="SodSul_2", size=SodSul_2['marker_size'], line_color='black', fill_color=SodSul_2['marker_color'])
figures[fig_n]['fig'].diamond(d34S32S_sam[SodSul_3['index']], d34S_vcdt[SodSul_3['index']],
                          legend_label="SodSul_3", size=SodSul_3['marker_size'], line_color='black', fill_color=SodSul_3['marker_color'])
figures[fig_n]['fig'].diamond(d34S32S_sam[SodSul_4['index']], d34S_vcdt[SodSul_4['index']],
                          legend_label="SodSul_4", size=SodSul_4['marker_size'], line_color='black', fill_color=SodSul_4['marker_color'])
figures[fig_n]['fig'].diamond(d34S32S_sam[SodSul_5['index']], d34S_vcdt[SodSul_5['index']],
                          legend_label="SodSul_5", size=SodSul_5['marker_size'], line_color='black', fill_color=SodSul_5['marker_color'])
figures[fig_n]['fig'].diamond(d34S32S_sam[SodSul_6['index']], d34S_vcdt[SodSul_6['index']],
                          legend_label="SodSul_6", size=SodSul_6['marker_size'], line_color='black', fill_color=SodSul_6['marker_color'])
figures[fig_n]['fig'].diamond(d34S32S_sam[SodSul_7['index']], d34S_vcdt[SodSul_7['index']],
                          legend_label="SodSul_7", size=SodSul_7['marker_size'], line_color='black', fill_color=SodSul_7['marker_color'])
figures[fig_n]['fig'].diamond(d34S32S_sam[SodSul_8['index']], d34S_vcdt[SodSul_8['index']],
                          legend_label="SodSul_8", size=SodSul_8['marker_size'], line_color='black', fill_color=SodSul_8['marker_color'])
figures[fig_n]['fig'].diamond(d34S32S_sam[samples['index']], d34S_vcdt[samples['index']],
                          legend_label="samples", size=12, line_color='black', fill_color='red')
figures[fig_n]['fig'].legend.location = 'left'


fig_n += 1


figures[fig_n] = {}
figures[fig_n]['cap'] = f"""Figure {fig_n}. blank corrected d66."""
figures[fig_n]['fig'] = figure(width=1800, height=700, x_axis_label="Area All (Vs)", y_axis_label="d66 blank corrected (permil)",
                           tools="pan, box_zoom, reset, save", active_drag="box_zoom")
figures[fig_n]['fig'].circle(AreaAll_sam[blank['index']], d66[blank['index']],
                         legend_label="Blank", size=blank['marker_size'], line_color='black', fill_color=blank['marker_color'])
figures[fig_n]['fig'].square(AreaAll_sam[SodSul_1['index']], d66_blank_corrected[SodSul_1['index']],
                         legend_label="SodSul_1 Microbalance", size=SodSul_1['marker_size'], line_color='black', fill_color=SodSul_1['marker_color'])
figures[fig_n]['fig'].triangle(AreaAll_sam[BaSO4['index']], d66_blank_corrected[BaSO4['index']],
                           legend_label="BaSO4", size=marker_sizes[1], line_color='black', fill_color=BaSO4['marker_color'])
figures[fig_n]['fig'].triangle(AreaAll_sam[Ag2S['index']], d66_blank_corrected[Ag2S['index']],
                           legend_label="Ag2S", size=marker_sizes[1], line_color='black', fill_color=Ag2S['marker_color'])
figures[fig_n]['fig'].triangle(AreaAll_sam[ZnS['index']], d66_blank_corrected[ZnS['index']],
                           legend_label="ZnS", size=marker_sizes[1], line_color='black', fill_color=ZnS['marker_color'])
figures[fig_n]['fig'].legend.location = 'left'




# -------------------- export summary data file --------------------
print(f'\n    Creating summary data file.')
summary_data_filename = f'shrekS_summary_data_MiraRoth.csv'
summary_data_file = os.path.join(project_directory, report_dir, 'data/', summary_data_filename)
summary_file_headers = ['Sample ID', 'Date', 'Analysis Number', 'Mass (mg)', 'Peak Area (Vs)', 'Sulfur amount (ug)', 'd34S vs working gas (permil)']
data_to_write = '[Identifier1[ii], Date[ii], int(Analysis[ii]), Amount[ii], AreaAll_sam_norm[ii], Sqty[ii], round(d34S32S_sam[ii], 2)]'
data_to_write = str(data_to_write).replace("'", "")
with open(summary_data_file, 'w', newline='') as csvfile:
    datawriter = csv.writer(csvfile, quoting=csv.QUOTE_MINIMAL)
    datawriter.writerow(summary_file_headers)
    for ii in all_stds:
        datawriter.writerow(eval(data_to_write))
    datawriter.writerow('\n')
    for ii in samples['index']:
        datawriter.writerow(eval(data_to_write))


# -------------------- make html summary report --------------------

# copy report files
shutil.copy2(os.path.join(python_directory, 'py_report_style.css'), os.path.join(project_directory, report_dir))
[shutil.copy2(os.path.join(python_directory, script), os.path.join(project_directory, report_dir, f"python/{script}_REPORT_COPY")) for script in python_scripts]


header = f"""
    <!DOCTYPE html>
    <html lang="en">
    <head>
        <!-- py by Andrew Schauer -->
        <meta http-equiv="Content-Type" content="text/html; charset=UTF-8" />
        <meta name="viewport" content="width=device-width,initial-scale=1">
        <link rel="stylesheet" type="text/css" href="py_report_style.css">
        <title>Shrek S Report</title>
    </head>
    <body id="top">\n

    <div class="created-date">Created - {str(dt.datetime.now())}</div>

    <h2>Shrek Sulfur Report - {Date[run_index_first_row][0][0:8]} to {Date[run_index_first_row][-1][0:8]}</h2>

    <h2>Introduction</h2>
    <div class="text-indent">
    <p>This report is meant to be a stand-alone collection of methods,
    data, scripts, and notes related to sulfur isotopic (<sup>34</sup>S)
    analysis using Shrek (a MAT 253 coupled to an EA and Conflo III).
    Analytical details are described in Jongebloed et al. (2022) as well as the lab
    overview and method (see <a href="#refs">References</a> below).</p>

    <p>The data and python scripts used to generate this page are linked
    and described in the <a href="#refs">References</a> section below. If you wish
        to save this report, download the <a href="report.zip">zip file</a> and all html, images, data files, and
        python scripts will be saved. <strong><a href="report.zip">Save a copy if you are finished
        analyzing your samples</a></strong>.</p></div>

    <h2>My data</h2>
    <div class="text-indent">
        <p>All this technical stuff is fine and all but where are <strong><a href="data/{summary_data_filename}">my data</a></strong>?
        This summary file contains sample IDs, date, unique analysis numbers, mass of material, sulfur amount, and the d34S vs VCDT.
        Each section of data is separated by an empty row. The first section of data are the trusted reference materials; the second
        section of data are trusted samples. If you are done analyzing samples, please save a copy of
        the entire report directory elsewhere, not just a copy of your data file.</p>
    </div>

    <h2>Inventory</h2>
    <div class="text-indent">
        <table>
            <tr><th>Category</th><th>n</th></tr>
            <tr><td>Total number of runs</td><td>{len(run_index_first_row)}</td></tr>
            <tr><td>Total number of analyses</td><td>{len(Analysis)}</td></tr>
            <tr><td>Total number of standards analyzed <sup>*</sup></td><td>{len(all_stds)}</td></tr>
            <tr><td>Total number of samples analyzed</td><td>{len(samples['index'])}</td></tr>
            <tr><td><br></td></tr>
            <tr><td>Number of <a href="#excluded">excluded analyses</a></td><td>{len(flag0_indices)}</td></tr>
        </table>
        <sup>*</sup><small> - Reference material analyses are accumulated from many individual runs and projects into a single reference
                              frame until the instrumentation shifts and a new reference frame must be considered.</small>
    </div>

    <h2>Figures</h2>
        <div class="text-indent">"""


figure_block = [f"""<div class="clear-both">{file_html(figures[i]['fig'], CDN)}{figures[i]['cap']}<hr></div>""" for i in figures.keys()]


excluded_analyses_block = str([f"<tr><td>{original_data['Analysis'][i]}</td><td>{original_data['Identifier1'][i]}</td><td>{original_data['notes'][i]}</td></tr>" for i in flag0_indices]).replace("[","").replace("'","").replace("]","").replace(", ","")

excluded_analysis = f"""</div>
    \n<h2 id="excluded">Excluded analyses</h2>
    <div class="text-indent"><p>This is a list of analyses that have been excluded from further data processing showing
    the sample ID, the analysis number, and the reason for exclusion. If you want to included select analyses back into
    data processing, open the log file, find the sample of interest, and change flag to 1, then rerun the script.</p>

    <table>
        <tr><th>Sample ID</th><th>Analysis Number</th><th>Reason for excluding</th></tr>
        {excluded_analyses_block}
    </table></div>"""

python_scripts_string = [f'<a href="python/{script}_REPORT_COPY">{script}</a>, ' for script in python_scripts]
python_scripts_string = str(python_scripts_string)
python_scripts_string = python_scripts_string.replace("'","").replace("[","").replace("]","").replace(", , ", ", ")

footer = f"""
    \n<h2 id="refs">References</h2>
    <div class="references">
    <ul>
        <li>Jongebloed et al. 2022 submitted</li>
        <li>Python scripts - {python_scripts_string}</li>
        <li>Data files - <a href="data/{log_file_name}_REPORT_COPY">{log_file_name}</a>, <a href="data/{summary_data_filename}">{summary_data_filename}</a></li>
        <li><a href="https://isolab.ess.washington.edu/laboratory/solid-S.php">IsoLab's sulfur isotope analysis overiew.</a></li>
        <li><a href="https://isolab.ess.washington.edu/SOPs/shrek-s.php">IsoLab's sulfur isotope analysis method.</a></li>
        <li><a href="report.zip">Zip file of entire report directory.</a></strong>.</li>
    </ul>
    </div>
    </body></html>"""


log_summary_page = os.path.join(project_directory, report_dir, 'report.html')
with open(log_summary_page, 'w') as html_page:

    html_page.write(header)

    [html_page.write(i) for i in figure_block]

    html_page.write(excluded_analysis)

    html_page.write(footer)


webbrowser.open(log_summary_page)


# create zip file of entire report directory
shutil.make_archive('report', 'zip', os.path.join(project_directory, report_dir))
if os.path.exists(os.path.join(project_directory, report_dir, 'report.zip')):
    os.remove(os.path.join(project_directory, report_dir, 'report.zip'))
shutil.move(os.path.join('report.zip'), os.path.join(project_directory, report_dir))
