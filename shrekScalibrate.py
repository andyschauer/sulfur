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
"""

__authors__ = "Andy Schauer, Ursula Jongebloed"
__email__ = "aschauer@uw.edu"
__last_modified__ = "2023-12-04"
__version__ = "3.0"
__copyright__ = "Copyright 2023, Andy Schauer"
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
    src = os.path.join(home_dir, project_dir, log_file_name)
    dst = os.path.join(home_dir, project_dir, 'archive/', archive_log_file_name)
    shutil.copy2(src, dst)
    os.remove(os.path.join(home_dir, project_dir, log_file_name))

    # use the existing accepted log file headers (defined in shrekS_standards.py) to construct the export data
    data_to_write = [f"original_data['{i}'][ii]" for i in log_file_headers]
    data_to_write = str(data_to_write).replace('"', '')

    # write data to csv log file
    with open(os.path.join(home_dir, project_dir, log_file_name), 'w', newline='') as csvfile:
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
python_scripts = ['shrekS.py', 'shrekScalibrate.py', 'shrekS_standards.py']

if lab.where_am_i() == 'poppy':
    home_dir = '/home/aschauer/'
    python_dir = '/home/aschauer/python/pybob/'
else:
    home_dir = 'S:/Data/'
    python_dir = 'S:/Data/python/'

project_dir = 'projects/shrekS/'
report_dir = 'report/'
fig_dir = 'figures/'

log_file_name = 'shrekS_analysis_log.csv'
# log_file_name = 'shrekS_analysis_log_BlanksEmptyColumn.csv'

if os.path.isdir(os.path.join(home_dir, project_dir)) is False:
    print('\n    The project directory does not exist...exiting....')
    sys.exit()

# copy log file to report directory
shutil.copy2(os.path.join(home_dir, project_dir, log_file_name), os.path.join(home_dir, project_dir, report_dir, f"data/{log_file_name}_REPORT_COPY"))

# load reference material information
with open('/home/aschauer/ReferenceMaterials/reference_materials.json', 'r') as f:
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
headers, data = lab.read_file(os.path.join(home_dir, project_dir, log_file_name), ',')


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
qtycal['index'] = [i for i, e in enumerate(Identifier1) if str(e).lower() in (name.lower() for name in qtycal['names'])]
zero['index'] = [i for i, e in enumerate(Identifier1) if str(e).lower() in (name.lower() for name in zero['names'])]

Ag2S['index'] = [i for i, e in enumerate(Identifier1) if str(e).lower() in (name.lower() for name in Ag2S['names'])]
BaSO4['index'] = [i for i, e in enumerate(Identifier1) if str(e).lower() in (name.lower() for name in BaSO4['names'])]
Na2SO4['index'] = [i for i, e in enumerate(Identifier1) if str(e).lower() in (name.lower() for name in Na2SO4['names'])]
ZnS['index'] = [i for i, e in enumerate(Identifier1) if str(e).lower() in (name.lower() for name in ZnS['names'])]
all_stds = Ag2S['index'] + BaSO4['index'] + Na2SO4['index'] + ZnS['index'] + qtycal['index']

Ag2SO4['index'] = [i for i, e in enumerate(Identifier1) if str(e).lower() in (name.lower() for name in Ag2SO4['names'])]

samples = {'index': list(set(all_analyses) - set(Ag2S['index']) - set(BaSO4['index']) - set(Na2SO4['index']) - set(ZnS['index']) - set(void['index']) - set(blank['index']) - set(qtycal['index']) - set(zero['index']))}
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


S5p0ug = np.where(Sqty[qtycal['index']]==5.0)[0]
S1p0ug = np.where(Sqty[qtycal['index']]==1.0)[0]
S0p5ug = np.where(Sqty[qtycal['index']]==0.5)[0]
S0p1ug = np.where(Sqty[qtycal['index']]==0.1)[0]


# # -------------------- Calculations --------------------
print('    Maths...')
# calculate how long each sample takes
analysis_date_time = [i + ' ' + j for i,j in zip(Date, Time)]
timediff = [(dt.datetime.strptime(i, '%m/%d/%y %H:%M:%S').timestamp() - dt.datetime.strptime(j, '%m/%d/%y %H:%M:%S').timestamp()) for i, j in zip(analysis_date_time[1:], analysis_date_time[:-1]) if len(i)==17 and len(j)==17]
timediff = [i for i in timediff if i>2000 and i<4000]
meantimediff = np.mean(timediff)

# # Normalize peak area to reference peak fluctuations
AreaAll_wg_norm_ratio = AreaAll_wg / np.nanmean(AreaAll_wg)
AreaAll_sam_norm = AreaAll_sam * AreaAll_wg_norm_ratio
AreaAll_sam_norm = AreaAll_sam

# blank correction
blank['mean_peak_area'] = np.mean(AreaAll_sam_norm[blank['index']])
blank['mean_d34S'] = np.mean(d34S32S_sam[blank['index']])

d34S_blank_corr = ((d34S32S_sam * AreaAll_sam_norm) - (blank['mean_d34S'] * blank['mean_peak_area'])) / AreaAll_sam_norm

AreaAll_sam_norm_blank_corr = AreaAll_sam_norm - blank['mean_peak_area']


# -------------------- Sqty vs Peak Area least squares fit --------------------
Sfit = lab.linreg(Sqty[qtycal['index']], AreaAll_sam_norm[qtycal['index']])  # slope, intercept
Sfit_eq_str = f"y = {round(Sfit[0], 1)} * x + {round(Sfit[1], 1)}"
Sqty_calc = (AreaAll_sam_norm - Sfit[1]) / Sfit[0]
qtycal['area_sort_i'] = np.argsort(AreaAll_sam_norm[qtycal['index']])
Sfit_residual = AreaAll_sam_norm[all_stds] - (Sfit[0] * Sqty[all_stds] + Sfit[1])



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
# #    Currently correcting to VCDT in two distinct size bins
# ag2s['ilt20'] = [i for i in ag2s['index'] if Sqty_calc[i]<20]  # index less than 20 ug of S
# ag2s['igt20'] = [i for i in ag2s['index'] if Sqty_calc[i]>=20]  # index greater than or equal to 20 ug of S
# ZnS['ilt20'] = [i for i in ZnS['index'] if Sqty_calc[i]<20]  # index less than 20 ug of S
# ZnS['igt20'] = [i for i in ZnS['index'] if Sqty_calc[i]>=20]  # index greater than or equal to 20 ug of S
# ilt20 = [i for i in all_analyses if Sqty_calc[i]<20]
# igt20 = [i for i in all_analyses if Sqty_calc[i]>=20]

# d34S_vcdt_fits = {'lt20': lab.linreg([np.mean(d34S32S_sam[ag2s['ilt20']]), np.mean(d34S32S_sam[ZnS['ilt20']])], [ag2s['d34Sacc'], ZnS['d34Sacc']]),
                  # 'gt20': lab.linreg([np.mean(d34S32S_sam[ag2s['igt20']]), np.mean(d34S32S_sam[ZnS['igt20']])], [ag2s['d34Sacc'], ZnS['d34Sacc']])}

# d34S_vcdt = np.zeros(len(all_analyses))-999
# d34S_vcdt[ilt20] = d34S_vcdt_fits['lt20'][0] * d34S32S_sam[ilt20] + d34S_vcdt_fits['lt20'][1]
# d34S_vcdt[igt20] = d34S_vcdt_fits['gt20'][0] * d34S32S_sam[igt20] + d34S_vcdt_fits['gt20'][1]

# d34S32S_sam = S_sam_d34S32S
# d34S_vcdt = S_sam_d34S32S

# d34S residual
d34S_residual = [d34S32S_sam[Na2SO4['index']] - np.mean(d34S32S_sam[Na2SO4['index']])]
d34S_residual = np.append(d34S_residual, [d34S32S_sam[qtycal['index']] - np.mean(d34S32S_sam[qtycal['index']])])


# -------------------- figures --------------------
print('    Making figures...')

qtycal['S_quantities'] = np.unique(Sqty[qtycal['index']])
qtycal['palette'] = palettes.Reds[len(qtycal['S_quantities'])]
qtycal['marker'] = [10, 14, 18, 22]
qtycal['marker_color'] = [None for i in range(1,len(qtycal['index'])+1)]
qtycal['marker_size'] = [None for i in range(1,len(qtycal['index'])+1)]

for i, S in enumerate(Sqty[qtycal['index']]):
    for j, e in enumerate(qtycal['S_quantities']):
        if S==e:
            qtycal['marker_color'][i] = qtycal['palette'][j]
            qtycal['marker_size'][i] = qtycal['marker'][j]

figures = {}

figures[1] = {}
figures[1]['fig'] = figure(width=1800, height=700, x_axis_label="Analysis", y_axis_label="Normalized Peak Area (Vs)", tools="pan, box_zoom, reset, save", active_drag="box_zoom")
figures[1]['fig'].circle(Analysis[void['index']], AreaAll_sam_norm[void['index']], legend_label="Void", size=8, line_color='black', fill_color='cyan')
figures[1]['fig'].circle(Analysis[blank['index']], AreaAll_sam_norm[blank['index']], legend_label="Blank", size=12, line_color='black', fill_color='black')
figures[1]['fig'].triangle(Analysis[qtycal['index']], AreaAll_sam_norm[qtycal['index']], legend_label="Na2SO4", size=qtycal['marker_size'], line_color='black', fill_color=qtycal['marker_color'])
[figures[1]['fig'].line([Analysis[i], Analysis[i]], [np.nanmin(AreaAll_sam_norm), np.nanmax(AreaAll_sam_norm)]) for i in run_index_first_row]
figures[1]['cap'] = f"""Figure 1. Normalized peak area versus analysis number. Peak Area is normalized by the working gas to account for fluctuations
                        in the mass spec sensitivity."""

figures[2] = {}
figures[2]['fig'] = figure(width=1800, height=700, x_axis_label="Target sulfur quantity (ug)", y_axis_label="Normalized Peak Area (Vs)", tools="pan, box_zoom, reset, save", active_drag="box_zoom")
figures[2]['fig'].triangle(Sqty[qtycal['index']], AreaAll_sam_norm[qtycal['index']], legend_label="Na2SO4", size=qtycal['marker_size'], line_color='black', fill_color=qtycal['marker_color'])
figures[2]['fig'].line(Sqty_calc[qtycal['index']][qtycal['area_sort_i']], AreaAll_sam_norm[qtycal['index']][qtycal['area_sort_i']])
figures[2]['cap'] = f"""Figure 2. Normalized peak area versus target sulfur quantity. Peak Area is normalized by the working gas to account for fluctuations
                        in the mass spec sensitivity. The line of best fit is {Sfit_eq_str}."""

figures[3] = {}
figures[3]['fig'] = figure(width=1800, height=700, x_axis_label="Calculated sulfur quantity (ug)", y_axis_label="d34S vs working gas (permil)", tools="pan, box_zoom, reset, save", active_drag="box_zoom")
figures[3]['fig'].triangle(Sqty_calc[qtycal['index']], d34S32S_sam[qtycal['index']], legend_label="Na2SO4", size=qtycal['marker_size'], line_color='black', fill_color=qtycal['marker_color'])
figures[3]['cap'] = f"""Figure 3. d34S versus calculated sulfur quantity."""

figures[4] = {}
figures[4]['fig'] = figure(width=1800, height=700, x_axis_label="Analysis", y_axis_label="d34S vs working gas (permil)", tools="pan, box_zoom, reset, save", active_drag="box_zoom")
figures[4]['fig'].triangle(Analysis[qtycal['index']], d34S32S_sam[qtycal['index']], legend_label="Na2SO4", size=qtycal['marker_size'], line_color='black', fill_color=qtycal['marker_color'])
figures[4]['cap'] = f"""Figure 4. d34S versus analysis number. """


# -------------------- export summary data file --------------------
print(f'\n    Creating summary data file.')
summary_data_filename = f'shrekS_summary_data_MiraRoth.csv'
summary_data_file = os.path.join(home_dir, project_dir, report_dir, 'data/', summary_data_filename)
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
shutil.copy2(os.path.join(python_dir, 'py_report_style.css'), os.path.join(home_dir, project_dir, report_dir))
[shutil.copy2(os.path.join(python_dir, script), os.path.join(home_dir, project_dir, report_dir, f"python/{script}_REPORT_COPY")) for script in python_scripts]


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


log_summary_page = os.path.join(home_dir, project_dir, report_dir, 'report.html')
with open(log_summary_page, 'w') as html_page:

    html_page.write(header)

    [html_page.write(i) for i in figure_block]

    html_page.write(excluded_analysis)

    html_page.write(footer)


webbrowser.open(log_summary_page)


# create zip file of entire report directory
shutil.make_archive('report', 'zip', os.path.join(home_dir, project_dir, report_dir))
if os.path.exists(os.path.join(home_dir, project_dir, report_dir, 'report.zip')):
    os.remove(os.path.join(home_dir, project_dir, report_dir, 'report.zip'))
shutil.move(os.path.join('report.zip'), os.path.join(home_dir, project_dir, report_dir))
