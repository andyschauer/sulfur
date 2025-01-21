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
    version 4.1 - 2024.05.11 - changed lab import to isolab_lib, replaced old DIY least squares regression with numpy polyfit
    version 4.2 - 2024.05.12 - now using dateutil.parser because open .csv file in excel and saving changes the date and time format
    version 4.3 - 2024.06.09 - changed flag to trust, tried to incorporate CN code updates, 
    version 4.4 - 2024.07.03 - added Amplitude figure
"""

__authors__ = "Andy Schauer, Ursula Jongebloed"
__email__ = "aschauer@uw.edu"
__last_modified__ = "2024.06.09"
__version__ = "4.3"
__copyright__ = "Copyright 2025, Andy Schauer"
__license__ = "Apache 2.0"
__acknowledgements__ = "Alli Moon, Drew Pronovost"



# -------------------- imports --------------------
import argparse
from bokeh.embed import file_html
from bokeh.models import NumeralTickFormatter
from bokeh.plotting import figure, show
from bokeh.resources import CDN, INLINE
import csv
import datetime as dt
import dateutil.parser
import isolab_lib
import json
import matplotlib.pyplot as pplt
import numpy as np
import os
from shrekS_lib import *
import shutil
import sys
import time
import webbrowser



# ---------- FUNCTIONS ---------- 
def add_calculation_note(note):
    calculation_notes.append(note)



# -------------------- setup --------------------
print('Setup...')

# verbose mode
parser = argparse.ArgumentParser()
parser.add_argument("--verbose", help="Include exhaustive diagnostic information and figures in report.", action="store_true")
args = parser.parse_args()
if args.verbose:
    verbose = True
    print('\n    ****    Verbose mode on    ****    ')
else:
    verbose = False

version = os.path.basename(__file__) + ' - ' + time.ctime(os.path.getctime(__file__))

python_directory = isolab_lib.get_path("shrekS", "python")
method_directory = isolab_lib.get_path("shrekS", "project")
new_data_directory = 'rawdata_new'
archive_data_directory = 'rawdata_archive'
junk_data_directory = 'rawdata_junk'

python_scripts = {'shrekS_lib.py': '', 'shrekS.py': '', 'shrekScalibrate.py': ''}
python_scripts = {key: (time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(os.path.getmtime(f'{python_directory}{key}')))) for key, value in python_scripts.items()}

S_log_file_list = isolab_lib.make_file_list(method_directory, '_analysis_log.csv')

print('\nWhat analysis log file to you wish to process?\n')

[print(f'    {i}') for i in S_log_file_list]
identified_file = 0
while identified_file == 0:
    S_log_file_search = input('\nEnter a project analysis log file from above: ')
    isfile = [S_log_file_search[0: len(S_log_file_search)] in x for x in S_log_file_list]
    if len(np.where(isfile)[0]) == 1:
        identified_file = 1
        log_file_name = S_log_file_list[np.where(isfile)[0][0]]
        print(f'    Processing S log file {log_file_name}...')
    else:
        print('\n** More than one file found. **\n')

if os.path.isdir(method_directory) is False:
    print('Method directory does not exist...exiting....')
    sys.exit()

project_root_name = log_file_name[:-17]
project_directory = os.path.join(method_directory, f'{project_root_name}/')
if os.path.exists(project_directory)==False:
    os.mkdir(project_directory)

archive_directory = os.path.join(project_directory, 'archive/')
if os.path.exists(archive_directory)==False:
    os.mkdir(archive_directory)

report_directory = os.path.join(project_directory, "report/")
if os.path.exists(report_directory):
    shutil.move(report_directory, os.path.join(archive_directory, f"report_archive_{int(dt.datetime.utcnow().timestamp())}"))
os.mkdir(report_directory)
os.mkdir(os.path.join(report_directory, "data/"))
os.mkdir(os.path.join(report_directory, "python/"))
shutil.copy2(os.path.join(python_directory, 'py_report_style.css'), report_directory)
[shutil.copy2(os.path.join(python_directory, script), os.path.join(report_directory, f"python/{script}_REPORT_COPY")) for script in python_scripts]
shutil.copy2(os.path.join(method_directory, log_file_name), os.path.join(report_directory, 'data/'))
report_page = os.path.join(report_directory, f'{project_root_name}_calibration_summary.html')

# copy log file to report directory
shutil.copy2(os.path.join(method_directory, log_file_name), os.path.join(method_directory, report_directory, f"data/{log_file_name}_REPORT_COPY"))


# --------------- load reference material information ------------------------
with open(isolab_lib.get_path("shrekS", "standards"), 'r') as f:
    refmat = json.load(f)

refmat_keys = refmat['sulfates'].keys()
for i in refmat_keys:
    globals()[i] = refmat['sulfates'][i]
    globals()[i]['index'] = np.empty(0, dtype="int16")

refmat_keys = refmat['sulfides'].keys()
for i in refmat_keys:
    globals()[i] = refmat['sulfides'][i]
    globals()[i]['index'] = np.empty(0, dtype="int16")

refmat_keys = refmat['organics'].keys()
for i in refmat_keys:
    globals()[i] = refmat['organics'][i]
    globals()[i]['index'] = np.empty(0, dtype="int16")


# -------------------- get data --------------------
print('    Importing data...')
headers, data = isolab_lib.read_file(os.path.join(method_directory, log_file_name), ',')


# -------------------- organize data --------------------
print('    Organizing data...')

data['empty'] = ['' for i in data['Analysis']]  # currently a carriage return is probably written into the "empty" column, which I don't want

# remove trust 0 analyses from data
original_data = data.copy()
trust0_indices = [i for i, e in enumerate(data['trust']) if int(e) == 0]
trust1_indices = [i for i, e in enumerate(data['trust']) if int(e) == 1]
for header in headers[:-1]:
    data[header] = [data[header][index] for index in trust1_indices]

strlist = set(headers) - set(numlist)

for i in numlist:
    globals()[i] = np.asarray(data[i], dtype=float)

for i in strlist:
    globals()[i] = np.asarray(data[i])

run_index_first_row = [np.where(file == j)[0][0] for j in sorted(set(file))]


# -------------------- categorize all Identifier 1 (Identifier1) analyses into appropriate material type --------------------
calculation_notes = []

knowns_indices = []
for i in knowns_list:
    temp_indices = [j for j, e in enumerate(Identifier1) if str(e).lower() in (name.lower() for name in eval(i)['names'])]
    eval(i)['index'] = temp_indices
    knowns_indices.extend(eval(i)['index'])

included_isotope_standards = list(set([i for i in Identifier1 if i in refmat_list]))

sample_indices = list(set(range(0,len(Analysis))) - set(knowns_indices))

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
timediff = [(dateutil.parser.parse(i).timestamp() - dateutil.parser.parse(j).timestamp()) for i, j in zip(analysis_date_time[1:], analysis_date_time[:-1])]
timediff = [i for i in timediff if i>2000 and i<6000]
meantimediff = np.mean(timediff)

# Normalize peak area to reference peak fluctuations
AreaAll_wg_norm_ratio = AreaAll_wg / np.nanmean(AreaAll_wg)
AreaAll_sam_norm = AreaAll_sam / AreaAll_wg_norm_ratio
add_calculation_note("sample peak area corrected for working gas fluctuations")


# blank correction
blank = void  # in case no blanks are analyzed, treat voids as blanks
add_calculation_note("no blanks (empty tins) so using voids instead")
blank['mean_peak_area'] = np.mean(AreaAll_sam_norm[blank['index']])
blank['mean_d34S'] = np.mean(d34S32S_sam[blank['index']])

# (meas_d34S * meas_size) = (sample_d34S * sample_size) + (blank_d34S * blank_size)
# (meas_d34S * meas_size) - (blank_d34S * blank_size) = (sample_d34S * sample_size) 
# ((meas_d34S * meas_size) - (blank_d34S * blank_size)) / sample_size = sample_d34S 
# sample_size = meas_size - blank_size

d34S_blank_corr = ((d34S32S_sam * AreaAll_sam_norm) - (blank['mean_d34S'] * blank['mean_peak_area'])) / (AreaAll_sam_norm - blank['mean_peak_area'])
AreaAll_sam_norm_blank_corr = AreaAll_sam_norm - blank['mean_peak_area']
add_calculation_note("peak area and d34S were blank corrected")


# -------------------- Sqty vs Peak Area least squares fit --------------------
if qtycal_BaSO4['index']:
    qtycal_BaSO4['Sqty'] = Sqty[qtycal_BaSO4['index']]
    qtycal_BaSO4['Sfit'] = np.polyfit(qtycal_BaSO4['Sqty'], AreaAll_sam_norm[qtycal_BaSO4['index']], 1)

    # trying to find a better way to fit these data but still end up with something reasonable. 
    #    Here, a slope 
    A = qtycal_BaSO4['Sqty'][:,np.newaxis]
    slope, _, _, _ = np.linalg.lstsq(A, AreaAll_sam_norm_blank_corr[qtycal_BaSO4['index']], rcond=None)
    Sqty2_pred = AreaAll_sam_norm / slope

else:
    qtycal_BaSO4['Sqty'] = Sqty[qtycal_BaSO4['index']]
    qtycal_BaSO4['Sfit'] = [1,1]
    A = qtycal_BaSO4['Sqty'][:,np.newaxis]
    slope = 1
    Sqty2_pred = AreaAll_sam_norm / slope




qtycal_BaSO4['Sfit_eq_str'] = f"y = {round(qtycal_BaSO4['Sfit'][0], 1)} * x + {round(qtycal_BaSO4['Sfit'][1], 1)}"
Sqty_pred = (AreaAll_sam_norm - qtycal_BaSO4['Sfit'][1]) / qtycal_BaSO4['Sfit'][0]
Sqty_pred[blank['index']] = AreaAll_sam_norm[blank['index']] / qtycal_BaSO4['Sfit'][0]
qtycal_BaSO4['area_sort_i'] = np.argsort(AreaAll_sam_norm[qtycal_BaSO4['index']])

# blank size analysis
blank['Sqty_pred'] = AreaAll_sam_norm[blank['index']] / qtycal_BaSO4['Sfit'][0]


# ---------- Isotope Calibration Setup ----------
#     pick extreme isotope standards for maximum calibration range
d34S_std_1_index = np.argmin([eval(i)['d34S'] for i in included_isotope_standards])
d34S_std_1 = included_isotope_standards[d34S_std_1_index]
d34S_std_2_index = np.argmax([eval(i)['d34S'] for i in included_isotope_standards])
d34S_std_2 = included_isotope_standards[d34S_std_2_index]
d34S_std_3 = list(set(included_isotope_standards) - set([d34S_std_1, d34S_std_2]))[0]

for i in included_isotope_standards:
    eval(i)['purpose'] = ''

eval(d34S_std_1)['purpose'] += 'd34S calibration; '
eval(d34S_std_2)['purpose'] += 'd34S calibration; '
eval(d34S_std_3)['purpose'] += 'd34S quality control; '


# ---------- Residuals ----------
qtycal_BaSO4['Sresidual'] = AreaAll_sam_norm[qtycal_BaSO4['index']] - (qtycal_BaSO4['Sfit'][0] * qtycal_BaSO4['Sqty'] + qtycal_BaSO4['Sfit'][1])

Nqty_residual = np.concatenate(((Sqty[eval(d34S_std_1)['index']] - Amount[eval(d34S_std_1)['index']] * eval(d34S_std_1)['fractionS']) * 1000,
                                (Sqty[eval(d34S_std_2)['index']] - Amount[eval(d34S_std_2)['index']] * eval(d34S_std_2)['fractionS']) * 1000,
                                (Sqty[eval(d34S_std_3)['index']] - Amount[eval(d34S_std_3)['index']] * eval(d34S_std_3)['fractionS']) * 1000))

eval(d34S_std_1)['d34S_residual'] = d34S_blank_corr[eval(d34S_std_1)['index']] - np.nanmean(d34S_blank_corr[eval(d34S_std_1)['index']])
eval(d34S_std_2)['d34S_residual'] = d34S_blank_corr[eval(d34S_std_2)['index']] - np.nanmean(d34S_blank_corr[eval(d34S_std_2)['index']])
eval(d34S_std_3)['d34S_residual'] = d34S_blank_corr[eval(d34S_std_3)['index']] - np.nanmean(d34S_blank_corr[eval(d34S_std_3)['index']])
d34S_residual_std = np.std(np.concatenate([eval(d34S_std_1)['d34S_residual'], eval(d34S_std_2)['d34S_residual'], eval(d34S_std_3)['d34S_residual']]))


# ---------- Isotope Drift Calculation ----------
# Sdrift_fit = np.polyfit(np.concatenate([Analysis[eval(d34S_std_1)['index']], Analysis[eval(d34S_std_2)['index']]]),
                        # np.concatenate([eval(d34S_std_1)['d34S_residual'], eval(d34S_std_2)['d34S_residual']]),
                        # 1)
# Sdrift_corrfac = Sdrift_fit[0] * Analysis + Sdrift_fit[1]
# d34S_drift_corr = d34S_blank_corr - Sdrift_corrfac

d34S_drift_corr = d34S32S_sam
add_calculation_note("d34S was NOT corrected for drift")



# ---------- Isotope Calibration ----------
d34S_VCDT_fit = np.polyfit([np.nanmean(d34S_drift_corr[eval(d34S_std_1)['index']]), np.nanmean(d34S_drift_corr[eval(d34S_std_2)['index']])], [eval(d34S_std_1)['d34S'], eval(d34S_std_2)['d34S']], 1)
d34S_VCDT = d34S_VCDT_fit[0] * d34S_drift_corr + d34S_VCDT_fit[1]

add_calculation_note("d34S was normalized to VCDT, respectively")



# ---------- Post-normalization residual calculation ----------
eval(d34S_std_1)['d34S_VCDT_residual'] = d34S_VCDT[eval(d34S_std_1)['index']] - np.nanmean(d34S_VCDT[eval(d34S_std_1)['index']])
eval(d34S_std_2)['d34S_VCDT_residual'] = d34S_VCDT[eval(d34S_std_2)['index']] - np.nanmean(d34S_VCDT[eval(d34S_std_2)['index']])
eval(d34S_std_3)['d34S_VCDT_residual'] = d34S_VCDT[eval(d34S_std_3)['index']] - np.nanmean(d34S_VCDT[eval(d34S_std_3)['index']])
d34S_VCDT_residual_std = np.std(np.concatenate([eval(d34S_std_1)['d34S_VCDT_residual'], eval(d34S_std_2)['d34S_VCDT_residual'], eval(d34S_std_3)['d34S_VCDT_residual']]))


# -------------------- figures --------------------
print('    Making figures...')

marker_sizes = [10, 15, 20, 25, 30]
blank['marker_size'] = marker_sizes[0]
blank['marker_color'] = 'black'
qtycal_BaSO4['marker_size'] = marker_sizes[3]
qtycal_BaSO4['marker_color'] = 'yellow'
BaSO4['marker_size'] = marker_sizes[1]
BaSO4['marker_color'] = 'white'
Ag2S['marker_size'] = marker_sizes[1]
Ag2S['marker_color'] = 'black'
ZnS['marker_size'] = marker_sizes[1]
ZnS['marker_color'] = 'tan'
SodSul_1['marker_size'] = marker_sizes[2]
SodSul_1['marker_color'] = 'hotpink'
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
font_size = "24pt"
fig_n = 1

if verbose:

    figures[fig_n] = {}
    figures[fig_n]['cap'] = f"""Figure {fig_n}. Working gas peak amplitude versus analysis number. """
    figures[fig_n]['fig'] = figure(width=1800, height=700, x_axis_label="Analysis", y_axis_label="Working Gas Peak Area (Vs)", tools="pan, box_zoom, reset, save", active_drag="box_zoom")
    figures[fig_n]['fig'].square(Analysis, Ampl64_wg, size=marker_sizes[0], line_color='black', fill_color='black')
    figures[fig_n]['fig'].xaxis.formatter = NumeralTickFormatter(format="0")
    figures[fig_n]['fig'].xaxis.axis_label_text_font_size = font_size
    figures[fig_n]['fig'].yaxis.axis_label_text_font_size = font_size
    figures[fig_n]['fig'].xaxis.major_label_text_font_size = font_size
    figures[fig_n]['fig'].yaxis.major_label_text_font_size = font_size
    figures[fig_n]['fig'].title.text_font_size = font_size
    
    fig_n += 1
    
    figures[fig_n] = {}
    figures[fig_n]['cap'] = f"""Figure {fig_n}. Working gas peak area versus analysis number. If you see fluctuations, this is the reason to normalize sample peak area. Sample peak area is normalized by this working gas peak area to account for fluctuations in the mass spec sensitivity."""
    figures[fig_n]['fig'] = figure(width=1800, height=700, x_axis_label="Analysis", y_axis_label="Working Gas Peak Area (Vs)", tools="pan, box_zoom, reset, save", active_drag="box_zoom")
    figures[fig_n]['fig'].square(Analysis, AreaAll_wg, size=marker_sizes[0], line_color='black', fill_color='black')
    figures[fig_n]['fig'].xaxis.formatter = NumeralTickFormatter(format="0")
    figures[fig_n]['fig'].xaxis.axis_label_text_font_size = font_size
    figures[fig_n]['fig'].yaxis.axis_label_text_font_size = font_size
    figures[fig_n]['fig'].xaxis.major_label_text_font_size = font_size
    figures[fig_n]['fig'].yaxis.major_label_text_font_size = font_size
    figures[fig_n]['fig'].title.text_font_size = font_size

    fig_n += 1

    figures[fig_n] = {}
    figures[fig_n]['cap'] = f"""Figure {fig_n}. Normalized peak area versus analysis number. Peak Area is normalized by the working gas to account for fluctuations in the mass spec sensitivity."""
    figures[fig_n]['fig'] = figure(width=1800, height=700, x_axis_label="Analysis", y_axis_label="Normalized Peak Area (Vs)", tools="pan, box_zoom, reset, save", active_drag="box_zoom")
    figures[fig_n]['fig'].circle(Analysis[blank['index']], AreaAll_sam_norm[blank['index']], legend_label="Blank", size=blank['marker_size'], line_color='black', fill_color=blank['marker_color'])
    figures[fig_n]['fig'].triangle(Analysis[BaSO4['index']], AreaAll_sam_norm[BaSO4['index']], legend_label="BaSO4", size=BaSO4['marker_size'], line_color='black', fill_color=BaSO4['marker_color'])
    figures[fig_n]['fig'].triangle(Analysis[Ag2S['index']], AreaAll_sam_norm[Ag2S['index']], legend_label="Ag2S", size=Ag2S['marker_size'], line_color='black', fill_color=Ag2S['marker_color'])
    figures[fig_n]['fig'].triangle(Analysis[ZnS['index']], AreaAll_sam_norm[ZnS['index']], legend_label="ZnS", size=ZnS['marker_size'], line_color='black', fill_color=ZnS['marker_color'])
    # figures[fig_n]['fig'].square(Analysis[SodSul_1['index']], AreaAll_sam_norm[SodSul_1['index']], legend_label="SodSul_1", size=SodSul_1['marker_size'], line_color='black', fill_color=SodSul_1['marker_color'])
    # figures[fig_n]['fig'].diamond(Analysis[SodSul_2['index']], AreaAll_sam_norm[SodSul_2['index']], legend_label="SodSul_2", size=SodSul_2['marker_size'], line_color='black', fill_color=SodSul_2['marker_color'])
    # figures[fig_n]['fig'].diamond(Analysis[SodSul_3['index']], AreaAll_sam_norm[SodSul_3['index']], legend_label="SodSul_3", size=SodSul_3['marker_size'], line_color='black', fill_color=SodSul_3['marker_color'])
    # figures[fig_n]['fig'].diamond(Analysis[SodSul_4['index']], AreaAll_sam_norm[SodSul_4['index']], legend_label="SodSul_4", size=SodSul_4['marker_size'], line_color='black', fill_color=SodSul_4['marker_color'])
    # figures[fig_n]['fig'].diamond(Analysis[SodSul_5['index']], AreaAll_sam_norm[SodSul_5['index']], legend_label="SodSul_5", size=SodSul_5['marker_size'], line_color='black', fill_color=SodSul_5['marker_color'])
    # figures[fig_n]['fig'].diamond(Analysis[SodSul_6['index']], AreaAll_sam_norm[SodSul_6['index']], legend_label="SodSul_6", size=SodSul_6['marker_size'], line_color='black', fill_color=SodSul_6['marker_color'])
    # figures[fig_n]['fig'].diamond(Analysis[SodSul_7['index']], AreaAll_sam_norm[SodSul_7['index']], legend_label="SodSul_7", size=SodSul_7['marker_size'], line_color='black', fill_color=SodSul_7['marker_color'])
    # figures[fig_n]['fig'].diamond(Analysis[SodSul_8['index']], AreaAll_sam_norm[SodSul_8['index']], legend_label="SodSul_8", size=SodSul_8['marker_size'], line_color='black', fill_color=SodSul_8['marker_color'])
    figures[fig_n]['fig'].diamond(Analysis[sample_indices], AreaAll_sam_norm[sample_indices], legend_label="samples", size=12, line_color='black', fill_color='red')
    [figures[fig_n]['fig'].line([Analysis[i], Analysis[i]], [np.nanmin(AreaAll_sam_norm), np.nanmax(AreaAll_sam_norm)]) for i in run_index_first_row]
    figures[fig_n]['fig'].xaxis.formatter = NumeralTickFormatter(format="0")
    figures[fig_n]['fig'].xaxis.axis_label_text_font_size = font_size
    figures[fig_n]['fig'].yaxis.axis_label_text_font_size = font_size
    figures[fig_n]['fig'].xaxis.major_label_text_font_size = font_size
    figures[fig_n]['fig'].yaxis.major_label_text_font_size = font_size
    figures[fig_n]['fig'].title.text_font_size = font_size
    figures[fig_n]['fig'].legend.label_text_font_size = font_size


fig_n += 1


figures[fig_n] = {}
figures[fig_n]['cap'] = f"""Figure {fig_n}. Normalized peak area versus target sulfur quantity. Peak Area is normalized by the working gas to account for fluctuations in the mass spec sensitivity. The line of best fit is {qtycal_BaSO4['Sfit_eq_str']}."""
figures[fig_n]['fig'] = figure(width=1800, height=700, x_axis_label="Target sulfur quantity (ug)", y_axis_label="Normalized Peak Area (Vs)", tools="pan, box_zoom, reset, save", active_drag="box_zoom")
figures[fig_n]['fig'].circle(Sqty[blank['index']], AreaAll_sam_norm[blank['index']], legend_label="Blank", size=blank['marker_size'], line_color='black', fill_color=blank['marker_color'])
figures[fig_n]['fig'].circle(Sqty[qtycal_BaSO4['index']], AreaAll_sam_norm[qtycal_BaSO4['index']], legend_label="qtycal", size=qtycal_BaSO4['marker_size'], line_color='black', fill_color=qtycal_BaSO4['marker_color'])
figures[fig_n]['fig'].triangle(Sqty[BaSO4['index']], AreaAll_sam_norm[BaSO4['index']], legend_label="BaSO4", size=BaSO4['marker_size'], line_color='black', fill_color=BaSO4['marker_color'])
figures[fig_n]['fig'].triangle(Sqty[Ag2S['index']], AreaAll_sam_norm[Ag2S['index']], legend_label="Ag2S", size=Ag2S['marker_size'], line_color='black', fill_color=Ag2S['marker_color'])
figures[fig_n]['fig'].triangle(Sqty[ZnS['index']], AreaAll_sam_norm[ZnS['index']], legend_label="ZnS", size=ZnS['marker_size'], line_color='black', fill_color=ZnS['marker_color'])
# figures[fig_n]['fig'].square(Sqty[SodSul_1['index']], AreaAll_sam_norm[SodSul_1['index']], legend_label="SodSul_1", size=SodSul_1['marker_size'], line_color='black', fill_color=SodSul_1['marker_color'])
# figures[fig_n]['fig'].diamond(Sqty[SodSul_2['index']], AreaAll_sam_norm[SodSul_2['index']], legend_label="SodSul_2", size=SodSul_2['marker_size'], line_color='black', fill_color=SodSul_2['marker_color'])
# figures[fig_n]['fig'].diamond(Sqty[SodSul_3['index']], AreaAll_sam_norm[SodSul_3['index']], legend_label="SodSul_3", size=SodSul_3['marker_size'], line_color='black', fill_color=SodSul_3['marker_color'])
# figures[fig_n]['fig'].diamond(Sqty[SodSul_4['index']], AreaAll_sam_norm[SodSul_4['index']], legend_label="SodSul_4", size=SodSul_4['marker_size'], line_color='black', fill_color=SodSul_4['marker_color'])
# figures[fig_n]['fig'].diamond(Sqty[SodSul_5['index']], AreaAll_sam_norm[SodSul_5['index']], legend_label="SodSul_5", size=SodSul_5['marker_size'], line_color='black', fill_color=SodSul_5['marker_color'])
# figures[fig_n]['fig'].diamond(Sqty[SodSul_6['index']], AreaAll_sam_norm[SodSul_6['index']], legend_label="SodSul_6", size=SodSul_6['marker_size'], line_color='black', fill_color=SodSul_6['marker_color'])
# figures[fig_n]['fig'].diamond(Sqty[SodSul_7['index']], AreaAll_sam_norm[SodSul_7['index']], legend_label="SodSul_7", size=SodSul_7['marker_size'], line_color='black', fill_color=SodSul_7['marker_color'])
# figures[fig_n]['fig'].diamond(Sqty[SodSul_8['index']], AreaAll_sam_norm[SodSul_8['index']], legend_label="SodSul_8", size=SodSul_8['marker_size'], line_color='black', fill_color=SodSul_8['marker_color'])
figures[fig_n]['fig'].diamond(Sqty[sample_indices], AreaAll_sam_norm[sample_indices], legend_label="samples", size=12, line_color='black', fill_color='red')
figures[fig_n]['fig'].line(Sqty_pred[qtycal_BaSO4['index']][qtycal_BaSO4['area_sort_i']], AreaAll_sam_norm[qtycal_BaSO4['index']][qtycal_BaSO4['area_sort_i']])
figures[fig_n]['fig'].line(Sqty2_pred[qtycal_BaSO4['index']][qtycal_BaSO4['area_sort_i']], AreaAll_sam_norm[qtycal_BaSO4['index']][qtycal_BaSO4['area_sort_i']])
figures[fig_n]['fig'].xaxis.formatter = NumeralTickFormatter(format="0")
figures[fig_n]['fig'].xaxis.axis_label_text_font_size = font_size
figures[fig_n]['fig'].yaxis.axis_label_text_font_size = font_size
figures[fig_n]['fig'].xaxis.major_label_text_font_size = font_size
figures[fig_n]['fig'].yaxis.major_label_text_font_size = font_size
figures[fig_n]['fig'].title.text_font_size = font_size
figures[fig_n]['fig'].legend.label_text_font_size = font_size


fig_n += 1


figures[fig_n] = {}
figures[fig_n]['cap'] = f"""Figure {fig_n}. Calculated sulfur quantity versus analysis number."""
figures[fig_n]['fig'] = figure(width=1800, height=700, x_axis_label="Analysis", y_axis_label="Sulfur quantity (ug)", tools="pan, box_zoom, reset, save", active_drag="box_zoom")
figures[fig_n]['fig'].circle(Analysis[blank['index']], Sqty_pred[blank['index']], legend_label="Blank", size=blank['marker_size'], line_color='black', fill_color=blank['marker_color'])
figures[fig_n]['fig'].triangle(Analysis[BaSO4['index']], Sqty_pred[BaSO4['index']], legend_label="BaSO4", size=BaSO4['marker_size'], line_color='black', fill_color=BaSO4['marker_color'])
figures[fig_n]['fig'].triangle(Analysis[Ag2S['index']], Sqty_pred[Ag2S['index']], legend_label="Ag2S", size=Ag2S['marker_size'], line_color='black', fill_color=Ag2S['marker_color'])
figures[fig_n]['fig'].triangle(Analysis[ZnS['index']], Sqty_pred[ZnS['index']], legend_label="ZnS", size=ZnS['marker_size'], line_color='black', fill_color=ZnS['marker_color'])
# figures[fig_n]['fig'].square(Analysis[SodSul_1['index']], Sqty_pred[SodSul_1['index']], legend_label="SodSul_1", size=SodSul_1['marker_size'], line_color='black', fill_color=SodSul_1['marker_color'])
# figures[fig_n]['fig'].diamond(Analysis[SodSul_2['index']], Sqty_pred[SodSul_2['index']], legend_label="SodSul_2", size=SodSul_2['marker_size'], line_color='black', fill_color=SodSul_2['marker_color'])
# figures[fig_n]['fig'].diamond(Analysis[SodSul_3['index']], Sqty_pred[SodSul_3['index']], legend_label="SodSul_3", size=SodSul_3['marker_size'], line_color='black', fill_color=SodSul_3['marker_color'])
# figures[fig_n]['fig'].diamond(Analysis[SodSul_4['index']], Sqty_pred[SodSul_4['index']], legend_label="SodSul_4", size=SodSul_4['marker_size'], line_color='black', fill_color=SodSul_4['marker_color'])
# figures[fig_n]['fig'].diamond(Analysis[SodSul_5['index']], Sqty_pred[SodSul_5['index']], legend_label="SodSul_5", size=SodSul_5['marker_size'], line_color='black', fill_color=SodSul_5['marker_color'])
# figures[fig_n]['fig'].diamond(Analysis[SodSul_6['index']], Sqty_pred[SodSul_6['index']], legend_label="SodSul_6", size=SodSul_6['marker_size'], line_color='black', fill_color=SodSul_6['marker_color'])
# figures[fig_n]['fig'].diamond(Analysis[SodSul_7['index']], Sqty_pred[SodSul_7['index']], legend_label="SodSul_7", size=SodSul_7['marker_size'], line_color='black', fill_color=SodSul_7['marker_color'])
# figures[fig_n]['fig'].diamond(Analysis[SodSul_8['index']], Sqty_pred[SodSul_8['index']], legend_label="SodSul_8", size=SodSul_8['marker_size'], line_color='black', fill_color=SodSul_8['marker_color'])
figures[fig_n]['fig'].diamond(Analysis[sample_indices], Sqty_pred[sample_indices], legend_label="samples", size=12, line_color='black', fill_color='red')
[figures[fig_n]['fig'].line([Analysis[i], Analysis[i]], [np.nanmin(Sqty_pred), np.nanmax(Sqty_pred)]) for i in run_index_first_row]
figures[fig_n]['fig'].xaxis.formatter = NumeralTickFormatter(format="0")
figures[fig_n]['fig'].xaxis.axis_label_text_font_size = font_size
figures[fig_n]['fig'].yaxis.axis_label_text_font_size = font_size
figures[fig_n]['fig'].xaxis.major_label_text_font_size = font_size
figures[fig_n]['fig'].yaxis.major_label_text_font_size = font_size
figures[fig_n]['fig'].title.text_font_size = font_size
figures[fig_n]['fig'].legend.label_text_font_size = font_size


fig_n += 1


figures[fig_n] = {}
figures[fig_n]['cap'] = f"""Figure {fig_n}. Yield."""
figures[fig_n]['fig'] = figure(width=1800, height=700, x_axis_label="Analysis", y_axis_label="Yield (Vs / ug)", tools="pan, box_zoom, reset, save", active_drag="box_zoom")
figures[fig_n]['fig'].triangle(Analysis[BaSO4['index']], AreaAll_sam_norm[BaSO4['index']]/Sqty[BaSO4['index']], legend_label="BaSO4", size=BaSO4['marker_size'], line_color='black', fill_color=BaSO4['marker_color'])
figures[fig_n]['fig'].triangle(Analysis[Ag2S['index']], AreaAll_sam_norm[Ag2S['index']]/Sqty[Ag2S['index']], legend_label="Ag2S", size=Ag2S['marker_size'], line_color='black', fill_color=Ag2S['marker_color'])
figures[fig_n]['fig'].triangle(Analysis[ZnS['index']], AreaAll_sam_norm[ZnS['index']]/Sqty[ZnS['index']], legend_label="ZnS", size=ZnS['marker_size'], line_color='black', fill_color=ZnS['marker_color'])
# figures[fig_n]['fig'].square(Analysis[SodSul_1['index']], AreaAll_sam_norm[SodSul_1['index']]/Sqty[SodSul_1['index']], legend_label="SodSul_1", size=SodSul_1['marker_size'], line_color='black', fill_color=SodSul_1['marker_color'])
# figures[fig_n]['fig'].diamond(Analysis[SodSul_2['index']], AreaAll_sam_norm[SodSul_2['index']]/Sqty[SodSul_2['index']], legend_label="SodSul_2", size=SodSul_2['marker_size'], line_color='black', fill_color=SodSul_2['marker_color'])
# figures[fig_n]['fig'].diamond(Analysis[SodSul_3['index']], AreaAll_sam_norm[SodSul_3['index']]/Sqty[SodSul_3['index']], legend_label="SodSul_3", size=SodSul_3['marker_size'], line_color='black', fill_color=SodSul_3['marker_color'])
# figures[fig_n]['fig'].diamond(Analysis[SodSul_4['index']], AreaAll_sam_norm[SodSul_4['index']]/Sqty[SodSul_4['index']], legend_label="SodSul_4", size=SodSul_4['marker_size'], line_color='black', fill_color=SodSul_4['marker_color'])
# figures[fig_n]['fig'].diamond(Analysis[SodSul_5['index']], AreaAll_sam_norm[SodSul_5['index']]/Sqty[SodSul_5['index']], legend_label="SodSul_5", size=SodSul_5['marker_size'], line_color='black', fill_color=SodSul_5['marker_color'])
# figures[fig_n]['fig'].diamond(Analysis[SodSul_6['index']], AreaAll_sam_norm[SodSul_6['index']]/Sqty[SodSul_6['index']], legend_label="SodSul_6", size=SodSul_6['marker_size'], line_color='black', fill_color=SodSul_6['marker_color'])
# figures[fig_n]['fig'].diamond(Analysis[SodSul_7['index']], AreaAll_sam_norm[SodSul_7['index']]/Sqty[SodSul_7['index']], legend_label="SodSul_7", size=SodSul_7['marker_size'], line_color='black', fill_color=SodSul_7['marker_color'])
# figures[fig_n]['fig'].diamond(Analysis[SodSul_8['index']], AreaAll_sam_norm[SodSul_8['index']]/Sqty[SodSul_8['index']], legend_label="SodSul_8", size=SodSul_8['marker_size'], line_color='black', fill_color=SodSul_8['marker_color'])
[figures[fig_n]['fig'].line([Analysis[i], Analysis[i]], [0, 200]) for i in run_index_first_row]
figures[fig_n]['fig'].xaxis.formatter = NumeralTickFormatter(format="0")
figures[fig_n]['fig'].xaxis.axis_label_text_font_size = font_size
figures[fig_n]['fig'].yaxis.axis_label_text_font_size = font_size
figures[fig_n]['fig'].xaxis.major_label_text_font_size = font_size
figures[fig_n]['fig'].yaxis.major_label_text_font_size = font_size
figures[fig_n]['fig'].title.text_font_size = font_size
figures[fig_n]['fig'].legend.label_text_font_size = font_size


fig_n += 1


figures[fig_n] = {}
figures[fig_n]['cap'] = f"""Figure {fig_n}. Measured d34S vs working gas."""
figures[fig_n]['fig'] = figure(width=1800, height=700, x_axis_label="Analysis", y_axis_label="d34S measured (permil)", tools="pan, box_zoom, reset, save", active_drag="box_zoom")
figures[fig_n]['fig'].circle(Analysis[blank['index']], d34S32S_sam[blank['index']], legend_label="Blank", size=blank['marker_size'], line_color='black', fill_color=blank['marker_color'])
figures[fig_n]['fig'].triangle(Analysis[BaSO4['index']], d34S32S_sam[BaSO4['index']], legend_label="BaSO4", size=BaSO4['marker_size'], line_color='black', fill_color=BaSO4['marker_color'])
figures[fig_n]['fig'].triangle(Analysis[Ag2S['index']], d34S32S_sam[Ag2S['index']], legend_label="Ag2S", size=Ag2S['marker_size'], line_color='black', fill_color=Ag2S['marker_color'])
figures[fig_n]['fig'].triangle(Analysis[ZnS['index']], d34S32S_sam[ZnS['index']], legend_label="ZnS", size=ZnS['marker_size'], line_color='black', fill_color=ZnS['marker_color'])
# figures[fig_n]['fig'].square(Analysis[SodSul_1['index']], d34S32S_sam[SodSul_1['index']], legend_label="SodSul_1", size=SodSul_1['marker_size'], line_color='black', fill_color=SodSul_1['marker_color'])
# figures[fig_n]['fig'].diamond(Analysis[SodSul_2['index']], d34S32S_sam[SodSul_2['index']], legend_label="SodSul_2", size=SodSul_2['marker_size'], line_color='black', fill_color=SodSul_2['marker_color'])
# figures[fig_n]['fig'].diamond(Analysis[SodSul_3['index']], d34S32S_sam[SodSul_3['index']], legend_label="SodSul_3", size=SodSul_3['marker_size'], line_color='black', fill_color=SodSul_3['marker_color'])
# figures[fig_n]['fig'].diamond(Analysis[SodSul_4['index']], d34S32S_sam[SodSul_4['index']], legend_label="SodSul_4", size=SodSul_4['marker_size'], line_color='black', fill_color=SodSul_4['marker_color'])
# figures[fig_n]['fig'].diamond(Analysis[SodSul_5['index']], d34S32S_sam[SodSul_5['index']], legend_label="SodSul_5", size=SodSul_5['marker_size'], line_color='black', fill_color=SodSul_5['marker_color'])
# figures[fig_n]['fig'].diamond(Analysis[SodSul_6['index']], d34S32S_sam[SodSul_6['index']], legend_label="SodSul_6", size=SodSul_6['marker_size'], line_color='black', fill_color=SodSul_6['marker_color'])
# figures[fig_n]['fig'].diamond(Analysis[SodSul_7['index']], d34S32S_sam[SodSul_7['index']], legend_label="SodSul_7", size=SodSul_7['marker_size'], line_color='black', fill_color=SodSul_7['marker_color'])
# figures[fig_n]['fig'].diamond(Analysis[SodSul_8['index']], d34S32S_sam[SodSul_8['index']], legend_label="SodSul_8", size=SodSul_8['marker_size'], line_color='black', fill_color=SodSul_8['marker_color'])
figures[fig_n]['fig'].diamond(Analysis[sample_indices], d34S32S_sam[sample_indices], legend_label="samples", size=12, line_color='black', fill_color='red')
[figures[fig_n]['fig'].line([Analysis[i], Analysis[i]], [np.nanmin(d34S32S_sam), np.nanmax(d34S32S_sam)]) for i in run_index_first_row]
figures[fig_n]['fig'].xaxis.formatter = NumeralTickFormatter(format="0")
figures[fig_n]['fig'].xaxis.axis_label_text_font_size = font_size
figures[fig_n]['fig'].yaxis.axis_label_text_font_size = font_size
figures[fig_n]['fig'].xaxis.major_label_text_font_size = font_size
figures[fig_n]['fig'].yaxis.major_label_text_font_size = font_size
figures[fig_n]['fig'].title.text_font_size = font_size
figures[fig_n]['fig'].legend.label_text_font_size = font_size



fig_n += 1



figures[fig_n] = {}
figures[fig_n]['cap'] = f"""Figure {fig_n}. d34S residual vs calculated Sulfur quantity."""
figures[fig_n]['fig'] = figure(width=1800, height=700, x_axis_label="Sqty_pred", y_axis_label="residual d34S (permil)", tools="pan, box_zoom, reset, save", active_drag="box_zoom")
figures[fig_n]['fig'].triangle(Sqty_pred[BaSO4['index']], d34S32S_sam[BaSO4['index']] - np.nanmean(d34S32S_sam[BaSO4['index']]), legend_label="BaSO4", size=BaSO4['marker_size'], line_color='black', fill_color=BaSO4['marker_color'])
figures[fig_n]['fig'].triangle(Sqty_pred[Ag2S['index']], d34S32S_sam[Ag2S['index']] - np.nanmean(d34S32S_sam[Ag2S['index']]), legend_label="Ag2S", size=Ag2S['marker_size'], line_color='black', fill_color=Ag2S['marker_color'])
figures[fig_n]['fig'].triangle(Sqty_pred[ZnS['index']], d34S32S_sam[ZnS['index']] - np.nanmean(d34S32S_sam[ZnS['index']]), legend_label="ZnS", size=ZnS['marker_size'], line_color='black', fill_color=ZnS['marker_color'])
# figures[fig_n]['fig'].square(Sqty_pred[SodSul_1['index']], d34S32S_sam[SodSul_1['index']] - np.nanmean(d34S32S_sam[SodSul_1['index']]), legend_label="SodSul_1", size=SodSul_1['marker_size'], line_color='black', fill_color=SodSul_1['marker_color'])
# figures[fig_n]['fig'].diamond(Sqty_pred[SodSul_2['index']], d34S32S_sam[SodSul_2['index']] - np.nanmean(d34S32S_sam[SodSul_2['index']]), legend_label="SodSul_2", size=SodSul_2['marker_size'], line_color='black', fill_color=SodSul_2['marker_color'])
# figures[fig_n]['fig'].diamond(Sqty_pred[SodSul_3['index']], d34S32S_sam[SodSul_3['index']] - np.nanmean(d34S32S_sam[SodSul_3['index']]), legend_label="SodSul_3", size=SodSul_3['marker_size'], line_color='black', fill_color=SodSul_3['marker_color'])
# figures[fig_n]['fig'].diamond(Sqty_pred[SodSul_4['index']], d34S32S_sam[SodSul_4['index']] - np.nanmean(d34S32S_sam[SodSul_4['index']]), legend_label="SodSul_4", size=SodSul_4['marker_size'], line_color='black', fill_color=SodSul_4['marker_color'])
# figures[fig_n]['fig'].diamond(Sqty_pred[SodSul_5['index']], d34S32S_sam[SodSul_5['index']] - np.nanmean(d34S32S_sam[SodSul_5['index']]), legend_label="SodSul_5", size=SodSul_5['marker_size'], line_color='black', fill_color=SodSul_5['marker_color'])
# figures[fig_n]['fig'].diamond(Sqty_pred[SodSul_6['index']], d34S32S_sam[SodSul_6['index']] - np.nanmean(d34S32S_sam[SodSul_6['index']]), legend_label="SodSul_6", size=SodSul_6['marker_size'], line_color='black', fill_color=SodSul_6['marker_color'])
# figures[fig_n]['fig'].diamond(Sqty_pred[SodSul_7['index']], d34S32S_sam[SodSul_7['index']] - np.nanmean(d34S32S_sam[SodSul_7['index']]), legend_label="SodSul_7", size=SodSul_7['marker_size'], line_color='black', fill_color=SodSul_7['marker_color'])
# figures[fig_n]['fig'].diamond(Sqty_pred[SodSul_8['index']], d34S32S_sam[SodSul_8['index']] - np.nanmean(d34S32S_sam[SodSul_8['index']]), legend_label="SodSul_8", size=SodSul_8['marker_size'], line_color='black', fill_color=SodSul_8['marker_color'])
figures[fig_n]['fig'].xaxis.formatter = NumeralTickFormatter(format="0")
figures[fig_n]['fig'].xaxis.axis_label_text_font_size = font_size
figures[fig_n]['fig'].yaxis.axis_label_text_font_size = font_size
figures[fig_n]['fig'].xaxis.major_label_text_font_size = font_size
figures[fig_n]['fig'].yaxis.major_label_text_font_size = font_size
figures[fig_n]['fig'].title.text_font_size = font_size
figures[fig_n]['fig'].legend.label_text_font_size = font_size


fig_n += 1



figures[fig_n] = {}
figures[fig_n]['cap'] = f"""Figure {fig_n}. d34S residual vs Peak Area."""
figures[fig_n]['fig'] = figure(width=1800, height=700, x_axis_label="AreaAll_sam", y_axis_label="residual d34S (permil)", tools="pan, box_zoom, reset, save", active_drag="box_zoom")
figures[fig_n]['fig'].triangle(AreaAll_sam[BaSO4['index']], d34S32S_sam[BaSO4['index']] - np.nanmean(d34S32S_sam[BaSO4['index']]), legend_label="BaSO4", size=BaSO4['marker_size'], line_color='black', fill_color=BaSO4['marker_color'])
figures[fig_n]['fig'].triangle(AreaAll_sam[Ag2S['index']], d34S32S_sam[Ag2S['index']] - np.nanmean(d34S32S_sam[Ag2S['index']]), legend_label="Ag2S", size=Ag2S['marker_size'], line_color='black', fill_color=Ag2S['marker_color'])
figures[fig_n]['fig'].triangle(AreaAll_sam[ZnS['index']], d34S32S_sam[ZnS['index']] - np.nanmean(d34S32S_sam[ZnS['index']]), legend_label="ZnS", size=ZnS['marker_size'], line_color='black', fill_color=ZnS['marker_color'])
# figures[fig_n]['fig'].square(Sqty_pred[SodSul_1['index']], d34S32S_sam[SodSul_1['index']] - np.nanmean(d34S32S_sam[SodSul_1['index']]), legend_label="SodSul_1", size=SodSul_1['marker_size'], line_color='black', fill_color=SodSul_1['marker_color'])
# figures[fig_n]['fig'].diamond(Sqty_pred[SodSul_2['index']], d34S32S_sam[SodSul_2['index']] - np.nanmean(d34S32S_sam[SodSul_2['index']]), legend_label="SodSul_2", size=SodSul_2['marker_size'], line_color='black', fill_color=SodSul_2['marker_color'])
# figures[fig_n]['fig'].diamond(Sqty_pred[SodSul_3['index']], d34S32S_sam[SodSul_3['index']] - np.nanmean(d34S32S_sam[SodSul_3['index']]), legend_label="SodSul_3", size=SodSul_3['marker_size'], line_color='black', fill_color=SodSul_3['marker_color'])
# figures[fig_n]['fig'].diamond(Sqty_pred[SodSul_4['index']], d34S32S_sam[SodSul_4['index']] - np.nanmean(d34S32S_sam[SodSul_4['index']]), legend_label="SodSul_4", size=SodSul_4['marker_size'], line_color='black', fill_color=SodSul_4['marker_color'])
# figures[fig_n]['fig'].diamond(Sqty_pred[SodSul_5['index']], d34S32S_sam[SodSul_5['index']] - np.nanmean(d34S32S_sam[SodSul_5['index']]), legend_label="SodSul_5", size=SodSul_5['marker_size'], line_color='black', fill_color=SodSul_5['marker_color'])
# figures[fig_n]['fig'].diamond(Sqty_pred[SodSul_6['index']], d34S32S_sam[SodSul_6['index']] - np.nanmean(d34S32S_sam[SodSul_6['index']]), legend_label="SodSul_6", size=SodSul_6['marker_size'], line_color='black', fill_color=SodSul_6['marker_color'])
# figures[fig_n]['fig'].diamond(Sqty_pred[SodSul_7['index']], d34S32S_sam[SodSul_7['index']] - np.nanmean(d34S32S_sam[SodSul_7['index']]), legend_label="SodSul_7", size=SodSul_7['marker_size'], line_color='black', fill_color=SodSul_7['marker_color'])
# figures[fig_n]['fig'].diamond(Sqty_pred[SodSul_8['index']], d34S32S_sam[SodSul_8['index']] - np.nanmean(d34S32S_sam[SodSul_8['index']]), legend_label="SodSul_8", size=SodSul_8['marker_size'], line_color='black', fill_color=SodSul_8['marker_color'])
figures[fig_n]['fig'].xaxis.formatter = NumeralTickFormatter(format="0")
figures[fig_n]['fig'].xaxis.axis_label_text_font_size = font_size
figures[fig_n]['fig'].yaxis.axis_label_text_font_size = font_size
figures[fig_n]['fig'].xaxis.major_label_text_font_size = font_size
figures[fig_n]['fig'].yaxis.major_label_text_font_size = font_size
figures[fig_n]['fig'].title.text_font_size = font_size
figures[fig_n]['fig'].legend.label_text_font_size = font_size

fig_n += 1


figures[fig_n] = {}
figures[fig_n]['cap'] = f"""Figure {fig_n}. d34S vcdt vs d34S measured."""
figures[fig_n]['fig'] = figure(width=1800, height=700, x_axis_label="d34S measured", y_axis_label="d34S VCDT (permil)", tools="pan, box_zoom, reset, save", active_drag="box_zoom")
figures[fig_n]['fig'].triangle(d34S32S_sam[BaSO4['index']], d34S_VCDT[BaSO4['index']], legend_label="BaSO4", size=BaSO4['marker_size'], line_color='black', fill_color=BaSO4['marker_color'])
figures[fig_n]['fig'].triangle(d34S32S_sam[Ag2S['index']], d34S_VCDT[Ag2S['index']], legend_label="Ag2S", size=Ag2S['marker_size'], line_color='black', fill_color=Ag2S['marker_color'])
figures[fig_n]['fig'].triangle(d34S32S_sam[ZnS['index']], d34S_VCDT[ZnS['index']], legend_label="ZnS", size=ZnS['marker_size'], line_color='black', fill_color=ZnS['marker_color'])
# figures[fig_n]['fig'].square(d34S32S_sam[SodSul_1['index']], d34S_VCDT[SodSul_1['index']], legend_label="SodSul_1", size=SodSul_1['marker_size'], line_color='black', fill_color=SodSul_1['marker_color'])
# figures[fig_n]['fig'].diamond(d34S32S_sam[SodSul_2['index']], d34S_VCDT[SodSul_2['index']], legend_label="SodSul_2", size=SodSul_2['marker_size'], line_color='black', fill_color=SodSul_2['marker_color'])
# figures[fig_n]['fig'].diamond(d34S32S_sam[SodSul_3['index']], d34S_VCDT[SodSul_3['index']], legend_label="SodSul_3", size=SodSul_3['marker_size'], line_color='black', fill_color=SodSul_3['marker_color'])
# figures[fig_n]['fig'].diamond(d34S32S_sam[SodSul_4['index']], d34S_VCDT[SodSul_4['index']], legend_label="SodSul_4", size=SodSul_4['marker_size'], line_color='black', fill_color=SodSul_4['marker_color'])
# figures[fig_n]['fig'].diamond(d34S32S_sam[SodSul_5['index']], d34S_VCDT[SodSul_5['index']], legend_label="SodSul_5", size=SodSul_5['marker_size'], line_color='black', fill_color=SodSul_5['marker_color'])
# figures[fig_n]['fig'].diamond(d34S32S_sam[SodSul_6['index']], d34S_VCDT[SodSul_6['index']], legend_label="SodSul_6", size=SodSul_6['marker_size'], line_color='black', fill_color=SodSul_6['marker_color'])
# figures[fig_n]['fig'].diamond(d34S32S_sam[SodSul_7['index']], d34S_VCDT[SodSul_7['index']], legend_label="SodSul_7", size=SodSul_7['marker_size'], line_color='black', fill_color=SodSul_7['marker_color'])
# figures[fig_n]['fig'].diamond(d34S32S_sam[SodSul_8['index']], d34S_VCDT[SodSul_8['index']], legend_label="SodSul_8", size=SodSul_8['marker_size'], line_color='black', fill_color=SodSul_8['marker_color'])
figures[fig_n]['fig'].diamond(d34S32S_sam[sample_indices], d34S_VCDT[sample_indices], legend_label="samples", size=12, line_color='black', fill_color='red')
figures[fig_n]['fig'].xaxis.formatter = NumeralTickFormatter(format="0")
figures[fig_n]['fig'].legend.location = 'left'
figures[fig_n]['fig'].xaxis.axis_label_text_font_size = font_size
figures[fig_n]['fig'].yaxis.axis_label_text_font_size = font_size
figures[fig_n]['fig'].xaxis.major_label_text_font_size = font_size
figures[fig_n]['fig'].yaxis.major_label_text_font_size = font_size
figures[fig_n]['fig'].title.text_font_size = font_size
figures[fig_n]['fig'].legend.label_text_font_size = font_size


fig_n += 1






# -------------------- export summary data file --------------------
print(f'\n    Creating summary data file.')
summary_data_filename = f'shrekS_summary_data.csv'
summary_data_file = os.path.join(method_directory, report_directory, 'data/', summary_data_filename)
summary_file_headers = ['Sample ID', 'Date', 'Analysis Number', 'Mass (mg)', 'Peak Area (Vs)', 'Sulfur amount (ug)', 'd34S vs VCDT (permil)']
data_to_write = '[Identifier1[ii], Date[ii], int(Analysis[ii]), Amount[ii], AreaAll_sam_norm[ii], Sqty[ii], round(d34S_VCDT[ii], 2)]'
data_to_write = str(data_to_write).replace("'", "")
with open(summary_data_file, 'w', newline='') as csvfile:
    datawriter = csv.writer(csvfile, quoting=csv.QUOTE_MINIMAL)
    datawriter.writerow(summary_file_headers)
    for ii in knowns_indices:
        datawriter.writerow(eval(data_to_write))
    datawriter.writerow('\n')
    for ii in sample_indices:
        datawriter.writerow(eval(data_to_write))


# -------------------- make html summary report --------------------

# copy report files
shutil.copy2(os.path.join(python_directory, 'py_report_style.css'), os.path.join(method_directory, report_directory))
[shutil.copy2(os.path.join(python_directory, script), os.path.join(method_directory, report_directory, f"python/{script}_REPORT_COPY")) for script in python_scripts]


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
    <body class="entire_page">\n

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
            <tr><td>Total number of standards analyzed <sup>*</sup></td><td>{len(knowns_indices)}</td></tr>
            <tr><td>Total number of samples analyzed</td><td>{len(sample_indices)}</td></tr>
            <tr><td><br></td></tr>
            <tr><td>Number of <a href="#excluded">excluded analyses</a></td><td>{len(trust0_indices)}</td></tr>
        </table>
        <sup>*</sup><small> - Reference material analyses are accumulated from many individual runs and projects into a single reference
                              frame until the instrumentation shifts and a new reference frame must be considered.</small>
    </div>

    <h2>Figures</h2>
        <div class="text-indent">"""


figure_block = [f"""<div class="clear-both">{file_html(figures[i]['fig'], CDN)}{figures[i]['cap']}<hr></div>""" for i in figures.keys()]


excluded_analyses_block = str([f"<tr><td>{original_data['Analysis'][i]}</td><td>{original_data['Identifier1'][i]}</td><td>{original_data['notes'][i]}</td></tr>" for i in trust0_indices]).replace("[","").replace("'","").replace("]","").replace(", ","")

excluded_analysis = f"""</div>
    \n<h2 id="excluded">Excluded analyses</h2>
    <div class="text-indent"><p>This is a list of analyses that have been excluded from further data processing showing
    the sample ID, the analysis number, and the reason for exclusion. If you want to included select analyses back into
    data processing, open the log file, find the sample of interest, and change trust to 1, then rerun the script.</p>

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


log_summary_page = os.path.join(method_directory, report_directory, 'report.html')
with open(log_summary_page, 'w') as html_page:

    html_page.write(header)

    [html_page.write(i) for i in figure_block]

    html_page.write(excluded_analysis)

    html_page.write(footer)




# create zip file of entire report directory
# try:
    # shutil.make_archive('report', 'zip', os.path.join(method_directory, report_directory))
# except PermissionError as e:
    # print(f"PermissionError: {e}")
# except Exception as e:
    # print(f"An error occurred: {e}")
    
# if os.path.exists(os.path.join(method_directory, report_directory, 'report.zip')):
    # os.remove(os.path.join(method_directory, report_directory, 'report.zip'))
# shutil.move(os.path.join('report.zip'), os.path.join(method_directory, report_directory))



webbrowser.open(log_summary_page)
