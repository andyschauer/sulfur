#!/usr/bin/env python3
"""
Library of functions used by the IsoLab shrekS_* suite of python scripts.

    version 3.0 - 2023.12.04 - Starting off with version 3 to match the other shrekS scripts. This library did not exist
        prior to this version.
"""

__authors__ = "Andy Schauer, Ursula Jongebloed"
__email__ = "aschauer@uw.edu"
__last_modified__ = "2023-12-04"
__version__ = "3.0"
__copyright__ = "Copyright 2023, Andy Schauer"
__license__ = "Apache 2.0"
__acknowledgements__ = "Alli Moon, Drew Pronovost"


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
