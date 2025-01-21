#!/usr/bin/env python3
"""
Library of functions used by the IsoLab shrekS_* suite of python scripts.

    version 3.0 - 2023.12.04 - Starting off with version 3 to match the other shrekS scripts. This library did not exist
        prior to this version.
    version 4.0 - 2024.03.11 - added get_path function
    version 4.1 - 2024.05.11 - removed get path function in favor of the one in isolab_lib
    version 4.2 - 2024.06.09 - changed flag to trust, cleaned up and made similar to CN version of this library
"""


__authors__ = "Andy Schauer, Ursula Jongebloed"
__email__ = "aschauer@uw.edu"
__last_modified__ = "2024-06-09"
__version__ = "4.2"
__copyright__ = "Copyright 2025, Andy Schauer"
__license__ = "Apache 2.0"
__acknowledgements__ = "Alli Moon, Drew Pronovost"


refmat_list = ['Ag2S', 'BaSO4', 'ZnS']
knowns_list = refmat_list[:] + ['blank', 'qtycal_BaSO4', 'qtycal_Na2SO4', 'void', 'zero', 'Ag2SO4', 'NIST1547', 'SodSul_1', 'SodSul_2', 'SodSul_3', 'SodSul_4', 'SodSul_5', 'SodSul_6', 'SodSul_7', 'SodSul_8']


meta_headers = ['Identifier1', 'Analysis', 'Amount', 'Date', 'Time', 'Comment', 'Identifier2', 'Preparation', 'Method']
S_headers = ['Ampl64', 'Ampl66', 'AreaAll', 'Area64', 'Area66', 'BGD64', 'BGD66', 'R34S32S', 'R66SO264SO2', 'd34S32S',
             'PeakNr', 'Start', 'Width', 'End']
supp_headers = ['Sqty', 'file', 'trust', 'notes', 'peak_center', 'rows_per_sample', 'pyversions', 'empty']


numlist = ['Analysis', 'Amount', 'Ampl64_wg', 'Ampl66_wg', 'AreaAll_wg', 'Area64_wg', 'Area66_wg',
           'BGD64_wg', 'BGD66_wg', 'R34S32S_wg', 'R66SO264SO2_wg', 'd34S32S_wg', 'PeakNr_wg', 'Start_wg',
           'Width_wg', 'End_wg', 'Ampl64_sam', 'Ampl66_sam', 'AreaAll_sam', 'Area64_sam', 'Area66_sam',
           'BGD64_sam', 'BGD66_sam', 'R34S32S_sam', 'R66SO264SO2_sam', 'd34S32S_sam', 'PeakNr_sam',
           'Start_sam', 'Width_sam', 'End_sam', 'Sqty', 'trust', 'peak_center', 'rows_per_sample']


blank = {'names': ['blank'],
         'material': 'a tin capsule packed with WO3 and Sn powder'}
qtycal_Na2SO4 = {'names': ['qtycal_Na2SO4', 'qtycal'],
          'material': 'Na2SO4',
          'fractionS': 0.22574,
          'notes': 'sodium sulfate solution based amounts'}
qtycal_BaSO4 = {'names': ['qtycal_BaSO4'],
          'material': 'BaSO4',
          'fractionS': 0.13739,
          'notes': 'barium sulfate solid material'}
void = {'names': ['void'],
         'material': None,
         'notes': 'no material dropped into EA'}
zero = {'names': ['zero'],
        'material': 'reference gas peaks treated as unknowns'}


meta_data = {}
S_wg_data = {}
S_sam_data = {}
supp_data = {}

shrekS_analysis_log_headers = []
shrekS_analysis_log_headers.extend(meta_headers)
shrekS_analysis_log_headers.extend([f"{i}_wg" for i in S_headers])
shrekS_analysis_log_headers.extend([f"{i}_sam" for i in S_headers])
shrekS_analysis_log_headers.extend(supp_headers)

data_to_write = []
data_to_write.extend([f"meta_data['{i}'][ii]" for i in meta_headers])
data_to_write.extend([f"S_wg_data['{i}'][ii]" for i in S_headers])
data_to_write.extend([f"S_sam_data['{i}'][ii]" for i in S_headers])
data_to_write.extend([f"supp_data['{i}'][ii]" for i in supp_headers])
data_to_write = str(data_to_write).replace("\"", "")
