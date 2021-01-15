#!/usr/bin/env python

import matlab.engine
import os
import re
import argparse


"""
Python wrapper for passing command line arguments to the matlab changepoint detection script
"""

parser = argparse.ArgumentParser()
parser.add_argument('--path', type=str, help='path to a directory of per base coverage files')
parser.add_argument('--var_start', type=int, help='variant start locus, only base value not including chromosome')
parser.add_argument('--var_stop', type=int, help='variant stop locus, only base value not including chromosome')
parser.add_argument('--extension_flank', type=float, help='decimal value to multiply the variant length by and flank the start and stop by this amount. The changepoint detection script needs some flanking region to work correctly. Example: enter 0.5 and putativeVarStart will be var_start-0.5*var_len and putativeVarStop will be var_stop+0.5*var_len') 
parser.add_argument('--output_dir', type=str, help='output for matlab script')
args = parser.parse_args()

path, var_start, var_stop, extension_flank, output_dir = args.path, args.var_start, args.var_stop, args.extension_flank, args.output_dir 

eng = matlab.engine.start_matlab()

var_len = abs(var_stop-var_start)
flankingStart = int( var_start -(extension_flank*var_len) )
flankingStop = int( var_stop + (extension_flank*var_len) )

eng.SignalDetection_CL_args(flankingStart, flankingStop, path+'/*', output_dir, nargout=0)

eng.quit()

