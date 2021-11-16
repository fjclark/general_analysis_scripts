#!/usr/bin/env python
# coding: utf-8

import sys

# First argument is file, follwed by desired number of lambda windows

#generate list of desired lambdas
no_windows = int(sys.argv[2])
lam_vals = []
for i in range(no_windows):
    lam_vals.append(i*(1/(no_windows-1)))

#generate lambda line in correct format
lambda_line = 'lambda array = '
for lam_val in lam_vals:
    lambda_line += f' {lam_val:.3f},'
lambda_line = lambda_line[:-1]
lambda_line +='\n'

print(f'\n#############      Writing lambda line      ###############\n{lambda_line}')

#read from file
file = sys.argv[1]
lam_line_idx=0
with open(file,'r') as config_file:
    lines = config_file.readlines()
    for i, l in enumerate(lines):
        if l.startswith('lambda array ='):
            lam_line_idx=i

lines[lam_line_idx] = lambda_line

#write lines
with open(file,'w') as config_file:
    config_file.writelines(lines)


