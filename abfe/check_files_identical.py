#!/usr/bin/python

import sys

with open(sys.argv[1],'r') as buffer:
    file1=buffer.readlines()

with open(sys.argv[2],'r') as buffer:
    file2=buffer.readlines()

for i in range(len(file1)):
    print(file1[i]==file2[i])
