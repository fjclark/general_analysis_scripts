#!/usr/bin/env python
# coding: utf-8

# Takes desired run number as string (1-100) and returns all energy
# contributions to the free energy of binding along with associated error calculated
# by propogation of the MBAR errors (a large underestimation of the uncertainty
# calculated from multiple runs). Must be run in the directory containing the free 
# and bound directories.

# In[1]:


import numpy as np
import sys
import re
import os


# ### Arguments

# In[2]:


# 01 to allow for up to 99 runs

run_no = '0'


# input desired run number as a string
input_no = sys.argv[1]

if int(input_no) < 10:
    run_no += input_no
else:
    run_no = input_no


# ### Functions to Read Files and Extract Energies

# In[3]:


def getEnergyMBAR(input_file):
    '''
    Finds MBAR free energy from .dat input file.
    '''
    with open(input_file,'r') as file:
        read_file = file.readlines()
    
    energy =''
    standard_dev =''

    for i, line in enumerate(read_file):
        if line.startswith('#MBAR free energy difference in kcal/mol:'):
            energy=read_file[i+1].split()[0][:-1]
            standard_dev=read_file[i+1].split()[1][:-1]
        
    return (energy, standard_dev)


# In[4]:


def getEnergyTI(input_file):
    '''
    Finds TI free energy from .dat input file.
    '''
    with open(input_file,'r') as file:
        read_file = file.readlines()
    
    energy =''

    for i, line in enumerate(read_file):
        if line.startswith('#TI free energy difference in kcal/mol:'):
            energy=read_file[i+1].split()[0][:-1]
        
    return energy


# In[5]:


def getLJcorr(input_file):
    '''Finds LJ correction from .dat file'''
    
    with open(input_file,'r') as file:
        read_file = file.readlines()
    
    correction = read_file[0].split()[2]
    standard_dev= read_file[0].split()[4]
        
    return (correction,standard_dev)


# In[6]:


def getStdState(input_file):
    '''Finds standard state correction from .dat file and returns correction as string'''
    
    with open(input_file,'r') as file:
        read_file = file.readlines()
        
    energy=''
    
    for i, line in enumerate(read_file):
        if line.startswith('Free energy change upon removing the restraint and applying standard state conditions'):
            energy=read_file[i].split()[-2]
            
    return energy


# ### Paths to Files

# Find dir names

dir_names = os.listdir("bound")
r = re.compile(f"run0{run_no}")
dir_name = list(filter(r.match, dir_names))[0]

energy_paths = {'free_discharge':f'free/{dir_name}/discharge/output/freenrg-MBAR-p-83-overlap.dat',
         'free_vanish':f'free/{dir_name}/vanish/output/freenrg-MBAR-p-83-overlap.dat',
         'bound_discharge':f'bound/{dir_name}/discharge/output/freenrg-MBAR-p-83-overlap.dat',
         'bound_vanish':f'bound/{dir_name}/vanish/output/freenrg-MBAR-p-83-overlap.dat'}


# ### Get Energies

# In[10]:


energies_mbar = {}

for section in energy_paths:
    energies_mbar[section]=float(getEnergyMBAR(energy_paths[section])[0])


# In[11]:


sd_mbar = {}

for section in energy_paths:
    sd_mbar[section]=float(getEnergyMBAR(energy_paths[section])[1])


# In[12]:


energies_ti = {}

for section in energy_paths:
    energies_ti[section]=float(getEnergyTI(energy_paths[section]))


# ### Get Corrections

# In[13]:


free_lj_corr_str, free_sd_lj_corr_str = getLJcorr(f'free/{dir_name}/vanish/output/freenrg-LJCOR.dat')[0:2]
free_lj_corr = float(free_lj_corr_str)
free_sd_lj_corr = float(free_sd_lj_corr_str)


# In[14]:


bound_lj_corr_str, bound_sd_lj_corr_str = getLJcorr(f'bound/{dir_name}/vanish/output/freenrg-LJCOR.dat')[0:2]
bound_lj_corr = float(bound_lj_corr_str)
bound_sd_lj_corr = float(bound_sd_lj_corr_str)


# In[15]:


std_state_corr = float(getStdState(f'bound/{dir_name}/vanish/output/standard-state-s-1-b-4-d-0.25-o-6.dat'))


# ### Calculate Output

# In[16]:


dg_bind_mbar=energies_mbar['free_discharge']+energies_mbar['free_vanish']+free_lj_corr-energies_mbar['bound_discharge']-energies_mbar['bound_vanish']-bound_lj_corr-std_state_corr


# In[17]:


dg_bind_ti=energies_ti['free_discharge']+energies_ti['free_vanish']+free_lj_corr-energies_ti['bound_discharge']-energies_ti['bound_vanish']-bound_lj_corr-std_state_corr


# In[19]:


var = 0

for stage in sd_mbar:
    var+= sd_mbar[stage]**2
    
var+= free_sd_lj_corr**2
var+= bound_sd_lj_corr**2

sd_dg_bind_mbar = np.sqrt(var)


dg_elec_mbar=energies_mbar['free_discharge']-energies_mbar['bound_discharge']

sd_dg_elec_mbar=np.sqrt(sd_mbar['bound_discharge']**2+sd_mbar['free_discharge']**2)

dg_elec_ti=energies_ti['free_discharge']-energies_ti['bound_discharge']

dg_lj_mbar=energies_mbar['free_vanish']+free_lj_corr-energies_mbar['bound_vanish']-bound_lj_corr

sd_dg_lj_mbar=np.sqrt(sd_mbar['bound_vanish']**2+sd_mbar['free_vanish']**2+free_sd_lj_corr**2+bound_sd_lj_corr**2)

dg_lj_ti=energies_ti['free_vanish']+free_lj_corr-energies_ti['bound_vanish']-bound_lj_corr


# ### Write Output

# In[22]:


for section in energy_paths:
    if section == 'bound_vanish':
        print('Section: ',section)
        print('MBAR estimate:', energies_mbar[section],'+/-',sd_mbar[section],'kcal/mol')
        print('TI estimate:', energies_ti[section],'kcal/mol')
        print('LJ correction: ', bound_lj_corr,'+/-',bound_sd_lj_corr,'kcal/mol')
        print('Standard state correction ', std_state_corr, 'kcal/mol','\n')
        
    elif section == 'free_vanish':
        print('Section: ',section)
        print('MBAR estimate:', energies_mbar[section],'+/-',sd_mbar[section],'kcal/mol')
        print('TI estimate:', energies_ti[section],'kcal/mol')
        print('LJ correction: ', free_lj_corr,'+/-',free_sd_lj_corr,'kcal/mol','\n')

    else:
        print('Section: ',section)
        print('MBAR estimate:', energies_mbar[section],'+/-',sd_mbar[section],'kcal/mol')
        print('TI estimate:', energies_ti[section],'kcal/mol','\n')
    
    
print('###########################################################################')

print('dg_electrostatic from MBAR: ',dg_elec_mbar,' +/-', sd_dg_elec_mbar,'kcal/mol')
print('dg_electrostatic from TI: ',dg_elec_ti,'kcal/mol')    

print('dg_LJ from MBAR: ',dg_lj_mbar,' +/-', sd_dg_elec_mbar,'kcal/mol')
print('dg_LJ from TI: ',dg_lj_ti,'kcal/mol','\n')  
    
print('###########################################################################')

print('dg_bind from MBAR: ',dg_bind_mbar,' +/-', sd_dg_bind_mbar,'kcal/mol')
print('dg_bind from TI: ',dg_bind_ti,'kcal/mol','\n')

print('###########################################################################\n')

