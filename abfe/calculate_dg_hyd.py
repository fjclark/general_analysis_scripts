#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np


# ### Functions to Read Files and Extract Energies

# In[2]:


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


# In[3]:


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


# In[4]:


def getLJcorr(input_file):
    '''Finds LJ correction from .dat file'''
    
    with open(input_file,'r') as file:
        read_file = file.readlines()
    
    correction = read_file[0].split()[2]
    standard_dev= read_file[0].split()[4]
        
    return (correction,standard_dev)


# ### Paths to Files

# In[5]:


paths = {'vacuum_discharge':'vacuum/discharge/output/freenrg-MBAR-p-83-overlap.dat',
         'vacuum_vanish':'vacuum/vanish/output/freenrg-MBAR-p-83-overlap.dat',
         'free_discharge':'free/discharge/output/freenrg-MBAR-p-83-overlap.dat',
         'free_vanish':'free/vanish/output/freenrg-MBAR-p-83-overlap.dat'}


# ### Get Energies

# In[6]:


energies_mbar = {}

for section in paths:
    energies_mbar[section]=float(getEnergyMBAR(paths[section])[0])


# In[7]:


sd_mbar = {}

for section in paths:
    sd_mbar[section]=float(getEnergyMBAR(paths[section])[1])


# In[8]:


energies_ti = {}

for section in paths:
    energies_ti[section]=float(getEnergyTI(paths[section]))


# In[9]:


lj_corr = float(getLJcorr('free/vanish/output/freenrg-LJCOR.dat')[0])
sd_lj_corr = float(getLJcorr('free/vanish/output/freenrg-LJCOR.dat')[1])


# In[10]:


dg_hyd_mbar=energies_mbar['vacuum_discharge']+energies_mbar['vacuum_vanish']-energies_mbar['free_discharge']-energies_mbar['free_vanish']-lj_corr


# In[11]:


dg_hyd_ti= energies_ti['vacuum_discharge']+energies_ti['vacuum_vanish']-energies_ti['free_discharge']-energies_ti['free_vanish']-lj_corr


# In[12]:


sd_dg_hyd_mbar=np.sqrt(sd_mbar['vacuum_discharge']**2+sd_mbar['vacuum_vanish']**2+sd_mbar['free_discharge']**2+sd_mbar['free_vanish']**2+sd_lj_corr**2)


# ### Write Output

# In[16]:


for section in paths:
    if section != 'free_vanish':
        print('Section: ',section)
        print('MBAR estimate:', energies_mbar[section],'+/-',sd_mbar[section],'kcal/mol')
        print('TI estimate:', energies_ti[section],'kcal/mol','\n')
    else:
        print('Section: ',section)
        print('MBAR estimate:', energies_mbar[section],'+/-',sd_mbar[section],'kcal/mol')
        print('TI estimate:', energies_ti[section],'kcal/mol')
        print('LJ correction: ', lj_corr,'+/-',sd_lj_corr,'kcal/mol','\n')
    
print('###########################################################################')

print('dg_hyd from MBAR: ',dg_hyd_mbar,' +/-', sd_dg_hyd_mbar,'kcal/mol')
print('dg_hyd from TI: ',dg_hyd_ti,'kcal/mol')

print('###########################################################################','\n')


# In[ ]:




