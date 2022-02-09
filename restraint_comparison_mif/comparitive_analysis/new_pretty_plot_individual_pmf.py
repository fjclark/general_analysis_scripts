#!/usr/bin/env python
# coding: utf-8
import sys
import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as colors

# Modified new_pretty_plot.py to produce one figure showing convergence plots
# for all stages.

# Limitations - assumes 1 epoch

nepochs = 1
       
#####################convergence-epoch-{nepochs}-{leg}-{stage}.dat

def plot_stage(datafile, ax, leg, stage):
    istream = open(datafile,'r')
    buffer = istream.readlines()
    istream.close()
    nruns = len(buffer[1].split())-3
    x_vals = []
    y_vals = []
    for y in range(0,1):
        y_vals.append([])
    y_avg = []
    error_vals = []
    
    epoch_step = 9
    epoch_counter = 1
    bars = []
    for line in buffer:
        if line.startswith("#"):
            continue
        elems = line.split()
        ctime = float(elems[0][:-1])
        x_vals.append(ctime)
        epoch_counter += 1
        if (epoch_counter > epoch_step):
            bars.append(ctime)
            #print ctime
            epoch_counter = 0
        for runid in range(0,1):
            y = float(elems[runid+1])
            y_vals[runid].append(y)
        DG_avg = float(elems[-2])
        CI = float(elems[-1])
        y_avg.append(DG_avg)
        error_vals.append(CI)
    # Convergence plot
    ax.plot(x_vals,y_avg)
    for entry in y_vals:
        ax.plot(x_vals,entry, linestyle='dashed')
    ax.set_title(f'{leg} {stage}')
    ax.set_ylabel('$\Delta \it{G}$ / kcal.mol$^-$$^1$')
    ax.set_xlabel('Cumulative sampling time / ns')
    ax.fill_between(x_vals, np.array(y_avg)-np.array(error_vals), np.array(y_avg)+np.array(error_vals), alpha=0.5, facecolor='#ffa500' )

def plotPMF(datafile, ax, leg, stage):
    '''Plots PMFs for each cumulative sampling time and 
    returns mapper which can be used for colorbar'''
    ofile = np.loadtxt(datafile, skiprows=1) 
    x_vals = ofile[:, 0]
    # create extra dimension for cumulative sampling time
    c = []
    c_incr = 0
    for y in range(len(ofile[0,1:])):
        c_incr += 0.05
        c.append(c_incr)
    # create mapper to map sampling time to colour
    norm = colors.Normalize(vmin=c[0], vmax=c[-1], clip=True)
    mapper = cm.ScalarMappable(norm=norm, cmap=cm.brg)
    # plot PMF for each time
    for y in range(len(ofile[0,1:])):
        ax.plot(x_vals, ofile[:, y+1], c=mapper.to_rgba(c[y]), alpha=0.2, lw=0.5)
    ax.set_title(f'{leg} {stage}')
    ax.set_xlabel(r'$\lambda$')
    ax.set_ylabel('$\Delta \it{G}$ / kcal.mol$^-$$^1$')

    return(mapper)

fig1 = plt.figure(figsize=(12,12))
#plot bound on top row
for i, stage in enumerate(['discharge','vanish']):
    ax = fig1.add_subplot(2,2,i+1)
    plot_stage(f'convergence-epoch-{nepochs}-bound-{stage}.dat', ax, 'bound', stage)
# plot free on bottom row
for i, stage in enumerate(['discharge','vanish']):
    ax = fig1.add_subplot(2,2,i+3)
    plot_stage(f'convergence-epoch-{nepochs}-free-{stage}.dat', ax, 'free', stage)
    
fig1.tight_layout()
fig1.savefig(f'convergence-epoch-{nepochs}_individual.png')
fig1.show()

fig2 = plt.figure(figsize=(12,12))
mapper = None
#plot bound on top row and get mapper
for i, stage in enumerate(['discharge','vanish']):
    ax = fig2.add_subplot(2,2,i+1)
    mapper = plotPMF(f'convergence-epoch-{nepochs}-bound-{stage}_pmf.dat', ax, 'bound', stage)
# plot free on bottom row
for i, stage in enumerate(['discharge','vanish']):
    ax = fig2.add_subplot(2,2,i+3)
    plotPMF(f'convergence-epoch-{nepochs}-free-{stage}_pmf.dat', ax, 'free', stage)
fig2.colorbar(mapper).set_label('Cumulative Sampling Time / ns')
fig2.tight_layout()
fig2.savefig(f'convergence-epoch-{nepochs}_individual_pmf.png')
fig2.show()
