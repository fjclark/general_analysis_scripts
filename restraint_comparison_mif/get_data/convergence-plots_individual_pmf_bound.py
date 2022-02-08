#!/usr/bin/env python
# coding: utf-8

# Modified from convergence-plots_debugged.py to produce plots for each stage
# and to include MBAR PMF outputs.

# 1) Work out binding free energies using all available sampling in the subfolders
# 2) Work out list of windows whose stdev exceeds threshold 
# 3) Generate submission scripts to carry out one round of adaptive sampling
import glob,sys,os
import copy
import numpy as np
import math
import matplotlib.pyplot as plt
import pickle

#debug
import pdb

### ADAPTIVE SAMPLING ROUTE
NRUNS = 1
EPOCHTIME = 11 # ns  
DISCARD = 6 # ns
CHUNKSIZE = 0.05 # ns
TIMESTEP = 0.000004 #ns
NRGFREQ = 100

RUNFORPMF=1
### HELPER FUNCTIONS
def doMBAR(simfiles, startime, endtime):
    cmd = "rm -rf tmp ; mkdir tmp" 
    os.system(cmd)
    cumtime = 0.0
    for simfile in simfiles:
        #print (simfile, max_time)
        # Ugly
        lamval = os.path.split(simfile)[-2].split("/")[-1]
        #print (lamval)
        istream = open(simfile,'r')
        ostream = open("tmp/%s.dat" % lamval,'w')
        ifile = istream.readlines()
        for line in ifile:
            if line.startswith("#"):
                ostream.write(line)
                continue
            elems = line.split()
            time = (float(elems[0])+NRGFREQ)*TIMESTEP
            if (time < startime):
                continue
            if (time > endtime):
                break
            ostream.write(line)
        cumtime += time#(time - startime)
        istream.close()
        ostream.close()
    # TODO ADD OVERLAP MATRIX TO OUTPUT
    cmd = "/home/finlayclark/anaconda3/envs/biosimspace-dev/bin/analyse_freenrg mbar -i tmp/*.dat -o tmp/mbar.dat -p 100"
    os.system(cmd)
    istream = open("tmp/mbar.dat","r")
    results = {}
    ofile = istream.readlines()
    pmf_mbar = np.empty([3])
    no_lam_vals = len(simfiles)
    inDGsection = False
    for l, line in enumerate(ofile):
        if line.startswith("#DG from neighbouring"):
            inDGsection = True
            continue
        elif line.startswith("#PMF from MBAR"):
            inDGsection = False
            for i in range(l+1, l+no_lam_vals+1):
                pmf_mbar = np.vstack([pmf_mbar, [float(i) for i in ofile[i].strip().split(' ')]])
        elif inDGsection:
            elems = line.split()
            lami = float(elems[0])
            lamj = float(elems[1])
            DGij = float(elems[2])
            DGij_sig = float(elems[3])
            results["%.5f-%.5f" % (lami,lamj)] = (DGij,DGij_sig)
    totline = ofile[-3].split()
    mbar = float(totline[0].strip(","))
    mbar_err = float(totline[1])
    results['DGtot'] = (mbar, mbar_err, cumtime, pmf_mbar[1:, ])
        
    return results


def calcEnergies(energies, basefolder, simprofile, startime, endtime):
    for run in range(1,NRUNS+1):
        leg = "bound"
        for stage in ("discharge","vanish"):
            rundir = "%s/%s/run00%d/%s/output/lambda-*/simfile.dat" % (basefolder,leg,run,stage)
            simfiles = glob.glob(rundir)
            # use analyse_freenrg to get free energy of chunks
            # time left after subtracting DISCARD should be multiple of CHUNKSIZE
            nchunks = int((endtime-DISCARD)/CHUNKSIZE)
            for chunk in range(nchunks):
                try:
                    energies[run]
                except KeyError:
                    energies[run] = {}
                try:
                    energies[run][leg]
                except KeyError:
                    energies[run][leg] = {}
                try:
                    energies[run][leg][stage] 
                except KeyError:
                    energies[run][leg][stage] = {}
                startime = DISCARD
                chunkendtime = ((chunk+1)*CHUNKSIZE)+DISCARD
                chunkendtime = round(chunkendtime,2)
                try:
                    DG = energies[run][leg][stage][(startime,chunkendtime)]
                    print ("USING cached MBAR results for %s, start %s, end %s, DG %s " % (rundir,startime,chunkendtime,DG))
                except KeyError:
                    # Check if all windows have a maxendtime < endtime, if so use previous endtime result
                    maxsimfiletime = -1
                    for simfile in simfiles:
                        simfiletime = simprofile[run][leg][stage][simfile]
                        simfiletime = round(simfiletime,2)
                        if maxsimfiletime < simfiletime:
                            maxsimfiletime = simfiletime
                    if maxsimfiletime < chunkendtime:
                        DG = energies[run][leg][stage][(startime,maxsimfiletime)]
                        print ("No additional data to process over this time interval (%s-%s)" % (startime, chunkendtime))
                        energies[run][leg][stage][(startime,chunkendtime)] = DG
                        print ("reusing cached MBAR results %s, start %s, end %s, DG %s " % (rundir,startime,chunkendtime,DG))
                    else:
                        # Otherwise do MBAR from scratch
                        DG = doMBAR(simfiles, startime, chunkendtime)
                        print ("DOING NEW MBAR on %s, start %s, end %s, DG %s " % (rundir,startime,chunkendtime,DG))
                        energies[run][leg][stage][(startime,chunkendtime)] = DG


def plotDGbind(energies, nepochs=0):
    runs = list(energies.keys())
    runs.sort()
    chunks = list(energies[runs[0]]['bound']['discharge'].keys())
    chunks.sort()
    x_vals = []
    y_vals = []
    error_vals = []
    y_runs = []
    for x in range(0,NRUNS):
        y_runs.append([])
    for chunk in chunks:
        DGbinds = []
        cumtime = 0.0  
        for run in runs:
            DG_bound_discharge = energies[run]['bound']['discharge'][chunk]['DGtot'][0]
            DG_bound_discharge_time = energies[run]['bound']['discharge'][chunk]['DGtot'][2]
            DG_bound_vanish = energies[run]['bound']['vanish'][chunk]['DGtot'][0]
            DG_bound_vanish_time = energies[run]['bound']['vanish'][chunk]['DGtot'][2]
            DGbind = - (DG_bound_discharge + DG_bound_vanish)
            DGbinds.append(DGbind)
            chunktime = DG_bound_discharge_time + DG_bound_vanish_time
            cumtime += chunktime
            y_runs[run-1].append(DGbind)
            #print (chunk,DGbind)
        DGbind_avg = np.array(DGbinds).mean()
        DGbind_ste = np.array(DGbinds).std()/math.sqrt(len(DGbinds))*1.96# 95% CI
        print (chunk[1], DGbind_avg,DGbind_ste)
        x_vals.append(cumtime)
        y_vals.append(DGbind_avg)
        error_vals.append(DGbind_ste)
    
    ostream = open("convergence-epoch-%s.dat" % nepochs,"w")
    header = "# cumulative_time DGrunX... avgDG 95CI\n"
    ostream.write(header)
    for x in range(0,len(x_vals)):
        line = " %8.5f " % x_vals[x]
        for entry in y_runs:
            line += " %8.5f " % entry[x]
        line += " %8.5f %8.5f \n" % (y_vals[x],error_vals[x])
        ostream.write(line)
    ostream.close()

    #print (x_vals)
    # Convergence plot
# =============================================================================
#     plt.plot(x_vals,y_vals)
#     for entry in y_runs:
#         plt.plot(x_vals,entry)
#     plt.ylabel('DG / kcal.mol-1')
#     plt.xlabel('Cumulative sampling time / ns')
#     plt.fill_between(x_vals, np.array(y_vals)-np.array(error_vals), np.array(y_vals)+np.array(error_vals), alpha=0.5, face#color='#ffa500' )
#     plt.show()
# =============================================================================


# New function to plot convergence for each individual stage 
def plotDGbindIndividual(energies, leg, stage, nepochs=0):
    runs = list(energies.keys())
    runs.sort()
    chunks = list(energies[runs[0]]['bound']['discharge'].keys())
    chunks.sort()
    x_vals = []
    y_vals = []
    error_vals = []
    y_runs = []
    for x in range(0,NRUNS):
        y_runs.append([])
    for chunk in chunks:
        DGs = []
        cumtime = 0.0  
        for run in runs:
            DG = energies[run][leg][stage][chunk]['DGtot'][0]
            DGs.append(DG)
            cumtime += energies[run][leg][stage][chunk]['DGtot'][2]
            y_runs[run-1].append(DG)
            #print (chunk,DGbind)
        DG_avg = np.array(DGs).mean()
        DG_ste = np.array(DGs).std()/math.sqrt(len(DGs))*1.96# 95% CI
        print (chunk[1], DG_avg,DG_ste)
        x_vals.append(cumtime)
        y_vals.append(DG_avg)
        error_vals.append(DG_ste)
    
    ostream = open(f"convergence-epoch-{nepochs}-{leg}-{stage}.dat", "w")
    header = "# cumulative_time DGrunX... avgDG 95CI\n"
    ostream.write(header)
    for x in range(0,len(x_vals)):
        line = " %8.5f " % x_vals[x]
        for entry in y_runs:
            line += " %8.5f " % entry[x]
        line += " %8.5f %8.5f \n" % (y_vals[x],error_vals[x])
        ostream.write(line)
    ostream.close()
    

def plotPMFbyChunk(run, energies, leg, stage, nepochs=0):
    '''
    Plot PMFs for each stage at each chunk for one run

    Parameters
    ----------
    run : str
    energies : dict
    leg : str
    stage : str
    nepochs : int, optional
        DESCRIPTION. The default is 0.

    Returns
    -------
    None.

    '''
    chunks = list(energies[run]['bound']['discharge'].keys())
    chunks.sort()
    x_vals = energies[run][leg][stage][chunks[0]]['DGtot'][3][:,0] #gets lam vals
    y_vals = []
    for chunk in chunks:
        y_vals.append(energies[run][leg][stage][chunk]['DGtot'][3][:,1])
    
    ostream = open(f"convergence-epoch-{nepochs}-{leg}-{stage}_pmf.dat", "w")
    header = "# PMF data. First column is lambda values, following columns give PMF for each value of cumulative sampling\n"
    ostream.write(header)
    for x in range(len(x_vals)):
        line = " %8.5f " % x_vals[x]
        for y_array in y_vals:
            line += " %8.5f " % y_array[x]
        line += "\n"
        ostream.write(line)
    ostream.close()


#####################
basefolder = "./"
simprofile={}
# scan data to work out in which epoch we are
maxendtime = -1
for run in range(1,NRUNS+1):
    simprofile[run] = {}
    leg = "bound"
    simprofile[run][leg] = {}
    for stage in ("discharge","vanish"):
        simprofile[run][leg][stage] = {}
        rundir = "%s/%s/run00%d/%s/output/lambda-*/simfile.dat" % (basefolder,leg,run,stage)
        simfiles = glob.glob(rundir)
        for simfile in simfiles:
            #Find how much sampling was done
            istream = open(simfile,'r')
            buffer = istream.readlines()
            istream.close()
            nsteps = int(buffer[-1].split()[0])
            endtime = (nsteps+NRGFREQ)*TIMESTEP
            if endtime > maxendtime:
                maxendtime = endtime
            simprofile[run][leg][stage][simfile] = endtime
#print (simprofile)
print ("MAXIMUM SAMPLING TIME IS %s ns " % maxendtime)
epoch = int(maxendtime/EPOCHTIME)
print ("WE ARE IN EPOCH %s " % epoch)
# Load cached energies (if any)
if os.path.exists("energies.pickle"):
    energies = pickle.load(open("energies.pickle","rb"))
else:
    energies = {}
#sys.exit(-1)
# Now compute energies for each sim using variable chunks
calcEnergies(energies, basefolder, simprofile, DISCARD, maxendtime)
plotDGbind(energies, nepochs=epoch)
# Now compute energies for each individual stage
for leg in ['bound']:
    for stage in ['discharge','vanish']:
        plotDGbindIndividual(energies, leg, stage, nepochs=epoch)
# Now calculate PMFs for each stage
for leg in ['bound']:
    for stage in ['discharge','vanish']:
        plotPMFbyChunk(RUNFORPMF, energies, leg, stage, nepochs=epoch)
# Save energies processed 
dumpme = open("energies.pickle","wb")
pickle.dump(energies, dumpme)
dumpme.close()

