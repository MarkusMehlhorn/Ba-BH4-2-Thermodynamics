import os
import numpy as np
import shutil as sh


# create volume points
# must be the same as in the varVol directory
VolRange = np.linspace(0.95,1.05,21)
#print(VolRange)

# split volumes in two groups
run1Vols = np.empty(11)
for i in range(len(VolRange)):
    if divmod(i, 2)[1] != 1:
        run1Vols[int(i/2)]=VolRange[i]
#print(run1Vols)
run2Vols = np.empty(10)
for i in range(len(VolRange)):
    if divmod(i, 2)[1] != 0:
        run2Vols[int(i/2)]=VolRange[i]
#print(run2Vols)

# select the volumes to be calcumlated
runVols = run2Vols
#print(runVols)

# go to the varVol directory
os.chdir('../varVol_tr2')

data = {}
for vol in runVols:
    # go to the volume directories and read the outputsfiles of the respective geometry optimization
    print(str(round(vol, 3))+":")
    os.chdir(str(round(vol, 3)))
    with open('qe-pw.out', 'r') as o:
        rawOut = o.readlines()
    # parse the output files
    for i in range(len(rawOut)):
        if "new unit-cell volume" in rawOut[i]:
            volume = float(rawOut[i].split()[-3])
        elif "celldm(1) =" in rawOut[i]:
            celldm1 = float(rawOut[i].split()[-1])*0.5291772105638411 # *Bohr Umrechnung in Angstrom
            celldm2 = float(rawOut[i+1].split()[-1])
            celldm3 = float(rawOut[i+2].split()[-1])
        elif 'number of atoms/cell' in rawOut[i]:
            nat = int(rawOut[i].split()[-1])
        elif "Begin final coordinates" in rawOut[i]:
            ifc=i

    geometry=[]
    for i in range(ifc+9,ifc+10+nat):
        geometry.append(rawOut[i])
    # transfer variables to the dictionary and calculate the cellconstants from the celldms
    data[vol]={}
    data[vol]['cellconstants']=[]
    data[vol]['cellconstants'].append(celldm1) # A
    data[vol]['cellconstants'].append(celldm2*celldm1) # B
    data[vol]['cellconstants'].append(celldm3*celldm1) # C
    data[vol]['cellvolume'] = volume
    data[vol]['geometry'] = geometry
    os.chdir('..')
    print('analysed')

# go to the template directory
os.chdir('../phonon_tr2/template')
with open('qe-ph+.job','r') as t:
    jobTemplate=t.read()
with open('run-pw.inp', 'r') as t:
    runTemplate=t.readlines()
os.chdir('..')

for vol in runVols:
    # create volume directories and put run files for ph, q2r and matdyn into it
    os.mkdir(str(round(vol, 3)))
    sh.copy('./template/run-ph.inp', './'+str(round(vol, 3)))
    sh.copy('./template/run-q2r.inp', './' + str(round(vol, 3)))
    sh.copy('./template/run-matdyn.inp', './' + str(round(vol, 3)))
    os.chdir(str(round(vol, 3)))
    # name the jobs and write the job file for the complete run
    job = jobTemplate.replace("phdos", "dos"+str(round(vol, 3)))
    with open('qe-ph+.job','w') as t:
        t.write(job)
    # write a runfile with the respective geometry
    prerun = []
    for i in range(len(runTemplate)-nat-1):
            if 'A =' in runTemplate[i]:
                prerun.append(' A = ' + str(data[vol]['cellconstants'][0])+',\n')
            elif 'B =' in runTemplate[i]:
                prerun.append(' B = ' + str(data[vol]['cellconstants'][1])+',\n')
            elif 'C =' in runTemplate[i]:
                prerun.append(' C = ' + str(data[vol]['cellconstants'][2])+',\n')
            else:
                prerun.append(runTemplate[i])
    for i in range(len(data[vol]['geometry'])):
        prerun.append(data[vol]['geometry'][i])
    with open('run-pw.inp', 'w') as t:
        t.writelines(prerun)
    # start the job
    print("calculating volume "+str(round(vol, 3)))
    os.system('sbatch qe-ph+.job')
    os.chdir('..')
