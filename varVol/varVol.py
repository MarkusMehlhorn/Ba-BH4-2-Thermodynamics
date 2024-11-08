import os
import numpy as np

#defining the volume variation range
VolRange = np.array([1.055, 1.065, 1.075, 1.085, 1.095, 1.105, 1.115, 1.125, 1.135, 1.145])
print(VolRange)

# select the volumes

runVols = VolRange

#go to the template directory and read the job file and the run files
os.chdir('template')

with open('qe-pw.job', 'r') as j:
    jobtemplate=j.read()
with open('run.inp', 'r') as r:
    runtemplate=r.readlines()

#search the run template for the atomic coordinates to remove them later
# and read the number of atoms
for i in range(len(runtemplate)):
    if 'nat =' in runtemplate[i]:
        nat = int(runtemplate[i].split()[2])
    elif 'ATOMIC_POSITIONS' in runtemplate[i]:
        positionLine = i
with open('qe-pw.out', 'r') as o:
    output=o.readlines()

# read the cell constants of the equilibrium structure
positions=[]
for i in range(len(output)):
    if 'End final coordinates' in output[i]:
        coordinatesPosition=i
    elif 'lattice parameter (alat)' in output[i]:
        alat=float(output[i].split()[4])

A = float(output[coordinatesPosition-nat-5].split()[0])*alat*0.5291772105638411
B = float(output[coordinatesPosition-nat-4].split()[1])*alat*0.5291772105638411
C = float(output[coordinatesPosition-nat-3].split()[2])*alat*0.5291772105638411 #(Bohr) (no ase on dyson )


# write the runs
prerun=[]
for i in range(len(runtemplate)):
    if i not in range(positionLine, positionLine+nat+1):
        if 'A =' in runtemplate[i]:
            prerun.append(' A = '+str(A))
        elif 'B =' in runtemplate[i]:
            prerun.append(' B = '+str(B))
        elif 'C =' in runtemplate[i]:
            prerun.append(' C = '+str(C))
        else:
            prerun.append(runtemplate[i])
# write the equilibrium atomic coordinates to the prerun
for i in range(coordinatesPosition-nat-1, coordinatesPosition):
    prerun.append(output[i])

# print the prerun
for i in range(len(prerun)):
    print(prerun[i])
os.chdir('..')


for vol in runVols:
    # create a directory for every volume
    os.mkdir(str(round(vol, 3)))
    os.chdir(str(round(vol, 3)))
    # change the name of the job
    prejob=jobtemplate.replace("varVol", str(round(vol, 3)))
    # write job script in the directory
    with open('qe-pw.job', 'w') as j:
        j.write(prejob)
    run = []
    # adjust the cell constants for the certain volume
    for i in range(len(prerun)):
        if 'A =' in prerun[i]:
            run.append('  A = '+str(A*(vol**(1/3)))+',\n')
        elif 'B =' in prerun[i]:
            run.append('  B = '+str(B*(vol**(1/3)))+',\n')
        elif 'C =' in prerun[i]:
            run.append('  C = '+str(C*(vol**(1/3)))+',\n')
        else:
            run.append(prerun[i])
    # write the run file for the certain volume
    with open('run.inp', 'w') as r:
        r.writelines(run)
    # start the job
    os.system('sbatch qe-pw.job')
    os.chdir('..')