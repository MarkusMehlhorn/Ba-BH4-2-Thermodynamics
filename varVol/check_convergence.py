import os
import numpy as np

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

for vol in runVols:
    # go to the volume directories and read the output files of the respective geometry optimization
    print(str(round(vol, 3))+":")
    os.chdir(str(round(vol, 3)))
    with open('qe-pw.out', 'r') as o:
        rawOut = o.readlines()
    # parse the output files
    convergence = False
    for i in range(len(rawOut)):
        if "Begin final coordinates" in rawOut[i]:
            convergence=True
    if convergence == True:
        print("converged")
    else:
        print("not converged")
    os.chdir('..')
