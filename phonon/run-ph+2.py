import os
import numpy as np
import shutil as sh


runVols=np.array([0.955, 0.965, 0.975, 0.985, 0.995, 1.005, 1.015, 1.025, 1.035, 1.045])

data = {}

# go to the template directory
os.chdir('../phonon_tr2/template')
with open('qe-ph+2.job','r') as t:
    jobTemplate=t.read()
os.chdir('..')

for vol in runVols:
    # create volume directories and put run files for q2r and matdyn into it
    sh.copy('./template/run-matdyn2.inp', './' + str(round(vol, 3)))
    os.chdir(str(round(vol, 3)))
    # name the jobs and write the job file for the complete run
    job = jobTemplate.replace("phdos", "dos"+str(round(vol, 3)))
    with open('qe-ph+2.job','w') as t:
        t.write(job)
    # write a runfile with the respective geometry
    prerun = []
    # start the job
    print("calculating volume "+str(round(vol, 3)))
    os.system('sbatch qe-ph+2.job')
    os.chdir('..')
