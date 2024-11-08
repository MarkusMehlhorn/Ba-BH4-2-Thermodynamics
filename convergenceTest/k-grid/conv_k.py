import os

#defining tested energy cutoffs
ecutrho_multiplyer=4
ecut=100 #Ry

Ns=list(range(1,5))

#go to the template directory and read the job and run files
os.chdir('template')
with open('qe-pw.job', 'r') as j:
    jobtemplate=j.read()
with open('run.inp', 'r') as r:
    runtemplate=r.read()
os.chdir('..')
for N in Ns:
    os.mkdir(str(N))
    os.chdir(str(N))
    with open('qe-pw.job', 'w') as j:
        j.write(jobtemplate)
    runtemplate_ecut=runtemplate.replace('MecutM', str(ecut))
    runtemplate_ecutrho=runtemplate_ecut.replace('MecutrhoM', str(ecutrho_multiplyer*ecut))
    runtemplate_N=runtemplate_ecutrho.replace('MNM', str(N))
    with open('run.inp', 'w') as r:
        r.write(runtemplate_N)
    os.system('sbatch qe-pw.job')
    os.chdir('..')
