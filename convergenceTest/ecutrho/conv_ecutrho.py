import os

#defining tested energy cutoffs
ecutrho_multiplyers=list(range(3,12))
ecut=100 #Ry
print(ecutrho_multiplyers)

#go to the template directory and read the job and run files
os.chdir('template')
with open('qe-pw.job', 'r') as j:
    jobtemplate=j.read()
with open('run.inp', 'r') as r:
    runtemplate=r.read()
os.chdir('..')
for ecutrho_multiplyer in ecutrho_multiplyers:
    os.mkdir(str(ecutrho_multiplyer))
    os.chdir(str(ecutrho_multiplyer))
    with open('qe-pw.job', 'w') as j:
        j.write(jobtemplate)
    runtemplate_ecutrho=runtemplate.replace('MecutM', str(ecut))
    runtemplate_ecutrho=runtemplate_ecutrho.replace('MecutrhoM', str(ecutrho_multiplyer*ecut))
    with open('run.inp', 'w') as r:
        r.write(runtemplate_ecutrho)
    os.system('sbatch qe-pw.job')
    os.chdir('..')
