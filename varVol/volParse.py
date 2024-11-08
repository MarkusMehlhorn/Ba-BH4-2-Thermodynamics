import os
import numpy as np
from pwtools import parse
from matplotlib import pyplot as plt

volDirs = [ f.name for f in os.scandir() if f.is_dir() ]
volDirs.remove('.idea')
volDirs.remove('template')

data={'total energy': [], 'unitcell volume': [], 'A' :[], 'B': [], 'C': []}
volDirs.sort()
print(volDirs)

for volDir in volDirs:
    os.chdir(volDir)
    out=parse.PwMDOutputFile('qe-pw.out')
    data['total energy'].append(out.get_etot()[-1])
    # parse volumes
    v=1

    print(volDir)
    for i in range(3):
        v = v*out.get_cell()[-1][i][i]
        if i == 0 :
            data['A'].append(out.get_cell()[-1][i][i])
        elif i == 1 :
            data['B'].append(out.get_cell()[-1][i][i])
        elif i == 2 :
            data['C'].append(out.get_cell()[-1][i][i])
    data['unitcell volume'].append(v*0.14818471118724036) # *Bohr^3 Umrechnung in Angstrom^3
    os.chdir('..')


intDir=np.array(volDirs)

print(data)
print(os.getcwd())

plt.scatter(data['unitcell volume'],data['total energy'] )
#plt.scatter(intDir ,data['total energy'] )
plt.xlabel('unit cell volume [Angstrom$^3$]',fontsize=18)
plt.ylabel('$E_{tot}$ [Ry]',fontsize=18)
plt.axis([200, 235, -156.476, -156.469])
plt.subplots_adjust(left=0.2)
#plt.figure(figsize=[6.4, 4.8])
plt.savefig('Etot-V-plot.png')
'''

plt.scatter(data['unitcell volume'],data['A'] )
plt.scatter(data['unitcell volume'],data['B'] )
plt.scatter(data['unitcell volume'],data['C'] )
plt.xlabel('unit cell volume [Angstrom$^3$]', fontsize=18)
plt.ylabel('A',fontsize=18)
plt.savefig('A-V-plot.png')
'''