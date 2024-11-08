import os
import numpy as np
from pwtools import parse
from matplotlib import pyplot as plt

volDirs = [ f.name for f in os.scandir() if f.is_dir() ]
#volDirs.remove('.idea')
volDirs.remove('template')

data={'total energy': [], 'unitcell volume':[], 'A': []}

for volDir in volDirs:
    os.chdir(volDir)
    out=parse.PwSCFOutputFile('qe-pw.out')
    data['total energy'].append(out.get_etot())
    # parse volumes
    v=1
    print(out.get_cell())
    for i in range(3):
        v = v*out.get_cell()[i][i]
        print(out.get_cell()[i][i])
        if i == 0 :
            data['A'].append(out.get_cell()[i][i])

    data['unitcell volume'].append(v*0.14818471118724036) # *Bohr^3 Umrechnung in Angstrom^3
    os.chdir('..')


print(data)
print(os.getcwd())

plt.scatter(data['unitcell volume'],data['total energy'] )
plt.xlabel('unit cell volume [Angstrom$^3$]', fontsize=18)
plt.ylabel('$E_{tot}$ [Ry]',fontsize=18)
plt.axis([200, 220, -156.475, -156.473])
plt.subplots_adjust(left=0.2)
#plt.figure(figsize=[6.4, 4.8])
plt.savefig('Etot-V-plot.png')
'''

plt.scatter(data['unitcell volume'],data['A'] )
plt.xlabel('unit cell volume [Angstrom$^3$]', fontsize=18)
plt.ylabel('A',fontsize=18)
plt.savefig('A-V-plot.png')
'''