import glob
import os

tests = [ f.name for f in os.scandir() if f.is_dir() ]
tests.remove('.idea')
tests.remove('template')
for test in tests:
    os.chdir(test)
    os.chdir('pwscf.save')
    file = glob.glob('./*')
    for f in file:
        os.remove(f)
    os.chdir('..')
    os.rmdir('pwscf.save')
    file=glob.glob('./*')
    for f in file:
        os.remove(f)
    os.chdir('..')
    os.rmdir(test)