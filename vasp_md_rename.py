import glob
import os
OSZICARs = glob.glob('OSZICAR_*')
OUTCARs = glob.glob('OUTCAR_*')
XDATCARs = glob.glob('XDATCAR_*')
if len(OSZICARs)>0:
    os.system('cp CONTCAR POSCAR')
    OSZICAR_number = max([int(s.replace('OSZICAR_','')) for s in OSZICARs])
    os.system('cp OSZICAR OSZICAR_'+str(OSZICAR_number+1))
    os.system('cp OUTCAR OUTCAR_'+str(OSZICAR_number+1))
    os.system('cp XDATCAR XDATCAR_'+str(OSZICAR_number+1))
else:
    os.system('cp CONTCAR POSCAR')
    os.system('cp OSZICAR OSZICAR_1')
    os.system('cp OUTCAR OUTCAR_1')
    os.system('cp XDATCAR XDATCAR_1')