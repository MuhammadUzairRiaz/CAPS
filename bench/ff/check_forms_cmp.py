import re
c=open('caps.txt').read().split('\n'); ce=[float(v) for v in c[0].split()]
cf={int(l.split()[0]):[float(v) for v in l.split()[1:]] for l in c[1:] if l.strip()}
log=open('log.lammps').read(); m=re.search(r'E_bond\s+E_angle.*\n\s*(.*)\n',log); le=[float(v) for v in m.group(1).split()]
worst=max(abs(a-b) for l in open('lmp.dump').read().split('ITEM: ATOMS')[1].split('\n')[1:] if l.strip() for a,b in zip(cf[int(l.split()[0])],[float(v) for v in l.split()[1:]]))
print(' '.join(f'{n}:{a:.6f}/{b:.6f}' for n,a,b in zip(['bond','angle','dih','imp','vdw'],ce,le) if abs(a)+abs(b)>0), f'maxdF {worst:.2e}', 'PASS' if worst<1e-4 and all(abs(a-b)<=1e-5*max(1,abs(b)) for a,b in zip(ce,le)) else 'FAIL')
