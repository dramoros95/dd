# TEST 4 : optimized code (iwf=2, izv=0) vs the original program with the Z(r)
# term removed (qxx=0) and the 4p coefficient 0.03172, on a reduced p grid,
# 3 angles. Usage: python3 t_regression.py  (run in kr/tests)
# Expected: identical TDCS up to the 4p normalisation factor 1/norm^2 = 1.0000025
import subprocess
o=open('../original/DI-Kr-4p-BBK-zv-3P.for').read()
o=o.replace('0.003172','0.03172').replace("c\tprint*,qxx1,qxx2\n","\tqxx1=(0.d0,0.d0)\n\tqxx2=(0.d0,0.d0)\n",1)
for a,b in [('np1=20','np1=4'),('np2=120','np2=8'),('np3=40','np3=4'),('np4=5','np4=2'),('nt=80','nt=12'),('nh=80','nh=12'),('do 800 i1=1,n1+1','do 800 i1=1,3')]:
    o=o.replace(a,b,1)
open('kr_reg.for','w').write(o)
n=open('../kr3e_opt.for').read()
for a,b in [('      iwf=1\n','      iwf=2\n'),('      izv=1\n','      izv=0\n'),('np1=20 ; np2=120 ; np3=40 ; np4=5','np1=4 ; np2=8 ; np3=4 ; np4=2'),('      nt=80\n','      nt=12\n'),('      nh=80\n','      nh=12\n'),('      acon=2.d0/27.21165d0','      acon=2./27.21165')]:
    n=n.replace(a,b,1)
open('new_reg.for','w').write(n)
f='-O2 -fopenmp -ffixed-line-length-132 -w'.split()
subprocess.run(['gfortran','-O2','-ffixed-line-length-132','-w','kr_reg.for','-o','kr_reg'],check=True)  # no OpenMP: its 960 MB arrays must stay static
subprocess.run(['gfortran']+f+['new_reg.for','-o','new_reg'],check=True)
a=[l.split() for l in subprocess.run(['./kr_reg'],capture_output=True,text=True).stdout.splitlines() if 'attention' not in l]
subprocess.run(['./new_reg'],input='1 3\n',capture_output=True,text=True)
b=[l.split() for l in open('kr3e_3P.dat')]
for x,y in zip(a,b): print('theta1=%6.1f original %.10e  optimized %.10e  ratio %.8f'%(float(x[0]),float(x[1]),float(y[1]),float(y[1])/float(x[1])))
