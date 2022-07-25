%oldchk=cas.chk
%chk=new
#P CAS(6,5, StateAverage, NRoot=3, FullDiag) guess=(read) NoSymm ChkBasis iop(4/21=10)

#DISP#

0 1
#XYZ#

0.333 0.333 0.333

