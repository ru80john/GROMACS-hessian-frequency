%nprocshared=12
%mem=50000MB
%chk=BEN_hessian.chk
#Opt Freq wb97xD  6-311++G(d, p) pop=(CM5, ESP, NBOREAD) 

BEN

0 1
  C     0.000000     0.695000     1.205000
  C     0.000000    -0.695000     1.205000
  C     0.000000     1.391000     0.000000
  H     0.000000     1.238000     2.144000
  C     0.000000    -1.391000     0.000000
  H     0.000000    -1.238000     2.144000
  C     0.000000    -0.695000    -1.205000
  H     0.000000    -2.475000     0.000000
  C     0.000000     0.695000    -1.205000
  H     0.000000    -1.238000    -2.144000
  H     0.000000     1.238000    -2.144000
  H     0.000000     2.475000     0.000000

$nbo BNDIDX $end

