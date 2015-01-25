module m1

using JFinEALE

t0 = time()

rho=1.21*1e-9;# mass density
c =345.0*1000;# millimeters per second
bulk= c^2*rho;
Lx=1900.0;# length of the box, millimeters
Ly=800.0; # length of the box, millimeters

fens,fes = Q4block(Lx,Ly,3,2); # Mesh
show(fes.conn)

end
