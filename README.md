The program gibbs is written in Fortran-77. 
It calculates the equilibrium volume and density of two phases of a liquid
as desribed by the Panagiotopiulos paper.
The Mone Carlo sampling is enhanced by cavity biasing for particle exchange and 
virial biasing for volume exchange.

The array dimensions are represented by special symbols in the file gibbs.for and 
creates the file gibbs.f.
Edit the script as follows to thange the default values in pre.csh.
# maxmol: maximum number of molecules
set maxmol = 1000
# maxrgrd (#RG) : Maximum no of grid points for g(r)'s 
set maxrgrd = 210
# maxslv (#SV) : Maximum number of atoms per molecule
set maxslv = 10
# #GX : Maximum no of grid points for cavity grids
set maxcgrd = 100 
# maxcav (#CV)  : Maximum no of cavities
set maxcav = 1000000
Ececute the sript as 
pre.csh gibbs
and compile the filew gibbs.f with a Fortran compiler.
