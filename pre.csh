#/bin/csh
# This c-shell script dimensions the Gibbs ens. Monte Carlo program. It assumes
# an input file with .for extension and prepares a new one with .f extension.
echo  "Arguments: file [, max molecules, max RDF grids]"
set file = $argv[1]
set old = $argv[1].for
set new = $file.f
set bin = $file.bin
echo Old file: $old
echo New file: $new
echo New executable: $bin
# maxmol (#MO) : Maximum no of molecules
set maxmol = 1000
if ( $#argv  > 1 ) then
  set maxmol = $argv[2]
endif
# maxdist (#IJ) : Maximum no of entries in the distance matrix
@ maxdist = ( $maxmol * ( $maxmol - 1 ) ) / 2
# maxrgrd (#RG) : Maximum no of grid points for g(r)'s 
set maxrgrd = 210
if ( $#argv > 2 ) then
  set maxrgrd = $argv[3]
endif
# maxslv (#SV) : Maximum number of atoms per molecule
set maxslv = 10
# #NA: maximum number of atoms
@ maxat = $maxmol * $maxslv
# #GX : Maximum no of grid points for cavity grids
set maxcgrd = 100 
# maxcav (#CV)  : Maximum no of cavities
set maxcav = 1000000
echo Maximum number of molecules: $maxmol
echo Maximum number of RDF grids: $maxrgrd
echo Maximum number of atoms per molecule: $maxslv
echo Maximum number of cavities: $maxcav
echo Maximum number of cavity grids: $maxcgrd 
set sub1 = 's/#MO/'$maxmol'/g'
set sub2 = 's/#IJ/'$maxdist'/g'
set sub3 = 's/#RG/'$maxrgrd'/g'
set sub4 = 's/#NA/'$maxat'/g'
set sub5 = 's/#SV/'$maxslv'/g'
set sub6 = 's/#CV/'$maxcav'/g'
set sub7 = 's/#GX/'$maxcgrd'/g'
echo "Start substitution"
cat $old | sed $sub1 | sed $sub2 | sed $sub3 | sed $sub4 | sed \
               $sub5 | sed $sub6 | sed $sub7 > $new
echo "Ready to compile"
