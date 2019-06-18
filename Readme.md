# Clas12 Analysis tool for Hipo4 format
## Requisites
- hipo library
  https://github.com/gavalian/Clas12Tool/tree/hipo4
```
git clone -b hipo4 --recurse-submodules https://github.com/gavalian/Clas12Tool.git
cd Clas12Tool
mkdir lib # could be included in Makefile...
make
setenv CLAS12TOOL `pwd` # or you can add it to your login script
```
  
## Compiling
- define the environment variables HIPO4DIR, LZ4DIR and CLAS12ANAHIPO4 to the root installation of the different packages (the last one to this folder)
- run make on this folder
- summary:
```
git clone -b hipo4 https://github.com/orsosa/Clas12Ana.git
setenv HIPO4DIR <path_to_hipo>
setenv LZ4DIR <path_to_lz4>
cd Clas12Ana
setenv CLAS12ANAHIPO4 `pwd`
make
```
- don't forget to include the CLAS12ANA/shlib to $LD_LIBRARY_PATH or $DYLD_LIBRARY_PATH

## Two examples included:

- get_simple_tuple uses the library to read hipo files and produce root ntuples.
- particle_mix uses the ntuple produced by the get_simple_tuple and make combination of particles (reconstruction of short living particles, pions correlation, etc.)


- examples not running because shared library libTIdentificatorClas12.so not found. Try adding the ${CLAS12ANAHIPO4}/shlib at the begining of your LD_LIBRARY_PATH
