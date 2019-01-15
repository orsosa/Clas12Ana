# Clas12 Analysis tool
## Compiling
- define the environment variables HIPODIR, LZ4DIR and CLAS12ANA to the root installation of the different packages (the last one to this folder)
- run scons on this folder
- summary:
```
git clone https://github.com/orsosa/Clas12Ana.git
setenv HIPODIR <path_to_hipo>
setenv LZ4DIR <path_to_lz4>
cd Clas12Ana
setenv CLAS12ANA `pwd`
scons
```

## Two examples included:

- get_simple_tuple uses the library to read hipo files and produce root ntuples.
- particle_mix uses the ntuple produced by the get_simple_tuple and make combination of particles (reconstruction of short living particles, pions correlation, etc.)
