# -*- python -*-
from subprocess import PIPE,Popen
import os
import sys
import glob

env = Environment();
ENV = os.environ
env.Append(ENV = ENV)

if "HIPODIR" in ENV.keys():
 HIPODIR =   ENV["HIPODIR"]
else:
 print "you must define hipo root installation directory"


if "LZ4DIR" in ENV.keys():
 LZ4DIR =   ENV["LZ4DIR"]
else:
 print "not using lz4, define LZ4DIR variable if you have it."
 Exit(1)

if "CLAS12ANA" in ENV.keys():
 CLAS12ANA =   ENV["CLAS12ANA"]
else:
 print "you must define CLAS12ANA root directory, may be: export CLAS12ANA=`pwd`"
 Exit(1)	    

####### ROOT ENVIRONMENT #####
out = Popen(["root-config","--glibs"],stdout=PIPE)
glibs_out = out.communicate()[0].replace("\n","").split()
glibs = map(lambda x: x.lstrip('-l'),filter(lambda x: '-l' in x, glibs_out))
libdir =  map(lambda x: x.lstrip("-L"),filter(lambda x: '-L' in x, glibs_out))
ldflags=  filter(lambda x: '-l' not in x and '-L' not in x, glibs_out)

out = Popen(["root-config","--cflags"],stdout=PIPE)
cflags_out = out.communicate()[0].replace("\n","").split()
cflags = filter(lambda x: "-I" not in x,cflags_out)
incdir = map(lambda x: x.lstrip("-I"),filter(lambda x: "-I" in x,cflags_out))

######################
LZ4DIR=LZ4DIR.rstrip("/")
HIPODIR=HIPODIR.rstrip("/")

SRCDIR ="src"
INCDIR ="include"
SHLIBDIR = CLAS12ANA + "/shlib"

print SHLIBDIR
env.Append(CPPPATH=[INCDIR,SRCDIR,"/usr/include","/usr/local/include","/opt/local/include",LZ4DIR + "/lib",HIPODIR + "/libcpp"])
env.Append(CPPPATH=incdir)

env.Append(CCFLAGS=["-O2","-fPIC","-m64","-fmessage-length=0","-g"])
env.Append(CCFLAGS=cflags)

env.Append(LINKFLAGS=ldflags)
env.Append(LIBPATH=["/opt/local/lib","/usr/lib","/usr/local/lib",LZ4DIR + "/lib","lib",HIPODIR + "/lib"])
env.Append(LIBPATH=libdir)

env.Append(CONFIGUREDIR=[LZ4DIR + "/lib",HIPODIR + "/lib"])

env.Append(LIBS=["Gui" ,"Core" ,"Imt" ,"RIO" ,"Net", "Hist", "Graf" ,"Graf3d" ,"Gpad" ,"Tree" ,"TreePlayer" ,"Rint" ,"Postscript" ,"Matrix" ,"Physics" ,"MathCore","Thread", "MultiProc" ,"m", "dl","EG"])

lib_target = "TIdentificatorCLAS12"
libsrc = ["Categorize.cxx","TIdentificatorCLAS12.cxx","node_assignment.cxx"]

libsrc = map(lambda x: SRCDIR + "/" + x,libsrc)
lib_obj = map(lambda x: x.replace(".cxx",".o"),libsrc)

conf = Configure(env)
if conf.CheckLib('libhipo'):
   print '\n\033[32m[**] >>>>> found library : HIPO'
   print ''
   env.Append(CCFLAGS="-D__HIPO__")
    
if conf.CheckLib('liblz4'):
   print '\n\033[32m[**] >>>>> found library : LZ4'
   print '[**] >>>>> enabling lz4 compression. \033[0m'
   print ''
   env.Append(CCFLAGS="-D__LZ4__")

if conf.CheckLib('libz'):
   print '\n\033[32m[**] >>>>> found library : libz'
   print '[**] >>>>> enabling gzip compression. \033[0m'
   print ''
   env.Append(CCFLAGS="-D__LIBZ__")

#object = env.Object(target="TIdentificatorCLAS12.o",source = libsrc)
shlib = env.SharedLibrary(target =  "TIdentificatorCLAS12", source = libsrc)
env.Command(SHLIBDIR + "/${SOURCE.file}" , shlib, Copy("$TARGET", "$SOURCE"))

Export('env shlib SHLIBDIR')
SConscript("examples/SConscript")
