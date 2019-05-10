#!/usr/bin/env python
from sys import argv
from subprocess import Popen, PIPE
import shlex
src_file="/home/orsosa/Clas12Tool/Utils/writeROOT.cc"
res0 = Popen(shlex.split('grep "hipo::node" ' + src_file),stdout=PIPE)
r1 = Popen(shlex.split('grep "reader.getBranch" ' + src_file),stdout=PIPE)
res1 = Popen(shlex.split("sed 's/reader\./fReader->/g'"),stdin=r1.stdout,stdout=PIPE)
r1.stdout.close()
fout0 = open("include/node_declaration.h","w")
fout1 = open("src/node_assignment.cxx","w")
fout1.write("//// File automatically produced by format_hiponodes.py do not make changes here!!\n")
fout1.write('#include "TIdentificatorCLAS12.h"\n')
fout1.write("int TIdentificatorCLAS12::InitNodes()\n")
fout1.write("{\n")
fout0.write("//// File automatically produced by format_hiponodes.py do not make changes here!!\n")


for line in res0.communicate():
    if line is not None: fout0.write(line)

for line in res1.communicate():
    if line is not None: fout1.write(line)
fout1.write("}\n")
