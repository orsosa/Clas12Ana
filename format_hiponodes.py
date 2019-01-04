#!/usr/bin/env python
from sys import argv
fin = open("include/hiponodes.h")
fout0 = open("include/node_declaration.h","w")
fout1 = open("src/node_assignment.cxx","w")
fout1.write("int TIdentificatorCLAS12::InitNodes()\n")
fout1.write("{\n")
for line in fin:
    linearray = line.split("=")
    fout0.write(linearray[0] + ";\n")
    fout1.write("  " + linearray[0].split("*")[1] + " = " + linearray[1])
fout1.write("}\n")
