#!/usr/bin/env python
from sys import argv
fin = open("schemas.txt")
fout0 = open("include/node_declaration.h","w")
fout1 = open("src/node_assignment.cxx","w")
fout1.write("//// File automatically produced by format_hiponodes.py do not make changes here!!\n")
fout1.write('#include "TIdentificatorCLAS12.h"\n')
fout1.write("\n")
fout0.write("//// File automatically produced by create_hiponodes.py do not make changes here!!\n")

for line in fin:
    if 'schema : ' not in line: continue
    
    bname = line.strip().split("schema : {")[-1].split("}{")[0].replace("}","").split("/")[0].replace(":","_")

    branches = line.strip().split("schema : {")[-1].split("}{")[-1].replace("}","").split(",")

    fout0.write("hipo::bank *" + bname + " ;\n")
    fout0.write("int get_"+ bname + '(int row);\n')
    fout1.write("int TIdentificatorCLAS12::get_"+ bname + '(int row){\n')
    for br in branches:
        typ =  br.split("/")[1]
        get = "get"
        if typ == "F":
            typ = "float "
            get += "Float("
        if typ == "S":
            typ = "short "
            get += "Short("
        if typ == "B":
            typ = "short "
            get += "Byte("
        if typ == "D":
            typ = "double "
            get += "Double("
        if typ == "I":
            typ = "int "
            get += "Int("
        if typ == "L":
            typ = "long "
            get += "Long("
            
        br = br.split("/")[0]        
        vname = bname + "_" + br
        fout0.write(typ + " " + vname + " ;\n")

        fout1.write("\t" + vname + " = "+ bname + "->" + get + '"' + br + '",row);\n')
    fout1.write('\treturn 0;\n')
    fout1.write("} \n\n")

    
fin = open("schemas.txt")
fout1.write("int TIdentificatorCLAS12::InitBanks(){\n")
for line in fin:
    if 'schema : ' not in line: continue
    bname = line.strip().split("schema : {")[-1].split("}{")[0].replace("}","").split("/")[0]
    fout1.write('\tif (fFactory->hasSchema("'+ bname + '"))\n');
    fout1.write("\t\t" + bname.replace(":","_") + ' = new hipo::bank(fFactory->getSchema("'+ bname + '"));\n')

fout1.write("return 0;\n}\n\n")

fin = open("schemas.txt")
fout1.write("int TIdentificatorCLAS12::FillBanks(){\n")
for line in fin:
    if 'schema : ' not in line: continue
    bname = line.strip().split("schema : {")[-1].split("}{")[0].replace("}","").split("/")[0]
    fout1.write('\tif (fFactory->hasSchema("'+ bname + '"))\n');
    fout1.write('\t\t fEvent->getStructure(*'+ bname.replace(":","_") + ');\n')

fout1.write("return 0;\n}\n\n")
    
        

