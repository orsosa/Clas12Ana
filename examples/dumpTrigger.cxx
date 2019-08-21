//******************************************************************
//*  ██╗  ██╗██╗██████╗  ██████╗     ██╗  ██╗    ██████╗
//*  ██║  ██║██║██╔══██╗██╔═══██╗    ██║  ██║   ██╔═████╗
//*  ███████║██║██████╔╝██║   ██║    ███████║   ██║██╔██║
//*  ██╔══██║██║██╔═══╝ ██║   ██║    ╚════██║   ████╔╝██║
//*  ██║  ██║██║██║     ╚██████╔╝         ██║██╗╚██████╔╝
//*  ╚═╝  ╚═╝╚═╝╚═╝      ╚═════╝          ╚═╝╚═╝ ╚═════╝
//************************ Jefferson National Lab (2017) ***********
//******************************************************************
//* Program for reading HIPO-4 Schema from file.
//* Reads the file and output the schema.
//*--
//* Author: O. Soto
//*

#include <cstdlib>
#include <iostream>
#include "reader.h"
#include "TH1F.h"
#include "TFile.h"


int main(int argc, char** argv) {

   std::cout << " reading file (HIPO) "  << __cplusplus << std::endl;

   char inputFile[256];

   if(argc>1) {
      sprintf(inputFile,"%s",argv[1]);
      //sprintf(outputFile,"%s",argv[2]);
   } else {
      std::cout << " *** please provide a file name..." << std::endl;
     exit(0);
   }

   hipo::reader  reader;
   reader.open(inputFile);
   hipo::dictionary  factory;
   reader.readDictionary(factory);
   // getchar();
   hipo::bank  part(factory.getSchema("REC::Particle"));
   hipo::bank  rconfig(factory.getSchema("RUN::config"));
   hipo::bank  recevent(factory.getSchema("REC::Event"));

   hipo::event event;

   int counter = 0;
   float Ne=0,Npim=0,Npip=0;     
   int Nn=0,Np=0,Npip_s=0,Npim_s=0,Nd=0,Nkp=0,Nkm=0;
   while(reader.next()==true){
     reader.read(event);
     event.getStructure(recevent);
     event.getStructure(part);
     event.getStructure(rconfig);

     int en, rn;
     int nrows =0;

     if (rconfig.getRows()>0){
       en = rconfig.getInt("event",0);
       rn = rconfig.getInt("run",0);
     }
     else 
       continue;

     nrows = part.getRows();
     if (nrows >=0){
       int pid = part.getInt("pid",0);   
     if (pid==11)
       Ne++;
     if (pid==-211)
       Npim++;
     if (pid==211)
       Npip++;
     }

     for (int k = 1; k<nrows;k++){
       int pid  = part.getInt("pid",0);
       if(pid==211) Npip_s++;
       if(pid==-211) Npim_s++;
       if(pid==2112) Nn++;
       if(pid==2212) Np++;
       if(pid==321) Nkp++;
       if(pid==-321) Nkm++;
       if(pid==45) Nd++;
     }

     counter++;
     
     //     if (nrows>0) printf("\n");
   }
   std::cout<<"events Ne/Npim/Npip ---- "<<Ne<<"/"<<Npim<<"/"<<Npip<<std::endl;
   std::cout<<"fract. Ne/Npim/Npip ---- "<<Ne/counter<<"/"<<Npim/counter<<"/"<<Npip/counter<<std::endl;
   std::cout<<"second events Npip/Npim/Nn/Np/Nkp/Nkm/Nd ---- "<<Npip_s<<"/"<<Npim_s<<"/"<<Nn<<"/"<<Np<<"/"<<Nn<<"/"<<Nkp<<"/"<<Nkm<<"/"<<Nd<<std::endl;
   std::cout<<"second frac. Npip/Npim/Nn/Np/Nkp/Nkm/Nd ---- "<<Npip_s/counter<<"/"<<Npim_s/counter<<"/"<<Nn/counter<<"/"<<Np/counter<<"/"<<Nn/counter<<"/"<<Nkp/counter<<"/"<<Nkm/counter<<"/"<<Nd/counter<<std::endl;
   printf("processed events = %d\n",counter);
}
