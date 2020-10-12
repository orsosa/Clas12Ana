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
   hipo::bank  rconfig(factory.getSchema("RUN::config"));
   hipo::bank  lund(factory.getSchema("MC::Lund"));

   hipo::event event;

   int counter = 0;
   std::cout<<"detector:layer \n\n";
     
   while(reader.next()==true){
     reader.read(event);
     event.getStructure(rconfig);
     event.getStructure(lund);
	       
     int en, rn;
     int nrows =0;

     if (rconfig.getRows()>0){
       en = rconfig.getInt("event",0);
       rn = rconfig.getInt("run",0);
     }
     else 
       continue;

     std::cout<<"run: "<<rn<<" -- event: "<<en<<std::endl;
     nrows = lund.getRows();
     std::cout<<"MC::Lund (row, type, pid, e, px, py, pz)\n";
     for(int row = 0; row < nrows; row++){	
       int   pid = lund.getInt("pid",row);
       int   type = lund.getInt("type",row);
       float px = lund.getFloat("px",row);
       float py = lund.getFloat("py",row);
       float pz = lund.getFloat("pz",row);
       float e = lund.getFloat("energy",row);
       std::cout<<"("<<row<<", "<<type<<", "<<pid<<", "<<e<<", "<<px<<", "<<py<<", "<<pz<<") ";
     }
     std::cout<<std::endl;
     counter++;
     
     //     if (nrows>0) printf("\n");
   }

   printf("processed events = %d\n",counter);
}
