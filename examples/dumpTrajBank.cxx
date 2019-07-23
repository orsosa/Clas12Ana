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
   hipo::bank  traj(factory.getSchema("REC::Traj"));
   hipo::bank  rconfig(factory.getSchema("RUN::config"));
   hipo::bank  recevent(factory.getSchema("REC::Event"));

   hipo::event event;

   int counter = 0;
   std::cout<<"detector:layer \n\n";
     
   while(reader.next()==true){
     reader.read(event);
     event.getStructure(recevent);
     event.getStructure(traj);
     event.getStructure(rconfig);

     int en, rn;
     int nrows =0;

     if (rconfig.getRows()>0){
       en = rconfig.getInt("event",0);
       rn = rconfig.getInt("run",0);
     }
     else 
       continue;

     std::cout<<"run: "<<rn<<" -- event: "<<en<<std::endl;
     nrows = traj.getRows();
     for(int row = 0; row < nrows; row++){	
       int   det = traj.getInt("detector",row);
       int   lay = traj.getInt("layer",row);
       std::cout<<det<<":"<<lay<<" // ";
     }
     std::cout<<std::endl;
     counter++;
     
     //     if (nrows>0) printf("\n");
   }

   printf("processed events = %d\n",counter);
}
