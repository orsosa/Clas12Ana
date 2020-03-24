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
   hipo::bank  fmt_ht(factory.getSchema("FMTRec::Hits"));
   hipo::bank  fmt_cl(factory.getSchema("FMTRec::Clusters"));
   hipo::bank  fmt_cr(factory.getSchema("FMTRec::Crosses"));

   hipo::event event;

   int counter = 0;
   std::cout<<"detector:layer \n\n";
     
   while(reader.next()==true){
     reader.read(event);
     event.getStructure(rconfig);
     event.getStructure(fmt_ht);
     event.getStructure(fmt_cl);
     event.getStructure(fmt_cr);
	       
     int en, rn;
     int nrows =0;

     if (rconfig.getRows()>0){
       en = rconfig.getInt("event",0);
       rn = rconfig.getInt("run",0);
     }
     else 
       continue;

     std::cout<<"run: "<<rn<<" -- event: "<<en<<std::endl;
     // FMTRec::Hits, some of the variables
     nrows = fmt_ht.getRows();
     std::cout<<"FMTRec::Hits: (ID,fitResidual)\n";
     for(int row = 0; row < nrows; row++){	
       int   id = fmt_ht.getShort("ID",row);
       float fRes = fmt_ht.getFloat("fitResidual",row);
       std::cout<<"("<<id<<", "<<fRes<<") ";
     }
     std::cout<<std::endl;
     // FMTRec::Clusters, some of the variables
     nrows = fmt_cl.getRows();
     std::cout<<"FMTRec::Clusters: (ID,ETot)\n";
     for(int row = 0; row < nrows; row++){	
       int   id = fmt_cl.getShort("ID",row);
       float Etot = fmt_cl.getFloat("ETot",row);
       std::cout<<"("<<id<<", "<<Etot<<") ";
     }
     std::cout<<std::endl;
     // FMTRec::Crosses, some of the variables
     nrows = fmt_cr.getRows();
     std::cout<<"FMTRec::Crosses: (ID,x,y,z)\n";
     for(int row = 0; row < nrows; row++){	
       int   id = fmt_cr.getShort("ID",row);
       float x = fmt_cr.getFloat("x",row);
       float y = fmt_cr.getFloat("y",row);
       float z = fmt_cr.getFloat("z",row);
       std::cout<<"("<<id<<", "<<x<<", "<<y<<", "<<z<<") ";
     }
     std::cout<<std::endl;

     counter++;
     
     //     if (nrows>0) printf("\n");
   }

   printf("processed events = %d\n",counter);
}
