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
#include "TMath.h"


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
   hipo::bank  rec(factory.getSchema("REC::Particle"));

   hipo::event event;

   int counter = 0;
   std::cout<<"detector:layer \n\n";
   int Ne = 0;     
   float th;
   float Q2;
   float W;
   float p, px, py, pz;
   while(reader.next()==true){
     reader.read(event);
     event.getStructure(rconfig);
     event.getStructure(rec);	       
     int en, rn;
     int nrows = 0;
     int pid = 0;
     if (rconfig.getRows()>0){
       en = rconfig.getInt("event",0);
       rn = rconfig.getInt("run",0);
     }       
     else 
       continue;
     if (rec.getRows()>0){
       nrows = rec.getRows();
       px =  rec.getFloat("px",0);
       py =  rec.getFloat("py",0);
       pz =  rec.getFloat("pz",0);
       pid = rec.getInt("pid",0);
     }
     else 
       continue;
     if (pid!=11) continue; 
     
     p = sqrt(px*px + py*py + pz*pz);
     th = acos(pz/p);
     Q2 = 4*10.6*p*(sin((th/2))*sin((th/2)));
     W = -Q2 + 0.93827*0.93827 + 2*(10.6 - p)*0.93827;

     th  = th*TMath::RadToDeg();
     if (15<th&&th<25&&3<p&&p<6&&Q2>1&&W>2)
       std::cout<<en<<" "<<Q2<<" "<<W<<" "<<p<<" "<<th<<" "<<nrows<<std::endl;
   }
   //std::cout<<"Ne: "<<Ne<<std::endl;
   //revent         Q2          W         Pe       th_e 

}
