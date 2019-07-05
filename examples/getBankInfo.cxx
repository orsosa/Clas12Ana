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
   getchar();
   hipo::bank  honline(factory.getSchema("HEL::online"));
   hipo::bank  hflip(factory.getSchema("HEL::flip"));
   hipo::bank  rconfig(factory.getSchema("RUN::config"));

   hipo::event event;

   int counter = 0;

   TH1F *hhon = new TH1F("hhon","HEL::online helicity",30,-1.5,1.5);
   TH1F *hhonRaw = new TH1F("hhonRaw","HEL::online helicityRaw",30,-1.5,1.5);      
   TH1F *hf = new TH1F("hf","HEL::flip helicity",30,-1.5,1.5);
   TH1F *hfRaw = new TH1F("hfRaw","HEL::flip helicityRaw",30,-1.5,1.5);
   TH1F *hfpatt = new TH1F("hpatt","HEL::flip pattern",30,-1.5,1.5);
   TH1F *hfpair = new TH1F("hpair","HEL::flip pair",30,-1.5,1.5);
   TH1F *hfstat = new TH1F("hstat","HEL::flip status",30,-1.5,1.5);
   TH1F *hDt = new TH1F("hDt","Flip #Deltat change ms",200,33.4,34.4);

   while(reader.next()==true){
     reader.read(event);
     
     event.getStructure(honline);
     event.getStructure(hflip);
     event.getStructure(rconfig);
     
     int en = rconfig.getInt("event",0);
     int rn = rconfig.getInt("run",0);
     long ts = rconfig.getLong("timestamp",0);

    
     int nrows =0;
     
     nrows = honline.getRows();
     for(int row = 0; row < nrows; row++){	
       int   hel = honline.getByte("helicity",row);
       int   helRaw = honline.getByte("helicityRaw",row); 
       hhon->Fill(hel);
       hhonRaw->Fill(helRaw);
       //  printf("%d,%d: honline %d, honlineRaw %d",rn,en,hel,helRaw);
     }
     //printf("\n");
     
     long hts_ms_prev = 0;
     nrows = hflip.getRows();
     for(int row = 0; row < nrows; row++){	
       int   hel = hflip.getByte("helicity",row);
       int   helRaw = hflip.getByte("helicityRaw",row);
       int   patt = hflip.getByte("pattern",row);
       int   pair = hflip.getByte("pair",row);
       int   hstat = hflip.getByte("status",row);
       int   hev = hflip.getInt("event",row);
       int   hrun = hflip.getInt("run",row);
       long   hts = hflip.getLong("timestamp",row);
       if (counter>1)
	 hDt->Fill(hts/4e6 - hts_ms_prev);
       hts_ms_prev =  hts/4e6;
       hf->Fill(hel);
       hfRaw->Fill(helRaw);
       hfpatt->Fill(patt);
       hfpair->Fill(pair);
       hfstat->Fill(hstat);
       //    printf("%d,%d: hel %d, helRaw %d, patt %d, pair %d, hrun %d, hev %d, hstat %d, hts %ld, ts %ld",rn,en,hel,helRaw,patt,pair,hrun,hev,hstat,hts,ts);
     }
     counter++;
     //     if (nrows>0) printf("\n");
   }
   TFile fout("hist.root","recreate");
   hhon->Write();
   hhonRaw->Write();
   hf->Write();
   hfRaw->Write();
   hfpatt->Write();
   hfpair->Write();
   hfstat->Write();
   hDt->Write();
   printf("processed events = %d\n",counter);
}
