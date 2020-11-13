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
   int Npip = 0;
   int Npim = 0;     
   float th, Q2, W, Nu;
   
   float pe, pex, pey, pez, q, qx, qy, qz, p, px, py, pz;;

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
       pex =  rec.getFloat("px",0);
       pey =  rec.getFloat("py",0);
       pez =  rec.getFloat("pz",0);
       pid = rec.getInt("pid",0);
       if (pid!=11) continue;
       pe = sqrt(pex*pex + pey*pey + pez*pez);
       th = acos(pez/pe);
       Q2 = 4*10.6*pe*(sin((th/2))*sin((th/2)));
       W = sqrt(-Q2 + 0.93827*0.93827 + 2*(10.6 - pe)*0.93827);
       Nu = 10.6 - pe;
       th  = th*TMath::RadToDeg();
       if (15<th&&th<25&&3<pe&&pe<6&&Q2>1&&W>2)
	 Ne++;
       else
	 continue;
       //       std::cout<<en<<" "<<Q2<<" "<<W<<" "<<pe<<" "<<th<<" "<<nrows<<std::endl;
       qx = -pex;
       qy = -pey;
       qz = 10.6-pez;
       for (int k=1;k<nrows; k++){
	 px = rec.getFloat("px",k);
	 py = rec.getFloat("py",k);
	 pz = rec.getFloat("pz",k);
	 p = sqrt(px*px + py*py + pz*pz);
	 float mx2 = W*W + 0.13957*0.13957 - 2*((Nu+0.93827)*sqrt(p*p + 0.13957*0.13957) - (qx*px + qy*py + qz*pz));
	 //          W*W + 0.13957*0.13957 - 2*((Nu+0.93827)*sqrt(P*P + 0.13957*0.13957) - (-Pex*Px -Pey*Py + (10.6-Pez)*Pz));
	 pid = rec.getInt("pid",k);
	 if (0<mx2&&mx2<10&&(pid==211 || pid==-211)){
	   if (pid==211) {Npip++;

	   }
	   else if (pid==-211){ Npim++;
	     std::cout<<en<<" "<<Q2<<" "<<W<<" "<<pe<<" "<<th<<" "<<nrows<<" "<<mx2<<std::endl;
	   }
	 }
	 
       }
       std::cout<<"\n";
     }
     else continue;




     
   }
   std::cout<<"Ne/Npip/Npim: "<<Ne<<"/"<<Npip<<"/"<<Npim<<std::endl;
   

}
