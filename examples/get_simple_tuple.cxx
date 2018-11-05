#include "Riostream.h"
#include "TApplication.h"
#include "TBenchmark.h"
#include "TROOT.h"
#include "TFile.h"
#include "TNtuple.h"
#include "TIdentificatorCLAS12.h"
#include "TMath.h"
#include "TDatabasePDG.h"
#include "TParticlePDG.h"
using namespace std;

int main(int argc, char **argv)
{
  //gROOT->Reset();
  TBenchmark bm;
  bm.Start(argv[0]);
  bool simul_key = 0;

  if (argc!=2)
  {
    std::cout<<"you must supply one hipo file.\n";
    exit(1);
  }
  
  TString fname = argv[1];
  TDatabasePDG pdg;
  Double_t kMe =pdg.GetParticle(11)->Mass();
  const char* NtupleName;
  TString  VarList = "TargType:Q2:Nu:Xb:W:SectorEl:ThetaPQ:PhiPQ:Zh:Pt:W2p:Xf:T:P:T4:deltaZ:E:Ee:Pe:Ect:Sct:Ecr:Scr:evnt:Px:Py:Pz:Xe:Ye:Ze:Xec:Yec:Zec:TEc:ECX:ECY:ECZ:Pex:Pey:Pez:Ein:Eout:Eine:Eoute:pid:Beta:vxh:vyh:vzh:npheltcc:nphehtcc:e_npheltcc:e_nphehtcc:e_chi2pid:chi2pid:e_Epcal:Epcal:e_sector_ltcc:e_sector_htcc:e_sector_ecal:sector_ltcc:sector_htcc:sector_ecal:helic:e_pcal_lu:e_pcal_lv:e_pcal_lw:e_ecin_lu:e_ecin_lv:e_ecin_lw:e_ecout_lu:e_ecout_lv:e_ecout_lw:pcal_lu:pcal_lv:pcal_lw:ecin_lu:ecin_lv:ecin_lw:ecout_lu:ecout_lv:ecout_lw:e_pcal_hx:e_pcal_hy:e_pcal_hz:e_ecin_hx:e_ecin_hy:e_ecin_hz:e_ecout_hx:e_ecout_hy:e_ecout_hz";
  Int_t Nvar = VarList.CountChar(':')+1;
 
  Float_t *vars = new Float_t[Nvar];
  TVector3 *vert;
  TIdentificatorCLAS12 *t = new TIdentificatorCLAS12(fname);
  
  //  Long_t nEntries = 0;
  //t->GetEntries();

  TFile *output;
  
  if(simul_key == 0) {
    NtupleName = "ntuple_data";
    output = new TFile("outfiles/prune_data_test.root", "RECREATE", "Data of particles");
  } else { 
    NtupleName = "ntuple_accept";
    output = new TFile("outfiles/prune_simul.root", "RECREATE", "Data of particles");
  }


  
  TNtuple *tElec = new TNtuple("e_rec","All Electrons","Q2:W:Nu:vxec:vyec:vzec:vxe:vye:vze:Pex:Pey:Pez:event:P:E:Ein:Eout:Epcal:npheltcc:nphehtcc:helic:e_pcal_lu:e_pcal_lv:e_pcal_lw:e_ecin_lu:e_ecin_lv:e_ecin_lw:e_ecout_lu:e_ecout_lv:e_ecout_lw:e_pcal_hx:e_pcal_hy:e_pcal_hz:e_ecin_hx:e_ecin_hy:e_ecin_hz:e_ecout_hx:e_ecout_hy:e_ecout_hz");
  Float_t DataElec[tElec->GetNvar()];

  TNtuple *ntuple = new TNtuple(NtupleName,"stable particles",VarList);
  TNtuple *ntuple_thrown = 0;
  TNtuple *e_thrown=0;
  if(simul_key == 1) {
    ntuple_thrown = new TNtuple("ntuple_thrown","particles pluses",VarList);
    e_thrown = new TNtuple("e_thrown","All Electrons","Q2:W:Nu:vxec:vyec:vzec:vxe:vye:vze:Pex:Pey:Pez:event");
}

//  TH1F *ht = new TH1F("ht","tdiff",1000,-15,15); 
  cout.width(4);
  Int_t event=0;
  while (t->Next()==true)
  {
      //    t->PrintMaps();
    Int_t nRows = t->GetNRows();
    //    const char * tt = "C";
    //if(nRows>0 && (t->GetCategorization(0,tt)) == "electron" && t -> Q2() > 1. && t -> W() > 2. && t -> Nu() / 5.015 < 0.85)
    if(nRows>0 && (t->GetCategorization(0)) == "electron")  
    {
      //Q2:W:Nu:vxec:vyec:vzec:vxe:vye:vze:Pex:Pey:Pez:event
      DataElec[0] = t -> Q2();
      DataElec[1] = t -> W();
      DataElec[2] = t -> Nu();
      vert = t->GetCorrectedVert();
      Float_t vxec=vert->X(); 
      Float_t vyec=vert->Y(); 
      Float_t vzec=vert->Z(); 
      DataElec[3] = vxec; 
      DataElec[4] = vyec; 
      DataElec[5] = vzec;

      DataElec[6] = t->X(0);
      DataElec[7] = t->Y(0);
      DataElec[8] = t->Z(0);
      DataElec[9] = t -> Px(0);
      DataElec[10] = t -> Py(0);
      DataElec[11] = t -> Pz(0);
      DataElec[12] = event;

      DataElec[13] = t->Momentum(0);
      DataElec[14] = t->Etot(0);
      DataElec[15] = t->Ein(0);
      DataElec[16] = t->Eout(0);
      DataElec[17] = t->Epcal(0);
      DataElec[18] = t->NpheLTCC(0);
      DataElec[19] = t->NpheHTCC(0);
      DataElec[20] = t->Helic();

      DataElec[21] = t->LU_PCAL();
      DataElec[22] = t->LV_PCAL();
      DataElec[23] = t->LW_PCAL();
      DataElec[24] = t->LU_ECIN();
      DataElec[25] = t->LV_ECIN();
      DataElec[26] = t->LW_ECIN();
      DataElec[27] = t->LU_ECOUT();
      DataElec[28] = t->LV_ECOUT();
      DataElec[29] = t->LW_ECOUT();

      DataElec[30] = t->HX_PCAL();
      DataElec[31] = t->HY_PCAL();
      DataElec[32] = t->HZ_PCAL();
      DataElec[33] = t->HX_ECIN();
      DataElec[34] = t->HY_ECIN();
      DataElec[35] = t->HZ_ECIN();
      DataElec[36] = t->HX_ECOUT();
      DataElec[37] = t->HY_ECOUT();
      DataElec[38] = t->HZ_ECOUT();

      
      tElec->Fill(DataElec);

      Int_t NmbPion = 0;
      for (Int_t i = 1; i < nRows; i++) 
      {
	//	std::cout<<"\tnr: "<<i<<std::endl; 
      	TString category = t->GetCategorization(i);

	//      	if (category == "gamma" || category == "pi-" || category == "high energy pion +" || category == "low energy pion +" || category == "s_electron" || category == "positron") 
      	if (category == "gamma" || category == "pi-" || category == "pi+")

	{
	  vars[0] = 0;//t -> ElecVertTarg();
	  vars[1] = t -> Q2();
	  vars[2] = t -> Nu();
	  vars[3] = t -> Xb();
	  vars[4] = t -> W();
	  vars[5] = t -> Sector(0);
	  vars[6] = t -> ThetaPQ(i);
	  vars[7] = t -> PhiPQ(i);
	  vars[8] = t -> Zh(i);
	  vars[9] = TMath::Sqrt(t -> Pt2(i));
	  vars[10] = t -> Mx2(i);
	  vars[11] = t -> Xf(i);
	  vars[12] = t -> T(i);
	  vars[13] = t -> Momentum(i);
	  vars[14] = 0;//t -> TimeCorr4(0.139570,i);
	  vars[15] = (t -> Z(i)) - (t -> Z(0));
	  vars[16] = t->Etot(i);//,t->Ein(i)+t->Eout(i));
	  vars[17] = t->Etot(0);//,t->Ein(0)+t->Eout(0));
          vars[18] = t->Momentum(0);
          vars[19] = 0;//t->TimeEC(0);
          vars[20] = 0;//t->TimeSC(0);
          vars[21] = 0;//t->PathEC(0);
          vars[22] = 0;//t->PathSC(0);
          vars[23] = event;
          vars[24] = t->Px(i);
          vars[25] = t->Py(i);
          vars[26] = t->Pz(i);
          vars[27] = t->X(0);
          vars[28] = t->Y(0);
          vars[29] = t->Z(0);
          vert = t->GetCorrectedVert();
          vars[30] = vert->X(); 
          vars[31] = vert->Y(); 
          vars[32] = vert->Z(); 
          vars[33] = 0;//t->TimeEC(i);
          vars[34] = 0;//t->XEC(i);
          vars[35] = 0;//t->YEC(i);
          vars[36] = 0;//t->ZEC(i);
          vars[37] = t->Px(0);
          vars[38] = t->Py(0);
          vars[39] = t->Pz(0);

          vars[40] = t->Ein(i);
          vars[41] = t->Eout(i);
          vars[42] = t->Ein(0);
          vars[43] = t->Eout(0);
	  vars[44] = t->Pid(i);
	  vars[45] = t->Beta(i);
          vars[46] = t->X(i);
          vars[47] = t->Y(i);
          vars[48] = t->Z(i);
	  vars[49] = t->NpheLTCC(i);
	  vars[50] = t->NpheHTCC(i);
	  vars[51] = t->NpheLTCC(0);
	  vars[52] = t->NpheHTCC(0);
	  vars[53] = t->Chi2pid(0);
	  vars[54] = t->Chi2pid(i);
          vars[55] = t->Epcal(0);
	  vars[56] = t->Epcal(i);
	  vars[57] = t->SectorLTCC(0);
	  vars[58] = t->SectorHTCC(0);
	  vars[59] = t->SectorECAL(0);
	  vars[60] = t->SectorLTCC(i);
	  vars[61] = t->SectorHTCC(i);
	  vars[62] = t->SectorECAL(i);
	  vars[63] = t->Helic();

	  vars[64] = t->LU_PCAL();
	  vars[65] = t->LV_PCAL();
	  vars[66] = t->LW_PCAL();
	  vars[67] = t->LU_ECIN();
	  vars[68] = t->LV_ECIN();
	  vars[69] = t->LW_ECIN();
	  vars[70] = t->LU_ECOUT();
	  vars[71] = t->LV_ECOUT();
	  vars[72] = t->LW_ECOUT();

	  vars[73] = t->LU_PCAL(i);
	  vars[74] = t->LV_PCAL(i);
	  vars[75] = t->LW_PCAL(i);
	  vars[76] = t->LU_ECIN(i);
	  vars[77] = t->LV_ECIN(i);
	  vars[78] = t->LW_ECIN(i);
	  vars[79] = t->LU_ECOUT(i);
	  vars[80] = t->LV_ECOUT(i);
	  vars[81] = t->LW_ECOUT(i);
  
	  vars[82] = t->HX_PCAL();
	  vars[83] = t->HY_PCAL();
	  vars[84] = t->HZ_PCAL();
	  vars[85] = t->HX_ECIN();
	  vars[86] = t->HY_ECIN();
	  vars[87] = t->HZ_ECIN();
	  vars[88] = t->HX_ECOUT();
	  vars[89] = t->HY_ECOUT();
	  vars[90] = t->HZ_ECOUT();

	  ntuple->Fill(vars);
	}
      }
    }

    
    /*
    if(false)//simul_key == 1 && t -> Id(0,1)==11) 
    {
      
      
      DataElec[0] = t -> Q2(1);
      DataElec[1] = t -> W(1);
      DataElec[2] = t -> Nu(1);
      DataElec[3] = t -> Z(0,1); 
      DataElec[4] = t -> Px(0,1);
      DataElec[5] = t -> Py(0,1);
      DataElec[6] = t -> Pz(0,1);
      DataElec[7] = 0;

      e_thrown->Fill(DataElec);

      //      std::cout<<"got electron gsim"<<std::endl;
      Int_t NmbPion = 0;
      for(t->Next()) 

      {
	
      	if(t -> Id(i,1)==22 || t -> Id(i,1)==-211 || t -> Id(i,1)==211 ) //gamma: 1/22, pi0,+,-: 7/111,8/211,9 (Geant3/pdg)
        {
	        vars[0] = t -> ElecVertTarg(1);
	        vars[1] = t -> Q2(1);
	        vars[2] = t -> Nu(1);
	        vars[3] = t -> Xb(1);
	        vars[4] = t -> W(1);
	        vars[5] = t -> Sector(0,1);
	        vars[6] = t -> ThetaPQ(i,1);
	        vars[7] = t -> PhiPQ(i,1);
	        vars[8] = t -> Zh(i,1);
	        vars[9] = TMath::Sqrt(t -> Pt2(i,1));
	        vars[10] = t -> Mx2(i,1);
	        vars[11] = t -> Xf(i,1);
	        vars[12] = t -> T(i,1);
	        vars[13] = t -> Momentum(i,1);
	        vars[14] = 0;//t -> TimeCorr4(0.139570,i);
	        vars[15] = (t -> Z(i,1)) - (t -> Z(0,1));
	        vars[16] = t->Momentum(i,1);//TMath::Max(t->Etot(i),t->Ein(i)+t->Eout(i));;
		vars[17] = TMath::Sqrt(t->Momentum(0,1)*t->Momentum(0,1)+kMe*kMe); //TMath::Max(t->Etot(0),t->Ein(0)+t->Eout(0));
		vars[18] =t->Momentum(0,1);
          vars[19] = 0;//t->TimeEC(0);
          vars[20] = 0;//t->TimeSC(0);
          vars[21] = 0;//t->PathEC(0);
          vars[22] = 0;//t->PathSC(0);
          vars[23] = k;
          vars[24] = t->Px(i,1);
          vars[25] = t->Py(i,1);
          vars[26] = t->Pz(i,1);
          vars[27] = t->X(0,1);
          vars[28] = t->Y(0,1);
          vars[29] = t->Z(0,1);
          //vert = t->GetCorrectedVert();
          vars[30] = t->X(0,1);//vert->X(); 
          vars[31] = t->Y(0,1);//vert->Y(); 
          vars[32] = t->Z(0,1);//vert->Z(); 
          vars[33] = 0;//t->TimeEC(i);
          vars[34] = 0;//t->XEC(i);
          vars[35] = 0;//t->YEC(i);
          vars[36] = 0;//t->ZEC(i);
          vars[37] = t->Px(0,1);
          vars[38] = t->Py(0,1);
          vars[39] = t->Pz(0,1);

	  vars[40] = 0;
	  vars[41] = 0;
	  vars[42] = 0;
	  vars[43] = 0;
	  vars[44] = t -> Id(i,1);
	  vars[45] = t->Betta(i,1);
          vars[46] = t->X(i,1);
          vars[47] = t->Y(i,1);
          vars[48] = t->Z(i,1);
    
      ntuple_thrown->Fill(vars);
      
	}
      }
    }*/
    //cout<<std::right<<float(k+1)/nEntries*100<<"%\r";
    //cout.flush();
    //    cout<<std::right<<float(k+1)/nEntries*100<<"%\n";
    
    cout<<std::right<<event++<<"\r";
    cout.flush();
  }
  
  cout<<std::right<<event++<<"\n";
    
  output->Write();
  output->Close();
  cout << "Done." << endl;
  bm.Show(argv[0]);
  return 0;
}
