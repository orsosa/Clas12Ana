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

int rotate_dcxy(Float_t dcx,Float_t dcy,Float_t &dcx_rot,Float_t &dcy_rot);

int main(int argc, char **argv)
{
  //gROOT->Reset();
  TBenchmark bm;
  bm.Start(argv[0]);
  bool simul_key = 0;

  if (argc<2)
  {
    std::cout<<"you must supply one hipo file.\n";
    exit(1);
  }
  TString fname = argv[1];
  for (int k=2;k<argc;k++){fname=fname + " " + argv[k];}
  TDatabasePDG pdg;
  Double_t kMe =pdg.GetParticle(11)->Mass();
  const char* NtupleName;

  TString  VarList = "TargType:Q2:Nu:Xb:W:SectorEl:ThetaPQ:PhiPQ:Zh:Pt:W2p:Xf:T:P:T4:deltaZ:E:Ee:Pe:Ect:Sct:Ecr:Scr:evnt:Px:Py:Pz:Xe:Ye:Ze:Xec:Yec:Zec:TEc:DCX:DCY:DCZ:Pex:Pey:Pez:Ein:Eout:Eine:Eoute:pid:Beta:vxh:vyh:vzh:npheltcc:nphehtcc:e_npheltcc:e_nphehtcc:e_chi2pid:chi2pid:e_Epcal:Epcal:e_sector_ltcc:e_sector_htcc:e_sector_ecal:sector_ltcc:sector_htcc:sector_ecal:helic:e_pcal_lu:e_pcal_lv:e_pcal_lw:e_ecin_lu:e_ecin_lv:e_ecin_lw:e_ecout_lu:e_ecout_lv:e_ecout_lw:pcal_lu:pcal_lv:pcal_lw:ecin_lu:ecin_lv:ecin_lw:ecout_lu:ecout_lv:ecout_lw:e_pcal_hx:e_pcal_hy:e_pcal_hz:e_ecin_hx:e_ecin_hy:e_ecin_hz:e_ecout_hx:e_ecout_hy:e_ecout_hz:sector_dc:statPart:e_statPart:e_DCPx:e_DCPy:e_DCPz:DCPx:DCPy:DCPz:trajx_sl0:trajx_sl1:trajx_sl2:trajx_sl3:trajx_sl4:trajx_sl5:trajy_sl0:trajy_sl1:trajy_sl2:trajy_sl3:trajy_sl4:trajy_sl5:trajz_sl0:trajz_sl1:trajz_sl2:trajz_sl3:trajz_sl4:trajz_sl5:trajdcxr0:trajdcxr1:trajdcxr2:trajdcyr0:trajdcyr1:trajdcyr2:trajdczr0:trajdczr1:trajdczr2:e_trajdcxr0:e_trajdcxr1:e_trajdcxr2:e_trajdcyr0:e_trajdcyr1:e_trajdcyr2:e_trajdczr0:e_trajdczr1:e_trajdczr2:e_pathtof:e_timetof:pathtof:timetof:e_sector_tof:sector_tof:e_Beta:STTime:RFTime:e_dcx_rot_0:e_dcy_rot_0:e_dcx_rot_1:e_dcy_rot_1:e_dcx_rot_2:e_dcy_rot_2:dcx_rot_0:dcy_rot_0:dcx_rot_1:dcy_rot_1:dcx_rot_2:dcy_rot_2:rich_h_x:rich_h_y:rich_h_z:rich_h_t:rich_c_x:rich_c_y:rich_c_z:rich_c_t:rich_rr_x:rich_rr_y:rich_rr_z:rich_rr_hx:rich_rr_hy:rich_rr_hz";

  Int_t Nvar = VarList.CountChar(':')+1;
 
  Float_t *vars = new Float_t[Nvar];
  TVector3 *vert;
  TIdentificatorCLAS12 *t = new TIdentificatorCLAS12(fname);
  
  TFile *output;

  if(simul_key == 0) {
    NtupleName = "ntuple_data";
    output = new TFile("outfiles/prune_data_test.root", "RECREATE", "Data of particles");
  } else { 
    NtupleName = "ntuple_accept";
    output = new TFile("outfiles/prune_simul.root", "RECREATE", "Data of particles");
  }

  TNtuple *tElec = new TNtuple("e_rec","All Electrons","Q2:W:Nu:vxec:vyec:vzec:vxe:vye:vze:Pex:Pey:Pez:event:P:E:Ein:Eout:Epcal:npheltcc:nphehtcc:helic:e_pcal_lu:e_pcal_lv:e_pcal_lw:e_ecin_lu:e_ecin_lv:e_ecin_lw:e_ecout_lu:e_ecout_lv:e_ecout_lw:e_pcal_hx:e_pcal_hy:e_pcal_hz:e_ecin_hx:e_ecin_hy:e_ecin_hz:e_ecout_hx:e_ecout_hy:e_ecout_hz:e_trajdcxr0:e_trajdcxr1:e_trajdcxr2:e_trajdcyr0:e_trajdcyr1:e_trajdcyr2:e_trajdczr0:e_trajdczr1:e_trajdczr2:e_pathtof:e_timetof:e_sector_tof:e_Beta:STTime:RFTime:e_dcx_rot_0:e_dcy_rot_0:e_dcx_rot_1:e_dcy_rot_1:e_dcx_rot_2:e_dcy_rot_2:e_sector_ltcc:e_sector_htcc:e_sector_ecal:rich_h_x:rich_h_y:rich_h_z:rich_h_t:rich_c_x:rich_c_y:rich_c_z:rich_c_t:rich_rr_x:rich_rr_y:rich_rr_z:rich_rr_hx:rich_rr_hy:rich_rr_hz");

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

  while (t->Next())
  {
      //    t->PrintMaps();
    Int_t nRows = t->GetNRows();
    //    const char * tt = "C";
    //if(nRows>0 && (t->GetCategorization(0,tt)) == "electron" && t -> Q2() > 1. && t -> W() > 2. && t -> Nu() / 5.015 < 0.85)
    if(nRows>0 && (t->GetCategorization(0)) == "pi-")  
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
      
      DataElec[39] = t->TrajDCX(0,0);
      DataElec[40] = t->TrajDCX(0,1);
      DataElec[41] = t->TrajDCX(0,2);
      DataElec[42] = t->TrajDCY(0,0);
      DataElec[43] = t->TrajDCY(0,1);
      DataElec[44] = t->TrajDCY(0,2);
      DataElec[45] = t->TrajDCZ(0,0);
      DataElec[46] = t->TrajDCZ(0,1);
      DataElec[47] = t->TrajDCZ(0,2);

      DataElec[48] = t->PathTOF(0);
      DataElec[49] = t->TimeTOF(0);
      DataElec[50] = t->SectorTOF(0);
      DataElec[51] = t->Beta(0);
      DataElec[52] = t->STTime();
      DataElec[53] = t->RFTime();

      
      Float_t dcx,dcy,dcx_rot,dcy_rot,dcth,DCsec;
      dcx   = t->TrajDCX(0,0);
      dcy   = t->TrajDCY(0,0);
      rotate_dcxy(dcx,dcy,dcx_rot,dcy_rot);
      DataElec[54] = dcx_rot;//region 0
      DataElec[55] = dcy_rot;//region 0
      dcx   = t->TrajDCX(0,1);
      dcy   = t->TrajDCY(0,1);
      rotate_dcxy(dcx,dcy,dcx_rot,dcy_rot);
      DataElec[56] = dcx_rot;//region 1
      DataElec[57] = dcy_rot;//region 1
      dcx   = t->TrajDCX(0,2);
      dcy   = t->TrajDCY(0,2);
      rotate_dcxy(dcx,dcy,dcx_rot,dcy_rot);
      DataElec[58] = dcx_rot;//region 2
      DataElec[59] = dcy_rot;//region 2

      DataElec[60] = t->SectorLTCC(0);
      DataElec[61] = t->SectorHTCC(0);
      DataElec[62] = t->SectorECAL(0);

      DataElec[63] = t->RICH_HAD_X();
      DataElec[64] = t->RICH_HAD_Y();
      DataElec[65] = t->RICH_HAD_Z();
      DataElec[66] = t->RICH_HAD_T();
      DataElec[67] = t->RICH_CLUSTER_X();
      DataElec[68] = t->RICH_CLUSTER_Y();
      DataElec[69] = t->RICH_CLUSTER_Z();
      DataElec[70] = t->RICH_CLUSTER_T();

      DataElec[71] = t->RICH_RR_X();
      DataElec[72] = t->RICH_RR_Y();
      DataElec[73] = t->RICH_RR_Z();
      DataElec[74] = t->RICH_RR_HX();
      DataElec[75] = t->RICH_RR_HY();
      DataElec[76] = t->RICH_RR_HZ();
      
      
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
          vars[34] = t->VX_DC(i);//t->XEC(i);
          vars[35] = t->VY_DC(i);//t->YEC(i);
          vars[36] = t->VZ_DC(i);//t->ZEC(i);
          vars[37] = t->Px(0);
          vars[38] = t->Py(0);
          vars[39] = t->Pz(0);

          vars[40] = t->Ein(i);
          vars[41] = t->Eout(i);
          vars[42] = t->Ein(0);
          vars[43] = t->Eout(0);
	  vars[44] = t->Pid(i) ;
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
	  //vars[91] = t->SectorDC(i);  // Warning seg. fault. on sim files, must be fixed!!!!!!!!!1
	  vars[92] = t->Status(i);
	  vars[93] = t->Status(0);

	  vars[94] = t->Px_DC(0);
	  vars[95] = t->Py_DC(0);
	  vars[96] = t->Pz_DC(0);
	  vars[97] = t->Px_DC(i);
	  vars[98] = t->Py_DC(i);
	  vars[99] = t->Pz_DC(i);

	  vars[100] = t->TrajX(i,0);
	  vars[101] = t->TrajX(i,1);
	  vars[102] = t->TrajX(i,2);
	  vars[103] = t->TrajX(i,3);
	  vars[104] = t->TrajX(i,4);
	  vars[105] = t->TrajX(i,5);
	  vars[106] = t->TrajY(i,0);
	  vars[107] = t->TrajY(i,1);
	  vars[108] = t->TrajY(i,2);
	  vars[109] = t->TrajY(i,3);
	  vars[110] = t->TrajY(i,4);
	  vars[111] = t->TrajY(i,5);
	  vars[112] = t->TrajZ(i,0);
	  vars[113] = t->TrajZ(i,1);
	  vars[114] = t->TrajZ(i,2);
	  vars[115] = t->TrajZ(i,3);
	  vars[116] = t->TrajZ(i,4);
	  vars[117] = t->TrajZ(i,5);

	  vars[118] = t->TrajDCX(i,0);
	  vars[119] = t->TrajDCX(i,1);
	  vars[120] = t->TrajDCX(i,2);
	  vars[121] = t->TrajDCY(i,0);
	  vars[122] = t->TrajDCY(i,1);
	  vars[123] = t->TrajDCY(i,2);
	  vars[124] = t->TrajDCZ(i,0);
	  vars[125] = t->TrajDCZ(i,1);
	  vars[126] = t->TrajDCZ(i,2);

	  vars[127] = t->TrajDCX(0,0);
	  vars[128] = t->TrajDCX(0,1);
	  vars[129] = t->TrajDCX(0,2);
	  vars[130] = t->TrajDCY(0,0);
	  vars[131] = t->TrajDCY(0,1);
	  vars[132] = t->TrajDCY(0,2);
	  vars[133] = t->TrajDCZ(0,0);
	  vars[134] = t->TrajDCZ(0,1);
	  vars[135] = t->TrajDCZ(0,2);

	  vars[136] = t->PathTOF(0);
	  vars[137] = t->TimeTOF(0);
	  vars[138] = t->PathTOF(i);
	  vars[139] = t->TimeTOF(i);

	  vars[140] = t->SectorTOF(0);
	  vars[141] = t->SectorTOF(i);

	  vars[142] = t->Beta(0);
	  vars[143] = t->STTime();
	  vars[144] = t->RFTime();


	  Float_t dcx,dcy,dcx_rot,dcy_rot,dcth,DCsec;
	  dcx   = t->TrajDCX(0,0);
	  dcy   = t->TrajDCY(0,0);
	  rotate_dcxy(dcx,dcy,dcx_rot,dcy_rot);
	  vars[145] = dcx_rot;//region 0
	  vars[146] = dcy_rot;//region 0
	  dcx   = t->TrajDCX(0,1);
	  dcy   = t->TrajDCY(0,1);
	  rotate_dcxy(dcx,dcy,dcx_rot,dcy_rot);
	  vars[147] = dcx_rot;//region 1
	  vars[148] = dcy_rot;//region 1
	  dcx   = t->TrajDCX(0,2);
	  dcy   = t->TrajDCY(0,2);
	  rotate_dcxy(dcx,dcy,dcx_rot,dcy_rot);
	  vars[149] = dcx_rot;//region 2
	  vars[150] = dcy_rot;//region 2

	  dcx   = t->TrajDCX(i,0);
	  dcy   = t->TrajDCY(i,0);
	  rotate_dcxy(dcx,dcy,dcx_rot,dcy_rot);
	  vars[151] = dcx_rot;//region 0
	  vars[152] = dcy_rot;//region 0
	  dcx   = t->TrajDCX(i,1);
	  dcy   = t->TrajDCY(i,1);
	  rotate_dcxy(dcx,dcy,dcx_rot,dcy_rot);
	  vars[153] = dcx_rot;//region 1
	  vars[154] = dcy_rot;//region 1
	  dcx   = t->TrajDCX(i,2);
	  dcy   = t->TrajDCY(i,2);
	  rotate_dcxy(dcx,dcy,dcx_rot,dcy_rot);
	  vars[155] = dcx_rot;//region 2
	  vars[156] = dcy_rot;//region 2


	  vars[157] = t->RICH_HAD_X(i);
	  vars[158] = t->RICH_HAD_Y(i);
	  vars[159] = t->RICH_HAD_Z(i);
	  vars[160] = t->RICH_HAD_T(i);
	  vars[161] = t->RICH_CLUSTER_X(i);
	  vars[162] = t->RICH_CLUSTER_Y(i);
	  vars[163] = t->RICH_CLUSTER_Z(i);
	  vars[164] = t->RICH_CLUSTER_T(i);


	  vars[165] = t->RICH_RR_X(i);
	  vars[166] = t->RICH_RR_Y(i);
	  vars[167] = t->RICH_RR_Z(i);
	  vars[168] = t->RICH_RR_HX(i);
	  vars[169] = t->RICH_RR_HY(i);
	  vars[170] = t->RICH_RR_HZ(i);

	  
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

int rotate_dcxy(Float_t dcx,Float_t dcy,Float_t &dcx_rot,Float_t &dcy_rot)
{
  Float_t dcth  = atan2(dcy,dcx)*TMath::RadToDeg();
  Int_t DCsec = (-30<dcth&&dcth<30)*1
	+ (30<dcth&&dcth<90)*2
	+ (90<dcth&&dcth<150)*3
	+ (150<dcth||dcth<-150)*4
	+ (-150<dcth&&dcth<-90)*5
	+ (-90<dcth&&dcth<-30)*6;

  dcx_rot = (DCsec==1)*(cos(TMath::Pi()/2)*dcx - sin(TMath::Pi()/2)*dcy)
  + (DCsec==2)*(cos(TMath::Pi()/6)*dcx - sin(TMath::Pi()/6)*dcy)
  + (DCsec==3)*(cos(-TMath::Pi()/6)*dcx - sin(-TMath::Pi()/6)*dcy)
  + (DCsec==4)*(cos(-TMath::Pi()/2)*dcx - sin(-TMath::Pi()/2)*dcy)
  + (DCsec==5)*(cos(-5*TMath::Pi()/6)*dcx - sin(-5*TMath::Pi()/6)*dcy)
  + (DCsec==6)*(cos(-7*TMath::Pi()/6)*dcx - sin(-7*TMath::Pi()/6)*dcy);
  
  dcy_rot = (DCsec==1)*(cos(TMath::Pi()/2)*dcy + sin(TMath::Pi()/2)*dcx)
  + (DCsec==2)*(cos(TMath::Pi()/6)*dcy + sin(TMath::Pi()/6)*dcx)
  + (DCsec==3)*(cos(-TMath::Pi()/6)*dcy + sin(-TMath::Pi()/6)*dcx)
  + (DCsec==4)*(cos(-TMath::Pi()/2)*dcy + sin(-TMath::Pi()/2)*dcx)
  + (DCsec==5)*(cos(-5*TMath::Pi()/6)*dcy + sin(-5*TMath::Pi()/6)*dcx)
  + (DCsec==6)*(cos(-7*TMath::Pi()/6)*dcy + sin(-7*TMath::Pi()/6)*dcx);
  return 0;
}
