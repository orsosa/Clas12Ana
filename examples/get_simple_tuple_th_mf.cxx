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
#include "TThread.h"
#include "TSemaphore.h"
#include "TMutex.h"
#include "TList.h"

using namespace std;

float EBEAM=10.6;
TString fname="";
Int_t Nvar, NvarElec;
Int_t Nth = 5;
Long_t Ntotal=0;
Long_t TH_MAX = 10000;
Long_t Nsum=0;
Float_t *progress_th;
Bool_t *slotAvailable, QUIET = false;
Int_t NthActive = 0;
TCondition cond(NULL);
TCondition slotCond(NULL);
TSemaphore sem(0);
TSemaphore initsem(0);
TSemaphore endsem(1);

TMutex *fileMutex;
TMutex *pdgMutex;
TMutex *coutMutex;


std::map <Int_t,TString> thn_ind;

bool simul_key = 0;

Int_t Nt = 0;
TString  varList = "TargType:Q2:Nu:Xb:W:SectorEl:ThetaPQ:PhiPQ:Zh:Pt2:Mx2:Xf:T:P:T4:deltaZ:E:Ee:Pe:Ect:Sct:Ecr:Scr:evnt:Px:Py:Pz:Xe:Ye:Ze:Xec:Yec:Zec:TEc:DCX:DCY:DCZ:Pex:Pey:Pez:Ein:Eout:Eine:Eoute:pid:Beta:vxh:vyh:vzh:npheltcc:nphehtcc:e_npheltcc:e_nphehtcc:e_chi2pid:chi2pid:e_Epcal:Epcal:e_sector_ltcc:e_sector_htcc:e_sector_ecal:sector_ltcc:sector_htcc:sector_ecal:helic:e_pcal_lu:e_pcal_lv:e_pcal_lw:e_ecin_lu:e_ecin_lv:e_ecin_lw:e_ecout_lu:e_ecout_lv:e_ecout_lw:pcal_lu:pcal_lv:pcal_lw:ecin_lu:ecin_lv:ecin_lw:ecout_lu:ecout_lv:ecout_lw:e_pcal_hx:e_pcal_hy:e_pcal_hz:e_ecin_hx:e_ecin_hy:e_ecin_hz:e_ecout_hx:e_ecout_hy:e_ecout_hz:sector_dc:statPart:e_statPart:e_DCPx:e_DCPy:e_DCPz:DCPx:DCPy:DCPz:trajx_sl0:trajx_sl1:trajx_sl2:trajx_sl3:trajx_sl4:trajx_sl5:trajy_sl0:trajy_sl1:trajy_sl2:trajy_sl3:trajy_sl4:trajy_sl5:trajz_sl0:trajz_sl1:trajz_sl2:trajz_sl3:trajz_sl4:trajz_sl5:trajdcxr0:trajdcxr1:trajdcxr2:trajdcyr0:trajdcyr1:trajdcyr2:trajdczr0:trajdczr1:trajdczr2:e_trajdcxr0:e_trajdcxr1:e_trajdcxr2:e_trajdcyr0:e_trajdcyr1:e_trajdcyr2:e_trajdczr0:e_trajdczr1:e_trajdczr2:e_pathtof:e_timetof:pathtof:timetof:e_sector_tof:sector_tof:e_Beta:STTime:RFTime:e_dcx_rot_0:e_dcy_rot_0:e_dcx_rot_1:e_dcy_rot_1:e_dcx_rot_2:e_dcy_rot_2:dcx_rot_0:dcy_rot_0:dcx_rot_1:dcy_rot_1:dcx_rot_2:dcy_rot_2:mcmass:revent:Npip_rec:Npim_rec:Npip_mc:Npim_mc:rec_elec:dc_chi2:ftof1ax:ftof1ay:ftof1az:pcalx:pcaly:pcalz:ecalx:ecaly:ecalz:ltccx:ltccy:ltccz:htccx:htccy:htccz:e_dc_chi2:e_ftof1ax:e_ftof1ay:e_ftof1az:e_pcalx:e_pcaly:e_pcalz:e_ecalx:e_ecaly:e_ecalz:e_ltccx:e_ltccy:e_ltccz:e_htccx:e_htccy:e_htccz:ftof1bx:ftof1by:ftof1bz:ftof2x:ftof2y:ftof2z:e_ftof1bx:e_ftof1by:e_ftof1bz:e_ftof2x:e_ftof2y:e_ftof2z:helonline_hel:helonline_helRaw:helflip_hel:helflip_helRaw:helflip_event:y:th_e:helicRaw:dc_status:dc_ndf:e_dc_status:e_dc_ndf";

TString varListElec = "Q2:W:Nu:vxec:vyec:vzec:vxe:vye:vze:Pex:Pey:Pez:event:P:E:Ein:Eout:Epcal:npheltcc:nphehtcc:helic:e_pcal_lu:e_pcal_lv:e_pcal_lw:e_ecin_lu:e_ecin_lv:e_ecin_lw:e_ecout_lu:e_ecout_lv:e_ecout_lw:e_pcal_hx:e_pcal_hy:e_pcal_hz:e_ecin_hx:e_ecin_hy:e_ecin_hz:e_ecout_hx:e_ecout_hy:e_ecout_hz:e_trajdcxr0:e_trajdcxr1:e_trajdcxr2:e_trajdcyr0:e_trajdcyr1:e_trajdcyr2:e_trajdczr0:e_trajdczr1:e_trajdczr2:e_pathtof:e_timetof:e_sector_tof:e_Beta:STTime:RFTime:e_dcx_rot_0:e_dcy_rot_0:e_dcx_rot_1:e_dcy_rot_1:e_dcx_rot_2:e_dcy_rot_2:e_sector_ltcc:e_sector_htcc:e_sector_ecal:revent:Npip_rec:Npim_rec:Npip_mc:Npim_mc:rec_elec:e_dc_chi2:helonline_hel:helonline_helRaw:helflip_hel:helflip_helRaw:helflip_event:helicRaw:e_dc_status:e_dc_ndf:e_sector_dc";


TDatabasePDG *pdg;


int rotate_dcxy(Float_t dcx,Float_t dcy,Float_t &dcx_rot,Float_t &dcy_rot);
void cleanUp(void *arg);
void filter(void *arg);
void create_threads(void *arg);

int main(int argc, char **argv)
{
  //gROOT->Reset();
  TBenchmark bm;
  bm.Start(argv[0]);

  
  if (argc<2)
  {
    std::cout<<"you must supply at least one hipo file.\n";
    
    exit(1);
  }

  for (int k=1;k<argc;k++)
  {
    if (strcmp(argv[k],"-s")==0) // simulation analysis
      simul_key = 1;
    else if (strcmp(argv[k],"-e")==0){ //energy set
      EBEAM = atof(argv[++k]);
    }
    else if (strcmp(argv[k],"-nth")==0){ //n threads
      Nth = atof(argv[++k]);
    }
    else if (strcmp(argv[k],"-nmax")==0){ //max entries read
      Ntotal = atol(argv[++k]);
    }
    else if (strcmp(argv[k],"-q")==0){ //quiet mode
      QUIET = true;
    }

    else
      fname = fname + " " + argv[k];
  }
  if (!QUIET) std::cout<<"File list: "<<fname<<std::endl;
  if (!QUIET) std::cout<<"Processing "<<((simul_key==1)?"simulations":"data")<<std::endl;
  if (!QUIET) std::cout<<"Beam energy: "<< EBEAM<<" GeV"<<std::endl;
  pdg = new TDatabasePDG();
  
  fileMutex = new TMutex();
  pdgMutex = new TMutex();
  coutMutex = new TMutex();

  Nvar = varList.CountChar(':')+1;
 
  //  cout.width(100);
  coutMutex->Lock();
  TIdentificatorCLAS12 *t = new TIdentificatorCLAS12(fname,EBEAM,true); // March - 19 cooking
  coutMutex->UnLock();
  if (!Ntotal)  Ntotal = t->getNevents();
  if (!QUIET) std::cout<<Ntotal<<std::endl;

  progress_th = new Float_t[Nth];
  memset(progress_th,0,Nth*sizeof(float));
  slotAvailable = new Bool_t[Nth];
  for (int k = 0; k<Nth;k++) slotAvailable[k]=kTRUE;

  Nsum=0;
  TThread *th = new TThread("threads_creator",create_threads,(void *)0);
  th->Run();
  initsem.Wait();
  
  while (NthActive){
    //Bool_t exitFlag = kTRUE;
    //    cout<<"Nth "<<NthActive<<" --- ";
    coutMutex->Lock();
    if (!QUIET) std::cout<<"Nth "<<TThread::Exists()<<" --- ";
    coutMutex->UnLock();
    for (int k = 0; k<Nth;k++){
      coutMutex->Lock();
      if (!QUIET) std::cout<<progress_th[k]/TH_MAX*100<<"\% // ";
      coutMutex->UnLock();
      //      exitFlag &= progress_th[k] == TH_MAX;
    }
    coutMutex->Lock();
    if (!QUIET) std::cout<<"  --> "<<(float)Nsum/Ntotal*100<<"\% \r";
    cout.flush();
    coutMutex->UnLock();
    //if (exitFlag) break;
  }
  endsem.Wait();
  th->Delete();
  delete t;
  delete pdg;
  cout << "\nDone." << endl;
  bm.Show(argv[0]);
  return 0;
}
////////// Threads Creator ///////////
void create_threads(void *arg){
  TThread *th;
  Long_t start_ind[2];
  
  start_ind[0] = 0;  start_ind[1] = 0;
  int k = 0;
  while (start_ind[0]<Ntotal){
    while (NthActive == Nth){
      sleep(1);
    }
    for (k=0;k<Nth;k++){
      if (slotAvailable[k]) break;
    }
    if (k==Nth) continue; 
    //TThread::Ps();
    
    //    th = TThread::GetThread(thn_ind[k]);
    //TThread::Delete(th);

    start_ind[1] = k;
    TThread::Lock();
    progress_th[k] = 0;
    NthActive++;
    Nt++;
    TThread::UnLock();
    if (NthActive==1) initsem.Post();
    thn_ind[k]=Form("th_%d",Nt);
    th = new TThread(thn_ind[k],filter,(void *)Form("%ld %ld",start_ind[0], start_ind[1]) );
    endsem.TryWait();
    th->Run();
    start_ind[0] += TH_MAX;

    sem.Wait();
  }
  //  TThread::Printf("############################################ k:%d",k);
  //  th = TThread::Self();
  //  th->Delete();
  //TThread::Exit();
}
////////////
/// cleanup
void cleanUp(void *arg) {
      // here the user cleanup is done
  TThread::Printf("Cleaned Up");

}
///////// Filter thread ///////////
void filter(void *arg)
{
  //sleep(1);
  TThread::CleanUpPush((void*)&cleanUp, (void*) NULL);
  TVector3 *vert;
  TFile *output;
  Int_t Nt_local=0;
  TString args = ((char*)arg);
  Ssiz_t stt=0;
  TString ss;
  args.Tokenize(ss,stt," ");
  Long_t startcnt = atol(ss.Data());
  args.Tokenize(ss,stt," ");
  Int_t ind = atoi(ss.Data());
  TThread::Lock();
  slotAvailable[ind]=kFALSE;
  TThread::UnLock();
  Nt_local=Nt;
  sem.Post();

  TThread::SetCancelOn(); // enable thread canceling

  TString NtupleName;

  fileMutex->Lock();
  if(simul_key == 0) {
    NtupleName = "ntuple_data";
    output = new TFile(Form("outfiles/pruned_dataH4_%d.root",Nt_local), "RECREATE", "Data of particles");
  } else { 
    NtupleName = "ntuple_accept";
    output = new TFile(Form("outfiles/pruned_simulH4_%d.root",Nt_local), "RECREATE", "Data of particles");
  }
  fileMutex->UnLock();
  //  TObjArray *objArr = new TObjArray(5);
  
  TNtuple *tElec = new TNtuple("e_rec","All Electrons",varListElec);
  //objArr->Add(tElec);
  
  Int_t NvarElec = tElec->GetNvar();
  Int_t NvarElecTh = 0;

  Float_t *DataElecTh = 0;
  
  TNtuple *ntuple = new TNtuple(NtupleName,"stable particles",varList);
  //objArr->Add(ntuple);
  TNtuple *ntuple_thrown = 0;
  TNtuple *e_thrown = 0;
  if(simul_key == 1) {
    ntuple_thrown = new TNtuple("ntuple_thrown","particles pluses",varList);
    //objArr->Add(ntuple_thrown);
    e_thrown = new TNtuple("e_thrown","All Electrons","Q2:W:Nu:Pex:Pey:Pez:vxe:vye:vze:mcmass:Npip_rec:Npim_rec:Npip_mc:Npim_mc:rec_elec:event");
    //objArr->Add(e_thrown);
    NvarElecTh = e_thrown->GetNvar();
    DataElecTh = new Float_t[NvarElecTh];
  }

  Float_t *vars = new Float_t[Nvar];
  Float_t *DataElec = new Float_t[NvarElec];
  coutMutex->Lock();
  TIdentificatorCLAS12 *t = new TIdentificatorCLAS12(fname,EBEAM,true); // March - 19 cooking
  coutMutex->UnLock();
  TThread::Printf("###################################### on thread %s, startcnt: %ld, progress_ind: %d",TThread::Self()->GetName(),startcnt,ind);
  Int_t event=0;
  Int_t rec_elec;
  Int_t Npip_mc = -1;
  Int_t Npim_mc = -1;
  Int_t Npip_rec = -1;
  Int_t Npim_rec = -1;
  for (int k=0;k<startcnt;k++) t->Next();
  while (t->Next())
  {
    Int_t nRows = t->GetNRows();

    if(simul_key == 1) {
      Npip_mc = t->GetNPart(211,1);
      Npim_mc = t->GetNPart(-211,1);
    }
    Npip_rec = t->GetNPart(211);
    Npim_rec = t->GetNPart(-211);
    rec_elec = 0;
    
    //if(false)
    if(nRows>0 && (t->GetCategorization(0)) == "electron")  
    {
      
      rec_elec=1;
      DataElec[0] = t -> Q2();
      DataElec[1] = t -> W();
      DataElec[2] = t -> Nu();
      vert = t->GetCorrectedVert();
      Float_t vxec=vert->X(); 
      Float_t vyec=vert->Y(); 
      Float_t vzec=vert->Z();
      delete vert;
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

      DataElec[39] = t->TrajDetIdX(0,"DC","DCSL1");
      DataElec[40] = t->TrajDetIdX(0,"DC","DCSL3");
      DataElec[41] = t->TrajDetIdX(0,"DC","DCSL5");
      DataElec[42] = t->TrajDetIdY(0,"DC","DCSL1");
      DataElec[43] = t->TrajDetIdY(0,"DC","DCSL3");
      DataElec[44] = t->TrajDetIdY(0,"DC","DCSL5");
      DataElec[45] = t->TrajDetIdZ(0,"DC","DCSL1");
      DataElec[46] = t->TrajDetIdZ(0,"DC","DCSL3");
      DataElec[47] = t->TrajDetIdZ(0,"DC","DCSL5");

      DataElec[48] = t->PathTOF(0);
      DataElec[49] = t->TimeTOF(0);
      DataElec[50] = t->SectorTOF(0);
      DataElec[51] = t->Beta(0);
      DataElec[52] = t->STTime();
      DataElec[53] = t->RFTime();
     
      Float_t dcx,dcy,dcx_rot,dcy_rot;

      dcx = t->TrajDetIdX(0,"DC","DCSL1");
      dcy = t->TrajDetIdY(0,"DC","DCSL1");
      rotate_dcxy(dcx,dcy,dcx_rot,dcy_rot);
      DataElec[54] = dcx_rot;//region 0
      DataElec[55] = dcy_rot;//region 0
      dcx = t->TrajDetIdX(0,"DC","DCSL3");
      dcy = t->TrajDetIdY(0,"DC","DCSL3");
      rotate_dcxy(dcx,dcy,dcx_rot,dcy_rot);
      DataElec[56] = dcx_rot;//region 1
      DataElec[57] = dcy_rot;//region 1
      dcx = t->TrajDetIdX(0,"DC","DCSL5");
      dcy = t->TrajDetIdY(0,"DC","DCSL5");
      rotate_dcxy(dcx,dcy,dcx_rot,dcy_rot);
      DataElec[58] = dcx_rot;//region 2
      DataElec[59] = dcy_rot;//region 2

      DataElec[60] = t->SectorLTCC(0);
      DataElec[61] = t->SectorHTCC(0);
      DataElec[62] = t->SectorECAL(0);
      DataElec[63] = t->Event();
      DataElec[64] = Npip_rec;
      DataElec[65] = Npim_rec;
      DataElec[66] = Npip_mc;
      DataElec[67] = Npim_mc;

      DataElec[68] = rec_elec;
      DataElec[69] = t->DCChi2(0);

      DataElec[70] = t->HelicOnline();
      DataElec[71] = t->HelicOnlineRaw();

      DataElec[72] = t->HelicFlip();
      DataElec[73] = t->HelicFlipRaw();

      DataElec[74] = t->HelicFlipEvent();
      DataElec[75] = t->HelicRaw();

      DataElec[76] = t->DCStatus();
      DataElec[77] = t->DCNDF();

      DataElec[78] = t->SectorDC();
	    
      tElec->Fill(DataElec);
      
      for (Int_t i = 1; i < nRows; i++) 
      {
	
      	TString category = t->GetCategorization(i);
	if (category == "pi-" || category == "pi+" || category == "gamma" || category == "proton" || category == "deuteron" || category=="K+" || category=="K-" || category=="neutron")
	//	if (category == "pi-" || category == "pi+" || category == "gamma")
	{

	  vars[0] = 0;//t -> ElecVertTarg();
	  vars[1] = t -> Q2();
	  vars[2] = t -> Nu();
	  vars[3] = t -> Xb();
	  vars[4] = t -> W();
	  vars[5] = t -> Sector(0);
	  vars[6] = t -> ThetaPQ(i);
	  vars[7] = t -> PhiPQ(i);
	  Float_t pid = t->Pid(i);
	  pdgMutex->Lock();
	  Float_t mass = ( pid == 45 )?1.85756:pdg->GetParticle(pid)->Mass();
	  pdgMutex->UnLock();
	  
	  vars[8] = t -> Zh(i,mass);
	  
	  vars[9] = t -> Pt2(i);
	  vars[10] = t -> Mx2(i,mass);
	  vars[11] = t -> Xf(i,mass);
	  vars[12] = t -> T(i,mass);
	  vars[13] = t -> Momentum(i);
	  vars[14] = 0;//t -> TimeCorr4(mass,i);
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
	  delete vert;
          vars[33] = 0;//t->TimeEC(i);
          vars[34] = 0;//t->VX_DC(i);//t->XEC(i);
          vars[35] = 0;//t->VY_DC(i);//t->YEC(i);
          vars[36] = 0;//t->VZ_DC(i);//t->ZEC(i);
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
	  vars[91] = t->SectorDC(i);  // Warning seg. fault. on sim files, must be fixed!!!!!!!!!
	  vars[92] = t->Status(i);
	  vars[93] = t->Status(0);
	  /*
	  vars[94] = t->Px_DC(0);
	  vars[95] = t->Py_DC(0);
	  vars[96] = t->Pz_DC(0);
	  vars[97] = t->Px_DC(i);
	  vars[98] = t->Py_DC(i);
	  vars[99] = t->Pz_DC(i);
	  */
	  
	  vars[100] = t->TrajDetIdX(i,"DC","DCSL1");
	  vars[101] = t->TrajDetIdX(i,"DC","DCSL2");
	  vars[102] = t->TrajDetIdX(i,"DC","DCSL3");
	  vars[103] = t->TrajDetIdX(i,"DC","DCSL4");
	  vars[104] = t->TrajDetIdX(i,"DC","DCSL5");
	  vars[105] = t->TrajDetIdX(i,"DC","DCSL6");
	  vars[106] = t->TrajDetIdY(i,"DC","DCSL1");
	  vars[107] = t->TrajDetIdY(i,"DC","DCSL2");
	  vars[108] = t->TrajDetIdY(i,"DC","DCSL3");
	  vars[109] = t->TrajDetIdY(i,"DC","DCSL4");
	  vars[110] = t->TrajDetIdY(i,"DC","DCSL5");
	  vars[111] = t->TrajDetIdY(i,"DC","DCSL6");
	  vars[112] = t->TrajDetIdZ(i,"DC","DCSL1");
	  vars[113] = t->TrajDetIdZ(i,"DC","DCSL2");
	  vars[114] = t->TrajDetIdZ(i,"DC","DCSL3");
	  vars[115] = t->TrajDetIdZ(i,"DC","DCSL4");
	  vars[116] = t->TrajDetIdZ(i,"DC","DCSL5");
	  vars[117] = t->TrajDetIdZ(i,"DC","DCSL6");

	  vars[118] = t->TrajDetIdX(i,"DC","DCSL1");
	  vars[119] = t->TrajDetIdX(i,"DC","DCSL3");
	  vars[120] = t->TrajDetIdX(i,"DC","DCSL5");
	  vars[121] = t->TrajDetIdY(i,"DC","DCSL1");
	  vars[122] = t->TrajDetIdY(i,"DC","DCSL3");
	  vars[123] = t->TrajDetIdY(i,"DC","DCSL5");
	  vars[124] = t->TrajDetIdZ(i,"DC","DCSL1");
	  vars[125] = t->TrajDetIdZ(i,"DC","DCSL3");
	  vars[126] = t->TrajDetIdZ(i,"DC","DCSL5");

	  vars[127] = t->TrajDetIdX(0,"DC","DCSL1");
	  vars[128] = t->TrajDetIdX(0,"DC","DCSL3");
	  vars[129] = t->TrajDetIdX(0,"DC","DCSL5");
	  vars[130] = t->TrajDetIdY(0,"DC","DCSL1");
	  vars[131] = t->TrajDetIdY(0,"DC","DCSL3");
	  vars[132] = t->TrajDetIdY(0,"DC","DCSL5");
	  vars[133] = t->TrajDetIdZ(0,"DC","DCSL1");
	  vars[134] = t->TrajDetIdZ(0,"DC","DCSL3");
	  vars[135] = t->TrajDetIdZ(0,"DC","DCSL5");
	  
	  vars[136] = t->PathTOF(0);
	  vars[137] = t->TimeTOF(0);
	  vars[138] = t->PathTOF(i);
	  vars[139] = t->TimeTOF(i);

	  vars[140] = t->SectorTOF(0);
	  vars[141] = t->SectorTOF(i);

	  vars[142] = t->Beta(0);
	  vars[143] = t->STTime();
	  vars[144] = t->RFTime();


	  Float_t dcx,dcy,dcx_rot,dcy_rot;
	  
	  dcx   = t->TrajDetIdX(0,"DC","DCSL1");
	  dcy   = t->TrajDetIdY(0,"DC","DCSL1");
	  rotate_dcxy(dcx,dcy,dcx_rot,dcy_rot);
	  vars[145] = dcx_rot;//region 0
	  vars[146] = dcy_rot;//region 0
	  dcx   = t->TrajDetIdX(0,"DC","DCSL3");
	  dcy   = t->TrajDetIdY(0,"DC","DCSL3");
	  rotate_dcxy(dcx,dcy,dcx_rot,dcy_rot);
	  vars[147] = dcx_rot;//region 1
	  vars[148] = dcy_rot;//region 1
	  dcx   = t->TrajDetIdX(0,"DC","DCSL5");
	  dcy   = t->TrajDetIdY(0,"DC","DCSL5");
	  rotate_dcxy(dcx,dcy,dcx_rot,dcy_rot);
	  vars[149] = dcx_rot;//region 2
	  vars[150] = dcy_rot;//region 2

	  dcx   = t->TrajDetIdX(i,"DC","DCSL1");
	  dcy   = t->TrajDetIdY(i,"DC","DCSL1");
	  rotate_dcxy(dcx,dcy,dcx_rot,dcy_rot);
	  vars[151] = dcx_rot;//region 0
	  vars[152] = dcy_rot;//region 0
	  dcx   = t->TrajDetIdX(i,"DC","DCSL3");
	  dcy   = t->TrajDetIdY(i,"DC","DCSL3");
	  rotate_dcxy(dcx,dcy,dcx_rot,dcy_rot);
	  vars[153] = dcx_rot;//region 1
	  vars[154] = dcy_rot;//region 1
	  dcx   = t->TrajDetIdX(i,"DC","DCSL5");
	  dcy   = t->TrajDetIdY(i,"DC","DCSL5");
	  rotate_dcxy(dcx,dcy,dcx_rot,dcy_rot);
	  vars[155] = dcx_rot;//region 2
	  vars[156] = dcy_rot;//region 2
	  vars[157] = 0;// mc mass

	  vars[158] = t->Event();// event number

	  vars[159] = Npip_rec;
	  vars[160] = Npim_rec;
	  vars[161] = Npip_mc;
	  vars[162] = Npim_mc;

	  vars[163] = rec_elec;
	  vars[164] = t->DCChi2(i);
	  vars[165] = t->TrajFTOF1AX(i);
	  vars[166] = t->TrajFTOF1AY(i);
	  vars[167] = t->TrajFTOF1AZ(i);
   	  vars[168] = t->TrajPCALX(i);
	  vars[169] = t->TrajPCALY(i);
	  vars[170] = t->TrajPCALZ(i);
	  vars[171] = t->TrajECINX(i);
	  vars[172] = t->TrajECINY(i);
	  vars[173] = t->TrajECINZ(i);
	  vars[174] = t->TrajLTCCX(i);
	  vars[175] = t->TrajLTCCY(i);
	  vars[176] = t->TrajLTCCZ(i);
	  vars[177] = t->TrajHTCCX(i);
	  vars[178] = t->TrajHTCCY(i);
	  vars[179] = t->TrajHTCCZ(i);
	  vars[180] = t->DCChi2(0);
	  vars[181] = t->TrajFTOF1AX(0);
	  vars[182] = t->TrajFTOF1AY(0);
	  vars[183] = t->TrajFTOF1AZ(0);
   	  vars[184] = t->TrajPCALX(0);
	  vars[185] = t->TrajPCALY(0);
	  vars[186] = t->TrajPCALZ(0);
	  vars[187] = t->TrajECINX(0);
	  vars[188] = t->TrajECINY(0);
	  vars[189] = t->TrajECINZ(0);
	  vars[190] = t->TrajLTCCX(0);
	  vars[191] = t->TrajLTCCY(0);
	  vars[192] = t->TrajLTCCZ(0);
	  vars[193] = t->TrajHTCCX(0);
	  vars[194] = t->TrajHTCCY(0);
	  vars[195] = t->TrajHTCCZ(0);
	  vars[196] = t->TrajFTOF1BX(i);
	  vars[197] = t->TrajFTOF1BY(i);
	  vars[198] = t->TrajFTOF1BZ(i);
	  vars[199] = t->TrajFTOF2X(i);
	  vars[200] = t->TrajFTOF2Y(i);
	  vars[201] = t->TrajFTOF2Z(i);
	  vars[202] = t->TrajFTOF1BX(0);
	  vars[203] = t->TrajFTOF1BY(0);
	  vars[204] = t->TrajFTOF1BZ(0);
	  vars[205] = t->TrajFTOF2X(0);
	  vars[206] = t->TrajFTOF2Y(0);
	  vars[207] = t->TrajFTOF2Z(0);

	  vars[208] = t->HelicOnline();
	  vars[209] = t->HelicOnlineRaw();
	  vars[210] = t->HelicFlip();
	  vars[211] = t->HelicFlipRaw();
	  vars[212] = t->HelicFlipEvent();

	  vars[213] = t->Nu()/EBEAM;
	  vars[214] = TMath::ACos(t->Pz(0)/t->Momentum(0))*TMath::RadToDeg();
	  vars[215] = t->HelicRaw();
	  vars[216] = t->DCStatus(i);
	  vars[217] = t->DCNDF(i);
	  vars[218] = t->DCStatus(0);
	  vars[219] = t->DCNDF(0);

	  ntuple->Fill(vars);

	}
      }
    }
    nRows = t->GetMCNRows();
    Int_t ind_first = 3;
    if(nRows>3 && (simul_key == 1 && t -> Pid(ind_first,1)==11 && t -> LundType(ind_first)==1) ) 
    {

      DataElecTh[0] = t -> Q2(1);
      DataElecTh[1] = t -> W(1);
      DataElecTh[2] = t -> Nu(1);
      DataElecTh[3] = t -> Px(ind_first,1);
      DataElecTh[4] = t -> Py(ind_first,1);
      DataElecTh[5] = t -> Pz(ind_first,1);
      DataElecTh[6] = t -> X(ind_first,1);
      DataElecTh[7] = t -> Y(ind_first,1);
      DataElecTh[8] = t -> Z(ind_first,1);
      DataElecTh[9] = t->MCMass(ind_first);

      DataElecTh[10] = Npip_rec;
      DataElecTh[11] = Npim_rec;
      DataElecTh[12] = Npip_mc;
      DataElecTh[13] = Npim_mc;

      DataElecTh[14] = rec_elec;
      DataElecTh[15] = event;

      e_thrown->Fill(DataElecTh);
    
      int npip=0,npim=0;
      for(int i=ind_first + 1; i<nRows; i++) 
      {    
      	if((t -> LundType(i) == 1) && (t -> Pid(i,1)==-211 || t -> Pid(i,1)==211|| t -> Pid(i,1) == 22 || t -> Pid(i,1) == 2212 || t -> Pid(i,1) == 45 || t -> Pid(i,1) == 321 || t -> Pid(i,1) == -321 || t -> Pid(i,1) == 2112 )) //gamma: 1/22, pi0,+,-: 7/111,8/211,9/-211 (Geant3/pdg)
        {
	  npip += (t -> Pid(i,1)==211)?1:0;
	  npim += (t -> Pid(i,1)==-211)?1:0;
	  vars[0] = 0;//t -> ElecVertTarg();
	  vars[1] = t -> Q2(1);
	  vars[2] = t -> Nu(1);
	  vars[3] = t -> Xb(1);
	  vars[4] = t -> W(1);
	  vars[5] = t -> Sector(ind_first,1);
	  vars[6] = t -> ThetaPQ(i,1);
	  vars[7] = t -> PhiPQ(i,1);
	  Float_t pid = t->Pid(i,1);
	  pdgMutex->Lock();
	  Float_t mass = ( pid == 45 )?1.85756:pdg->GetParticle(pid)->Mass();
	  pdgMutex->UnLock();
	  vars[8] = t -> Zh(i,mass,1);
	  vars[9] = t -> Pt2(i,1);
	  vars[10] = t -> Mx2(i,mass,1);
	  vars[11] = t -> Xf(i,mass,1);
	  vars[12] = t -> T(i,mass,1);
	  vars[13] = t -> Momentum(i,1);
	  vars[14] = 0;//t -> TimeCorr4(0.139570,i);
	  vars[15] = (t -> Z(i,1)) - (t -> Z(ind_first,1));
	  vars[16] = sqrt(t -> Momentum(i,1)*t -> Momentum(i,1) + t->MCMass(i)*t->MCMass(i));//,t->Ein(i)+t->Eout(i));
	  vars[17] = sqrt(t -> Momentum(ind_first,1)*t -> Momentum(ind_first,1) + t->MCMass(ind_first)*t->MCMass(ind_first));//,t->Ein(0)+t->Eout(0));
          vars[18] = t->Momentum(ind_first,1);
          vars[19] = 0;//t->TimeEC(0);
          vars[20] = 0;//t->TimeSC(0);
          vars[21] = 0;//t->PathEC(0);
          vars[22] = 0;//t->PathSC(0);
          vars[23] = event;
          vars[24] = t->Px(i,1);
          vars[25] = t->Py(i,1);
          vars[26] = t->Pz(i,1);
          vars[27] = t->X(ind_first,1);
          vars[28] = t->Y(ind_first,1);
          vars[29] = t->Z(ind_first,1);
          vars[30] = 0;//vert->X(); 
          vars[31] = 0;//vert->Y(); 
          vars[32] = 0;//vert->Z(); 
          vars[33] = 0;//t->TimeEC(i);
          vars[34] = 0;//t->VX_DC(i);//t->XEC(i);
          vars[35] = 0;//t->VY_DC(i);//t->YEC(i);
          vars[36] = 0;//t->VZ_DC(i);//t->ZEC(i);
          vars[37] = t->Px(ind_first,1);
          vars[38] = t->Py(ind_first,1);
          vars[39] = t->Pz(ind_first,1);
          vars[40] = 0;//t->Ein(i);
          vars[41] = 0;//t->Eout(i);
          vars[42] = 0;//t->Ein(0);
          vars[43] = 0;//t->Eout(0);
	  vars[44] = t->Pid(i,1);
	  vars[45] = t->Momentum(i,1)/t->Etot(i,1);
          vars[46] = t->X(i,1);
          vars[47] = t->Y(i,1);
          vars[48] = t->Z(i,1);
	  vars[49] = 0;//t->NpheLTCC(i);
	  vars[50] = 0;//t->NpheHTCC(i);
	  vars[51] = 0;//t->NpheLTCC(0);
	  vars[52] = 0;//t->NpheHTCC(0);
	  vars[53] = 0;//t->Chi2pid(0);
	  vars[54] = 0;//t->Chi2pid(i);
          vars[55] = 0;//t->Epcal(0);
	  vars[56] = 0;//t->Epcal(i);
	  vars[57] = 0;//t->SectorLTCC(0);
	  vars[58] = 0;//t->SectorHTCC(0);
	  vars[59] = 0;//t->SectorECAL(0);
	  vars[60] = 0;//t->SectorLTCC(i);
	  vars[61] = 0;//t->SectorHTCC(i);
	  vars[62] = 0;//t->SectorECAL(i);
	  vars[63] = 0;//t->Helic();

	  vars[64] = 0;//t->LU_PCAL();
	  vars[65] = 0;//t->LV_PCAL();
	  vars[66] = 0;//t->LW_PCAL();
	  vars[67] = 0;//t->LU_ECIN();
	  vars[68] = 0;//t->LV_ECIN();
	  vars[69] = 0;//t->LW_ECIN();
	  vars[70] = 0;//t->LU_ECOUT();
	  vars[71] = 0;//t->LV_ECOUT();
	  vars[72] = 0;//t->LW_ECOUT();

	  vars[73] = 0;//t->LU_PCAL(i);
	  vars[74] = 0;//t->LV_PCAL(i);
	  vars[75] = 0;//t->LW_PCAL(i);
	  vars[76] = 0;//t->LU_ECIN(i);
	  vars[77] = 0;//t->LV_ECIN(i);
	  vars[78] = 0;//t->LW_ECIN(i);
	  vars[79] = 0;//t->LU_ECOUT(i);
	  vars[80] = 0;//t->LV_ECOUT(i);
	  vars[81] = 0;//t->LW_ECOUT(i);

	  vars[82] = 0;//t->HX_PCAL();
	  vars[83] = 0;//t->HY_PCAL();
	  vars[84] = 0;//t->HZ_PCAL();
	  vars[85] = 0;//t->HX_ECIN();
	  vars[86] = 0;//t->HY_ECIN();
	  vars[87] = 0;//t->HZ_ECIN();
	  vars[88] = 0;//t->HX_ECOUT();
	  vars[89] = 0;//t->HY_ECOUT();
	  vars[90] = 0;//t->HZ_ECOUT();
	  //vars[91] = 0;//t->SectorDC(i);  // Warning seg. fault. on sim files, must be fixed!!!!!!!!!1
	  vars[92] = 0;//t->Status(i);
	  vars[93] = 0;//t->Status(0);
	  vars[94] = 0;//t->Px_DC(0);
	  vars[95] = 0;//t->Py_DC(0);
	  vars[96] = 0;//t->Pz_DC(0);
	  vars[97] = 0;//t->Px_DC(i);
	  vars[98] = 0;//t->Py_DC(i);
	  vars[99] = 0;//t->Pz_DC(i);
	  vars[100] = 0;//t->TrajX(i,0);
	  vars[101] = 0;//t->TrajX(i,1);
	  vars[102] = 0;//t->TrajX(i,2);
	  vars[103] = 0;//t->TrajX(i,3);
	  vars[104] = 0;//t->TrajX(i,4);
	  vars[105] = 0;//t->TrajX(i,5);
	  vars[106] = 0;//t->TrajY(i,0);
	  vars[107] = 0;//t->TrajY(i,1);
	  vars[108] = 0;//t->TrajY(i,2);
	  vars[109] = 0;//t->TrajY(i,3);
	  vars[110] = 0;//t->TrajY(i,4);
	  vars[111] = 0;//t->TrajY(i,5);
	  vars[112] = 0;//t->TrajZ(i,0);
	  vars[113] = 0;//t->TrajZ(i,1);
	  vars[114] = 0;//t->TrajZ(i,2);
	  vars[115] = 0;//t->TrajZ(i,3);
	  vars[116] = 0;//t->TrajZ(i,4);
	  vars[117] = 0;//t->TrajZ(i,5);

	  vars[118] = 0;//t->TrajDCX(i,0);
	  vars[119] = 0;//t->TrajDCX(i,1);
	  vars[120] = 0;//t->TrajDCX(i,2);
	  vars[121] = 0;//t->TrajDCY(i,0);
	  vars[122] = 0;//t->TrajDCY(i,1);
	  vars[123] = 0;//t->TrajDCY(i,2);
	  vars[124] = 0;//t->TrajDCZ(i,0);
	  vars[125] = 0;//t->TrajDCZ(i,1);
	  vars[126] = 0;//t->TrajDCZ(i,2);
	  vars[127] = 0;//t->TrajDCX(0,0);
	  vars[128] = 0;//t->TrajDCX(0,1);
	  vars[129] = 0;//t->TrajDCX(0,2);
	  vars[130] = 0;//t->TrajDCY(0,0);
	  vars[131] = 0;//t->TrajDCY(0,1);
	  vars[132] = 0;//t->TrajDCY(0,2);
	  vars[133] = 0;//t->TrajDCZ(0,0);
	  vars[134] = 0;//t->TrajDCZ(0,1);
	  vars[135] = 0;//t->TrajDCZ(0,2);

	  vars[136] = 0;//t->PathTOF(0);
	  vars[137] = 0;//t->TimeTOF(0);
	  vars[138] = 0;//t->PathTOF(i);
	  vars[139] = 0;//t->TimeTOF(i);

	  vars[140] = 0;//t->SectorTOF(0);
	  vars[141] = 0;//t->SectorTOF(i);

	  vars[142] = t->Momentum(ind_first,1)/t->Etot(ind_first,1);
	  vars[143] = 0;//t->STTime();
	  vars[144] = 0;//t->RFTime();

	  vars[145] = 0;//dcx_rot;//region 0
	  vars[146] = 0;//dcy_rot;//region 0
	  vars[147] = 0;//dcx_rot;//region 1
	  vars[148] = 0;//dcy_rot;//region 1
	  vars[149] = 0;//dcx_rot;//region 2
	  vars[150] = 0;//dcy_rot;//region 2
	  vars[151] = 0;//dcx_rot;//region 0
	  vars[152] = 0;//dcy_rot;//region 0
	  vars[153] = 0;//dcx_rot;//region 1
	  vars[154] = 0;//dcy_rot;//region 1
	  vars[155] = 0;//dcx_rot;//region 2
	  vars[156] = 0;//dcy_rot;//region 2

	  vars[157] = t->MCMass(i);//dcy_rot;//region 2
	  vars[158] = t->Event();//dcy_rot;//region 2

	  vars[159] = Npip_rec;
	  vars[160] = Npim_rec;
	  vars[161] = Npip_mc;
	  vars[162] = Npim_mc;

	  vars[163] = rec_elec;
	  vars[164] = 0;
	  vars[165] = 0;
	  vars[166] = 0;
	  vars[167] = 0;

	  vars[168] = 0;
	  vars[169] = 0;
	  vars[170] = 0;
	  vars[171] = 0;
	  vars[172] = 0;
	  vars[173] = 0;
	  vars[174] = 0;
	  vars[175] = 0;
	  vars[176] = 0;
	  vars[177] = 0;
	  vars[178] = 0;
	  vars[179] = 0;
	  vars[180] = 0;
	  vars[181] = 0;
	  vars[182] = 0;
	  vars[183] = 0;
   	  vars[184] = 0;
	  vars[185] = 0;
	  vars[186] = 0;
	  vars[187] = 0;
	  vars[188] = 0;
	  vars[189] = 0;
	  vars[190] = 0;
	  vars[191] = 0;
	  vars[192] = 0;
	  vars[193] = 0;
	  vars[194] = 0;
	  vars[195] = 0;

	  vars[196] = 0;
	  vars[197] = 0;
	  vars[198] = 0;
	  vars[199] = 0;
	  vars[200] = 0;
	  vars[201] = 0;
	  vars[202] = 0;
	  vars[203] = 0;
	  vars[204] = 0;
	  vars[205] = 0;
	  vars[206] = 0;
	  vars[207] = 0;

	  vars[208] = 0;
	  vars[209] = 0;
	  vars[210] = 0;
	  vars[211] = 0;
	  vars[212] = 0;
	  
	  vars[213] = t->Nu()/EBEAM;
	  vars[214] = TMath::ACos(t->Pz(0)/t->Momentum(0))*TMath::RadToDeg();
	  vars[215] = 0;
	  vars[216] = t->DCStatus(i);
	  vars[217] = t->DCNDF(i);
	  vars[218] = t->DCStatus(0);
	  vars[219] = t->DCNDF(0);

	  
	  ntuple_thrown->Fill(vars);

	}
      }
      if (npip != Npip_mc || npim != Npim_mc)
      {
	/*	TThread::Printf("");
	TThread::Printf("WRONG Npid on mc %d",Npip_mc);
	TThread::Printf("%d :: %d ",npim,Npim_mc);*/
      }
    }

    TThread::Lock();
    progress_th[ind] = ++event;
    Nsum++;
    TThread::UnLock();
    if ((event==TH_MAX) || (Nsum==Ntotal)) break;
  }

  fileMutex->Lock();
  output->Write();
  //  objArr->Write("",TObject::kOverwrite);
  output->Close();
  fileMutex->UnLock();
  
  TThread::Printf("%d",__LINE__);

  if(DataElecTh) delete[] DataElecTh;
  delete[] DataElec;
  delete[] vars;
  delete t;
  fileMutex->Lock();
  delete output;
  fileMutex->UnLock();
  
  TThread::Lock();
  NthActive--;
  slotAvailable[ind]=kTRUE;
  TThread::UnLock();
  endsem.Post();
      
  TThread *th =   TThread::Self();
  th->Delete();
  //TThread::Exit();
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
