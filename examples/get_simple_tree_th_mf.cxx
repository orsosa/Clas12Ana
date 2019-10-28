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
#include "data_struct.h"


using namespace std;

float EBEAM=10.6;
TString fname="";
Int_t Nvar, NvarElec;
Int_t Nth = 1;
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
 
  //  cout.width(100);
  coutMutex->Lock();
  TIdentificatorCLAS12 *t = new TIdentificatorCLAS12(fname,EBEAM,true); // March - 19 cooking
  coutMutex->UnLock();
  if (!Ntotal)  Ntotal = t->getNevents();
  if (!QUIET) std::cout<<Ntotal<<std::endl;

  TH_MAX = Ntotal/Nth+1;// max, Nth entries more in last thread, will end reaching Ntotal.
  
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
    
    //th = TThread::GetThread(thn_ind[k]);
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

  TString treeName = "evData";
  fileMutex->Lock();
  if(simul_key == 0) {
    output = new TFile(Form("outfiles/pruned_dataH4_%d.root",Nt_local), "RECREATE", "Data of particles");
  } else { 
    output = new TFile(Form("outfiles/pruned_simulH4_%d.root",Nt_local), "RECREATE", "Data of particles");
  }
  fileMutex->UnLock();

  //TObjArray *objArr = new TObjArray(5);
  //objArr->Add(tElec);
  
  TTree *evTree = new TTree(treeName,"reconstructed particles");
  evTree->SetMaxVirtualSize(100e6); // 100MB max
  
  DATA Evnt;

  initTree(evTree,&Evnt);

  resetDATA(&Evnt);
  
  coutMutex->Lock();
  TIdentificatorCLAS12 *t = new TIdentificatorCLAS12(fname,EBEAM,true); // March - 19 cooking
  coutMutex->UnLock();
  TThread::Printf("###################################### on thread %s, startcnt: %ld, progress_ind: %d",TThread::Self()->GetName(),startcnt,ind);
  Int_t event = 0;
  Int_t npart = 0;
  Float_t revent = 0;
  for (int k=0;k<startcnt;k++) t->Next();
  while (t->Next())
  {
    Int_t nRows = t->GetNRows();
    /*
    if(simul_key == 1) {
      Npip_mc = t->GetNPart(211,1);
      Npim_mc = t->GetNPart(-211,1);
    }
    Npip_rec = t->GetNPart(211);
    Npim_rec = t->GetNPart(-211);
    rec_elec = 0;
    */
    //if(false)
    npart = 0;
    resetDATA(&Evnt);
    revent = t->Event();
    if(nRows>0 && (t->GetCategorization(0)) == "electron")  
    {

      Evnt.Q2 = t->Q2();
      Evnt.W = t->W();
      Evnt.Nu = t->Nu();
      Evnt.Xb = t->Xb();
      vert = t->GetCorrectedVert();
      Float_t vxec=vert->X(); 
      Float_t vyec=vert->Y(); 
      Float_t vzec=vert->Z();
      delete vert;
      Evnt.vxec = vxec;
      Evnt.vyec = vyec;
      Evnt.vzec = vzec;
      Evnt.vxe = t->X(0);
      Evnt.vye = t->Y(0);
      Evnt.vze = t->Z(0);
      Evnt.Pex = t->Px(0);
      Evnt.Pey = t->Py(0);
      Evnt.Pez = t->Pz(0);
      Evnt.event = event; 
      Evnt.Pe = t->Momentum(0);
      Evnt.Ee = t->Etot(0);
      Evnt.e_Ein =  t->Ein(0);
      Evnt.e_Eout = t->Eout(0);
      Evnt.e_Epcal =  t->Epcal(0);
      Evnt.e_npheltcc =  t->NpheLTCC(0);
      Evnt.e_nphehtcc =  t->NpheHTCC(0);
      Evnt.helic =  t->Helic();
      Evnt.e_chi2pid = t->Chi2pid(0);
      Evnt.e_pcal_lu =  t->LU_PCAL();
      Evnt.e_pcal_lv =  t->LV_PCAL();
      Evnt.e_pcal_lw =  t->LW_PCAL();
      Evnt.e_ecin_lu =  t->LU_ECIN();
      Evnt.e_ecin_lv =  t->LV_ECIN();
      Evnt.e_ecin_lw =  t->LW_ECIN();
      Evnt.e_ecout_lu =  t->LU_ECOUT();
      Evnt.e_ecout_lv =  t->LV_ECOUT();
      Evnt.e_ecout_lw =  t->LW_ECOUT();
      Evnt.e_pcal_hx =  t->HX_PCAL();
      Evnt.e_pcal_hy =  t->HY_PCAL();
      Evnt.e_pcal_hz =  t->HZ_PCAL();
      Evnt.e_ecin_hx = t->HX_ECIN();
      Evnt.e_ecin_hy = t->HY_ECIN();
      Evnt.e_ecin_hz = t->HZ_ECIN();
      Evnt.e_ecout_hx = t->HX_ECOUT();
      Evnt.e_ecout_hy = t->HY_ECOUT();
      Evnt.e_ecout_hz = t->HZ_ECOUT();
      Evnt.e_trajx_sl0 = t->TrajDetIdX(0,"DC","DCSL1");
      Evnt.e_trajx_sl1 = t->TrajDetIdX(0,"DC","DCSL2");
      Evnt.e_trajx_sl2 = t->TrajDetIdX(0,"DC","DCSL3");
      Evnt.e_trajx_sl3 = t->TrajDetIdX(0,"DC","DCSL4");
      Evnt.e_trajx_sl4 = t->TrajDetIdX(0,"DC","DCSL5");
      Evnt.e_trajx_sl5 = t->TrajDetIdX(0,"DC","DCSL6");
      Evnt.e_trajy_sl0 = t->TrajDetIdY(0,"DC","DCSL1");
      Evnt.e_trajy_sl1 = t->TrajDetIdY(0,"DC","DCSL2");
      Evnt.e_trajy_sl2 = t->TrajDetIdY(0,"DC","DCSL3");
      Evnt.e_trajy_sl3 = t->TrajDetIdY(0,"DC","DCSL4");
      Evnt.e_trajy_sl4 = t->TrajDetIdY(0,"DC","DCSL5");
      Evnt.e_trajy_sl5 = t->TrajDetIdY(0,"DC","DCSL6");
      Evnt.e_trajz_sl0 = t->TrajDetIdZ(0,"DC","DCSL1");
      Evnt.e_trajz_sl1 = t->TrajDetIdZ(0,"DC","DCSL2");
      Evnt.e_trajz_sl2 = t->TrajDetIdZ(0,"DC","DCSL3");
      Evnt.e_trajz_sl3 = t->TrajDetIdZ(0,"DC","DCSL4");
      Evnt.e_trajz_sl4 = t->TrajDetIdZ(0,"DC","DCSL5");
      Evnt.e_trajz_sl5 = t->TrajDetIdZ(0,"DC","DCSL6");
      Evnt.e_pathtof =    t->PathTOF(0);
      Evnt.e_timetof =    t->TimeTOF(0);
      Evnt.e_sector_tof = t->SectorTOF(0);
      Evnt.e_Beta =       t->Beta(0);
      Evnt.STTime =       t->STTime();
      Evnt.RFTime =       t->RFTime();

      Float_t dcx,dcy,dcx_rot,dcy_rot;

      dcx = t->TrajDetIdX(0,"DC","DCSL1");
      dcy = t->TrajDetIdY(0,"DC","DCSL1");
      rotate_dcxy(dcx,dcy,dcx_rot,dcy_rot);
      Evnt.e_dcx_rot_0 =  dcx_rot;
      Evnt.e_dcy_rot_0 = dcy_rot;

      dcx = t->TrajDetIdX(0,"DC","DCSL3");
      dcy = t->TrajDetIdY(0,"DC","DCSL3");
      rotate_dcxy(dcx,dcy,dcx_rot,dcy_rot);
      Evnt.e_dcx_rot_1 =  dcx_rot;
      Evnt.e_dcy_rot_1 =  dcy_rot;

      dcx = t->TrajDetIdX(0,"DC","DCSL5");
      dcy = t->TrajDetIdY(0,"DC","DCSL5");
      rotate_dcxy(dcx,dcy,dcx_rot,dcy_rot);
      Evnt.e_dcx_rot_2 = dcx_rot;
      Evnt.e_dcy_rot_2 = dcy_rot;
      
      Evnt.e_sector_ltcc =  t->SectorLTCC(0);
      Evnt.e_sector_htcc = t->SectorHTCC(0);
      Evnt.e_sector_ecal = t->SectorECAL(0);
      Evnt.revent = t->Event();
      Evnt.e_dc_chi2 = t->DCChi2(0);
      Evnt.e_ftof1ax = t->TrajFTOF1AX(0);
      Evnt.e_ftof1ay = t->TrajFTOF1AY(0);
      Evnt.e_ftof1az = t->TrajFTOF1AZ(0);
      Evnt.e_pcalx =  t->TrajPCALX(0);
      Evnt.e_pcaly =  t->TrajPCALY(0);
      Evnt.e_pcalz =  t->TrajPCALZ(0);
      Evnt.e_ecalx =  t->TrajECINX(0);
      Evnt.e_ecaly =  t->TrajECINY(0);
      Evnt.e_ecalz =  t->TrajECINZ(0);
      Evnt.e_ltccx =  t->TrajLTCCX(0);
      Evnt.e_ltccy =  t->TrajLTCCY(0);
      Evnt.e_ltccz =  t->TrajLTCCZ(0);
      Evnt.e_htccx =  t->TrajHTCCX(0);
      Evnt.e_htccy =  t->TrajHTCCY(0);
      Evnt.e_htccz =  t->TrajHTCCZ(0);
      Evnt.e_ftof1bx = t->TrajFTOF1BX(0);
      Evnt.e_ftof1by = t->TrajFTOF1BY(0);
      Evnt.e_ftof1bz = t->TrajFTOF1BZ(0);
      Evnt.e_ftof2x =  t->TrajFTOF2X(0);
      Evnt.e_ftof2y =  t->TrajFTOF2Y(0);
      Evnt.e_ftof2z =  t->TrajFTOF2Z(0);
      Evnt.helonline_hel = t->HelicOnline();
      Evnt.helonline_helRaw = t->HelicOnlineRaw();
      Evnt.helflip_hel =  t->HelicFlip();
      Evnt.helflip_helRaw = t->HelicFlipRaw();
      Evnt.helflip_event = t->HelicFlipEvent();
      Evnt.helicRaw = t->HelicRaw();
      Evnt.e_dc_status = t->DCStatus(0);
      Evnt.e_dc_ndf = t->DCNDF(0);
      Evnt.e_sector_dc = t->DCStatus(0);
      Evnt.e_statPart = t->Status(0);
      Evnt.e_DCPx = 0; //t->Px_DC(0);
      Evnt.e_DCPy = 0; //t->Py_DC(0);
      Evnt.e_DCPz = 0; //t->Pz_DC(0);
      Evnt.revent = revent;
      Evnt.y = t->Nu()/EBEAM;
      Evnt.th_e = TMath::ACos(t->Pz(0)/t->Momentum(0))*TMath::RadToDeg();
      Evnt.helicRaw = t->HelicRaw();
          
      //      tElec->Fill(DataElec);
      for (Int_t i = 1; i < nRows; i++) 
      {
	
      	TString category = t->GetCategorization(i);
	if (category == "pi-" || category == "pi+" || category == "gamma" || category == "proton" || category == "deuteron" || category=="K+" || category=="K-" || category=="neutron")
	//	if (category == "pi-" || category == "pi+" || category == "gamma")
	{

	  Evnt.ThetaPQ[npart] = t -> ThetaPQ(i);
	  Evnt.PhiPQ[npart] = t -> PhiPQ(i);
	  Float_t pid = t->Pid(i);
	  pdgMutex->Lock();
	  Float_t mass = ( pid == 45 )?1.85756:pdg->GetParticle(pid)->Mass();
	  pdgMutex->UnLock();
	  Evnt.Zh[npart] = t -> Zh(i,mass);
	  Evnt.Pt2[npart] =  t -> Pt2(i);
	  Evnt.Mx2[npart] =  t -> Mx2(i,mass);
	  Evnt.Xf[npart] = t -> Xf(i,mass);
	  Evnt.T[npart] =  t -> T(i,mass);
	  Evnt.P[npart] =  t -> Momentum(i);
	  Evnt.deltaZ[npart] =  (t -> Z(i)) - (t -> Z(0));
	  Evnt.E[npart] = t->Etot(i);
	  Evnt.Px[npart] =  t->Px(i);
	  Evnt.Py[npart] =  t->Py(i);
	  Evnt.Pz[npart] =  t->Pz(i);
	  Evnt.Ein[npart] =  t->Ein(i);
	  Evnt.Eout[npart] = t->Eout(i);
	  Evnt.pid[npart] =  pid;
	  Evnt.Beta[npart] = t->Beta(i);
	  Evnt.vxh[npart] =  t->X(i);
	  Evnt.vyh[npart] =  t->Y(i);
	  Evnt.vzh[npart] =  t->Z(i);
	  Evnt.npheltcc[npart] =  t->NpheLTCC(i);
	  Evnt.nphehtcc[npart] =  t->NpheHTCC(i);
	  Evnt.chi2pid[npart] =  t->Chi2pid(i);
	  Evnt.Epcal[npart] =  t->Epcal(i);
	  Evnt.sector_ltcc[npart] = t->SectorLTCC(i);
	  Evnt.sector_htcc[npart] = t->SectorHTCC(i);
	  Evnt.sector_ecal[npart] = t->SectorECAL(i);
	  Evnt.pcal_lu[npart] =  t->LU_PCAL(i);
	  Evnt.pcal_lv[npart] =  t->LV_PCAL(i);
	  Evnt.pcal_lw[npart] =  t->LW_PCAL(i);
	  Evnt.ecin_lu[npart] =  t->LU_ECIN(i);
	  Evnt.ecin_lv[npart] =  t->LV_ECIN(i);
	  Evnt.ecin_lw[npart] =  t->LW_ECIN(i);
	  Evnt.ecout_lu[npart] =  t->LU_ECOUT(i);
	  Evnt.ecout_lv[npart] =  t->LV_ECOUT(i);
	  Evnt.ecout_lw[npart] =  t->LW_ECOUT(i);
	  Evnt.sector_dc[npart] = t->SectorDC(i);
	  Evnt.statPart[npart] = t->Status(i);
	  Evnt.DCPx[npart] =  0;//t->Px_DC(i);
	  Evnt.DCPy[npart] =  0;//t->Py_DC(i);
	  Evnt.DCPz[npart] =  0;//t->Pz_DC(i);
	  Evnt.trajx_sl0[npart] = t->TrajDetIdX(i,"DC","DCSL1");
	  Evnt.trajx_sl1[npart] = t->TrajDetIdX(i,"DC","DCSL2");
	  Evnt.trajx_sl2[npart] = t->TrajDetIdX(i,"DC","DCSL3");
	  Evnt.trajx_sl3[npart] = t->TrajDetIdX(i,"DC","DCSL4");
	  Evnt.trajx_sl4[npart] = t->TrajDetIdX(i,"DC","DCSL5");
	  Evnt.trajx_sl5[npart] = t->TrajDetIdX(i,"DC","DCSL6");
	  Evnt.trajy_sl0[npart] = t->TrajDetIdY(i,"DC","DCSL1");
	  Evnt.trajy_sl1[npart] = t->TrajDetIdY(i,"DC","DCSL2");
	  Evnt.trajy_sl2[npart] = t->TrajDetIdY(i,"DC","DCSL3");
	  Evnt.trajy_sl3[npart] = t->TrajDetIdY(i,"DC","DCSL4");
	  Evnt.trajy_sl4[npart] = t->TrajDetIdY(i,"DC","DCSL5");
	  Evnt.trajy_sl5[npart] = t->TrajDetIdY(i,"DC","DCSL6");
	  Evnt.trajz_sl0[npart] = t->TrajDetIdZ(i,"DC","DCSL1");
	  Evnt.trajz_sl1[npart] = t->TrajDetIdZ(i,"DC","DCSL2");
	  Evnt.trajz_sl2[npart] = t->TrajDetIdZ(i,"DC","DCSL3");
	  Evnt.trajz_sl3[npart] = t->TrajDetIdZ(i,"DC","DCSL4");
	  Evnt.trajz_sl4[npart] = t->TrajDetIdZ(i,"DC","DCSL5");
	  Evnt.trajz_sl5[npart] = t->TrajDetIdZ(i,"DC","DCSL6");
	  Evnt.pathtof[npart] = t->PathTOF(i);
	  Evnt.timetof[npart] = t->TimeTOF(i);
	  Evnt.sector_tof[npart] = t->SectorTOF(i);

	  Float_t dcx,dcy,dcx_rot,dcy_rot;

	  dcx   = t->TrajDetIdX(i,"DC","DCSL1");
	  dcy   = t->TrajDetIdY(i,"DC","DCSL1");
	  rotate_dcxy(dcx,dcy,dcx_rot,dcy_rot);
	  Evnt.dcx_rot_0[npart] = dcx_rot;
	  Evnt.dcy_rot_0[npart] = dcy_rot;

	  dcx   = t->TrajDetIdX(i,"DC","DCSL3");
	  dcy   = t->TrajDetIdY(i,"DC","DCSL3");
	  rotate_dcxy(dcx,dcy,dcx_rot,dcy_rot);
	  Evnt.dcx_rot_1[npart] = dcx_rot;
	  Evnt.dcy_rot_1[npart] = dcy_rot;

	  dcx   = t->TrajDetIdX(i,"DC","DCSL5");
	  dcy   = t->TrajDetIdY(i,"DC","DCSL5");
	  rotate_dcxy(dcx,dcy,dcx_rot,dcy_rot);
	  Evnt.dcx_rot_2[npart] = dcx_rot;
	  Evnt.dcy_rot_2[npart] = dcy_rot;

	  Evnt.dc_chi2[npart] = t->DCChi2(i);
	  Evnt.ftof1ax[npart] =  t->TrajFTOF1AX(i);
	  Evnt.ftof1ay[npart] =  t->TrajFTOF1AY(i);
	  Evnt.ftof1az[npart] =  t->TrajFTOF1AZ(i);
	  Evnt.pcalx[npart] = t->TrajPCALX(i);
	  Evnt.pcaly[npart] = t->TrajPCALY(i);
	  Evnt.pcalz[npart] = t->TrajPCALZ(i);
	  Evnt.ecalx[npart] =  t->TrajECINX(i);
	  Evnt.ecaly[npart] =  t->TrajECINY(i);
	  Evnt.ecalz[npart] =  t->TrajECINZ(i);
	  Evnt.ltccx[npart] = t->TrajLTCCX(i);
	  Evnt.ltccy[npart] = t->TrajLTCCY(i);
	  Evnt.ltccz[npart] = t->TrajLTCCZ(i);
	  Evnt.htccx[npart] =  t->TrajHTCCX(i);
	  Evnt.htccy[npart] =  t->TrajHTCCY(i);
	  Evnt.htccz[npart] =  t->TrajHTCCZ(i);
	  Evnt.ftof1bx[npart] = t->TrajFTOF1BX(i);
	  Evnt.ftof1by[npart] = t->TrajFTOF1BY(i);
	  Evnt.ftof1bz[npart] = t->TrajFTOF1BZ(i);
	  Evnt.ftof2x[npart] = t->TrajFTOF2X(i);
	  Evnt.ftof2y[npart] = t->TrajFTOF2Y(i);
	  Evnt.ftof2z[npart] = t->TrajFTOF2Z(i);
	  Evnt.dc_status[npart] =  t->DCStatus(i);
	  Evnt.dc_ndf[npart] = t->DCNDF(i);
	  npart++;
	}
      }
    }
    Evnt.npart = npart;
    npart = 0;
    nRows = t->GetMCNRows();
    Int_t ind_first = 3;
    if(nRows>3 && (simul_key == 1 && t -> Pid(ind_first,1)==11 && t -> LundType(ind_first)==1) ) 
    {
      Evnt.mc_Q2 =  t -> Q2(1);
      Evnt.mc_W = t -> W(1);
      Evnt.mc_Nu = t -> Nu(1);
      Evnt.mc_Xb = t -> Xb(1);
      Evnt.mc_vxe = t -> X(ind_first,1);
      Evnt.mc_vye = t -> Y(ind_first,1);
      Evnt.mc_vze = t -> Z(ind_first,1);
      Evnt.mc_Pex = t -> Px(ind_first,1);
      Evnt.mc_Pey = t -> Py(ind_first,1);
      Evnt.mc_Pez = t -> Pz(ind_first,1);
      Evnt.mc_event = event;
      Evnt.e_mcmass =  t->MCMass(ind_first);
      Evnt.mc_Pe = t->Momentum(ind_first,1);
      Evnt.mc_Ee =  sqrt(t -> Momentum(ind_first,1)*t -> Momentum(ind_first,1) + t->MCMass(ind_first)*t->MCMass(ind_first));
      Evnt.mc_revent = revent;
      Evnt.mc_y = t->Nu()/EBEAM;
      Evnt.mc_th_e =TMath::ACos(t->Pz(ind_first,1)/t->Momentum(ind_first,1))*TMath::RadToDeg();
      Evnt.mc_e_Beta =  t->Momentum(ind_first,1)/t->Etot(ind_first,1);

      for(int i=ind_first + 1; i<nRows; i++) 
      {    
      	if((t -> LundType(i) == 1) && (t -> Pid(i,1)==-211 || t -> Pid(i,1)==211|| t -> Pid(i,1) == 22 || t -> Pid(i,1) == 2212 || t -> Pid(i,1) == 45 || t -> Pid(i,1) == 321 || t -> Pid(i,1) == -321 || t -> Pid(i,1) == 2112 )) //gamma: 1/22, pi0,+,-: 7/111,8/211,9/-211 (Geant3/pdg)
        {
	  Evnt.mc_ThetaPQ[npart] = t -> ThetaPQ(i,1);
	  Evnt.mc_PhiPQ[npart] = t -> PhiPQ(i,1);
	  Float_t pid = t->Pid(i,1);
	  pdgMutex->Lock();
	  Float_t mass = ( pid == 45 )?1.85756:pdg->GetParticle(pid)->Mass();
	  pdgMutex->UnLock();

	  Evnt.mc_Zh[npart] = t -> Zh(i,mass,1);
	  Evnt.mc_Pt2[npart] =  t -> Pt2(i,1);
	  Evnt.mc_Mx2[npart] =  t -> Mx2(i,mass,1);
	  Evnt.mc_Xf[npart] = t -> Xf(i,mass,1);
	  Evnt.mc_T[npart] = t -> T(i,mass,1);
	  Evnt.mc_P[npart] =  t -> Momentum(i,1);
	  Evnt.mc_deltaZ[npart] = (t -> Z(i,1)) - (t -> Z(ind_first,1));
	  Evnt.mc_E[npart] =  sqrt(t -> Momentum(i,1)*t -> Momentum(i,1) + t->MCMass(i)*t->MCMass(i));
	  Evnt.mc_Px[npart] = t->Px(i,1);
	  Evnt.mc_Py[npart] = t->Py(i,1);
	  Evnt.mc_Pz[npart] = t->Py(i,1);
	  Evnt.mc_pid[npart] = pid;
	  Evnt.mc_Beta[npart] =  t->Momentum(i,1)/t->Etot(i,1);
	  Evnt.mc_vxh[npart] = t->X(i,1);
	  Evnt.mc_vyh[npart] = t->Y(i,1);
	  Evnt.mc_vzh[npart] = t->Z(i,1);
	  Evnt.mcmass[npart] =  t->MCMass(i);
	  npart++;

	}
      }
     
    }
    Evnt.mc_npart = npart;
    evTree->Fill();

    
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
