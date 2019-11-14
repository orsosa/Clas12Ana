#include "TSpectrum.h"
#include "Riostream.h"
#include "TApplication.h"
#include "TROOT.h"
#include "TFile.h"
#include "TNtuple.h"
#include "TMath.h"
#include "TBenchmark.h"
#include "TLorentzVector.h"
#include "TCanvas.h"
#include "TH2F.h"
#include "TLatex.h"
#include "TDatabasePDG.h"
#include "TParticlePDG.h"
#include "TMatrixD.h"
#include "TClonesArray.h"
#include "TChain.h"
#include <vector>
#include <iostream>
#include <iomanip>      // std::setw
#include <sys/stat.h>
#include <cstdarg>
#include <algorithm>
#include <map>
#include <string.h>
#include "TRandom3.h"
#include "particle_mix_tree.h"
#include "data_struct.h"
#include "TRegexp.h"
#include "PARTDATA.h"
#include "DETDATA.h"


bool DEBUG=false;
Float_t HELIC=-111;
TRandom3 *rndm;
TString INFILE="data.root", INDIR="",REACTION=":pi+_pi-",OFILE="outfiles/pippim_all.root";
TH1F *hW;
TH1F *hW2;
TH1F *hWmb;
TH1F *hW2mb;
TH1F *hT;
//TH1F *hcfm;
TH2F *hEpi0_th;
Float_t kMpi = TDatabasePDG::Instance()->GetParticle("pi-")->Mass();
Float_t kMpi0 = TDatabasePDG::Instance()->GetParticle("pi0")->Mass();

//Float_t kSpi0=1.94034e-2;//Got from sim.
Float_t kMh =  TDatabasePDG::Instance()->GetParticle("pi-")->Mass();
//Float_t kMPi0=5.39609e-01;
//Float_t kSPi0=5.98542e-02;
bool GSIM=false,EFLAG=false;
Float_t kPt2,kEvent; 
TNtuple *tuple;
//TNtuple *tuple_sim;
TNtuple *tuplemb;
TNtuple *tuplePi0_gamma, *tupleGamma;

DATA Event;
Float_t kEbeam=-1;

long Ne = -1,START = 0;
char st[3]= "C"; // solid target: C Fe Pb
char tt[3] = "C"; // cut on solid target or Deuterium : (st) or D.

Float_t kMprt=0.938272, kMn =0.939565;// kMprt to avoid replace confusion with kMpi
TClonesArray *P4Arr;

class Particle: public TLorentzVector
{
public:
  Float_t vx, vy, vz, pid, time, beta, dcx, dcy, dcz, statPart, dc_chi2, dc_ndf, pcal_lu, pcal_lv, pcal_lw, phiHs, helic002_phiH, helic005_phiH, helic010_phiH, helic020_phiH, m2b;
  inline Double_t P2() const {return P()*P();}
  //TParticlePDG *info;
  Particle() : TLorentzVector(), vx(0),vy(0),vz(0),pid(0),time(0),beta(0),dcx(0),dcy(0),dcz(0),statPart(-1),dc_chi2(-111),dc_ndf(-111),pcal_lu(-111),pcal_lv(-111),pcal_lw(-111), phiHs(-1111), helic002_phiH(-111), helic005_phiH(-111), helic010_phiH(-111), helic020_phiH(-111) {m2b = (P2() - beta*beta*P2())/(beta*beta);}
  Particle(Float_t px,Float_t py, Float_t pz, Float_t e, Float_t x, Float_t y, Float_t z, Float_t pid=0, Float_t t=0, Float_t b =0, Float_t dx=0, Float_t dy=0, Float_t dz=0, Float_t sP=-1, Float_t chi2=-111, Float_t ndf=-111, Float_t plu=-111, Float_t plv =-111, Float_t plw =-111, Float_t phihs = -1111, Float_t h002_phiH = -111, Float_t h005_phiH = -111, Float_t h010_phiH = -111, Float_t h020_phiH = -111): TLorentzVector(px,py,pz,e),vx(x),vy(y),vz(z),pid(pid),time(t),beta(b),dcx(dx),dcy(dy),dcz(dz),statPart(sP),dc_chi2(chi2),dc_ndf(ndf),pcal_lu(plu),pcal_lv(plv),pcal_lw(plw),phiHs(phihs),helic002_phiH(h002_phiH), helic005_phiH(h005_phiH), helic010_phiH(h010_phiH), helic020_phiH(h020_phiH) {m2b = (P2() - beta*beta*P2())/(beta*beta);}
  Particle(TLorentzVector lv, Float_t x=0, Float_t y=0, Float_t z=0, Float_t pid=0, Float_t t=0, Float_t b =0, Float_t dx=0, Float_t dy=0, Float_t dz=0, Float_t sP=-1, Float_t chi2=-111, Float_t ndf=-111, Float_t plu=-111, Float_t plv =-111, Float_t plw =-111, Float_t phihs = -1111, Float_t h002_phiH = -111, Float_t h005_phiH = -111, Float_t h010_phiH = -111, Float_t h020_phiH = -111): TLorentzVector(lv),vx(x),vy(y),vz(z),pid(pid),time(t),beta(b),dcx(dx),dcy(dy),dcz(dz),statPart(sP),dc_chi2(chi2),dc_ndf(ndf),pcal_lu(plu),pcal_lv(plv),pcal_lw(plw),phiHs(phihs),helic002_phiH(h002_phiH), helic005_phiH(h005_phiH), helic010_phiH(h010_phiH), helic020_phiH(h020_phiH){m2b = (P2() - beta*beta*P2())/(beta*beta);}
  Particle(Particle &p):vx(p.vx),vy(p.vy),vz(p.vz),pid(p.pid),time(p.time),beta(p.beta),dcx(p.dcx),dcy(p.dcy),dcz(p.dcz),statPart(p.statPart),dc_chi2(p.dc_chi2),dc_ndf(p.dc_ndf),pcal_lu(p.pcal_lu),pcal_lv(p.pcal_lv),pcal_lw(p.pcal_lw),phiHs(p.phiHs),helic002_phiH(p.helic002_phiH), helic005_phiH(p.helic005_phiH), helic010_phiH(p.helic010_phiH), helic020_phiH(p.helic020_phiH), m2b(p.m2b) {SetVect(p.Vect()); SetT(p.T());}

  inline Particle operator + (const Particle & q) const //const: the object that owns the method will not be modified by this method
  {
    Particle p;
    p.SetVect(Vect()+q.Vect());
    p.SetT(E()+q.T());
    return p;
  }

  inline Bool_t checkFiducial()
  {
    // Converts x,y,z EC hit in CLAS coordinate system
    // into u,v,w distances of the EC hit, and then test the fiducial cut. (Phonetic names)
    Bool_t test=kFALSE;

    Float_t ex=0., wy=0., zd=0., yu=0., ve=0.,  wu=0., xi=0., yi=0., zi=0., ec_phy = 0., phy = 0., rot[3][3];

    // Parameters
    Float_t ec_the = 0.4363323;
    Float_t ylow = -182.974;
    Float_t yhi = 189.956;
    Float_t tgrho = 1.95325; 
    Float_t sinrho = 0.8901256; 
    Float_t cosrho = 0.455715;
    
    // Variables
    ex = vx;
    wy = vy;
    zd = vz;
    
    phy = TMath::ATan2(wy,ex)*57.29578;
    if(phy<0.){phy = phy + 360;}
    phy = phy+30.;
    if(phy>360.){phy = phy-360.;}
    
    ec_phy = ((Int_t) (phy/60.))*1.0471975;
    
    rot[0][0] = TMath::Cos(ec_the)*TMath::Cos(ec_phy);
    rot[0][1] = -TMath::Sin(ec_phy);
    rot[0][2] = TMath::Sin(ec_the)*TMath::Cos(ec_phy);
    rot[1][0] = TMath::Cos(ec_the)*TMath::Sin(ec_phy);
    rot[1][1] = TMath::Cos(ec_phy);
    rot[1][2] = TMath::Sin(ec_the)*TMath::Sin(ec_phy);
    rot[2][0] = -TMath::Sin(ec_the);
    rot[2][1] = 0.;
    rot[2][2] = TMath::Cos(ec_the);
    
    yi = ex*rot[0][0]+wy*rot[1][0]+zd*rot[2][0];
    xi = ex*rot[0][1]+wy*rot[1][1]+zd*rot[2][1];
    zi = ex*rot[0][2]+wy*rot[1][2]+zd*rot[2][2];
    zi = zi-510.32 ;
    
    yu = (yi-ylow)/sinrho;
    ve = (yhi-ylow)/tgrho - xi + (yhi-yi)/tgrho;
    wu = ((yhi-ylow)/tgrho + xi + (yhi-yi)/tgrho)/2./cosrho;

    //U in ]40, 410[ , V in [0,370[ and W in [0,410[.)
    if ((40<yu&&yu<410) && (ve<370)  && (wu<410))
      test=kTRUE;
    //    TVector3 * result3= new TVector3(yu,ve,wu);
  
    return test;
  }

  inline Particle operator += (const Particle & q) 
  {
    SetVect(Vect()+q.Vect());
    SetT(E()+q.T());
    return *this;
  }
  const char *getName()
  {
    return ((TDatabasePDG::Instance())->GetParticle(pid))->GetName(); 
  }
};

class Combo
{
public:
  std::vector<Particle*> kParticles;
  static constexpr Float_t kQ2Tol=0.3,kNuTol=0.3;
  int Npart;
  Float_t lastEvent,kQ2,kNu;
  TLorentzVector *q4; // virtual photon 4th vector.
  Combo(): Npart(0),lastEvent(0),kQ2(0),kNu(0),q4(0){}

  Combo(Combo &c)
  {
    Npart=0;
    for (int k=0;k<c.Npart;k++)
    {
      addParticle(c.kParticles[k]);
    }
    lastEvent=c.lastEvent;
    kQ2=c.kQ2;
    kNu=c.kNu;
  }

  inline  Int_t size(){return kParticles.size();}
  inline  Double_t Px(){return getSum().Px();}
  inline  Double_t Py(){return getSum().Py();}
  inline  Double_t Pz(){return getSum().Pz();}
  inline  Double_t P2(){return getSum().P2();}
  inline  Double_t P(){return getSum().P();}
  inline  Double_t E(){return getSum().E();}
  inline  Double_t M(){return getSum().M();}
  inline  Double_t M2(){return getSum().M2();}



  //  Combo(Particle *&p): {}
  ~Combo()
  {
    //    std::cout<<"Combo destructor called: "<<kParticles.size()<<"\n";
    //for (int k=0;k< kParticles.size();k++)
    // delete kParticles[k];
    delete q4;
    clear();
  }
  void clear(){ kParticles.clear();Npart=0;}

  void boost() //Go to CM frame.
  {
    //Particle *p=&getSum();
    for (int k =0;k<Npart;k++)
      kParticles[k]->Boost(-getSum().BoostVector());
  }

  int addParticle(Particle *p,Float_t ev=0,Bool_t rotfirst=kFALSE, Bool_t fid=kFALSE)
  {
    lastEvent=ev;
    if (Npart==0)
    {
      kQ2= Event.Q2; // using global variable, must be changed!.
      kNu= Event.Nu; // using global variable, must be changed!.
      q4 = new TLorentzVector(-Event.Pex,-Event.Pey,kEbeam-Event.Pez,kEbeam-Event.Ee);
    }
    else if (rotfirst)
    {
      TLorentzVector *q4n = new TLorentzVector(-Event.Pex,-Event.Pey,kEbeam-Event.Pez,kEbeam-Event.Ee);
      // Double_t Dth = q4->Theta() - q4n->Theta(); //no theta rotation, beam direction.
      Double_t Dphi = q4->Phi() - q4n->Phi();
      //p->SetTheta(p->Theta()+Dth);
      p->SetPhi(p->Phi()+Dphi);
    }
    if (fid && !p->checkFiducial()) return Npart;
    // std::cout<<__FILE__<<"::"<<__LINE__<<std::endl;
    kParticles.push_back(p);
    Npart++;

    return Npart;
  }

  Particle getSum()
  {
    //    Particle *p= new Particle();
    Particle p;
      
    for (int k=0;k<Npart;k++)
    {
      p+=*kParticles[k];
    }
    return p;
  }

  Int_t findPid(Int_t pid)
  {
    Int_t ret =0;
    for (int k=0;k<Npart;k++)
      if (pid==kParticles[k]->pid)
	ret++;
    return ret;
  }

  Bool_t isCompatible()
  {
    /////// using global variables, must be changed!
    if ( !( (-kQ2Tol< (kQ2 - Event.Q2) && (kQ2 - Event.Q2) < kQ2Tol) && (-kNuTol< (kNu - Event.Nu) && (kNu - Event.Nu) < kNuTol) ) )
      return kFALSE;

    return kTRUE;
  }


  
  inline Particle* operator [] (const int & i) const
  {
    if (i>=Npart||i<0)
    {
      std::cout << "Index out of bounds" <<std::endl; 
      exit(1);
    }
    return kParticles[i];
  } 


  inline Combo operator + (const Combo & c) const //const: the object that owns the method will not be modified by this method
  {
    
    Combo r;
    for (int k=0;k<c.Npart;k++)
    {
      r.addParticle(c[k]);
    }
    for (int k=0;k<this->Npart;k++)
    {
      r.addParticle((*this)[k]);
    }

    return r;
  }

  inline Combo operator += (const Combo &q) 
  {
    for (int k =0;k<q.Npart;k++)
    {
      addParticle( q[k] );
    }
    return *this;
  }
  
  void print()
  {
    
    for (int k=0;k<Npart;k++)
      std::cout<<kParticles[k]->getName()<<" ";
    std::cout << std::endl;
  }

};


class Reaction
{
public:
  char name[50];
  char filename[50];
  Combo *kPrimary; //primary
  
  
  TParticlePDG *kPdgInfo;
  std::vector<Particle*> *kSecondary; //all secondary

  //std::vector<Particle*> *kCombo; //partial combinations.
  std::vector <Combo *> *kCombo; //partial combinations.
  std::vector <Combo *> kBkgnd; //background combinations.

  static const Int_t kBUFFLIMIT=100;

  std::vector<TH1F*> hSPid;
  TFile *kOutFile;

  DATAMIX mixEvent;
  DATAMIX mixEvent_bkg;

  TClonesArray *DET_CA;
  TClonesArray *PDATA_CA;
  TClonesArray *MC_PDATA_CA;

  
  TTree *kOutData, *kOutBkgnd;

  //  TTree *kP4Tree;
  bool fEMatch;
  
  int kPPid;
  int kNSecondary;
  std::vector<int>::iterator kSIt;
  std::vector<int> kSPid; //Amount of secondary particles for a given pid.
  std::map<int,int> kNSPid;//Amount of secondary pid particles required.
  Reaction(){strcpy(name,"eta -> pi+ pi- a"),strcpy(filename,"test_eta_pippima.root");init();}
  Reaction(const char *n,const char *fn,bool fEMatch=false): fEMatch(fEMatch) {strcpy(name,n); strcpy(filename,fn); init();}

  int store()
  {
    for (int k=0;k<(int)kSPid.size();k++)
    {
      hSPid[k]->Write("",TObject::kOverwrite);
    }
    if (kOutBkgnd) kOutBkgnd->Write("",TObject::kOverwrite);

    std::cout<<"Writing root file: "<<filename<<std::endl;
    return kOutData->Write("",TObject::kOverwrite);
  }
  
  ~Reaction()
  {
    clear();
    delete[] kSecondary;
    delete[] kCombo;
    delete PDATA_CA;
    delete DET_CA;
    kOutFile->Close();
    delete kOutFile;
    //    delete kPdgInfo;
    hSPid.clear();
    std::cout<<"All deleted"<<std::endl;
  }
  
  void clear()
  {
    
    for (int k=0;k<(int)kSPid.size();k++)
    {
      kSecondary[k].clear();
      kCombo[k].clear();
    }
    
  }

  //// To store different information modify init() and fill()
  int init()
  {
    kNSecondary=0;
    kOutFile = new TFile(filename,"recreate");
    //kOutData=new TNtuple("outdata",Form("%s",name),"M:Phx:Phy:Phz:Nu:Q2:Z:Cospq:Pt2:Event:M2_01:M2_02:M_c:Phx_c:Phy_c:Phz_c:Z_c:Cospq_c:Pt2_c:Chi2:qx1:qy1:qz1:qx2:qy2:qz2");
    //kOutData = new TNtuple("outdata",Form("%s",name),

    DET_CA = new TClonesArray("DETDATA");
    PDATA_CA = new TClonesArray("PARTDATA");
    MC_PDATA_CA = new TClonesArray("PARTDATA");
  
    kOutData = new TTree("outdata",Form("%s",name));
    kOutData->SetMaxVirtualSize(100e6); // 100MB max

    kOutBkgnd = 0;//new TTree("outbkgnd",Form("%s",name));
    initMixTree(kOutData,&mixEvent,DET_CA,PDATA_CA,MC_PDATA_CA);

    /*
    ////// DETECTOR DATA
    kOutData->Branch("npart",&detData.npart,"npart/I");
    kOutData->Branch("beta",detData.beta,"beta[npart]/F");
    kOutData->Branch("m2b",detData.m2b,"m2b[npart]/F");
    kOutData->Branch("vx",detData.vx,"vx[npart]/F");
    kOutData->Branch("vy",detData.vy,"vy[npart]/F");
    kOutData->Branch("vz",detData.vz,"vz[npart]/F");
    kOutData->Branch("dcx",detData.dcx,"dcx[npart]/F");
    kOutData->Branch("dcy",detData.dcy,"dcy[npart]/F");
    kOutData->Branch("dcz",detData.dcz,"dcz[npart]/F");
    kOutData->Branch("statPart",detData.statPart,"statPart[npart]/F");
    kOutData->Branch("dc_chi2",detData.dc_chi2,"dc_chi2[npart]/F");
    kOutData->Branch("dc_ndf",detData.dc_ndf,"dc_ndf[npart]/F");
    kOutData->Branch("pcal_lu",detData.pcal_lu,"pcal_lu[npart]/F");
    kOutData->Branch("pcal_lv",detData.pcal_lv,"pcal_lv[npart]/F");
    kOutData->Branch("pcal_lw",detData.pcal_lw,"pcal_lw[npart]/F");
    ////////////  
    kOutBkgnd = kOutData->CloneTree(0);
    kOutBkgnd->SetName("outbkgnd");
    */
    
    return 0;
  }

  int fillTree()
  {
    kOutData->Fill();
    PDATA_CA->Clear();
    DET_CA->Clear();
    return 0;
  }
 
  int setVars(int k,bool isMC = false)
  {
    setVars(kPrimary, &mixEvent, k, isMC);
    return 0;
  }

  int setElecVars(bool isMC = false){
    
    if (isMC){
      mixEvent.mc_Q2 = Event.mc_Q2;
      mixEvent.mc_W = Event.mc_W;
      mixEvent.mc_Nu = Event.mc_Nu;
      mixEvent.mc_Xb = Event.mc_Xb;
      mixEvent.mc_vxe = Event.mc_vxe;
      mixEvent.mc_vye = Event.mc_vye;
      mixEvent.mc_vze = Event.mc_vze;
      mixEvent.mc_Pex = Event.mc_Pex;
      mixEvent.mc_Pey = Event.mc_Pey;
      mixEvent.mc_Pez = Event.mc_Pez;
      mixEvent.mc_event = Event.mc_event;
      mixEvent.e_mcmass = Event.e_mcmass;
      mixEvent.mc_Pe = Event.mc_Pe;
      mixEvent.mc_Ee = Event.mc_Ee;
      mixEvent.mc_revent = Event.mc_revent;
      mixEvent.mc_y = Event.mc_y;
      mixEvent.mc_th_e = Event.mc_th_e;
      mixEvent.mc_phi_e = Event.mc_phi_e;
      mixEvent.mc_e_Beta = Event.mc_e_Beta;
      mixEvent.mc_helic = 0;
      Float_t y = Event.mc_y;
      Float_t g =  2*kMprt*Event.mc_Xb/sqrt(Event.mc_Q2);
      Float_t epsilon =  (1 - y - 1/4.*g*g*y*y )/ (1 - y + 1./2*y*y + 1./4*g*g*y*y);
      mixEvent.mc_gamm = g;
      mixEvent.mc_epsilon = epsilon;
      Float_t A = y*y/(2*(1-epsilon));
      mixEvent.mc_fA = A;
      mixEvent.mc_fB = A*epsilon;
      mixEvent.mc_fC = A*sqrt(1-epsilon*epsilon) ;
      mixEvent.mc_fV = A*sqrt(2*epsilon*(1+epsilon)); 
      mixEvent.mc_fW = A*sqrt(2*epsilon*(1-epsilon));
    }
    else {
      mixEvent.Q2 = Event.Q2;
      mixEvent.W = Event.W;
      mixEvent.Nu = Event.Nu;
      mixEvent.Xb = Event.Xb;
      mixEvent.vxec = Event.vxec;
      mixEvent.vyec = Event.vyec;
      mixEvent.vzec = Event.vzec;
      mixEvent.vxe = Event.vxe;
      mixEvent.vye = Event.vye;
      mixEvent.vze = Event.vze;
      mixEvent.Pex = Event.Pex;
      mixEvent.Pey = Event.Pey;
      mixEvent.Pez = Event.Pez;
      mixEvent.event = Event.event;
      mixEvent.Pe = Event.Pe;
      mixEvent.Ee = Event.Ee;
      mixEvent.e_Ein = Event.e_Ein;
      mixEvent.e_Eout = Event.e_Eout;
      mixEvent.e_Epcal = Event.e_Epcal;
      mixEvent.e_npheltcc = Event.e_npheltcc;
      mixEvent.e_nphehtcc = Event.e_nphehtcc;
      mixEvent.helic = Event.helic;
      mixEvent.e_chi2pid = Event.e_chi2pid;
      mixEvent.e_pcal_lu = Event.e_pcal_lu;
      mixEvent.e_pcal_lv = Event.e_pcal_lv;
      mixEvent.e_pcal_lw = Event.e_pcal_lw;
      mixEvent.e_ecin_lu = Event.e_ecin_lu;
      mixEvent.e_ecin_lv = Event.e_ecin_lv;
      mixEvent.e_ecin_lw = Event.e_ecin_lw;
      mixEvent.e_ecout_lu = Event.e_ecout_lu;
      mixEvent.e_ecout_lv = Event.e_ecout_lv;
      mixEvent.e_ecout_lw = Event.e_ecout_lw;
      mixEvent.e_pcal_hx = Event.e_pcal_hx;
      mixEvent.e_pcal_hy = Event.e_pcal_hy;
      mixEvent.e_pcal_hz = Event.e_pcal_hz;
      mixEvent.e_ecin_hx = Event.e_ecin_hx;
      mixEvent.e_ecin_hy = Event.e_ecin_hy;
      mixEvent.e_ecin_hz = Event.e_ecin_hz;
      mixEvent.e_ecout_hx = Event.e_ecout_hx;
      mixEvent.e_ecout_hy = Event.e_ecout_hy;
      mixEvent.e_ecout_hz = Event.e_ecout_hz;
      mixEvent.e_trajx_sl0 = Event.e_trajx_sl0;
      mixEvent.e_trajx_sl1 = Event.e_trajx_sl1;
      mixEvent.e_trajx_sl2 = Event.e_trajx_sl2;
      mixEvent.e_trajx_sl3 = Event.e_trajx_sl3;
      mixEvent.e_trajx_sl4 = Event.e_trajx_sl4;
      mixEvent.e_trajx_sl5 = Event.e_trajx_sl5;
      mixEvent.e_trajy_sl0 = Event.e_trajy_sl0;
      mixEvent.e_trajy_sl1 = Event.e_trajy_sl1;
      mixEvent.e_trajy_sl2 = Event.e_trajy_sl2;
      mixEvent.e_trajy_sl3 = Event.e_trajy_sl3;
      mixEvent.e_trajy_sl4 = Event.e_trajy_sl4;
      mixEvent.e_trajy_sl5 = Event.e_trajy_sl5;
      mixEvent.e_trajz_sl0 = Event.e_trajz_sl0;
      mixEvent.e_trajz_sl1 = Event.e_trajz_sl1;
      mixEvent.e_trajz_sl2 = Event.e_trajz_sl2;
      mixEvent.e_trajz_sl3 = Event.e_trajz_sl3;
      mixEvent.e_trajz_sl4 = Event.e_trajz_sl4;
      mixEvent.e_trajz_sl5 = Event.e_trajz_sl5;
      mixEvent.e_pathtof = Event.e_pathtof;
      mixEvent.e_timetof = Event.e_timetof;
      mixEvent.e_sector_tof = Event.e_sector_tof;
      mixEvent.e_Beta = Event.e_Beta;
      mixEvent.STTime = Event.STTime;
      mixEvent.RFTime = Event.RFTime;
      mixEvent.e_dcx_rot_0 = Event.e_dcx_rot_0;
      mixEvent.e_dcy_rot_0 = Event.e_dcy_rot_0;
      mixEvent.e_dcx_rot_1 = Event.e_dcx_rot_1;
      mixEvent.e_dcy_rot_1 = Event.e_dcy_rot_1;
      mixEvent.e_dcx_rot_2 = Event.e_dcx_rot_2;
      mixEvent.e_dcy_rot_2 = Event.e_dcy_rot_2;
      mixEvent.e_sector_ltcc = Event.e_sector_ltcc;
      mixEvent.e_sector_htcc = Event.e_sector_htcc;
      mixEvent.e_sector_ecal = Event.e_sector_ecal;
      mixEvent.e_dc_chi2 = Event.e_dc_chi2;
      mixEvent.e_ftof1ax = Event.e_ftof1ax;
      mixEvent.e_ftof1ay = Event.e_ftof1ay;
      mixEvent.e_ftof1az = Event.e_ftof1az;
      mixEvent.e_pcalx = Event.e_pcalx;
      mixEvent.e_pcaly = Event.e_pcaly;
      mixEvent.e_pcalz = Event.e_pcalz;
      mixEvent.e_ecalx = Event.e_ecalx;
      mixEvent.e_ecaly = Event.e_ecaly;
      mixEvent.e_ecalz = Event.e_ecalz;
      mixEvent.e_ltccx = Event.e_ltccx;
      mixEvent.e_ltccy = Event.e_ltccy;
      mixEvent.e_ltccz = Event.e_ltccz;
      mixEvent.e_htccx = Event.e_htccx;
      mixEvent.e_htccy = Event.e_htccy;
      mixEvent.e_htccz = Event.e_htccz;
      mixEvent.e_ftof1bx = Event.e_ftof1bx;
      mixEvent.e_ftof1by = Event.e_ftof1by;
      mixEvent.e_ftof1bz = Event.e_ftof1bz;
      mixEvent.e_ftof2x = Event.e_ftof2x;
      mixEvent.e_ftof2y = Event.e_ftof2y;
      mixEvent.e_ftof2z = Event.e_ftof2z;
      mixEvent.helonline_hel = Event.helonline_hel;
      mixEvent.helonline_helRaw = Event.helonline_helRaw;
      mixEvent.helflip_hel = Event.helflip_hel;
      mixEvent.helflip_helRaw = Event.helflip_helRaw;
      mixEvent.helflip_event = Event.helflip_event;
      mixEvent.e_dc_status = Event.e_dc_status;
      mixEvent.e_dc_ndf = Event.e_dc_ndf;
      mixEvent.e_sector_dc = Event.e_sector_dc;
      mixEvent.e_statPart = Event.e_statPart;
      mixEvent.e_DCPx = Event.e_DCPx;
      mixEvent.e_DCPy = Event.e_DCPy;
      mixEvent.e_DCPz = Event.e_DCPz;
      mixEvent.revent = Event.revent;
      mixEvent.y = Event.y;
      mixEvent.th_e = Event.th_e;
      mixEvent.phi_e = Event.phi_e;
      mixEvent.helicRaw = Event.helicRaw;
      Float_t y = Event.y;
      Float_t g =  2*kMprt*Event.Xb/sqrt(Event.Q2);
      Float_t epsilon =  (1 - y - 1/4.*g*g*y*y )/ (1 - y + 1./2*y*y + 1./4*g*g*y*y);
      mixEvent.gamm = g;
      mixEvent.epsilon = epsilon;
      Float_t A = y*y/(2*(1-epsilon));
      mixEvent.fA = A;
      mixEvent.fB = A*epsilon;
      mixEvent.fC = A*sqrt(1-epsilon*epsilon) ;
      mixEvent.fV = A*sqrt(2*epsilon*(1+epsilon)); 
      mixEvent.fW = A*sqrt(2*epsilon*(1-epsilon));
    }
    return 0;
  }

  int setVars(Combo *comb,DATAMIX *evnt,int k,bool isMC = false){
    if (DEBUG ) std::cout<<"##### combo size: "<<comb->size()<<" #######"<<std::endl;
    Float_t Pex,Pez,Pey,Q2,Nu,W,event;
    if (!isMC){
      
      Pex = Event.Pex;
      Pey = Event.Pey;
      Pez = Event.Pez;
      Nu  = Event.Nu;
      Q2  = Event.Q2;
      W   = Event.W;
      event = Event.event;
    }
    else{
      Pex = Event.mc_Pex;
      Pey = Event.mc_Pey;
      Pez = Event.mc_Pez;
      Nu  = Event.mc_Nu;
      Q2  = Event.mc_Q2;
      W   = Event.mc_W;
      event = Event.mc_event;
    }

    Double_t Px = comb->Px();
    Double_t Py = comb->Py();
    Double_t Pz = comb->Pz();
    Double_t E = comb->E();
    Double_t P2 = comb->P2();
    Double_t M2 = comb->M2();
    Double_t M =  (M2>=0)?TMath::Sqrt(M2):-1.0;
    //Float_t theta=Pz/sqrt(P2);
    Float_t theta_0 = (*comb)[0]->Theta();
    Float_t theta_1 = (*comb)[1]->Theta();

    Float_t cospq = ((kEbeam-Pez)*Pz - Pex*Px - Pey*Py)/( sqrt((Q2 + Nu*Nu)*P2) );

    Float_t cospq0 = ((kEbeam-Pez)*(*comb)[0]->Pz() - Pex*(*comb)[0]->Px() - Pey*(*comb)[0]->Py() )/( sqrt((Q2 + Nu*Nu))*(*comb)[0]->P());

    Float_t cospq1 = ((kEbeam-Pez)*(*comb)[1]->Pz() - Pex*(*comb)[1]->Px() - Pey*(*comb)[1]->Py())/( sqrt((Q2 + Nu*Nu))*(*comb)[1]->P());
 
    Float_t Pt2 = P2*(1-cospq*cospq);
    Float_t Pl2 = P2*cospq*cospq;

    Float_t Pl = sqrt(P2)*cospq;
    Float_t Pl0 = (*comb)[0]->P()*cospq0;
    Float_t Pl1 = (*comb)[1]->P()*cospq1;

    Double_t phi_pq;
    TVector3 Vhad(Px,Py,Pz);
    TVector3 Vvirt(-Pex,-Pey,kEbeam-Pez);
    Double_t phi_z = TMath::Pi()-Vvirt.Phi();
    Vvirt.RotateZ(phi_z);
    Vhad.RotateZ(phi_z);
    TVector3 Vhelp(0.,0.,1.);
    Double_t phi_y = Vvirt.Angle(Vhelp);
    Vvirt.RotateY(phi_y);
    Vhad.RotateY(phi_y);
    phi_pq=Vhad.Phi() * 180./(TMath::Pi());

    TLorentzVector q_lv(-Pex,-Pey,kEbeam-Pez,Nu); // virtual photon 4vec
    TLorentzVector P_lv(0,0,0,kMprt); // Nucleon 4vec
    TLorentzVector Ptot_lv = q_lv+P_lv; // total 4vec

    TLorentzVector Ph(Px,Py,Pz,E);
    TLorentzVector P0((*comb)[0]->Px(),(*comb)[0]->Py(),(*comb)[0]->Pz(),(*comb)[0]->E());
    TLorentzVector P1((*comb)[1]->Px(),(*comb)[1]->Py(),(*comb)[1]->Pz(),(*comb)[1]->E());
    TLorentzVector k_in(0,0,kEbeam,kEbeam);

    TLorentzVector P0_dicm = P0;
    TLorentzVector P1_dicm = P1;
    
   /// Boost to di-hadron cm frame///
    P0_dicm.Boost(-Ph.BoostVector());
    P1_dicm.Boost(-Ph.BoostVector());
    /////////////////////////////////

    /*
      TLorentzVector P0_BF(0.1939,0.0749,0.3729,0.8912),
      P1_BF(0.4010,0.1279,0.4797,0.6522),
      q_BF(0.6616,0.8348,2.5974,2.5268),
      P_BF(0,0,0,kMprt),
      qin_BF(0,0,q_BF.M(),0);
      TLorentzVector PX_BF = q_BF + P_BF - P0_BF - P1_BF;


      TLorentzVector p_BF_boost(q_BF-qin_BF );
      q_BF.Boost ( -p_BF_boost.BoostVector() );
      q_BF.Print();

    */

    //// pseudorapidity in Breit frame (light-front coor.)
    TLorentzVector P0_BF = P0;
    TLorentzVector P1_BF = P1;
    TLorentzVector q_BF = q_lv;
    TLorentzVector P_BF = P_lv;
    TLorentzVector PX_BF = q_BF + P_BF - P0_BF - P1_BF;
    
    Float_t xbj = Q2/2/Nu/kMprt;
    //Nachtmann variable
    Float_t xn = 2*xbj/(1 + sqrt(1 + 4*xbj*xbj*kMprt*kMprt/Q2));
    Float_t Q = sqrt(Q2);

    Float_t z0 = (*comb)[0]->E()/Nu;
    Float_t z1 = (*comb)[1]->E()/Nu;
    Float_t ztot_inv = 1./(z0+z1);
    
    TVector3 p0Lv_BF = (P0_BF.Vect()*q_BF.Vect().Unit())*q_BF.Vect().Unit();
    TVector3 p1Lv_BF = (P1_BF.Vect()*q_BF.Vect().Unit())*q_BF.Vect().Unit();
    TVector3 p0Tv_BF = P0_BF.Vect() - p0Lv_BF;
    TVector3 p1Tv_BF = P1_BF.Vect() - p1Lv_BF;

    // pip /////    
    Float_t MhT20 = p0Tv_BF.Mag2() + kMh*kMh;
    Float_t etaBF_0plus = Q*z0*(Q2 - xn*xn*kMprt*kMprt)/
      (2*xn*xn*kMprt*kMprt*sqrt(MhT20));

    Float_t etaBF_0minus = Q*z0*(Q2 - xn*xn*kMprt*kMprt)/
      (2*xn*xn*kMprt*kMprt*sqrt(MhT20));

    etaBF_0plus += Q/(kMprt*xn)*sqrt( z0*z0*TMath::Power((Q2 - xn*xn*kMprt*kMprt),2)/(4*xn*xn*kMprt*kMprt*MhT20) -1);

    etaBF_0minus -= Q/(kMprt*xn)*sqrt( z0*z0*TMath::Power((Q2 - xn*xn*kMprt*kMprt),2)/(4*xn*xn*kMprt*kMprt*MhT20) -1); 

    etaBF_0plus = TMath::Log(etaBF_0plus);
    etaBF_0minus = TMath::Log(etaBF_0minus);

    //////// end pip ////////

    // pim /////    

    Float_t MhT21 = p1Tv_BF.Mag2() + kMh*kMh;
    Float_t etaBF_1plus = Q*z1*(Q2 - xn*xn*kMprt*kMprt)/
      (2*xn*xn*kMprt*kMprt*sqrt(MhT21));

    Float_t etaBF_1minus = Q*z1*(Q2 - xn*xn*kMprt*kMprt)/
      (2*xn*xn*kMprt*kMprt*sqrt(MhT21));

    etaBF_1plus += Q/(xn*kMprt)*sqrt( z1*z1*TMath::Power((Q2 - xn*xn*kMprt*kMprt),2)/(4*xn*xn*kMprt*kMprt*MhT21) -1);

    etaBF_1minus -= Q/(xn*kMprt)*sqrt( z1*z1*TMath::Power((Q2 - xn*xn*kMprt*kMprt),2)/(4*xn*xn*kMprt*kMprt*MhT21) -1); 

    etaBF_1plus = log(etaBF_1plus);
    etaBF_1minus = log(etaBF_1minus);

    //////// end pim ////////
    
    Float_t cos_theta_P0cm = P0_dicm.Vect()*Ph.Vect()/P0_dicm.Vect().Mag()/Ph.Vect().Mag();

    Float_t sin_theta_P0cm = P0_dicm.Vect().Cross(Ph.Vect()).Mag()/P0_dicm.Vect().Mag()/Ph.Vect().Mag();

    // -Ptot_lv.BoostVector() // to the a-P c.m.
    // -Ph.BoostVector() // to the di hadron c.m.
    
    //// Boost to cm virtual photon - nucleon.
    q_lv.Boost(-Ptot_lv.BoostVector());
    P_lv.Boost(-Ptot_lv.BoostVector());
    Ph.Boost(-Ptot_lv.BoostVector());
    P0.Boost(-Ptot_lv.BoostVector());
    P1.Boost(-Ptot_lv.BoostVector());
    k_in.Boost(-Ptot_lv.BoostVector());

    Float_t Plt = Ph.Vect()*q_lv.Vect().Unit();
    Float_t Pl0m = P0.Vect()*q_lv.Vect().Unit();
    Float_t Pl1m = P1.Vect()*q_lv.Vect().Unit();

    //// rapidity in CM
    Float_t etaCM_0 = 0.5*TMath::Log( (P0.E() + Pl0m) / (P0.E() - Pl0m) );
    Float_t etaCM_1 = 0.5*TMath::Log( (P1.E() + Pl1m) / (P1.E() - Pl1m) );

    //// rapidity BF harut
    //    Float_t g2    = 4*xbj*xbj*kMprt*kMprt/Q2;
    //Float_t g24   = xn*xn*kMprt*kMprt/Q2;
    Float_t yhc1  = xn*xn*kMprt*kMprt+xn*Q2;
    Float_t yhc2  = (1.0-xn)*Q2;
    Float_t yhc   = TMath::Log(sqrt(yhc1/yhc2));
    Float_t etaBF0 = -etaCM_0 - yhc;
    Float_t etaBF1 = -etaCM_1 - yhc;
    
    //// phiR //////////////////
    TVector3 Ph_u = Ph.Vect().Unit();
    TVector3 R = P0.Vect() - P1.Vect();
    R=R*0.5;
    TVector3 RT = R-(R*Ph_u)*Ph_u;

    Float_t qxkRT_sign = q_lv.Vect().Cross(k_in.Vect())*RT;
    qxkRT_sign /= TMath::Abs(qxkRT_sign);

    // Float_t qxkST_sign = q_lv.Vect().Cross(k_in)*ST;
    // qxkST_sign /= TMath::Abs(qxkST_sign);

    Float_t qxkqxRT = (q_lv.Vect().Cross(k_in.Vect()))*(q_lv.Vect().Cross(RT));
    Float_t qxkqxRT_max = (q_lv.Vect().Cross(k_in.Vect())).Mag()*(q_lv.Vect().Cross(RT)).Mag();
    
    Float_t phiR =  qxkRT_sign*TMath::ACos(qxkqxRT/qxkqxRT_max)*TMath::RadToDeg();
    phiR=phiR<0?phiR+360:phiR;
    ////////////////

    TVector3 phTv = Ph.Vect() - (Ph.Vect()*q_lv.Vect().Unit())*q_lv.Vect().Unit();
    //// phiR_cov //////////////////
    TVector3 p0Lv = (P0.Vect()*q_lv.Vect().Unit())*q_lv.Vect().Unit();
    TVector3 p1Lv = (P1.Vect()*q_lv.Vect().Unit())*q_lv.Vect().Unit();
    TVector3 p0Tv = P0.Vect() - p0Lv;
    TVector3 p1Tv = P1.Vect() - p1Lv; 

    TVector3 RT_cov = (z1*p0Tv - z0*p1Tv)*ztot_inv;

    Float_t qxkRT_cov_sign = q_lv.Vect().Cross(k_in.Vect())*RT_cov;
    qxkRT_cov_sign /= TMath::Abs(qxkRT_cov_sign);

    Float_t qxkqxRT_cov = (q_lv.Vect().Cross(k_in.Vect()))*(q_lv.Vect().Cross(RT_cov));
    Float_t qxkqxRT_cov_max = (q_lv.Vect().Cross(k_in.Vect())).Mag()*(q_lv.Vect().Cross(RT_cov)).Mag();
    
    Float_t phiR_cov =  qxkRT_cov_sign*TMath::ACos(qxkqxRT_cov/qxkqxRT_cov_max)*TMath::RadToDeg();
    phiR_cov=phiR_cov<0?phiR_cov+360:phiR_cov;
    ////////////////

    /////// phiR_ha transverse to q ///

    TVector3 RTha = 0.5*(p0Tv - p1Tv);

    Float_t qxkRTha_sign = q_lv.Vect().Cross(k_in.Vect())*RTha;
    qxkRTha_sign /= TMath::Abs(qxkRTha_sign);

    // Float_t qxkST_sign = q_lv.Vect().Cross(k_in)*ST;
    // qxkST_sign /= TMath::Abs(qxkST_sign);

    Float_t qxkqxRTha = (q_lv.Vect().Cross(k_in.Vect()))*(q_lv.Vect().Cross(RTha));
    Float_t qxkqxRTha_max = (q_lv.Vect().Cross(k_in.Vect())).Mag()*(q_lv.Vect().Cross(RTha)).Mag();
    
    Float_t phiR_ha =  qxkRTha_sign*TMath::ACos(qxkqxRTha/qxkqxRTha_max)*TMath::RadToDeg();
    phiR_ha=phiR_ha<0?phiR_ha+360:phiR_ha;
    ////////////////

    
    //// phiR_covH Transverse to Hadron //////////////////
    TVector3 p0LvH = (P0.Vect()*Ph_u)*Ph_u;
    TVector3 p1LvH = (P1.Vect()*Ph_u)*Ph_u;
    TVector3 p0TvH = P0.Vect() - p0LvH;
    TVector3 p1TvH = P1.Vect() - p1LvH; 

    TVector3 RTH_cov = (z1*p0TvH - z0*p1TvH)*ztot_inv;

    Float_t qxkRTH_cov_sign = q_lv.Vect().Cross(k_in.Vect())*RTH_cov;
    qxkRTH_cov_sign /= TMath::Abs(qxkRTH_cov_sign);

    Float_t qxkqxRTH_cov = (q_lv.Vect().Cross(k_in.Vect()))*(q_lv.Vect().Cross(RTH_cov));
    Float_t qxkqxRTH_cov_max = (q_lv.Vect().Cross(k_in.Vect())).Mag()*(q_lv.Vect().Cross(RTH_cov)).Mag();
    
    Float_t phiR_covH =  qxkRTH_cov_sign*TMath::ACos(qxkqxRTH_cov/qxkqxRTH_cov_max)*TMath::RadToDeg();
    phiR_covH=phiR_covH<0?phiR_covH+360:phiR_covH;
    ////////////////

    //// phiH //////////////////
    TVector3 Phv = Ph.Vect();
    Float_t qxkPhv_sign = q_lv.Vect().Cross(k_in.Vect())*Phv;
    qxkPhv_sign /= TMath::Abs(qxkPhv_sign);

    Float_t qxkqxPhv = (q_lv.Vect().Cross(k_in.Vect()))*(q_lv.Vect().Cross(Phv));
    Float_t qxkqxPhv_max = (q_lv.Vect().Cross(k_in.Vect())).Mag()*(q_lv.Vect().Cross(Phv)).Mag();
    
    Float_t phiH =  qxkPhv_sign*TMath::ACos(qxkqxPhv/qxkqxPhv_max)*TMath::RadToDeg();
    phiH=phiH<0?phiH+360:phiH;
    ////////////////
    //// phiH0 /// TRENTO ////
    TVector3 P0v = P0.Vect();
    Float_t qxkP0v_sign = q_lv.Vect().Cross(k_in.Vect())*P0v;
    qxkP0v_sign /= TMath::Abs(qxkP0v_sign);
    
    Float_t qxkqxP0v = (q_lv.Vect().Cross(k_in.Vect()))*(q_lv.Vect().Cross(P0v));
    Float_t qxkqxP0v_max = (q_lv.Vect().Cross(k_in.Vect())).Mag()*(q_lv.Vect().Cross(P0v)).Mag();
    Float_t phiH0 =  qxkP0v_sign*TMath::ACos(qxkqxP0v/qxkqxP0v_max)*TMath::RadToDeg();
    phiH0=phiH0<0?phiH0+360:phiH0;
    //////////////////

    //// phiH1 /// TRENTO ////
    TVector3 P1v = P1.Vect();
    Float_t qxkP1v_sign = q_lv.Vect().Cross(k_in.Vect())*P1v;
    qxkP1v_sign /= TMath::Abs(qxkP1v_sign);
    
    Float_t qxkqxP1v = (q_lv.Vect().Cross(k_in.Vect()))*(q_lv.Vect().Cross(P1v));
    Float_t qxkqxP1v_max = (q_lv.Vect().Cross(k_in.Vect())).Mag()*(q_lv.Vect().Cross(P1v)).Mag();
    Float_t phiH1 =  qxkP1v_sign*TMath::ACos(qxkqxP1v/qxkqxP1v_max)*TMath::RadToDeg();
    phiH1=phiH1<0?phiH1+360:phiH1;
    //////////////////
  
    Float_t Mx2 = W*W + M*M - 2*( (Nu+kMprt)*E - sqrt((Q2 + Nu*Nu)*Pl2));

    Float_t Phmax = TMath::Sqrt( TMath::Power( W*W + M*M - kMprt*kMprt, 2 ) - 4*W*W*M*M )/(2* W);

    Float_t Ph0max = TMath::Sqrt( TMath::Power( W*W + kMpi*kMpi - kMprt*kMprt - 2*W*(*comb)[1]->E(), 2 ) - 4*W*W*kMpi*kMpi )/(2* W);
    Float_t Ph1max = TMath::Sqrt( TMath::Power( W*W + kMpi*kMpi - kMprt*kMprt - 2*W*(*comb)[0]->E(), 2 ) - 4*W*W*kMpi*kMpi )/(2* W);

    ////// LORENTZ BOOST //////////
    Float_t b=TMath::Sqrt(Q2 + Nu*Nu)/(Nu+kMprt);
    Float_t g=(Nu+kMprt)/W;

    Float_t PlCM = g*(Pl - b*E);
    Float_t Pl0CM = g*( Pl0 - b*(*comb)[0]->E() );
    Float_t Pl1CM = g*( Pl1 - b*(*comb)[1]->E() );
      
    Float_t xF = PlCM/Phmax;
    Float_t xF0 = Pl0CM/Ph0max;
    Float_t xF1 = Pl1CM/Ph1max;

    Float_t xFm = 2*PlCM/W;
    Float_t xFm0 = 2*Pl0CM/W;
    Float_t xFm1 = 2*Pl1CM/W;

    Float_t xFo = 2*Plt/W;
    Float_t xF0o = 2*Pl0m/W;
    Float_t xF1o = 2*Pl1m/W;

    /*** set MC data ***/
    if (isMC){
      evnt->mc_M[k] = M;
      evnt->mc_Phx[k] = Px;
      evnt->mc_Phy[k] = Py;
      evnt->mc_Phz[k] = Pz;
      evnt->mc_Z[k] = E/Nu;
      evnt->mc_Cospq[k] = cospq;
      evnt->mc_Pt2[k] = Pt2;
      evnt->mc_Event[k] = event;
      evnt->mc_M2_01[k] = ((comb->Npart==3)? ( *(*comb)[0] + *(*comb)[1]).M2() : 0);
      evnt->mc_M2_02[k] = ((comb->Npart==3)? ( *(*comb)[0] + *(*comb)[2]).M2() : 0);
      evnt->mc_phiH[k] = phiH;
      evnt->mc_phiR[k] = phiR;

      evnt->mc_helic002_phiR[k] = get_helicity(phiR,0.02);
      evnt->mc_helic005_phiR[k] = get_helicity(phiR,0.05);
      evnt->mc_helic010_phiR[k] = get_helicity(phiR,0.1);
      evnt->mc_helic020_phiR[k] = get_helicity(phiR,0.2);

      evnt->mc_helic002_phiRst[k] = get_helicity(phiR,0.02,sin_theta_P0cm);
      evnt->mc_helic005_phiRst[k] = get_helicity(phiR,0.05,sin_theta_P0cm);
      evnt->mc_helic010_phiRst[k] = get_helicity(phiR,0.1,sin_theta_P0cm);
      evnt->mc_helic020_phiRst[k] = get_helicity(phiR,0.2,sin_theta_P0cm);
      
      evnt->mc_wUxS_phiR[k] = get_UxS(phiR,0.2,0.15);
      evnt->mc_Mx2[k] = Mx2;
      evnt->mc_xF[k] = xF;
      evnt->mc_xF0[k] = xF0;
      evnt->mc_xF1[k] = xF1;
      evnt->mc_plcm[k] = PlCM;
      evnt->mc_plcm0[k] = Pl0CM;
      evnt->mc_plcm1[k] = Pl1CM;
      evnt->mc_Eh[k] = E;
      evnt->mc_xFm[k] = xFm;
      evnt->mc_xFm0[k] = xFm0;
      evnt->mc_xFm1[k] = xFm1;
      evnt->mc_theta0[k] = theta_0*TMath::RadToDeg();
      evnt->mc_theta1[k] = theta_1*TMath::RadToDeg();
      evnt->mc_cos_theta_P0cm[k] = cos_theta_P0cm;
      evnt->mc_sin_theta_P0cm[k] = sin_theta_P0cm;
      evnt->mc_xFo[k]  = xFo;
      evnt->mc_xFo0[k] = xF0o;
      evnt->mc_xFo1[k] = xF1o;

      Float_t dphi = (phiH-phiR);
      while (dphi<0)
	dphi+=360.;
      while (dphi>360)
	dphi-=360.;
      evnt->mc_phiH_phiR[k] = dphi;

      evnt->mc_phiR_cov[k] = phiR_cov;
      evnt->mc_p0T2[k] = p0Tv.Mag2();
      evnt->mc_p1T2[k] = p1Tv.Mag2();
      evnt->mc_phipq[k] = phi_pq;
      evnt->mc_phT2[k] = phTv.Mag2();
    
      evnt->mc_etaCM0[k] = etaCM_0;
      evnt->mc_etaCM1[k] = etaCM_1;
      evnt->mc_etaBF0p[k] = etaBF_0plus;
      evnt->mc_etaBF1p[k] = etaBF_1plus;
      evnt->mc_etaBF0m[k] = etaBF_0minus;
      evnt->mc_etaBF1m[k] = etaBF_1minus;

      evnt->mc_etaBF0[k] = etaBF0;
      evnt->mc_etaBF1[k] = etaBF1;
    
      evnt->mc_phiR_ha[k] = phiR_ha;
      evnt->mc_plcm0_r[k] = Pl0m;
      evnt->mc_plcm1_r[k] = Pl1m;
      evnt->mc_phiR_covH[k] = phiR_covH;

      evnt->mc_E0_phcm[k]  = P0_dicm.E();
      evnt->mc_E1_phcm[k] = P1_dicm.E();

      evnt->mc_R[k] = 0.5*sqrt(M*M - 4*kMpi*kMpi);
      evnt->mc_KF[k] =  mixEvent.mc_fW/mixEvent.mc_fA*kMprt/sqrt(Event.mc_Q2)*mixEvent.mc_R[k]/M;
      
      evnt->mc_mix_npart = comb->Npart;
      PARTDATA 	*pdata = (PARTDATA *)MC_PDATA_CA->ConstructedAt(k);
      for (int i =0; i<evnt->mc_mix_npart; i++){
	pdata->e[i] = (*comb)[i]->E();
	pdata->px[i] = (*comb)[i]->Px();
	pdata->py[i] = (*comb)[i]->Py();
	pdata->pz[i] = (*comb)[i]->Pz();
	pdata->helic002_phiH[i] = (*comb)[i]->helic002_phiH;
	pdata->helic005_phiH[i] = (*comb)[i]->helic005_phiH;
	pdata->helic010_phiH[i] = (*comb)[i]->helic010_phiH;
	pdata->helic020_phiH[i] = (*comb)[i]->helic020_phiH;
	Float_t phiHp = (*comb)[i]->phiHs;
	phiHp = (phiHp<0?phiHp+360:(phiHp>360?phiHp-360:phiHp));
	pdata->phiHs[i] = phiHp;
	
      }
      /*
      for (int i =0; i<evnt->mc_mix_npart; i++){
	Double_t px = (*comb)[i]->Px();
	Double_t py = (*comb)[i]->Py();
	Double_t pz = (*comb)[i]->Pz();
	Double_t e = (*comb)[i]->E();
	evnt->mc_e[k][i] = e;
	evnt->mc_px[k][i] = px;
	evnt->mc_py[k][i] = py;
	evnt->mc_pz[k][i] = pz;
	evnt->mc_helic002_phiH[k][i] = (*comb)[i]->helic002_phiH;
	evnt->mc_helic005_phiH[k][i] = (*comb)[i]->helic005_phiH;
	evnt->mc_helic010_phiH[k][i] = (*comb)[i]->helic010_phiH;
	evnt->mc_helic020_phiH[k][i] = (*comb)[i]->helic020_phiH;
	Float_t phiHp = (*comb)[i]->phiHs;
	phiHp = (phiHp<0?phiHp+360:(phiHp>360?phiHp-360:phiHp));
	evnt->mc_phiHs[k][i] = phiHp;
      }
      */
    }
    /*** end set MC data ***/
    /*** set rec data ***/
    else {
      evnt->M[k] = M;
      evnt->Phx[k] = Px;
      evnt->Phy[k] = Py;
      evnt->Phz[k] = Pz;
      evnt->Z[k] = (Float_t)E/Nu;
      evnt->Cospq[k] = cospq;
      evnt->Pt2[k] = Pt2;
      evnt->Event[k] = event;

      evnt->M2_01[k] = ((comb->Npart==3)? ( *(*comb)[0] + *(*comb)[1]).M2() : 0);
      evnt->M2_02[k] = ((comb->Npart==3)? ( *(*comb)[0] + *(*comb)[2]).M2() : 0);

      evnt->phiH[k] = phiH;
      evnt->phiR[k] = phiR;
      evnt->helic002_phiR[k] = evnt->mc_helic002_phiR[k]; 
      evnt->helic005_phiR[k] = evnt->mc_helic005_phiR[k];
      evnt->helic010_phiR[k] = evnt->mc_helic010_phiR[k];
      evnt->helic020_phiR[k] = evnt->mc_helic020_phiR[k];
      evnt->helic002_phiRst[k] = evnt->mc_helic002_phiRst[k]; 
      evnt->helic005_phiRst[k] = evnt->mc_helic005_phiRst[k];
      evnt->helic010_phiRst[k] = evnt->mc_helic010_phiRst[k];
      evnt->helic020_phiRst[k] = evnt->mc_helic020_phiRst[k];

      evnt->wUxS_phiR[k] = evnt->mc_wUxS_phiR[k];
      evnt->Mx2[k] = Mx2;
      evnt->xF[k] = xF;
      evnt->xF0[k] = xF0;
      evnt->xF1[k] = xF1;
      evnt->plcm[k] = PlCM;
      evnt->plcm0[k] = Pl0CM;
      evnt->plcm1[k] = Pl1CM;
      evnt->Eh[k] = E;
      evnt->xFm[k] = xFm;
      evnt->xFm0[k] = xFm0;
      evnt->xFm1[k] = xFm1;
      evnt->theta0[k] = theta_0*TMath::RadToDeg();
      evnt->theta1[k] = theta_1*TMath::RadToDeg();
      evnt->cos_theta_P0cm[k] = cos_theta_P0cm;
      evnt->sin_theta_P0cm[k] = sin_theta_P0cm;
      evnt->xFo[k]  = xFo;
      evnt->xFo0[k] = xF0o;
      evnt->xFo1[k] = xF1o;

      Float_t dphi = (phiH-phiR);
      while (dphi<0)
	dphi+=360.;
      while (dphi>360)
	dphi-=360.;
      evnt->phiH_phiR[k] = dphi;

      evnt->phiR_cov[k] = phiR_cov;
      evnt->p0T2[k] = p0Tv.Mag2();
      evnt->p1T2[k] = p1Tv.Mag2();
      evnt->phipq[k] = phi_pq;
      evnt->phT2[k] = phTv.Mag2();
    
      evnt->etaCM0[k] = etaCM_0;
      evnt->etaCM1[k] = etaCM_1;
      evnt->etaBF0p[k] = etaBF_0plus;
      evnt->etaBF1p[k] = etaBF_1plus;
      evnt->etaBF0m[k] = etaBF_0minus;
      evnt->etaBF1m[k] = etaBF_1minus;

      evnt->etaBF0[k] = etaBF0;
      evnt->etaBF1[k] = etaBF1;
    
      evnt->phiR_ha[k] = phiR_ha;
      evnt->plcm0_r[k] = Pl0m;
      evnt->plcm1_r[k] = Pl1m;
      evnt->phiR_covH[k] = phiR_covH;

      evnt->E0_phcm[k]  = P0_dicm.E();
      evnt->E1_phcm[k] = P1_dicm.E();
      
      evnt->R[k] = 0.5*sqrt(M*M - 4*kMpi*kMpi);
      evnt->KF[k] =  mixEvent.fW/mixEvent.fA*kMprt/sqrt(Event.Q2)*mixEvent.R[k]/M;
      
      evnt->mix_npart = comb->Npart;
      PARTDATA *pdata = (PARTDATA *)PDATA_CA->ConstructedAt(k);
      DETDATA *det = (DETDATA *)DET_CA->ConstructedAt(k);
      for (int i = 0; i < evnt->mix_npart; i++){
	pdata->e[i] = (*comb)[i]->E();
	pdata->px[i] = (*comb)[i]->Px();
	pdata->py[i] = (*comb)[i]->Py();
	pdata->pz[i] = (*comb)[i]->Pz();
	pdata->helic002_phiH[i] = (*comb)[i]->helic002_phiH;
	pdata->helic005_phiH[i] = (*comb)[i]->helic005_phiH;
	pdata->helic010_phiH[i] = (*comb)[i]->helic010_phiH;
	pdata->helic020_phiH[i] = (*comb)[i]->helic020_phiH;
	Float_t phiHp = (*comb)[i]->phiHs;
	phiHp = (phiHp<0?phiHp+360:(phiHp>360?phiHp-360:phiHp));
	pdata->phiHs[i] = phiHp;
	/////// 
	det->beta[i] = (*comb)[i]->beta;
	det->m2b[i] =  (*comb)[i]->m2b;
	det->vx[i] = (*comb)[i]->vx;
	det->vy[i] = (*comb)[i]->vy;
	det->vz[i] = (*comb)[i]->vz;
	det->dcx[i] = (*comb)[i]->dcx;
	det->dcy[i] = (*comb)[i]->dcy;
	det->dcz[i] = (*comb)[i]->dcz;
	det->statPart[i] = (*comb)[i]->statPart;
	det->dc_chi2[i] = (*comb)[i]->dc_chi2;
	det->dc_ndf[i] = (*comb)[i]->dc_ndf;
	det->pcal_lu[i] = (*comb)[i]->pcal_lu;
	det->pcal_lv[i] = (*comb)[i]->pcal_lv;
	det->pcal_lw[i] = (*comb)[i]->pcal_lw;
      }
      /*
      for (int i =0; i<evnt->mix_npart; i++){
	evnt->beta[k][i] = (*comb)[i]->beta;
	evnt->m2b[k][i] =  (*comb)[i]->m2b;
	evnt->vx[k][i] = (*comb)[i]->vx;
	evnt->vy[k][i] = (*comb)[i]->vy;
	evnt->vz[k][i] = (*comb)[i]->vz;
	evnt->dcx[k][i] = (*comb)[i]->dcx;
	evnt->dcy[k][i] = (*comb)[i]->dcy;
	evnt->dcz[k][i] = (*comb)[i]->dcz;
	evnt->statPart[k][i] = (*comb)[i]->statPart;
	evnt->dc_chi2[k][i] = (*comb)[i]->dc_chi2;
	evnt->dc_ndf[k][i] = (*comb)[i]->dc_ndf;
	evnt->pcal_lu[k][i] = (*comb)[i]->pcal_lu;
	evnt->pcal_lv[k][i] = (*comb)[i]->pcal_lv;
	evnt->pcal_lw[k][i] = (*comb)[i]->pcal_lw;
	Double_t px = (*comb)[i]->Px();
	Double_t py = (*comb)[i]->Py();
	Double_t pz = (*comb)[i]->Pz();
	Double_t e = (*comb)[i]->E();
	evnt->e[k][i] = e;
	evnt->px[k][i] = px;
	evnt->py[k][i] = py;
	evnt->pz[k][i] = pz;
	evnt->helic002_phiH[k][i] = evnt->mc_helic002_phiH[k][i];
	evnt->helic005_phiH[k][i] = evnt->mc_helic005_phiH[k][i];
	evnt->helic010_phiH[k][i] = evnt->mc_helic010_phiH[k][i];
	evnt->helic020_phiH[k][i] = evnt->mc_helic020_phiH[k][i];
	Float_t phiHp = (*comb)[i]->phiHs;
	phiHp = (phiHp<0?phiHp+360:(phiHp>360?phiHp-360:phiHp));
	evnt->phiHs[k][i] = phiHp;
      }
      */
    }
    /*** end set rec data ***/
    //    std::cout<<"### "<<evnt->beta[k][0]<<std::endl<<std::endl;
    return 0;
  }
  
  Bool_t FidCheck(int pid)
  {
    return true;
  }
  //////////////////////////
  inline   Double_t kinFit(Double_t *W, Double_t *Wa, TMatrixD &V ){
    bool db=false;
    if (db)std::cout<<"V Matrix on  "<<__LINE__<<std::endl;
    if (db)V.Print();


    //    TMatrixD V(4,4) = *Vm;//covariance matrix
    TMatrixD V_new(4,4);//covariance matrix updated
    TMatrixD VD(1,1);//auxiliary matrix
    TMatrixD D(1,4);//derivatives of restriction on xa
    TMatrixD d(1,1);//restriction on xa
    TMatrixD x(4,1);//corrected point
    TMatrixD xa(4,1);//linearization point
    TMatrixD x0(4,1);// initial point
    TMatrixD lam(1,1);//constrain multiplier
    Double_t Chi2 =0;//chi2 value
    
   
    D(0,0) = -2*Wa[0]; D(0,1) = -2*Wa[1]; D(0,2) = -2*Wa[2]; D(0,3) = 2*Wa[3];
    if (db)std::cout<<"D Matrix on  "<<__LINE__<<std::endl;
    if (db)D.Print();
    //    x(0,0) = W[0]; x(1,0) = W[1]; x(2,0) = W[2]; x(3,0) = W[3];
    xa(0,0) = Wa[0]; xa(1,0) = Wa[1]; xa(2,0) = Wa[2]; xa(3,0) = Wa[3];
    x0(0,0) = W[0]; x0(1,0) = W[1]; x0(2,0) = W[2]; x0(3,0) = W[3];
    
    //lam(0,0) = 0;
    d(0,0)=Wa[3]*Wa[3] - Wa[0]*Wa[0] - Wa[1]*Wa[1] - Wa[2]*Wa[2] - kPdgInfo->Mass()*kPdgInfo->Mass();
    
    TMatrixD DT =(D).Transpose(D);
    if (db)std::cout<<"DT Matrix on  "<<__LINE__<<std::endl;
    if (db) DT.Print();
    D.Transpose(D);
    
    TMatrixD VDL = (D*V);
    if (db)std::cout<<"VDL Matrix on  "<<__LINE__<<std::endl;
    if (db)VDL.Print();
    VD=VDL*DT;
    if (db)std::cout<<"VD Matrix on  "<<__LINE__<<std::endl;
    if (db)VD.Print();
    VD.Invert();
    if (db)std::cout<<"VD Matrix on  "<<__LINE__<<std::endl;
    if (db)VD.Print();
    lam = VD*(D*(x0-xa) + d);
    
    TMatrixD lamT = lam.Transpose(lam);
    lam.Transpose(lam);
    
    x = x0 + (V*DT)*lam;
    W[0]=x(0,0);W[1]=x(1,0);W[2]=x(2,0);W[3]=x(3,0);
    V_new = V - (V*DT)*VD*(D*V);
    V=V_new;
    TMatrixD aux=(lamT*(D*(x0 - xa) + d));
    //aux.Print();
    Chi2=aux(0,0);
    
    return Chi2;
}
  int pushSecondary(int pid)
  {
    kNSecondary++;
    if (!isPid(pid))
    {
      kSPid.push_back(pid);
      kNSPid[pid]=1;
      int k = kSPid.size();
      hSPid.push_back(new TH1F(Form("hNPart%d",k), TDatabasePDG::Instance()->GetParticle(pid)->GetName(),20,0,20) );
    }
    else
      kNSPid[pid]++;


    return pid;
  }

  int addSecondary(int pid)
  {
    pid = TDatabasePDG::Instance()->GetParticle(pid)->PdgCode();
    pushSecondary(pid);
    return pid;
  }

  int addSecondary(const char *name)
  {
    int pid  = TDatabasePDG::Instance()->GetParticle(name)->PdgCode();
    pushSecondary(pid);
    return pid;
  }


  int addPrimary(int pid){kPdgInfo = TDatabasePDG::Instance()->GetParticle(pid); kPPid=kPdgInfo->PdgCode(); return kPPid;}
  int addPrimary(const char *name){kPdgInfo = TDatabasePDG::Instance()->GetParticle(name); kPPid=kPdgInfo->PdgCode(); return kPPid;}


  bool isPid(int pid)
  {
    kSIt = std::find (kSPid.begin(), kSPid.end(), pid);
    if (kSIt!=kSPid.end())
      return true;
    else
      return false;
  }
  bool checkMinPart()
  {
    bool minPart=true;
    for (int k=0;k<(int)kSPid.size();k++)
    {
      if (fEMatch)
	minPart=minPart && ((int)kSecondary[k].size() == kNSPid[kSPid[k]]);
      else
	minPart=minPart && ((int)kSecondary[k].size() >= kNSPid[kSPid[k]]);

      //      std::cout<<__LINE__<<": minPart "<<minPart<<" : "<<kSecondary[k].size()<<std::endl;
      hSPid[k]->Fill(kSecondary[k].size());
    }
    return minPart;
  }

  bool checkMinPart(Combo *c)
  {
    bool minPart=true;
    for (int k=0;k<(int)kSPid.size();k++)
    {
	minPart=minPart && (c->findPid(kSPid[k]) == kNSPid[kSPid[k]]);
    }
    return minPart;
  }


  //int takeN(int N,int kspid, int pos=0,Particle p=Particle(),int count=0)
  //  int takeN(int N,int kspid, int pos=0,Combo *c= new Combo() ,int count=0)
  int takeN(int N, int kspid, int pos = 0, Combo *c = 0, int count = 0)
  {
    Combo *c_new;
    if (DEBUG) std::cout<<"### take N:  "<<N<<"##############"<<std::endl; 
    if (N<1) return -1;
    if (N!=1)
    {
      for (int k =pos;k<(int)kSecondary[kspid].size()-N+1;k++)
      {
	if (c==0) c_new = new Combo();
	else c_new = new Combo(*c);

	c_new->addParticle(kSecondary[kspid][pos]);
	count=takeN(N-1,kspid,k+1,c_new,count);
      }
      //if (c!=0) delete c;
    }
    
    else
    {
      for (int k=pos;k<(int)kSecondary[kspid].size();k++)
      {
	if (c==0) c_new = new Combo();
	else c_new = new Combo(*c);

	c_new->addParticle(kSecondary[kspid][k]);
	
	if (DEBUG) std::cout<<"############ Npart from takeN: "<<c_new->Npart<<"#### pid: "<<kSPid[kspid]<<"#############"<<std::endl;
	//	kCombo[kspid].push_back(new Combo(*c) );
	kCombo[kspid].push_back(c_new);
	//kParticles.push_back(new Particle(*kSecondary[kspid][k]));
	//std::cout<<__LINE__<<" "<< kCombo[kspid].back()->M()<<std::endl;
	count++;

      }
    }
    if (c!=0) delete c;
    return count;
  }
  int findSecondary(Particle *p)
  {
    Bool_t notused = kTRUE;
    int count = 0;
    for (int k =0;k<(int)kSPid.size();k++)
    {
      //if(pid == kSPid[k]&&checkDCFidPhi())
      if(p->pid == kSPid[k])
      {
	kSecondary[k].push_back(p);
	//std::cout<<dcx_rot_0<<" / "<<dcy_rot_0<<" / "<<dcz_r0<<std::endl;
	count++;
	notused = kFALSE;
	break;
      }
    }
    if (notused) delete p;
    return count;
  }
  /*
  int correct_momentum()
  {
    Float_t Rt = TMath::Sqrt( ECX*ECX + ECY*ECY );
    Float_t R = TMath::Sqrt( ECX*ECX + ECY*ECY + (ECZ-Ze)*(ECZ -Ze) );
    Float_t theta_gam=TMath::ASin(Rt/R);
    Float_t phi_gam = TMath::ATan2(ECY,ECX);
    Px=Ep*TMath::Sin(theta_gam)*TMath::Cos(phi_gam);
    Py=Ep*TMath::Sin(theta_gam)*TMath::Sin(phi_gam);
    Pz=Ep*TMath::Cos(theta_gam);
    
    return 0;
  }
  
  Bool_t checkDCFidPhi()
  {
    Bool_t ret=false;
    Float_t phi = atan2(DCPy,DCPx);
    Float_t rho = sqrt( DCPx*DCPx + DCPy*DCPy );
    for (int k=0;k<6;k++)
    {
      ret = ret || (-2.8 + k*1.04 <phi&&phi<-2 + k*1.04);
    }
    ret = ret && (rho>0.1);
    return ret;
  }
  */
  Float_t getE(Float_t Es=0, Float_t p=0, Float_t pid=211)
  {
    //Float_t pip[4] = {-0.001114,0.1146,1.103e-5,-0.01259};
    //Float_t pim[4] = {0.007386,0.09922,-0.001244,-0.01057};

    Float_t energy =-1;

    if (pid == 45)
      energy = sqrt(p*p + TMath::Power(1.8756,2) );
    else
      energy = sqrt(p*p + TMath::Power(TDatabasePDG::Instance()->GetParticle(pid)->Mass(),2) );
    //    Float_t sf = ((pid==211)?pip[0] + pip[1]/x + pip[2]*x + pip[3]/x/x:((pid==-211)?pim[0] + pim[1]/x + pim[2]*x + pim[3]/x/x:0)) ;
    // Float_t energy = Es/sf;

    return energy;
  }
  
  int setOutVars(bool isMC = false){
    Float_t pid,E,P,Px,Py,Pz,vxh,vyh,vzh,beta,dcx_rot_0,dcy_rot_0,trajz_sl0,statPart,dc_chi2,dc_ndf,pcal_lu,pcal_lv,pcal_lw,phiHs;
    Int_t npart = 0;
    Float_t Ep = 0;
    
    if (!isMC)
      npart = Event.npart;
    else
      npart = Event.mc_npart;
    /*** rec particles ***/
    for ( int k = 0; k < npart; k++){
      if (!isMC){
	pid = Event.pid[k];
	E = Event.E[k];
	P = Event.P[k];
	Px = Event.Px[k];
	Py = Event.Py[k];
	Pz = Event.Pz[k];
	vxh = Event.vxh[k];
	vyh = Event.vyh[k];
	vzh = Event.vzh[k];
	beta = Event.Beta[k];
	dcx_rot_0 = Event.dcx_rot_0[k];
	dcy_rot_0 = Event.dcy_rot_0[k];
	trajz_sl0 = Event.trajz_sl0[k];
	statPart = Event.statPart[k];
	dc_chi2 = Event.dc_chi2[k];
	dc_ndf = Event.dc_ndf[k];
	pcal_lu = Event.pcal_lu[k];
	pcal_lv = Event.pcal_lv[k];
	pcal_lw = Event.pcal_lw[k];
	phiHs = Event.PhiPQ[k];
      }
      else{
	pid = Event.mc_pid[k];
	E = Event.mc_E[k];
	P = Event.mc_P[k];
	Px = Event.mc_Px[k];
	Py = Event.mc_Py[k];
	Pz = Event.mc_Pz[k];
	vxh = Event.mc_vxh[k];
	vyh = Event.mc_vyh[k];
	vzh = Event.mc_vzh[k];
	beta = Event.mc_Beta[k];
	dcx_rot_0 = 0;
	dcy_rot_0 = 0;
	trajz_sl0 = 0;
	statPart = 0;
	dc_chi2 = 0;
	dc_ndf = 0;
	pcal_lu = 0;
	pcal_lv = 0;
	pcal_lw = 0;
	phiHs = Event.mc_PhiPQ[k];
      }

      if (isMC)
	Ep = E; //gsim
      else
	Ep = (pid==22)? (E/0.23):getE(E,P,pid);
      
      if (FidCheck(11)){// no fid
	Particle *p = new Particle (Px,Py,Pz,Ep,vxh,vyh,vzh,pid,0,beta,dcx_rot_0,dcy_rot_0,trajz_sl0,statPart,dc_chi2,dc_ndf,pcal_lu,pcal_lv,pcal_lw,phiHs,get_helicity(phiHs,0.02),get_helicity(phiHs,0.05),get_helicity(phiHs,0.1),get_helicity(phiHs,0.2));
	//if (!isMC)push_bkgnd(p);
	findSecondary(new Particle(*p));
	delete p;
      }
    }

    /*** end adding particles to the set***/
    
    /*** get comb ***/
    //    if (!isMC) pop_bkgnd();
    if (isMC) mixEvent.mc_npart = 0;
    else mixEvent.npart = 0;
      
    if (checkMinPart()){
      int Npart = 1;
	
      for (int k =0; k<(int)kSPid.size(); k++){
	if (DEBUG) std::cout<<"############ Nsecondary of pid: "<<kSPid[k]<<" ::: "<<kSecondary[k].size()<<"###########"<<std::endl; 
	takeN(kNSPid[kSPid[k]] ,k);
	  
	Npart *= kCombo[k].size();
      }
      if (DEBUG) std::cout<<"############ N candidates for primary: "<<Npart<<"####################\n##########################"<<std::endl;
      if (isMC) mixEvent.mc_npart = Npart;
      else mixEvent.npart = Npart;
      
      for(int k=0; k<Npart; k++){
	kPrimary = new Combo();
	int div=1;
	for (int l=0; l<(int)kSPid.size(); l++){
	  int size = kCombo[l].size();
	  //std::cout<<__LINE__<<": \n";
	  *kPrimary+=*kCombo[l][ (k/div)%size ];
	  //std::cout<<__LINE__<<": \n";
	    
	  div*=size;
	}
	//  if (DEBUG) std::cout<<"############ Npart from primary: "<<kPrimary->Npart<<"####################"<<std::endl;
	setVars(k,isMC);
	delete kPrimary;
      }
      setElecVars(isMC); 
      //if (mixEvent.npart>0) std::cout<<"#### "<<mixEvent.beta[0][0]<<" ####"<<std::endl<<std::endl;
    }
    
    clear();
    return 0;
  }
  
  int getCombinations(TChain *t, long st = 0){
    kSecondary = new std::vector<Particle*> [kSPid.size()];
    kCombo = new std::vector<Combo*> [kSPid.size()];

    std::cout<<"processing...\n";
    std::cout.fill('.');
    
    std::cout<<"start / Ne: "<<st<<" / "<<Ne;
    
    for ( Long64_t i = st; i < Ne; i++){
      t->GetEntry(i);
      if (Event.mc_npart>0) {
	kEbeam = Event.mc_Nu + Event.mc_Pe;
      }
      else  {
	kEbeam = Event.Nu + Event.Pe;
      }
      resetDATAMIX(&mixEvent);
      setOutVars(true); // MC
      setOutVars();
      /*** end rec particles ***/
      if (mixEvent.npart>0 || mixEvent.mc_npart>0)
	fillTree();

      //	if (data_type<2&&FidCheck(pid))
      std::cout<<std::setw(15)<<float(i+1)/Ne*100<<" %"<<"\r";
      std::cout.flush();
      
    }
    return 0;
  }
  int push_bkgnd(Particle *p)
  {
    if (isPid(p->pid) &&  kBkgnd.size() < kBUFFLIMIT)
    { 
      if (!kBkgnd.empty())
      {
	int i;
	for (i =0;i<(int)kBkgnd.size();i++)
	{
	  if ( (kBkgnd[i]->findPid(p->pid)<kNSPid[p->pid] ) && (kBkgnd[i]->lastEvent!=Event.event) && (kBkgnd[i]->isCompatible()) )
	  {
	    kBkgnd[i]->addParticle(p,Event.event,kTRUE);
	    break;
	  }
	  
	}
	if (i==(int)kBkgnd.size())
	{
	  Combo *c = new Combo();
	  c->addParticle(p,Event.event);
	  kBkgnd.push_back(c);
	}
      }
      else 
      {
	Combo *c = new Combo();
	c->addParticle(p,Event.event);
	kBkgnd.push_back(c);
      }
    }
    return 0;
  }

  bool pop_bkgnd(){
    bool ret = false;
    if (!kBkgnd.empty()){
      int cnt=0;
      for (int i =0;i<(int)kBkgnd.size();i++){
	if (checkMinPart(kBkgnd[i])){
	  //	  std::cout<<__FILE__<<"::"<<__LINE__<<std::endl;
	  setVars(kBkgnd[i],&mixEvent_bkg,cnt++);
	  delete kBkgnd[i];
	  kBkgnd.erase(kBkgnd.begin()+i);
	  ret = true;
	}
      } 
    } 
    return ret;
  }

  /*
  int getBkgnd(TChain *t)
  {
    
    t->GetEntry(0);
    setElectVar();
    Event.event=Event.event;
    std::cout<<"processing...\n";
    std::cout.fill('.');
    for ( int i = 0; i < Ne; i++)
    {
      t->GetEntry(i);
      Ep = (pid==22)? (E/0.23):sqrt(P*P+ TMath::Power(TDatabasePDG::Instance()->GetParticle(pid)->Mass(),2));
      
      if (pid==22)
	correct_momentum();

      if (Event.event==Event.event)
      {
	//std::cout<<__LINE__<<" "<<findSecondary()<<std::endl;
	Particle *p = new Particle(Px,Py,Pz,Ep,vx,vy,vz,pid);
	push_bkgnd(p);
      }
      else
      {
	pop_bkgnd();	
	setElectVar();
	//

	Particle *p = new Particle(Px,Py,Pz,Ep,vx,vy,vz,pid);
	push_bkgnd(p);
      }
      std::cout<<std::setw(15)<<float(i+1)/Ne*100<<" %"<<"\r";
      std::cout.flush();
    } 
    return 0;
  }*/
};



void check_dir(const char *outdir)
{
  struct stat sb;
  if (stat(outdir, &sb) != 0)
  {
      system(Form("mkdir %s",outdir));
  }
  
}

int main(int argc, char *argv[])
{
  TBenchmark *bm = new TBenchmark();
  bm->Start("main_program");
  parseopt(argc,argv);

  time_t seconds;
  seconds = time(NULL);
  rndm = new TRandom3(seconds);

  TChain *t = new TChain();

  TString PATH="";
  if (!INDIR.IsNull()) PATH = INDIR + "/*.root";
  else PATH = INFILE;

  std::cout<<"processing: "<<PATH<<std::endl;
  std::cout<<"output file name: "<<OFILE<<std::endl;
  std::cout<<"Exact match: "<<EFLAG<<std::endl;
  std::cout<<"Reaction: "<<REACTION<<std::endl;
  
  t->Add(PATH + "/evData");

  setTreeAddress(t,&Event);
  if (Ne==-1) Ne = t->GetEntries();
  else if ( t->GetEntries() < (START+Ne) ) Ne = t->GetEntries();
  else Ne += START;

  std::cout<<"Number of entries to be processed: "<<Ne<<std::endl;
  std::cout<<"Starting event number: "<<START<<std::endl;

 
  if (!check_reaction())
    return 1;
  
  TString rname = REACTION + (EFLAG?", exact-match":", all combinations");
  Reaction r(rname,OFILE,EFLAG);
  TString primary = get_primary(),secondary="";
  Ssiz_t start= REACTION.Index(":")+1;
  if (primary!="")
    r.addPrimary(primary);
  std::cout<<primary<<" -> ";
  
  while ((secondary=pop_secondary(start)) != ""){
    std::cout<<secondary<<" ";
    r.addSecondary(secondary);
  }
  std::cout<<std::endl;

  r.getCombinations(t,START);
  r.store();
  std::cout<<"\n";
  r.kOutData->Print();
  //  r.kOutBkgnd->Print();
  //  corrfile->Close();
  bm->Show("main_program");
  delete bm;
  return 0;  

}
