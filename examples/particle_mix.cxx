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
#include "particle_mix.h"
#define MAXPART 10

bool DEBUG=false;
Float_t HELIC=-111;
TString INFILE="data.root", INDIR="";
TH1F *hW;
TH1F *hW2;
TH1F *hWmb;
TH1F *hW2mb;
TH1F *hT;
//TH1F *hcfm;
TH2F *hEpi0_th;
Float_t kMpi = TDatabasePDG::Instance()->GetParticle("pi-")->Mass();
Float_t kMPi0=1.33196e-1;//Got from sim.
Float_t kSPi0=1.94034e-2;//Got from sim.
Float_t kMh =  TDatabasePDG::Instance()->GetParticle("pi-")->Mass();
//Float_t kMPi0=5.39609e-01;
//Float_t kSPi0=5.98542e-02;
bool GSIM=false;
int data_type=0;
Float_t kPt2,kEvent; 
TNtuple *tuple;
//TNtuple *tuple_sim;
TNtuple *tuplemb;
TNtuple *tuplePi0_gamma, *tupleGamma;
Float_t kEbeam=10.2,E,Ee,Ee_prev,Beta,Ep,P,Px,Py,Pz,evnt,evnt_prev,revent,revent_prev,Ze,Ze_prev,Ye,Ye_prev,Xe,Xe_prev,TEc,Q2,Q2_prev,W,W_prev,Nu,Nu_prev,helic,helic_prev,Pex,Pex_prev,Pey,Pey_prev,Pez,Pez_prev,TargType,TargType_prev,TargTypeO=0,TargTypeO_prev=0,pid,vx,vy,vz,DCX,DCY,DCZ,ECX,ECY,ECZ,DCPx,DCPy,DCPz,dcx_r0,dcy_r0,dcz_r0;

typedef struct {
  Int_t npart;
  Float_t beta[MAXPART];
  Float_t m2b[MAXPART];
  Float_t vx[MAXPART];
  Float_t vy[MAXPART];
  Float_t vz[MAXPART];
  Float_t dcx[MAXPART];
  Float_t dcy[MAXPART];
  Float_t dcz[MAXPART];
} detData_t;


long Ne = -1;
char st[3]= "C"; // solid target: C Fe Pb
char tt[3] = "C"; // cut on solid target or Deuterium : (st) or D.

Float_t kMprt=0.938272, kMn =0.939565;// kMprt to avoid replace confusion with kMpi
TClonesArray *P4Arr;


class Particle: public TLorentzVector
{
public:
  Float_t vx,vy,vz,pid,time,beta,dcx,dcy,dcz,m2b;
  inline Double_t P2() const {return P()*P();}
  //TParticlePDG *info;
  Particle() : TLorentzVector(), vx(0),vy(0),vz(0),pid(0),time(0),beta(0),dcx(0),dcy(0),dcz(0){m2b = (P2() - beta*beta*P2())/(beta*beta);}
  Particle(Float_t px,Float_t py, Float_t pz, Float_t e, Float_t x, Float_t y, Float_t z, Float_t pid=0, Float_t t=0, Float_t b =0, Float_t dx=0, Float_t dy=0, Float_t dz=0): TLorentzVector(px,py,pz,e),vx(x),vy(y),vz(z),pid(pid),time(t),beta(b),dcx(dx),dcy(dy),dcz(dz){m2b = (P2() - beta*beta*P2())/(beta*beta);}
  Particle(TLorentzVector lv, Float_t x=0, Float_t y=0, Float_t z=0, Float_t pid=0, Float_t t=0, Float_t b =0, Float_t dx=0, Float_t dy=0, Float_t dz=0): TLorentzVector(lv),vx(x),vy(y),vz(z),pid(pid),time(t),beta(b),dcx(dx),dcy(dy),dcz(dz){m2b = (P2() - beta*beta*P2())/(beta*beta);}
  Particle(Particle &p):vx(p.vx),vy(p.vy),vz(p.vz),pid(p.pid),time(p.time),beta(p.beta),dcx(p.dcx),dcy(p.dcy),dcz(p.dcz),m2b(p.m2b) {SetVect(p.Vect()); SetT(p.T());}

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
      kQ2= Q2; // using global variable, must be changed!.
      kNu= Nu; // using global variable, must be changed!.
      q4 = new TLorentzVector(-Pex,-Pey,kEbeam-Pez,kEbeam-Ee);
    }
    else if (rotfirst)
    {
      TLorentzVector *q4n = new TLorentzVector(-Pex,-Pey,kEbeam-Pez,kEbeam-Ee);
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
    if ( !( (-kQ2Tol< (kQ2 - Q2) && (kQ2 - Q2) < kQ2Tol) && (-kNuTol< (kNu - Nu) && (kNu - Nu) < kNuTol) ) )
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
  detData_t detData;

  TParticlePDG *kPdgInfo;
  std::vector<Particle*> *kSecondary; //all secondary

  //std::vector<Particle*> *kCombo; //partial combinations.
  std::vector <Combo *> *kCombo;//partial combinations.
  std::vector <Combo *> kBkgnd;//background combinations.

  static const Int_t kBUFFLIMIT=100;

  std::vector<TH1F*> hSPid;
  TFile *kOutFile;
  Float_t *kData;
  Float_t *keData;
  
  TNtuple  *kElecData;

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
    for (int k=0;k<kSPid.size();k++)
    {
      hSPid[k]->Write("",TObject::kOverwrite);
    }
    kOutBkgnd->Write("",TObject::kOverwrite);
    kElecData->Write("",TObject::kOverwrite);
    //kP4Tree->Write("",TObject::kOverwrite);
    return kOutData->Write("",TObject::kOverwrite);
  }
  ~Reaction()
  {
    clear();
    delete[] kSecondary;
    delete[] kCombo;
    delete kOutData;
    delete kOutBkgnd;
    delete[] kData;
    delete[] keData;
    delete P4Arr;
    kOutFile->Close();
    delete kOutFile;
    delete kPdgInfo;
    hSPid.clear();
  }
  void clear()
  {
    for (int k=0;k<kSPid.size();k++)
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
    TString varlist="M:Phx:Phy:Phz:Nu:Q2:Z:Cospq:Pt2:Event:M2_01:M2_02:vzec:z1:z2:z3:W:vxec:vyec:qx1:qy1:qz1:qx2:qy2:qz2:E1:E2:E1c:E2c:x1:y1:x2:y2:TargType:TargTypeO:phiH:phiR:Mx2:xF:xF0:xF1:plcm:plcm0:plcm1:Eh:xFm:xFm0:xFm1:helic:theta0:theta1:cos_theta_P0cm:xFo:xF0o:xF1o:event:phiH_phiR:Ee:phiR_cov:p0T2:p1T2:phipq:phT2:revent:etaCM0:etaCM1:etaBF0p:etaBF1p:etaBF0m:etaBF1m:etaBF0:etaBF1:phiH0:phiH1:qx:qy:qz:phiR_ha:plcm0_r:plcm1_r:phiR_covH:E0_phcm:E1_phcm:y:th_e";
    kElecData = new TNtuple("ElecData",Form("%s",name),"Nu:Q2:Event:vze:Ee:Pex:Pey:Pez:W");


    kData = new Float_t[varlist.CountChar(':')+1];
    keData = new Float_t[kElecData->GetNvar()];
    P4Arr = new TClonesArray("TLorentzVector");
    kOutData = new TTree("outdata",Form("%s",name));
    kOutData->Branch("P4",&P4Arr,6);
    kOutData->Branch("primary",kData,varlist);

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
    ////////////
    
    kOutBkgnd = kOutData->CloneTree(0);
    kOutBkgnd->SetName("outbkgnd");
    
    
    return 0;
  }

  int fill()
  {
    fill(kPrimary,kOutData);
  }

  int fill_elec()
  {

    keData[0] = Nu_prev;
    keData[1] = Q2_prev;
    keData[2] = evnt_prev;
    keData[3] = Ze_prev;
    keData[4] = Ee_prev;
    keData[5] = Pex_prev;
    keData[6] = Pey_prev;
    keData[7] = Pez_prev;
    keData[8] = W_prev;
    kElecData->Fill(keData);
  }


  int fill(Combo *comb,TTree *ttree)
  {
    if (DEBUG )std::cout<<"##### combo size: "<<comb->size()<<" #######"<<std::endl;
    
    Double_t Px = comb->Px();
    Double_t Py = comb->Py();
    Double_t Pz = comb->Pz();
    Double_t E = comb->E();
    Double_t P2 = comb->P2();
    Double_t M2 = comb->M2();
    Double_t M =  (M2>=0)?TMath::Sqrt(M2):-1.0;
    Float_t theta=Pz/sqrt(P2);
    Float_t theta_0 = (*comb)[0]->Theta();
    Float_t theta_1 = (*comb)[1]->Theta();

    Float_t cospq = ((kEbeam-Pez_prev)*Pz - Pex_prev*Px - Pey_prev*Py)/( sqrt((Q2_prev + Nu_prev*Nu_prev)*P2) );

    Float_t cospq0 = ((kEbeam-Pez_prev)*(*comb)[0]->Pz() - Pex_prev*(*comb)[0]->Px() - Pey_prev*(*comb)[0]->Py() )/( sqrt((Q2_prev + Nu_prev*Nu_prev))*(*comb)[0]->P());
    Float_t cospq1 = ((kEbeam-Pez_prev)*(*comb)[1]->Pz() - Pex_prev*(*comb)[1]->Px() - Pey_prev*(*comb)[1]->Py())/( sqrt((Q2_prev + Nu_prev*Nu_prev))*(*comb)[1]->P());
 
    Float_t Pt2 = P2*(1-cospq*cospq);
    Float_t Pl2 = P2*cospq*cospq;

    Float_t Pl = sqrt(P2)*cospq;
    Float_t Pl0 = (*comb)[0]->P()*cospq0;
    Float_t Pl1 = (*comb)[1]->P()*cospq1;

    Double_t phi_pq;
    TVector3 Vhad(Px,Py,Pz);
    TVector3 Vvirt(-Pex_prev,-Pey_prev,kEbeam-Pez_prev);
    Double_t phi_z = TMath::Pi()-Vvirt.Phi();
    Vvirt.RotateZ(phi_z);
    Vhad.RotateZ(phi_z);
    TVector3 Vhelp(0.,0.,1.);
    Double_t phi_y = Vvirt.Angle(Vhelp);
    Vvirt.RotateY(phi_y);
    Vhad.RotateY(phi_y);
    phi_pq=Vhad.Phi() * 180./(TMath::Pi());

    TLorentzVector q_lv(-Pex_prev,-Pey_prev,kEbeam-Pez_prev,Nu_prev); // virtual photon 4vec
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
    
    Float_t xbj = Q2_prev/2/Nu_prev/kMprt;
    //Nachtmann variable
    Float_t xn = 2*xbj/(1 + sqrt(1 + 4*xbj*xbj*kMprt*kMprt/Q2_prev));
    Float_t Q = sqrt(Q2_prev);

    Float_t z0 = (*comb)[0]->E()/Nu_prev;
    Float_t z1 = (*comb)[1]->E()/Nu_prev;
    Float_t ztot_inv = 1./(z0+z1);
    
    TVector3 p0Lv_BF = (P0_BF.Vect()*q_BF.Vect().Unit())*q_BF.Vect().Unit();
    TVector3 p1Lv_BF = (P1_BF.Vect()*q_BF.Vect().Unit())*q_BF.Vect().Unit();
    TVector3 p0Tv_BF = P0_BF.Vect() - p0Lv_BF;
    TVector3 p1Tv_BF = P1_BF.Vect() - p1Lv_BF;

    // pip /////    
    Float_t MhT20 = p0Tv_BF.Mag2() + kMh*kMh;
    Float_t etaBF_0plus = Q*z0*(Q2_prev - xn*xn*kMprt*kMprt)/
      (2*xn*xn*kMprt*kMprt*sqrt(MhT20));

    Float_t etaBF_0minus = Q*z0*(Q2_prev - xn*xn*kMprt*kMprt)/
      (2*xn*xn*kMprt*kMprt*sqrt(MhT20));

    etaBF_0plus += Q/(kMprt*xn)*sqrt( z0*z0*TMath::Power((Q2_prev - xn*xn*kMprt*kMprt),2)/(4*xn*xn*kMprt*kMprt*MhT20) -1);

    etaBF_0minus -= Q/(kMprt*xn)*sqrt( z0*z0*TMath::Power((Q2_prev - xn*xn*kMprt*kMprt),2)/(4*xn*xn*kMprt*kMprt*MhT20) -1); 

    etaBF_0plus = TMath::Log(etaBF_0plus);
    etaBF_0minus = TMath::Log(etaBF_0minus);

    //////// end pip ////////

    // pim /////    

    Float_t MhT21 = p1Tv_BF.Mag2() + kMh*kMh;
    Float_t etaBF_1plus = Q*z1*(Q2_prev - xn*xn*kMprt*kMprt)/
      (2*xn*xn*kMprt*kMprt*sqrt(MhT21));

    Float_t etaBF_1minus = Q*z1*(Q2_prev - xn*xn*kMprt*kMprt)/
      (2*xn*xn*kMprt*kMprt*sqrt(MhT21));

    etaBF_1plus += Q/(xn*kMprt)*sqrt( z1*z1*TMath::Power((Q2_prev - xn*xn*kMprt*kMprt),2)/(4*xn*xn*kMprt*kMprt*MhT21) -1);

    etaBF_1minus -= Q/(xn*kMprt)*sqrt( z1*z1*TMath::Power((Q2_prev - xn*xn*kMprt*kMprt),2)/(4*xn*xn*kMprt*kMprt*MhT21) -1); 

    etaBF_1plus = log(etaBF_1plus);
    etaBF_1minus = log(etaBF_1minus);

    //////// end pim ////////
    
    
    
    Float_t cos_theta_P0cm = P0_dicm.Vect()*Ph.Vect()/P0_dicm.Vect().Mag()/Ph.Vect().Mag();

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
    Float_t g2    = 4*xbj*xbj*kMprt*kMprt/Q2_prev;
    Float_t g24   = xn*xn*kMprt*kMprt/Q2_prev;
    Float_t yhc1  = xn*xn*kMprt*kMprt+xn*Q2_prev;
    Float_t yhc2  = (1.0-xn)*Q2_prev;
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
  
    Float_t Mx2 = W_prev*W_prev + M*M - 2*( (Nu_prev+kMprt)*E - sqrt((Q2_prev + Nu_prev*Nu_prev)*Pl2));

    Float_t Phmax = TMath::Sqrt( TMath::Power( W_prev*W_prev + M*M - kMprt*kMprt, 2 ) - 4*W_prev*W_prev*M*M )/(2* W_prev);

    Float_t Ph0max = TMath::Sqrt( TMath::Power( W_prev*W_prev + kMpi*kMpi - kMprt*kMprt - 2*W_prev*(*comb)[1]->E(), 2 ) - 4*W_prev*W_prev*kMpi*kMpi )/(2* W_prev);
    Float_t Ph1max = TMath::Sqrt( TMath::Power( W_prev*W_prev + kMpi*kMpi - kMprt*kMprt - 2*W_prev*(*comb)[0]->E(), 2 ) - 4*W_prev*W_prev*kMpi*kMpi )/(2* W_prev);

    ////// LORENTZ BOOST //////////
    Float_t b=TMath::Sqrt(Q2_prev + Nu_prev*Nu_prev)/(Nu_prev+kMprt);
    Float_t g=(Nu_prev+kMprt)/W_prev;

    Float_t PlCM = g*(Pl - b*E);
    Float_t Pl0CM = g*( Pl0 - b*(*comb)[0]->E() );
    Float_t Pl1CM = g*( Pl1 - b*(*comb)[1]->E() );
      
    Float_t xF = PlCM/Phmax;
    Float_t xF0 = Pl0CM/Ph0max;
    Float_t xF1 = Pl1CM/Ph1max;

    Float_t xFm = 2*PlCM/W_prev;
    Float_t xFm0 = 2*Pl0CM/W_prev;
    Float_t xFm1 = 2*Pl1CM/W_prev;

    Float_t xFo = 2*Plt/W_prev;
    Float_t xF0o = 2*Pl0m/W_prev;
    Float_t xF1o = 2*Pl1m/W_prev;
    
    kData[0] = M;
    kData[1] = Px;
    kData[2] = Py;
    kData[3] = Pz;
    kData[4] = Nu_prev;
    kData[5] = Q2_prev;
    kData[6] = (Float_t)E/Nu_prev;
    kData[7] = cospq;
    kData[8] = Pt2;
    kData[9] = evnt_prev;

    //    if (DEBUG) std::cout<<std::endl<<__LINE__<<": Npart:"<<comb->kParticles.size()<<std::endl;

    kData[10]=((comb->Npart==3)? ( *(*comb)[0] + *(*comb)[1]).M2() : 0);
    kData[11]=((comb->Npart==3)? ( *(*comb)[0] + *(*comb)[2]).M2() : 0);
    kData[12] = Ze_prev; // z electron corrected.

    kData[13] = (*comb)[0]->vz;
    kData[14] = ((comb->Npart>1)?(*comb)[1]->vz:0);
    kData[15] = ((comb->Npart>2)?(*comb)[2]->vz:0);
    kData[16] = W_prev;
    kData[17] = Xe_prev;
    kData[18] = Ye_prev;

    Double_t qx1,qy1,qz1,qx2,qy2,qz2;
    qx1=(*comb)[0]->Px();
    qy1=(*comb)[0]->Py();
    qz1=(*comb)[0]->Pz();
    qx2=((comb->Npart>1)?(*comb)[1]->Px():0);
    qy2=((comb->Npart>1)?(*comb)[1]->Py():0);
    qz2=((comb->Npart>1)?(*comb)[1]->Pz():0);
    Float_t E1=(*comb)[0]->E(); 
    Float_t E2=((comb->Npart>1)?(*comb)[1]->E():0); 
    kData[19] = qx1;
    kData[20] = qy1;
    kData[21] = qz1;
    kData[22] = qx2;
    kData[23] = qy2;
    kData[24] = qz2;
    kData[25] = E1;
    kData[26] = E2;
    kData[27] = E1;
    kData[28] = E2;
    //if (0.35<E1&&E1<1.2) kData[27] /= hcfm->GetBinContent(hcfm->FindBin(E1));
    //if (0.35<E2&&E2<1.2) kData[28] /= hcfm->GetBinContent(hcfm->FindBin(E2));
    kData[29] =(*comb)[0]->vx;
    kData[30] =(*comb)[0]->vy;
    kData[31] =((comb->Npart>1)?(*comb)[1]->vx:0);
    kData[32] =((comb->Npart>1)?(*comb)[1]->vy:0);
    kData[33] = TargType_prev;
    kData[34] = TargTypeO_prev;
    kData[35] = phiH;
    kData[36] = phiR;
    kData[37] = Mx2;
    kData[38] = xF;
    kData[39] = xF0;
    kData[40] = xF1;
    kData[41] = PlCM;
    kData[42] = Pl0CM;
    kData[43] = Pl1CM;
    kData[44] = E;
    kData[45] = xFm;
    kData[46] = xFm0;
    kData[47] = xFm1;
    kData[48] = (HELIC==-111)?helic_prev:HELIC;// if not given take it from file.
    kData[49] = theta_0*TMath::RadToDeg();
    kData[50] = theta_1*TMath::RadToDeg();
    kData[51] = cos_theta_P0cm;
    kData[52] = xFo;
    kData[53] = xF0o;
    kData[54] = xF1o;
    kData[55] = evnt_prev;
    Float_t dphi = (phiH-phiR);
    while (dphi<0)
      dphi+=360.;
    while (dphi>360)
      dphi-=360.;
    kData[56] = dphi;
    kData[57] = Ee_prev;

    kData[58] = phiR_cov;
    kData[59] = p0Tv.Mag2();
    kData[60] = p1Tv.Mag2();
    kData[61] = phi_pq;
    kData[62] = phTv.Mag2();
    kData[63] = revent_prev;

    kData[64] = etaCM_0;
    kData[65] = etaCM_1;
    kData[66] = etaBF_0plus;
    kData[67] = etaBF_1plus;
    kData[68] = etaBF_0minus;
    kData[69] = etaBF_1minus;

    kData[70] = etaBF0;
    kData[71] = etaBF1;
    
    kData[72] = phiH0;
    kData[73] = phiH1;

    kData[74] = -Pex_prev;
    kData[75] = -Pey_prev;
    kData[76] = kEbeam -Pez_prev;

    kData[77] = phiR_ha;
    kData[78] = Pl0m;
    kData[79] = Pl1m;
    kData[80] = phiR_covH;


    kData[81] = P0_dicm.E();
    kData[82] = P1_dicm.E();
    kData[83] = Nu_prev/kEbeam;
    kData[84] = acos(Pez_prev/sqrt(Pex_prev*Pex_prev + Pey_prev*Pey_prev + Pez_prev*Pez_prev))*TMath::RadToDeg();

    detData.npart = comb->Npart;
    for (int k =0 ;k<detData.npart;k++){
      detData.beta[k] = (*comb)[k]->beta;
      detData.m2b[k] =  (*comb)[k]->m2b;
      detData.vx[k] = (*comb)[k]->vx;
      detData.vy[k] = (*comb)[k]->vy;
      detData.vz[k] = (*comb)[k]->vz;
      detData.dcx[k] = (*comb)[k]->dcx;
      detData.dcy[k] = (*comb)[k]->dcy;
      detData.dcz[k] = (*comb)[k]->dcz;
    }

    /*  
    Double_t *W = new Double_t[4];
    Double_t *Wa = new Double_t[4];
    TMatrixD V(4,4);

    V(0,0) = 3.5*0.103/TMath::Sqrt(E);
    V(1,1) = 3.5*0.103/TMath::Sqrt(E);
    V(2,2) = 3.5*0.103/TMath::Sqrt(E);
    V(3,3) = 3.5*0.103/TMath::Sqrt(E);


    W[0]=Px;W[1]=Py;W[2]=Pz;W[3]=E;
    Wa[0]=Px;Wa[1]=Py;Wa[2]=Pz;Wa[3]=E;

    Double_t Chi2 =0;

    for (int k=0;k<1;k++){
          Chi2 = kinFit(W,Wa,V);
      
    }
    Double_t Phx_c = W[0], Phy_c= W[1], Phz_c= W[2], E_c = W[3], Z_c,Cospq_c, Pt2_c,M_c,M2_c,P2_c;
    P2_c =  Phx_c* Phx_c + Phy_c* Phy_c + Phz_c* Phz_c ; 
    M2_c = E_c*E_c -  P2_c;
    M_c=(M2_c>0)?TMath::Sqrt(E_c*E_c -  P2_c):M2_c;
    Cospq_c = ((kEbeam-Pez_prev)*Phz_c -Pex_prev*Phx_c - Pey_prev*Phy_c)/( sqrt((Q2_prev + Nu_prev*Nu_prev)*P2_c) );
    Pt2_c = P2_c*(1-Cospq_c*Cospq_c);
    Z_c = E_c/Nu_prev;

    delete[] W;
    delete[] Wa;

    kData[12] = M_c;
    kData[13] = Phx_c;
    kData[14] = Phy_c;
    kData[15] = Phz_c;
    kData[16] = Z_c;
    kData[17] = Cospq_c;
    kData[18] = Pt2_c; 
    kData[19] = Chi2;
    Double_t qx1,qy1,qz1,qx2,qy2,qz2;
    qx1=(*comb)[0]->Px();
    qy1=(*comb)[0]->Py();
    qz1=(*comb)[0]->Pz();
    qx2=(*comb)[1]->Px();
    qy2=(*comb)[1]->Py();
    qz2=(*comb)[1]->Pz();
  
    kData[20] = qx1;
    kData[21] = qy1;
    kData[22] = qz1;
    kData[23] = qx2;
    kData[24] = qy2;
    kData[25] = qz2;
*/
    /*
    ////// 3pi0-> 6a ///////////////////////////////////////////
    Double_t dist=1e7;
    Double_t M2_etaN=0;
    Double_t M2_eta0=2*0.135*0.135;
    Int_t order=4321;
    Int_t indexes[15]={4321, 4231, 4132};
    Int_t kk=0;
    for (int j=0;j<3;j++)
    {
      kk=indexes[j];
      M2_etaN=(*(*comb)[kk/int(1e3)-1] + *(*comb)[kk%int(1e3)/int(1e2)-1]).M2() \
            + (*(*comb)[kk%int(1e2)/int(1e1)-1] + *(*comb)[kk%int(1e1)-1]).M2();
      if (TMath::Abs(M2_etaN - M2_eta0)<dist)
      {
	order=indexes[j];
	dist=TMath::Abs(M2_etaN - M2_eta0);
      }
    }
    
    //Int_t ind[4]={order/int(1e3), order%int(1e3)/int(1e2), order%int(1e2)/int(1e1), order%int(1e1)};
    */
    ////////////////////////////////////////////////////////

    for (int k=0;k<kNSecondary;k++)
    {
      
      Double_t px= (*comb)[k]->Px();
      Double_t py= (*comb)[k]->Py();
      Double_t pz= (*comb)[k]->Pz();
      Double_t e= (*comb)[k]->E();
      
      new ((*P4Arr)[k]) TLorentzVector(px,py,pz,e);
    }


    /*
    Int_t ind[4]={1,2,3,4};
    for (int k=0;k<kNSecondary;k++)
    {
      
      Double_t px= (*comb)[ind[k]-1]->Px();
      Double_t py= (*comb)[ind[k]-1]->Py();
      Double_t pz= (*comb)[ind[k]-1]->Pz();
      Double_t e= (*comb)[ind[k]-1]->E();
      
      new ((*P4Arr)[k]) TLorentzVector(px,py,pz,e);
    }
    */


    return ttree->Fill();

  }
  Bool_t FidCheck(int pid)
  {
    Bool_t ret = false;
    Float_t pip_dcymin_r0 = 25,pip_dcymax_r0 = 135, pim_dcymin_r0 = 32,pim_dcymax_r0 = 145;
    Float_t dcy0,dcy_min,dcy_max,dcth_min,dcth_max;
    if (pid == 211)
    {
      dcy0 = 10;
      dcy_min = 25;
      dcy_max = 135;
      dcth_min = 64;
      dcth_max = 116;
    }
    else if (pid == -211)
    {
      
      /*// Sep-18 cooking
      dcy0 = 15;
      dcy_min = 32;
      dcy_max = 145;
      dcth_min = 65;
      dcth_max = 115;*/

      // March-19 cooking
      dcy0 = 25;
      dcy_min = 32;
      dcy_max = 145;
      dcth_min = 58;
      dcth_max = 122;
    }
    else return true;
    
    Float_t th_0 = TMath::ATan2(dcy_r0-dcy0,dcx_r0)*TMath::RadToDeg();
    
    if (dcy_min<dcy_r0&&dcy_r0<dcy_max
	&&dcth_min<th_0&&th_0<dcth_max)
      ret = true;
    
    return ret;

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
    for (int k=0;k<kSPid.size();k++)
    {
      if (fEMatch)
	minPart=minPart && (kSecondary[k].size() == kNSPid[kSPid[k]]);
      else
	minPart=minPart && (kSecondary[k].size() >= kNSPid[kSPid[k]]);

      //      std::cout<<__LINE__<<": minPart "<<minPart<<" : "<<kSecondary[k].size()<<std::endl;
      hSPid[k]->Fill(kSecondary[k].size());
    }
    return minPart;
  }

  bool checkMinPart(Combo *c)
  {
    bool minPart=true;
    for (int k=0;k<kSPid.size();k++)
    {
	minPart=minPart && (c->findPid(kSPid[k]) == kNSPid[kSPid[k]]);
    }
    return minPart;
  }


  void setElectVar()
  {
    Nu_prev = Nu;
    Q2_prev = Q2;
    Ze_prev = Ze;
    Xe_prev = Xe;
    Ye_prev = Ye;
    W_prev = W;
    Pex_prev = Pex;
    Pey_prev = Pey;
    Pez_prev = Pez;
    Ee_prev = Ee;
    TargType_prev = TargType;
    TargTypeO_prev = TargTypeO;
    helic_prev = helic;
    revent_prev = revent;
  }

  //int takeN(int N,int kspid, int pos=0,Particle p=Particle(),int count=0)
  //  int takeN(int N,int kspid, int pos=0,Combo *c= new Combo() ,int count=0)
  int takeN(int N,int kspid, int pos=0,Combo *c= 0 ,int count=0)
  {
    Combo *c_new;
    if (DEBUG) std::cout<<"### take N:  "<<N<<"##############"<<std::endl; 
    if (N<1) return -1;
    if (N!=1)
    {
      for (int k =pos;k<kSecondary[kspid].size()-N+1;k++)
      {
	if (c==0) c_new = new Combo();
	else c_new = new Combo(*c);

	c_new->addParticle(kSecondary[kspid][pos]);
	count=takeN(N-1,kspid,k+1,c_new,count);
      }
    }

    else
   {
     
      for (int k=pos;k<kSecondary[kspid].size();k++)
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
    return count;
  }
  int findSecondary()
  { 
    int count=0;
    for (int k =0;k<kSPid.size();k++)
    {
      //if(pid == kSPid[k]&&checkDCFidPhi())
      if(pid == kSPid[k])
      {
	kSecondary[k].push_back(new Particle(Px,Py,Pz,Ep,vx,vy,vz,pid,0,Beta,dcx_r0,dcy_r0,dcz_r0));
	//std::cout<<dcx_r0<<" / "<<dcy_r0<<" / "<<dcz_r0<<std::endl;
	count++;
      }
    }
    return count;
  }
  
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

  Float_t getE(Float_t Es=0, Float_t p=0, Float_t pid=211)
  {
    Float_t pip[4] = {-0.001114,0.1146,1.103e-5,-0.01259};
    Float_t pim[4] = {0.007386,0.09922,-0.001244,-0.01057};

    Float_t energy ;
    if (pid == 45)
      energy == sqrt(p*p + TMath::Power(1.8756,2) );
    else
      energy = sqrt(p*p + TMath::Power(TDatabasePDG::Instance()->GetParticle(pid)->Mass(),2) );
    //    Float_t sf = ((pid==211)?pip[0] + pip[1]/x + pip[2]*x + pip[3]/x/x:((pid==-211)?pim[0] + pim[1]/x + pim[2]*x + pim[3]/x/x:0)) ;
    // Float_t energy = Es/sf;

    return energy;
  }
  
  int getCombinations(TChain *t)
  {
    kSecondary=new std::vector<Particle*> [kSPid.size()];
    kCombo = new std::vector<Combo*> [kSPid.size()];

    t->GetEntry(0);
    setElectVar();
    evnt_prev=evnt;
    std::cout<<"processing...\n";
    std::cout.fill('.');
    for ( int i = 0; i < Ne; i++)
    {
      t->GetEntry(i);
      Ep=E;
      if (data_type<2)// not gsim
      {
	Ep = (pid==22)? (E/0.23):getE(E,P,pid);
	//  if (pid==22)
	//  correct_momentum();
	
      }
      if (evnt==evnt_prev)
      {
	//std::cout<<__LINE__<<" "<<findSecondary()<<std::endl;
	//if (data_type<2&&FidCheck(pid))
	if ((data_type==2)||FidCheck(11))// no fid
	{
	  Particle *p = new Particle(Px,Py,Pz,Ep,vx,vy,vz,pid,0,Beta,dcx_r0,dcy_r0,dcz_r0);
	 
	  push_bkgnd(p);
	  
	  findSecondary();
	}
      }
      else
      {
	pop_bkgnd();
	if (checkMinPart())
	{
	  int Npart = 1;

	  for (int k =0;k<kSPid.size();k++)
	  {
	    if (DEBUG) std::cout<<"############ Nsecondary of pid: "<<kSPid[k]<<" ::: "<<kSecondary[k].size()<<"###########"<<std::endl; 
	    takeN(kNSPid[kSPid[k]] ,k);

	    Npart*=kCombo[k].size();
	  }
	  if (DEBUG) std::cout<<"############ N candidates for primary: "<<Npart<<"####################\n##########################"<<std::endl;
	  for(int k=0;k<Npart;k++)
	  {
	    kPrimary = new Combo();
	    int div=1;
	    for (int l =0;l<kSPid.size();l++)
	    {
	      int size = kCombo[l].size();
	      //std::cout<<__LINE__<<": \n";
	      *kPrimary+=*kCombo[l][ (k/div)%size ];
	      //std::cout<<__LINE__<<": \n";

	      div*=size;
	    }
	    //  if (DEBUG) std::cout<<"############ Npart from primary: "<<kPrimary->Npart<<"####################"<<std::endl;
	    fill_elec();
  	    fill();
	    delete kPrimary;
	  }
	}
	clear();
	setElectVar();
	evnt_prev=evnt;
	//	if (data_type<2&&FidCheck(pid))

	if ((data_type==2)||FidCheck(11))// no fid
	{
	  Particle *p =new Particle(Px,Py,Pz,Ep,vx,vy,vz,pid,0,Beta,dcx_r0,dcy_r0,dcz_r0);
	  push_bkgnd(p);
	  findSecondary();
	}
      }
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
	for (i =0;i<kBkgnd.size();i++)
	{
	  if ( (kBkgnd[i]->findPid(p->pid)<kNSPid[p->pid] ) && (kBkgnd[i]->lastEvent!=evnt) && (kBkgnd[i]->isCompatible()) )
	  {
	    kBkgnd[i]->addParticle(p,evnt,kTRUE);
	    break;
	  }
	  
	}
	if (i==kBkgnd.size())
	{
	  Combo *c = new Combo();
	  c->addParticle(p,evnt);
	  kBkgnd.push_back(c);
	}
      }
      else 
      {
	Combo *c = new Combo();
	c->addParticle(p,evnt);
	kBkgnd.push_back(c);
      }
    }
    return 0;
  }

  bool pop_bkgnd()
  {
    bool ret = false;
    if (!kBkgnd.empty())
    {
     for (int i =0;i<kBkgnd.size();i++)
      {
	if (checkMinPart(kBkgnd[i]))
	{
	  //	  std::cout<<__FILE__<<"::"<<__LINE__<<std::endl;
	  fill(kBkgnd[i],kOutBkgnd);
	  delete kBkgnd[i];
	  kBkgnd.erase(kBkgnd.begin()+i);
	  ret = true;
	}
      } 
    } 
    return ret;
  }

  int getBkgnd(TChain *t)
  {
    
    t->GetEntry(0);
    setElectVar();
    evnt_prev=evnt;
    std::cout<<"processing...\n";
    std::cout.fill('.');
    for ( int i = 0; i < Ne; i++)
    {
      t->GetEntry(i);
      Ep = (pid==22)? (E/0.23):sqrt(P*P+ TMath::Power(TDatabasePDG::Instance()->GetParticle(pid)->Mass(),2));
      
      if (pid==22)
	correct_momentum();

      if (evnt==evnt_prev)
      {
	//std::cout<<__LINE__<<" "<<findSecondary()<<std::endl;
	Particle *p = new Particle(Px,Py,Pz,Ep,vx,vy,vz,pid);
	push_bkgnd(p);
      }
      else
      {
	pop_bkgnd();	
	setElectVar();
	evnt_prev=evnt;

	Particle *p = new Particle(Px,Py,Pz,Ep,vx,vy,vz,pid);
	push_bkgnd(p);
      }
      std::cout<<std::setw(15)<<float(i+1)/Ne*100<<" %"<<"\r";
      std::cout.flush();
    } 
    return 0;
  }
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
  bm->Start("get_eta");
  parseopt(argc,argv);
  //TFile * corrfile = new TFile("gammECorr.root","read");
  //hcfm = (TH1F*)corrfile->Get("hcfm");
  //char outdir[50];
  //strcpy(outdir,Form("test_eta_%sD_%s",st,tt));

  TChain *t = new TChain();

  TString PATH="";
  if (!INDIR.IsNull()) PATH = INDIR + "/*.root";
  else PATH = INFILE;

  std::cout<<"processing: "<<PATH<<std::endl;
  std::cout<<"data_type: "<<data_type<<std::endl;
  
  if (data_type==1) t->Add(PATH + "/ntuple_accept");
  else if (data_type==2) t->Add(PATH + "/ntuple_thrown");
  else 
    t->Add(PATH + "/ntuple_data");
  
  //  t->Add("/user/o/orsosa/osoto_ana/local/prune_simul.root/ntuple_accept"); //simrec test
  //t->Add("/data/atlas/users/orsosa/eg2_sim_pruned/C/pruned_simul_*.root/ntuple_accept"); //C  sim rec 
  //  t->Add("/data/atlas/users/orsosa/eg2_sim_pruned/C/pruned_simul_*.root/ntuple_thrown"); //C sim gsim 
  //  t->Add("/data/atlas/users/orsosa/eg2_sim_pruned/C_D/pruned*.root/ntuple_accept"); //C D sim rec
  //  t->Add("/data/atlas/users/orsosa/eg2_sim_pruned/C_D/pruned*.root/ntuple_thrown"); //C D sim rec

  //  t->Add("/data/atlas/users/orsosa/eg2_data_pruned/C-thickD2/pruned*.root/ntuple_data"); //data C and D
  //t->Add("/data/atlas/users/orsosa/eg2_data_pruned/Pb-thinD2/pruned*.root/ntuple_data"); //data Pb and D
  //t->Add("/data/atlas/users/orsosa/eg2_data_pruned/Fe-thickD2/pruned*.root/ntuple_data"); //data Fe and D


  //check_dir(outdir);


  t->SetBranchStatus("*",0);
  t->SetBranchStatus("E",1);
  t->SetBranchStatus("Ee",1);
  t->SetBranchStatus("P",1);
  t->SetBranchStatus("Px",1);
  t->SetBranchStatus("Py",1);
  t->SetBranchStatus("Pz",1);
  t->SetBranchStatus("evnt",1);
  t->SetBranchStatus("revent",1);
  t->SetBranchStatus("Ze",1);
  t->SetBranchStatus("Ye",1);
  t->SetBranchStatus("Xe",1);
  t->SetBranchStatus("TEc",1);
  t->SetBranchStatus("Q2",1);
  t->SetBranchStatus("Nu",1);
  t->SetBranchStatus("W",1);
  t->SetBranchStatus("Pex",1);
  t->SetBranchStatus("Pey",1);
  t->SetBranchStatus("Pez",1);
  t->SetBranchStatus("pid",1);
  t->SetBranchStatus("vxh",1);
  t->SetBranchStatus("vyh",1);
  t->SetBranchStatus("vzh",1);
  t->SetBranchStatus("TargType",1);
  t->SetBranchStatus("DCX",1);
  t->SetBranchStatus("DCY",1);
  t->SetBranchStatus("DCZ",1);
  t->SetBranchStatus("helic",1);
  t->SetBranchStatus("DCPx",1);
  t->SetBranchStatus("DCPy",1);
  t->SetBranchStatus("DCPz",1);
  t->SetBranchStatus("dcx_rot_0",1);
  t->SetBranchStatus("dcy_rot_0",1);
  t->SetBranchStatus("trajdczr0",1);
  t->SetBranchStatus("Beta",1);
  

  //  t->SetBranchStatus("TargTypeO",1);

  t->SetBranchAddress("E",&E);
  t->SetBranchAddress("Ee",&Ee);
  t->SetBranchAddress("P",&P);
  t->SetBranchAddress("Px",&Px);
  t->SetBranchAddress("Py",&Py);
  t->SetBranchAddress("Pz",&Pz);
  t->SetBranchAddress("Ze",&Ze);
  t->SetBranchAddress("Ye",&Ye);
  t->SetBranchAddress("Xe",&Xe);
  t->SetBranchAddress("evnt",&evnt);
  t->SetBranchAddress("revent",&revent);
  t->SetBranchAddress("TEc",&TEc);
  t->SetBranchAddress("Q2",&Q2);
  t->SetBranchAddress("Nu",&Nu);
  t->SetBranchAddress("W",&W);
  t->SetBranchAddress("Pex",&Pex);
  t->SetBranchAddress("Pey",&Pey);
  t->SetBranchAddress("Pez",&Pez);
  t->SetBranchAddress("pid",&pid);
  t->SetBranchAddress("vxh",&vx);
  t->SetBranchAddress("vyh",&vy);
  t->SetBranchAddress("vzh",&vz);
  t->SetBranchAddress("TargType",&TargType);
  t->SetBranchAddress("DCX",&DCX);
  t->SetBranchAddress("DCY",&DCY);
  t->SetBranchAddress("DCZ",&DCZ);
  t->SetBranchAddress("helic",&helic);
  t->SetBranchAddress("DCPx",&DCPx);
  t->SetBranchAddress("DCPy",&DCPy);
  t->SetBranchAddress("DCPz",&DCPz);
  t->SetBranchAddress("dcx_rot_0",&dcx_r0);
  t->SetBranchAddress("dcy_rot_0",&dcy_r0);
  t->SetBranchAddress("trajdczr0",&dcz_r0);
  t->SetBranchAddress("Beta",&Beta);

  //t->SetBranchAddress("TargTypeO",&TargTypeO);

  //  t->SetMaxEntryLoop(1e4);
  //  t->SetAlias("Eh",Form("(pid==22)*E/0.272 + (pid!=22)*(TMath::Sqrt(P**2 + %f**2)",TDatabasePDG::Instance()->GetParticle("pi-")->Mass() ));
  if (Ne==-1)  Ne = t->GetEntries();
  std::cout<<"Number of entries to be processed: "<<Ne<<std::endl;

  /*
  /// eta -> pi+ pi- a
  Reaction r("eta -> pi+ pi- a","test_pipapimOnly.root",true); 
  r.addPrimary("eta");
  r.addSecondary("pi+");
  r.addSecondary("gamma");
  r.addSecondary("pi-");
  */

  /*    
  // K+ -> pi+ pi+ pi-
  Reaction r("K+ -> pi- pi+ pi+","test_pimpippipOnlyC_sim.root",true);
  r.addPrimary("K+");
  r.addSecondary("pi-");
  r.addSecondary("pi+");
  r.addSecondary("pi+");
  */

  /*
  // K- -> pi- pi- pi+
  Reaction r("K- -> pi- pi- pi+","test_pimpimpip.root");
  r.addPrimary("K-");
  r.addSecondary("pi-");
  r.addSecondary("pi-");
  r.addSecondary("pi+");
  */ 

  /*  
  // eta -> pi0 pi0 pi0 -> 6a
  Reaction r("eta -> 6a","etaout_6a_all.root",false); //
  r.addPrimary("eta");
  r.addSecondary("gamma");
  r.addSecondary("gamma");
  r.addSecondary("gamma");
  r.addSecondary("gamma");
  r.addSecondary("gamma");
  r.addSecondary("gamma");
  */

  /*  
  // pi0 -> e- e+ a
  Reaction r("pi0 -> e- e+ a","pi0_CD_e-e+a.root",true);
  r.addPrimary("pi0");
  r.addSecondary("e-");
  r.addSecondary("e+");  
  r.addSecondary("gamma");
  */
  /*  
  
  //pi0 -> a a
  Reaction r("pi0 -> a a","outfiles/aa_all.root",false);
  r.addPrimary("pi0");
  r.addSecondary("gamma");
  r.addSecondary("gamma");
*/
  /*
  // Sigma0 -> L(p pi-) a
  Reaction r("Sigma0 -> p pi- a","outfiles/ppima_all.root",false);
  r.addPrimary("Xi0");
  r.addSecondary("proton");
  r.addSecondary("pi-");
  r.addSecondary("gamma");
  */

  
  /*  
  // Lambda0 -> p pi-
  Reaction r("Lambda0 -> p pi-","outfiles/ppim_all.root",false);
  r.addPrimary("Lambda0");
  r.addSecondary("proton");
  r.addSecondary("pi-");
  */


  // K0 -> pi+ pi-
  Reaction r("K0 -> pi+ pi-","outfiles/pippim_all.root",false);
  r.addPrimary("K0");
  r.addSecondary("pi+");
  r.addSecondary("pi-");

  
  /*
  // K0 -> pi+ pi-
  Reaction r("K0 -> pi+ pi-","outfiles/pippim_only.root",true);
  r.addPrimary("K0");
  r.addSecondary("pi+");
  r.addSecondary("pi-");
  */
  

  /*
  // Xi0 -> L(p pi-) pi0(a a)
  Reaction r("Xi0 -> p pi- a a","outfiles/ppimaa_all.root",false);
  r.addPrimary("Xi0");
  r.addSecondary("proton");
  r.addSecondary("pi-");
  r.addSecondary("gamma");
  r.addSecondary("gamma");
  */

  /*
  // Xi- -> L(p pi-) pi-
  Reaction r("Sigma0 -> p pi- pi-","outfiles/ppimpim_all.root",false);
  r.addPrimary("Xi-");
  r.addSecondary("proton");
  r.addSecondary("pi-");
  r.addSecondary("pi-");
  */
  
  /*    
  // Sigma0 -> L(p pi-) a
  Reaction r("Sigma0 -> p pi- a","outfiles/ppima_all.root",false);
  r.addPrimary("Xi0");
  r.addSecondary("proton");
  r.addSecondary("pi-");
  r.addSecondary("gamma");
  */
 
  
  /*  
  //eta -> a a
  Reaction r("eta -> a a","etaout_aa_only.root",true);
  r.addPrimary("eta");
  r.addSecondary("gamma");
  r.addSecondary("gamma");  
  
  */
  


  /*
  // w -> pi+ pi- a a
  Reaction r("w -> pi+ pi- a a","wout.root");
  r.addPrimary("omega");
  r.addSecondary("gamma");
  r.addSecondary("gamma");
  r.addSecondary("pi+");
  r.addSecondary("pi-");
  */

  /*        
  //eta -> a a pi+ pi-
  Reaction r("eta -> a a pi+ pi-","etaout_pippimaa_all_eta_exact.root",true);
  r.addPrimary("eta");
  r.addSecondary("gamma");
  r.addSecondary("gamma");
  r.addSecondary("pi+");
  r.addSecondary("pi-");
  */
  
  /*
  //eta -> a a
  Reaction r("eta -> a a","etaout_aa_only.root",true);
  r.addPrimary("eta");
  r.addSecondary("gamma");
  r.addSecondary("gamma");  
  */

  /*
  //eta -> a a
  Reaction r("eta -> a a","etaout_aa_all_bk.root",false);
  r.addPrimary("eta");
  r.addSecondary("gamma");
  r.addSecondary("gamma");  
  */

  /*
 //eta -> 3pi0 -> 6a 
  Reaction r("eta ->  3pi0 -> 6a","etaout6a.root",true);
  r.addPrimary("eta");
  r.addSecondary("gamma");
  r.addSecondary("gamma");
  r.addSecondary("gamma");
  r.addSecondary("gamma");
  r.addSecondary("gamma");
  r.addSecondary("gamma");
  */

  /*
  //eta -> 2pi0 -> 4a 
  Reaction r("eta ->  2pi0 -> 4a","etaout4a.root",true);
  r.addPrimary("eta");
  r.addSecondary("gamma");
  r.addSecondary("gamma");
  r.addSecondary("gamma");
  r.addSecondary("gamma");
  */


  /*  
  //corr -> pi+ pi+ pi-
  Reaction r("K0 -> pi+ pi+ pi-","pip_corr.root",true);
  r.addPrimary("K0");
  r.addSecondary("pi+");
  r.addSecondary("pi+"); 
  r.addSecondary("pi-");
  */


  r.getCombinations(t);
  r.store();
  std::cout<<"\n";
  r.kOutData->Print();
  r.kOutBkgnd->Print();
  //  corrfile->Close();
  bm->Show("get_pi0");
  return 0;  

}
