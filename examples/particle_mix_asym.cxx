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
#include "TRandom3.h"
#include <vector>
#include <iostream>
#include <iomanip>      // std::setw
#include <sys/stat.h>
#include <cstdarg>
#include <algorithm>
#include <map>
#include <string.h>
#include <time.h>

#include "particle_mix_asym.h"

bool DEBUG=false;
Float_t HELIC=-111;
TRandom3 *rndm;
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
//Float_t kMPi0=5.39609e-01;
//Float_t kSPi0=5.98542e-02;
bool GSIM=false;
int data_type=0;
Float_t kPt2,kEvent; 
TNtuple *tuple;
//TNtuple *tuple_sim;
TNtuple *tuplemb;
TNtuple *tuplePi0_gamma, *tupleGamma;
Float_t kEbeam=10.6,E,Ee,Ee_prev,Ep,P,Px,Py,Pz,evnt,evnt_prev,Ze,Ze_prev,Ye,Ye_prev,Xe,Xe_prev,TEc,Q2,Q2_prev,W,W_prev,Nu,Nu_prev,helic,helic_prev,Pex,Pex_prev,Pey,Pey_prev,Pez,Pez_prev,TargType,TargType_prev,TargTypeO=0,TargTypeO_prev=0,pid,vx,vy,vz,DCX,DCY,DCZ,ECX,ECY,ECZ,DCPx,DCPy,DCPz,Npip_rec_prev,Npim_rec_prev,Npip_mc_prev,Npim_mc_prev,rec_elec_prev,Npip_rec,Npim_rec,Npip_mc,Npim_mc,rec_elec;
long Ne = -1;
char st[3]= "C"; // solid target: C Fe Pb
char tt[3] = "C"; // cut on solid target or Deuterium : (st) or D.

Float_t kMprt=0.938272, kMntr =0.939565;
TClonesArray *P4Arr;

class Particle: public TLorentzVector
{
public:
  Float_t vx,vy,vz,pid,time;
  //TParticlePDG *info;
  Particle() : TLorentzVector(), vx(0),vy(0),vz(0),pid(0),time(0){}
  Particle(Float_t px,Float_t py, Float_t pz, Float_t e, Float_t x, Float_t y, Float_t z, Float_t pid=0, Float_t t=0): TLorentzVector(px,py,pz,e),vx(x),vy(y),vz(z),pid(pid),time(t){}
  Particle(TLorentzVector lv, Float_t x=0, Float_t y=0, Float_t z=0, Float_t pid=0, Float_t t=0): TLorentzVector(lv),vx(x),vy(y),vz(z),pid(pid),time(time){}
  Particle(Particle &p):vx(p.vx),vy(p.vy),vz(p.vz),pid(p.pid),time(p.time) {SetVect(p.Vect()); SetT(p.T());}
  inline Double_t P2() const {return P()*P();}
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

  TTree *kOutMCData, *kOutBkgnd, *kOutRSData;

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
    kOutRSData->Write("",TObject::kOverwrite);
    
    return kOutMCData->Write("",TObject::kOverwrite);
  }
  ~Reaction()
  {
    clear();
    delete[] kSecondary;
    delete[] kCombo;
    delete kOutMCData;
    delete kOutRSData;
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
    kNSecondary = 0;
    kSecondary = 0;
    kCombo = 0;
    kOutFile = new TFile(filename,"recreate");
    TString varlist="M:Phx:Phy:Phz:Nu:Q2:Z:Cospq:Pt2:Event:M2_01:M2_02:vzec:z1:z2:z3:W:vxec:vyec:qx1:qy1:qz1:qx2:qy2:qz2:E1:E2:E1c:E2c:x1:y1:x2:y2:TargType:TargTypeO:phiH:phiR:Mx2:xF:xF0:xF1:plcm:plcm0:plcm1:Eh:xFm:xFm0:xFm1:helic:theta0:theta1:cos_theta_P0cm:xFo:xF0o:xF1o:event:asym:Npip_rec:Npim_rec:Npip_mc:Npim_mc:rec_elec";
    kElecData = new TNtuple("ElecData",Form("%s",name),"Nu:Q2:Event:vze:Ee:Pex:Pey:Pez:W");

    kData = new Float_t[varlist.CountChar(':')+1];
    keData = new Float_t[kElecData->GetNvar()];
    P4Arr = new TClonesArray("TLorentzVector");
    kOutMCData = new TTree("outdata_mc",Form("%s",name));
    kOutMCData->Branch("P4",&P4Arr,6);
    kOutMCData->Branch("primary",kData,varlist);

    kOutRSData = new TTree("outdata_rs",Form("%s",name));
    kOutRSData->Branch("P4",&P4Arr,6);
    kOutRSData->Branch("primary",kData,varlist);
    
    kOutBkgnd = kOutMCData->CloneTree(0);
    kOutBkgnd->SetName("outbkgnd");
    
    return 0;
  }

  int fill(TTree *ttree,Int_t &asym,Bool_t is_mc = false)
  {
    fill(kPrimary,ttree,asym, is_mc);
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


  int fill(Combo *comb,TTree *ttree,Int_t &asym, Bool_t is_mc)
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
    Float_t phiH = phi_pq;
    
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
    //// /////////////////////////////////

    Float_t Plt = Ph.Vect()*q_lv.Vect().Unit();
    Float_t Pl0m = P0.Vect()*q_lv.Vect().Unit();
    Float_t Pl1m = P1.Vect()*q_lv.Vect().Unit();
    
 
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
    kData[35] = phi_pq;
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
   

    ///// test to include the asymmetry. ////////
    if (is_mc)
    {
      Float_t uval = rndm->Uniform(); // 0-1 uniform distribution.
      Float_t aval = asym_n(E/Nu_prev,phiR,HELIC);
      if (uval < aval)
	asym = 1;
      else
	asym = 0;
    }
      
    kData[56] = asym;
    kData[57] = Npip_rec_prev;
    kData[58] = Npim_rec_prev;
    kData[59] = Npip_mc_prev;
    kData[60] = Npim_mc_prev;
    kData[61] = rec_elec_prev;
    



    
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
  //////////////////////////

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
    Nu_prev=Nu;
    Q2_prev=Q2;
    Ze_prev=Ze;
    Xe_prev=Xe;
    Ye_prev=Ye;
    W_prev = W;
    Pex_prev=Pex;
    Pey_prev=Pey;
    Pez_prev=Pez;
    Ee_prev = Ee;
    TargType_prev = TargType;
    TargTypeO_prev = TargTypeO;
    helic_prev = helic;
    Npip_rec_prev = Npip_rec;
    Npim_rec_prev = Npim_rec;
    Npip_mc_prev = Npip_mc;
    Npim_mc_prev = Npim_mc;
    rec_elec_prev = rec_elec;
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
	kSecondary[k].push_back(new Particle(Px,Py,Pz,Ep,vx,vy,vz,pid));
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
    Float_t x=sqrt(p*p+0.139*0.139);
    Float_t sf = ((pid==211)?pip[0] + pip[1]/x + pip[2]*x + pip[3]/x/x:((pid==-211)?pim[0] + pim[1]/x + pim[2]*x + pim[3]/x/x:0)) ;
    // Float_t energy = Es/sf;
    Float_t energy = x;
   

    return energy;
  }
  
  int next_event(TChain *t,Int_t start_ind, Int_t &end_ind)
  {
    Int_t Ne = t->GetEntries();
    t->GetEntry(start_ind);
    evnt_prev=evnt;
    for ( int i = start_ind; i < Ne; i++)
    {
      t->GetEntry(i);
      if  (evnt_prev!=evnt)
      {
	end_ind = i;
	return evnt;
      }
    }
    return -1;
    
  }
  int getCombinations(TChain *t,Int_t start_ind, Int_t &end_ind,Int_t &asym, Bool_t is_mc)
  {
    Int_t Ne = t->GetEntries();
    if (!kSecondary)
      kSecondary=new std::vector<Particle*> [kSPid.size()];
    if (!kCombo)
      kCombo = new std::vector<Combo*> [kSPid.size()];

    
    t->GetEntry(start_ind);
    setElectVar();
    evnt_prev=evnt;
    for ( int i = start_ind; i < Ne; i++)
    {
      t->GetEntry(i);
      Ep=E;
      if (!is_mc)// not gsim
      {
	Ep = (pid==22)? (E/0.23):( (pid==211 || pid==-211)?(getE(E,P,pid) ):E);
 
	//  if (pid==22)
	//  correct_momentum();
	
      }

      if (evnt==evnt_prev)
      {
	//std::cout<<__LINE__<<" "<<findSecondary()<<std::endl;
	Particle *p = new Particle(Px,Py,Pz,Ep,vx,vy,vz,pid);
	
	push_bkgnd(p);
	
	findSecondary();
      }
      else
      {
	end_ind = i;
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
  	    if (is_mc) fill(kOutMCData,asym,is_mc);
	    else fill(kOutRSData,asym,is_mc);
	    
	    delete kPrimary;
	  }
	}
	clear();
	/*	setElectVar();
	evnt_prev=evnt;
	Particle *p =new Particle(Px,Py,Pz,Ep,vx,vy,vz,pid);
	push_bkgnd(p);
	findSecondary();*/
	return 0;
      }
      
    } 
    return 1;
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
	  int asym_bg;
	  fill(kBkgnd[i],kOutBkgnd,asym_bg,false);
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
      Ep = (pid==22)? (E/0.23):sqrt(P*P+ TMath::Power(TDatabasePDG::Instance()->GetParticle("pi-")->Mass(),2));
      
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
  time_t seconds;
  seconds = time(NULL);
  rndm = new TRandom3(seconds);
  TChain *t_mc = new TChain();
  TChain *t_rs = new TChain();
  t_rs->Add("data.root/ntuple_accept");
  t_mc->Add("data.root/ntuple_thrown");
 
  /////// MC Tuple ////////
  t_mc->SetBranchStatus("*",0);
  t_mc->SetBranchStatus("E",1);
  t_mc->SetBranchStatus("Ee",1);
  t_mc->SetBranchStatus("P",1);
  t_mc->SetBranchStatus("Px",1);
  t_mc->SetBranchStatus("Py",1);
  t_mc->SetBranchStatus("Pz",1);
  t_mc->SetBranchStatus("evnt",1);
  t_mc->SetBranchStatus("Ze",1);
  t_mc->SetBranchStatus("Ye",1);
  t_mc->SetBranchStatus("Xe",1);
  t_mc->SetBranchStatus("TEc",1);
  t_mc->SetBranchStatus("Q2",1);
  t_mc->SetBranchStatus("Nu",1);
  t_mc->SetBranchStatus("W",1);
  t_mc->SetBranchStatus("Pex",1);
  t_mc->SetBranchStatus("Pey",1);
  t_mc->SetBranchStatus("Pez",1);
  t_mc->SetBranchStatus("pid",1);
  t_mc->SetBranchStatus("vxh",1);
  t_mc->SetBranchStatus("vyh",1);
  t_mc->SetBranchStatus("vzh",1);
  t_mc->SetBranchStatus("TargType",1);
  t_mc->SetBranchStatus("DCX",1);
  t_mc->SetBranchStatus("DCY",1);
  t_mc->SetBranchStatus("DCZ",1);
  t_mc->SetBranchStatus("helic",1);
  t_mc->SetBranchStatus("DCPx",1);
  t_mc->SetBranchStatus("DCPy",1);
  t_mc->SetBranchStatus("DCPz",1);

  t_mc->SetBranchStatus("Npip_rec",1);
  t_mc->SetBranchStatus("Npim_rec",1);
  t_mc->SetBranchStatus("Npip_mc",1);
  t_mc->SetBranchStatus("Npim_mc",1);
  t_mc->SetBranchStatus("rec_elec",1);
  
  t_mc->SetBranchAddress("E",&E);
  t_mc->SetBranchAddress("Ee",&Ee);
  t_mc->SetBranchAddress("P",&P);
  t_mc->SetBranchAddress("Px",&Px);
  t_mc->SetBranchAddress("Py",&Py);
  t_mc->SetBranchAddress("Pz",&Pz);
  t_mc->SetBranchAddress("Ze",&Ze);
  t_mc->SetBranchAddress("Ye",&Ye);
  t_mc->SetBranchAddress("Xe",&Xe);
  t_mc->SetBranchAddress("evnt",&evnt);
  t_mc->SetBranchAddress("TEc",&TEc);
  t_mc->SetBranchAddress("Q2",&Q2);
  t_mc->SetBranchAddress("Nu",&Nu);
  t_mc->SetBranchAddress("W",&W);
  t_mc->SetBranchAddress("Pex",&Pex);
  t_mc->SetBranchAddress("Pey",&Pey);
  t_mc->SetBranchAddress("Pez",&Pez);
  t_mc->SetBranchAddress("pid",&pid);
  t_mc->SetBranchAddress("vxh",&vx);
  t_mc->SetBranchAddress("vyh",&vy);
  t_mc->SetBranchAddress("vzh",&vz);
  t_mc->SetBranchAddress("TargType",&TargType);
  t_mc->SetBranchAddress("DCX",&DCX);
  t_mc->SetBranchAddress("DCY",&DCY);
  t_mc->SetBranchAddress("DCZ",&DCZ);
  t_mc->SetBranchAddress("helic",&helic);
  t_mc->SetBranchAddress("DCPx",&DCPx);
  t_mc->SetBranchAddress("DCPy",&DCPy);
  t_mc->SetBranchAddress("DCPz",&DCPz);

  t_mc->SetBranchAddress("Npip_rec",&Npip_rec);
  t_mc->SetBranchAddress("Npim_rec",&Npim_rec);
  t_mc->SetBranchAddress("Npip_mc",&Npip_mc);
  t_mc->SetBranchAddress("Npim_mc",&Npim_mc);
  t_mc->SetBranchAddress("rec_elec",&rec_elec);
  
  ///////// end mc tuple //////////
  
  //////// RS Tuple /////
  t_rs->SetBranchStatus("*",0);
  t_rs->SetBranchStatus("E",1);
  t_rs->SetBranchStatus("Ee",1);
  t_rs->SetBranchStatus("P",1);
  t_rs->SetBranchStatus("Px",1);
  t_rs->SetBranchStatus("Py",1);
  t_rs->SetBranchStatus("Pz",1);
  t_rs->SetBranchStatus("evnt",1);
  t_rs->SetBranchStatus("Ze",1);
  t_rs->SetBranchStatus("Ye",1);
  t_rs->SetBranchStatus("Xe",1);
  t_rs->SetBranchStatus("TEc",1);
  t_rs->SetBranchStatus("Q2",1);
  t_rs->SetBranchStatus("Nu",1);
  t_rs->SetBranchStatus("W",1);
  t_rs->SetBranchStatus("Pex",1);
  t_rs->SetBranchStatus("Pey",1);
  t_rs->SetBranchStatus("Pez",1);
  t_rs->SetBranchStatus("pid",1);
  t_rs->SetBranchStatus("vxh",1);
  t_rs->SetBranchStatus("vyh",1);
  t_rs->SetBranchStatus("vzh",1);
  t_rs->SetBranchStatus("TargType",1);
  t_rs->SetBranchStatus("DCX",1);
  t_rs->SetBranchStatus("DCY",1);
  t_rs->SetBranchStatus("DCZ",1);
  t_rs->SetBranchStatus("helic",1);
  t_rs->SetBranchStatus("DCPx",1);
  t_rs->SetBranchStatus("DCPy",1);
  t_rs->SetBranchStatus("DCPz",1);
  t_rs->SetBranchStatus("Npip_rec",1);
  t_rs->SetBranchStatus("Npim_rec",1);
  t_rs->SetBranchStatus("Npip_mc",1);
  t_rs->SetBranchStatus("Npim_mc",1);
  t_rs->SetBranchStatus("rec_elec",1);

  
  t_rs->SetBranchAddress("E",&E);
  t_rs->SetBranchAddress("Ee",&Ee);
  t_rs->SetBranchAddress("P",&P);
  t_rs->SetBranchAddress("Px",&Px);
  t_rs->SetBranchAddress("Py",&Py);
  t_rs->SetBranchAddress("Pz",&Pz);
  t_rs->SetBranchAddress("Ze",&Ze);
  t_rs->SetBranchAddress("Ye",&Ye);
  t_rs->SetBranchAddress("Xe",&Xe);
  t_rs->SetBranchAddress("evnt",&evnt);
  t_rs->SetBranchAddress("TEc",&TEc);
  t_rs->SetBranchAddress("Q2",&Q2);
  t_rs->SetBranchAddress("Nu",&Nu);
  t_rs->SetBranchAddress("W",&W);
  t_rs->SetBranchAddress("Pex",&Pex);
  t_rs->SetBranchAddress("Pey",&Pey);
  t_rs->SetBranchAddress("Pez",&Pez);
  t_rs->SetBranchAddress("pid",&pid);
  t_rs->SetBranchAddress("vxh",&vx);
  t_rs->SetBranchAddress("vyh",&vy);
  t_rs->SetBranchAddress("vzh",&vz);
  t_rs->SetBranchAddress("TargType",&TargType);
  t_rs->SetBranchAddress("DCX",&DCX);
  t_rs->SetBranchAddress("DCY",&DCY);
  t_rs->SetBranchAddress("DCZ",&DCZ);
  t_rs->SetBranchAddress("helic",&helic);
  t_rs->SetBranchAddress("DCPx",&DCPx);
  t_rs->SetBranchAddress("DCPy",&DCPy);
  t_rs->SetBranchAddress("DCPz",&DCPz);
  t_rs->SetBranchAddress("Npip_rec",&Npip_rec);
  t_rs->SetBranchAddress("Npim_rec",&Npim_rec);
  t_rs->SetBranchAddress("Npip_mc",&Npip_mc);
  t_rs->SetBranchAddress("Npim_mc",&Npim_mc);
  t_rs->SetBranchAddress("rec_elec",&rec_elec);

  /////////// end acc tuple ////


  if (Ne==-1)  Ne = t_mc->GetEntries();
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

  // pi0 -> a a
  Reaction r("pi0 -> a a","pi0_CD_2aonly_c.root");
  r.addPrimary("pi0");
  r.addSecondary("gamma");
  r.addSecondary("gamma");
  */

  /*
  // pi0 -> a a
  Reaction r("pi0 -> a a","pi0out.root");
  r.addPrimary("pi0");
  r.addSecondary("gamma");
  r.addSecondary("gamma");
  */

  
  // K0 -> pi+ pi-
  Reaction r("K0 -> pi+ pi-","pippim_exact_asym.root",true);
  r.addPrimary("K0");
  r.addSecondary("pi+");
  r.addSecondary("pi-");
  

  
  /*  
  //eta -> a a
  Reaction r("eta -> a a","etaout_aa_only.root",true);
  r.addPrimary("eta");
  r.addSecondary("gamma");
  r.addSecondary("gamma");  
  */

  
  
  /*
    // K0 -> pi+ pi-
  Reaction r("K0 -> pi+ pi-","pippim_all.root",false);
  r.addPrimary("K0");
  r.addSecondary("pi+");
  r.addSecondary("pi-");
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
  Int_t end_ind_rs=0, start_ind_rs = 0, asym, end_ind_mc, mc_ev, rs_ev;
  std::cout<<"processing...\n";
  std::cout.fill('.');
  int ret;
  for (int k = 0; k<Ne;k++)
  {

    HELIC = rndm->Uniform()>0.5?1:0;
    t_mc->GetEntry(k);
    mc_ev = evnt;
    t_rs->GetEntry(start_ind_rs);
    rs_ev = evnt;

    if (mc_ev > rs_ev) // no MC on record.
    {
      rs_ev = r.next_event(t_rs, start_ind_rs, end_ind_rs);
      start_ind_rs = end_ind_rs;
    }

    
    ret = r.getCombinations(t_mc, k, end_ind_mc, asym, true);
    if (ret) break;

    if (rs_ev==mc_ev)
    {
      //      std::cout<<__LINE__<<std::endl<<std::endl;
      r.getCombinations(t_rs, start_ind_rs, end_ind_rs, asym, false);
      
    }
    start_ind_rs = end_ind_rs;
    k = end_ind_mc;
    
    std::cout<<std::setw(15)<<float(k+1)/Ne*100<<" %"<<"\r";
    std::cout.flush();

  }


  r.store();
  std::cout<<"\n";
  r.kOutMCData->Print();
  r.kOutRSData->Print();
  //  corrfile->Close();
  bm->Show("get_pi0");
  return 0;  

}
