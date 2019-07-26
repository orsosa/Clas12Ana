#ifndef TIDENTIFICATOR_H
#define TIDENTIFICATOR_H

#include "TObject.h"
#include "reader.h"
#include "TString.h"
#include "TVector3.h"

#define  MAXNFILES 2000

class TIdentificatorCLAS12 {
public:
  TIdentificatorCLAS12();
  explicit TIdentificatorCLAS12(hipo::reader *reader = 0,Double_t beamE=10.6);
  explicit TIdentificatorCLAS12(TString fname, Double_t beamE=10.6, Bool_t mcf=false);
  ~TIdentificatorCLAS12();
    
  Bool_t Next(); // explore the hipo file opened

    // HEADER bank
  Float_t NEvent();                             // inline

    // REC::Particle
    Double_t Beta(Int_t k);         // inline
    Double_t Charge(Int_t k, Bool_t kind = 0);    // inline
    Double_t Px(Int_t k, Bool_t kind = 0);        // inline
    Double_t Py(Int_t k, Bool_t kind = 0);        // inline
    Double_t Pz(Int_t k, Bool_t kind = 0);        // inline
    Double_t X(Int_t k, Bool_t kind = 0);         // inline
    Double_t Y(Int_t k, Bool_t kind = 0);         // inline
    Double_t Z(Int_t k, Bool_t kind = 0);         // inline
    Int_t Helic();                       // inline
    Int_t HelicRaw();                       // inline
    Int_t HelicOnline();                       // inline
    Int_t HelicOnlineRaw();                       // inline
    Int_t HelicFlip();                       // inline
    Int_t HelicFlipRaw();                       // inline
    Int_t HelicFlipEvent();                       // inline
    Int_t Event();                       // inline
    Int_t RunN();                       // inline
    Float_t Torus();
    Float_t Solenoid();
    Float_t STTime(Int_t k=0);
    Float_t RFTime(Int_t k=0);
    Int_t GetNRows();                             // inline
    Int_t GetMCNRows();                           // inline
    Float_t MCMass(Int_t k=0);                    // inline
    Float_t LundType(Int_t k=0);                  // inline
    Int_t GetNPart(Int_t pid = 211, Bool_t kind = 0);  // inline
    
    Int_t StatCC(Int_t k);                        // inline
    Int_t StatSC(Int_t k);                        // inline
    Int_t StatDC(Int_t k);                        // inline
    Int_t StatEC(Int_t k);                        // inline
    Int_t Status(Int_t k);                        // inline


    // CC
    Double_t Nphe(Int_t k);                       // inline
    Double_t NpheLTCC(Int_t k);                   // inline
    Double_t NpheHTCC(Int_t k);                   // inline
    Double_t Chi2pid(Int_t k);                    // inline
    Double_t Chi2CC(Int_t k);                     // inline

    Double_t CCStatus(Int_t k);                   // inline

    // DCPB
    Double_t DCChi2(Int_t k=0);                     // inline

    // ECPB
    Double_t Etot(Int_t k,Bool_t kind=0);         // inline
    Double_t Ein(Int_t k);                        // inline
    Double_t Eout(Int_t k);                       // inline
    Double_t Epcal(Int_t k);                      // inline
    Double_t ECStatus(Int_t k);                   // inline
    Float_t XEC(Int_t k);                         // inline 
    Float_t YEC(Int_t k);                         // inline 
    Float_t ZEC(Int_t k);                         // inline 
    Float_t TimeEC(Int_t k);                      // inline 
    Float_t PathEC(Int_t k);                      // inline 

    Float_t LU_PCAL(Int_t k=0);                   // inline 
    Float_t LV_PCAL(Int_t k=0);                   // inline 
    Float_t LW_PCAL(Int_t k=0);                   // inline 

    Float_t LU_ECIN(Int_t k=0);                   // inline 
    Float_t LV_ECIN(Int_t k=0);                   // inline 
    Float_t LW_ECIN(Int_t k=0);                   // inline 

    Float_t LU_ECOUT(Int_t k=0);                  // inline 
    Float_t LV_ECOUT(Int_t k=0);                  // inline 
    Float_t LW_ECOUT(Int_t k=0);                  // inline 


    Float_t HX_PCAL(Int_t k=0);                   // inline 
    Float_t HY_PCAL(Int_t k=0);                   // inline 
    Float_t HZ_PCAL(Int_t k=0);                   // inline 

    Float_t HX_ECIN(Int_t k=0);                         // inline 
    Float_t HY_ECIN(Int_t k=0);                         // inline 
    Float_t HZ_ECIN(Int_t k=0);                         // inline 

    Float_t HX_ECOUT(Int_t k=0);                         // inline 
    Float_t HY_ECOUT(Int_t k=0);                         // inline 
    Float_t HZ_ECOUT(Int_t k=0);                         // inline 
    
    Float_t Px_DC(Int_t k=0);                         // inline 
    Float_t Py_DC(Int_t k=0);                         // inline 
    Float_t Pz_DC(Int_t k=0);                         // inline 

    Float_t VX_DC(Int_t k=0);                         // inline 
    Float_t VY_DC(Int_t k=0);                         // inline 
    Float_t VZ_DC(Int_t k=0);                         // inline 

    Float_t RICH_HAD_X(Int_t k=0); //inline
    Float_t RICH_HAD_Y(Int_t k=0); //inline
    Float_t RICH_HAD_Z(Int_t k=0); //inline
    Float_t RICH_HAD_T(Int_t k=0); //inline

    Float_t RICH_HAD_WX(Int_t k=0); //inline
    Float_t RICH_HAD_WY(Int_t k=0); //inline
    Float_t RICH_HAD_WZ(Int_t k=0); //inline
    Float_t RICH_HAD_WT(Int_t k=0); //inline

    Float_t RICH_CLUSTER_X(Int_t k=0); //inline
    Float_t RICH_CLUSTER_Y(Int_t k=0); //inline
    Float_t RICH_CLUSTER_Z(Int_t k=0); //inline
    Float_t RICH_CLUSTER_T(Int_t k=0); //inline

    Float_t RICH_CLUSTER_WX(Int_t k=0); //inline
    Float_t RICH_CLUSTER_WY(Int_t k=0); //inline
    Float_t RICH_CLUSTER_WZ(Int_t k=0); //inline
    Float_t RICH_CLUSTER_WT(Int_t k=0); //inline

    Float_t RICH_RR_X(Int_t k=0); //inline
    Float_t RICH_RR_Y(Int_t k=0); //inline
    Float_t RICH_RR_Z(Int_t k=0); //inline
    Float_t RICH_RR_T(Int_t k=0); //inline

    Float_t RICH_RR_HX(Int_t k=0); //inline
    Float_t RICH_RR_HY(Int_t k=0); //inline
    Float_t RICH_RR_HZ(Int_t k=0); //inline
    Float_t RICH_RR_HT(Int_t k=0); //inline

    Float_t RICH_PMT(Int_t k=0); //inline
    Float_t RICH_ANODE(Int_t k=0); //inline
    Float_t RICH_HITS_X(Int_t k=0); //inline
    Float_t RICH_HITS_Y(Int_t k=0); //inline
    Float_t RICH_HITS_Z(Int_t k=0); //inline

    
    // SCPB
    Double_t PathSC(Int_t k);                     // inline
    Double_t TimeSC(Int_t k);                     // inline
    Double_t EdepSC(Int_t k);                     // inline
    Double_t SCStatus(Int_t k);                   // inline
    

    // Derived observables
    Double_t Momentum(Int_t k, Bool_t kind = 0);
    Double_t Mass2(Int_t k, Bool_t kind = 0);
    Double_t ThetaLab(Int_t k, Bool_t kind = 0);
    Double_t PhiLab(Int_t k, Bool_t kind = 0);
    Double_t ThetaVirtLab(Bool_t kind = 0);
    Double_t PhiVirtLab(Bool_t kind = 0);
    Double_t ThetaPQ(Int_t k, Bool_t kind = 0);
    Double_t PhiPQ(Int_t k, Bool_t kind = 0);
    Double_t CosThetaPQ(Int_t k, Bool_t kind = 0);
    Double_t PTrans2PQ(Int_t k, Bool_t kind = 0);
    Double_t PLong2PQ(Int_t k, Bool_t kind = 0);
    Int_t Sector(Int_t k, Bool_t kind = 0);
    Int_t SectorLTCC(Int_t k=0);
    Int_t SectorHTCC(Int_t k=0);
    Int_t SectorECAL(Int_t k=0);
    Int_t SectorDC(Int_t k=0);
    Int_t TrajDetId(Int_t k=0);
    Float_t TrajX(Int_t k=0,Int_t sl=0);
    Float_t TrajY(Int_t k=0,Int_t sl=0);
    Float_t TrajZ(Int_t k=0,Int_t sl=0);

    Float_t TrajCX(Int_t k=0,Int_t reg=0);
    Float_t TrajCY(Int_t k=0,Int_t reg=0);
    Float_t TrajCZ(Int_t k=0,Int_t reg=0);

    Float_t TrajDCX(Int_t k=0,Int_t reg=0);
    Float_t TrajDCY(Int_t k=0,Int_t reg=0);
    Float_t TrajDCZ(Int_t k=0,Int_t reg=0);
    
    Float_t TrajDetIdX(Int_t k=0, TString dname="", TString layer="");
    Float_t TrajDetIdY(Int_t k=0, TString dname="", TString layer="");
    Float_t TrajDetIdZ(Int_t k=0, TString dname="", TString layer="");

    Float_t TrajDetId_CX(Int_t k=0, TString dname="", TString layer="");
    Float_t TrajDetId_CY(Int_t k=0, TString dname="", TString layer="");
    Float_t TrajDetId_CZ(Int_t k=0, TString dname="", TString layer="");
    
    Float_t TrajECINX(Int_t k=0);
    Float_t TrajECINY(Int_t k=0);
    Float_t TrajECINZ(Int_t k=0);

    Float_t TrajECOUTX(Int_t k=0);
    Float_t TrajECOUTY(Int_t k=0);
    Float_t TrajECOUTZ(Int_t k=0);
    
    Float_t TrajPCALX(Int_t k=0);
    Float_t TrajPCALY(Int_t k=0);
    Float_t TrajPCALZ(Int_t k=0);

    Float_t TrajLTCCX(Int_t k=0);
    Float_t TrajLTCCY(Int_t k=0);
    Float_t TrajLTCCZ(Int_t k=0);

    Float_t TrajFTOF1AX(Int_t k=0);
    Float_t TrajFTOF1AY(Int_t k=0);
    Float_t TrajFTOF1AZ(Int_t k=0);

    Float_t TrajFTOF1BX(Int_t k=0);
    Float_t TrajFTOF1BY(Int_t k=0);
    Float_t TrajFTOF1BZ(Int_t k=0);

    Float_t TrajFTOF2X(Int_t k=0);
    Float_t TrajFTOF2Y(Int_t k=0);
    Float_t TrajFTOF2Z(Int_t k=0);

    Float_t Traj101X(Int_t k=0);
    Float_t Traj101Y(Int_t k=0);
    Float_t Traj101Z(Int_t k=0);

    Float_t Traj102X(Int_t k=0);
    Float_t Traj102Y(Int_t k=0);
    Float_t Traj102Z(Int_t k=0);

    Float_t TrajHTCCX(Int_t k=0);
    Float_t TrajHTCCY(Int_t k=0);
    Float_t TrajHTCCZ(Int_t k=0);

    Float_t TrajFMT1X(Int_t k=0);
    Float_t TrajFMT1Y(Int_t k=0);
    Float_t TrajFMT1Z(Int_t k=0);

    Float_t TrajFMT2X(Int_t k=0);
    Float_t TrajFMT2Y(Int_t k=0);
    Float_t TrajFMT2Z(Int_t k=0);

    Float_t TrajFMT3X(Int_t k=0);
    Float_t TrajFMT3Y(Int_t k=0);
    Float_t TrajFMT3Z(Int_t k=0);

    Float_t TrajFMT4X(Int_t k=0);
    Float_t TrajFMT4Y(Int_t k=0);
    Float_t TrajFMT4Z(Int_t k=0);

    Float_t TrajFMT5X(Int_t k=0);
    Float_t TrajFMT5Y(Int_t k=0);
    Float_t TrajFMT5Z(Int_t k=0);

    Float_t TrajFMT6X(Int_t k=0);
    Float_t TrajFMT6Y(Int_t k=0);
    Float_t TrajFMT6Z(Int_t k=0);

    ///CXCYCZ
    Float_t TrajECIN_CX(Int_t k=0);
    Float_t TrajECIN_CY(Int_t k=0);
    Float_t TrajECIN_CZ(Int_t k=0);

    Float_t TrajECOUT_CX(Int_t k=0);
    Float_t TrajECOUT_CY(Int_t k=0);
    Float_t TrajECOUT_CZ(Int_t k=0);
    
    Float_t TrajPCAL_CX(Int_t k=0);
    Float_t TrajPCAL_CY(Int_t k=0);
    Float_t TrajPCAL_CZ(Int_t k=0);

    Float_t TrajLTCC_CX(Int_t k=0);
    Float_t TrajLTCC_CY(Int_t k=0);
    Float_t TrajLTCC_CZ(Int_t k=0);

    Float_t TrajFTOF1A_CX(Int_t k=0);
    Float_t TrajFTOF1A_CY(Int_t k=0);
    Float_t TrajFTOF1A_CZ(Int_t k=0);

    Float_t TrajFTOF1B_CX(Int_t k=0);
    Float_t TrajFTOF1B_CY(Int_t k=0);
    Float_t TrajFTOF1B_CZ(Int_t k=0);

    Float_t TrajFTOF2_CX(Int_t k=0);
    Float_t TrajFTOF2_CY(Int_t k=0);
    Float_t TrajFTOF2_CZ(Int_t k=0);

    Float_t Traj101_CX(Int_t k=0);
    Float_t Traj101_CY(Int_t k=0);
    Float_t Traj101_CZ(Int_t k=0);

    Float_t Traj102_CX(Int_t k=0);
    Float_t Traj102_CY(Int_t k=0);
    Float_t Traj102_CZ(Int_t k=0);

    Float_t TrajHTCC_CX(Int_t k=0);
    Float_t TrajHTCC_CY(Int_t k=0);
    Float_t TrajHTCC_CZ(Int_t k=0);

    Float_t TrajFMT1_CX(Int_t k=0);
    Float_t TrajFMT1_CY(Int_t k=0);
    Float_t TrajFMT1_CZ(Int_t k=0);

    Float_t TrajFMT2_CX(Int_t k=0);
    Float_t TrajFMT2_CY(Int_t k=0);
    Float_t TrajFMT2_CZ(Int_t k=0);

    Float_t TrajFMT3_CX(Int_t k=0);
    Float_t TrajFMT3_CY(Int_t k=0);
    Float_t TrajFMT3_CZ(Int_t k=0);

    Float_t TrajFMT4_CX(Int_t k=0);
    Float_t TrajFMT4_CY(Int_t k=0);
    Float_t TrajFMT4_CZ(Int_t k=0);

    Float_t TrajFMT5_CX(Int_t k=0);
    Float_t TrajFMT5_CY(Int_t k=0);
    Float_t TrajFMT5_CZ(Int_t k=0);

    Float_t TrajFMT6_CX(Int_t k=0);
    Float_t TrajFMT6_CY(Int_t k=0);
    Float_t TrajFMT6_CZ(Int_t k=0);


    ///
    Float_t PathTOF(Int_t k);
    Float_t TimeTOF(Int_t k);
    Int_t SectorTOF(Int_t k);
    
    
    //Added in hayk's code
    Double_t Pt2(Int_t, Bool_t = 0);
    Double_t Pl2(Int_t, Bool_t = 0);
    Double_t PlCM(Int_t k, Double_t mp = 0.13957, Bool_t kind = 0);
    Double_t PmaxCM(Bool_t = 0);

    // Kinematic variables
    Double_t Q2(Bool_t kind = 0);
    Double_t W(Bool_t kind = 0);
    Double_t Nu(Bool_t kind = 0);
    Double_t Xb(Bool_t kind = 0);
    Double_t Yb(Bool_t kind = 0);
    Double_t ZhPi(Int_t k, Double_t Mass, Bool_t kind = 0);

    //Added in hayk's code
    Double_t Zh(Int_t k, Double_t mp = 0.13957, Bool_t kind = 0);
    Double_t Xf(Int_t k, Double_t mp = 0.13957, Bool_t kind = 0);
    Double_t Mx2(Int_t k, Float_t mp = 0.13957, Bool_t kind = 0);
    Double_t T(Int_t k, Double_t mp = 0.13957, Bool_t kind = 0);

    // Correction functions
    Double_t TimeCorr4(Double_t mass, Int_t k);
    TVector3 *GetCorrectedVert();

    
    //Added in hayk's code
    //Int_t ElecVertTarg(Bool_t = 0);

    // Particle Identification cuts
    //TString GetCategorization(Int_t k, const char*tt = "",bool mflag=false);
    TString GetCategorization(Int_t k);
    TString GetCategorizationOld(Int_t k);
    TString GetCategorizationMin(Int_t k); 
    TString* GetCategorization();
    int Pid(int k=0,Bool_t kind=0);
    void PrintCategorization();
    void PrintCategorization(TString* partIds);

    // Fiducial Cut
    Double_t FidTheta(Int_t kind = 0);
    Double_t FidThetaMin();
    Double_t FidThetaMinPiPlus(Int_t k);
    Double_t FidFunc(Int_t side, Int_t param);
    Double_t FidFuncPiPlus(Int_t side, Int_t param, Int_t k);
    Double_t FidPhi(Int_t k, Bool_t kind = 0);
    Double_t FidPhiMin();
    Double_t FidPhiMax();
    Double_t FidPhiMinPiPlus(Int_t k);
    Double_t FidPhiMaxPiPlus(Int_t k);
    Bool_t FidCheckCut();
    Bool_t FidCheckCutPiPlus(Int_t k);
    Int_t FidSector(Int_t k, Bool_t kind = 0);

    // Target methods.
    Int_t ElecVertTarg();
    Int_t ElecVertTarg(Bool_t kind);
    Bool_t PionVertTarg(Int_t k);

    // Other methods.
    TVector3 *XYZToUVW(TVector3 *xyz);       
    bool SampFracCheck(const char* tt = "");                         // inline
    // ClassDef(TIdentificatorCLAS12,1)


    int FillResponseMaps();
    int FillMap(hipo::bank *bank, std::map <int,std::vector<int>> &mp, TString piname = "pindex");
    int FillMapRev(hipo::bank *bank, std::map <int,std::vector<int>> &mp, TString piname = "pindex");
    int ClearMaps();
    int PrintMaps();
    int PrintMap(std::map <int,std::vector<int>> &mp);
    long getMapAddr(){return (long)&calorimeterMap;}

private:
#include "node_declaration.h" // got from hipo-io --code file.hipo all and format_nodes.py

    const Double_t kEbeam;    // The energy of incoming electron beam
    const Double_t kMpi;      // The mass of the pion
    const Double_t kMprt;     // The mass of the proton
    const Double_t kMntr;     // The mass of the neutron
    const Double_t kGOOD;     // The key for the exceptions (should be improved to avoid it at all !!!)
    const Bool_t kMCFlag;     //Marco contalbrigo flag, traj dc detId.
    Int_t kIndLundFirst; // index of the first particle comming out in the Lund simulation (the electron)
    hipo::reader *fReader;
    Int_t Nfiles;
    Int_t kCurrentFileIndex;
    std::vector<TString> flist;
    hipo::event *fEvent;
    hipo::dictionary *fFactory;

    /*
    TClasTool *fCT;           // Pointer to the main ClasTool object
    TEVNTClass *fEVNT;        // Pointer to the EVNT object
    TGSIMClass *fGSIM;        // Pointer to the GSIM object
    TCCPBClass *fCCPB;        // Pointer to the CCPB object
    TECPBClass *fECPB;        // Pointer to the ECPB object
    TSCPBClass *fSCPB;        // Pointer to the SCPB object
    TDCPBClass *fDCPB;        // Pointer to the DCPB object
    */
    TString* fPartIds;        // Array with the categories of the particles belonging to an event.
    std::map<TString,int> detectorType;
    std::map<TString,int> layerType;
    std::map<int,int> DCSuperLayer;
    std::map<TString,int> trajDetId;

    int InitBanks();
    int FillBanks();
    int InitDetectorMap();
    int InitLayerMap();
    int InitDCSuperLayerMap();
    int InitTrajDetId();

    std::map <int,std::vector<int>> calorimeterMap;
    std::map <int,std::vector<int>> cherenkovMap;
    std::map <int,std::vector<int>> scintillatorMap;
    std::map <int,std::vector<int>> trackMap;
    std::map <int,std::vector<int>> trajMap;
    std::map <int,std::vector<int>> covMatMap;
    std::map <int,std::vector<int>> richHadPartMap;
    std::map <int,std::vector<int>> richHadClusterMap;
    std::map <int,std::vector<int>> richRRPartMap;
    std::map <int,std::vector<int>> richRRClusterMap;
    std::map <int,std::vector<int>> richHitClusterMap;

    
};

#include "TIdentificatorCLAS12_inlines.h"

#endif

