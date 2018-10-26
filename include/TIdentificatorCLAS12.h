#ifndef TIDENTIFICATOR_H
#define TIDENTIFICATOR_H

#include "TObject.h"
#include "reader.h"
#include "node.h"
#include "TString.h"
#include "TVector3.h"

class TIdentificatorCLAS12 {
public:
  TIdentificatorCLAS12();
  explicit TIdentificatorCLAS12(hipo::reader *reader = 0,Double_t beamE=10.6);
  explicit TIdentificatorCLAS12(TString fname, Double_t beamE=10.6);
    ~TIdentificatorCLAS12();
    
    Bool_t Next(); // explore the hipo file opened

    // HEADER bank
    Float_t NEvent();                             // inline

    // REC::Particle
    Double_t Beta(Int_t k, Bool_t kind = 0);         // inline
    Double_t Id(Int_t k, Bool_t kind = 0);        // inline
    Double_t Charge(Int_t k, Bool_t kind = 0);    // inline
    Double_t Px(Int_t k, Bool_t kind = 0);        // inline
    Double_t Py(Int_t k, Bool_t kind = 0);        // inline
    Double_t Pz(Int_t k, Bool_t kind = 0);        // inline
    Double_t X(Int_t k, Bool_t kind = 0);         // inline
    Double_t Y(Int_t k, Bool_t kind = 0);         // inline
    Double_t Z(Int_t k, Bool_t kind = 0);         // inline
    Int_t Helic(Int_t k);                         //inline
    Int_t GetNRows();                             //inline
     
    Int_t StatCC(Int_t k);                        // inline
    Int_t StatSC(Int_t k);                        // inline
    Int_t StatDC(Int_t k);                        // inline
    Int_t StatEC(Int_t k);                        // inline
    Int_t Status(Int_t k);                     // inline


    
    // CCPB
    Double_t Nphe(Int_t k);                       // inline
    Double_t NpheLTCC(Int_t k);                       // inline
    Double_t NpheHTCC(Int_t k);                       // inline
    Double_t Chi2pid(Int_t k);                     // inline
    Double_t Chi2CC(Int_t k);                     // inline

    Double_t CCStatus(Int_t k);                   // inline

    // DCPB
    Double_t DCStatus(Int_t k);                   // inline

    // ECPB
    Double_t Etot(Int_t k);                       // inline
    Double_t Ein(Int_t k);                        // inline
    Double_t Eout(Int_t k);                       // inline
    Double_t Epcal(Int_t k);                       // inline
    Double_t ECStatus(Int_t k);                   // inline
    Float_t XEC(Int_t k);                         // inline // osoto add
    Float_t YEC(Int_t k);                         // inline // osoto add
    Float_t ZEC(Int_t k);                         // inline // osoto add
    Float_t TimeEC(Int_t k);                      // inline // osoto add
    Float_t PathEC(Int_t k);                      // inline // osoto add

    // SCPB
    Double_t PathSC(Int_t k);                     // inline
    Double_t TimeSC(Int_t k);                     // inline
    Double_t EdepSC(Int_t k);                     // inline
    Double_t SCStatus(Int_t k);                   // inline

    // Derived observables
    Double_t Momentum(Int_t k, Bool_t kind = 0);
    Double_t Moment(Int_t k, Bool_t kind = 0);    // Deprecated
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
	 //Added in hayk's code
    Double_t Pt2(Int_t, Bool_t = 0);
    Double_t Pl2(Int_t, Bool_t = 0);
    Double_t PlCM(Int_t, Bool_t = 0);
    Double_t PmaxCM(Bool_t = 0);

    // Kinematic variables
    Double_t Q2(Bool_t kind = 0);
    Double_t W(Bool_t kind = 0);
    Double_t Nu(Bool_t kind = 0);
    Double_t Xb(Bool_t kind = 0);
    Double_t Yb(Bool_t kind = 0);
    Double_t ZhPi(Int_t k, Double_t Mass, Bool_t kind = 0);

	 //Added in hayk's code
    Double_t Zh(Int_t, Bool_t = 0, Double_t = 0.139);
    Double_t Xf(Int_t, Bool_t = 0);
    Double_t Mx2(Int_t, Bool_t = 0);
    Double_t T(Int_t, Bool_t = 0);


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
    int Pid(int k=0);
    void PrintCategorization();
    void PrintCategorization(TString* partIds);

    // Fiducial Cut
    Double_t FidTheta(Int_t k, Bool_t kind = 0);
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

    //Target methods.
    Int_t ElecVertTarg();
    Int_t ElecVertTarg(Bool_t kind);
    Bool_t PionVertTarg(Int_t k);

    //Other methods.
    TVector3 *XYZToUVW(TVector3 *xyz);            //osoto add
    bool SampFracCheck(const char* tt = "");                         // inline //osoto add 
    //    ClassDef(TIdentificatorCLAS12,1)


    int FillResponseMaps();
    int FillMap(hipo::node<int16_t> *pi,std::map <int,std::vector<int>> &mp);
    int ClearMaps();
    int PrintMaps();
    int PrintMap(std::map <int,std::vector<int>> &mp);


private:
#include "node_declaration.h" // got from hipo-io --code file.hipo all and format_nodes.py

    const Double_t kEbeam;    // The energy of incoming electron beam
    const Double_t kMpi;      // The mass of the pion
    const Double_t kMprt;     // The mass of the proton
    const Double_t kMntr;     // The mass of the neutron
    const Double_t kGOOD;     // The key for the exceptions (should be improved to avoid it at all !!!)
    hipo::reader *fReader;
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
    int InitNodes();
    int InitDetectorMap();
    int InitLayerMap();

    
    std::map <int,std::vector<int>> calorimeterMap;
    std::map <int,std::vector<int>> cherenkovMap;
    std::map <int,std::vector<int>> scintillatorMap;
    std::map <int,std::vector<int>> trackMap;
    std::map <int,std::vector<int>> trajMap;
    std::map <int,std::vector<int>> covMatMap;

    
};

#include "TIdentificatorCLAS12_inlines.h"

#endif

