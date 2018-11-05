#ifndef TIDENTIFICATOR_INLINES_H
#define TIDENTIFICATOR_INLINES_H
extern const Double_t kMeanFePar[6][3];
extern const Double_t kSigmaFePar[6][2];

extern const Double_t kMeanCPar[6][3];
extern const Double_t kSigmaCPar[6][2];

extern const Double_t kMeanPbPar[6][3];
extern const Double_t kSigmaPbPar[6][2];
/*
inline Float_t TIdentificatorCLAS12::NEvent()
{
    // Return the event number from HEAD bank.

    return ((THEADERClass *) fCT->GetHEADER())->GetNEvent();
}
*/



inline Double_t TIdentificatorCLAS12::Beta(Int_t k, Bool_t kind)
{
  if (kind == 0){
    return REC__Particle_beta->getValue(k);
  }  else{
    return -111; // to be implemented
  }
}
/*
inline Double_t TIdentificatorCLAS12::Id(Int_t k, Bool_t kind)
{
    // Return the ID of the particle in the row k, from Particle Data Group
    // (PDG), estimated from some preliminary analysis during data
    // calibration.
    //
    // If kind is zero, the EVNT bank is used. If not, the GSIM bank is used
    // instead.

    if (kind == 0) {
        fEVNT = (TEVNTClass *) fCT->GetBankRow("EVNT", k);
        return fEVNT->Id;
    } else {                            // Fix this in case kind != 1
        fGSIM = (TGSIMClass *) fCT->GetBankRow("GSIM", k);
        return fGSIM->Id;
    }
}

*/

inline Double_t TIdentificatorCLAS12::Charge(Int_t k, Bool_t kind)
{
    // Return the electrical charge for the particle in the row k of the EVNT
    // bank.

    if (kind == 0) {
      return REC__Particle_charge->getValue(k);
    } else {                           
  
      return -111; // to be implemented
    }
}



inline Double_t TIdentificatorCLAS12::Px(Int_t k, Bool_t kind)
{

    if (kind == 0) {
      return REC__Particle_px->getValue(k);
    } else {                            // Fix this in case kind != 1
      return -111;
    }
}



inline Double_t TIdentificatorCLAS12::Py(Int_t k, Bool_t kind)
{
    if (kind == 0) {
	return REC__Particle_py->getValue(k);
    } else {                            // Fix this in case kind != 1
        return -111;
    }
}



inline Double_t TIdentificatorCLAS12::Pz(Int_t k, Bool_t kind)
{
    if (kind == 0) {
	return REC__Particle_pz->getValue(k);
    } else {                            // Fix this in case kind != 1
        return -111;
    }
}



inline Double_t TIdentificatorCLAS12::X(Int_t k, Bool_t kind)
{
    if (kind == 0) {
      return REC__Particle_vx->getValue(k);
    } else { 
      return -111;
    }
}


inline Double_t TIdentificatorCLAS12::Y(Int_t k, Bool_t kind)
{
    if (kind == 0) {
      return REC__Particle_vy->getValue(k);
    } else { 
      return -111;
    }
}


inline Double_t TIdentificatorCLAS12::Z(Int_t k, Bool_t kind)
{
    if (kind == 0) {
      return REC__Particle_vz->getValue(k);
    } else {                            // Fix this
        return -111;
    }
}

inline Int_t TIdentificatorCLAS12::Helic(Int_t k) /// from event
{
  return (Int_t)REC__Event_Helic->getValue(k);
}

inline Int_t TIdentificatorCLAS12::Status(Int_t k)
{
  return (Int_t)REC__Particle_status->getValue(k);
}

/*

inline Int_t TIdentificatorCLAS12::StatSC(Int_t k)
{
    // Return the SCPB key for the particle in the row k of the EVNT bank.

    fEVNT = (TEVNTClass *) fCT->GetBankRow("EVNT", k);
    return fEVNT->Scstat;
}



inline Int_t TIdentificatorCLAS12::StatDC(Int_t k)
{
    // Return the DCPB key for the particle in the row k of the EVNT bank.

    fEVNT = (TEVNTClass *) fCT->GetBankRow("EVNT", k);
    return fEVNT->Dcstat;
}



inline Int_t TIdentificatorCLAS12::StatEC(Int_t k)
{
    // Return the ECPB key for the particle in the row k of the EVNT bank.

    fEVNT = (TEVNTClass *) fCT->GetBankRow("EVNT", k);
    return fEVNT->Ecstat;
}



*/

inline Double_t TIdentificatorCLAS12::NpheLTCC(Int_t k)
{
  int N = cherenkovMap[k].size();
  int index=-1;
  for (int i=0;i<N;i++) // get index of the correct detector
  {
    if (REC__Cherenkov_detector->getValue(cherenkovMap[k][i]) == detectorType["LTCC"])
    {
      index=i;
      break;
    }
    
  }
  if (index>=0)
    return REC__Cherenkov_nphe->getValue(index);
  else
    return -111;
}

inline Double_t TIdentificatorCLAS12::NpheHTCC(Int_t k)
{
  int N = cherenkovMap[k].size();
  int index=-1;
  for (int i=0;i<N;i++) // get index of the correct detector
  {
    if (REC__Cherenkov_detector->getValue(cherenkovMap[k][i]) == detectorType["HTCC"])
    {
      index=i;
      break;
    }
    
  }
  if (index>=0)
    return REC__Cherenkov_nphe->getValue(index);
  else
    return -111;
}



inline Double_t TIdentificatorCLAS12::Chi2pid(Int_t k)
{
  return REC__Particle_chi2pid->getValue(k);
}

inline Int_t TIdentificatorCLAS12::GetNRows()
{
  return  REC__Particle_pid->getLength();
}
/*

inline Double_t TIdentificatorCLAS12::CCStatus(Int_t k)
{
    // Return the signal goodness in the CCPB for the particle in the row k of
    // the EVNT bank.

    fEVNT = (TEVNTClass *) fCT->GetBankRow("EVNT", k);

    if (fEVNT->Ccstat >= 1) {
        Int_t cc_stat = fEVNT->Ccstat - 1;
        fCCPB = (TCCPBClass *) fCT->GetBankRow("CCPB", cc_stat);
        return fCCPB->Status;
    } else {
        return kGOOD;
    }
}



inline Double_t TIdentificatorCLAS12::DCStatus(Int_t k)
{
    // Return the signal goodness in the DCPB for the particle in the row k of
    // the EVNT bank.

    fEVNT = (TEVNTClass *) fCT->GetBankRow("EVNT", k);
    if (fEVNT->Dcstat >= 1) {
        Int_t dc_stat = fEVNT->Dcstat - 1;
        fDCPB = (TDCPBClass *) fCT->GetBankRow("DCPB", dc_stat);
        return fDCPB->Status;
    } else {
        return kGOOD;
    }
}

*/

inline Double_t TIdentificatorCLAS12::Etot(Int_t k)
{
    // Return total energy deposited in the calorimeter for the particle in
    // the row k of the EVNT bank.

  return Ein(k) + Eout(k) + Epcal(k);
}



inline Double_t TIdentificatorCLAS12::Ein(Int_t k)
{
  
  Int_t N = calorimeterMap[k].size();
  Int_t index=-1;
  for (int i=0;i<N;i++)
  {
    if (REC__Calorimeter_detector->getValue(calorimeterMap[k][i])==detectorType["ECAL"]
	&&REC__Calorimeter_layer->getValue(calorimeterMap[k][i])==layerType["EC_Inner"] )
    {
      index = calorimeterMap[k][i];
      break;
    }
  }
  if (index>=0)
    return REC__Calorimeter_energy->getValue(index);
  else 
    return 0;
}



inline Double_t TIdentificatorCLAS12::Eout(Int_t k)
{
  
  Int_t N = calorimeterMap[k].size();
  Int_t index=-1;
  for (int i=0;i<N;i++)
  {
    if (REC__Calorimeter_detector->getValue(calorimeterMap[k][i])==detectorType["ECAL"]
	&&REC__Calorimeter_layer->getValue(calorimeterMap[k][i])==layerType["EC_Outer"] )
    {
      index = calorimeterMap[k][i];
      break;
    }
  }
  if (index>=0)
    return REC__Calorimeter_energy->getValue(index);
  else 
    return 0;
}


inline Double_t TIdentificatorCLAS12::Epcal(Int_t k)
{
  
  Int_t N = calorimeterMap[k].size();
  Int_t index=-1;
  for (int i=0;i<N;i++)
  {
    if (REC__Calorimeter_detector->getValue(calorimeterMap[k][i])==detectorType["ECAL"]
	&&REC__Calorimeter_layer->getValue(calorimeterMap[k][i])==layerType["PCAL"] )
    {
      index = calorimeterMap[k][i];
      break;
    }
  }
  if (index>=0)
    return REC__Calorimeter_energy->getValue(index);
  else 
    return 0;
}


inline Int_t TIdentificatorCLAS12::SectorLTCC(Int_t k,Bool_t kind)
{
  int N = cherenkovMap[k].size();
  int index=-1;
  for (int i=0;i<N;i++) // get index of the correct detector
  {
    if (REC__Cherenkov_detector->getValue(cherenkovMap[k][i]) == detectorType["LTCC"])
    {
      index=i;
      break;
    }
    
  }
  if (index>=0)
    return REC__Cherenkov_sector->getValue(index);
  else
    return -111;
}


inline Int_t TIdentificatorCLAS12::SectorHTCC(Int_t k,Bool_t kind)
{
  int N = cherenkovMap[k].size();
  int index=-1;
  for (int i=0;i<N;i++) // get index of the correct detector
  {
    if (REC__Cherenkov_detector->getValue(cherenkovMap[k][i]) == detectorType["HTCC"])
    {
      index=i;
      break;
    }
    
  }
  if (index>=0)
    return REC__Cherenkov_sector->getValue(index);
  else
    return -111;
}

inline Int_t TIdentificatorCLAS12::SectorECAL(Int_t k,Bool_t kind)
{
  
  Int_t N = calorimeterMap[k].size();
  Int_t index=-1;
  for (int i=0;i<N;i++)
  {
    if (REC__Calorimeter_detector->getValue(calorimeterMap[k][i])==detectorType["ECAL"])
    {
      index = calorimeterMap[k][i];
      break;
    }
  }
  if (index>=0)
    return REC__Calorimeter_sector->getValue(index);
  else 
    return 0;
}


   
inline Float_t TIdentificatorCLAS12::LU_PCAL(Int_t k)
{
  
  Int_t N = calorimeterMap[k].size();
  Int_t index=-1;
  for (int i=0;i<N;i++)
  {
    if (REC__Calorimeter_detector->getValue(calorimeterMap[k][i])==detectorType["ECAL"]
	&&REC__Calorimeter_layer->getValue(calorimeterMap[k][i])==layerType["PCAL"] )
    {
      index = calorimeterMap[k][i];
      break;
    }
  }
  if (index>=0)
    return REC__Calorimeter_lu->getValue(index);
  else 
    return -111;
}

inline Float_t TIdentificatorCLAS12::LV_PCAL(Int_t k)
{
  
  Int_t N = calorimeterMap[k].size();
  Int_t index=-1;
  for (int i=0;i<N;i++)
  {
    if (REC__Calorimeter_detector->getValue(calorimeterMap[k][i])==detectorType["ECAL"]
	&&REC__Calorimeter_layer->getValue(calorimeterMap[k][i])==layerType["PCAL"] )
    {
      index = calorimeterMap[k][i];
      break;
    }
  }
  if (index>=0)
    return REC__Calorimeter_lv->getValue(index);
  else 
    return -111;
}

inline Float_t TIdentificatorCLAS12::LW_PCAL(Int_t k)
{
  
  Int_t N = calorimeterMap[k].size();
  Int_t index=-1;
  for (int i=0;i<N;i++)
  {
    if (REC__Calorimeter_detector->getValue(calorimeterMap[k][i])==detectorType["ECAL"]
	&&REC__Calorimeter_layer->getValue(calorimeterMap[k][i])==layerType["PCAL"] )
    {
      index = calorimeterMap[k][i];
      break;
    }
  }
  if (index>=0)
    return REC__Calorimeter_lw->getValue(index);
  else 
    return -111;
}


   
inline Float_t TIdentificatorCLAS12::LU_ECIN(Int_t k)
{
  
  Int_t N = calorimeterMap[k].size();
  Int_t index=-1;
  for (int i=0;i<N;i++)
  {
    if (REC__Calorimeter_detector->getValue(calorimeterMap[k][i])==detectorType["ECAL"]
	&&REC__Calorimeter_layer->getValue(calorimeterMap[k][i])==layerType["EC_Inner"] )
    {
      index = calorimeterMap[k][i];
      break;
    }
  }
  if (index>=0)
    return REC__Calorimeter_lu->getValue(index);
  else 
    return -111;
}

inline Float_t TIdentificatorCLAS12::LV_ECIN(Int_t k)
{
  
  Int_t N = calorimeterMap[k].size();
  Int_t index=-1;
  for (int i=0;i<N;i++)
  {
    if (REC__Calorimeter_detector->getValue(calorimeterMap[k][i])==detectorType["ECAL"]
	&&REC__Calorimeter_layer->getValue(calorimeterMap[k][i])==layerType["EC_Inner"] )
    {
      index = calorimeterMap[k][i];
      break;
    }
  }
  if (index>=0)
    return REC__Calorimeter_lv->getValue(index);
  else 
    return -111;
}

inline Float_t TIdentificatorCLAS12::LW_ECIN(Int_t k)
{
  
  Int_t N = calorimeterMap[k].size();
  Int_t index=-1;
  for (int i=0;i<N;i++)
  {
    if (REC__Calorimeter_detector->getValue(calorimeterMap[k][i])==detectorType["ECAL"]
	&&REC__Calorimeter_layer->getValue(calorimeterMap[k][i])==layerType["EC_Inner"] )
    {
      index = calorimeterMap[k][i];
      break;
    }
  }
  if (index>=0)
    return REC__Calorimeter_lw->getValue(index);
  else 
    return -111;
}

   
inline Float_t TIdentificatorCLAS12::LU_ECOUT(Int_t k)
{
  
  Int_t N = calorimeterMap[k].size();
  Int_t index=-1;
  for (int i=0;i<N;i++)
  {
    if (REC__Calorimeter_detector->getValue(calorimeterMap[k][i])==detectorType["ECAL"]
	&&REC__Calorimeter_layer->getValue(calorimeterMap[k][i])==layerType["EC_Outer"] )
    {
      index = calorimeterMap[k][i];
      break;
    }
  }
  if (index>=0)
    return REC__Calorimeter_lu->getValue(index);
  else 
    return -111;
}

inline Float_t TIdentificatorCLAS12::LV_ECOUT(Int_t k)
{
  
  Int_t N = calorimeterMap[k].size();
  Int_t index=-1;
  for (int i=0;i<N;i++)
  {
    if (REC__Calorimeter_detector->getValue(calorimeterMap[k][i])==detectorType["ECAL"]
	&&REC__Calorimeter_layer->getValue(calorimeterMap[k][i])==layerType["EC_Outer"] )
    {
      index = calorimeterMap[k][i];
      break;
    }
  }
  if (index>=0)
    return REC__Calorimeter_lv->getValue(index);
  else 
    return -111;
}

inline Float_t TIdentificatorCLAS12::LW_ECOUT(Int_t k)
{
  
  Int_t N = calorimeterMap[k].size();
  Int_t index=-1;
  for (int i=0;i<N;i++)
  {
    if (REC__Calorimeter_detector->getValue(calorimeterMap[k][i])==detectorType["ECAL"]
	&&REC__Calorimeter_layer->getValue(calorimeterMap[k][i])==layerType["EC_Outer"] )
    {
      index = calorimeterMap[k][i];
      break;
    }
  }
  if (index>=0)
    return REC__Calorimeter_lw->getValue(index);
  else 
    return -111;
}



inline Float_t TIdentificatorCLAS12::HX_PCAL(Int_t k)
{
  
  Int_t N = calorimeterMap[k].size();
  Int_t index=-1;
  for (int i=0;i<N;i++)
  {
    if (REC__Calorimeter_detector->getValue(calorimeterMap[k][i])==detectorType["ECAL"]
	&&REC__Calorimeter_layer->getValue(calorimeterMap[k][i])==layerType["PCAL"] )
    {
      index = calorimeterMap[k][i];
      break;
    }
  }
  if (index>=0)
    return REC__Calorimeter_hx->getValue(index);
  else 
    return -111;
}

inline Float_t TIdentificatorCLAS12::HY_PCAL(Int_t k)
{
  
  Int_t N = calorimeterMap[k].size();
  Int_t index=-1;
  for (int i=0;i<N;i++)
  {
    if (REC__Calorimeter_detector->getValue(calorimeterMap[k][i])==detectorType["ECAL"]
	&&REC__Calorimeter_layer->getValue(calorimeterMap[k][i])==layerType["PCAL"] )
    {
      index = calorimeterMap[k][i];
      break;
    }
  }
  if (index>=0)
    return REC__Calorimeter_hy->getValue(index);
  else 
    return -111;
}

inline Float_t TIdentificatorCLAS12::HZ_PCAL(Int_t k)
{
  
  Int_t N = calorimeterMap[k].size();
  Int_t index=-1;
  for (int i=0;i<N;i++)
  {
    if (REC__Calorimeter_detector->getValue(calorimeterMap[k][i])==detectorType["ECAL"]
	&&REC__Calorimeter_layer->getValue(calorimeterMap[k][i])==layerType["PCAL"] )
    {
      index = calorimeterMap[k][i];
      break;
    }
  }
  if (index>=0)
    return REC__Calorimeter_hz->getValue(index);
  else 
    return -111;
}

inline Float_t TIdentificatorCLAS12::HX_ECIN(Int_t k)
{
  
  Int_t N = calorimeterMap[k].size();
  Int_t index=-1;
  for (int i=0;i<N;i++)
  {
    if (REC__Calorimeter_detector->getValue(calorimeterMap[k][i])==detectorType["ECAL"]
	&&REC__Calorimeter_layer->getValue(calorimeterMap[k][i])==layerType["EC_Inner"] )
    {
      index = calorimeterMap[k][i];
      break;
    }
  }
  if (index>=0)
    return REC__Calorimeter_hx->getValue(index);
  else 
    return -111;
}

inline Float_t TIdentificatorCLAS12::HY_ECIN(Int_t k)
{
  
  Int_t N = calorimeterMap[k].size();
  Int_t index=-1;
  for (int i=0;i<N;i++)
  {
    if (REC__Calorimeter_detector->getValue(calorimeterMap[k][i])==detectorType["ECAL"]
	&&REC__Calorimeter_layer->getValue(calorimeterMap[k][i])==layerType["EC_Inner"] )
    {
      index = calorimeterMap[k][i];
      break;
    }
  }
  if (index>=0)
    return REC__Calorimeter_hy->getValue(index);
  else 
    return -111;
}

inline Float_t TIdentificatorCLAS12::HZ_ECIN(Int_t k)
{
  
  Int_t N = calorimeterMap[k].size();
  Int_t index=-1;
  for (int i=0;i<N;i++)
  {
    if (REC__Calorimeter_detector->getValue(calorimeterMap[k][i])==detectorType["ECAL"]
	&&REC__Calorimeter_layer->getValue(calorimeterMap[k][i])==layerType["EC_Inner"] )
    {
      index = calorimeterMap[k][i];
      break;
    }
  }
  if (index>=0)
    return REC__Calorimeter_hz->getValue(index);
  else 
    return -111;
}

inline Float_t TIdentificatorCLAS12::HX_ECOUT(Int_t k)
{
  
  Int_t N = calorimeterMap[k].size();
  Int_t index=-1;
  for (int i=0;i<N;i++)
  {
    if (REC__Calorimeter_detector->getValue(calorimeterMap[k][i])==detectorType["ECAL"]
	&&REC__Calorimeter_layer->getValue(calorimeterMap[k][i])==layerType["EC_Outer"] )
    {
      index = calorimeterMap[k][i];
      break;
    }
  }
  if (index>=0)
    return REC__Calorimeter_hx->getValue(index);
  else 
    return -111;
}

inline Float_t TIdentificatorCLAS12::HY_ECOUT(Int_t k)
{
  
  Int_t N = calorimeterMap[k].size();
  Int_t index=-1;
  for (int i=0;i<N;i++)
  {
    if (REC__Calorimeter_detector->getValue(calorimeterMap[k][i])==detectorType["ECAL"]
	&&REC__Calorimeter_layer->getValue(calorimeterMap[k][i])==layerType["EC_Outer"] )
    {
      index = calorimeterMap[k][i];
      break;
    }
  }
  if (index>=0)
    return REC__Calorimeter_hy->getValue(index);
  else 
    return -111;
}

inline Float_t TIdentificatorCLAS12::HZ_ECOUT(Int_t k)
{
  
  Int_t N = calorimeterMap[k].size();
  Int_t index=-1;
  for (int i=0;i<N;i++)
  {
    if (REC__Calorimeter_detector->getValue(calorimeterMap[k][i])==detectorType["ECAL"]
	&&REC__Calorimeter_layer->getValue(calorimeterMap[k][i])==layerType["EC_Outer"] )
    {
      index = calorimeterMap[k][i];
      break;
    }
  }
  if (index>=0)
    return REC__Calorimeter_hz->getValue(index);
  else 
    return -111;
}




/*

inline Double_t TIdentificatorCLAS12::ECStatus(Int_t k)
{
    // Return the signal goodness in the ECPB for the particle in the row k of
    // the EVNT bank.

    fEVNT = (TEVNTClass *) fCT->GetBankRow("EVNT", k);
    if (fEVNT->Ecstat >= 1) {
        Int_t ec_stat = fEVNT->Ecstat - 1;
        fECPB = (TECPBClass *) fCT->GetBankRow("ECPB", ec_stat);
        return fECPB->Status;
    } else {
        return kGOOD;
    }
}



inline Double_t TIdentificatorCLAS12::PathSC(Int_t k)
{
    // Return the path length from target, in the SCPB, for the particle in
    // the row k of the EVNT bank.

    fEVNT = (TEVNTClass *) fCT->GetBankRow("EVNT", k);

    if (fEVNT->Scstat >= 1) {
        Int_t sc_stat = fEVNT->Scstat - 1;
        fSCPB = (TSCPBClass *) fCT->GetBankRow("SCPB", sc_stat);
        return fSCPB->Path;
    } else {
        return kGOOD;
    }
}



inline Double_t TIdentificatorCLAS12::TimeSC(Int_t k)
{
    // Return the flight time relative to the event start time, in the SCPB,
    // for the particle in the row k of the EVNT bank.

    fEVNT = (TEVNTClass *) fCT->GetBankRow("EVNT", k);

    if (fEVNT->Scstat >= 1) {
        Int_t sc_stat = fEVNT->Scstat - 1;
        fSCPB = (TSCPBClass *) fCT->GetBankRow("SCPB", sc_stat);
        return fSCPB->Time;
    } else {
        return kGOOD;
    }
}



inline Double_t TIdentificatorCLAS12::EdepSC(Int_t k)
{
    // Return the deposited energy (dE/dX) in the SCPB, for the particle in
    // the row k of the EVNT bank.

    fEVNT = (TEVNTClass *) fCT->GetBankRow("EVNT", k);

    if (fEVNT->Scstat >= 1) {
        Int_t sc_stat = fEVNT->Scstat - 1;
        fSCPB = (TSCPBClass *) fCT->GetBankRow("SCPB", sc_stat);
        return fSCPB->Edep;
    } else {
        return kGOOD;
    }
}



inline Double_t TIdentificatorCLAS12::SCStatus(Int_t k)
{
    // Return the signal goodness in the SCPB for the particle in the row k of
    // the EVNT bank.

    fEVNT = (TEVNTClass *) fCT->GetBankRow("EVNT", k);

    if (fEVNT->Scstat >= 1) {
        Int_t sc_stat = fEVNT->Scstat -1;
        fSCPB = (TSCPBClass *) fCT->GetBankRow("SCPB", sc_stat);
        return fSCPB->Status;
    } else {
        return kGOOD;
    }
}

inline Float_t TIdentificatorCLAS12::XEC(Int_t k)
{
    fEVNT = (TEVNTClass *) fCT->GetBankRow("EVNT", k);

    if (fEVNT->Ecstat >= 1) {
        Int_t ec_stat= fEVNT->Ecstat - 1;
        fECPB = (TECPBClass *) fCT->GetBankRow("ECPB", ec_stat);
        return fECPB->X;
    } else {
        return kGOOD;
    }
}

inline Float_t TIdentificatorCLAS12::YEC(Int_t k)
{
    fEVNT = (TEVNTClass *) fCT->GetBankRow("EVNT", k);

    if (fEVNT->Ecstat >= 1) {
        Int_t ec_stat= fEVNT->Ecstat - 1;
        fECPB = (TECPBClass *) fCT->GetBankRow("ECPB", ec_stat);
        return fECPB->Y;
    } else {
        return kGOOD;
    }
}

inline Float_t TIdentificatorCLAS12::ZEC(Int_t k)
{
    fEVNT = (TEVNTClass *) fCT->GetBankRow("EVNT", k);

    if (fEVNT->Ecstat >= 1) {
        Int_t ec_stat= fEVNT->Ecstat - 1;
        fECPB = (TECPBClass *) fCT->GetBankRow("ECPB", ec_stat);
        return fECPB->Z;
    } else {
        return kGOOD;
    }
}

inline Float_t TIdentificatorCLAS12::TimeEC(Int_t k)
{
    fEVNT = (TEVNTClass *) fCT->GetBankRow("EVNT", k);

    if (fEVNT->Ecstat >= 1) {
        Int_t ec_stat= fEVNT->Ecstat - 1;
        fECPB = (TECPBClass *) fCT->GetBankRow("ECPB", ec_stat);
        return fECPB->Time;
    } else {
        return kGOOD;
    }
}

inline Float_t TIdentificatorCLAS12::PathEC(Int_t k)
{
    fEVNT = (TEVNTClass *) fCT->GetBankRow("EVNT", k);

    if (fEVNT->Ecstat >= 1) {
        Int_t ec_stat= fEVNT->Ecstat - 1;
        fECPB = (TECPBClass *) fCT->GetBankRow("ECPB", ec_stat);
        return fECPB->Path;
    } else {
        return kGOOD;
    }
}

inline bool TIdentificatorCLAS12::SampFracCheck(const char *tt)
{
  if ( !strcmp(tt,"") ) return true;
  bool ret=false;
  Double_t E=TMath::Max(Etot(0),Ein(0)+Eout(0));
  Double_t p=Momentum(0);
  Int_t s = Sector(0);
  Double_t MeanE;
  Double_t SigmaE;
  if (!strcmp(tt,"Fe") ) 
  {
     MeanE = kMeanFePar[s][0] + kMeanFePar[s][1]*p + kMeanFePar[s][2]*p*p;
     SigmaE = TMath::Sqrt(kSigmaFePar[s][0]*kSigmaFePar[s][0] + kSigmaFePar[s][1]*kSigmaFePar[s][1]/p);
  }
  else if (!strcmp(tt,"C") ) 
  {
     MeanE = kMeanCPar[s][0] + kMeanCPar[s][1]*p + kMeanCPar[s][2]*p*p;
     SigmaE = TMath::Sqrt(kSigmaCPar[s][0]*kSigmaCPar[s][0] + kSigmaCPar[s][1]*kSigmaCPar[s][1]/p);
  }
  else if (!strcmp(tt,"Pb") ) 
  {
     MeanE = kMeanPbPar[s][0] + kMeanPbPar[s][1]*p + kMeanPbPar[s][2]*p*p;
     SigmaE = TMath::Sqrt(kSigmaPbPar[s][0]*kSigmaPbPar[s][0] + kSigmaPbPar[s][1]*kSigmaPbPar[s][1]/p);
  }
  if (TMath::Abs( E/p - MeanE )< 2.5*SigmaE)
    ret=true;
  return ret;
}
*/
#endif
