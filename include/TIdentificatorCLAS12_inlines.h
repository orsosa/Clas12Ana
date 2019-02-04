
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



inline Double_t TIdentificatorCLAS12::Beta(Int_t k)
{
    return REC__Particle_beta->getValue(k);
}

inline Double_t TIdentificatorCLAS12::Charge(Int_t k, Bool_t kind)
{
    // Return the electrical charge for the particle in the row k of the EVNT
    // bank.

    if (kind == 0) {
      return REC__Particle_charge->getValue(k);
    } else {
      return -111; // Not found in MC__xxx banks
    }
}



inline Double_t TIdentificatorCLAS12::Px(Int_t k, Bool_t kind)
{

    if (kind == 0) {
      return REC__Particle_px->getValue(k);
    } else {                
      return MC__Lund_px->getValue(k);
    }
}



inline Double_t TIdentificatorCLAS12::Py(Int_t k, Bool_t kind)
{
    if (kind == 0) {
      return REC__Particle_py->getValue(k);
    } else {                
      return MC__Lund_py->getValue(k);
    }
}



inline Double_t TIdentificatorCLAS12::Pz(Int_t k, Bool_t kind)
{
    if (kind == 0) {
      return REC__Particle_pz->getValue(k);
    } else {                
      return MC__Lund_pz->getValue(k);
    }
}



inline Double_t TIdentificatorCLAS12::X(Int_t k, Bool_t kind)
{
    if (kind == 0) {
      return REC__Particle_vx->getValue(k);
    } else { 
      return MC__Lund_vx->getValue(k);
    }
}


inline Double_t TIdentificatorCLAS12::Y(Int_t k, Bool_t kind)
{
    if (kind == 0) {
      return REC__Particle_vy->getValue(k);
    } else { 
      return MC__Lund_vy->getValue(k);
    }
}


inline Double_t TIdentificatorCLAS12::Z(Int_t k, Bool_t kind)
{
    if (kind == 0) {
      return REC__Particle_vz->getValue(k);
    } else {                            // Fix this
      return MC__Lund_vz->getValue(k);
    }
}

inline Int_t TIdentificatorCLAS12::Helic(Int_t k) /// from event
{
  return (Int_t)REC__Event_Helic->getValue(k);
}

inline Int_t TIdentificatorCLAS12::Event(Int_t k) /// from event
{
  return (Int_t)RUN__config_event->getValue(k);
}

inline Float_t TIdentificatorCLAS12::STTime(Int_t k) /// from event
{
  return (Int_t)REC__Event_STTime->getValue(k);
}

inline Float_t TIdentificatorCLAS12::RFTime(Int_t k) /// from event
{
  return (Int_t)REC__Event_RFTime->getValue(k);
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

inline Int_t TIdentificatorCLAS12::GetMCNRows()
{
  return  MC__Lund_pid->getLength();
}

inline Float_t TIdentificatorCLAS12::MCMass(Int_t k)
{
  return MC__Lund_mass->getValue(k);
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

inline Double_t TIdentificatorCLAS12::Etot(Int_t k,Bool_t kind)
{
    // Return total energy deposited in the calorimeter for the particle in
    // the row k of the EVNT bank.
  if (kind==0)
    return Ein(k) + Eout(k) + Epcal(k);
  else
    return MC__Lund_E->getValue(k);
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


inline Int_t TIdentificatorCLAS12::SectorLTCC(Int_t k)
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


inline Int_t TIdentificatorCLAS12::SectorHTCC(Int_t k)
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

inline Int_t TIdentificatorCLAS12::SectorECAL(Int_t k)
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



inline Float_t TIdentificatorCLAS12::VX_DC(Int_t k)
{
  
  Int_t N = trackMap[k].size();
  Int_t index=-1;
  for (int i=0;i<N;i++)
  {
    if (REC__Track_detector->getValue(trackMap[k][i])==detectorType["DC"])
    {
      index = trackMap[k][i];
      break;
    }
  }
  if (index>=0)
    return REC__Track_vx_nomm->getValue(index);
  else 
    return -1111;
}

inline Float_t TIdentificatorCLAS12::VY_DC(Int_t k)
{
  
  Int_t N = trackMap[k].size();
  Int_t index=-1;
  for (int i=0;i<N;i++)
  {
    if (REC__Track_detector->getValue(trackMap[k][i])==detectorType["DC"])
    {
      index = trackMap[k][i];
      break;
    }
  }
  if (index>=0)
    return REC__Track_vy_nomm->getValue(index);
  else 
    return -1111;
}

inline Float_t TIdentificatorCLAS12::VZ_DC(Int_t k)
{
  
  Int_t N = trackMap[k].size();
  Int_t index=-1;
  for (int i=0;i<N;i++)
  {
    if (REC__Track_detector->getValue(trackMap[k][i])==detectorType["DC"])
    {
      index = trackMap[k][i];
      break;
    }
  }
  if (index>=0)
    return REC__Track_vz_nomm->getValue(index);
  else 
    return -1111;
}

inline Int_t TIdentificatorCLAS12::SectorDC(Int_t k)
{
  Int_t N = trackMap[k].size();
  Int_t index=-1;
  for (int i=0;i<N;i++)
  {  
    if (REC__Track_detector->getValue(trackMap[k][i])==detectorType["DC"])
    {
      index = trackMap[k][i];
      break;
    }
  }

  if (index>=0){
    REC__Track_sector->getValue(index);
    return REC__Track_sector->getValue(index);}
  else 
    return -111;
}

inline Float_t TIdentificatorCLAS12::Px_DC(Int_t k)
{
  
  Int_t N = trackMap[k].size();
  Int_t index=-1;
  for (int i=0;i<N;i++)
  {
    if (REC__Track_detector->getValue(trackMap[k][i])==detectorType["DC"])
    {
      index = trackMap[k][i];
      break;
    }
  }
  if (index>=0)
    return REC__Track_px_nomm->getValue(index);
  else 
    return -111;
}

inline Float_t TIdentificatorCLAS12::Py_DC(Int_t k)
{
  
  Int_t N = trackMap[k].size();
  Int_t index=-1;
  for (int i=0;i<N;i++)
  {
    if (REC__Track_detector->getValue(trackMap[k][i])==detectorType["DC"])
    {
      index = trackMap[k][i];
      break;
    }
  }
  if (index>=0)
    return REC__Track_py_nomm->getValue(index);
  else 
    return -111;
}

inline Float_t TIdentificatorCLAS12::Pz_DC(Int_t k)
{
  
  Int_t N = trackMap[k].size();
  Int_t index=-1;
  for (int i=0;i<N;i++)
  {
    if (REC__Track_detector->getValue(trackMap[k][i])==detectorType["DC"])
    {
      index = trackMap[k][i];
      break;
    }
  }
  if (index>=0)
    return REC__Track_pz_nomm->getValue(index);
  else 
    return -111;
}


inline Float_t TIdentificatorCLAS12::TrajX(Int_t k,Int_t sl)
{
  
  Int_t N = trajMap[k].size();
  Int_t index=-1;
  for (int i=0;i<N;i++)
  {
    if (REC__Traj_detId->getValue(trajMap[k][i])==DCSuperLayer[sl])
    {
      index = trajMap[k][i];
      break;
    }
  }
  if (index>=0)
    return REC__Traj_x->getValue(index);
  else 
    return -1111;
}

inline Float_t TIdentificatorCLAS12::TrajY(Int_t k,Int_t sl)
{
  
  Int_t N = trajMap[k].size();
  Int_t index=-1;
  for (int i=0;i<N;i++)
  {
    if (REC__Traj_detId->getValue(trajMap[k][i])==DCSuperLayer[sl])
    {
      index = trajMap[k][i];
      break;
    }
  }
  if (index>=0)
    return REC__Traj_y->getValue(index);
  else 
    return -1111;
}

inline Float_t TIdentificatorCLAS12::TrajZ(Int_t k,Int_t sl)
{
  
  Int_t N = trajMap[k].size();
  Int_t index=-1;
  for (int i=0;i<N;i++)
  {
    if (REC__Traj_detId->getValue(trajMap[k][i])==DCSuperLayer[sl])
    {
      index = trajMap[k][i];
      break;
    }
  }
  if (index>=0)
    return REC__Traj_z->getValue(index);
  else 
    return -1111;
}

/////// traj dc cx,cy,cz on each super layer

inline Float_t TIdentificatorCLAS12::TrajCX(Int_t k,Int_t sl)
{
  
  Int_t N = trajMap[k].size();
  Int_t index=-1;
  for (int i=0;i<N;i++)
  {
    if (REC__Traj_detId->getValue(trajMap[k][i])==DCSuperLayer[sl])
    {
      index = trajMap[k][i];
      break;
    }
  }
  if (index>=0)
    return REC__Traj_cx->getValue(index);
  else 
    return -1111;
}

inline Float_t TIdentificatorCLAS12::TrajCY(Int_t k,Int_t sl)
{
  
  Int_t N = trajMap[k].size();
  Int_t index=-1;
  for (int i=0;i<N;i++)
  {
    if (REC__Traj_detId->getValue(trajMap[k][i])==DCSuperLayer[sl])
    {
      index = trajMap[k][i];
      break;
    }
  }
  if (index>=0)
    return REC__Traj_cy->getValue(index);
  else 
    return -1111;
}

inline Float_t TIdentificatorCLAS12::TrajCZ(Int_t k,Int_t sl)
{
  
  Int_t N = trajMap[k].size();
  Int_t index=-1;
  for (int i=0;i<N;i++)
  {
    if (REC__Traj_detId->getValue(trajMap[k][i])==DCSuperLayer[sl])
    {
      index = trajMap[k][i];
      break;
    }
  }
  if (index>=0)
    return REC__Traj_cz->getValue(index);
  else 
    return -1111;
}

////
/// Traj detID XYZ ///
inline Float_t TIdentificatorCLAS12::TrajDetIdX(Int_t k,TString dname)
{
  
  Int_t N = trajMap[k].size();
  Int_t index=-1;
  for (int i=0;i<N;i++)
  {
    if (REC__Traj_detId->getValue(trajMap[k][i])==trajDetId[dname])
    {
      index = trajMap[k][i];
      break;
    }
  }
  if (index>=0)
    return REC__Traj_x->getValue(index);
  else 
    return -1111;
}

inline Float_t TIdentificatorCLAS12::TrajDetIdY(Int_t k,TString dname)
{
  
  Int_t N = trajMap[k].size();
  Int_t index=-1;
  for (int i=0;i<N;i++)
  {
    if (REC__Traj_detId->getValue(trajMap[k][i])==trajDetId[dname])
    {
      index = trajMap[k][i];
      break;
    }
  }
  if (index>=0)
    return REC__Traj_y->getValue(index);
  else 
    return -1111;
}

inline Float_t TIdentificatorCLAS12::TrajDetIdZ(Int_t k,TString dname)
{
  
  Int_t N = trajMap[k].size();
  Int_t index=-1;
  for (int i=0;i<N;i++)
  {
    if (REC__Traj_detId->getValue(trajMap[k][i])==trajDetId[dname])
    {
      index = trajMap[k][i];
      break;
    }
  }
  if (index>=0)
    return REC__Traj_z->getValue(index);
  else 
    return -1111;
}

//////////
/// Traj XYZ EC
inline Float_t TIdentificatorCLAS12::TrajECX(Int_t k)
{
  return TrajDetIdX(k,"EC");
}
inline Float_t TIdentificatorCLAS12::TrajECY(Int_t k)
{
  return TrajDetIdY(k,"EC");
}
inline Float_t TIdentificatorCLAS12::TrajECZ(Int_t k)
{
  return TrajDetIdZ(k,"EC");
}
/////////////
////// Traj XYZ PCAL
inline Float_t TIdentificatorCLAS12::TrajPCALX(Int_t k)
{
  return TrajDetIdX(k,"PCAL");
}
inline Float_t TIdentificatorCLAS12::TrajPCALY(Int_t k)
{
  return TrajDetIdY(k,"PCAL");
}
inline Float_t TIdentificatorCLAS12::TrajPCALZ(Int_t k)
{
  return TrajDetIdZ(k,"PCAL");
}
//////////
////// Traj XYZ LTCC
inline Float_t TIdentificatorCLAS12::TrajLTCCX(Int_t k)
{
  return TrajDetIdX(k,"LTCC");
}
inline Float_t TIdentificatorCLAS12::TrajLTCCY(Int_t k)
{
  return TrajDetIdY(k,"LTCC");
}
inline Float_t TIdentificatorCLAS12::TrajLTCCZ(Int_t k)
{
  return TrajDetIdZ(k,"LTCC");
}
////////
////// Traj XYZ FTOF1A
inline Float_t TIdentificatorCLAS12::TrajFTOF1AX(Int_t k)
{
  return TrajDetIdX(k,"FTOT1A");
}
inline Float_t TIdentificatorCLAS12::TrajFTOF1AY(Int_t k)
{
  return TrajDetIdY(k,"FTOT1A");
}
inline Float_t TIdentificatorCLAS12::TrajFTOF1AZ(Int_t k)
{
  return TrajDetIdZ(k,"FTOT1A");
}
////////
////// Traj XYZ FTOF1B
inline Float_t TIdentificatorCLAS12::TrajFTOF1BX(Int_t k)
{
  return TrajDetIdX(k,"FTOF1B");
}
inline Float_t TIdentificatorCLAS12::TrajFTOF1BY(Int_t k)
{
  return TrajDetIdY(k,"FTOF1B");
}
inline Float_t TIdentificatorCLAS12::TrajFTOF1BZ(Int_t k)
{
  return TrajDetIdZ(k,"FTOF1B");
}
////////
////// Traj XYZ FTOF2
inline Float_t TIdentificatorCLAS12::TrajFTOF2X(Int_t k)
{
  return TrajDetIdX(k,"FTOF2");
}
inline Float_t TIdentificatorCLAS12::TrajFTOF2Y(Int_t k)
{
  return TrajDetIdY(k,"FTOF2");
}
inline Float_t TIdentificatorCLAS12::TrajFTOF2Z(Int_t k)
{
  return TrajDetIdZ(k,"FTOF2");
}
////////
////// Traj XYZ 101
inline Float_t TIdentificatorCLAS12::Traj101X(Int_t k)
{
  return TrajDetIdX(k,"101");
}
inline Float_t TIdentificatorCLAS12::Traj101Y(Int_t k)
{
  return TrajDetIdY(k,"101");
}
inline Float_t TIdentificatorCLAS12::Traj101Z(Int_t k)
{
  return TrajDetIdZ(k,"101");
}
////////
////// Traj XYZ 102
inline Float_t TIdentificatorCLAS12::Traj102X(Int_t k)
{
  return TrajDetIdX(k,"102");
}
inline Float_t TIdentificatorCLAS12::Traj102Y(Int_t k)
{
  return TrajDetIdY(k,"102");
}
inline Float_t TIdentificatorCLAS12::Traj102Z(Int_t k)
{
  return TrajDetIdZ(k,"102");
}
///////////
////// Traj XYZ HTCC
inline Float_t TIdentificatorCLAS12::TrajHTCCX(Int_t k)
{
  return TrajDetIdX(k,"HTCC");
}
inline Float_t TIdentificatorCLAS12::TrajHTCCY(Int_t k)
{
  return TrajDetIdY(k,"HTCC");
}
inline Float_t TIdentificatorCLAS12::TrajHTCCZ(Int_t k)
{
  return TrajDetIdZ(k,"HTCC");
}
///////////
////// Traj XYZ FMT1
inline Float_t TIdentificatorCLAS12::TrajFMT1X(Int_t k)
{
  return TrajDetIdX(k,"FMT1");
}
inline Float_t TIdentificatorCLAS12::TrajFMT1Y(Int_t k)
{
  return TrajDetIdY(k,"FMT1");
}
inline Float_t TIdentificatorCLAS12::TrajFMT1Z(Int_t k)
{
  return TrajDetIdZ(k,"FMT1");
}
///////////
////// Traj XYZ FMT2
inline Float_t TIdentificatorCLAS12::TrajFMT2X(Int_t k)
{
  return TrajDetIdX(k,"FMT2");
}
inline Float_t TIdentificatorCLAS12::TrajFMT2Y(Int_t k)
{
  return TrajDetIdY(k,"FMT2");
}
inline Float_t TIdentificatorCLAS12::TrajFMT2Z(Int_t k)
{
  return TrajDetIdZ(k,"FMT2");
}
///////////
////// Traj XYZ FMT3
inline Float_t TIdentificatorCLAS12::TrajFMT3X(Int_t k)
{
  return TrajDetIdX(k,"FMT3");
}
inline Float_t TIdentificatorCLAS12::TrajFMT3Y(Int_t k)
{
  return TrajDetIdY(k,"FMT3");
}
inline Float_t TIdentificatorCLAS12::TrajFMT3Z(Int_t k)
{
  return TrajDetIdZ(k,"FMT3");
}
///////////
////// Traj XYZ FMT4
inline Float_t TIdentificatorCLAS12::TrajFMT4X(Int_t k)
{
  return TrajDetIdX(k,"FMT4");
}
inline Float_t TIdentificatorCLAS12::TrajFMT4Y(Int_t k)
{
  return TrajDetIdY(k,"FMT4");
}
inline Float_t TIdentificatorCLAS12::TrajFMT4Z(Int_t k)
{
  return TrajDetIdZ(k,"FMT4");
}
///////////
////// Traj XYZ FMT5
inline Float_t TIdentificatorCLAS12::TrajFMT5X(Int_t k)
{
  return TrajDetIdX(k,"FMT5");
}
inline Float_t TIdentificatorCLAS12::TrajFMT5Y(Int_t k)
{
  return TrajDetIdY(k,"FMT5");
}
inline Float_t TIdentificatorCLAS12::TrajFMT5Z(Int_t k)
{
  return TrajDetIdZ(k,"FMT5");
}
///////////
////// Traj XYZ FMT6
inline Float_t TIdentificatorCLAS12::TrajFMT6X(Int_t k)
{
  return TrajDetIdX(k,"FMT6");
}
inline Float_t TIdentificatorCLAS12::TrajFMT6Y(Int_t k)
{
  return TrajDetIdY(k,"FMT6");
}
inline Float_t TIdentificatorCLAS12::TrajFMT6Z(Int_t k)
{
  return TrajDetIdZ(k,"FMT6");
}
///////////


inline Float_t TIdentificatorCLAS12::TrajDCX(Int_t k,Int_t reg)
{
  Int_t N = trajMap[k].size();
  Int_t index=-1;
  int sl=2*reg;
  for (int i=0;i<N;i++)
  {
    if (REC__Traj_detId->getValue(trajMap[k][i])==DCSuperLayer[sl] || REC__Traj_detId->getValue(trajMap[k][i])==DCSuperLayer[sl+1] )
    {
      index = trajMap[k][i];
      break;
    }
  }
  if (index>=0)
    return REC__Traj_x->getValue(index);
  else 
    return -1111;
}

inline Float_t TIdentificatorCLAS12::TrajDCY(Int_t k,Int_t reg)
{
  
  Int_t N = trajMap[k].size();
  Int_t index=-1;
  int sl=2*reg;
  for (int i=0;i<N;i++)
  {
    if (REC__Traj_detId->getValue(trajMap[k][i])==DCSuperLayer[sl] || REC__Traj_detId->getValue(trajMap[k][i])==DCSuperLayer[sl+1] )
    {
      index = trajMap[k][i];
      break;
    }
  }
  if (index>=0)
    return REC__Traj_y->getValue(index);
  else 
    return -1111;
}

inline Float_t TIdentificatorCLAS12::TrajDCZ(Int_t k,Int_t reg)
{
  
  Int_t N = trajMap[k].size();
  Int_t index=-1;
  int sl=2*reg;
  for (int i=0;i<N;i++)
  {
    if (REC__Traj_detId->getValue(trajMap[k][i])==DCSuperLayer[sl] || REC__Traj_detId->getValue(trajMap[k][i])==DCSuperLayer[sl+1] )
    {
      index = trajMap[k][i];
      break;
    }
  }
  if (index>=0)
    return REC__Traj_z->getValue(index);
  else 
    return -1111;
}


inline Float_t TIdentificatorCLAS12::PathTOF(Int_t k)
{
  Int_t N = scintillatorMap[k].size();
  Int_t index=-1;
  for (int i=0;i<N;i++)
  {
    if (REC__Scintillator_detector->getValue(scintillatorMap[k][i])==detectorType["FTOF"])
    {
      if (REC__Scintillator_layer->getValue(scintillatorMap[k][i])==layerType["FTOF_1A"] )
      {
	index = scintillatorMap[k][i];
	break;
      }
      else if(REC__Scintillator_layer->getValue(scintillatorMap[k][i])==layerType["FTOF_1B"])
      {
	index = scintillatorMap[k][i];
	break;
      }
      else if(REC__Scintillator_layer->getValue(scintillatorMap[k][i])==layerType["FTOF_2"])
      {
	index = scintillatorMap[k][i];
	break;
      }
    }
    
  }
  if (index>=0)
    return REC__Scintillator_path->getValue(index);
  else 
    return -111;
}

inline Float_t TIdentificatorCLAS12::TimeTOF(Int_t k)
{
  Int_t N = scintillatorMap[k].size();
  Int_t index=-1;
  for (int i=0;i<N;i++)
  {
    if (REC__Scintillator_detector->getValue(scintillatorMap[k][i])==detectorType["FTOF"])
    {
      if (REC__Scintillator_layer->getValue(scintillatorMap[k][i])==layerType["FTOF_1A"] )
      {
	index = scintillatorMap[k][i];
	break;
      }
      else if(REC__Scintillator_layer->getValue(scintillatorMap[k][i])==layerType["FTOF_1B"])
      {
	index = scintillatorMap[k][i];
	break;
      }
      else if(REC__Scintillator_layer->getValue(scintillatorMap[k][i])==layerType["FTOF_2"])
      {
	index = scintillatorMap[k][i];
	break;
      }
    }
    
  }
  if (index>=0)
    return REC__Scintillator_time->getValue(index);
  else 
    return -111;
}

inline Int_t TIdentificatorCLAS12::SectorTOF(Int_t k)
{
  Int_t N = scintillatorMap[k].size();
  Int_t index=-1;
  for (int i=0;i<N;i++)
  {
    if (REC__Scintillator_detector->getValue(scintillatorMap[k][i])==detectorType["FTOF"])
    {
      if (REC__Scintillator_layer->getValue(scintillatorMap[k][i])==layerType["FTOF_1A"] )
      {
	index = scintillatorMap[k][i];
	break;
      }
      else if(REC__Scintillator_layer->getValue(scintillatorMap[k][i])==layerType["FTOF_1B"])
      {
	index = scintillatorMap[k][i];
	break;
      }
      else if(REC__Scintillator_layer->getValue(scintillatorMap[k][i])==layerType["FTOF_2"])
      {
	index = scintillatorMap[k][i];
	break;
      }
    }
    
  }
  if (index>=0)
    return REC__Scintillator_sector->getValue(index);
  else 
    return -111;
}




inline Float_t TIdentificatorCLAS12::RICH_HAD_X(Int_t k)
{
  
  Int_t N = richHadPartMap[k].size();
  Int_t index=-1;
  if (N>0) index=richHadPartMap[k][0];
  if (N>1) std::cout<<"Warning, more than one hit per event!! taking the dirst one. hadrons\n";
  
  if (index>=0)
    return RICH__hadrons_traced_hitx->getValue(index);
  else 
    return -1111;
}

inline Float_t TIdentificatorCLAS12::RICH_HAD_Y(Int_t k)
{
  
  Int_t N = richHadPartMap[k].size();
  Int_t index=-1;
  if (N>0) index=richHadPartMap[k][0];
  if (N>1) std::cout<<"Warning, more than one hit per event!! taking the dirst one. hadrons\n";
  
  if (index>=0)
    return RICH__hadrons_traced_hity->getValue(index);
  else 
    return -1111;
}

inline Float_t TIdentificatorCLAS12::RICH_HAD_Z(Int_t k)
{
  
  Int_t N = richHadPartMap[k].size();
  Int_t index=-1;
  if (N>0) index=richHadPartMap[k][0];
  if (N>1) std::cout<<"Warning, more than one hit per event!! taking the dirst one. hadrons\n";
  
  if (index>=0)
    return RICH__hadrons_traced_hitz->getValue(index);
  else 
    return -1111;
}

inline Float_t TIdentificatorCLAS12::RICH_HAD_T(Int_t k)
{
  
  Int_t N = richHadPartMap[k].size();
  Int_t index=-1;
  if (N>0) index=richHadPartMap[k][0];
  if (N>1) std::cout<<"Warning, more than one hit per event!! taking the first one. hadrons\n";
  
  if (index>=0)
    return RICH__hadrons_traced_time->getValue(index);
  else 
    return -1111;
}


inline Float_t TIdentificatorCLAS12::RICH_CLUSTER_X(Int_t k)
{
  
  Int_t N = richHadPartMap[k].size();
  Int_t index=-1;
  if (N>0)
  {
    index=richHadClusterMap[richHadPartMap[k][0]][0];
  }
  if (N>1) std::cout<<"Warning, more than one hit per event!! taking the first one. NEWCLUSTER\n";
  
  if (index>=0)
    return RICH__newclusters_x->getValue(index);
  else 
    return -1111;
}

inline Float_t TIdentificatorCLAS12::RICH_CLUSTER_Y(Int_t k)
{
  
  Int_t N = richHadPartMap[k].size();
  Int_t index=-1;
  if (N>0)
  {
    index=richHadClusterMap[richHadPartMap[k][0]][0];
  }
  if (N>1) std::cout<<"Warning, more than one hit per event!! taking the first one. NEWCLUSTER\n";
  
  if (index>=0)
    return RICH__newclusters_y->getValue(index);
  else 
    return -1111;
}

inline Float_t TIdentificatorCLAS12::RICH_CLUSTER_Z(Int_t k)
{
  
  Int_t N = richHadPartMap[k].size();
  Int_t index=-1;
  if (N>0)
  {
    index=richHadClusterMap[richHadPartMap[k][0]][0];
  }
  if (N>1) std::cout<<"Warning, more than one hit per event!! taking the first one. NEWCLUSTER\n";
  
  if (index>=0)
    return RICH__newclusters_z->getValue(index);
  else 
    return -1111;
}

inline Float_t TIdentificatorCLAS12::RICH_CLUSTER_T(Int_t k)
{
  
  Int_t N = richHadPartMap[k].size();
  Int_t index=-1;
  if (N>0)
  {
    index=richHadClusterMap[richHadPartMap[k][0]][0];
  }
  if (N>1) std::cout<<"Warning, more than one hit per event!! taking the first one. NEWCLUSTER\n";
  
  if (index>=0)
    return RICH__newclusters_time->getValue(index);
  else 
    return -1111;
}

///// CLUSTER WX,WY,WZ,WT
inline Float_t TIdentificatorCLAS12::RICH_CLUSTER_WX(Int_t k)
{
  
  Int_t N = richHadPartMap[k].size();
  Int_t index=-1;
  if (N>0)
  {
    index=richHadClusterMap[richHadPartMap[k][0]][0];
  }
  if (N>1) std::cout<<"Warning, more than one hit per event!! taking the first one. NEWCLUSTER\n";
  
  if (index>=0)
    return RICH__newclusters_wx->getValue(index);
  else 
    return -1111;
}

inline Float_t TIdentificatorCLAS12::RICH_CLUSTER_WY(Int_t k)
{
  
  Int_t N = richHadPartMap[k].size();
  Int_t index=-1;
  if (N>0)
  {
    index=richHadClusterMap[richHadPartMap[k][0]][0];
  }
  if (N>1) std::cout<<"Warning, more than one hit per event!! taking the first one. NEWCLUSTER\n";
  
  if (index>=0)
    return RICH__newclusters_wy->getValue(index);
  else 
    return -1111;
}

inline Float_t TIdentificatorCLAS12::RICH_CLUSTER_WZ(Int_t k)
{
  
  Int_t N = richHadPartMap[k].size();
  Int_t index=-1;
  if (N>0)
  {
    index=richHadClusterMap[richHadPartMap[k][0]][0];
  }
  if (N>1) std::cout<<"Warning, more than one hit per event!! taking the first one. NEWCLUSTER\n";
  
  if (index>=0)
    return RICH__newclusters_wz->getValue(index);
  else 
    return -1111;
}

inline Float_t TIdentificatorCLAS12::RICH_CLUSTER_WT(Int_t k)
{
  
  Int_t N = richHadPartMap[k].size();
  Int_t index=-1;
  if (N>0)
  {
    index=richHadClusterMap[richHadPartMap[k][0]][0];
  }
  if (N>1) std::cout<<"Warning, more than one hit per event!! taking the first one. NEWCLUSTER\n";
  
  if (index>=0)
    return RICH__newclusters_wtime->getValue(index);
  else 
    return -1111;
}

/////////

inline Float_t TIdentificatorCLAS12::RICH_RR_X(Int_t k)
{
  // REC::RICH
  
  Int_t N = richRRPartMap[k].size();
  Int_t index=-1;
  if (N>0) index=richRRPartMap[k][0];
  if (N>1) std::cout<<"Warning, more than one hit per event!! taking the first one. REC::RICH X\n";
  
  if (index>=0)
    return REC__RICH_x->getValue(index);
  else 
    return -1111;
}

inline Float_t TIdentificatorCLAS12::RICH_RR_Y(Int_t k)
{
  // REC::RICH
  
  Int_t N = richRRPartMap[k].size();
  Int_t index=-1;
  if (N>0) index=richRRPartMap[k][0];
  if (N>1) std::cout<<"Warning, more than one hit per event!! taking the first one. REC::RICH Y\n";
  
  if (index>=0)
    return REC__RICH_y->getValue(index);
  else 
    return -1111;
}

inline Float_t TIdentificatorCLAS12::RICH_RR_Z(Int_t k)
{
  // REC::RICH
  
  Int_t N = richRRPartMap[k].size();
  Int_t index=-1;
  if (N>0) index=richRRPartMap[k][0];
  if (N>1) std::cout<<"Warning, more than one hit per event!! taking the first one. REC::RICH Z\n";
  
  if (index>=0)
    return REC__RICH_z->getValue(index);
  else 
    return -1111;
}


inline Float_t TIdentificatorCLAS12::RICH_RR_HX(Int_t k)
{
  // REC::RICH
  
  Int_t N = richRRPartMap[k].size();
  Int_t index=-1;
  if (N>0) index=richRRPartMap[k][0];
  if (N>1) std::cout<<"Warning, more than one hit per event!! taking the first one. REC::RICH HX\n";
  
  if (index>=0)
    return REC__RICH_hx->getValue(index);
  else 
    return -1111;
}

inline Float_t TIdentificatorCLAS12::RICH_RR_HY(Int_t k)
{
  // REC::RICH
  
  Int_t N = richRRPartMap[k].size();
  Int_t index=-1;
  if (N>0) index=richRRPartMap[k][0];
  if (N>1) std::cout<<"Warning, more than one hit per event!! taking the first one. REC::RICH HY\n";
  
  if (index>=0)
    return REC__RICH_hy->getValue(index);
  else 
    return -1111;
}

inline Float_t TIdentificatorCLAS12::RICH_RR_HZ(Int_t k)
{
  // REC::RICH
  
  Int_t N = richRRPartMap[k].size();
  Int_t index=-1;
  if (N>0) index=richRRPartMap[k][0];
  if (N>1) std::cout<<"Warning, more than one hit per event!! taking the first one. REC::RICH HZ\n";
  
  if (index>=0)
    return REC__RICH_hz->getValue(index);
  else 
    return -1111;
}


inline Float_t TIdentificatorCLAS12::RICH_PMT(Int_t k)
{
  
  Int_t N = richRRPartMap[k].size();
  Int_t index=-1;
  if (N>0)
  {
    int RRrow = richRRPartMap[k][0];// Map
    int CLrow = richRRClusterMap[RRrow][0]; //rev Map
    index=richHitClusterMap[CLrow-1][0]; // Map //hit id = row+1
  }
  if (N>1) std::cout<<"Warning, more than one hit per event!! taking the first one. NEWCLUSTER\n";
  
  if (index>=0)
    return RICH__newhits_pmt->getValue(index);
  else 
    return -1111;
}

inline Float_t TIdentificatorCLAS12::RICH_ANODE(Int_t k)
{
  
  Int_t N = richRRPartMap[k].size();
  Int_t index=-1;
  if (N>0)
  {
    int RRrow = richRRPartMap[k][0];// Map
    int CLrow = richRRClusterMap[RRrow][0]; //rev Map
    index=richHitClusterMap[CLrow-1][0]; // Map //hit id = row+1
  }
  if (N>1) std::cout<<"Warning, more than one hit per event!! taking the first one. NEWCLUSTER\n";
  
  if (index>=0)
    return RICH__newhits_anode->getValue(index);
  else 
    return -1111;
}

inline Float_t TIdentificatorCLAS12::RICH_HITS_X(Int_t k)
{
  
  Int_t N = richRRPartMap[k].size();
  Int_t index=-1;
  if (N>0)
  {
    int RRrow = richRRPartMap[k][0];// Map
    int CLrow = richRRClusterMap[RRrow][0]; //rev Map
    index=richHitClusterMap[CLrow-1][0]; // Map //hit id = row+1
  }
  if (N>1) std::cout<<"Warning, more than one hit per event!! taking the first one. NEWCLUSTER\n";
  
  if (index>=0)
    return RICH__newhits_x->getValue(index);
  else 
    return -1111;
}

inline Float_t TIdentificatorCLAS12::RICH_HITS_Y(Int_t k)
{
  
  Int_t N = richRRPartMap[k].size();
  Int_t index=-1;
  if (N>0)
  {
    int RRrow = richRRPartMap[k][0];// Map
    int CLrow = richRRClusterMap[RRrow][0]; //rev Map
    index=richHitClusterMap[CLrow-1][0]; // Map //hit id = row+1
  }
  if (N>1) std::cout<<"Warning, more than one hit per event!! taking the first one. NEWCLUSTER\n";
  
  if (index>=0)
    return RICH__newhits_y->getValue(index);
  else 
    return -1111;
}

inline Float_t TIdentificatorCLAS12::RICH_HITS_Z(Int_t k)
{
  
  Int_t N = richRRPartMap[k].size();
  Int_t index=-1;
  if (N>0)
  {
    int RRrow = richRRPartMap[k][0];// Map
    int CLrow = richRRClusterMap[RRrow][0]; //rev Map
    index=richHitClusterMap[CLrow-1][0]; // Map //cluster id = row+1
  }
  if (N>1) std::cout<<"Warning, more than one hit per event!! taking the first one. NEWCLUSTER\n";
  
  if (index>=0)
    return RICH__newhits_z->getValue(index);
  else 
    return -1111;
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
