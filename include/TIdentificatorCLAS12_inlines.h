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
  get_REC__Particle(k);
  return REC__Particle_beta;
}

inline Double_t TIdentificatorCLAS12::Charge(Int_t k, Bool_t kind)
{
    // Return the electrical charge for the particle in the row k of the EVNT
    // bank.

  if (kind == 0) {
    get_REC__Particle(k);
    return REC__Particle_charge;
  } else {
    return -111; // Not found in MC__xxx banks
  }
}



inline Double_t TIdentificatorCLAS12::Px(Int_t k, Bool_t kind)
{

    if (kind == 0) {
      get_REC__Particle(k);
      return REC__Particle_px;
    } else {
      get_MC__Lund(k);
      return MC__Lund_px;
    }
}



inline Double_t TIdentificatorCLAS12::Py(Int_t k, Bool_t kind)
{
    if (kind == 0) {
      get_REC__Particle(k);
      return REC__Particle_py;
    } else {
      get_MC__Lund(k);
      return MC__Lund_py;
    }
}



inline Double_t TIdentificatorCLAS12::Pz(Int_t k, Bool_t kind)
{
    if (kind == 0) {
      get_REC__Particle(k);
      return REC__Particle_pz;
    } else {
      get_MC__Lund(k);
      return MC__Lund_pz;
    }
}



inline Double_t TIdentificatorCLAS12::X(Int_t k, Bool_t kind)
{
    if (kind == 0) {
      get_REC__Particle(k);
      return REC__Particle_vx;
    } else {
      get_MC__Lund(k);
      return MC__Lund_vx;
    }
}


inline Double_t TIdentificatorCLAS12::Y(Int_t k, Bool_t kind)
{
    if (kind == 0) {
      get_REC__Particle(k);
      return REC__Particle_vy;
    } else {
      get_MC__Lund(k);
      return MC__Lund_vy;
    }
}


inline Double_t TIdentificatorCLAS12::Z(Int_t k, Bool_t kind)
{
    if (kind == 0) {
      get_REC__Particle(k);
      return REC__Particle_vz;
    } else {
      get_MC__Lund(k);
      return MC__Lund_vz;
    }
}

inline Int_t TIdentificatorCLAS12::Helic() /// from event
{
  get_REC__Event(0);
  return (Int_t)REC__Event_helicity;
}

inline Int_t TIdentificatorCLAS12::HelicOnline() /// from event
{
  get_HEL__online(0);
  return HEL__online_helicity;
}

inline Int_t TIdentificatorCLAS12::HelicOnlineRaw() /// from event
{
  get_HEL__online(0);
  return HEL__online_helicityRaw;
}


inline Int_t TIdentificatorCLAS12::HelicFlip() /// from event
{
  get_HEL__flip(0);
  return HEL__flip_helicity;
}

inline Int_t TIdentificatorCLAS12::HelicFlipRaw() /// from event
{
  get_HEL__flip(0);
  return HEL__flip_helicityRaw;
}

inline Int_t TIdentificatorCLAS12::HelicFlipEvent() /// from event
{
  get_HEL__flip(0);
  return HEL__flip_event;
}


inline Int_t TIdentificatorCLAS12::Event() /// from event
{
  get_RUN__config(0);
  return (Int_t)RUN__config_event;
}

inline Int_t TIdentificatorCLAS12::RunN() /// from event
{
  get_RUN__config(0);
  return (Int_t)RUN__config_run;
}

inline Float_t TIdentificatorCLAS12::Torus() /// from event
{
  get_RUN__config(0);
  return RUN__config_torus;
}

inline Float_t TIdentificatorCLAS12::Solenoid() /// from event
{
  get_RUN__config(0);
  return RUN__config_solenoid;
}

inline Float_t TIdentificatorCLAS12::STTime(Int_t k) /// from event
{
  get_REC__Event(k);
  return REC__Event_startTime;
}

inline Float_t TIdentificatorCLAS12::RFTime(Int_t k) /// from event
{
  get_REC__Event(k);
  return REC__Event_RFTime;
}


inline Int_t TIdentificatorCLAS12::Status(Int_t k)
{
  get_REC__Particle(k);
  return (Int_t)REC__Particle_status;
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
    get_REC__Cherenkov(cherenkovMap[k][i]);
    if (REC__Cherenkov_detector == detectorType["LTCC"])
    {
      index=i;
      break;
    }
    
  }
  if (index>=0){
    get_REC__Cherenkov(index);
    return REC__Cherenkov_nphe;
  }
  else
    return -111;
}

inline Double_t TIdentificatorCLAS12::NpheHTCC(Int_t k)
{
  int N = cherenkovMap[k].size();
  int index=-1;
  for (int i=0;i<N;i++) // get index of the correct detector
  {
    get_REC__Cherenkov(cherenkovMap[k][i]);
    if (REC__Cherenkov_detector == detectorType["HTCC"])
    {
      index=i;
      break;
    }
    
  }
  if (index>=0){
    get_REC__Cherenkov(index);
    return REC__Cherenkov_nphe;
  }
  else
    return -111;
}



inline Double_t TIdentificatorCLAS12::Chi2pid(Int_t k)
{
  get_REC__Particle(k);
  return REC__Particle_chi2pid;
}

inline Int_t TIdentificatorCLAS12::GetNRows()
{

  return  REC__Particle->getRows();
}

inline Int_t TIdentificatorCLAS12::GetMCNRows()
{
  return  MC__Lund->getRows();
}

inline Float_t TIdentificatorCLAS12::MCMass(Int_t k)
{
  get_MC__Lund(k);
  return MC__Lund_mass;
}

inline Float_t TIdentificatorCLAS12::LundType(Int_t k)
{
  get_MC__Lund(k);
  return  MC__Lund_type;
}

inline Int_t TIdentificatorCLAS12::GetNPart(Int_t pid, Bool_t kind)
{
  int Npart = 0;
  if (!kind)
  {
    int n = REC__Particle->getRows();
    for (int k=0;k<n;k++){
      get_REC__Particle(k);
      if (REC__Particle_pid==pid)
	Npart++;
    }
  }
  else
  {
    int n = MC__Lund->getRows();
    for (int k=kIndLundFirst;k<n;k++){
      get_MC__Lund(k);
      if (MC__Lund_pid == pid && MC__Lund_type == 1)
	Npart++;
    }
  }
  
  return Npart;
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



*/

inline Double_t TIdentificatorCLAS12::DCChi2(Int_t k)
{

  Int_t N = trackMap[k].size();
  Int_t index=-1;
  for (int i=0;i<N;i++)
  {
    get_REC__Track(trackMap[k][i]);
    if (REC__Track_detector==detectorType["DC"])
    {
      index = trackMap[k][i];
      break;
    }
  }
  if (index>=0){
    get_REC__Track(index);
    return REC__Track_chi2;
  }
  else 
    return -1111;

}


inline Double_t TIdentificatorCLAS12::Etot(Int_t k,Bool_t kind)
{
    // Return total energy deposited in the calorimeter for the particle in
    // the row k of the EVNT bank.
  if (kind==0){
    return Ein(k) + Eout(k) + Epcal(k);
  }
  else{
    get_MC__Lund(k);
    return MC__Lund_energy;
  }
}



inline Double_t TIdentificatorCLAS12::Ein(Int_t k)
{
  
  Int_t N = calorimeterMap[k].size();
  Int_t index=-1;
  for (int i=0;i<N;i++)
  {
    get_REC__Calorimeter(calorimeterMap[k][i]);
    if (REC__Calorimeter_detector==detectorType["ECAL"]
	&&REC__Calorimeter_layer==layerType["EC_Inner"] )
    {
      index = calorimeterMap[k][i];
      break;
    }
  }
  if (index>=0){
    get_REC__Calorimeter(index);
    return REC__Calorimeter_energy;
  }
  else 
    return 0;
}



inline Double_t TIdentificatorCLAS12::Eout(Int_t k)
{
  
  Int_t N = calorimeterMap[k].size();
  Int_t index=-1;
  for (int i=0;i<N;i++)
  {
    get_REC__Calorimeter(calorimeterMap[k][i]);
    if (REC__Calorimeter_detector==detectorType["ECAL"]
	&&REC__Calorimeter_layer==layerType["EC_Outer"] )
    {
      index = calorimeterMap[k][i];
      break;
    }
  }
  if (index>=0){
    get_REC__Calorimeter(index);
    return REC__Calorimeter_energy;
  }
  else 
    return 0;
}


inline Double_t TIdentificatorCLAS12::Epcal(Int_t k)
{
  
  Int_t N = calorimeterMap[k].size();
  Int_t index=-1;
  for (int i=0;i<N;i++)
  {
    get_REC__Calorimeter(calorimeterMap[k][i]);
    if (REC__Calorimeter_detector==detectorType["ECAL"]
	&&REC__Calorimeter_layer==layerType["PCAL"] )
    {
      index = calorimeterMap[k][i];
      break;
    }
  }
  if (index>=0){
    get_REC__Calorimeter(index);
    return REC__Calorimeter_energy;
  }
  else 
    return 0;
}


inline Int_t TIdentificatorCLAS12::SectorLTCC(Int_t k)
{
  int N = cherenkovMap[k].size();
  int index=-1;
  for (int i=0;i<N;i++) // get index of the correct detector
  {
    get_REC__Cherenkov(cherenkovMap[k][i]);
    if (REC__Cherenkov_detector == detectorType["LTCC"])
    {
      index=i;
      break;
    }
    
  }
  if (index>=0){
    get_REC__Cherenkov(index);
    return REC__Cherenkov_sector;
  }
  else
    return -111;
}


inline Int_t TIdentificatorCLAS12::SectorHTCC(Int_t k)
{
  int N = cherenkovMap[k].size();
  int index=-1;
  for (int i=0;i<N;i++) // get index of the correct detector
  {
    get_REC__Cherenkov(cherenkovMap[k][i]);
    if (REC__Cherenkov_detector == detectorType["HTCC"])
    {
      index=i;
      break;
    }
    
  }
  if (index>=0){
    get_REC__Cherenkov(index);
    return REC__Cherenkov_sector;
  }
  else
    return -111;
}

inline Int_t TIdentificatorCLAS12::SectorECAL(Int_t k)
{
  
  Int_t N = calorimeterMap[k].size();
  Int_t index=-1;
  for (int i=0;i<N;i++)
  {
    get_REC__Calorimeter(calorimeterMap[k][i]);
    if (REC__Calorimeter_detector==detectorType["ECAL"])
    {
      index = calorimeterMap[k][i];
      break;
    }
  }
  if (index>=0){
    get_REC__Calorimeter(index);
    return REC__Calorimeter_sector;
  }
  else 
    return 0;
}

   
inline Float_t TIdentificatorCLAS12::LU_PCAL(Int_t k)
{
  
  Int_t N = calorimeterMap[k].size();
  Int_t index=-1;
  for (int i=0;i<N;i++)
  {
    get_REC__Calorimeter(calorimeterMap[k][i]);
    if (REC__Calorimeter_detector==detectorType["ECAL"]
	&&REC__Calorimeter_layer==layerType["PCAL"] )
    {
      index = calorimeterMap[k][i];
      break;
    }
  }
  if (index>=0){
    get_REC__Calorimeter(index);
    return REC__Calorimeter_lu;
  }
  else 
    return -111;
}

inline Float_t TIdentificatorCLAS12::LV_PCAL(Int_t k)
{
  
  Int_t N = calorimeterMap[k].size();
  Int_t index=-1;
  for (int i=0;i<N;i++)
  {
    get_REC__Calorimeter(calorimeterMap[k][i]);
    if (REC__Calorimeter_detector==detectorType["ECAL"]
	&&REC__Calorimeter_layer==layerType["PCAL"] )
    {
      index = calorimeterMap[k][i];
      break;
    }
  }
  if (index>=0){
    get_REC__Calorimeter(index);
    return REC__Calorimeter_lv;
  }
  else 
    return -111;
}

inline Float_t TIdentificatorCLAS12::LW_PCAL(Int_t k)
{
  
  Int_t N = calorimeterMap[k].size();
  Int_t index=-1;
  for (int i=0;i<N;i++)
  {
    get_REC__Calorimeter(calorimeterMap[k][i]);
    if (REC__Calorimeter_detector==detectorType["ECAL"]
	&&REC__Calorimeter_layer==layerType["PCAL"] )
    {
      index = calorimeterMap[k][i];
      break;
    }
  }
  if (index>=0){
    get_REC__Calorimeter(index);
    return REC__Calorimeter_lw;
  }
  else 
    return -111;
}


   
inline Float_t TIdentificatorCLAS12::LU_ECIN(Int_t k)
{
  
  Int_t N = calorimeterMap[k].size();
  Int_t index=-1;
  for (int i=0;i<N;i++)
  {
    get_REC__Calorimeter(calorimeterMap[k][i]);
    if (REC__Calorimeter_detector==detectorType["ECAL"]
	&&REC__Calorimeter_layer==layerType["EC_Inner"] )
    {
      index = calorimeterMap[k][i];
      break;
    }
  }
  if (index>=0){
    get_REC__Calorimeter(index);
    return REC__Calorimeter_lu;
  }
  else 
    return -111;
}

inline Float_t TIdentificatorCLAS12::LV_ECIN(Int_t k)
{
  
  Int_t N = calorimeterMap[k].size();
  Int_t index=-1;
  for (int i=0;i<N;i++)
  {
    get_REC__Calorimeter(calorimeterMap[k][i]);
    if (REC__Calorimeter_detector==detectorType["ECAL"]
	&&REC__Calorimeter_layer==layerType["EC_Inner"] )
    {
      index = calorimeterMap[k][i];
      break;
    }
  }
  if (index>=0){
    get_REC__Calorimeter(index);
    return REC__Calorimeter_lv;
  }
  else 
    return -111;
}

inline Float_t TIdentificatorCLAS12::LW_ECIN(Int_t k)
{
  
  Int_t N = calorimeterMap[k].size();
  Int_t index=-1;
  for (int i=0;i<N;i++)
  {
    get_REC__Calorimeter(calorimeterMap[k][i]);
    if (REC__Calorimeter_detector==detectorType["ECAL"]
	&&REC__Calorimeter_layer==layerType["EC_Inner"] )
    {
      index = calorimeterMap[k][i];
      break;
    }
  }
  if (index>=0){
    get_REC__Calorimeter(index);
    return REC__Calorimeter_lw;
  }
  else 
    return -111;
}

   
inline Float_t TIdentificatorCLAS12::LU_ECOUT(Int_t k)
{
  
  Int_t N = calorimeterMap[k].size();
  Int_t index=-1;
  for (int i=0;i<N;i++)
  {
    get_REC__Calorimeter(calorimeterMap[k][i]);
    if (REC__Calorimeter_detector==detectorType["ECAL"]
	&&REC__Calorimeter_layer==layerType["EC_Outer"] )
    {
      index = calorimeterMap[k][i];
      break;
    }
  }
  if (index>=0){
    get_REC__Calorimeter(index);
    return REC__Calorimeter_lu;
  }
  else 
    return -111;
}

inline Float_t TIdentificatorCLAS12::LV_ECOUT(Int_t k)
{
  
  Int_t N = calorimeterMap[k].size();
  Int_t index=-1;
  for (int i=0;i<N;i++)
  {
    get_REC__Calorimeter(calorimeterMap[k][i]);
    if (REC__Calorimeter_detector==detectorType["ECAL"]
	&&REC__Calorimeter_layer==layerType["EC_Outer"] )
    {
      index = calorimeterMap[k][i];
      break;
    }
  }
  if (index>=0){
    get_REC__Calorimeter(index);
    return REC__Calorimeter_lv;
  }
  else 
    return -111;
}

inline Float_t TIdentificatorCLAS12::LW_ECOUT(Int_t k)
{
  
  Int_t N = calorimeterMap[k].size();
  Int_t index=-1;
  for (int i=0;i<N;i++)
  {
    get_REC__Calorimeter(calorimeterMap[k][i]);
    if (REC__Calorimeter_detector==detectorType["ECAL"]
	&&REC__Calorimeter_layer==layerType["EC_Outer"] )
    {
      index = calorimeterMap[k][i];
      break;
    }
  }
  if (index>=0){
    get_REC__Calorimeter(index);
    return REC__Calorimeter_lw;
  }
  else 
    return -111;
}



inline Float_t TIdentificatorCLAS12::HX_PCAL(Int_t k)
{
  
  Int_t N = calorimeterMap[k].size();
  Int_t index=-1;
  for (int i=0;i<N;i++)
  {
    get_REC__Calorimeter(calorimeterMap[k][i]);
    if (REC__Calorimeter_detector==detectorType["ECAL"]
	&&REC__Calorimeter_layer==layerType["PCAL"] )
    {
      index = calorimeterMap[k][i];
      break;
    }
  }
  if (index>=0){
    get_REC__Calorimeter(index);
    return REC__Calorimeter_hx;
  }
  else 
    return -111;
}

inline Float_t TIdentificatorCLAS12::HY_PCAL(Int_t k)
{
  
  Int_t N = calorimeterMap[k].size();
  Int_t index=-1;
  for (int i=0;i<N;i++)
  {
    get_REC__Calorimeter(calorimeterMap[k][i]);
    if (REC__Calorimeter_detector==detectorType["ECAL"]
	&&REC__Calorimeter_layer==layerType["PCAL"] )
    {
      index = calorimeterMap[k][i];
      break;
    }
  }
  if (index>=0){
    get_REC__Calorimeter(index);
    return REC__Calorimeter_hy;
  }
  else 
    return -111;
}

inline Float_t TIdentificatorCLAS12::HZ_PCAL(Int_t k)
{
  
  Int_t N = calorimeterMap[k].size();
  Int_t index=-1;
  for (int i=0;i<N;i++)
  {
    get_REC__Calorimeter(calorimeterMap[k][i]);
    if (REC__Calorimeter_detector==detectorType["ECAL"]
	&&REC__Calorimeter_layer==layerType["PCAL"] )
    {
      index = calorimeterMap[k][i];
      break;
    }
  }
  if (index>=0){
    get_REC__Calorimeter(index);
    return REC__Calorimeter_hz;
  }
  else 
    return -111;
}

inline Float_t TIdentificatorCLAS12::HX_ECIN(Int_t k)
{
  Int_t N = calorimeterMap[k].size();
  Int_t index=-1;
  for (int i=0;i<N;i++)
  {
    get_REC__Calorimeter(calorimeterMap[k][i]);
    if (REC__Calorimeter_detector==detectorType["ECAL"]
	&&REC__Calorimeter_layer==layerType["EC_Inner"] )
    {
      index = calorimeterMap[k][i];
      break;
    }
  }
  if (index>=0){
    get_REC__Calorimeter(index);
    return REC__Calorimeter_hx;
  }
  else 
    return -111;
}

inline Float_t TIdentificatorCLAS12::HY_ECIN(Int_t k)
{
  
  Int_t N = calorimeterMap[k].size();
  Int_t index=-1;
  for (int i=0;i<N;i++)
  {
    get_REC__Calorimeter(calorimeterMap[k][i]);
    if (REC__Calorimeter_detector==detectorType["ECAL"]
	&&REC__Calorimeter_layer==layerType["EC_Inner"] )
    {
      index = calorimeterMap[k][i];
      break;
    }
  }
  if (index>=0){
    get_REC__Calorimeter(index);
    return REC__Calorimeter_hy;
  }
  else 
    return -111;
}

inline Float_t TIdentificatorCLAS12::HZ_ECIN(Int_t k)
{
  
  Int_t N = calorimeterMap[k].size();
  Int_t index=-1;
  for (int i=0;i<N;i++)
  {
    get_REC__Calorimeter(calorimeterMap[k][i]);
    if (REC__Calorimeter_detector==detectorType["ECAL"]
	&&REC__Calorimeter_layer==layerType["EC_Inner"] )
    {
      index = calorimeterMap[k][i];
      break;
    }
  }
  if (index>=0){
    get_REC__Calorimeter(index);
    return REC__Calorimeter_hz;
  }
  else 
    return -111;
}

inline Float_t TIdentificatorCLAS12::HX_ECOUT(Int_t k)
{
  
  Int_t N = calorimeterMap[k].size();
  Int_t index=-1;
  for (int i=0;i<N;i++)
  {
    get_REC__Calorimeter(calorimeterMap[k][i]);
    if (REC__Calorimeter_detector==detectorType["ECAL"]
	&&REC__Calorimeter_layer==layerType["EC_Outer"] )
    {
      index = calorimeterMap[k][i];
      break;
    }
  }
  if (index>=0){
    get_REC__Calorimeter(index);
    return REC__Calorimeter_hx;
  }
  else 
    return -111;
}

inline Float_t TIdentificatorCLAS12::HY_ECOUT(Int_t k)
{
  
  Int_t N = calorimeterMap[k].size();
  Int_t index=-1;
  for (int i=0;i<N;i++)
  {
    get_REC__Calorimeter(calorimeterMap[k][i]);
    if (REC__Calorimeter_detector==detectorType["ECAL"]
	&&REC__Calorimeter_layer==layerType["EC_Outer"] )
    {
      index = calorimeterMap[k][i];
      break;
    }
  }
  if (index>=0){
    get_REC__Calorimeter(index);
    return REC__Calorimeter_hy;
  }
  else 
    return -111;
}

inline Float_t TIdentificatorCLAS12::HZ_ECOUT(Int_t k)
{
  
  Int_t N = calorimeterMap[k].size();
  Int_t index=-1;
  for (int i=0;i<N;i++)
  {
    get_REC__Calorimeter(calorimeterMap[k][i]);
    if (REC__Calorimeter_detector==detectorType["ECAL"]
	&&REC__Calorimeter_layer==layerType["EC_Outer"] )
    {
      index = calorimeterMap[k][i];
      break;
    }
  }
  if (index>=0){
    get_REC__Calorimeter(index);
    return REC__Calorimeter_hz;
  }
  else 
    return -111;
}


/*
inline Float_t TIdentificatorCLAS12::VX_DC(Int_t k)
{
  
  Int_t N = trackMap[k].size();
  Int_t index=-1;
  for (int i=0;i<N;i++)
  {
    get_REC__Track(trackMap[k][i]);
    if (REC__Track_detector==detectorType["DC"])
    {
      index = trackMap[k][i];
      break;
    }
  }
  if (index>=0){
    get_REC__Track(index);
    return REC__Track_vx_nomm;
  }
  else 
    return -1111;
}

inline Float_t TIdentificatorCLAS12::VY_DC(Int_t k)
{
  
  Int_t N = trackMap[k].size();
  Int_t index=-1;
  for (int i=0;i<N;i++)
  {
    get_REC__Track(trackMap[k][i]);
    if (REC__Track_detector==detectorType["DC"])
    {
      index = trackMap[k][i];
      break;
    }
  }
  if (index>=0){
    get_REC__Track(index);
    return REC__Track_vy_nomm;
  }
  else 
    return -1111;
}

inline Float_t TIdentificatorCLAS12::VZ_DC(Int_t k)
{
  
  Int_t N = trackMap[k].size();
  Int_t index=-1;
  for (int i=0;i<N;i++)
  {
    get_REC__Track(trackMap[k][i]);
    if (REC__Track_detector==detectorType["DC"])
    {
      index = trackMap[k][i];
      break;
    }
  }
  if (index>=0){
    get_REC__Track(index);
    return REC__Track_vz_nomm;
  }
  else 
    return -1111;
}
*/

inline Int_t TIdentificatorCLAS12::SectorDC(Int_t k)
{
  Int_t N = trackMap[k].size();
  Int_t index=-1;
  for (int i=0;i<N;i++)
  {
    get_REC__Track(trackMap[k][i]);
    if (REC__Track_detector==detectorType["DC"])
    {
      index = trackMap[k][i];
      break;
    }
  }

  if (index>=0){
    get_REC__Track(index);
    return REC__Track_sector;
  }
  else 
    return -111;
}
/*
inline Float_t TIdentificatorCLAS12::Px_DC(Int_t k)
{
  
  Int_t N = trackMap[k].size();
  Int_t index=-1;
  for (int i=0;i<N;i++)
  {
    get_REC__Track(trackMap[k][i]);
    if (REC__Track_detector==detectorType["DC"])
    {
      index = trackMap[k][i];
      break;
    }
  }
  if (index>=0){
    get_REC__Track(index);
    return REC__Track_px_nomm;
  }
  else 
    return -111;
}

inline Float_t TIdentificatorCLAS12::Py_DC(Int_t k)
{
  
  Int_t N = trackMap[k].size();
  Int_t index=-1;
  for (int i=0;i<N;i++)
  {
    get_REC__Track(trackMap[k][i]);
    if (REC__Track_detector==detectorType["DC"])
    {
      index = trackMap[k][i];
      break;
    }
  }
  if (index>=0){
    get_REC__Track(index);
    return REC__Track_py_nomm;
  }
  else 
    return -111;
}

inline Float_t TIdentificatorCLAS12::Pz_DC(Int_t k)
{
  
  Int_t N = trackMap[k].size();
  Int_t index=-1;
  for (int i=0;i<N;i++)
  {
    get_REC__Track(trackMap[k][i]);
    if (REC__Track_detector==detectorType["DC"])
    {
      index = trackMap[k][i];
      break;
    }
  }
  if (index>=0){
    get_REC__Track(index);
    return REC__Track_pz_nomm;
  }
  else 
    return -111;
}


///////////// Trajectory bank

inline Float_t TIdentificatorCLAS12::TrajX(Int_t k,Int_t sl)
{
  
  Int_t N = trajMap[k].size();
  Int_t index=-1;
  for (int i=0;i<N;i++)
  {
    get_REC__Traj(trajMap[k][i]);
    if (REC__Traj_detId==DCSuperLayer[sl])
    {
      index = trajMap[k][i];
      break;
    }
  }
  if (index>=0){
    get_REC__Traj(index);
    return REC__Traj_x;
  }
  else 
    return -1111;
}

inline Float_t TIdentificatorCLAS12::TrajY(Int_t k,Int_t sl)
{
  
  Int_t N = trajMap[k].size();
  Int_t index=-1;
  for (int i=0;i<N;i++)
  {
    get_REC__Traj(trajMap[k][i]);
    if (REC__Traj_detId==DCSuperLayer[sl])
    {
      index = trajMap[k][i];
      break;
    }
  }
  if (index>=0){
    get_REC__Traj(index);
    return REC__Traj_y;
  }
  else 
    return -1111;
}

inline Float_t TIdentificatorCLAS12::TrajZ(Int_t k,Int_t sl)
{
  
  Int_t N = trajMap[k].size();
  Int_t index=-1;
  for (int i=0;i<N;i++)
  {
    get_REC__Traj(trajMap[k][i]);
    if (REC__Traj_detId==DCSuperLayer[sl])
    {
      index = trajMap[k][i];
      break;
    }
  }
  if (index>=0){
    get_REC__Traj(index);
    return REC__Traj_z;
  }
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
    get_REC__Traj(trajMap[k][i]);
    if (REC__Traj_detId==DCSuperLayer[sl])
    {
      index = trajMap[k][i];
      break;
    }
  }
  if (index>=0){
    get_REC__Traj(index);
    return REC__Traj_cx;
  }
  else 
    return -1111;
}

inline Float_t TIdentificatorCLAS12::TrajCY(Int_t k,Int_t sl)
{
  
  Int_t N = trajMap[k].size();
  Int_t index=-1;
  for (int i=0;i<N;i++)
  {
    get_REC__Traj(trajMap[k][i]);
    if (REC__Traj_detId==DCSuperLayer[sl])
    {
      index = trajMap[k][i];
      break;
    }
  }
  if (index>=0){
    get_REC__Traj(index);
    return REC__Traj_cy;
  }
  else 
    return -1111;
}

inline Float_t TIdentificatorCLAS12::TrajCZ(Int_t k,Int_t sl)
{
  
  Int_t N = trajMap[k].size();
  Int_t index=-1;
  for (int i=0;i<N;i++)
  {
    get_REC__Traj(trajMap[k][i]);
    if (REC__Traj_detId==DCSuperLayer[sl])
    {
      index = trajMap[k][i];
      break;
    }
  }
  if (index>=0){
    get_REC__Traj(index);
    return REC__Traj_cz;
  }
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
    get_REC__Traj(trajMap[k][i]);
    if (REC__Traj_detId==trajDetId[dname])
    {
      index = trajMap[k][i];
      break;
    }
  }
  if (index>=0){
    get_REC__Traj(index);
    return REC__Traj_x;
  }
  else 
    return -1111;
}

inline Float_t TIdentificatorCLAS12::TrajDetIdY(Int_t k,TString dname)
{
  
  Int_t N = trajMap[k].size();
  Int_t index=-1;
  for (int i=0;i<N;i++)
  {
    get_REC__Traj(trajMap[k][i]);
    if (REC__Traj_detId==trajDetId[dname])
    {
      index = trajMap[k][i];
      break;
    }
  }
  if (index>=0){
    get_REC__Traj(index);
    return REC__Traj_y;
  }
  else 
    return -1111;
}

inline Float_t TIdentificatorCLAS12::TrajDetIdZ(Int_t k,TString dname)
{
  
  Int_t N = trajMap[k].size();
  Int_t index=-1;
  for (int i=0;i<N;i++)
  {
    get_REC__Traj(trajMap[k][i]);
    if (REC__Traj_detId==trajDetId[dname])
    {
      index = trajMap[k][i];
      break;
    }
  }
  if (index>=0){
    get_REC__Traj(index);
    return REC__Traj_z;
  }
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
  return TrajDetIdX(k,"FTOF1A");
}
inline Float_t TIdentificatorCLAS12::TrajFTOF1AY(Int_t k)
{
  return TrajDetIdY(k,"FTOF1A");
}
inline Float_t TIdentificatorCLAS12::TrajFTOF1AZ(Int_t k)
{
  return TrajDetIdZ(k,"FTOF1A");
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
    get_REC__Traj(trajMap[k][i]);
    if (REC__Traj_detId==DCSuperLayer[sl] || REC__Traj_detId==DCSuperLayer[sl+1] )
    {
      index = trajMap[k][i];
      break;
    }
  }
  if (index>=0){
    get_REC__Traj(index);
    return REC__Traj_x;
  }
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
    get_REC__Traj(trajMap[k][i]);
    if (REC__Traj_detId==DCSuperLayer[sl] || REC__Traj_detId==DCSuperLayer[sl+1] )
    {
      index = trajMap[k][i];
      break;
    }
  }
  if (index>=0){
    get_REC__Traj(index);
    return REC__Traj_y;
  }
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
    get_REC__Traj(trajMap[k][i]);
    if (REC__Traj_detId==DCSuperLayer[sl] || REC__Traj_detId==DCSuperLayer[sl+1] )
    {
      index = trajMap[k][i];
      break;
    }
  }
  if (index>=0){
    get_REC__Traj(index);
    return REC__Traj_z;
  }
  else 
    return -1111;
}

*/
inline Float_t TIdentificatorCLAS12::PathTOF(Int_t k)
{
  Int_t N = scintillatorMap[k].size();
  Int_t index=-1;
  for (int i=0;i<N;i++)
  {
    get_REC__Scintillator(scintillatorMap[k][i]);
    if (REC__Scintillator_detector==detectorType["FTOF"])
    {
      if (REC__Scintillator_layer==layerType["FTOF_1A"] )
      {
	index = scintillatorMap[k][i];
	break;
      }
      else if(REC__Scintillator_layer==layerType["FTOF_1B"])
      {
	index = scintillatorMap[k][i];
	break;
      }
      else if(REC__Scintillator_layer==layerType["FTOF_2"])
      {
	index = scintillatorMap[k][i];
	break;
      }
    }
    
  }
  if (index>=0){
    get_REC__Scintillator(index);
    return REC__Scintillator_path;
  }
  else 
    return -111;
}

inline Float_t TIdentificatorCLAS12::TimeTOF(Int_t k)
{
  Int_t N = scintillatorMap[k].size();
  Int_t index=-1;
  for (int i=0;i<N;i++)
  {
    get_REC__Scintillator(scintillatorMap[k][i]);
    if (REC__Scintillator_detector == detectorType["FTOF"])
    {
      if (REC__Scintillator_layer == layerType["FTOF_1B"] )
      {
	index = scintillatorMap[k][i];
	break;
      }
      else if(REC__Scintillator_layer == layerType["FTOF_1A"])
      {
	index = scintillatorMap[k][i];
	break;
      }
      else if(REC__Scintillator_layer == layerType["FTOF_2"])
      {
	index = scintillatorMap[k][i];
	break;
      }
    }
    
  }
  if (index>=0){
    get_REC__Scintillator(index);
    return REC__Scintillator_time;
  }
  else 
    return -111;
}

inline Int_t TIdentificatorCLAS12::SectorTOF(Int_t k)
{
  Int_t N = scintillatorMap[k].size();
  Int_t index=-1;
  for (int i=0;i<N;i++)
  {
    get_REC__Scintillator(scintillatorMap[k][i]);
    if (REC__Scintillator_detector == detectorType["FTOF"])
    {
      if (REC__Scintillator_layer == layerType["FTOF_1A"] )
      {
	index = scintillatorMap[k][i];
	break;
      }
      else if(REC__Scintillator_layer == layerType["FTOF_1B"])
      {
	index = scintillatorMap[k][i];
	break;
      }
      else if(REC__Scintillator_layer == layerType["FTOF_2"])
      {
	index = scintillatorMap[k][i];
	break;
      }
    }
    
  }
  if (index>=0){
    get_REC__Scintillator(index);
    return REC__Scintillator_sector;
  }
  else 
    return -111;
}




inline Float_t TIdentificatorCLAS12::RICH_HAD_X(Int_t k)
{
  
  Int_t N = richHadPartMap[k].size();
  Int_t index=-1;
  if (N>0) index=richHadPartMap[k][0];
  if (N>1) std::cout<<"Warning, more than one hit per event!! taking the dirst one. hadrons\n";
  
  if (index>=0){
    get_RICH__hadrons(index);
    return RICH__hadrons_traced_hitx;
  }
  else 
    return -1111;
}

inline Float_t TIdentificatorCLAS12::RICH_HAD_Y(Int_t k)
{
  
  Int_t N = richHadPartMap[k].size();
  Int_t index=-1;
  if (N>0) index=richHadPartMap[k][0];
  if (N>1) std::cout<<"Warning, more than one hit per event!! taking the dirst one. hadrons\n";
  
  if (index>=0){
    get_RICH__hadrons(index);
    return RICH__hadrons_traced_hity;
  }
  else 
    return -1111;
}

inline Float_t TIdentificatorCLAS12::RICH_HAD_Z(Int_t k)
{
  
  Int_t N = richHadPartMap[k].size();
  Int_t index=-1;
  if (N>0) index=richHadPartMap[k][0];
  if (N>1) std::cout<<"Warning, more than one hit per event!! taking the dirst one. hadrons\n";
  
  if (index>=0){
    get_RICH__hadrons(index);
    return RICH__hadrons_traced_hitz;
  }
  else 
    return -1111;
}

inline Float_t TIdentificatorCLAS12::RICH_HAD_T(Int_t k)
{
  
  Int_t N = richHadPartMap[k].size();
  Int_t index=-1;
  if (N>0) index=richHadPartMap[k][0];
  if (N>1) std::cout<<"Warning, more than one hit per event!! taking the first one. hadrons\n";
  
  if (index>=0){
    get_RICH__hadrons(index);
    return RICH__hadrons_traced_time;
  }
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
  
  if (index>=0){
    get_RICH__clusters(index);
    return RICH__clusters_x;
  }
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
  
  if (index>=0){
    get_RICH__clusters(index);
    return RICH__clusters_y;
  }
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
  
  if (index>=0){
    get_RICH__clusters(index);
    return RICH__clusters_z;
  }
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
  
  if (index>=0){
    get_RICH__clusters(index);
    return RICH__clusters_time;
  }
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
  
  if (index>=0){
    get_RICH__clusters(index);
    return RICH__clusters_wx;
  }
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
  
  if (index>=0){
    get_RICH__clusters(index);
    return RICH__clusters_wy;
  }
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
  
  if (index>=0){
    get_RICH__clusters(index);
    return RICH__clusters_wz;
  }
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
  
  if (index>=0){
    get_RICH__clusters(index);
    return RICH__clusters_wtime;
  }
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
  
  if (index>=0){
    get_REC__RICH(index);
    return REC__RICH_x;
  }
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
  
  if (index>=0){
    get_REC__RICH(index);
    return REC__RICH_y;
  }
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
  
  if (index>=0){
    get_REC__RICH(index);
    return REC__RICH_z;
  }
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
  
  if (index>=0){
    get_REC__RICH(index);
    return REC__RICH_hx;
  }
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
  
  if (index>=0){
    get_REC__RICH(index);
    return REC__RICH_hy;
  }
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
  
  if (index>=0){
    get_REC__RICH(index);
    return REC__RICH_hz;
  }
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
  
  if (index>=0){
    get_RICH__hits(index);
    return RICH__hits_pmt;
  }
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
  
  if (index>=0){
    get_RICH__hits(index);
    return RICH__hits_anode;
  }
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
  
  if (index>=0){
    get_RICH__hits(index);
    return RICH__hits_x;
  }
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
  
  if (index>=0){
    get_RICH__hits(index);
    return RICH__hits_y;
  }
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
  
  if (index>=0){
    get_RICH__hits(index);
    return RICH__hits_z;
  }
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
