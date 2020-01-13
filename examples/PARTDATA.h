#ifndef __PARTDATA__
#define __PARTDATA__
#define MAXPART_MIX 6
#define DEFAULTVALUE -111111
#include "TObject.h"
/////////// particle data sucture ///////////
class PARTDATA : public TObject {
 public:
  Float_t e[MAXPART_MIX];
  Float_t px[MAXPART_MIX];
  Float_t py[MAXPART_MIX];
  Float_t pz[MAXPART_MIX];
  Float_t mpid[MAXPART_MIX];
  Float_t mind[MAXPART_MIX];
  Float_t helic002_phiH[MAXPART_MIX];
  Float_t helic005_phiH[MAXPART_MIX];
  Float_t helic010_phiH[MAXPART_MIX];
  Float_t helic020_phiH[MAXPART_MIX];
  Float_t phiHs[MAXPART_MIX];
  PARTDATA(){
    for (int k = 0; k<MAXPART_MIX; k++){
      e[k] = DEFAULTVALUE;
      px[k] = DEFAULTVALUE;
      py[k] = DEFAULTVALUE;
      pz[k] = DEFAULTVALUE;
      mpid[k] = DEFAULTVALUE;
      mind[k] = DEFAULTVALUE;      
      helic002_phiH[k] = DEFAULTVALUE;
      helic005_phiH[k] = DEFAULTVALUE;
      helic010_phiH[k] = DEFAULTVALUE;
      helic020_phiH[k] = DEFAULTVALUE;
      phiHs[k] = DEFAULTVALUE;    
    }
  };
  ~PARTDATA() {};

  ClassDef(PARTDATA,1);
};
//////////////////////////////////
#endif
