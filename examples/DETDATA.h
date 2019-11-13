#ifndef __DETDATA__
#define __DETDATA__
#define MAXPART_MIX 6
#define DEFAULTVALUE -111111
#include "TObject.h"
////////////// detector data sucture ///////////
// If you add a variable here, add it also to the Particle Class.
class DETDATA: public TObject {
 public:
  Float_t beta[MAXPART_MIX];
  Float_t m2b[MAXPART_MIX];
  Float_t vx[MAXPART_MIX];
  Float_t vy[MAXPART_MIX];
  Float_t vz[MAXPART_MIX];
  Float_t dcx[MAXPART_MIX];
  Float_t dcy[MAXPART_MIX];
  Float_t dcz[MAXPART_MIX];
  Float_t statPart[MAXPART_MIX];
  Float_t dc_chi2[MAXPART_MIX];
  Float_t dc_ndf[MAXPART_MIX];
  Float_t pcal_lu[MAXPART_MIX];
  Float_t pcal_lv[MAXPART_MIX];
  Float_t pcal_lw[MAXPART_MIX];

  DETDATA(){
    for (int k = 0; k<MAXPART_MIX; k++){
      beta[k] = DEFAULTVALUE;
      m2b[k] = DEFAULTVALUE;
      vx[k] = DEFAULTVALUE;
      vy[k] = DEFAULTVALUE;
      vz[k] = DEFAULTVALUE;
      dcx[k] = DEFAULTVALUE;
      dcy[k] = DEFAULTVALUE;
      dcz[k] = DEFAULTVALUE;
      statPart[k] = DEFAULTVALUE;
      dc_chi2[k] = DEFAULTVALUE;
      dc_ndf[k] = DEFAULTVALUE;
      pcal_lu[k] = DEFAULTVALUE;
      pcal_lv[k] = DEFAULTVALUE;
      pcal_lw[k] = DEFAULTVALUE;
    }
  };
  
  ~DETDATA(){ };

  ClassDef(DETDATA,1);
};
//////////////////////////////////////////
#endif
