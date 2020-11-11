
#ifndef __DETDATA__
#define __DETDATA__
#define MAXPART_MIX 8
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

  Float_t dcx0[MAXPART_MIX];
  Float_t dcy0[MAXPART_MIX];
  Float_t dcz0[MAXPART_MIX];

  Float_t dcx1[MAXPART_MIX];
  Float_t dcy1[MAXPART_MIX];
  Float_t dcz1[MAXPART_MIX];

  Float_t dcx2[MAXPART_MIX];
  Float_t dcy2[MAXPART_MIX];
  Float_t dcz2[MAXPART_MIX];

  Float_t statPart[MAXPART_MIX];
  Float_t dc_chi2[MAXPART_MIX];
  Float_t dc_ndf[MAXPART_MIX];
  Float_t pcal_lu[MAXPART_MIX];
  Float_t pcal_lv[MAXPART_MIX];
  Float_t pcal_lw[MAXPART_MIX];
  Float_t chi2pid[MAXPART_MIX];

  DETDATA(){
    for (int k = 0; k<MAXPART_MIX; k++){
      beta[k] = DEFAULTVALUE;
      m2b[k] = DEFAULTVALUE;
      vx[k] = DEFAULTVALUE;
      vy[k] = DEFAULTVALUE;
      vz[k] = DEFAULTVALUE;

      dcx0[k] = DEFAULTVALUE;// region 1, superlayer 1
      dcy0[k] = DEFAULTVALUE;
      dcz0[k] = DEFAULTVALUE;

      dcx1[k] = DEFAULTVALUE;// region 2, superlayer 3
      dcy1[k] = DEFAULTVALUE;
      dcz1[k] = DEFAULTVALUE;

      dcx2[k] = DEFAULTVALUE;// region 3, superlayer 6
      dcy2[k] = DEFAULTVALUE;
      dcz2[k] = DEFAULTVALUE;

      statPart[k] = DEFAULTVALUE;
      dc_chi2[k] = DEFAULTVALUE;
      dc_ndf[k] = DEFAULTVALUE;
      pcal_lu[k] = DEFAULTVALUE;
      pcal_lv[k] = DEFAULTVALUE;
      pcal_lw[k] = DEFAULTVALUE;
      dc_chi2[k] = DEFAULTVALUE;
    }
  };
  
  ~DETDATA(){ };

  ClassDef(DETDATA,1);
};
//////////////////////////////////////////
#endif
