#ifndef __DATA_STRUCT__
#define __DATA_STRUCT__
#define MAXPART 50
#define MAXPART_MIX 6

#define DEFAULT_VALUE -111111
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
#include "TRandom3.h"
#include "TRegexp.h"
#include "PARTDATA.h"
#include "DETDATA.h"

//class DETDATA;
//class PARTDATA;

typedef struct DATAstr{
  Int_t npart;
  Float_t Q2, W, Nu, Xb, vxec, vyec, vzec, vxe, vye, vze, Pex, Pey, Pez, event, Pe, Ee, e_Ein, e_Eout, e_Epcal, e_npheltcc, e_nphehtcc, helic, e_chi2pid, e_pcal_lu, e_pcal_lv, e_pcal_lw, e_ecin_lu, e_ecin_lv, e_ecin_lw, e_ecout_lu, e_ecout_lv, e_ecout_lw, e_pcal_hx, e_pcal_hy, e_pcal_hz, e_ecin_hx, e_ecin_hy, e_ecin_hz, e_ecout_hx, e_ecout_hy, e_ecout_hz, e_trajx_sl0, e_trajx_sl1, e_trajx_sl2, e_trajx_sl3, e_trajx_sl4, e_trajx_sl5, e_trajy_sl0, e_trajy_sl1, e_trajy_sl2, e_trajy_sl3, e_trajy_sl4, e_trajy_sl5, e_trajz_sl0, e_trajz_sl1, e_trajz_sl2, e_trajz_sl3, e_trajz_sl4, e_trajz_sl5, e_pathtof, e_timetof, e_sector_tof, e_Beta, STTime, RFTime, e_dcx_rot_0, e_dcy_rot_0, e_dcx_rot_1, e_dcy_rot_1, e_dcx_rot_2, e_dcy_rot_2, e_sector_ltcc, e_sector_htcc, e_sector_ecal, e_dc_chi2, e_ftof1ax, e_ftof1ay, e_ftof1az, e_pcalx, e_pcaly, e_pcalz, e_ecalx, e_ecaly, e_ecalz, e_ltccx, e_ltccy, e_ltccz, e_htccx, e_htccy, e_htccz, e_ftof1bx, e_ftof1by, e_ftof1bz, e_ftof2x, e_ftof2y, e_ftof2z, helonline_hel, helonline_helRaw, helflip_hel, helflip_helRaw, helflip_event, e_dc_status, e_dc_ndf, e_sector_dc, e_statPart, e_DCPx, e_DCPy, e_DCPz, revent, y, th_e, phi_e, helicRaw;

 
  Float_t  ThetaPQ[MAXPART], PhiPQ[MAXPART], Zh[MAXPART], Pt2[MAXPART], Mx2[MAXPART], Xf[MAXPART], T[MAXPART], P[MAXPART], deltaZ[MAXPART], E[MAXPART], Px[MAXPART], Py[MAXPART], Pz[MAXPART], Ein[MAXPART], Eout[MAXPART], pid[MAXPART], Beta[MAXPART], vxh[MAXPART], vyh[MAXPART], vzh[MAXPART], npheltcc[MAXPART], nphehtcc[MAXPART],  chi2pid[MAXPART], Epcal[MAXPART], sector_ltcc[MAXPART], sector_htcc[MAXPART], sector_ecal[MAXPART], pcal_lu[MAXPART], pcal_lv[MAXPART], pcal_lw[MAXPART], ecin_lu[MAXPART], ecin_lv[MAXPART], ecin_lw[MAXPART], ecout_lu[MAXPART], ecout_lv[MAXPART], ecout_lw[MAXPART], sector_dc[MAXPART], statPart[MAXPART],  DCPx[MAXPART], DCPy[MAXPART], DCPz[MAXPART], trajx_sl0[MAXPART], trajx_sl1[MAXPART], trajx_sl2[MAXPART], trajx_sl3[MAXPART], trajx_sl4[MAXPART], trajx_sl5[MAXPART], trajy_sl0[MAXPART], trajy_sl1[MAXPART], trajy_sl2[MAXPART], trajy_sl3[MAXPART], trajy_sl4[MAXPART], trajy_sl5[MAXPART], trajz_sl0[MAXPART], trajz_sl1[MAXPART], trajz_sl2[MAXPART], trajz_sl3[MAXPART], trajz_sl4[MAXPART], trajz_sl5[MAXPART], pathtof[MAXPART], timetof[MAXPART], sector_tof[MAXPART], dcx_rot_0[MAXPART], dcy_rot_0[MAXPART], dcx_rot_1[MAXPART], dcy_rot_1[MAXPART], dcx_rot_2[MAXPART], dcy_rot_2[MAXPART], dc_chi2[MAXPART], ftof1ax[MAXPART], ftof1ay[MAXPART], ftof1az[MAXPART], pcalx[MAXPART], pcaly[MAXPART], pcalz[MAXPART], ecalx[MAXPART], ecaly[MAXPART], ecalz[MAXPART], ltccx[MAXPART], ltccy[MAXPART], ltccz[MAXPART], htccx[MAXPART], htccy[MAXPART], htccz[MAXPART], ftof1bx[MAXPART], ftof1by[MAXPART], ftof1bz[MAXPART], ftof2x[MAXPART], ftof2y[MAXPART], ftof2z[MAXPART], dc_status[MAXPART], dc_ndf[MAXPART];

  //// MC data ///
  Int_t mc_npart;
  Float_t mc_Q2, mc_W, mc_Nu, mc_Xb, mc_vxe, mc_vye, mc_vze, mc_Pex, mc_Pey, mc_Pez, mc_event, e_mcmass, mc_Pe, mc_Ee, mc_revent, mc_y, mc_th_e, mc_phi_e, mc_e_Beta;
 
  Float_t  mc_ThetaPQ[MAXPART], mc_PhiPQ[MAXPART], mc_Zh[MAXPART], mc_Pt2[MAXPART], mc_Mx2[MAXPART], mc_Xf[MAXPART], mc_T[MAXPART], mc_P[MAXPART], mc_deltaZ[MAXPART], mc_E[MAXPART], mc_Px[MAXPART], mc_Py[MAXPART], mc_Pz[MAXPART], mc_pid[MAXPART], mc_Beta[MAXPART], mc_vxh[MAXPART], mc_vyh[MAXPART], mc_vzh[MAXPART], mcmass[MAXPART];

} DATA;

typedef struct DATAMIXstr{
  Int_t npart;
  Float_t Q2, W, Nu, Xb, vxec, vyec, vzec, vxe, vye, vze, Pex, Pey, Pez, event, Pe, Ee, e_Ein, e_Eout, e_Epcal, e_npheltcc, e_nphehtcc, helic, e_chi2pid, e_pcal_lu, e_pcal_lv, e_pcal_lw, e_ecin_lu, e_ecin_lv, e_ecin_lw, e_ecout_lu, e_ecout_lv, e_ecout_lw, e_pcal_hx, e_pcal_hy, e_pcal_hz, e_ecin_hx, e_ecin_hy, e_ecin_hz, e_ecout_hx, e_ecout_hy, e_ecout_hz, e_trajx_sl0, e_trajx_sl1, e_trajx_sl2, e_trajx_sl3, e_trajx_sl4, e_trajx_sl5, e_trajy_sl0, e_trajy_sl1, e_trajy_sl2, e_trajy_sl3, e_trajy_sl4, e_trajy_sl5, e_trajz_sl0, e_trajz_sl1, e_trajz_sl2, e_trajz_sl3, e_trajz_sl4, e_trajz_sl5, e_pathtof, e_timetof, e_sector_tof, e_Beta, STTime, RFTime, e_dcx_rot_0, e_dcy_rot_0, e_dcx_rot_1, e_dcy_rot_1, e_dcx_rot_2, e_dcy_rot_2, e_sector_ltcc, e_sector_htcc, e_sector_ecal, e_dc_chi2, e_ftof1ax, e_ftof1ay, e_ftof1az, e_pcalx, e_pcaly, e_pcalz, e_ecalx, e_ecaly, e_ecalz, e_ltccx, e_ltccy, e_ltccz, e_htccx, e_htccy, e_htccz, e_ftof1bx, e_ftof1by, e_ftof1bz, e_ftof2x, e_ftof2y, e_ftof2z, helonline_hel, helonline_helRaw, helflip_hel, helflip_helRaw, helflip_event, e_dc_status, e_dc_ndf, e_sector_dc, e_statPart, e_DCPx, e_DCPy, e_DCPz, revent, y, th_e, phi_e, helicRaw, epsilon, gamm, fA, fB, fC, fV, fW;


  Int_t mix_npart;
  Int_t mc_mix_npart;
  /*  
  //// mix data///

  Int_t mix_npart;
  ////////////// detector data sucture ///////////
  // If you add a variable here, add it also to the Particle Class.
  Float_t beta[MAXPART][MAXPART_MIX], m2b[MAXPART][MAXPART_MIX], vx[MAXPART][MAXPART_MIX], vy[MAXPART][MAXPART_MIX], vz[MAXPART][MAXPART_MIX], dcx[MAXPART][MAXPART_MIX], dcy[MAXPART][MAXPART_MIX], dcz[MAXPART][MAXPART_MIX], statPart[MAXPART][MAXPART_MIX], dc_chi2[MAXPART][MAXPART_MIX], dc_ndf[MAXPART][MAXPART_MIX], pcal_lu[MAXPART][MAXPART_MIX], pcal_lv[MAXPART][MAXPART_MIX], pcal_lw[MAXPART][MAXPART_MIX];
  //// end detData ////
  //// particle data ////
  Float_t e[MAXPART][MAXPART_MIX], px[MAXPART][MAXPART_MIX], py[MAXPART][MAXPART_MIX], pz[MAXPART][MAXPART_MIX], helic002_phiH[MAXPART][MAXPART_MIX], helic005_phiH[MAXPART][MAXPART_MIX], helic010_phiH[MAXPART][MAXPART_MIX], helic020_phiH[MAXPART][MAXPART_MIX],phiHs[MAXPART][MAXPART_MIX];
  //// mc ///
  Int_t mc_mix_npart;
  Float_t mc_e[MAXPART][MAXPART_MIX], mc_px[MAXPART][MAXPART_MIX], mc_py[MAXPART][MAXPART_MIX], mc_pz[MAXPART][MAXPART_MIX], mc_helic002_phiH[MAXPART][MAXPART_MIX], mc_helic005_phiH[MAXPART][MAXPART_MIX], mc_helic010_phiH[MAXPART][MAXPART_MIX], mc_helic020_phiH[MAXPART][MAXPART_MIX], mc_phiHs[MAXPART][MAXPART_MIX];

  //// end mix data  ///
  */
  
  Float_t M[MAXPART], Phx[MAXPART], Phy[MAXPART], Phz[MAXPART], Z[MAXPART], Cospq[MAXPART], Pt2[MAXPART], Event[MAXPART], M2_01[MAXPART], M2_02[MAXPART], phiH[MAXPART], phiR[MAXPART], Mx2[MAXPART], xF[MAXPART], xF0[MAXPART], xF1[MAXPART], plcm[MAXPART], plcm0[MAXPART], plcm1[MAXPART], Eh[MAXPART], xFm[MAXPART], xFm0[MAXPART], xFm1[MAXPART], theta0[MAXPART], theta1[MAXPART], cos_theta_P0cm[MAXPART], sin_theta_P0cm[MAXPART], xFo[MAXPART], xFo0[MAXPART], xFo1[MAXPART], phiH_phiR[MAXPART], phiR_cov[MAXPART], p0T2[MAXPART], p1T2[MAXPART], phipq[MAXPART], phT2[MAXPART], etaCM0[MAXPART], etaCM1[MAXPART], etaBF0p[MAXPART], etaBF1p[MAXPART], etaBF0m[MAXPART], etaBF1m[MAXPART], etaBF0[MAXPART], etaBF1[MAXPART], phiR_ha[MAXPART], plcm0_r[MAXPART], plcm1_r[MAXPART], phiR_covH[MAXPART], E0_phcm[MAXPART], E1_phcm[MAXPART], helic002_phiR[MAXPART], helic005_phiR[MAXPART], helic010_phiR[MAXPART], helic020_phiR[MAXPART], KF[MAXPART], R[MAXPART], helic002_phiRst[MAXPART], helic005_phiRst[MAXPART], helic010_phiRst[MAXPART], helic020_phiRst[MAXPART], wUxS_phiR[MAXPART];
  
  //// MC data ///
  Int_t mc_npart;
  Float_t mc_Q2, mc_W, mc_Nu, mc_Xb, mc_vxe, mc_vye, mc_vze, mc_Pex, mc_Pey, mc_Pez, mc_event, e_mcmass, mc_Pe, mc_Ee, mc_revent, mc_y, mc_th_e, mc_phi_e, mc_e_Beta, mc_helic, mc_epsilon, mc_gamm, mc_fA, mc_fB, mc_fC, mc_fV, mc_fW;

  Float_t mc_M[MAXPART], mc_Phx[MAXPART], mc_Phy[MAXPART], mc_Phz[MAXPART], mc_Z[MAXPART], mc_Cospq[MAXPART], mc_Pt2[MAXPART], mc_Event[MAXPART], mc_M2_01[MAXPART], mc_M2_02[MAXPART], mc_phiH[MAXPART], mc_phiR[MAXPART], mc_Mx2[MAXPART], mc_xF[MAXPART], mc_xF0[MAXPART], mc_xF1[MAXPART], mc_plcm[MAXPART], mc_plcm0[MAXPART], mc_plcm1[MAXPART], mc_Eh[MAXPART], mc_xFm[MAXPART], mc_xFm0[MAXPART], mc_xFm1[MAXPART], mc_theta0[MAXPART], mc_theta1[MAXPART], mc_cos_theta_P0cm[MAXPART], mc_sin_theta_P0cm[MAXPART], mc_xFo[MAXPART], mc_xFo0[MAXPART], mc_xFo1[MAXPART], mc_phiH_phiR[MAXPART], mc_phiR_cov[MAXPART], mc_p0T2[MAXPART], mc_p1T2[MAXPART], mc_phipq[MAXPART], mc_phT2[MAXPART], mc_etaCM0[MAXPART], mc_etaCM1[MAXPART], mc_etaBF0p[MAXPART], mc_etaBF1p[MAXPART], mc_etaBF0m[MAXPART], mc_etaBF1m[MAXPART], mc_etaBF0[MAXPART], mc_etaBF1[MAXPART], mc_phiR_ha[MAXPART], mc_plcm0_r[MAXPART], mc_plcm1_r[MAXPART], mc_phiR_covH[MAXPART], mc_E0_phcm[MAXPART], mc_E1_phcm[MAXPART], mc_helic002_phiR[MAXPART], mc_helic005_phiR[MAXPART], mc_helic010_phiR[MAXPART], mc_helic020_phiR[MAXPART], mc_helic002_phiRst[MAXPART], mc_helic005_phiRst[MAXPART], mc_helic010_phiRst[MAXPART], mc_helic020_phiRst[MAXPART], mc_wUxS_phiR[MAXPART], mc_KF[MAXPART], mc_R[MAXPART];

} DATAMIX;

//// init mix tree ///
int initMixTree(TTree *t, DATAMIX *evnt = 0, TClonesArray *det = 0, TClonesArray *pdata = 0, TClonesArray *mc_pdata = 0){
  t->Branch("npart",&evnt->npart,"npart/I");
  t->Branch("M",evnt->M,"M[npart]/F");
  t->Branch("Phx",evnt->Phx,"Phx[npart]/F");
  t->Branch("Phy",evnt->Phy,"Phy[npart]/F");
  t->Branch("Phz",evnt->Phz,"Phz[npart]/F");
  t->Branch("Z",evnt->Z,"Z[npart]/F");
  t->Branch("Cospq",evnt->Cospq,"Cospq[npart]/F");
  t->Branch("Pt2",evnt->Pt2,"Pt2[npart]/F");
  t->Branch("Event",evnt->Event,"Event[npart]/F");
  t->Branch("M2_01",evnt->M2_01,"M2_01[npart]/F");
  t->Branch("M2_02",evnt->M2_02,"M2_02[npart]/F");
  t->Branch("phiH",evnt->phiH,"phiH[npart]/F");
  t->Branch("phiR",evnt->phiR,"phiR[npart]/F");
  t->Branch("helic002_phiR",evnt->helic002_phiR,"helic002_phiR[npart]/F");
  t->Branch("helic005_phiR",evnt->helic005_phiR,"helic005_phiR[npart]/F");
  t->Branch("helic010_phiR",evnt->helic010_phiR,"helic010_phiR[npart]/F");
  t->Branch("helic020_phiR",evnt->helic020_phiR,"helic020_phiR[npart]/F");
  t->Branch("helic002_phiRst",evnt->helic002_phiRst,"helic002_phiRst[npart]/F");
  t->Branch("helic005_phiRst",evnt->helic005_phiRst,"helic005_phiRst[npart]/F");
  t->Branch("helic010_phiRst",evnt->helic010_phiRst,"helic010_phiRst[npart]/F");
  t->Branch("helic020_phiRst",evnt->helic020_phiRst,"helic020_phiRst[npart]/F");

  t->Branch("wUxS_phiR",evnt->wUxS_phiR,"wUxS_phiR[npart]/F");
  t->Branch("R",evnt->R,"R[npart]/F");
  t->Branch("KF",evnt->KF,"KF[npart]/F");
  t->Branch("Mx2",evnt->Mx2,"Mx2[npart]/F");
  t->Branch("xF",evnt->xF,"xF[npart]/F");
  t->Branch("xF0",evnt->xF0,"xF0[npart]/F");
  t->Branch("xF1",evnt->xF1,"xF1[npart]/F");
  t->Branch("plcm",evnt->plcm,"plcm[npart]/F");
  t->Branch("plcm0",evnt->plcm0,"plcm0[npart]/F");
  t->Branch("plcm1",evnt->plcm1,"plcm1[npart]/F");
  t->Branch("Eh",evnt->Eh,"Eh[npart]/F");
  t->Branch("xFm",evnt->xFm,"xFm[npart]/F");
  t->Branch("xFm0",evnt->xFm0,"xFm0[npart]/F");
  t->Branch("xFm1",evnt->xFm1,"xFm1[npart]/F");
  t->Branch("theta0",evnt->theta0,"theta0[npart]/F");
  t->Branch("theta1",evnt->theta1,"theta1[npart]/F");
  t->Branch("cos_theta_P0cm",evnt->cos_theta_P0cm,"cos_theta_P0cm[npart]/F");
  t->Branch("sin_theta_P0cm",evnt->sin_theta_P0cm,"sin_theta_P0cm[npart]/F");
  t->Branch("xFo",evnt->xFo,"xFo[npart]/F");
  t->Branch("xFo0",evnt->xFo0,"xFo0[npart]/F");
  t->Branch("xFo1",evnt->xFo1,"xFo1[npart]/F");
  t->Branch("phiH_phiR",evnt->phiH_phiR,"phiH_phiR[npart]/F");
  t->Branch("phiR_cov",evnt->phiR_cov,"phiR_cov[npart]/F");
  t->Branch("p0T2",evnt->p0T2,"p0T2[npart]/F");
  t->Branch("p1T2",evnt->p1T2,"p1T2[npart]/F");
  t->Branch("phipq",evnt->phipq,"phipq[npart]/F");
  t->Branch("phT2",evnt->phT2,"phT2[npart]/F");
  t->Branch("etaCM0",evnt->etaCM0,"etaCM0[npart]/F");
  t->Branch("etaCM1",evnt->etaCM1,"etaCM1[npart]/F");
  t->Branch("etaBF0p",evnt->etaBF0p,"etaBF0p[npart]/F");
  t->Branch("etaBF1p",evnt->etaBF1p,"etaBF1p[npart]/F");
  t->Branch("etaBF0m",evnt->etaBF0m,"etaBF0m[npart]/F");
  t->Branch("etaBF1m",evnt->etaBF1m,"etaBF1m[npart]/F");
  t->Branch("etaBF0",evnt->etaBF0,"etaBF0[npart]/F");
  t->Branch("etaBF1",evnt->etaBF1,"etaBF1[npart]/F");
  t->Branch("phiR_ha",evnt->phiR_ha,"phiR_ha[npart]/F");
  t->Branch("plcm0_r",evnt->plcm0_r,"plcm0_r[npart]/F");
  t->Branch("plcm1_r",evnt->plcm1_r,"plcm1_r[npart]/F");
  t->Branch("phiR_covH",evnt->phiR_covH,"phiR_covH[npart]/F");
  t->Branch("E0_phcm",evnt->E0_phcm,"E0_phcm[npart]/F");
  t->Branch("E1_phcm",evnt->E1_phcm,"E1_phcm[npart]/F");
  //// end FS rec///

  t->Branch("mix_npart",&evnt->mix_npart,"mix_npart/I");
  t->Branch("mc_mix_npart",&evnt->mc_mix_npart,"mc_mix_npart/I");
  
  t->Branch("det",&det,MAXPART);
  t->Branch("pdata",&pdata,MAXPART);
  t->Branch("mc_pdata",&mc_pdata,MAXPART);

  /*
  //// mix data ////
  //// det data ///
  t->Branch("beta",evnt->beta,"beta[npart][5]/F");
  t->Branch("m2b",evnt->m2b,"m2b[npart][5]/F");
  t->Branch("vx",evnt->vx,"vx[npart][5]/F");
  t->Branch("vy",evnt->vy,"vy[npart][5]/F");
  t->Branch("vz",evnt->vz,"vz[npart][5]/F");
  t->Branch("dcx",evnt->dcx,"dcx[npart][5]/F");
  t->Branch("dcy",evnt->dcy,"dcy[npart][5]/F");
  t->Branch("dcz",evnt->dcz,"dcz[npart][5]/F");
  t->Branch("statPart",evnt->statPart,"statPart[npart][5]/F");
  t->Branch("dc_chi2",evnt->dc_chi2,"dc_chi2[npart][5]/F");
  t->Branch("dc_ndf",evnt->dc_ndf,"dc_ndf[npart][5]/F");
  t->Branch("pcal_lu",evnt->pcal_lu,"pcal_lu[npart][5]/F");
  t->Branch("pcal_lv",evnt->pcal_lv,"pcal_lv[npart][5]/F");
  t->Branch("pcal_lw",evnt->pcal_lw,"pcal_lw[npart][5]/F");
  t->Branch("helic002_phiH",evnt->helic002_phiH,"helic002_phiH[npart][5]/F");
  t->Branch("helic005_phiH",evnt->helic005_phiH,"helic005_phiH[npart][5]/F");
  t->Branch("helic010_phiH",evnt->helic010_phiH,"helic010_phiH[npart][5]/F");
  t->Branch("helic020_phiH",evnt->helic020_phiH,"helic020_phiH[npart][5]/F");
  t->Branch("phiHs",evnt->phiHs,"phiHs[npart][5]/F");

  //// end det data///
  //// set particles ///
  t->Branch("e",evnt->e,"e[npart][5]/F");
  t->Branch("px",evnt->px,"px[npart][5]/F");
  t->Branch("py",evnt->py,"py[npart][5]/F");
  t->Branch("pz",evnt->pz,"pz[npart][5]/F");
  //// mc ///
  t->Branch("mc_mix_npart",&evnt->mc_mix_npart,"mc_mix_npart/I");
  t->Branch("mc_e",evnt->mc_e,"mc_e[mc_npart][5]/F");
  t->Branch("mc_px",evnt->mc_px,"mc_px[mc_npart][5]/F");
  t->Branch("mc_py",evnt->mc_py,"mc_py[mc_npart][5]/F");
  t->Branch("mc_pz",evnt->mc_pz,"mc_pz[mc_npart][5]/F");
  t->Branch("mc_helic002_phiH",evnt->mc_helic002_phiH,"mc_helic002_phiH[mc_npart][5]/F");
  t->Branch("mc_helic005_phiH",evnt->mc_helic005_phiH,"mc_helic005_phiH[mc_npart][5]/F");
  t->Branch("mc_helic010_phiH",evnt->mc_helic010_phiH,"mc_helic010_phiH[mc_npart][5]/F");
  t->Branch("mc_helic020_phiH",evnt->mc_helic020_phiH,"mc_helic020_phiH[mc_npart][5]/F");
  t->Branch("mc_phiHs",evnt->mc_phiHs,"mc_phiHs[mc_npart][5]/F");
  //// end particle ///
  //// mix data end ////
  */

  
  //// FS MC ///
  t->Branch("mc_npart",&evnt->mc_npart,"mc_npart/I");
  t->Branch("mc_M",evnt->mc_M,"mc_M[mc_npart]/F");
  t->Branch("mc_Phx",evnt->mc_Phx,"mc_Phx[mc_npart]/F");
  t->Branch("mc_Phy",evnt->mc_Phy,"mc_Phy[mc_npart]/F");
  t->Branch("mc_Phz",evnt->mc_Phz,"mc_Phz[mc_npart]/F");
  t->Branch("mc_Z",evnt->mc_Z,"mc_Z[mc_npart]/F");
  t->Branch("mc_Cospq",evnt->mc_Cospq,"mc_Cospq[mc_npart]/F");
  t->Branch("mc_Pt2",evnt->mc_Pt2,"mc_Pt2[mc_npart]/F");
  t->Branch("mc_Event",evnt->mc_Event,"mc_Event[mc_npart]/F");
  t->Branch("mc_M2_01",evnt->mc_M2_01,"mc_M2_01[mc_npart]/F");
  t->Branch("mc_M2_02",evnt->mc_M2_02,"mc_M2_02[mc_npart]/F");
  t->Branch("mc_phiH",evnt->mc_phiH,"mc_phiH[mc_npart]/F");
  t->Branch("mc_phiR",evnt->mc_phiR,"mc_phiR[mc_npart]/F");
  t->Branch("mc_helic002_phiR",evnt->mc_helic002_phiR,"mc_helic002_phiR[mc_npart]/F");
  t->Branch("mc_helic005_phiR",evnt->mc_helic005_phiR,"mc_helic005_phiR[mc_npart]/F");
  t->Branch("mc_helic010_phiR",evnt->mc_helic010_phiR,"mc_helic010_phiR[mc_npart]/F");
  t->Branch("mc_helic020_phiR",evnt->mc_helic020_phiR,"mc_helic020_phiR[mc_npart]/F");
  t->Branch("mc_helic002_phiRst",evnt->mc_helic002_phiRst,"mc_helic002_phiRst[mc_npart]/F");
  t->Branch("mc_helic005_phiRst",evnt->mc_helic005_phiRst,"mc_helic005_phiRst[mc_npart]/F");
  t->Branch("mc_helic010_phiRst",evnt->mc_helic010_phiRst,"mc_helic010_phiRst[mc_npart]/F");
  t->Branch("mc_helic020_phiRst",evnt->mc_helic020_phiRst,"mc_helic020_phiRst[mc_npart]/F");
  t->Branch("mc_wUxS_phiR",evnt->mc_wUxS_phiR,"mc_wUxS_phiR[mc_npart]/F");
  t->Branch("mc_Mx2",evnt->mc_Mx2,"mc_Mx2[mc_npart]/F");
  t->Branch("mc_xF",evnt->mc_xF,"mc_xF[mc_npart]/F");
  t->Branch("mc_xF0",evnt->mc_xF0,"mc_xF0[mc_npart]/F");
  t->Branch("mc_xF1",evnt->mc_xF1,"mc_xF1[mc_npart]/F");
  t->Branch("mc_plcm",evnt->mc_plcm,"mc_plcm[mc_npart]/F");
  t->Branch("mc_plcm0",evnt->mc_plcm0,"mc_plcm0[mc_npart]/F");
  t->Branch("mc_plcm1",evnt->mc_plcm1,"mc_plcm1[mc_npart]/F");
  t->Branch("mc_Eh",evnt->mc_Eh,"mc_Eh[mc_npart]/F");
  t->Branch("mc_xFm",evnt->mc_xFm,"mc_xFm[mc_npart]/F");
  t->Branch("mc_xFm0",evnt->mc_xFm0,"mc_xFm0[mc_npart]/F");
  t->Branch("mc_xFm1",evnt->mc_xFm1,"mc_xFm1[mc_npart]/F");
  t->Branch("mc_theta0",evnt->mc_theta0,"mc_theta0[mc_npart]/F");
  t->Branch("mc_theta1",evnt->mc_theta1,"mc_theta1[mc_npart]/F");
  t->Branch("mc_cos_theta_P0cm",evnt->mc_cos_theta_P0cm,"mc_cos_theta_P0cm[mc_npart]/F");
  t->Branch("mc_sin_theta_P0cm",evnt->mc_sin_theta_P0cm,"mc_sin_theta_P0cm[mc_npart]/F");
  t->Branch("mc_xFo",evnt->mc_xFo,"mc_xFo[mc_npart]/F");
  t->Branch("mc_xFo0",evnt->mc_xFo0,"mc_xFo0[mc_npart]/F");
  t->Branch("mc_xFo1",evnt->mc_xFo1,"mc_xFo1[mc_npart]/F");
  t->Branch("mc_phiH_phiR",evnt->mc_phiH_phiR,"mc_phiH_phiR[mc_npart]/F");
  t->Branch("mc_phiR_cov",evnt->mc_phiR_cov,"mc_phiR_cov[mc_npart]/F");
  t->Branch("mc_p0T2",evnt->mc_p0T2,"mc_p0T2[mc_npart]/F");
  t->Branch("mc_p1T2",evnt->mc_p1T2,"mc_p1T2[mc_npart]/F");
  t->Branch("mc_phipq",evnt->mc_phipq,"mc_phipq[mc_npart]/F");
  t->Branch("mc_phT2",evnt->mc_phT2,"mc_phT2[mc_npart]/F");
  t->Branch("mc_etaCM0",evnt->mc_etaCM0,"mc_etaCM0[mc_npart]/F");
  t->Branch("mc_etaCM1",evnt->mc_etaCM1,"mc_etaCM1[mc_npart]/F");
  t->Branch("mc_etaBF0p",evnt->mc_etaBF0p,"mc_etaBF0p[mc_npart]/F");
  t->Branch("mc_etaBF1p",evnt->mc_etaBF1p,"mc_etaBF1p[mc_npart]/F");
  t->Branch("mc_etaBF0m",evnt->mc_etaBF0m,"mc_etaBF0m[mc_npart]/F");
  t->Branch("mc_etaBF1m",evnt->mc_etaBF1m,"mc_etaBF1m[mc_npart]/F");
  t->Branch("mc_etaBF0",evnt->mc_etaBF0,"mc_etaBF0[mc_npart]/F");
  t->Branch("mc_etaBF1",evnt->mc_etaBF1,"mc_etaBF1[mc_npart]/F");
  t->Branch("mc_phiR_ha",evnt->mc_phiR_ha,"mc_phiR_ha[mc_npart]/F");
  t->Branch("mc_plcm0_r",evnt->mc_plcm0_r,"mc_plcm0_r[mc_npart]/F");
  t->Branch("mc_plcm1_r",evnt->mc_plcm1_r,"mc_plcm1_r[mc_npart]/F");
  t->Branch("mc_phiR_covH",evnt->mc_phiR_covH,"mc_phiR_covH[mc_npart]/F");
  t->Branch("mc_E0_phcm",evnt->mc_E0_phcm,"mc_E0_phcm[mc_npart]/F");
  t->Branch("mc_E1_phcm",evnt->mc_E1_phcm,"mc_E1_phcm[mc_npart]/F");
  t->Branch("mc_R",evnt->mc_R,"mc_R[mc_npart]/F");
  t->Branch("mc_KF",evnt->mc_KF,"mc_KF[mc_npart]/F");

  //// end FS MC ///
  
  //// Electron variables ///
  t->Branch("Q2",&evnt->Q2,"Q2/F");
  t->Branch("W",&evnt->W,"W/F");
  t->Branch("Nu",&evnt->Nu,"Nu/F");
  t->Branch("Xb",&evnt->Xb,"Xb/F");
  t->Branch("vxec",&evnt->vxec,"vxec/F");
  t->Branch("vyec",&evnt->vyec,"vyec/F");
  t->Branch("vzec",&evnt->vzec,"vzec/F");
  t->Branch("vxe",&evnt->vxe,"vxe/F");
  t->Branch("vye",&evnt->vye,"vye/F");
  t->Branch("vze",&evnt->vze,"vze/F");
  t->Branch("Pex",&evnt->Pex,"Pex/F");
  t->Branch("Pey",&evnt->Pey,"Pey/F");
  t->Branch("Pez",&evnt->Pez,"Pez/F");
  t->Branch("event",&evnt->event,"event/F");
  t->Branch("Pe",&evnt->Pe,"Pe/F");
  t->Branch("Ee",&evnt->Ee,"Ee/F");
  t->Branch("e_Ein",&evnt->e_Ein,"e_Ein/F");
  t->Branch("e_Eout",&evnt->e_Eout,"e_Eout/F");
  t->Branch("e_Epcal",&evnt->e_Epcal,"e_Epcal/F");
  t->Branch("e_npheltcc",&evnt->e_npheltcc,"e_npheltcc/F");
  t->Branch("e_nphehtcc",&evnt->e_nphehtcc,"e_nphehtcc/F");
  t->Branch("helic",&evnt->helic,"helic/F");
  t->Branch("e_chi2pid",&evnt->e_chi2pid,"e_chi2pid/F");
  t->Branch("e_pcal_lu",&evnt->e_pcal_lu,"e_pcal_lu/F");
  t->Branch("e_pcal_lv",&evnt->e_pcal_lv,"e_pcal_lv/F");
  t->Branch("e_pcal_lw",&evnt->e_pcal_lw,"e_pcal_lw/F");
  t->Branch("e_ecin_lu",&evnt->e_ecin_lu,"e_ecin_lu/F");
  t->Branch("e_ecin_lv",&evnt->e_ecin_lv,"e_ecin_lv/F");
  t->Branch("e_ecin_lw",&evnt->e_ecin_lw,"e_ecin_lw/F");
  t->Branch("e_ecout_lu",&evnt->e_ecout_lu,"e_ecout_lu/F");
  t->Branch("e_ecout_lv",&evnt->e_ecout_lv,"e_ecout_lv/F");
  t->Branch("e_ecout_lw",&evnt->e_ecout_lw,"e_ecout_lw/F");
  t->Branch("e_pcal_hx",&evnt->e_pcal_hx,"e_pcal_hx/F");
  t->Branch("e_pcal_hy",&evnt->e_pcal_hy,"e_pcal_hy/F");
  t->Branch("e_pcal_hz",&evnt->e_pcal_hz,"e_pcal_hz/F");
  t->Branch("e_ecin_hx",&evnt->e_ecin_hx,"e_ecin_hx/F");
  t->Branch("e_ecin_hy",&evnt->e_ecin_hy,"e_ecin_hy/F");
  t->Branch("e_ecin_hz",&evnt->e_ecin_hz,"e_ecin_hz/F");
  t->Branch("e_ecout_hx",&evnt->e_ecout_hx,"e_ecout_hx/F");
  t->Branch("e_ecout_hy",&evnt->e_ecout_hy,"e_ecout_hy/F");
  t->Branch("e_ecout_hz",&evnt->e_ecout_hz,"e_ecout_hz/F");
  t->Branch("e_trajx_sl0",&evnt->e_trajx_sl0,"e_trajx_sl0/F");
  t->Branch("e_trajx_sl1",&evnt->e_trajx_sl1,"e_trajx_sl1/F");
  t->Branch("e_trajx_sl2",&evnt->e_trajx_sl2,"e_trajx_sl2/F");
  t->Branch("e_trajx_sl3",&evnt->e_trajx_sl3,"e_trajx_sl3/F");
  t->Branch("e_trajx_sl4",&evnt->e_trajx_sl4,"e_trajx_sl4/F");
  t->Branch("e_trajx_sl5",&evnt->e_trajx_sl5,"e_trajx_sl5/F");
  t->Branch("e_trajy_sl0",&evnt->e_trajy_sl0,"e_trajy_sl0/F");
  t->Branch("e_trajy_sl1",&evnt->e_trajy_sl1,"e_trajy_sl1/F");
  t->Branch("e_trajy_sl2",&evnt->e_trajy_sl2,"e_trajy_sl2/F");
  t->Branch("e_trajy_sl3",&evnt->e_trajy_sl3,"e_trajy_sl3/F");
  t->Branch("e_trajy_sl4",&evnt->e_trajy_sl4,"e_trajy_sl4/F");
  t->Branch("e_trajy_sl5",&evnt->e_trajy_sl5,"e_trajy_sl5/F");
  t->Branch("e_trajz_sl0",&evnt->e_trajz_sl0,"e_trajz_sl0/F");
  t->Branch("e_trajz_sl1",&evnt->e_trajz_sl1,"e_trajz_sl1/F");
  t->Branch("e_trajz_sl2",&evnt->e_trajz_sl2,"e_trajz_sl2/F");
  t->Branch("e_trajz_sl3",&evnt->e_trajz_sl3,"e_trajz_sl3/F");
  t->Branch("e_trajz_sl4",&evnt->e_trajz_sl4,"e_trajz_sl4/F");
  t->Branch("e_trajz_sl5",&evnt->e_trajz_sl5,"e_trajz_sl5/F");
  t->Branch("e_pathtof",&evnt->e_pathtof,"e_pathtof/F");
  t->Branch("e_timetof",&evnt->e_timetof,"e_timetof/F");
  t->Branch("e_sector_tof",&evnt->e_sector_tof,"e_sector_tof/F");
  t->Branch("e_Beta",&evnt->e_Beta,"e_Beta/F");
  t->Branch("STTime",&evnt->STTime,"STTime/F");
  t->Branch("RFTime",&evnt->RFTime,"RFTime/F");
  t->Branch("e_dcx_rot_0",&evnt->e_dcx_rot_0,"e_dcx_rot_0/F");
  t->Branch("e_dcy_rot_0",&evnt->e_dcy_rot_0,"e_dcy_rot_0/F");
  t->Branch("e_dcx_rot_1",&evnt->e_dcx_rot_1,"e_dcx_rot_1/F");
  t->Branch("e_dcy_rot_1",&evnt->e_dcy_rot_1,"e_dcy_rot_1/F");
  t->Branch("e_dcx_rot_2",&evnt->e_dcx_rot_2,"e_dcx_rot_2/F");
  t->Branch("e_dcy_rot_2",&evnt->e_dcy_rot_2,"e_dcy_rot_2/F");
  t->Branch("e_sector_ltcc",&evnt->e_sector_ltcc,"e_sector_ltcc/F");
  t->Branch("e_sector_htcc",&evnt->e_sector_htcc,"e_sector_htcc/F");
  t->Branch("e_sector_ecal",&evnt->e_sector_ecal,"e_sector_ecal/F");
  t->Branch("revent",&evnt->revent,"revent/F");
  t->Branch("e_dc_chi2",&evnt->e_dc_chi2,"e_dc_chi2/F");
  t->Branch("e_ftof1ax",&evnt->e_ftof1ax,"e_ftof1ax/F");
  t->Branch("e_ftof1ay",&evnt->e_ftof1ay,"e_ftof1ay/F");
  t->Branch("e_ftof1az",&evnt->e_ftof1az,"e_ftof1az/F");
  t->Branch("e_pcalx",&evnt->e_pcalx,"e_pcalx/F");
  t->Branch("e_pcaly",&evnt->e_pcaly,"e_pcaly/F");
  t->Branch("e_pcalz",&evnt->e_pcalz,"e_pcalz/F");
  t->Branch("e_ecalx",&evnt->e_ecalx,"e_ecalx/F");
  t->Branch("e_ecaly",&evnt->e_ecaly,"e_ecaly/F");
  t->Branch("e_ecalz",&evnt->e_ecalz,"e_ecalz/F");
  t->Branch("e_ltccx",&evnt->e_ltccx,"e_ltccx/F");
  t->Branch("e_ltccy",&evnt->e_ltccy,"e_ltccy/F");
  t->Branch("e_ltccz",&evnt->e_ltccz,"e_ltccz/F");
  t->Branch("e_htccx",&evnt->e_htccx,"e_htccx/F");
  t->Branch("e_htccy",&evnt->e_htccy,"e_htccy/F");
  t->Branch("e_htccz",&evnt->e_htccz,"e_htccz/F");
  t->Branch("e_ftof1bx",&evnt->e_ftof1bx,"e_ftof1bx/F");
  t->Branch("e_ftof1by",&evnt->e_ftof1by,"e_ftof1by/F");
  t->Branch("e_ftof1bz",&evnt->e_ftof1bz,"e_ftof1bz/F");
  t->Branch("e_ftof2x",&evnt->e_ftof2x,"e_ftof2x/F");
  t->Branch("e_ftof2y",&evnt->e_ftof2y,"e_ftof2y/F");
  t->Branch("e_ftof2z",&evnt->e_ftof2z,"e_ftof2z/F");
  t->Branch("helonline_hel",&evnt->helonline_hel,"helonline_hel/F");
  t->Branch("helonline_helRaw",&evnt->helonline_helRaw,"helonline_helRaw/F");
  t->Branch("helflip_hel",&evnt->helflip_hel,"helflip_hel/F");
  t->Branch("helflip_helRaw",&evnt->helflip_helRaw,"helflip_helRaw/F");
  t->Branch("helflip_event",&evnt->helflip_event,"helflip_event/F");
  t->Branch("e_dc_status",&evnt->e_dc_status,"e_dc_status/F");
  t->Branch("e_dc_ndf",&evnt->e_dc_ndf,"e_dc_ndf/F");
  t->Branch("e_sector_dc",&evnt->e_sector_dc,"e_sector_dc/F");
  t->Branch("e_statPart",&evnt->e_statPart,"e_statPart/F");
  t->Branch("e_DCPx",&evnt->e_DCPx,"e_DCPx/F");
  t->Branch("e_DCPy",&evnt->e_DCPy,"e_DCPy/F");
  t->Branch("e_DCPz",&evnt->e_DCPz,"e_DCPz/F");
  t->Branch("y",&evnt->y,"y/F");
  t->Branch("th_e",&evnt->th_e,"th_e/F");
  t->Branch("phi_e",&evnt->phi_e,"phi_e/F");
  t->Branch("helicRaw",&evnt->helicRaw,"helicRaw/F");
  t->Branch("fA",&evnt->fA,"fA/F");
  t->Branch("fB",&evnt->fB,"fB/F");
  t->Branch("fC",&evnt->fC,"fC/F");
  t->Branch("fV",&evnt->fV,"fV/F");
  t->Branch("fW",&evnt->fW,"fW/F");

  ////  End electron variables ///
  //// Electrons MC///
  t->Branch("mc_Q2",&evnt->mc_Q2,"mc_Q2/F");
  t->Branch("mc_W",&evnt->mc_W,"mc_W/F");
  t->Branch("mc_Nu",&evnt->mc_Nu,"mc_Nu/F");
  t->Branch("mc_Xb",&evnt->mc_Xb,"mc_Xb/F");
  t->Branch("mc_vxe",&evnt->mc_vxe,"mc_vxe/F");
  t->Branch("mc_vye",&evnt->mc_vye,"mc_vye/F");
  t->Branch("mc_vze",&evnt->mc_vze,"mc_vze/F");
  t->Branch("mc_Pex",&evnt->mc_Pex,"mc_Pex/F");
  t->Branch("mc_Pey",&evnt->mc_Pey,"mc_Pey/F");
  t->Branch("mc_Pez",&evnt->mc_Pez,"mc_Pez/F");
  t->Branch("mc_event",&evnt->mc_event,"mc_event/F");
  t->Branch("e_mcmass",&evnt->e_mcmass,"e_mcmass/F");
  t->Branch("mc_Pe",&evnt->mc_Pe,"mc_Pe/F");
  t->Branch("mc_Ee",&evnt->mc_Ee,"mc_Ee/F");
  t->Branch("mc_revent",&evnt->mc_revent,"mc_revent/F");
  t->Branch("mc_y",&evnt->mc_y,"mc_y/F");
  t->Branch("mc_th_e",&evnt->mc_th_e,"mc_th_e/F");
  t->Branch("mc_phi_e",&evnt->mc_phi_e,"mc_phi_e/F");
  t->Branch("mc_e_Beta",&evnt->mc_e_Beta,"mc_e_Beta/F");
  t->Branch("mc_helic",&evnt->mc_helic,"mc_helic/F");
  t->Branch("mc_fA",&evnt->mc_fA,"mc_fA/F");
  t->Branch("mc_fB",&evnt->mc_fB,"mc_fB/F");
  t->Branch("mc_fC",&evnt->mc_fC,"mc_fC/F");
  t->Branch("mc_fV",&evnt->mc_fV,"mc_fV/F");
  t->Branch("mc_fW",&evnt->mc_fW,"mc_fW/F");

  ////  END electron variables MC ///

  return 0; 
}
//// end init mix tree///

int initTree(TTree *t, DATA* evnt = 0){
  //// Electron variables ///
  t->Branch("npart",&evnt->npart,"npart/I");
  t->Branch("Q2",&evnt->Q2,"Q2/F");
  t->Branch("W",&evnt->W,"W/F");
  t->Branch("Nu",&evnt->Nu,"Nu/F");
  t->Branch("Xb",&evnt->Xb,"Xb/F");
  t->Branch("vxec",&evnt->vxec,"vxec/F");
  t->Branch("vyec",&evnt->vyec,"vyec/F");
  t->Branch("vzec",&evnt->vzec,"vzec/F");
  t->Branch("vxe",&evnt->vxe,"vxe/F");
  t->Branch("vye",&evnt->vye,"vye/F");
  t->Branch("vze",&evnt->vze,"vze/F");
  t->Branch("Pex",&evnt->Pex,"Pex/F");
  t->Branch("Pey",&evnt->Pey,"Pey/F");
  t->Branch("Pez",&evnt->Pez,"Pez/F");
  t->Branch("event",&evnt->event,"event/F");
  t->Branch("Pe",&evnt->Pe,"Pe/F");
  t->Branch("Ee",&evnt->Ee,"Ee/F");
  t->Branch("e_Ein",&evnt->e_Ein,"e_Ein/F");
  t->Branch("e_Eout",&evnt->e_Eout,"e_Eout/F");
  t->Branch("e_Epcal",&evnt->e_Epcal,"e_Epcal/F");
  t->Branch("e_npheltcc",&evnt->e_npheltcc,"e_npheltcc/F");
  t->Branch("e_nphehtcc",&evnt->e_nphehtcc,"e_nphehtcc/F");
  t->Branch("helic",&evnt->helic,"helic/F");
  t->Branch("e_chi2pid",&evnt->e_chi2pid,"e_chi2pid/F");
  t->Branch("e_pcal_lu",&evnt->e_pcal_lu,"e_pcal_lu/F");
  t->Branch("e_pcal_lv",&evnt->e_pcal_lv,"e_pcal_lv/F");
  t->Branch("e_pcal_lw",&evnt->e_pcal_lw,"e_pcal_lw/F");
  t->Branch("e_ecin_lu",&evnt->e_ecin_lu,"e_ecin_lu/F");
  t->Branch("e_ecin_lv",&evnt->e_ecin_lv,"e_ecin_lv/F");
  t->Branch("e_ecin_lw",&evnt->e_ecin_lw,"e_ecin_lw/F");
  t->Branch("e_ecout_lu",&evnt->e_ecout_lu,"e_ecout_lu/F");
  t->Branch("e_ecout_lv",&evnt->e_ecout_lv,"e_ecout_lv/F");
  t->Branch("e_ecout_lw",&evnt->e_ecout_lw,"e_ecout_lw/F");
  t->Branch("e_pcal_hx",&evnt->e_pcal_hx,"e_pcal_hx/F");
  t->Branch("e_pcal_hy",&evnt->e_pcal_hy,"e_pcal_hy/F");
  t->Branch("e_pcal_hz",&evnt->e_pcal_hz,"e_pcal_hz/F");
  t->Branch("e_ecin_hx",&evnt->e_ecin_hx,"e_ecin_hx/F");
  t->Branch("e_ecin_hy",&evnt->e_ecin_hy,"e_ecin_hy/F");
  t->Branch("e_ecin_hz",&evnt->e_ecin_hz,"e_ecin_hz/F");
  t->Branch("e_ecout_hx",&evnt->e_ecout_hx,"e_ecout_hx/F");
  t->Branch("e_ecout_hy",&evnt->e_ecout_hy,"e_ecout_hy/F");
  t->Branch("e_ecout_hz",&evnt->e_ecout_hz,"e_ecout_hz/F");
  t->Branch("e_trajx_sl0",&evnt->e_trajx_sl0,"e_trajx_sl0/F");
  t->Branch("e_trajx_sl1",&evnt->e_trajx_sl1,"e_trajx_sl1/F");
  t->Branch("e_trajx_sl2",&evnt->e_trajx_sl2,"e_trajx_sl2/F");
  t->Branch("e_trajx_sl3",&evnt->e_trajx_sl3,"e_trajx_sl3/F");
  t->Branch("e_trajx_sl4",&evnt->e_trajx_sl4,"e_trajx_sl4/F");
  t->Branch("e_trajx_sl5",&evnt->e_trajx_sl5,"e_trajx_sl5/F");
  t->Branch("e_trajy_sl0",&evnt->e_trajy_sl0,"e_trajy_sl0/F");
  t->Branch("e_trajy_sl1",&evnt->e_trajy_sl1,"e_trajy_sl1/F");
  t->Branch("e_trajy_sl2",&evnt->e_trajy_sl2,"e_trajy_sl2/F");
  t->Branch("e_trajy_sl3",&evnt->e_trajy_sl3,"e_trajy_sl3/F");
  t->Branch("e_trajy_sl4",&evnt->e_trajy_sl4,"e_trajy_sl4/F");
  t->Branch("e_trajy_sl5",&evnt->e_trajy_sl5,"e_trajy_sl5/F");
  t->Branch("e_trajz_sl0",&evnt->e_trajz_sl0,"e_trajz_sl0/F");
  t->Branch("e_trajz_sl1",&evnt->e_trajz_sl1,"e_trajz_sl1/F");
  t->Branch("e_trajz_sl2",&evnt->e_trajz_sl2,"e_trajz_sl2/F");
  t->Branch("e_trajz_sl3",&evnt->e_trajz_sl3,"e_trajz_sl3/F");
  t->Branch("e_trajz_sl4",&evnt->e_trajz_sl4,"e_trajz_sl4/F");
  t->Branch("e_trajz_sl5",&evnt->e_trajz_sl5,"e_trajz_sl5/F");
  t->Branch("e_pathtof",&evnt->e_pathtof,"e_pathtof/F");
  t->Branch("e_timetof",&evnt->e_timetof,"e_timetof/F");
  t->Branch("e_sector_tof",&evnt->e_sector_tof,"e_sector_tof/F");
  t->Branch("e_Beta",&evnt->e_Beta,"e_Beta/F");
  t->Branch("STTime",&evnt->STTime,"STTime/F");
  t->Branch("RFTime",&evnt->RFTime,"RFTime/F");
  t->Branch("e_dcx_rot_0",&evnt->e_dcx_rot_0,"e_dcx_rot_0/F");
  t->Branch("e_dcy_rot_0",&evnt->e_dcy_rot_0,"e_dcy_rot_0/F");
  t->Branch("e_dcx_rot_1",&evnt->e_dcx_rot_1,"e_dcx_rot_1/F");
  t->Branch("e_dcy_rot_1",&evnt->e_dcy_rot_1,"e_dcy_rot_1/F");
  t->Branch("e_dcx_rot_2",&evnt->e_dcx_rot_2,"e_dcx_rot_2/F");
  t->Branch("e_dcy_rot_2",&evnt->e_dcy_rot_2,"e_dcy_rot_2/F");
  t->Branch("e_sector_ltcc",&evnt->e_sector_ltcc,"e_sector_ltcc/F");
  t->Branch("e_sector_htcc",&evnt->e_sector_htcc,"e_sector_htcc/F");
  t->Branch("e_sector_ecal",&evnt->e_sector_ecal,"e_sector_ecal/F");
  t->Branch("revent",&evnt->revent,"revent/F");
  t->Branch("e_dc_chi2",&evnt->e_dc_chi2,"e_dc_chi2/F");
  t->Branch("e_ftof1ax",&evnt->e_ftof1ax,"e_ftof1ax/F");
  t->Branch("e_ftof1ay",&evnt->e_ftof1ay,"e_ftof1ay/F");
  t->Branch("e_ftof1az",&evnt->e_ftof1az,"e_ftof1az/F");
  t->Branch("e_pcalx",&evnt->e_pcalx,"e_pcalx/F");
  t->Branch("e_pcaly",&evnt->e_pcaly,"e_pcaly/F");
  t->Branch("e_pcalz",&evnt->e_pcalz,"e_pcalz/F");
  t->Branch("e_ecalx",&evnt->e_ecalx,"e_ecalx/F");
  t->Branch("e_ecaly",&evnt->e_ecaly,"e_ecaly/F");
  t->Branch("e_ecalz",&evnt->e_ecalz,"e_ecalz/F");
  t->Branch("e_ltccx",&evnt->e_ltccx,"e_ltccx/F");
  t->Branch("e_ltccy",&evnt->e_ltccy,"e_ltccy/F");
  t->Branch("e_ltccz",&evnt->e_ltccz,"e_ltccz/F");
  t->Branch("e_htccx",&evnt->e_htccx,"e_htccx/F");
  t->Branch("e_htccy",&evnt->e_htccy,"e_htccy/F");
  t->Branch("e_htccz",&evnt->e_htccz,"e_htccz/F");
  t->Branch("e_ftof1bx",&evnt->e_ftof1bx,"e_ftof1bx/F");
  t->Branch("e_ftof1by",&evnt->e_ftof1by,"e_ftof1by/F");
  t->Branch("e_ftof1bz",&evnt->e_ftof1bz,"e_ftof1bz/F");
  t->Branch("e_ftof2x",&evnt->e_ftof2x,"e_ftof2x/F");
  t->Branch("e_ftof2y",&evnt->e_ftof2y,"e_ftof2y/F");
  t->Branch("e_ftof2z",&evnt->e_ftof2z,"e_ftof2z/F");
  t->Branch("helonline_hel",&evnt->helonline_hel,"helonline_hel/F");
  t->Branch("helonline_helRaw",&evnt->helonline_helRaw,"helonline_helRaw/F");
  t->Branch("helflip_hel",&evnt->helflip_hel,"helflip_hel/F");
  t->Branch("helflip_helRaw",&evnt->helflip_helRaw,"helflip_helRaw/F");
  t->Branch("helflip_event",&evnt->helflip_event,"helflip_event/F");
  t->Branch("e_dc_status",&evnt->e_dc_status,"e_dc_status/F");
  t->Branch("e_dc_ndf",&evnt->e_dc_ndf,"e_dc_ndf/F");
  t->Branch("e_sector_dc",&evnt->e_sector_dc,"e_sector_dc/F");
  t->Branch("e_statPart",&evnt->e_statPart,"e_statPart/F");
  t->Branch("e_DCPx",&evnt->e_DCPx,"e_DCPx/F");
  t->Branch("e_DCPy",&evnt->e_DCPy,"e_DCPy/F");
  t->Branch("e_DCPz",&evnt->e_DCPz,"e_DCPz/F");
  t->Branch("y",&evnt->y,"y/F");
  t->Branch("th_e",&evnt->th_e,"th_e/F");
  t->Branch("phi_e",&evnt->phi_e,"phi_e/F");
  t->Branch("helicRaw",&evnt->helicRaw,"helicRaw/F");
  //// End electron variables ///

  //// FS particles ///
  t->Branch("ThetaPQ",evnt->ThetaPQ,"ThetaPQ[npart]/F");
  t->Branch("PhiPQ",evnt->PhiPQ,"PhiPQ[npart]/F");
  t->Branch("Zh",evnt->Zh,"Zh[npart]/F");
  t->Branch("Pt2",evnt->Pt2,"Pt2[npart]/F");
  t->Branch("Mx2",evnt->Mx2,"Mx2[npart]/F");
  t->Branch("Xf",evnt->Xf,"Xf[npart]/F");
  t->Branch("T",evnt->T,"T[npart]/F");
  t->Branch("P",evnt->P,"P[npart]/F");
  t->Branch("deltaZ",evnt->deltaZ,"deltaZ[npart]/F");
  t->Branch("E",evnt->E,"E[npart]/F");
  t->Branch("Px",evnt->Px,"Px[npart]/F");
  t->Branch("Py",evnt->Py,"Py[npart]/F");
  t->Branch("Pz",evnt->Pz,"Pz[npart]/F");
  t->Branch("Ein",evnt->Ein,"Ein[npart]/F");
  t->Branch("Eout",evnt->Eout,"Eout[npart]/F");
  t->Branch("pid",evnt->pid,"pid[npart]/F");
  t->Branch("Beta",evnt->Beta,"Beta[npart]/F");
  t->Branch("vxh",evnt->vxh,"vxh[npart]/F");
  t->Branch("vyh",evnt->vyh,"vyh[npart]/F");
  t->Branch("vzh",evnt->vzh,"vzh[npart]/F");
  t->Branch("npheltcc",evnt->npheltcc,"npheltcc[npart]/F");
  t->Branch("nphehtcc",evnt->nphehtcc,"nphehtcc[npart]/F");
  t->Branch("chi2pid",evnt->chi2pid,"chi2pid[npart]/F");
  t->Branch("Epcal",evnt->Epcal,"Epcal[npart]/F");
  t->Branch("sector_ltcc",evnt->sector_ltcc,"sector_ltcc[npart]/F");
  t->Branch("sector_htcc",evnt->sector_htcc,"sector_htcc[npart]/F");
  t->Branch("sector_ecal",evnt->sector_ecal,"sector_ecal[npart]/F");
  t->Branch("pcal_lu",evnt->pcal_lu,"pcal_lu[npart]/F");
  t->Branch("pcal_lv",evnt->pcal_lv,"pcal_lv[npart]/F");
  t->Branch("pcal_lw",evnt->pcal_lw,"pcal_lw[npart]/F");
  t->Branch("ecin_lu",evnt->ecin_lu,"ecin_lu[npart]/F");
  t->Branch("ecin_lv",evnt->ecin_lv,"ecin_lv[npart]/F");
  t->Branch("ecin_lw",evnt->ecin_lw,"ecin_lw[npart]/F");
  t->Branch("ecout_lu",evnt->ecout_lu,"ecout_lu[npart]/F");
  t->Branch("ecout_lv",evnt->ecout_lv,"ecout_lv[npart]/F");
  t->Branch("ecout_lw",evnt->ecout_lw,"ecout_lw[npart]/F");
  t->Branch("sector_dc",evnt->sector_dc,"sector_dc[npart]/F");
  t->Branch("statPart",evnt->statPart,"statPart[npart]/F");
  t->Branch("DCPx",evnt->DCPx,"DCPx[npart]/F");
  t->Branch("DCPy",evnt->DCPy,"DCPy[npart]/F");
  t->Branch("DCPz",evnt->DCPz,"DCPz[npart]/F");
  t->Branch("trajx_sl0",evnt->trajx_sl0,"trajx_sl0[npart]/F");
  t->Branch("trajx_sl1",evnt->trajx_sl1,"trajx_sl1[npart]/F");
  t->Branch("trajx_sl2",evnt->trajx_sl2,"trajx_sl2[npart]/F");
  t->Branch("trajx_sl3",evnt->trajx_sl3,"trajx_sl3[npart]/F");
  t->Branch("trajx_sl4",evnt->trajx_sl4,"trajx_sl4[npart]/F");
  t->Branch("trajx_sl5",evnt->trajx_sl5,"trajx_sl5[npart]/F");
  t->Branch("trajy_sl0",evnt->trajy_sl0,"trajy_sl0[npart]/F");
  t->Branch("trajy_sl1",evnt->trajy_sl1,"trajy_sl1[npart]/F");
  t->Branch("trajy_sl2",evnt->trajy_sl2,"trajy_sl2[npart]/F");
  t->Branch("trajy_sl3",evnt->trajy_sl3,"trajy_sl3[npart]/F");
  t->Branch("trajy_sl4",evnt->trajy_sl4,"trajy_sl4[npart]/F");
  t->Branch("trajy_sl5",evnt->trajy_sl5,"trajy_sl5[npart]/F");
  t->Branch("trajz_sl0",evnt->trajz_sl0,"trajz_sl0[npart]/F");
  t->Branch("trajz_sl1",evnt->trajz_sl1,"trajz_sl1[npart]/F");
  t->Branch("trajz_sl2",evnt->trajz_sl2,"trajz_sl2[npart]/F");
  t->Branch("trajz_sl3",evnt->trajz_sl3,"trajz_sl3[npart]/F");
  t->Branch("trajz_sl4",evnt->trajz_sl4,"trajz_sl4[npart]/F");
  t->Branch("trajz_sl5",evnt->trajz_sl5,"trajz_sl5[npart]/F");
  t->Branch("pathtof",evnt->pathtof,"pathtof[npart]/F");
  t->Branch("timetof",evnt->timetof,"timetof[npart]/F");
  t->Branch("sector_tof",evnt->sector_tof,"sector_tof[npart]/F");
  t->Branch("dcx_rot_0",evnt->dcx_rot_0,"dcx_rot_0[npart]/F");
  t->Branch("dcy_rot_0",evnt->dcy_rot_0,"dcy_rot_0[npart]/F");
  t->Branch("dcx_rot_1",evnt->dcx_rot_1,"dcx_rot_1[npart]/F");
  t->Branch("dcy_rot_1",evnt->dcy_rot_1,"dcy_rot_1[npart]/F");
  t->Branch("dcx_rot_2",evnt->dcx_rot_2,"dcx_rot_2[npart]/F");
  t->Branch("dcy_rot_2",evnt->dcy_rot_2,"dcy_rot_2[npart]/F");
  t->Branch("dc_chi2",evnt->dc_chi2,"dc_chi2[npart]/F");
  t->Branch("ftof1ax",evnt->ftof1ax,"ftof1ax[npart]/F");
  t->Branch("ftof1ay",evnt->ftof1ay,"ftof1ay[npart]/F");
  t->Branch("ftof1az",evnt->ftof1az,"ftof1az[npart]/F");
  t->Branch("pcalx",evnt->pcalx,"pcalx[npart]/F");
  t->Branch("pcaly",evnt->pcaly,"pcaly[npart]/F");
  t->Branch("pcalz",evnt->pcalz,"pcalz[npart]/F");
  t->Branch("ecalx",evnt->ecalx,"ecalx[npart]/F");
  t->Branch("ecaly",evnt->ecaly,"ecaly[npart]/F");
  t->Branch("ecalz",evnt->ecalz,"ecalz[npart]/F");
  t->Branch("ltccx",evnt->ltccx,"ltccx[npart]/F");
  t->Branch("ltccy",evnt->ltccy,"ltccy[npart]/F");
  t->Branch("ltccz",evnt->ltccz,"ltccz[npart]/F");
  t->Branch("htccx",evnt->htccx,"htccx[npart]/F");
  t->Branch("htccy",evnt->htccy,"htccy[npart]/F");
  t->Branch("htccz",evnt->htccz,"htccz[npart]/F");
  t->Branch("ftof1bx",evnt->ftof1bx,"ftof1bx[npart]/F");
  t->Branch("ftof1by",evnt->ftof1by,"ftof1by[npart]/F");
  t->Branch("ftof1bz",evnt->ftof1bz,"ftof1bz[npart]/F");
  t->Branch("ftof2x",evnt->ftof2x,"ftof2x[npart]/F");
  t->Branch("ftof2y",evnt->ftof2y,"ftof2y[npart]/F");
  t->Branch("ftof2z",evnt->ftof2z,"ftof2z[npart]/F");
  t->Branch("dc_status",evnt->dc_status,"dc_status[npart]/F");
  t->Branch("dc_ndf",evnt->dc_ndf,"dc_ndf[npart]/F");

  //// End FS particles ///

  //// Electrons MC///
  t->Branch("mc_npart",&evnt->mc_npart,"mc_npart/I");
  t->Branch("mc_Q2",&evnt->mc_Q2,"mc_Q2/F");
  t->Branch("mc_W",&evnt->mc_W,"mc_W/F");
  t->Branch("mc_Nu",&evnt->mc_Nu,"mc_Nu/F");
  t->Branch("mc_Xb",&evnt->mc_Xb,"mc_Xb/F");
  t->Branch("mc_vxe",&evnt->mc_vxe,"mc_vxe/F");
  t->Branch("mc_vye",&evnt->mc_vye,"mc_vye/F");
  t->Branch("mc_vze",&evnt->mc_vze,"mc_vze/F");
  t->Branch("mc_Pex",&evnt->mc_Pex,"mc_Pex/F");
  t->Branch("mc_Pey",&evnt->mc_Pey,"mc_Pey/F");
  t->Branch("mc_Pez",&evnt->mc_Pez,"mc_Pez/F");
  t->Branch("mc_event",&evnt->mc_event,"mc_event/F");
  t->Branch("e_mcmass",&evnt->e_mcmass,"e_mcmass/F");
  t->Branch("mc_Pe",&evnt->mc_Pe,"mc_Pe/F");
  t->Branch("mc_Ee",&evnt->mc_Ee,"mc_Ee/F");
  t->Branch("mc_revent",&evnt->mc_revent,"mc_revent/F");
  t->Branch("mc_y",&evnt->mc_y,"mc_y/F");
  t->Branch("mc_th_e",&evnt->mc_th_e,"mc_th_e/F");
  t->Branch("mc_phi_e",&evnt->mc_phi_e,"mc_phi_e/F");
  t->Branch("mc_e_Beta",&evnt->mc_e_Beta,"mc_e_Beta/F");
  ////  END electron variables MC ///
  
  ////  FS particles MC///
  t->Branch("mc_ThetaPQ",evnt->mc_ThetaPQ,"mc_ThetaPQ[mc_npart]/F");
  t->Branch("mc_PhiPQ",evnt->mc_PhiPQ,"mc_PhiPQ[mc_npart]/F");
  t->Branch("mc_Zh",evnt->mc_Zh,"mc_Zh[mc_npart]/F");
  t->Branch("mc_Pt2",evnt->mc_Pt2,"mc_Pt2[mc_npart]/F");
  t->Branch("mc_Mx2",evnt->mc_Mx2,"mc_Mx2[mc_npart]/F");
  t->Branch("mc_Xf",evnt->mc_Xf,"mc_Xf[mc_npart]/F");
  t->Branch("mc_T",evnt->mc_T,"mc_T[mc_npart]/F");
  t->Branch("mc_P",evnt->mc_P,"mc_P[mc_npart]/F");
  t->Branch("mc_deltaZ",evnt->mc_deltaZ,"mc_deltaZ[mc_npart]/F");
  t->Branch("mc_E",evnt->mc_E,"mc_E[mc_npart]/F");
  t->Branch("mc_Px",evnt->mc_Px,"mc_Px[mc_npart]/F");
  t->Branch("mc_Py",evnt->mc_Py,"mc_Py[mc_npart]/F");
  t->Branch("mc_Pz",evnt->mc_Pz,"mc_Pz[mc_npart]/F");
  t->Branch("mc_pid",evnt->mc_pid,"mc_pid[mc_npart]/F");
  t->Branch("mc_Beta",evnt->mc_Beta,"mc_Beta[mc_npart]/F");
  t->Branch("mc_vxh",evnt->mc_vxh,"mc_vxh[mc_npart]/F");
  t->Branch("mc_vyh",evnt->mc_vyh,"mc_vyh[mc_npart]/F");
  t->Branch("mc_vzh",evnt->mc_vzh,"mc_vzh[mc_npart]/F");
  t->Branch("mcmass",evnt->mcmass,"mcmass[mc_npart]/F");
    ////  END FS particles MC ///
  return 0;
}

int resetDATAMIX(DATAMIX *evnt = 0){

  evnt->npart = DEFAULT_VALUE;
  evnt->Q2 = DEFAULT_VALUE;
  evnt->W = DEFAULT_VALUE;
  evnt->Nu = DEFAULT_VALUE;
  evnt->Xb = DEFAULT_VALUE;
  evnt->vxec = DEFAULT_VALUE;
  evnt->vyec = DEFAULT_VALUE;
  evnt->vzec = DEFAULT_VALUE;
  evnt->vxe = DEFAULT_VALUE;
  evnt->vye = DEFAULT_VALUE;
  evnt->vze = DEFAULT_VALUE;
  evnt->Pex = DEFAULT_VALUE;
  evnt->Pey = DEFAULT_VALUE;
  evnt->Pez = DEFAULT_VALUE;
  evnt->event = DEFAULT_VALUE;
  evnt->Pe = DEFAULT_VALUE;
  evnt->Ee = DEFAULT_VALUE;
  evnt->e_Ein = DEFAULT_VALUE;
  evnt->e_Eout = DEFAULT_VALUE;
  evnt->e_Epcal = DEFAULT_VALUE;
  evnt->e_npheltcc = DEFAULT_VALUE;
  evnt->e_nphehtcc = DEFAULT_VALUE;
  evnt->helic = DEFAULT_VALUE;
  evnt->e_chi2pid = DEFAULT_VALUE;
  evnt->e_pcal_lu = DEFAULT_VALUE;
  evnt->e_pcal_lv = DEFAULT_VALUE;
  evnt->e_pcal_lw = DEFAULT_VALUE;
  evnt->e_ecin_lu = DEFAULT_VALUE;
  evnt->e_ecin_lv = DEFAULT_VALUE;
  evnt->e_ecin_lw = DEFAULT_VALUE;
  evnt->e_ecout_lu = DEFAULT_VALUE;
  evnt->e_ecout_lv = DEFAULT_VALUE;
  evnt->e_ecout_lw = DEFAULT_VALUE;
  evnt->e_pcal_hx = DEFAULT_VALUE;
  evnt->e_pcal_hy = DEFAULT_VALUE;
  evnt->e_pcal_hz = DEFAULT_VALUE;
  evnt->e_ecin_hx = DEFAULT_VALUE;
  evnt->e_ecin_hy = DEFAULT_VALUE;
  evnt->e_ecin_hz = DEFAULT_VALUE;
  evnt->e_ecout_hx = DEFAULT_VALUE;
  evnt->e_ecout_hy = DEFAULT_VALUE;
  evnt->e_ecout_hz = DEFAULT_VALUE;
  evnt->e_trajx_sl0 = DEFAULT_VALUE;
  evnt->e_trajx_sl1 = DEFAULT_VALUE;
  evnt->e_trajx_sl2 = DEFAULT_VALUE;
  evnt->e_trajx_sl3 = DEFAULT_VALUE;
  evnt->e_trajx_sl4 = DEFAULT_VALUE;
  evnt->e_trajx_sl5 = DEFAULT_VALUE;
  evnt->e_trajy_sl0 = DEFAULT_VALUE;
  evnt->e_trajy_sl1 = DEFAULT_VALUE;
  evnt->e_trajy_sl2 = DEFAULT_VALUE;
  evnt->e_trajy_sl3 = DEFAULT_VALUE;
  evnt->e_trajy_sl4 = DEFAULT_VALUE;
  evnt->e_trajy_sl5 = DEFAULT_VALUE;
  evnt->e_trajz_sl0 = DEFAULT_VALUE;
  evnt->e_trajz_sl1 = DEFAULT_VALUE;
  evnt->e_trajz_sl2 = DEFAULT_VALUE;
  evnt->e_trajz_sl3 = DEFAULT_VALUE;
  evnt->e_trajz_sl4 = DEFAULT_VALUE;
  evnt->e_trajz_sl5 = DEFAULT_VALUE;
  evnt->e_pathtof = DEFAULT_VALUE;
  evnt->e_timetof = DEFAULT_VALUE;
  evnt->e_sector_tof = DEFAULT_VALUE;
  evnt->e_Beta = DEFAULT_VALUE;
  evnt->STTime = DEFAULT_VALUE;
  evnt->RFTime = DEFAULT_VALUE;
  evnt->e_dcx_rot_0 = DEFAULT_VALUE;
  evnt->e_dcy_rot_0 = DEFAULT_VALUE;
  evnt->e_dcx_rot_1 = DEFAULT_VALUE;
  evnt->e_dcy_rot_1 = DEFAULT_VALUE;
  evnt->e_dcx_rot_2 = DEFAULT_VALUE;
  evnt->e_dcy_rot_2 = DEFAULT_VALUE;
  evnt->e_sector_ltcc = DEFAULT_VALUE;
  evnt->e_sector_htcc = DEFAULT_VALUE;
  evnt->e_sector_ecal = DEFAULT_VALUE;
  evnt->e_dc_chi2 = DEFAULT_VALUE;
  evnt->e_ftof1ax = DEFAULT_VALUE;
  evnt->e_ftof1ay = DEFAULT_VALUE;
  evnt->e_ftof1az = DEFAULT_VALUE;
  evnt->e_pcalx = DEFAULT_VALUE;
  evnt->e_pcaly = DEFAULT_VALUE;
  evnt->e_pcalz = DEFAULT_VALUE;
  evnt->e_ecalx = DEFAULT_VALUE;
  evnt->e_ecaly = DEFAULT_VALUE;
  evnt->e_ecalz = DEFAULT_VALUE;
  evnt->e_ltccx = DEFAULT_VALUE;
  evnt->e_ltccy = DEFAULT_VALUE;
  evnt->e_ltccz = DEFAULT_VALUE;
  evnt->e_htccx = DEFAULT_VALUE;
  evnt->e_htccy = DEFAULT_VALUE;
  evnt->e_htccz = DEFAULT_VALUE;
  evnt->e_ftof1bx = DEFAULT_VALUE;
  evnt->e_ftof1by = DEFAULT_VALUE;
  evnt->e_ftof1bz = DEFAULT_VALUE;
  evnt->e_ftof2x = DEFAULT_VALUE;
  evnt->e_ftof2y = DEFAULT_VALUE;
  evnt->e_ftof2z = DEFAULT_VALUE;
  evnt->helonline_hel = DEFAULT_VALUE;
  evnt->helonline_helRaw = DEFAULT_VALUE;
  evnt->helflip_hel = DEFAULT_VALUE;
  evnt->helflip_helRaw = DEFAULT_VALUE;
  evnt->helflip_event = DEFAULT_VALUE;
  evnt->e_dc_status = DEFAULT_VALUE;
  evnt->e_dc_ndf = DEFAULT_VALUE;
  evnt->e_sector_dc = DEFAULT_VALUE;
  evnt->e_statPart = DEFAULT_VALUE;
  evnt->e_DCPx = DEFAULT_VALUE;
  evnt->e_DCPy = DEFAULT_VALUE;
  evnt->e_DCPz = DEFAULT_VALUE;
  evnt->revent = DEFAULT_VALUE;
  evnt->y = DEFAULT_VALUE;
  evnt->th_e = DEFAULT_VALUE;
  evnt->phi_e = DEFAULT_VALUE;
  evnt->helicRaw = DEFAULT_VALUE;

  evnt->fA = DEFAULT_VALUE;
  evnt->fB = DEFAULT_VALUE;
  evnt->fC = DEFAULT_VALUE;
  evnt->fV = DEFAULT_VALUE;
  evnt->fW = DEFAULT_VALUE;

  evnt->mix_npart = DEFAULT_VALUE;
  
  //// MC data ///
  evnt->mc_npart = DEFAULT_VALUE;
  evnt->mc_Q2 = DEFAULT_VALUE;
  evnt->mc_W = DEFAULT_VALUE;
  evnt->mc_Nu = DEFAULT_VALUE;
  evnt->mc_Xb = DEFAULT_VALUE;
  evnt->mc_vxe = DEFAULT_VALUE;
  evnt->mc_vye = DEFAULT_VALUE;
  evnt->mc_vze = DEFAULT_VALUE;
  evnt->mc_Pex = DEFAULT_VALUE;
  evnt->mc_Pey = DEFAULT_VALUE;
  evnt->mc_Pez = DEFAULT_VALUE;
  evnt->mc_event = DEFAULT_VALUE;
  evnt->e_mcmass = DEFAULT_VALUE;
  evnt->mc_Pe = DEFAULT_VALUE;
  evnt->mc_Ee = DEFAULT_VALUE;
  evnt->mc_revent = DEFAULT_VALUE;
  evnt->mc_y = DEFAULT_VALUE;
  evnt->mc_th_e = DEFAULT_VALUE;
  evnt->mc_phi_e = DEFAULT_VALUE;
  evnt->mc_e_Beta = DEFAULT_VALUE;
  evnt->mc_helic = DEFAULT_VALUE;
  evnt->mc_fA = DEFAULT_VALUE;
  evnt->mc_fB = DEFAULT_VALUE;
  evnt->mc_fC = DEFAULT_VALUE;
  evnt->mc_fV = DEFAULT_VALUE;
  evnt->mc_fW = DEFAULT_VALUE;

  
  evnt->mc_mix_npart = DEFAULT_VALUE;
  return 0;
}


int resetDATA(DATA *evnt = 0){
  evnt->Q2 = DEFAULT_VALUE;
  evnt->W = DEFAULT_VALUE;
  evnt->Nu = DEFAULT_VALUE;
  evnt->Xb = DEFAULT_VALUE;
  evnt->vxec = DEFAULT_VALUE;
  evnt->vyec = DEFAULT_VALUE;
  evnt->vzec = DEFAULT_VALUE;
  evnt->vxe = DEFAULT_VALUE;
  evnt->vye = DEFAULT_VALUE;
  evnt->vze = DEFAULT_VALUE;
  evnt->Pex = DEFAULT_VALUE;
  evnt->Pey = DEFAULT_VALUE;
  evnt->Pez = DEFAULT_VALUE;
  evnt->event = DEFAULT_VALUE;
  evnt->Pe = DEFAULT_VALUE;
  evnt->Ee = DEFAULT_VALUE;
  evnt->e_Ein = DEFAULT_VALUE;
  evnt->e_Eout = DEFAULT_VALUE;
  evnt->e_Epcal = DEFAULT_VALUE;
  evnt->e_npheltcc = DEFAULT_VALUE;
  evnt->e_nphehtcc = DEFAULT_VALUE;
  evnt->helic = DEFAULT_VALUE;
  evnt->e_chi2pid = DEFAULT_VALUE;
  evnt->e_pcal_lu = DEFAULT_VALUE;
  evnt->e_pcal_lv = DEFAULT_VALUE;
  evnt->e_pcal_lw = DEFAULT_VALUE;
  evnt->e_ecin_lu = DEFAULT_VALUE;
  evnt->e_ecin_lv = DEFAULT_VALUE;
  evnt->e_ecin_lw = DEFAULT_VALUE;
  evnt->e_ecout_lu = DEFAULT_VALUE;
  evnt->e_ecout_lv = DEFAULT_VALUE;
  evnt->e_ecout_lw = DEFAULT_VALUE;
  evnt->e_pcal_hx = DEFAULT_VALUE;
  evnt->e_pcal_hy = DEFAULT_VALUE;
  evnt->e_pcal_hz = DEFAULT_VALUE;
  evnt->e_ecin_hx = DEFAULT_VALUE;
  evnt->e_ecin_hy = DEFAULT_VALUE;
  evnt->e_ecin_hz = DEFAULT_VALUE;
  evnt->e_ecout_hx = DEFAULT_VALUE;
  evnt->e_ecout_hy = DEFAULT_VALUE;
  evnt->e_ecout_hz = DEFAULT_VALUE;
  evnt->e_trajx_sl0 = DEFAULT_VALUE;
  evnt->e_trajx_sl1 = DEFAULT_VALUE;
  evnt->e_trajx_sl2 = DEFAULT_VALUE;
  evnt->e_trajx_sl3 = DEFAULT_VALUE;
  evnt->e_trajx_sl4 = DEFAULT_VALUE;
  evnt->e_trajx_sl5 = DEFAULT_VALUE;
  evnt->e_trajy_sl0 = DEFAULT_VALUE;
  evnt->e_trajy_sl1 = DEFAULT_VALUE;
  evnt->e_trajy_sl2 = DEFAULT_VALUE;
  evnt->e_trajy_sl3 = DEFAULT_VALUE;
  evnt->e_trajy_sl4 = DEFAULT_VALUE;
  evnt->e_trajy_sl5 = DEFAULT_VALUE;
  evnt->e_trajz_sl0 = DEFAULT_VALUE;
  evnt->e_trajz_sl1 = DEFAULT_VALUE;
  evnt->e_trajz_sl2 = DEFAULT_VALUE;
  evnt->e_trajz_sl3 = DEFAULT_VALUE;
  evnt->e_trajz_sl4 = DEFAULT_VALUE;
  evnt->e_trajz_sl5 = DEFAULT_VALUE;
  evnt->e_pathtof = DEFAULT_VALUE;
  evnt->e_timetof = DEFAULT_VALUE;
  evnt->e_sector_tof = DEFAULT_VALUE;
  evnt->e_Beta = DEFAULT_VALUE;
  evnt->STTime = DEFAULT_VALUE;
  evnt->RFTime = DEFAULT_VALUE;
  evnt->e_dcx_rot_0 = DEFAULT_VALUE;
  evnt->e_dcy_rot_0 = DEFAULT_VALUE;
  evnt->e_dcx_rot_1 = DEFAULT_VALUE;
  evnt->e_dcy_rot_1 = DEFAULT_VALUE;
  evnt->e_dcx_rot_2 = DEFAULT_VALUE;
  evnt->e_dcy_rot_2 = DEFAULT_VALUE;
  evnt->e_sector_ltcc = DEFAULT_VALUE;
  evnt->e_sector_htcc = DEFAULT_VALUE;
  evnt->e_sector_ecal = DEFAULT_VALUE;
  evnt->e_dc_chi2 = DEFAULT_VALUE;
  evnt->e_ftof1ax = DEFAULT_VALUE;
  evnt->e_ftof1ay = DEFAULT_VALUE;
  evnt->e_ftof1az = DEFAULT_VALUE;
  evnt->e_pcalx = DEFAULT_VALUE;
  evnt->e_pcaly = DEFAULT_VALUE;
  evnt->e_pcalz = DEFAULT_VALUE;
  evnt->e_ecalx = DEFAULT_VALUE;
  evnt->e_ecaly = DEFAULT_VALUE;
  evnt->e_ecalz = DEFAULT_VALUE;
  evnt->e_ltccx = DEFAULT_VALUE;
  evnt->e_ltccy = DEFAULT_VALUE;
  evnt->e_ltccz = DEFAULT_VALUE;
  evnt->e_htccx = DEFAULT_VALUE;
  evnt->e_htccy = DEFAULT_VALUE;
  evnt->e_htccz = DEFAULT_VALUE;
  evnt->e_ftof1bx = DEFAULT_VALUE;
  evnt->e_ftof1by = DEFAULT_VALUE;
  evnt->e_ftof1bz = DEFAULT_VALUE;
  evnt->e_ftof2x = DEFAULT_VALUE;
  evnt->e_ftof2y = DEFAULT_VALUE;
  evnt->e_ftof2z = DEFAULT_VALUE;
  evnt->helonline_hel = DEFAULT_VALUE;
  evnt->helonline_helRaw = DEFAULT_VALUE;
  evnt->helflip_hel = DEFAULT_VALUE;
  evnt->helflip_helRaw = DEFAULT_VALUE;
  evnt->helflip_event = DEFAULT_VALUE;
  evnt->e_dc_status = DEFAULT_VALUE;
  evnt->e_dc_ndf = DEFAULT_VALUE;
  evnt->e_sector_dc = DEFAULT_VALUE;
  evnt->e_statPart = DEFAULT_VALUE;
  evnt->e_DCPx = DEFAULT_VALUE;
  evnt->e_DCPy = DEFAULT_VALUE;
  evnt->e_DCPz = DEFAULT_VALUE;
  evnt->revent = DEFAULT_VALUE;
  evnt->y = DEFAULT_VALUE;
  evnt->th_e = DEFAULT_VALUE;
  evnt->phi_e = DEFAULT_VALUE;
  evnt->helicRaw = DEFAULT_VALUE;
  //// MC variables ///
  evnt->mc_Q2 = DEFAULT_VALUE;
  evnt->mc_W = DEFAULT_VALUE;
  evnt->mc_Nu = DEFAULT_VALUE;
  evnt->mc_Xb = DEFAULT_VALUE;
  evnt->mc_vxe = DEFAULT_VALUE;
  evnt->mc_vye = DEFAULT_VALUE;
  evnt->mc_vze = DEFAULT_VALUE;
  evnt->mc_Pex = DEFAULT_VALUE;
  evnt->mc_Pey = DEFAULT_VALUE;
  evnt->mc_Pez = DEFAULT_VALUE;
  evnt->mc_event = DEFAULT_VALUE;
  evnt->e_mcmass = DEFAULT_VALUE;
  evnt->mc_Pe = DEFAULT_VALUE;
  evnt->mc_Ee = DEFAULT_VALUE;
  evnt->mc_revent = DEFAULT_VALUE;
  evnt->mc_y = DEFAULT_VALUE;
  evnt->mc_th_e = DEFAULT_VALUE;
  evnt->mc_phi_e = DEFAULT_VALUE;
  evnt->mc_e_Beta = DEFAULT_VALUE;
  return 0;
}

//// Setting address for reading ///
int setTreeAddress(TTree *t, DATA* evnt = 0){
  //// Electron variables ///
  t->SetBranchAddress("npart",&evnt->npart);
  t->SetBranchAddress("Q2",&evnt->Q2);
  t->SetBranchAddress("W",&evnt->W);
  t->SetBranchAddress("Nu",&evnt->Nu);
  t->SetBranchAddress("Xb",&evnt->Xb);
  t->SetBranchAddress("vxec",&evnt->vxec);
  t->SetBranchAddress("vyec",&evnt->vyec);
  t->SetBranchAddress("vzec",&evnt->vzec);
  t->SetBranchAddress("vxe",&evnt->vxe);
  t->SetBranchAddress("vye",&evnt->vye);
  t->SetBranchAddress("vze",&evnt->vze);
  t->SetBranchAddress("Pex",&evnt->Pex);
  t->SetBranchAddress("Pey",&evnt->Pey);
  t->SetBranchAddress("Pez",&evnt->Pez);
  t->SetBranchAddress("event",&evnt->event);
  t->SetBranchAddress("Pe",&evnt->Pe);
  t->SetBranchAddress("Ee",&evnt->Ee);
  t->SetBranchAddress("e_Ein",&evnt->e_Ein);
  t->SetBranchAddress("e_Eout",&evnt->e_Eout);
  t->SetBranchAddress("e_Epcal",&evnt->e_Epcal);
  t->SetBranchAddress("e_npheltcc",&evnt->e_npheltcc);
  t->SetBranchAddress("e_nphehtcc",&evnt->e_nphehtcc);
  t->SetBranchAddress("helic",&evnt->helic);
  t->SetBranchAddress("e_chi2pid",&evnt->e_chi2pid);
  t->SetBranchAddress("e_pcal_lu",&evnt->e_pcal_lu);
  t->SetBranchAddress("e_pcal_lv",&evnt->e_pcal_lv);
  t->SetBranchAddress("e_pcal_lw",&evnt->e_pcal_lw);
  t->SetBranchAddress("e_ecin_lu",&evnt->e_ecin_lu);
  t->SetBranchAddress("e_ecin_lv",&evnt->e_ecin_lv);
  t->SetBranchAddress("e_ecin_lw",&evnt->e_ecin_lw);
  t->SetBranchAddress("e_ecout_lu",&evnt->e_ecout_lu);
  t->SetBranchAddress("e_ecout_lv",&evnt->e_ecout_lv);
  t->SetBranchAddress("e_ecout_lw",&evnt->e_ecout_lw);
  t->SetBranchAddress("e_pcal_hx",&evnt->e_pcal_hx);
  t->SetBranchAddress("e_pcal_hy",&evnt->e_pcal_hy);
  t->SetBranchAddress("e_pcal_hz",&evnt->e_pcal_hz);
  t->SetBranchAddress("e_ecin_hx",&evnt->e_ecin_hx);
  t->SetBranchAddress("e_ecin_hy",&evnt->e_ecin_hy);
  t->SetBranchAddress("e_ecin_hz",&evnt->e_ecin_hz);
  t->SetBranchAddress("e_ecout_hx",&evnt->e_ecout_hx);
  t->SetBranchAddress("e_ecout_hy",&evnt->e_ecout_hy);
  t->SetBranchAddress("e_ecout_hz",&evnt->e_ecout_hz);
  t->SetBranchAddress("e_trajx_sl0",&evnt->e_trajx_sl0);
  t->SetBranchAddress("e_trajx_sl1",&evnt->e_trajx_sl1);
  t->SetBranchAddress("e_trajx_sl2",&evnt->e_trajx_sl2);
  t->SetBranchAddress("e_trajx_sl3",&evnt->e_trajx_sl3);
  t->SetBranchAddress("e_trajx_sl4",&evnt->e_trajx_sl4);
  t->SetBranchAddress("e_trajx_sl5",&evnt->e_trajx_sl5);
  t->SetBranchAddress("e_trajy_sl0",&evnt->e_trajy_sl0);
  t->SetBranchAddress("e_trajy_sl1",&evnt->e_trajy_sl1);
  t->SetBranchAddress("e_trajy_sl2",&evnt->e_trajy_sl2);
  t->SetBranchAddress("e_trajy_sl3",&evnt->e_trajy_sl3);
  t->SetBranchAddress("e_trajy_sl4",&evnt->e_trajy_sl4);
  t->SetBranchAddress("e_trajy_sl5",&evnt->e_trajy_sl5);
  t->SetBranchAddress("e_trajz_sl0",&evnt->e_trajz_sl0);
  t->SetBranchAddress("e_trajz_sl1",&evnt->e_trajz_sl1);
  t->SetBranchAddress("e_trajz_sl2",&evnt->e_trajz_sl2);
  t->SetBranchAddress("e_trajz_sl3",&evnt->e_trajz_sl3);
  t->SetBranchAddress("e_trajz_sl4",&evnt->e_trajz_sl4);
  t->SetBranchAddress("e_trajz_sl5",&evnt->e_trajz_sl5);
  t->SetBranchAddress("e_pathtof",&evnt->e_pathtof);
  t->SetBranchAddress("e_timetof",&evnt->e_timetof);
  t->SetBranchAddress("e_sector_tof",&evnt->e_sector_tof);
  t->SetBranchAddress("e_Beta",&evnt->e_Beta);
  t->SetBranchAddress("STTime",&evnt->STTime);
  t->SetBranchAddress("RFTime",&evnt->RFTime);
  t->SetBranchAddress("e_dcx_rot_0",&evnt->e_dcx_rot_0);
  t->SetBranchAddress("e_dcy_rot_0",&evnt->e_dcy_rot_0);
  t->SetBranchAddress("e_dcx_rot_1",&evnt->e_dcx_rot_1);
  t->SetBranchAddress("e_dcy_rot_1",&evnt->e_dcy_rot_1);
  t->SetBranchAddress("e_dcx_rot_2",&evnt->e_dcx_rot_2);
  t->SetBranchAddress("e_dcy_rot_2",&evnt->e_dcy_rot_2);
  t->SetBranchAddress("e_sector_ltcc",&evnt->e_sector_ltcc);
  t->SetBranchAddress("e_sector_htcc",&evnt->e_sector_htcc);
  t->SetBranchAddress("e_sector_ecal",&evnt->e_sector_ecal);
  t->SetBranchAddress("revent",&evnt->revent);
  t->SetBranchAddress("e_dc_chi2",&evnt->e_dc_chi2);
  t->SetBranchAddress("e_ftof1ax",&evnt->e_ftof1ax);
  t->SetBranchAddress("e_ftof1ay",&evnt->e_ftof1ay);
  t->SetBranchAddress("e_ftof1az",&evnt->e_ftof1az);
  t->SetBranchAddress("e_pcalx",&evnt->e_pcalx);
  t->SetBranchAddress("e_pcaly",&evnt->e_pcaly);
  t->SetBranchAddress("e_pcalz",&evnt->e_pcalz);
  t->SetBranchAddress("e_ecalx",&evnt->e_ecalx);
  t->SetBranchAddress("e_ecaly",&evnt->e_ecaly);
  t->SetBranchAddress("e_ecalz",&evnt->e_ecalz);
  t->SetBranchAddress("e_ltccx",&evnt->e_ltccx);
  t->SetBranchAddress("e_ltccy",&evnt->e_ltccy);
  t->SetBranchAddress("e_ltccz",&evnt->e_ltccz);
  t->SetBranchAddress("e_htccx",&evnt->e_htccx);
  t->SetBranchAddress("e_htccy",&evnt->e_htccy);
  t->SetBranchAddress("e_htccz",&evnt->e_htccz);
  t->SetBranchAddress("e_ftof1bx",&evnt->e_ftof1bx);
  t->SetBranchAddress("e_ftof1by",&evnt->e_ftof1by);
  t->SetBranchAddress("e_ftof1bz",&evnt->e_ftof1bz);
  t->SetBranchAddress("e_ftof2x",&evnt->e_ftof2x);
  t->SetBranchAddress("e_ftof2y",&evnt->e_ftof2y);
  t->SetBranchAddress("e_ftof2z",&evnt->e_ftof2z);
  t->SetBranchAddress("helonline_hel",&evnt->helonline_hel);
  t->SetBranchAddress("helonline_helRaw",&evnt->helonline_helRaw);
  t->SetBranchAddress("helflip_hel",&evnt->helflip_hel);
  t->SetBranchAddress("helflip_helRaw",&evnt->helflip_helRaw);
  t->SetBranchAddress("helflip_event",&evnt->helflip_event);
  t->SetBranchAddress("e_dc_status",&evnt->e_dc_status);
  t->SetBranchAddress("e_dc_ndf",&evnt->e_dc_ndf);
  t->SetBranchAddress("e_sector_dc",&evnt->e_sector_dc);
  t->SetBranchAddress("e_statPart",&evnt->e_statPart);
  t->SetBranchAddress("e_DCPx",&evnt->e_DCPx);
  t->SetBranchAddress("e_DCPy",&evnt->e_DCPy);
  t->SetBranchAddress("e_DCPz",&evnt->e_DCPz);
  t->SetBranchAddress("y",&evnt->y);
  t->SetBranchAddress("th_e",&evnt->th_e);
  t->SetBranchAddress("phi_e",&evnt->phi_e);
  t->SetBranchAddress("helicRaw",&evnt->helicRaw);
  ////  End electron variables ///

  /// FS particles ///
  t->SetBranchAddress("ThetaPQ",evnt->ThetaPQ);
  t->SetBranchAddress("PhiPQ",evnt->PhiPQ);
  t->SetBranchAddress("Zh",evnt->Zh);
  t->SetBranchAddress("Pt2",evnt->Pt2);
  t->SetBranchAddress("Mx2",evnt->Mx2);
  t->SetBranchAddress("Xf",evnt->Xf);
  t->SetBranchAddress("T",evnt->T);
  t->SetBranchAddress("P",evnt->P);
  t->SetBranchAddress("deltaZ",evnt->deltaZ);
  t->SetBranchAddress("E",evnt->E);
  t->SetBranchAddress("Px",evnt->Px);
  t->SetBranchAddress("Py",evnt->Py);
  t->SetBranchAddress("Pz",evnt->Pz);
  t->SetBranchAddress("Ein",evnt->Ein);
  t->SetBranchAddress("Eout",evnt->Eout);
  t->SetBranchAddress("pid",evnt->pid);
  t->SetBranchAddress("Beta",evnt->Beta);
  t->SetBranchAddress("vxh",evnt->vxh);
  t->SetBranchAddress("vyh",evnt->vyh);
  t->SetBranchAddress("vzh",evnt->vzh);
  t->SetBranchAddress("npheltcc",evnt->npheltcc);
  t->SetBranchAddress("nphehtcc",evnt->nphehtcc);
  t->SetBranchAddress("chi2pid",evnt->chi2pid);
  t->SetBranchAddress("Epcal",evnt->Epcal);
  t->SetBranchAddress("sector_ltcc",evnt->sector_ltcc);
  t->SetBranchAddress("sector_htcc",evnt->sector_htcc);
  t->SetBranchAddress("sector_ecal",evnt->sector_ecal);
  t->SetBranchAddress("pcal_lu",evnt->pcal_lu);
  t->SetBranchAddress("pcal_lv",evnt->pcal_lv);
  t->SetBranchAddress("pcal_lw",evnt->pcal_lw);
  t->SetBranchAddress("ecin_lu",evnt->ecin_lu);
  t->SetBranchAddress("ecin_lv",evnt->ecin_lv);
  t->SetBranchAddress("ecin_lw",evnt->ecin_lw);
  t->SetBranchAddress("ecout_lu",evnt->ecout_lu);
  t->SetBranchAddress("ecout_lv",evnt->ecout_lv);
  t->SetBranchAddress("ecout_lw",evnt->ecout_lw);
  t->SetBranchAddress("sector_dc",evnt->sector_dc);
  t->SetBranchAddress("statPart",evnt->statPart);
  t->SetBranchAddress("DCPx",evnt->DCPx);
  t->SetBranchAddress("DCPy",evnt->DCPy);
  t->SetBranchAddress("DCPz",evnt->DCPz);
  t->SetBranchAddress("trajx_sl0",evnt->trajx_sl0);
  t->SetBranchAddress("trajx_sl1",evnt->trajx_sl1);
  t->SetBranchAddress("trajx_sl2",evnt->trajx_sl2);
  t->SetBranchAddress("trajx_sl3",evnt->trajx_sl3);
  t->SetBranchAddress("trajx_sl4",evnt->trajx_sl4);
  t->SetBranchAddress("trajx_sl5",evnt->trajx_sl5);
  t->SetBranchAddress("trajy_sl0",evnt->trajy_sl0);
  t->SetBranchAddress("trajy_sl1",evnt->trajy_sl1);
  t->SetBranchAddress("trajy_sl2",evnt->trajy_sl2);
  t->SetBranchAddress("trajy_sl3",evnt->trajy_sl3);
  t->SetBranchAddress("trajy_sl4",evnt->trajy_sl4);
  t->SetBranchAddress("trajy_sl5",evnt->trajy_sl5);
  t->SetBranchAddress("trajz_sl0",evnt->trajz_sl0);
  t->SetBranchAddress("trajz_sl1",evnt->trajz_sl1);
  t->SetBranchAddress("trajz_sl2",evnt->trajz_sl2);
  t->SetBranchAddress("trajz_sl3",evnt->trajz_sl3);
  t->SetBranchAddress("trajz_sl4",evnt->trajz_sl4);
  t->SetBranchAddress("trajz_sl5",evnt->trajz_sl5);
  t->SetBranchAddress("pathtof",evnt->pathtof);
  t->SetBranchAddress("timetof",evnt->timetof);
  t->SetBranchAddress("sector_tof",evnt->sector_tof);
  t->SetBranchAddress("dcx_rot_0",evnt->dcx_rot_0);
  t->SetBranchAddress("dcy_rot_0",evnt->dcy_rot_0);
  t->SetBranchAddress("dcx_rot_1",evnt->dcx_rot_1);
  t->SetBranchAddress("dcy_rot_1",evnt->dcy_rot_1);
  t->SetBranchAddress("dcx_rot_2",evnt->dcx_rot_2);
  t->SetBranchAddress("dcy_rot_2",evnt->dcy_rot_2);
  t->SetBranchAddress("dc_chi2",evnt->dc_chi2);
  t->SetBranchAddress("ftof1ax",evnt->ftof1ax);
  t->SetBranchAddress("ftof1ay",evnt->ftof1ay);
  t->SetBranchAddress("ftof1az",evnt->ftof1az);
  t->SetBranchAddress("pcalx",evnt->pcalx);
  t->SetBranchAddress("pcaly",evnt->pcaly);
  t->SetBranchAddress("pcalz",evnt->pcalz);
  t->SetBranchAddress("ecalx",evnt->ecalx);
  t->SetBranchAddress("ecaly",evnt->ecaly);
  t->SetBranchAddress("ecalz",evnt->ecalz);
  t->SetBranchAddress("ltccx",evnt->ltccx);
  t->SetBranchAddress("ltccy",evnt->ltccy);
  t->SetBranchAddress("ltccz",evnt->ltccz);
  t->SetBranchAddress("htccx",evnt->htccx);
  t->SetBranchAddress("htccy",evnt->htccy);
  t->SetBranchAddress("htccz",evnt->htccz);
  t->SetBranchAddress("ftof1bx",evnt->ftof1bx);
  t->SetBranchAddress("ftof1by",evnt->ftof1by);
  t->SetBranchAddress("ftof1bz",evnt->ftof1bz);
  t->SetBranchAddress("ftof2x",evnt->ftof2x);
  t->SetBranchAddress("ftof2y",evnt->ftof2y);
  t->SetBranchAddress("ftof2z",evnt->ftof2z);
  t->SetBranchAddress("dc_status",evnt->dc_status);
  t->SetBranchAddress("dc_ndf",evnt->dc_ndf);

  //// End FS particles ///

  //// Electrons MC///
  t->SetBranchAddress("mc_npart",&evnt->mc_npart);
  t->SetBranchAddress("mc_Q2",&evnt->mc_Q2);
  t->SetBranchAddress("mc_W",&evnt->mc_W);
  t->SetBranchAddress("mc_Nu",&evnt->mc_Nu);
  t->SetBranchAddress("mc_Xb",&evnt->mc_Xb);
  t->SetBranchAddress("mc_vxe",&evnt->mc_vxe);
  t->SetBranchAddress("mc_vye",&evnt->mc_vye);
  t->SetBranchAddress("mc_vze",&evnt->mc_vze);
  t->SetBranchAddress("mc_Pex",&evnt->mc_Pex);
  t->SetBranchAddress("mc_Pey",&evnt->mc_Pey);
  t->SetBranchAddress("mc_Pez",&evnt->mc_Pez);
  t->SetBranchAddress("mc_event",&evnt->mc_event);
  t->SetBranchAddress("e_mcmass",&evnt->e_mcmass);
  t->SetBranchAddress("mc_Pe",&evnt->mc_Pe);
  t->SetBranchAddress("mc_Ee",&evnt->mc_Ee);
  t->SetBranchAddress("mc_revent",&evnt->mc_revent);
  t->SetBranchAddress("mc_y",&evnt->mc_y);
  t->SetBranchAddress("mc_th_e",&evnt->mc_th_e);
  t->SetBranchAddress("mc_phi_e",&evnt->mc_phi_e);
  t->SetBranchAddress("mc_e_Beta",&evnt->mc_e_Beta);
  ////  END electron variables MC ///
  
  ////  FS particles MC ///
  t->SetBranchAddress("mc_ThetaPQ",evnt->mc_ThetaPQ);
  t->SetBranchAddress("mc_PhiPQ",evnt->mc_PhiPQ);
  t->SetBranchAddress("mc_Zh",evnt->mc_Zh);
  t->SetBranchAddress("mc_Pt2",evnt->mc_Pt2);
  t->SetBranchAddress("mc_Mx2",evnt->mc_Mx2);
  t->SetBranchAddress("mc_Xf",evnt->mc_Xf);
  t->SetBranchAddress("mc_T",evnt->mc_T);
  t->SetBranchAddress("mc_P",evnt->mc_P);
  t->SetBranchAddress("mc_deltaZ",evnt->mc_deltaZ);
  t->SetBranchAddress("mc_E",evnt->mc_E);
  t->SetBranchAddress("mc_Px",evnt->mc_Px);
  t->SetBranchAddress("mc_Py",evnt->mc_Py);
  t->SetBranchAddress("mc_Pz",evnt->mc_Pz);
  t->SetBranchAddress("mc_pid",evnt->mc_pid);
  t->SetBranchAddress("mc_Beta",evnt->mc_Beta);
  t->SetBranchAddress("mc_vxh",evnt->mc_vxh);
  t->SetBranchAddress("mc_vyh",evnt->mc_vyh);
  t->SetBranchAddress("mc_vzh",evnt->mc_vzh);
  t->SetBranchAddress("mcmass",evnt->mcmass);
    ////  END FS particles MC ///
  return 0;
}
//// END setting address///


#endif
