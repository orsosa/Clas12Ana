//// File automatically produced by format_hiponodes.py do not make changes here!!
#include "TIdentificatorCLAS12.h"

int TIdentificatorCLAS12::get_BMTRec__Hits(int row){
	BMTRec__Hits_ID = BMTRec__Hits->getShort("ID",row);
	BMTRec__Hits_sector = BMTRec__Hits->getByte("sector",row);
	BMTRec__Hits_layer = BMTRec__Hits->getByte("layer",row);
	BMTRec__Hits_strip = BMTRec__Hits->getInt("strip",row);
	BMTRec__Hits_fitResidual = BMTRec__Hits->getFloat("fitResidual",row);
	BMTRec__Hits_trkingStat = BMTRec__Hits->getInt("trkingStat",row);
	BMTRec__Hits_clusterID = BMTRec__Hits->getShort("clusterID",row);
	BMTRec__Hits_trkID = BMTRec__Hits->getShort("trkID",row);
	return 0;
} 

int TIdentificatorCLAS12::get_RAW__adc(int row){
	RAW__adc_crate = RAW__adc->getByte("crate",row);
	RAW__adc_slot = RAW__adc->getByte("slot",row);
	RAW__adc_channel = RAW__adc->getShort("channel",row);
	RAW__adc_order = RAW__adc->getByte("order",row);
	RAW__adc_ADC = RAW__adc->getInt("ADC",row);
	RAW__adc_time = RAW__adc->getFloat("time",row);
	RAW__adc_ped = RAW__adc->getShort("ped",row);
	return 0;
} 

int TIdentificatorCLAS12::get_BAND__adc(int row){
	BAND__adc_sector = BAND__adc->getByte("sector",row);
	BAND__adc_layer = BAND__adc->getByte("layer",row);
	BAND__adc_component = BAND__adc->getShort("component",row);
	BAND__adc_order = BAND__adc->getByte("order",row);
	BAND__adc_ADC = BAND__adc->getInt("ADC",row);
	BAND__adc_time = BAND__adc->getFloat("time",row);
	BAND__adc_ped = BAND__adc->getShort("ped",row);
	return 0;
} 

int TIdentificatorCLAS12::get_RUN__config(int row){
	RUN__config_run = RUN__config->getInt("run",row);
	RUN__config_event = RUN__config->getInt("event",row);
	RUN__config_unixtime = RUN__config->getInt("unixtime",row);
	RUN__config_trigger = RUN__config->getLong("trigger",row);
	RUN__config_timestamp = RUN__config->getLong("timestamp",row);
	RUN__config_type = RUN__config->getByte("type",row);
	RUN__config_mode = RUN__config->getByte("mode",row);
	RUN__config_torus = RUN__config->getFloat("torus",row);
	RUN__config_solenoid = RUN__config->getFloat("solenoid",row);
	return 0;
} 

int TIdentificatorCLAS12::get_RICH__clusters(int row){
	RICH__clusters_id = RICH__clusters->getShort("id",row);
	RICH__clusters_size = RICH__clusters->getShort("size",row);
	RICH__clusters_sector = RICH__clusters->getShort("sector",row);
	RICH__clusters_tile = RICH__clusters->getShort("tile",row);
	RICH__clusters_pmt = RICH__clusters->getShort("pmt",row);
	RICH__clusters_charge = RICH__clusters->getFloat("charge",row);
	RICH__clusters_time = RICH__clusters->getFloat("time",row);
	RICH__clusters_rawtime = RICH__clusters->getFloat("rawtime",row);
	RICH__clusters_x = RICH__clusters->getFloat("x",row);
	RICH__clusters_y = RICH__clusters->getFloat("y",row);
	RICH__clusters_z = RICH__clusters->getFloat("z",row);
	RICH__clusters_wtime = RICH__clusters->getFloat("wtime",row);
	RICH__clusters_wx = RICH__clusters->getFloat("wx",row);
	RICH__clusters_wy = RICH__clusters->getFloat("wy",row);
	RICH__clusters_wz = RICH__clusters->getFloat("wz",row);
	return 0;
} 

int TIdentificatorCLAS12::get_RECHB__Scintillator(int row){
	RECHB__Scintillator_index = RECHB__Scintillator->getShort("index",row);
	RECHB__Scintillator_pindex = RECHB__Scintillator->getShort("pindex",row);
	RECHB__Scintillator_detector = RECHB__Scintillator->getByte("detector",row);
	RECHB__Scintillator_sector = RECHB__Scintillator->getByte("sector",row);
	RECHB__Scintillator_layer = RECHB__Scintillator->getByte("layer",row);
	RECHB__Scintillator_component = RECHB__Scintillator->getShort("component",row);
	RECHB__Scintillator_energy = RECHB__Scintillator->getFloat("energy",row);
	RECHB__Scintillator_time = RECHB__Scintillator->getFloat("time",row);
	RECHB__Scintillator_path = RECHB__Scintillator->getFloat("path",row);
	RECHB__Scintillator_chi2 = RECHB__Scintillator->getFloat("chi2",row);
	RECHB__Scintillator_x = RECHB__Scintillator->getFloat("x",row);
	RECHB__Scintillator_y = RECHB__Scintillator->getFloat("y",row);
	RECHB__Scintillator_z = RECHB__Scintillator->getFloat("z",row);
	RECHB__Scintillator_hx = RECHB__Scintillator->getFloat("hx",row);
	RECHB__Scintillator_hy = RECHB__Scintillator->getFloat("hy",row);
	RECHB__Scintillator_hz = RECHB__Scintillator->getFloat("hz",row);
	RECHB__Scintillator_status = RECHB__Scintillator->getShort("status",row);
	return 0;
} 

int TIdentificatorCLAS12::get_REC__RingCher(int row){
	REC__RingCher_id = REC__RingCher->getShort("id",row);
	REC__RingCher_hindex = REC__RingCher->getShort("hindex",row);
	REC__RingCher_pindex = REC__RingCher->getShort("pindex",row);
	REC__RingCher_apath = REC__RingCher->getFloat("apath",row);
	REC__RingCher_atime = REC__RingCher->getFloat("atime",row);
	REC__RingCher_aEtaC = REC__RingCher->getFloat("aEtaC",row);
	REC__RingCher_tpath = REC__RingCher->getFloat("tpath",row);
	REC__RingCher_ttime = REC__RingCher->getFloat("ttime",row);
	REC__RingCher_tEtaC = REC__RingCher->getFloat("tEtaC",row);
	return 0;
} 

int TIdentificatorCLAS12::get_BSTRec__LayerEffs(int row){
	BSTRec__LayerEffs_sector = BSTRec__LayerEffs->getByte("sector",row);
	BSTRec__LayerEffs_layer = BSTRec__LayerEffs->getByte("layer",row);
	BSTRec__LayerEffs_residual = BSTRec__LayerEffs->getFloat("residual",row);
	BSTRec__LayerEffs_status = BSTRec__LayerEffs->getByte("status",row);
	return 0;
} 

int TIdentificatorCLAS12::get_RTPC__pos(int row){
	RTPC__pos_step = RTPC__pos->getFloat("step",row);
	RTPC__pos_time = RTPC__pos->getFloat("time",row);
	RTPC__pos_energy = RTPC__pos->getFloat("energy",row);
	RTPC__pos_posx = RTPC__pos->getFloat("posx",row);
	RTPC__pos_posy = RTPC__pos->getFloat("posy",row);
	RTPC__pos_posz = RTPC__pos->getFloat("posz",row);
	RTPC__pos_phi = RTPC__pos->getFloat("phi",row);
	return 0;
} 

int TIdentificatorCLAS12::get_TimeBasedTrkg__TBCrosses(int row){
	TimeBasedTrkg__TBCrosses_id = TimeBasedTrkg__TBCrosses->getShort("id",row);
	TimeBasedTrkg__TBCrosses_status = TimeBasedTrkg__TBCrosses->getShort("status",row);
	TimeBasedTrkg__TBCrosses_sector = TimeBasedTrkg__TBCrosses->getByte("sector",row);
	TimeBasedTrkg__TBCrosses_region = TimeBasedTrkg__TBCrosses->getByte("region",row);
	TimeBasedTrkg__TBCrosses_x = TimeBasedTrkg__TBCrosses->getFloat("x",row);
	TimeBasedTrkg__TBCrosses_y = TimeBasedTrkg__TBCrosses->getFloat("y",row);
	TimeBasedTrkg__TBCrosses_z = TimeBasedTrkg__TBCrosses->getFloat("z",row);
	TimeBasedTrkg__TBCrosses_err_x = TimeBasedTrkg__TBCrosses->getFloat("err_x",row);
	TimeBasedTrkg__TBCrosses_err_y = TimeBasedTrkg__TBCrosses->getFloat("err_y",row);
	TimeBasedTrkg__TBCrosses_err_z = TimeBasedTrkg__TBCrosses->getFloat("err_z",row);
	TimeBasedTrkg__TBCrosses_ux = TimeBasedTrkg__TBCrosses->getFloat("ux",row);
	TimeBasedTrkg__TBCrosses_uy = TimeBasedTrkg__TBCrosses->getFloat("uy",row);
	TimeBasedTrkg__TBCrosses_uz = TimeBasedTrkg__TBCrosses->getFloat("uz",row);
	TimeBasedTrkg__TBCrosses_err_ux = TimeBasedTrkg__TBCrosses->getFloat("err_ux",row);
	TimeBasedTrkg__TBCrosses_err_uy = TimeBasedTrkg__TBCrosses->getFloat("err_uy",row);
	TimeBasedTrkg__TBCrosses_err_uz = TimeBasedTrkg__TBCrosses->getFloat("err_uz",row);
	TimeBasedTrkg__TBCrosses_Segment1_ID = TimeBasedTrkg__TBCrosses->getShort("Segment1_ID",row);
	TimeBasedTrkg__TBCrosses_Segment2_ID = TimeBasedTrkg__TBCrosses->getShort("Segment2_ID",row);
	return 0;
} 

int TIdentificatorCLAS12::get_HitBasedTrkg__HBTracks(int row){
	HitBasedTrkg__HBTracks_id = HitBasedTrkg__HBTracks->getShort("id",row);
	HitBasedTrkg__HBTracks_status = HitBasedTrkg__HBTracks->getShort("status",row);
	HitBasedTrkg__HBTracks_sector = HitBasedTrkg__HBTracks->getByte("sector",row);
	HitBasedTrkg__HBTracks_c1_x = HitBasedTrkg__HBTracks->getFloat("c1_x",row);
	HitBasedTrkg__HBTracks_c1_y = HitBasedTrkg__HBTracks->getFloat("c1_y",row);
	HitBasedTrkg__HBTracks_c1_z = HitBasedTrkg__HBTracks->getFloat("c1_z",row);
	HitBasedTrkg__HBTracks_c1_ux = HitBasedTrkg__HBTracks->getFloat("c1_ux",row);
	HitBasedTrkg__HBTracks_c1_uy = HitBasedTrkg__HBTracks->getFloat("c1_uy",row);
	HitBasedTrkg__HBTracks_c1_uz = HitBasedTrkg__HBTracks->getFloat("c1_uz",row);
	HitBasedTrkg__HBTracks_c3_x = HitBasedTrkg__HBTracks->getFloat("c3_x",row);
	HitBasedTrkg__HBTracks_c3_y = HitBasedTrkg__HBTracks->getFloat("c3_y",row);
	HitBasedTrkg__HBTracks_c3_z = HitBasedTrkg__HBTracks->getFloat("c3_z",row);
	HitBasedTrkg__HBTracks_c3_ux = HitBasedTrkg__HBTracks->getFloat("c3_ux",row);
	HitBasedTrkg__HBTracks_c3_uy = HitBasedTrkg__HBTracks->getFloat("c3_uy",row);
	HitBasedTrkg__HBTracks_c3_uz = HitBasedTrkg__HBTracks->getFloat("c3_uz",row);
	HitBasedTrkg__HBTracks_t1_x = HitBasedTrkg__HBTracks->getFloat("t1_x",row);
	HitBasedTrkg__HBTracks_t1_y = HitBasedTrkg__HBTracks->getFloat("t1_y",row);
	HitBasedTrkg__HBTracks_t1_z = HitBasedTrkg__HBTracks->getFloat("t1_z",row);
	HitBasedTrkg__HBTracks_t1_px = HitBasedTrkg__HBTracks->getFloat("t1_px",row);
	HitBasedTrkg__HBTracks_t1_py = HitBasedTrkg__HBTracks->getFloat("t1_py",row);
	HitBasedTrkg__HBTracks_t1_pz = HitBasedTrkg__HBTracks->getFloat("t1_pz",row);
	HitBasedTrkg__HBTracks_Vtx0_x = HitBasedTrkg__HBTracks->getFloat("Vtx0_x",row);
	HitBasedTrkg__HBTracks_Vtx0_y = HitBasedTrkg__HBTracks->getFloat("Vtx0_y",row);
	HitBasedTrkg__HBTracks_Vtx0_z = HitBasedTrkg__HBTracks->getFloat("Vtx0_z",row);
	HitBasedTrkg__HBTracks_p0_x = HitBasedTrkg__HBTracks->getFloat("p0_x",row);
	HitBasedTrkg__HBTracks_p0_y = HitBasedTrkg__HBTracks->getFloat("p0_y",row);
	HitBasedTrkg__HBTracks_p0_z = HitBasedTrkg__HBTracks->getFloat("p0_z",row);
	HitBasedTrkg__HBTracks_x = HitBasedTrkg__HBTracks->getFloat("x",row);
	HitBasedTrkg__HBTracks_y = HitBasedTrkg__HBTracks->getFloat("y",row);
	HitBasedTrkg__HBTracks_z = HitBasedTrkg__HBTracks->getFloat("z",row);
	HitBasedTrkg__HBTracks_tx = HitBasedTrkg__HBTracks->getFloat("tx",row);
	HitBasedTrkg__HBTracks_ty = HitBasedTrkg__HBTracks->getFloat("ty",row);
	HitBasedTrkg__HBTracks_Cross1_ID = HitBasedTrkg__HBTracks->getShort("Cross1_ID",row);
	HitBasedTrkg__HBTracks_Cross2_ID = HitBasedTrkg__HBTracks->getShort("Cross2_ID",row);
	HitBasedTrkg__HBTracks_Cross3_ID = HitBasedTrkg__HBTracks->getShort("Cross3_ID",row);
	HitBasedTrkg__HBTracks_q = HitBasedTrkg__HBTracks->getByte("q",row);
	HitBasedTrkg__HBTracks_pathlength = HitBasedTrkg__HBTracks->getFloat("pathlength",row);
	HitBasedTrkg__HBTracks_chi2 = HitBasedTrkg__HBTracks->getFloat("chi2",row);
	HitBasedTrkg__HBTracks_ndf = HitBasedTrkg__HBTracks->getShort("ndf",row);
	return 0;
} 

int TIdentificatorCLAS12::get_CVTRec__Cosmics(int row){
	CVTRec__Cosmics_ID = CVTRec__Cosmics->getShort("ID",row);
	CVTRec__Cosmics_trkline_yx_slope = CVTRec__Cosmics->getFloat("trkline_yx_slope",row);
	CVTRec__Cosmics_trkline_yx_interc = CVTRec__Cosmics->getFloat("trkline_yx_interc",row);
	CVTRec__Cosmics_trkline_yz_slope = CVTRec__Cosmics->getFloat("trkline_yz_slope",row);
	CVTRec__Cosmics_trkline_yz_interc = CVTRec__Cosmics->getFloat("trkline_yz_interc",row);
	CVTRec__Cosmics_theta = CVTRec__Cosmics->getFloat("theta",row);
	CVTRec__Cosmics_phi = CVTRec__Cosmics->getFloat("phi",row);
	CVTRec__Cosmics_chi2 = CVTRec__Cosmics->getFloat("chi2",row);
	CVTRec__Cosmics_ndf = CVTRec__Cosmics->getShort("ndf",row);
	CVTRec__Cosmics_Cross1_ID = CVTRec__Cosmics->getShort("Cross1_ID",row);
	CVTRec__Cosmics_Cross2_ID = CVTRec__Cosmics->getShort("Cross2_ID",row);
	CVTRec__Cosmics_Cross3_ID = CVTRec__Cosmics->getShort("Cross3_ID",row);
	CVTRec__Cosmics_Cross4_ID = CVTRec__Cosmics->getShort("Cross4_ID",row);
	CVTRec__Cosmics_Cross5_ID = CVTRec__Cosmics->getShort("Cross5_ID",row);
	CVTRec__Cosmics_Cross6_ID = CVTRec__Cosmics->getShort("Cross6_ID",row);
	CVTRec__Cosmics_Cross7_ID = CVTRec__Cosmics->getShort("Cross7_ID",row);
	CVTRec__Cosmics_Cross8_ID = CVTRec__Cosmics->getShort("Cross8_ID",row);
	CVTRec__Cosmics_Cross9_ID = CVTRec__Cosmics->getShort("Cross9_ID",row);
	CVTRec__Cosmics_Cross10_ID = CVTRec__Cosmics->getShort("Cross10_ID",row);
	CVTRec__Cosmics_Cross11_ID = CVTRec__Cosmics->getShort("Cross11_ID",row);
	CVTRec__Cosmics_Cross12_ID = CVTRec__Cosmics->getShort("Cross12_ID",row);
	CVTRec__Cosmics_Cross13_ID = CVTRec__Cosmics->getShort("Cross13_ID",row);
	CVTRec__Cosmics_Cross14_ID = CVTRec__Cosmics->getShort("Cross14_ID",row);
	CVTRec__Cosmics_Cross15_ID = CVTRec__Cosmics->getShort("Cross15_ID",row);
	CVTRec__Cosmics_Cross16_ID = CVTRec__Cosmics->getShort("Cross16_ID",row);
	CVTRec__Cosmics_Cross17_ID = CVTRec__Cosmics->getShort("Cross17_ID",row);
	CVTRec__Cosmics_Cross18_ID = CVTRec__Cosmics->getShort("Cross18_ID",row);
	return 0;
} 

int TIdentificatorCLAS12::get_RECHB__Cherenkov(int row){
	RECHB__Cherenkov_index = RECHB__Cherenkov->getShort("index",row);
	RECHB__Cherenkov_pindex = RECHB__Cherenkov->getShort("pindex",row);
	RECHB__Cherenkov_detector = RECHB__Cherenkov->getByte("detector",row);
	RECHB__Cherenkov_sector = RECHB__Cherenkov->getByte("sector",row);
	RECHB__Cherenkov_nphe = RECHB__Cherenkov->getFloat("nphe",row);
	RECHB__Cherenkov_time = RECHB__Cherenkov->getFloat("time",row);
	RECHB__Cherenkov_path = RECHB__Cherenkov->getFloat("path",row);
	RECHB__Cherenkov_chi2 = RECHB__Cherenkov->getFloat("chi2",row);
	RECHB__Cherenkov_x = RECHB__Cherenkov->getFloat("x",row);
	RECHB__Cherenkov_y = RECHB__Cherenkov->getFloat("y",row);
	RECHB__Cherenkov_z = RECHB__Cherenkov->getFloat("z",row);
	RECHB__Cherenkov_dtheta = RECHB__Cherenkov->getFloat("dtheta",row);
	RECHB__Cherenkov_dphi = RECHB__Cherenkov->getFloat("dphi",row);
	RECHB__Cherenkov_status = RECHB__Cherenkov->getShort("status",row);
	return 0;
} 

int TIdentificatorCLAS12::get_BSTRec__Clusters(int row){
	BSTRec__Clusters_ID = BSTRec__Clusters->getShort("ID",row);
	BSTRec__Clusters_sector = BSTRec__Clusters->getByte("sector",row);
	BSTRec__Clusters_layer = BSTRec__Clusters->getByte("layer",row);
	BSTRec__Clusters_size = BSTRec__Clusters->getShort("size",row);
	BSTRec__Clusters_ETot = BSTRec__Clusters->getFloat("ETot",row);
	BSTRec__Clusters_seedE = BSTRec__Clusters->getFloat("seedE",row);
	BSTRec__Clusters_seedStrip = BSTRec__Clusters->getInt("seedStrip",row);
	BSTRec__Clusters_centroid = BSTRec__Clusters->getFloat("centroid",row);
	BSTRec__Clusters_centroidResidual = BSTRec__Clusters->getFloat("centroidResidual",row);
	BSTRec__Clusters_seedResidual = BSTRec__Clusters->getFloat("seedResidual",row);
	BSTRec__Clusters_Hit1_ID = BSTRec__Clusters->getShort("Hit1_ID",row);
	BSTRec__Clusters_Hit2_ID = BSTRec__Clusters->getShort("Hit2_ID",row);
	BSTRec__Clusters_Hit3_ID = BSTRec__Clusters->getShort("Hit3_ID",row);
	BSTRec__Clusters_Hit4_ID = BSTRec__Clusters->getShort("Hit4_ID",row);
	BSTRec__Clusters_Hit5_ID = BSTRec__Clusters->getShort("Hit5_ID",row);
	BSTRec__Clusters_trkID = BSTRec__Clusters->getShort("trkID",row);
	return 0;
} 

int TIdentificatorCLAS12::get_CVTRec__Trajectory(int row){
	CVTRec__Trajectory_id = CVTRec__Trajectory->getShort("id",row);
	CVTRec__Trajectory_detector = CVTRec__Trajectory->getShort("detector",row);
	CVTRec__Trajectory_sector = CVTRec__Trajectory->getByte("sector",row);
	CVTRec__Trajectory_layer = CVTRec__Trajectory->getByte("layer",row);
	CVTRec__Trajectory_x = CVTRec__Trajectory->getFloat("x",row);
	CVTRec__Trajectory_y = CVTRec__Trajectory->getFloat("y",row);
	CVTRec__Trajectory_z = CVTRec__Trajectory->getFloat("z",row);
	CVTRec__Trajectory_phi = CVTRec__Trajectory->getFloat("phi",row);
	CVTRec__Trajectory_theta = CVTRec__Trajectory->getFloat("theta",row);
	CVTRec__Trajectory_langle = CVTRec__Trajectory->getFloat("langle",row);
	CVTRec__Trajectory_centroid = CVTRec__Trajectory->getFloat("centroid",row);
	CVTRec__Trajectory_path = CVTRec__Trajectory->getFloat("path",row);
	return 0;
} 

int TIdentificatorCLAS12::get_RECHB__Calorimeter(int row){
	RECHB__Calorimeter_index = RECHB__Calorimeter->getShort("index",row);
	RECHB__Calorimeter_pindex = RECHB__Calorimeter->getShort("pindex",row);
	RECHB__Calorimeter_detector = RECHB__Calorimeter->getByte("detector",row);
	RECHB__Calorimeter_sector = RECHB__Calorimeter->getByte("sector",row);
	RECHB__Calorimeter_layer = RECHB__Calorimeter->getByte("layer",row);
	RECHB__Calorimeter_energy = RECHB__Calorimeter->getFloat("energy",row);
	RECHB__Calorimeter_time = RECHB__Calorimeter->getFloat("time",row);
	RECHB__Calorimeter_path = RECHB__Calorimeter->getFloat("path",row);
	RECHB__Calorimeter_chi2 = RECHB__Calorimeter->getFloat("chi2",row);
	RECHB__Calorimeter_x = RECHB__Calorimeter->getFloat("x",row);
	RECHB__Calorimeter_y = RECHB__Calorimeter->getFloat("y",row);
	RECHB__Calorimeter_z = RECHB__Calorimeter->getFloat("z",row);
	RECHB__Calorimeter_hx = RECHB__Calorimeter->getFloat("hx",row);
	RECHB__Calorimeter_hy = RECHB__Calorimeter->getFloat("hy",row);
	RECHB__Calorimeter_hz = RECHB__Calorimeter->getFloat("hz",row);
	RECHB__Calorimeter_lu = RECHB__Calorimeter->getFloat("lu",row);
	RECHB__Calorimeter_lv = RECHB__Calorimeter->getFloat("lv",row);
	RECHB__Calorimeter_lw = RECHB__Calorimeter->getFloat("lw",row);
	RECHB__Calorimeter_du = RECHB__Calorimeter->getFloat("du",row);
	RECHB__Calorimeter_dv = RECHB__Calorimeter->getFloat("dv",row);
	RECHB__Calorimeter_dw = RECHB__Calorimeter->getFloat("dw",row);
	RECHB__Calorimeter_m2u = RECHB__Calorimeter->getFloat("m2u",row);
	RECHB__Calorimeter_m2v = RECHB__Calorimeter->getFloat("m2v",row);
	RECHB__Calorimeter_m2w = RECHB__Calorimeter->getFloat("m2w",row);
	RECHB__Calorimeter_m3u = RECHB__Calorimeter->getFloat("m3u",row);
	RECHB__Calorimeter_m3v = RECHB__Calorimeter->getFloat("m3v",row);
	RECHB__Calorimeter_m3w = RECHB__Calorimeter->getFloat("m3w",row);
	RECHB__Calorimeter_status = RECHB__Calorimeter->getShort("status",row);
	return 0;
} 

int TIdentificatorCLAS12::get_TimeBasedTrkg__TBSegmentTrajectory(int row){
	TimeBasedTrkg__TBSegmentTrajectory_segmentID = TimeBasedTrkg__TBSegmentTrajectory->getShort("segmentID",row);
	TimeBasedTrkg__TBSegmentTrajectory_sector = TimeBasedTrkg__TBSegmentTrajectory->getByte("sector",row);
	TimeBasedTrkg__TBSegmentTrajectory_superlayer = TimeBasedTrkg__TBSegmentTrajectory->getByte("superlayer",row);
	TimeBasedTrkg__TBSegmentTrajectory_layer = TimeBasedTrkg__TBSegmentTrajectory->getByte("layer",row);
	TimeBasedTrkg__TBSegmentTrajectory_matchedHitID = TimeBasedTrkg__TBSegmentTrajectory->getShort("matchedHitID",row);
	TimeBasedTrkg__TBSegmentTrajectory_trkDoca = TimeBasedTrkg__TBSegmentTrajectory->getFloat("trkDoca",row);
	return 0;
} 

int TIdentificatorCLAS12::get_CTOF__tdc(int row){
	CTOF__tdc_sector = CTOF__tdc->getByte("sector",row);
	CTOF__tdc_layer = CTOF__tdc->getByte("layer",row);
	CTOF__tdc_component = CTOF__tdc->getShort("component",row);
	CTOF__tdc_order = CTOF__tdc->getByte("order",row);
	CTOF__tdc_TDC = CTOF__tdc->getInt("TDC",row);
	return 0;
} 

int TIdentificatorCLAS12::get_REC__Cherenkov(int row){
	REC__Cherenkov_index = REC__Cherenkov->getShort("index",row);
	REC__Cherenkov_pindex = REC__Cherenkov->getShort("pindex",row);
	REC__Cherenkov_detector = REC__Cherenkov->getByte("detector",row);
	REC__Cherenkov_sector = REC__Cherenkov->getByte("sector",row);
	REC__Cherenkov_nphe = REC__Cherenkov->getFloat("nphe",row);
	REC__Cherenkov_time = REC__Cherenkov->getFloat("time",row);
	REC__Cherenkov_path = REC__Cherenkov->getFloat("path",row);
	REC__Cherenkov_chi2 = REC__Cherenkov->getFloat("chi2",row);
	REC__Cherenkov_x = REC__Cherenkov->getFloat("x",row);
	REC__Cherenkov_y = REC__Cherenkov->getFloat("y",row);
	REC__Cherenkov_z = REC__Cherenkov->getFloat("z",row);
	REC__Cherenkov_dtheta = REC__Cherenkov->getFloat("dtheta",row);
	REC__Cherenkov_dphi = REC__Cherenkov->getFloat("dphi",row);
	REC__Cherenkov_status = REC__Cherenkov->getShort("status",row);
	return 0;
} 

int TIdentificatorCLAS12::get_BMTRec__LayerEffs(int row){
	BMTRec__LayerEffs_sector = BMTRec__LayerEffs->getByte("sector",row);
	BMTRec__LayerEffs_layer = BMTRec__LayerEffs->getByte("layer",row);
	BMTRec__LayerEffs_residual = BMTRec__LayerEffs->getFloat("residual",row);
	BMTRec__LayerEffs_status = BMTRec__LayerEffs->getByte("status",row);
	return 0;
} 

int TIdentificatorCLAS12::get_FTOF__adc(int row){
	FTOF__adc_sector = FTOF__adc->getByte("sector",row);
	FTOF__adc_layer = FTOF__adc->getByte("layer",row);
	FTOF__adc_component = FTOF__adc->getShort("component",row);
	FTOF__adc_order = FTOF__adc->getByte("order",row);
	FTOF__adc_ADC = FTOF__adc->getInt("ADC",row);
	FTOF__adc_time = FTOF__adc->getFloat("time",row);
	FTOF__adc_ped = FTOF__adc->getShort("ped",row);
	return 0;
} 

int TIdentificatorCLAS12::get_MC__Lund(int row){
	MC__Lund_index = MC__Lund->getByte("index",row);
	MC__Lund_lifetime = MC__Lund->getFloat("lifetime",row);
	MC__Lund_type = MC__Lund->getByte("type",row);
	MC__Lund_pid = MC__Lund->getInt("pid",row);
	MC__Lund_parent = MC__Lund->getByte("parent",row);
	MC__Lund_daughter = MC__Lund->getByte("daughter",row);
	MC__Lund_px = MC__Lund->getFloat("px",row);
	MC__Lund_py = MC__Lund->getFloat("py",row);
	MC__Lund_pz = MC__Lund->getFloat("pz",row);
	MC__Lund_energy = MC__Lund->getFloat("energy",row);
	MC__Lund_mass = MC__Lund->getFloat("mass",row);
	MC__Lund_vx = MC__Lund->getFloat("vx",row);
	MC__Lund_vy = MC__Lund->getFloat("vy",row);
	MC__Lund_vz = MC__Lund->getFloat("vz",row);
	return 0;
} 

int TIdentificatorCLAS12::get_DETECTOR__lcpb(int row){
	DETECTOR__lcpb_sector = DETECTOR__lcpb->getByte("sector",row);
	DETECTOR__lcpb_etot = DETECTOR__lcpb->getFloat("etot",row);
	DETECTOR__lcpb_ein = DETECTOR__lcpb->getFloat("ein",row);
	DETECTOR__lcpb_time = DETECTOR__lcpb->getFloat("time",row);
	DETECTOR__lcpb_path = DETECTOR__lcpb->getFloat("path",row);
	DETECTOR__lcpb_x = DETECTOR__lcpb->getFloat("x",row);
	DETECTOR__lcpb_y = DETECTOR__lcpb->getFloat("y",row);
	DETECTOR__lcpb_z = DETECTOR__lcpb->getFloat("z",row);
	return 0;
} 

int TIdentificatorCLAS12::get_MC__Header(int row){
	MC__Header_run = MC__Header->getInt("run",row);
	MC__Header_event = MC__Header->getInt("event",row);
	MC__Header_type = MC__Header->getByte("type",row);
	MC__Header_helicity = MC__Header->getFloat("helicity",row);
	return 0;
} 

int TIdentificatorCLAS12::get_CND__clusters(int row){
	CND__clusters_id = CND__clusters->getShort("id",row);
	CND__clusters_sector = CND__clusters->getByte("sector",row);
	CND__clusters_layer = CND__clusters->getByte("layer",row);
	CND__clusters_component = CND__clusters->getShort("component",row);
	CND__clusters_nhits = CND__clusters->getShort("nhits",row);
	CND__clusters_energy = CND__clusters->getFloat("energy",row);
	CND__clusters_x = CND__clusters->getFloat("x",row);
	CND__clusters_y = CND__clusters->getFloat("y",row);
	CND__clusters_z = CND__clusters->getFloat("z",row);
	CND__clusters_time = CND__clusters->getFloat("time",row);
	CND__clusters_status = CND__clusters->getShort("status",row);
	return 0;
} 

int TIdentificatorCLAS12::get_TimeBasedTrkg__TBCovMat(int row){
	TimeBasedTrkg__TBCovMat_id = TimeBasedTrkg__TBCovMat->getShort("id",row);
	TimeBasedTrkg__TBCovMat_C11 = TimeBasedTrkg__TBCovMat->getFloat("C11",row);
	TimeBasedTrkg__TBCovMat_C12 = TimeBasedTrkg__TBCovMat->getFloat("C12",row);
	TimeBasedTrkg__TBCovMat_C13 = TimeBasedTrkg__TBCovMat->getFloat("C13",row);
	TimeBasedTrkg__TBCovMat_C14 = TimeBasedTrkg__TBCovMat->getFloat("C14",row);
	TimeBasedTrkg__TBCovMat_C15 = TimeBasedTrkg__TBCovMat->getFloat("C15",row);
	TimeBasedTrkg__TBCovMat_C21 = TimeBasedTrkg__TBCovMat->getFloat("C21",row);
	TimeBasedTrkg__TBCovMat_C22 = TimeBasedTrkg__TBCovMat->getFloat("C22",row);
	TimeBasedTrkg__TBCovMat_C23 = TimeBasedTrkg__TBCovMat->getFloat("C23",row);
	TimeBasedTrkg__TBCovMat_C24 = TimeBasedTrkg__TBCovMat->getFloat("C24",row);
	TimeBasedTrkg__TBCovMat_C25 = TimeBasedTrkg__TBCovMat->getFloat("C25",row);
	TimeBasedTrkg__TBCovMat_C31 = TimeBasedTrkg__TBCovMat->getFloat("C31",row);
	TimeBasedTrkg__TBCovMat_C32 = TimeBasedTrkg__TBCovMat->getFloat("C32",row);
	TimeBasedTrkg__TBCovMat_C33 = TimeBasedTrkg__TBCovMat->getFloat("C33",row);
	TimeBasedTrkg__TBCovMat_C34 = TimeBasedTrkg__TBCovMat->getFloat("C34",row);
	TimeBasedTrkg__TBCovMat_C35 = TimeBasedTrkg__TBCovMat->getFloat("C35",row);
	TimeBasedTrkg__TBCovMat_C41 = TimeBasedTrkg__TBCovMat->getFloat("C41",row);
	TimeBasedTrkg__TBCovMat_C42 = TimeBasedTrkg__TBCovMat->getFloat("C42",row);
	TimeBasedTrkg__TBCovMat_C43 = TimeBasedTrkg__TBCovMat->getFloat("C43",row);
	TimeBasedTrkg__TBCovMat_C44 = TimeBasedTrkg__TBCovMat->getFloat("C44",row);
	TimeBasedTrkg__TBCovMat_C45 = TimeBasedTrkg__TBCovMat->getFloat("C45",row);
	TimeBasedTrkg__TBCovMat_C51 = TimeBasedTrkg__TBCovMat->getFloat("C51",row);
	TimeBasedTrkg__TBCovMat_C52 = TimeBasedTrkg__TBCovMat->getFloat("C52",row);
	TimeBasedTrkg__TBCovMat_C53 = TimeBasedTrkg__TBCovMat->getFloat("C53",row);
	TimeBasedTrkg__TBCovMat_C54 = TimeBasedTrkg__TBCovMat->getFloat("C54",row);
	TimeBasedTrkg__TBCovMat_C55 = TimeBasedTrkg__TBCovMat->getFloat("C55",row);
	return 0;
} 

int TIdentificatorCLAS12::get_RICH__hits(int row){
	RICH__hits_id = RICH__hits->getShort("id",row);
	RICH__hits_sector = RICH__hits->getShort("sector",row);
	RICH__hits_tile = RICH__hits->getShort("tile",row);
	RICH__hits_pmt = RICH__hits->getShort("pmt",row);
	RICH__hits_anode = RICH__hits->getShort("anode",row);
	RICH__hits_x = RICH__hits->getFloat("x",row);
	RICH__hits_y = RICH__hits->getFloat("y",row);
	RICH__hits_z = RICH__hits->getFloat("z",row);
	RICH__hits_time = RICH__hits->getFloat("time",row);
	RICH__hits_rawtime = RICH__hits->getFloat("rawtime",row);
	RICH__hits_cluster = RICH__hits->getShort("cluster",row);
	RICH__hits_xtalk = RICH__hits->getShort("xtalk",row);
	RICH__hits_duration = RICH__hits->getShort("duration",row);
	return 0;
} 

int TIdentificatorCLAS12::get_RECHB__Track(int row){
	RECHB__Track_index = RECHB__Track->getShort("index",row);
	RECHB__Track_pindex = RECHB__Track->getShort("pindex",row);
	RECHB__Track_detector = RECHB__Track->getByte("detector",row);
	RECHB__Track_sector = RECHB__Track->getByte("sector",row);
	RECHB__Track_status = RECHB__Track->getShort("status",row);
	RECHB__Track_q = RECHB__Track->getByte("q",row);
	RECHB__Track_chi2 = RECHB__Track->getFloat("chi2",row);
	RECHB__Track_NDF = RECHB__Track->getShort("NDF",row);
	return 0;
} 

int TIdentificatorCLAS12::get_MC__True(int row){
	MC__True_detector = MC__True->getByte("detector",row);
	MC__True_pid = MC__True->getInt("pid",row);
	MC__True_mpid = MC__True->getInt("mpid",row);
	MC__True_tid = MC__True->getInt("tid",row);
	MC__True_mtid = MC__True->getInt("mtid",row);
	MC__True_otid = MC__True->getInt("otid",row);
	MC__True_trackE = MC__True->getFloat("trackE",row);
	MC__True_totEdep = MC__True->getFloat("totEdep",row);
	MC__True_avgX = MC__True->getFloat("avgX",row);
	MC__True_avgY = MC__True->getFloat("avgY",row);
	MC__True_avgZ = MC__True->getFloat("avgZ",row);
	MC__True_avgLx = MC__True->getFloat("avgLx",row);
	MC__True_avgLy = MC__True->getFloat("avgLy",row);
	MC__True_avgLz = MC__True->getFloat("avgLz",row);
	MC__True_px = MC__True->getFloat("px",row);
	MC__True_py = MC__True->getFloat("py",row);
	MC__True_pz = MC__True->getFloat("pz",row);
	MC__True_vx = MC__True->getFloat("vx",row);
	MC__True_vy = MC__True->getFloat("vy",row);
	MC__True_vz = MC__True->getFloat("vz",row);
	MC__True_mvx = MC__True->getFloat("mvx",row);
	MC__True_mvy = MC__True->getFloat("mvy",row);
	MC__True_mvz = MC__True->getFloat("mvz",row);
	MC__True_avgT = MC__True->getFloat("avgT",row);
	MC__True_nsteps = MC__True->getInt("nsteps",row);
	MC__True_procID = MC__True->getInt("procID",row);
	MC__True_hitn = MC__True->getInt("hitn",row);
	return 0;
} 

int TIdentificatorCLAS12::get_BST__adc(int row){
	BST__adc_sector = BST__adc->getByte("sector",row);
	BST__adc_layer = BST__adc->getByte("layer",row);
	BST__adc_component = BST__adc->getShort("component",row);
	BST__adc_order = BST__adc->getByte("order",row);
	BST__adc_ADC = BST__adc->getInt("ADC",row);
	BST__adc_time = BST__adc->getFloat("time",row);
	BST__adc_ped = BST__adc->getShort("ped",row);
	BST__adc_timestamp = BST__adc->getLong("timestamp",row);
	return 0;
} 

int TIdentificatorCLAS12::get_MC__Event(int row){
	MC__Event_npart = MC__Event->getShort("npart",row);
	MC__Event_atarget = MC__Event->getShort("atarget",row);
	MC__Event_ztarget = MC__Event->getShort("ztarget",row);
	MC__Event_ptarget = MC__Event->getFloat("ptarget",row);
	MC__Event_pbeam = MC__Event->getFloat("pbeam",row);
	MC__Event_btype = MC__Event->getShort("btype",row);
	MC__Event_ebeam = MC__Event->getFloat("ebeam",row);
	MC__Event_targetid = MC__Event->getShort("targetid",row);
	MC__Event_processid = MC__Event->getShort("processid",row);
	MC__Event_weight = MC__Event->getFloat("weight",row);
	return 0;
} 

int TIdentificatorCLAS12::get_HitBasedTrkg__HBCrosses(int row){
	HitBasedTrkg__HBCrosses_id = HitBasedTrkg__HBCrosses->getShort("id",row);
	HitBasedTrkg__HBCrosses_status = HitBasedTrkg__HBCrosses->getShort("status",row);
	HitBasedTrkg__HBCrosses_sector = HitBasedTrkg__HBCrosses->getByte("sector",row);
	HitBasedTrkg__HBCrosses_region = HitBasedTrkg__HBCrosses->getByte("region",row);
	HitBasedTrkg__HBCrosses_x = HitBasedTrkg__HBCrosses->getFloat("x",row);
	HitBasedTrkg__HBCrosses_y = HitBasedTrkg__HBCrosses->getFloat("y",row);
	HitBasedTrkg__HBCrosses_z = HitBasedTrkg__HBCrosses->getFloat("z",row);
	HitBasedTrkg__HBCrosses_err_x = HitBasedTrkg__HBCrosses->getFloat("err_x",row);
	HitBasedTrkg__HBCrosses_err_y = HitBasedTrkg__HBCrosses->getFloat("err_y",row);
	HitBasedTrkg__HBCrosses_err_z = HitBasedTrkg__HBCrosses->getFloat("err_z",row);
	HitBasedTrkg__HBCrosses_ux = HitBasedTrkg__HBCrosses->getFloat("ux",row);
	HitBasedTrkg__HBCrosses_uy = HitBasedTrkg__HBCrosses->getFloat("uy",row);
	HitBasedTrkg__HBCrosses_uz = HitBasedTrkg__HBCrosses->getFloat("uz",row);
	HitBasedTrkg__HBCrosses_err_ux = HitBasedTrkg__HBCrosses->getFloat("err_ux",row);
	HitBasedTrkg__HBCrosses_err_uy = HitBasedTrkg__HBCrosses->getFloat("err_uy",row);
	HitBasedTrkg__HBCrosses_err_uz = HitBasedTrkg__HBCrosses->getFloat("err_uz",row);
	HitBasedTrkg__HBCrosses_Segment1_ID = HitBasedTrkg__HBCrosses->getShort("Segment1_ID",row);
	HitBasedTrkg__HBCrosses_Segment2_ID = HitBasedTrkg__HBCrosses->getShort("Segment2_ID",row);
	return 0;
} 

int TIdentificatorCLAS12::get_FTOF__clusters(int row){
	FTOF__clusters_id = FTOF__clusters->getShort("id",row);
	FTOF__clusters_status = FTOF__clusters->getShort("status",row);
	FTOF__clusters_trackid = FTOF__clusters->getShort("trackid",row);
	FTOF__clusters_sector = FTOF__clusters->getByte("sector",row);
	FTOF__clusters_layer = FTOF__clusters->getByte("layer",row);
	FTOF__clusters_component = FTOF__clusters->getShort("component",row);
	FTOF__clusters_energy = FTOF__clusters->getFloat("energy",row);
	FTOF__clusters_time = FTOF__clusters->getFloat("time",row);
	FTOF__clusters_energy_unc = FTOF__clusters->getFloat("energy_unc",row);
	FTOF__clusters_time_unc = FTOF__clusters->getFloat("time_unc",row);
	FTOF__clusters_x = FTOF__clusters->getFloat("x",row);
	FTOF__clusters_y = FTOF__clusters->getFloat("y",row);
	FTOF__clusters_z = FTOF__clusters->getFloat("z",row);
	FTOF__clusters_x_unc = FTOF__clusters->getFloat("x_unc",row);
	FTOF__clusters_y_unc = FTOF__clusters->getFloat("y_unc",row);
	FTOF__clusters_z_unc = FTOF__clusters->getFloat("z_unc",row);
	return 0;
} 

int TIdentificatorCLAS12::get_REC__TrackCross(int row){
	REC__TrackCross_index = REC__TrackCross->getShort("index",row);
	REC__TrackCross_pindex = REC__TrackCross->getShort("pindex",row);
	REC__TrackCross_detector = REC__TrackCross->getByte("detector",row);
	REC__TrackCross_sector = REC__TrackCross->getByte("sector",row);
	REC__TrackCross_layer = REC__TrackCross->getByte("layer",row);
	REC__TrackCross_c_x = REC__TrackCross->getFloat("c_x",row);
	REC__TrackCross_c_y = REC__TrackCross->getFloat("c_y",row);
	REC__TrackCross_c_z = REC__TrackCross->getFloat("c_z",row);
	REC__TrackCross_c_ux = REC__TrackCross->getFloat("c_ux",row);
	REC__TrackCross_c_uy = REC__TrackCross->getFloat("c_uy",row);
	REC__TrackCross_c_uz = REC__TrackCross->getFloat("c_uz",row);
	REC__TrackCross_status = REC__TrackCross->getShort("status",row);
	return 0;
} 

int TIdentificatorCLAS12::get_REC__Scintillator(int row){
	REC__Scintillator_index = REC__Scintillator->getShort("index",row);
	REC__Scintillator_pindex = REC__Scintillator->getShort("pindex",row);
	REC__Scintillator_detector = REC__Scintillator->getByte("detector",row);
	REC__Scintillator_sector = REC__Scintillator->getByte("sector",row);
	REC__Scintillator_layer = REC__Scintillator->getByte("layer",row);
	REC__Scintillator_component = REC__Scintillator->getShort("component",row);
	REC__Scintillator_energy = REC__Scintillator->getFloat("energy",row);
	REC__Scintillator_time = REC__Scintillator->getFloat("time",row);
	REC__Scintillator_path = REC__Scintillator->getFloat("path",row);
	REC__Scintillator_chi2 = REC__Scintillator->getFloat("chi2",row);
	REC__Scintillator_x = REC__Scintillator->getFloat("x",row);
	REC__Scintillator_y = REC__Scintillator->getFloat("y",row);
	REC__Scintillator_z = REC__Scintillator->getFloat("z",row);
	REC__Scintillator_hx = REC__Scintillator->getFloat("hx",row);
	REC__Scintillator_hy = REC__Scintillator->getFloat("hy",row);
	REC__Scintillator_hz = REC__Scintillator->getFloat("hz",row);
	REC__Scintillator_status = REC__Scintillator->getShort("status",row);
	return 0;
} 

int TIdentificatorCLAS12::get_ECAL__peaks(int row){
	ECAL__peaks_id = ECAL__peaks->getShort("id",row);
	ECAL__peaks_status = ECAL__peaks->getShort("status",row);
	ECAL__peaks_sector = ECAL__peaks->getByte("sector",row);
	ECAL__peaks_layer = ECAL__peaks->getByte("layer",row);
	ECAL__peaks_energy = ECAL__peaks->getFloat("energy",row);
	ECAL__peaks_time = ECAL__peaks->getFloat("time",row);
	ECAL__peaks_xo = ECAL__peaks->getFloat("xo",row);
	ECAL__peaks_yo = ECAL__peaks->getFloat("yo",row);
	ECAL__peaks_zo = ECAL__peaks->getFloat("zo",row);
	ECAL__peaks_xe = ECAL__peaks->getFloat("xe",row);
	ECAL__peaks_ye = ECAL__peaks->getFloat("ye",row);
	ECAL__peaks_ze = ECAL__peaks->getFloat("ze",row);
	ECAL__peaks_width = ECAL__peaks->getFloat("width",row);
	return 0;
} 

int TIdentificatorCLAS12::get_TimeBasedTrkg__TBClusters(int row){
	TimeBasedTrkg__TBClusters_id = TimeBasedTrkg__TBClusters->getShort("id",row);
	TimeBasedTrkg__TBClusters_status = TimeBasedTrkg__TBClusters->getShort("status",row);
	TimeBasedTrkg__TBClusters_sector = TimeBasedTrkg__TBClusters->getByte("sector",row);
	TimeBasedTrkg__TBClusters_superlayer = TimeBasedTrkg__TBClusters->getByte("superlayer",row);
	TimeBasedTrkg__TBClusters_Hit1_ID = TimeBasedTrkg__TBClusters->getShort("Hit1_ID",row);
	TimeBasedTrkg__TBClusters_Hit2_ID = TimeBasedTrkg__TBClusters->getShort("Hit2_ID",row);
	TimeBasedTrkg__TBClusters_Hit3_ID = TimeBasedTrkg__TBClusters->getShort("Hit3_ID",row);
	TimeBasedTrkg__TBClusters_Hit4_ID = TimeBasedTrkg__TBClusters->getShort("Hit4_ID",row);
	TimeBasedTrkg__TBClusters_Hit5_ID = TimeBasedTrkg__TBClusters->getShort("Hit5_ID",row);
	TimeBasedTrkg__TBClusters_Hit6_ID = TimeBasedTrkg__TBClusters->getShort("Hit6_ID",row);
	TimeBasedTrkg__TBClusters_Hit7_ID = TimeBasedTrkg__TBClusters->getShort("Hit7_ID",row);
	TimeBasedTrkg__TBClusters_Hit8_ID = TimeBasedTrkg__TBClusters->getShort("Hit8_ID",row);
	TimeBasedTrkg__TBClusters_Hit9_ID = TimeBasedTrkg__TBClusters->getShort("Hit9_ID",row);
	TimeBasedTrkg__TBClusters_Hit10_ID = TimeBasedTrkg__TBClusters->getShort("Hit10_ID",row);
	TimeBasedTrkg__TBClusters_Hit11_ID = TimeBasedTrkg__TBClusters->getShort("Hit11_ID",row);
	TimeBasedTrkg__TBClusters_Hit12_ID = TimeBasedTrkg__TBClusters->getShort("Hit12_ID",row);
	TimeBasedTrkg__TBClusters_avgWire = TimeBasedTrkg__TBClusters->getFloat("avgWire",row);
	TimeBasedTrkg__TBClusters_fitChisqProb = TimeBasedTrkg__TBClusters->getFloat("fitChisqProb",row);
	TimeBasedTrkg__TBClusters_fitSlope = TimeBasedTrkg__TBClusters->getFloat("fitSlope",row);
	TimeBasedTrkg__TBClusters_fitSlopeErr = TimeBasedTrkg__TBClusters->getFloat("fitSlopeErr",row);
	TimeBasedTrkg__TBClusters_fitInterc = TimeBasedTrkg__TBClusters->getFloat("fitInterc",row);
	TimeBasedTrkg__TBClusters_fitIntercErr = TimeBasedTrkg__TBClusters->getFloat("fitIntercErr",row);
	TimeBasedTrkg__TBClusters_size = TimeBasedTrkg__TBClusters->getByte("size",row);
	return 0;
} 

int TIdentificatorCLAS12::get_TAGGER__tgpb(int row){
	TAGGER__tgpb_status = TAGGER__tgpb->getByte("status",row);
	TAGGER__tgpb_tid = TAGGER__tgpb->getShort("tid",row);
	TAGGER__tgpb_eid = TAGGER__tgpb->getShort("eid",row);
	TAGGER__tgpb_time = TAGGER__tgpb->getFloat("time",row);
	TAGGER__tgpb_energy = TAGGER__tgpb->getFloat("energy",row);
	return 0;
} 

int TIdentificatorCLAS12::get_LTCC__clusters(int row){
	LTCC__clusters_id = LTCC__clusters->getShort("id",row);
	LTCC__clusters_status = LTCC__clusters->getByte("status",row);
	LTCC__clusters_sector = LTCC__clusters->getByte("sector",row);
	LTCC__clusters_segment = LTCC__clusters->getShort("segment",row);
	LTCC__clusters_x = LTCC__clusters->getFloat("x",row);
	LTCC__clusters_y = LTCC__clusters->getFloat("y",row);
	LTCC__clusters_z = LTCC__clusters->getFloat("z",row);
	LTCC__clusters_nphe = LTCC__clusters->getFloat("nphe",row);
	LTCC__clusters_time = LTCC__clusters->getFloat("time",row);
	LTCC__clusters_nHits = LTCC__clusters->getShort("nHits",row);
	LTCC__clusters_minTheta = LTCC__clusters->getFloat("minTheta",row);
	LTCC__clusters_maxTheta = LTCC__clusters->getFloat("maxTheta",row);
	LTCC__clusters_minPhi = LTCC__clusters->getFloat("minPhi",row);
	LTCC__clusters_maxPhi = LTCC__clusters->getFloat("maxPhi",row);
	return 0;
} 

int TIdentificatorCLAS12::get_RICH__hadrons(int row){
	RICH__hadrons_id = RICH__hadrons->getShort("id",row);
	RICH__hadrons_hit_index = RICH__hadrons->getShort("hit_index",row);
	RICH__hadrons_particle_index = RICH__hadrons->getShort("particle_index",row);
	RICH__hadrons_traced_the = RICH__hadrons->getFloat("traced_the",row);
	RICH__hadrons_traced_phi = RICH__hadrons->getFloat("traced_phi",row);
	RICH__hadrons_traced_hitx = RICH__hadrons->getFloat("traced_hitx",row);
	RICH__hadrons_traced_hity = RICH__hadrons->getFloat("traced_hity",row);
	RICH__hadrons_traced_hitz = RICH__hadrons->getFloat("traced_hitz",row);
	RICH__hadrons_traced_time = RICH__hadrons->getFloat("traced_time",row);
	RICH__hadrons_traced_path = RICH__hadrons->getFloat("traced_path",row);
	RICH__hadrons_traced_ilay = RICH__hadrons->getShort("traced_ilay",row);
	RICH__hadrons_traced_ico = RICH__hadrons->getShort("traced_ico",row);
	RICH__hadrons_traced_emix = RICH__hadrons->getFloat("traced_emix",row);
	RICH__hadrons_traced_emiy = RICH__hadrons->getFloat("traced_emiy",row);
	RICH__hadrons_traced_emiz = RICH__hadrons->getFloat("traced_emiz",row);
	RICH__hadrons_EtaC_ele = RICH__hadrons->getFloat("EtaC_ele",row);
	RICH__hadrons_EtaC_pi = RICH__hadrons->getFloat("EtaC_pi",row);
	RICH__hadrons_EtaC_k = RICH__hadrons->getFloat("EtaC_k",row);
	RICH__hadrons_EtaC_pr = RICH__hadrons->getFloat("EtaC_pr",row);
	RICH__hadrons_EtaC_min = RICH__hadrons->getFloat("EtaC_min",row);
	RICH__hadrons_EtaC_max = RICH__hadrons->getFloat("EtaC_max",row);
	return 0;
} 

int TIdentificatorCLAS12::get_DC__tdc(int row){
	DC__tdc_sector = DC__tdc->getByte("sector",row);
	DC__tdc_layer = DC__tdc->getByte("layer",row);
	DC__tdc_component = DC__tdc->getShort("component",row);
	DC__tdc_order = DC__tdc->getByte("order",row);
	DC__tdc_TDC = DC__tdc->getInt("TDC",row);
	return 0;
} 

int TIdentificatorCLAS12::get_BSTRec__Crosses(int row){
	BSTRec__Crosses_ID = BSTRec__Crosses->getShort("ID",row);
	BSTRec__Crosses_sector = BSTRec__Crosses->getByte("sector",row);
	BSTRec__Crosses_region = BSTRec__Crosses->getByte("region",row);
	BSTRec__Crosses_x = BSTRec__Crosses->getFloat("x",row);
	BSTRec__Crosses_y = BSTRec__Crosses->getFloat("y",row);
	BSTRec__Crosses_z = BSTRec__Crosses->getFloat("z",row);
	BSTRec__Crosses_err_x = BSTRec__Crosses->getFloat("err_x",row);
	BSTRec__Crosses_err_y = BSTRec__Crosses->getFloat("err_y",row);
	BSTRec__Crosses_err_z = BSTRec__Crosses->getFloat("err_z",row);
	BSTRec__Crosses_ux = BSTRec__Crosses->getFloat("ux",row);
	BSTRec__Crosses_uy = BSTRec__Crosses->getFloat("uy",row);
	BSTRec__Crosses_uz = BSTRec__Crosses->getFloat("uz",row);
	BSTRec__Crosses_Cluster1_ID = BSTRec__Crosses->getShort("Cluster1_ID",row);
	BSTRec__Crosses_Cluster2_ID = BSTRec__Crosses->getShort("Cluster2_ID",row);
	BSTRec__Crosses_trkID = BSTRec__Crosses->getShort("trkID",row);
	return 0;
} 

int TIdentificatorCLAS12::get_HTCC__rec(int row){
	HTCC__rec_id = HTCC__rec->getShort("id",row);
	HTCC__rec_nhits = HTCC__rec->getShort("nhits",row);
	HTCC__rec_nphe = HTCC__rec->getFloat("nphe",row);
	HTCC__rec_ntheta = HTCC__rec->getShort("ntheta",row);
	HTCC__rec_nphi = HTCC__rec->getShort("nphi",row);
	HTCC__rec_mintheta = HTCC__rec->getShort("mintheta",row);
	HTCC__rec_maxtheta = HTCC__rec->getShort("maxtheta",row);
	HTCC__rec_minphi = HTCC__rec->getShort("minphi",row);
	HTCC__rec_maxphi = HTCC__rec->getShort("maxphi",row);
	HTCC__rec_time = HTCC__rec->getFloat("time",row);
	HTCC__rec_theta = HTCC__rec->getFloat("theta",row);
	HTCC__rec_dtheta = HTCC__rec->getFloat("dtheta",row);
	HTCC__rec_phi = HTCC__rec->getFloat("phi",row);
	HTCC__rec_dphi = HTCC__rec->getFloat("dphi",row);
	HTCC__rec_x = HTCC__rec->getFloat("x",row);
	HTCC__rec_y = HTCC__rec->getFloat("y",row);
	HTCC__rec_z = HTCC__rec->getFloat("z",row);
	return 0;
} 

int TIdentificatorCLAS12::get_FTHODO__adc(int row){
	FTHODO__adc_sector = FTHODO__adc->getByte("sector",row);
	FTHODO__adc_layer = FTHODO__adc->getByte("layer",row);
	FTHODO__adc_component = FTHODO__adc->getShort("component",row);
	FTHODO__adc_order = FTHODO__adc->getByte("order",row);
	FTHODO__adc_ADC = FTHODO__adc->getInt("ADC",row);
	FTHODO__adc_time = FTHODO__adc->getFloat("time",row);
	FTHODO__adc_ped = FTHODO__adc->getShort("ped",row);
	return 0;
} 

int TIdentificatorCLAS12::get_FTOF__matchedclusters(int row){
	FTOF__matchedclusters_sector = FTOF__matchedclusters->getByte("sector",row);
	FTOF__matchedclusters_paddle_id1A = FTOF__matchedclusters->getShort("paddle_id1A",row);
	FTOF__matchedclusters_paddle_id1B = FTOF__matchedclusters->getShort("paddle_id1B",row);
	FTOF__matchedclusters_clus_1Aid = FTOF__matchedclusters->getShort("clus_1Aid",row);
	FTOF__matchedclusters_clus_1Bid = FTOF__matchedclusters->getShort("clus_1Bid",row);
	FTOF__matchedclusters_clusSize_1A = FTOF__matchedclusters->getShort("clusSize_1A",row);
	FTOF__matchedclusters_clusSize_1B = FTOF__matchedclusters->getShort("clusSize_1B",row);
	FTOF__matchedclusters_tminAlgo_1B_tCorr = FTOF__matchedclusters->getFloat("tminAlgo_1B_tCorr",row);
	FTOF__matchedclusters_midbarAlgo_1B_tCorr = FTOF__matchedclusters->getFloat("midbarAlgo_1B_tCorr",row);
	FTOF__matchedclusters_EmaxAlgo_1B_tCorr = FTOF__matchedclusters->getFloat("EmaxAlgo_1B_tCorr",row);
	return 0;
} 

int TIdentificatorCLAS12::get_ECAL__adc(int row){
	ECAL__adc_sector = ECAL__adc->getByte("sector",row);
	ECAL__adc_layer = ECAL__adc->getByte("layer",row);
	ECAL__adc_component = ECAL__adc->getShort("component",row);
	ECAL__adc_order = ECAL__adc->getByte("order",row);
	ECAL__adc_ADC = ECAL__adc->getInt("ADC",row);
	ECAL__adc_time = ECAL__adc->getFloat("time",row);
	ECAL__adc_ped = ECAL__adc->getShort("ped",row);
	return 0;
} 

int TIdentificatorCLAS12::get_FMT__adc(int row){
	FMT__adc_sector = FMT__adc->getByte("sector",row);
	FMT__adc_layer = FMT__adc->getByte("layer",row);
	FMT__adc_component = FMT__adc->getShort("component",row);
	FMT__adc_order = FMT__adc->getByte("order",row);
	FMT__adc_ADC = FMT__adc->getInt("ADC",row);
	FMT__adc_time = FMT__adc->getFloat("time",row);
	FMT__adc_ped = FMT__adc->getShort("ped",row);
	FMT__adc_integral = FMT__adc->getInt("integral",row);
	FMT__adc_timestamp = FMT__adc->getLong("timestamp",row);
	return 0;
} 

int TIdentificatorCLAS12::get_BMT__adc(int row){
	BMT__adc_sector = BMT__adc->getByte("sector",row);
	BMT__adc_layer = BMT__adc->getByte("layer",row);
	BMT__adc_component = BMT__adc->getShort("component",row);
	BMT__adc_order = BMT__adc->getByte("order",row);
	BMT__adc_ADC = BMT__adc->getInt("ADC",row);
	BMT__adc_time = BMT__adc->getFloat("time",row);
	BMT__adc_ped = BMT__adc->getShort("ped",row);
	BMT__adc_integral = BMT__adc->getInt("integral",row);
	BMT__adc_timestamp = BMT__adc->getLong("timestamp",row);
	return 0;
} 

int TIdentificatorCLAS12::get_ECAL__calib(int row){
	ECAL__calib_sector = ECAL__calib->getByte("sector",row);
	ECAL__calib_layer = ECAL__calib->getByte("layer",row);
	ECAL__calib_energy = ECAL__calib->getFloat("energy",row);
	ECAL__calib_rawEU = ECAL__calib->getFloat("rawEU",row);
	ECAL__calib_rawEV = ECAL__calib->getFloat("rawEV",row);
	ECAL__calib_rawEW = ECAL__calib->getFloat("rawEW",row);
	ECAL__calib_recEU = ECAL__calib->getFloat("recEU",row);
	ECAL__calib_recEV = ECAL__calib->getFloat("recEV",row);
	ECAL__calib_recEW = ECAL__calib->getFloat("recEW",row);
	return 0;
} 

int TIdentificatorCLAS12::get_BMTRec__Crosses(int row){
	BMTRec__Crosses_ID = BMTRec__Crosses->getShort("ID",row);
	BMTRec__Crosses_sector = BMTRec__Crosses->getByte("sector",row);
	BMTRec__Crosses_region = BMTRec__Crosses->getByte("region",row);
	BMTRec__Crosses_x = BMTRec__Crosses->getFloat("x",row);
	BMTRec__Crosses_y = BMTRec__Crosses->getFloat("y",row);
	BMTRec__Crosses_z = BMTRec__Crosses->getFloat("z",row);
	BMTRec__Crosses_err_x = BMTRec__Crosses->getFloat("err_x",row);
	BMTRec__Crosses_err_y = BMTRec__Crosses->getFloat("err_y",row);
	BMTRec__Crosses_err_z = BMTRec__Crosses->getFloat("err_z",row);
	BMTRec__Crosses_ux = BMTRec__Crosses->getFloat("ux",row);
	BMTRec__Crosses_uy = BMTRec__Crosses->getFloat("uy",row);
	BMTRec__Crosses_uz = BMTRec__Crosses->getFloat("uz",row);
	BMTRec__Crosses_Cluster1_ID = BMTRec__Crosses->getShort("Cluster1_ID",row);
	BMTRec__Crosses_Cluster2_ID = BMTRec__Crosses->getShort("Cluster2_ID",row);
	BMTRec__Crosses_trkID = BMTRec__Crosses->getShort("trkID",row);
	return 0;
} 

int TIdentificatorCLAS12::get_RUN__rf(int row){
	RUN__rf_id = RUN__rf->getShort("id",row);
	RUN__rf_time = RUN__rf->getFloat("time",row);
	return 0;
} 

int TIdentificatorCLAS12::get_BSTRec__Hits(int row){
	BSTRec__Hits_ID = BSTRec__Hits->getShort("ID",row);
	BSTRec__Hits_sector = BSTRec__Hits->getByte("sector",row);
	BSTRec__Hits_layer = BSTRec__Hits->getByte("layer",row);
	BSTRec__Hits_strip = BSTRec__Hits->getInt("strip",row);
	BSTRec__Hits_fitResidual = BSTRec__Hits->getFloat("fitResidual",row);
	BSTRec__Hits_trkingStat = BSTRec__Hits->getInt("trkingStat",row);
	BSTRec__Hits_clusterID = BSTRec__Hits->getShort("clusterID",row);
	BSTRec__Hits_trkID = BSTRec__Hits->getShort("trkID",row);
	return 0;
} 

int TIdentificatorCLAS12::get_HitBasedTrkg__HBSegmentTrajectory(int row){
	HitBasedTrkg__HBSegmentTrajectory_segmentID = HitBasedTrkg__HBSegmentTrajectory->getShort("segmentID",row);
	HitBasedTrkg__HBSegmentTrajectory_sector = HitBasedTrkg__HBSegmentTrajectory->getByte("sector",row);
	HitBasedTrkg__HBSegmentTrajectory_superlayer = HitBasedTrkg__HBSegmentTrajectory->getByte("superlayer",row);
	HitBasedTrkg__HBSegmentTrajectory_layer = HitBasedTrkg__HBSegmentTrajectory->getByte("layer",row);
	HitBasedTrkg__HBSegmentTrajectory_matchedHitID = HitBasedTrkg__HBSegmentTrajectory->getShort("matchedHitID",row);
	HitBasedTrkg__HBSegmentTrajectory_trkDoca = HitBasedTrkg__HBSegmentTrajectory->getFloat("trkDoca",row);
	return 0;
} 

int TIdentificatorCLAS12::get_LTCC__adc(int row){
	LTCC__adc_sector = LTCC__adc->getByte("sector",row);
	LTCC__adc_layer = LTCC__adc->getByte("layer",row);
	LTCC__adc_component = LTCC__adc->getShort("component",row);
	LTCC__adc_order = LTCC__adc->getByte("order",row);
	LTCC__adc_ADC = LTCC__adc->getInt("ADC",row);
	LTCC__adc_time = LTCC__adc->getFloat("time",row);
	LTCC__adc_ped = LTCC__adc->getShort("ped",row);
	return 0;
} 

int TIdentificatorCLAS12::get_CND__tdc(int row){
	CND__tdc_sector = CND__tdc->getByte("sector",row);
	CND__tdc_layer = CND__tdc->getByte("layer",row);
	CND__tdc_component = CND__tdc->getShort("component",row);
	CND__tdc_order = CND__tdc->getByte("order",row);
	CND__tdc_TDC = CND__tdc->getInt("TDC",row);
	return 0;
} 

int TIdentificatorCLAS12::get_HTCC__adc(int row){
	HTCC__adc_sector = HTCC__adc->getByte("sector",row);
	HTCC__adc_layer = HTCC__adc->getByte("layer",row);
	HTCC__adc_component = HTCC__adc->getShort("component",row);
	HTCC__adc_order = HTCC__adc->getByte("order",row);
	HTCC__adc_ADC = HTCC__adc->getInt("ADC",row);
	HTCC__adc_time = HTCC__adc->getFloat("time",row);
	HTCC__adc_ped = HTCC__adc->getShort("ped",row);
	return 0;
} 

int TIdentificatorCLAS12::get_RAW__vtp(int row){
	RAW__vtp_crate = RAW__vtp->getByte("crate",row);
	RAW__vtp_word = RAW__vtp->getInt("word",row);
	return 0;
} 

int TIdentificatorCLAS12::get_CVTRec__Tracks(int row){
	CVTRec__Tracks_ID = CVTRec__Tracks->getShort("ID",row);
	CVTRec__Tracks_fittingMethod = CVTRec__Tracks->getByte("fittingMethod",row);
	CVTRec__Tracks_c_x = CVTRec__Tracks->getFloat("c_x",row);
	CVTRec__Tracks_c_y = CVTRec__Tracks->getFloat("c_y",row);
	CVTRec__Tracks_c_z = CVTRec__Tracks->getFloat("c_z",row);
	CVTRec__Tracks_c_ux = CVTRec__Tracks->getFloat("c_ux",row);
	CVTRec__Tracks_c_uy = CVTRec__Tracks->getFloat("c_uy",row);
	CVTRec__Tracks_c_uz = CVTRec__Tracks->getFloat("c_uz",row);
	CVTRec__Tracks_pathlength = CVTRec__Tracks->getFloat("pathlength",row);
	CVTRec__Tracks_q = CVTRec__Tracks->getByte("q",row);
	CVTRec__Tracks_p = CVTRec__Tracks->getFloat("p",row);
	CVTRec__Tracks_pt = CVTRec__Tracks->getFloat("pt",row);
	CVTRec__Tracks_phi0 = CVTRec__Tracks->getFloat("phi0",row);
	CVTRec__Tracks_tandip = CVTRec__Tracks->getFloat("tandip",row);
	CVTRec__Tracks_z0 = CVTRec__Tracks->getFloat("z0",row);
	CVTRec__Tracks_d0 = CVTRec__Tracks->getFloat("d0",row);
	CVTRec__Tracks_cov_d02 = CVTRec__Tracks->getFloat("cov_d02",row);
	CVTRec__Tracks_cov_d0phi0 = CVTRec__Tracks->getFloat("cov_d0phi0",row);
	CVTRec__Tracks_cov_d0rho = CVTRec__Tracks->getFloat("cov_d0rho",row);
	CVTRec__Tracks_cov_phi02 = CVTRec__Tracks->getFloat("cov_phi02",row);
	CVTRec__Tracks_cov_phi0rho = CVTRec__Tracks->getFloat("cov_phi0rho",row);
	CVTRec__Tracks_cov_rho2 = CVTRec__Tracks->getFloat("cov_rho2",row);
	CVTRec__Tracks_cov_z02 = CVTRec__Tracks->getFloat("cov_z02",row);
	CVTRec__Tracks_cov_tandip2 = CVTRec__Tracks->getFloat("cov_tandip2",row);
	CVTRec__Tracks_circlefit_chi2_per_ndf = CVTRec__Tracks->getFloat("circlefit_chi2_per_ndf",row);
	CVTRec__Tracks_linefit_chi2_per_ndf = CVTRec__Tracks->getFloat("linefit_chi2_per_ndf",row);
	CVTRec__Tracks_chi2 = CVTRec__Tracks->getFloat("chi2",row);
	CVTRec__Tracks_ndf = CVTRec__Tracks->getShort("ndf",row);
	CVTRec__Tracks_Cross1_ID = CVTRec__Tracks->getShort("Cross1_ID",row);
	CVTRec__Tracks_Cross2_ID = CVTRec__Tracks->getShort("Cross2_ID",row);
	CVTRec__Tracks_Cross3_ID = CVTRec__Tracks->getShort("Cross3_ID",row);
	CVTRec__Tracks_Cross4_ID = CVTRec__Tracks->getShort("Cross4_ID",row);
	CVTRec__Tracks_Cross5_ID = CVTRec__Tracks->getShort("Cross5_ID",row);
	CVTRec__Tracks_Cross6_ID = CVTRec__Tracks->getShort("Cross6_ID",row);
	CVTRec__Tracks_Cross7_ID = CVTRec__Tracks->getShort("Cross7_ID",row);
	CVTRec__Tracks_Cross8_ID = CVTRec__Tracks->getShort("Cross8_ID",row);
	CVTRec__Tracks_Cross9_ID = CVTRec__Tracks->getShort("Cross9_ID",row);
	return 0;
} 

int TIdentificatorCLAS12::get_RECFT__Particle(int row){
	RECFT__Particle_pid = RECFT__Particle->getInt("pid",row);
	RECFT__Particle_beta = RECFT__Particle->getFloat("beta",row);
	RECFT__Particle_chi2pid = RECFT__Particle->getFloat("chi2pid",row);
	RECFT__Particle_status = RECFT__Particle->getShort("status",row);
	return 0;
} 

int TIdentificatorCLAS12::get_TimeBasedTrkg__TBTracks(int row){
	TimeBasedTrkg__TBTracks_id = TimeBasedTrkg__TBTracks->getShort("id",row);
	TimeBasedTrkg__TBTracks_status = TimeBasedTrkg__TBTracks->getShort("status",row);
	TimeBasedTrkg__TBTracks_sector = TimeBasedTrkg__TBTracks->getByte("sector",row);
	TimeBasedTrkg__TBTracks_c1_x = TimeBasedTrkg__TBTracks->getFloat("c1_x",row);
	TimeBasedTrkg__TBTracks_c1_y = TimeBasedTrkg__TBTracks->getFloat("c1_y",row);
	TimeBasedTrkg__TBTracks_c1_z = TimeBasedTrkg__TBTracks->getFloat("c1_z",row);
	TimeBasedTrkg__TBTracks_c1_ux = TimeBasedTrkg__TBTracks->getFloat("c1_ux",row);
	TimeBasedTrkg__TBTracks_c1_uy = TimeBasedTrkg__TBTracks->getFloat("c1_uy",row);
	TimeBasedTrkg__TBTracks_c1_uz = TimeBasedTrkg__TBTracks->getFloat("c1_uz",row);
	TimeBasedTrkg__TBTracks_c3_x = TimeBasedTrkg__TBTracks->getFloat("c3_x",row);
	TimeBasedTrkg__TBTracks_c3_y = TimeBasedTrkg__TBTracks->getFloat("c3_y",row);
	TimeBasedTrkg__TBTracks_c3_z = TimeBasedTrkg__TBTracks->getFloat("c3_z",row);
	TimeBasedTrkg__TBTracks_c3_ux = TimeBasedTrkg__TBTracks->getFloat("c3_ux",row);
	TimeBasedTrkg__TBTracks_c3_uy = TimeBasedTrkg__TBTracks->getFloat("c3_uy",row);
	TimeBasedTrkg__TBTracks_c3_uz = TimeBasedTrkg__TBTracks->getFloat("c3_uz",row);
	TimeBasedTrkg__TBTracks_t1_x = TimeBasedTrkg__TBTracks->getFloat("t1_x",row);
	TimeBasedTrkg__TBTracks_t1_y = TimeBasedTrkg__TBTracks->getFloat("t1_y",row);
	TimeBasedTrkg__TBTracks_t1_z = TimeBasedTrkg__TBTracks->getFloat("t1_z",row);
	TimeBasedTrkg__TBTracks_t1_px = TimeBasedTrkg__TBTracks->getFloat("t1_px",row);
	TimeBasedTrkg__TBTracks_t1_py = TimeBasedTrkg__TBTracks->getFloat("t1_py",row);
	TimeBasedTrkg__TBTracks_t1_pz = TimeBasedTrkg__TBTracks->getFloat("t1_pz",row);
	TimeBasedTrkg__TBTracks_Vtx0_x = TimeBasedTrkg__TBTracks->getFloat("Vtx0_x",row);
	TimeBasedTrkg__TBTracks_Vtx0_y = TimeBasedTrkg__TBTracks->getFloat("Vtx0_y",row);
	TimeBasedTrkg__TBTracks_Vtx0_z = TimeBasedTrkg__TBTracks->getFloat("Vtx0_z",row);
	TimeBasedTrkg__TBTracks_p0_x = TimeBasedTrkg__TBTracks->getFloat("p0_x",row);
	TimeBasedTrkg__TBTracks_p0_y = TimeBasedTrkg__TBTracks->getFloat("p0_y",row);
	TimeBasedTrkg__TBTracks_p0_z = TimeBasedTrkg__TBTracks->getFloat("p0_z",row);
	TimeBasedTrkg__TBTracks_Cross1_ID = TimeBasedTrkg__TBTracks->getShort("Cross1_ID",row);
	TimeBasedTrkg__TBTracks_Cross2_ID = TimeBasedTrkg__TBTracks->getShort("Cross2_ID",row);
	TimeBasedTrkg__TBTracks_Cross3_ID = TimeBasedTrkg__TBTracks->getShort("Cross3_ID",row);
	TimeBasedTrkg__TBTracks_q = TimeBasedTrkg__TBTracks->getByte("q",row);
	TimeBasedTrkg__TBTracks_pathlength = TimeBasedTrkg__TBTracks->getFloat("pathlength",row);
	TimeBasedTrkg__TBTracks_chi2 = TimeBasedTrkg__TBTracks->getFloat("chi2",row);
	TimeBasedTrkg__TBTracks_ndf = TimeBasedTrkg__TBTracks->getShort("ndf",row);
	return 0;
} 

int TIdentificatorCLAS12::get_DETECTOR__icpb(int row){
	DETECTOR__icpb_etc = DETECTOR__icpb->getFloat("etc",row);
	DETECTOR__icpb_ecc = DETECTOR__icpb->getFloat("ecc",row);
	DETECTOR__icpb_tc = DETECTOR__icpb->getFloat("tc",row);
	DETECTOR__icpb_xc = DETECTOR__icpb->getFloat("xc",row);
	DETECTOR__icpb_yc = DETECTOR__icpb->getFloat("yc",row);
	return 0;
} 

int TIdentificatorCLAS12::get_RF__tdc(int row){
	RF__tdc_sector = RF__tdc->getByte("sector",row);
	RF__tdc_layer = RF__tdc->getByte("layer",row);
	RF__tdc_component = RF__tdc->getShort("component",row);
	RF__tdc_order = RF__tdc->getByte("order",row);
	RF__tdc_TDC = RF__tdc->getInt("TDC",row);
	return 0;
} 

int TIdentificatorCLAS12::get_HEADER__info(int row){
	HEADER__info_nrun = HEADER__info->getInt("nrun",row);
	HEADER__info_nevt = HEADER__info->getInt("nevt",row);
	HEADER__info_trigger = HEADER__info->getInt("trigger",row);
	HEADER__info_helicity = HEADER__info->getByte("helicity",row);
	HEADER__info_fc = HEADER__info->getFloat("fc",row);
	HEADER__info_fcg = HEADER__info->getFloat("fcg",row);
	HEADER__info_stt = HEADER__info->getFloat("stt",row);
	HEADER__info_rastr1 = HEADER__info->getShort("rastr1",row);
	HEADER__info_rastr2 = HEADER__info->getShort("rastr2",row);
	return 0;
} 

int TIdentificatorCLAS12::get_FTCAL__clusters(int row){
	FTCAL__clusters_size = FTCAL__clusters->getShort("size",row);
	FTCAL__clusters_id = FTCAL__clusters->getShort("id",row);
	FTCAL__clusters_x = FTCAL__clusters->getFloat("x",row);
	FTCAL__clusters_y = FTCAL__clusters->getFloat("y",row);
	FTCAL__clusters_z = FTCAL__clusters->getFloat("z",row);
	FTCAL__clusters_widthX = FTCAL__clusters->getFloat("widthX",row);
	FTCAL__clusters_widthY = FTCAL__clusters->getFloat("widthY",row);
	FTCAL__clusters_radius = FTCAL__clusters->getFloat("radius",row);
	FTCAL__clusters_time = FTCAL__clusters->getFloat("time",row);
	FTCAL__clusters_energy = FTCAL__clusters->getFloat("energy",row);
	FTCAL__clusters_maxEnergy = FTCAL__clusters->getFloat("maxEnergy",row);
	FTCAL__clusters_recEnergy = FTCAL__clusters->getFloat("recEnergy",row);
	return 0;
} 

int TIdentificatorCLAS12::get_RAW__scaler(int row){
	RAW__scaler_crate = RAW__scaler->getByte("crate",row);
	RAW__scaler_slot = RAW__scaler->getByte("slot",row);
	RAW__scaler_channel = RAW__scaler->getShort("channel",row);
	RAW__scaler_helicity = RAW__scaler->getByte("helicity",row);
	RAW__scaler_quartet = RAW__scaler->getByte("quartet",row);
	RAW__scaler_value = RAW__scaler->getLong("value",row);
	return 0;
} 

int TIdentificatorCLAS12::get_BMTRec__Clusters(int row){
	BMTRec__Clusters_ID = BMTRec__Clusters->getShort("ID",row);
	BMTRec__Clusters_sector = BMTRec__Clusters->getByte("sector",row);
	BMTRec__Clusters_layer = BMTRec__Clusters->getByte("layer",row);
	BMTRec__Clusters_size = BMTRec__Clusters->getShort("size",row);
	BMTRec__Clusters_ETot = BMTRec__Clusters->getFloat("ETot",row);
	BMTRec__Clusters_seedE = BMTRec__Clusters->getFloat("seedE",row);
	BMTRec__Clusters_seedStrip = BMTRec__Clusters->getInt("seedStrip",row);
	BMTRec__Clusters_centroid = BMTRec__Clusters->getFloat("centroid",row);
	BMTRec__Clusters_centroidResidual = BMTRec__Clusters->getFloat("centroidResidual",row);
	BMTRec__Clusters_seedResidual = BMTRec__Clusters->getFloat("seedResidual",row);
	BMTRec__Clusters_Hit1_ID = BMTRec__Clusters->getShort("Hit1_ID",row);
	BMTRec__Clusters_Hit2_ID = BMTRec__Clusters->getShort("Hit2_ID",row);
	BMTRec__Clusters_Hit3_ID = BMTRec__Clusters->getShort("Hit3_ID",row);
	BMTRec__Clusters_Hit4_ID = BMTRec__Clusters->getShort("Hit4_ID",row);
	BMTRec__Clusters_Hit5_ID = BMTRec__Clusters->getShort("Hit5_ID",row);
	BMTRec__Clusters_trkID = BMTRec__Clusters->getShort("trkID",row);
	return 0;
} 

int TIdentificatorCLAS12::get_TimeBasedTrkg__TBSegments(int row){
	TimeBasedTrkg__TBSegments_id = TimeBasedTrkg__TBSegments->getShort("id",row);
	TimeBasedTrkg__TBSegments_status = TimeBasedTrkg__TBSegments->getShort("status",row);
	TimeBasedTrkg__TBSegments_sector = TimeBasedTrkg__TBSegments->getByte("sector",row);
	TimeBasedTrkg__TBSegments_superlayer = TimeBasedTrkg__TBSegments->getByte("superlayer",row);
	TimeBasedTrkg__TBSegments_Cluster_ID = TimeBasedTrkg__TBSegments->getShort("Cluster_ID",row);
	TimeBasedTrkg__TBSegments_Hit1_ID = TimeBasedTrkg__TBSegments->getShort("Hit1_ID",row);
	TimeBasedTrkg__TBSegments_Hit2_ID = TimeBasedTrkg__TBSegments->getShort("Hit2_ID",row);
	TimeBasedTrkg__TBSegments_Hit3_ID = TimeBasedTrkg__TBSegments->getShort("Hit3_ID",row);
	TimeBasedTrkg__TBSegments_Hit4_ID = TimeBasedTrkg__TBSegments->getShort("Hit4_ID",row);
	TimeBasedTrkg__TBSegments_Hit5_ID = TimeBasedTrkg__TBSegments->getShort("Hit5_ID",row);
	TimeBasedTrkg__TBSegments_Hit6_ID = TimeBasedTrkg__TBSegments->getShort("Hit6_ID",row);
	TimeBasedTrkg__TBSegments_Hit7_ID = TimeBasedTrkg__TBSegments->getShort("Hit7_ID",row);
	TimeBasedTrkg__TBSegments_Hit8_ID = TimeBasedTrkg__TBSegments->getShort("Hit8_ID",row);
	TimeBasedTrkg__TBSegments_Hit9_ID = TimeBasedTrkg__TBSegments->getShort("Hit9_ID",row);
	TimeBasedTrkg__TBSegments_Hit10_ID = TimeBasedTrkg__TBSegments->getShort("Hit10_ID",row);
	TimeBasedTrkg__TBSegments_Hit11_ID = TimeBasedTrkg__TBSegments->getShort("Hit11_ID",row);
	TimeBasedTrkg__TBSegments_Hit12_ID = TimeBasedTrkg__TBSegments->getShort("Hit12_ID",row);
	TimeBasedTrkg__TBSegments_avgWire = TimeBasedTrkg__TBSegments->getFloat("avgWire",row);
	TimeBasedTrkg__TBSegments_fitChisqProb = TimeBasedTrkg__TBSegments->getFloat("fitChisqProb",row);
	TimeBasedTrkg__TBSegments_fitSlope = TimeBasedTrkg__TBSegments->getFloat("fitSlope",row);
	TimeBasedTrkg__TBSegments_fitSlopeErr = TimeBasedTrkg__TBSegments->getFloat("fitSlopeErr",row);
	TimeBasedTrkg__TBSegments_fitInterc = TimeBasedTrkg__TBSegments->getFloat("fitInterc",row);
	TimeBasedTrkg__TBSegments_fitIntercErr = TimeBasedTrkg__TBSegments->getFloat("fitIntercErr",row);
	TimeBasedTrkg__TBSegments_SegEndPoint1X = TimeBasedTrkg__TBSegments->getFloat("SegEndPoint1X",row);
	TimeBasedTrkg__TBSegments_SegEndPoint1Z = TimeBasedTrkg__TBSegments->getFloat("SegEndPoint1Z",row);
	TimeBasedTrkg__TBSegments_SegEndPoint2X = TimeBasedTrkg__TBSegments->getFloat("SegEndPoint2X",row);
	TimeBasedTrkg__TBSegments_SegEndPoint2Z = TimeBasedTrkg__TBSegments->getFloat("SegEndPoint2Z",row);
	TimeBasedTrkg__TBSegments_resiSum = TimeBasedTrkg__TBSegments->getFloat("resiSum",row);
	TimeBasedTrkg__TBSegments_timeSum = TimeBasedTrkg__TBSegments->getFloat("timeSum",row);
	TimeBasedTrkg__TBSegments_size = TimeBasedTrkg__TBSegments->getByte("size",row);
	return 0;
} 

int TIdentificatorCLAS12::get_RECHB__ForwardTagger(int row){
	RECHB__ForwardTagger_index = RECHB__ForwardTagger->getShort("index",row);
	RECHB__ForwardTagger_pindex = RECHB__ForwardTagger->getShort("pindex",row);
	RECHB__ForwardTagger_detector = RECHB__ForwardTagger->getByte("detector",row);
	RECHB__ForwardTagger_layer = RECHB__ForwardTagger->getByte("layer",row);
	RECHB__ForwardTagger_energy = RECHB__ForwardTagger->getFloat("energy",row);
	RECHB__ForwardTagger_time = RECHB__ForwardTagger->getFloat("time",row);
	RECHB__ForwardTagger_path = RECHB__ForwardTagger->getFloat("path",row);
	RECHB__ForwardTagger_chi2 = RECHB__ForwardTagger->getFloat("chi2",row);
	RECHB__ForwardTagger_x = RECHB__ForwardTagger->getFloat("x",row);
	RECHB__ForwardTagger_y = RECHB__ForwardTagger->getFloat("y",row);
	RECHB__ForwardTagger_z = RECHB__ForwardTagger->getFloat("z",row);
	RECHB__ForwardTagger_dx = RECHB__ForwardTagger->getFloat("dx",row);
	RECHB__ForwardTagger_dy = RECHB__ForwardTagger->getFloat("dy",row);
	RECHB__ForwardTagger_radius = RECHB__ForwardTagger->getFloat("radius",row);
	RECHB__ForwardTagger_size = RECHB__ForwardTagger->getShort("size",row);
	RECHB__ForwardTagger_status = RECHB__ForwardTagger->getShort("status",row);
	return 0;
} 

int TIdentificatorCLAS12::get_REC__Calorimeter(int row){
	REC__Calorimeter_index = REC__Calorimeter->getShort("index",row);
	REC__Calorimeter_pindex = REC__Calorimeter->getShort("pindex",row);
	REC__Calorimeter_detector = REC__Calorimeter->getByte("detector",row);
	REC__Calorimeter_sector = REC__Calorimeter->getByte("sector",row);
	REC__Calorimeter_layer = REC__Calorimeter->getByte("layer",row);
	REC__Calorimeter_energy = REC__Calorimeter->getFloat("energy",row);
	REC__Calorimeter_time = REC__Calorimeter->getFloat("time",row);
	REC__Calorimeter_path = REC__Calorimeter->getFloat("path",row);
	REC__Calorimeter_chi2 = REC__Calorimeter->getFloat("chi2",row);
	REC__Calorimeter_x = REC__Calorimeter->getFloat("x",row);
	REC__Calorimeter_y = REC__Calorimeter->getFloat("y",row);
	REC__Calorimeter_z = REC__Calorimeter->getFloat("z",row);
	REC__Calorimeter_hx = REC__Calorimeter->getFloat("hx",row);
	REC__Calorimeter_hy = REC__Calorimeter->getFloat("hy",row);
	REC__Calorimeter_hz = REC__Calorimeter->getFloat("hz",row);
	REC__Calorimeter_lu = REC__Calorimeter->getFloat("lu",row);
	REC__Calorimeter_lv = REC__Calorimeter->getFloat("lv",row);
	REC__Calorimeter_lw = REC__Calorimeter->getFloat("lw",row);
	REC__Calorimeter_du = REC__Calorimeter->getFloat("du",row);
	REC__Calorimeter_dv = REC__Calorimeter->getFloat("dv",row);
	REC__Calorimeter_dw = REC__Calorimeter->getFloat("dw",row);
	REC__Calorimeter_m2u = REC__Calorimeter->getFloat("m2u",row);
	REC__Calorimeter_m2v = REC__Calorimeter->getFloat("m2v",row);
	REC__Calorimeter_m2w = REC__Calorimeter->getFloat("m2w",row);
	REC__Calorimeter_m3u = REC__Calorimeter->getFloat("m3u",row);
	REC__Calorimeter_m3v = REC__Calorimeter->getFloat("m3v",row);
	REC__Calorimeter_m3w = REC__Calorimeter->getFloat("m3w",row);
	REC__Calorimeter_status = REC__Calorimeter->getShort("status",row);
	return 0;
} 

int TIdentificatorCLAS12::get_REC__CovMat(int row){
	REC__CovMat_index = REC__CovMat->getShort("index",row);
	REC__CovMat_pindex = REC__CovMat->getShort("pindex",row);
	REC__CovMat_C11 = REC__CovMat->getFloat("C11",row);
	REC__CovMat_C12 = REC__CovMat->getFloat("C12",row);
	REC__CovMat_C13 = REC__CovMat->getFloat("C13",row);
	REC__CovMat_C14 = REC__CovMat->getFloat("C14",row);
	REC__CovMat_C15 = REC__CovMat->getFloat("C15",row);
	REC__CovMat_C22 = REC__CovMat->getFloat("C22",row);
	REC__CovMat_C23 = REC__CovMat->getFloat("C23",row);
	REC__CovMat_C24 = REC__CovMat->getFloat("C24",row);
	REC__CovMat_C25 = REC__CovMat->getFloat("C25",row);
	REC__CovMat_C33 = REC__CovMat->getFloat("C33",row);
	REC__CovMat_C34 = REC__CovMat->getFloat("C34",row);
	REC__CovMat_C35 = REC__CovMat->getFloat("C35",row);
	REC__CovMat_C44 = REC__CovMat->getFloat("C44",row);
	REC__CovMat_C45 = REC__CovMat->getFloat("C45",row);
	REC__CovMat_C55 = REC__CovMat->getFloat("C55",row);
	return 0;
} 

int TIdentificatorCLAS12::get_RAW__epics(int row){
	RAW__epics_json = RAW__epics->getByte("json",row);
	return 0;
} 

int TIdentificatorCLAS12::get_REC__VertDoca(int row){
	REC__VertDoca_index1 = REC__VertDoca->getShort("index1",row);
	REC__VertDoca_index2 = REC__VertDoca->getShort("index2",row);
	REC__VertDoca_x = REC__VertDoca->getFloat("x",row);
	REC__VertDoca_y = REC__VertDoca->getFloat("y",row);
	REC__VertDoca_z = REC__VertDoca->getFloat("z",row);
	REC__VertDoca_x1 = REC__VertDoca->getFloat("x1",row);
	REC__VertDoca_y1 = REC__VertDoca->getFloat("y1",row);
	REC__VertDoca_z1 = REC__VertDoca->getFloat("z1",row);
	REC__VertDoca_cx1 = REC__VertDoca->getFloat("cx1",row);
	REC__VertDoca_cy1 = REC__VertDoca->getFloat("cy1",row);
	REC__VertDoca_cz1 = REC__VertDoca->getFloat("cz1",row);
	REC__VertDoca_x2 = REC__VertDoca->getFloat("x2",row);
	REC__VertDoca_y2 = REC__VertDoca->getFloat("y2",row);
	REC__VertDoca_z2 = REC__VertDoca->getFloat("z2",row);
	REC__VertDoca_cx2 = REC__VertDoca->getFloat("cx2",row);
	REC__VertDoca_cy2 = REC__VertDoca->getFloat("cy2",row);
	REC__VertDoca_cz2 = REC__VertDoca->getFloat("cz2",row);
	REC__VertDoca_r = REC__VertDoca->getFloat("r",row);
	return 0;
} 

int TIdentificatorCLAS12::get_FTHODO__clusters(int row){
	FTHODO__clusters_size = FTHODO__clusters->getShort("size",row);
	FTHODO__clusters_id = FTHODO__clusters->getShort("id",row);
	FTHODO__clusters_x = FTHODO__clusters->getFloat("x",row);
	FTHODO__clusters_y = FTHODO__clusters->getFloat("y",row);
	FTHODO__clusters_z = FTHODO__clusters->getFloat("z",row);
	FTHODO__clusters_widthX = FTHODO__clusters->getFloat("widthX",row);
	FTHODO__clusters_widthY = FTHODO__clusters->getFloat("widthY",row);
	FTHODO__clusters_radius = FTHODO__clusters->getFloat("radius",row);
	FTHODO__clusters_time = FTHODO__clusters->getFloat("time",row);
	FTHODO__clusters_energy = FTHODO__clusters->getFloat("energy",row);
	return 0;
} 

int TIdentificatorCLAS12::get_RUN__scaler(int row){
	RUN__scaler_fcupgated = RUN__scaler->getFloat("fcupgated",row);
	RUN__scaler_fcup = RUN__scaler->getFloat("fcup",row);
	RUN__scaler_livetime = RUN__scaler->getFloat("livetime",row);
	return 0;
} 

int TIdentificatorCLAS12::get_BAND__hits(int row){
	BAND__hits_id = BAND__hits->getShort("id",row);
	BAND__hits_sector = BAND__hits->getByte("sector",row);
	BAND__hits_layer = BAND__hits->getByte("layer",row);
	BAND__hits_component = BAND__hits->getShort("component",row);
	BAND__hits_meantimeTdc = BAND__hits->getFloat("meantimeTdc",row);
	BAND__hits_meantimeFadc = BAND__hits->getFloat("meantimeFadc",row);
	BAND__hits_difftimeTdc = BAND__hits->getFloat("difftimeTdc",row);
	BAND__hits_difftimeFadc = BAND__hits->getFloat("difftimeFadc",row);
	BAND__hits_adcLcorr = BAND__hits->getFloat("adcLcorr",row);
	BAND__hits_adcRcorr = BAND__hits->getFloat("adcRcorr",row);
	BAND__hits_tFadcLcorr = BAND__hits->getFloat("tFadcLcorr",row);
	BAND__hits_tFadcRcorr = BAND__hits->getFloat("tFadcRcorr",row);
	BAND__hits_tTdcLcorr = BAND__hits->getFloat("tTdcLcorr",row);
	BAND__hits_tTdcRcorr = BAND__hits->getFloat("tTdcRcorr",row);
	BAND__hits_x = BAND__hits->getFloat("x",row);
	BAND__hits_y = BAND__hits->getFloat("y",row);
	BAND__hits_z = BAND__hits->getFloat("z",row);
	BAND__hits_ux = BAND__hits->getFloat("ux",row);
	BAND__hits_uy = BAND__hits->getFloat("uy",row);
	BAND__hits_uz = BAND__hits->getFloat("uz",row);
	return 0;
} 

int TIdentificatorCLAS12::get_ECAL__hits(int row){
	ECAL__hits_id = ECAL__hits->getShort("id",row);
	ECAL__hits_status = ECAL__hits->getShort("status",row);
	ECAL__hits_sector = ECAL__hits->getByte("sector",row);
	ECAL__hits_layer = ECAL__hits->getByte("layer",row);
	ECAL__hits_strip = ECAL__hits->getByte("strip",row);
	ECAL__hits_peakid = ECAL__hits->getByte("peakid",row);
	ECAL__hits_energy = ECAL__hits->getFloat("energy",row);
	ECAL__hits_time = ECAL__hits->getFloat("time",row);
	return 0;
} 

int TIdentificatorCLAS12::get_ECAL__clusters(int row){
	ECAL__clusters_id = ECAL__clusters->getShort("id",row);
	ECAL__clusters_status = ECAL__clusters->getShort("status",row);
	ECAL__clusters_sector = ECAL__clusters->getByte("sector",row);
	ECAL__clusters_layer = ECAL__clusters->getByte("layer",row);
	ECAL__clusters_x = ECAL__clusters->getFloat("x",row);
	ECAL__clusters_y = ECAL__clusters->getFloat("y",row);
	ECAL__clusters_z = ECAL__clusters->getFloat("z",row);
	ECAL__clusters_energy = ECAL__clusters->getFloat("energy",row);
	ECAL__clusters_time = ECAL__clusters->getFloat("time",row);
	ECAL__clusters_widthU = ECAL__clusters->getFloat("widthU",row);
	ECAL__clusters_widthV = ECAL__clusters->getFloat("widthV",row);
	ECAL__clusters_widthW = ECAL__clusters->getFloat("widthW",row);
	ECAL__clusters_idU = ECAL__clusters->getByte("idU",row);
	ECAL__clusters_idV = ECAL__clusters->getByte("idV",row);
	ECAL__clusters_idW = ECAL__clusters->getByte("idW",row);
	ECAL__clusters_coordU = ECAL__clusters->getInt("coordU",row);
	ECAL__clusters_coordV = ECAL__clusters->getInt("coordV",row);
	ECAL__clusters_coordW = ECAL__clusters->getInt("coordW",row);
	return 0;
} 

int TIdentificatorCLAS12::get_FT__particles(int row){
	FT__particles_id = FT__particles->getShort("id",row);
	FT__particles_charge = FT__particles->getByte("charge",row);
	FT__particles_energy = FT__particles->getFloat("energy",row);
	FT__particles_cx = FT__particles->getFloat("cx",row);
	FT__particles_cy = FT__particles->getFloat("cy",row);
	FT__particles_cz = FT__particles->getFloat("cz",row);
	FT__particles_time = FT__particles->getFloat("time",row);
	FT__particles_calID = FT__particles->getShort("calID",row);
	FT__particles_hodoID = FT__particles->getShort("hodoID",row);
	FT__particles_trkID = FT__particles->getShort("trkID",row);
	return 0;
} 

int TIdentificatorCLAS12::get_HitBasedTrkg__HBClusters(int row){
	HitBasedTrkg__HBClusters_id = HitBasedTrkg__HBClusters->getShort("id",row);
	HitBasedTrkg__HBClusters_status = HitBasedTrkg__HBClusters->getShort("status",row);
	HitBasedTrkg__HBClusters_sector = HitBasedTrkg__HBClusters->getByte("sector",row);
	HitBasedTrkg__HBClusters_superlayer = HitBasedTrkg__HBClusters->getByte("superlayer",row);
	HitBasedTrkg__HBClusters_Hit1_ID = HitBasedTrkg__HBClusters->getShort("Hit1_ID",row);
	HitBasedTrkg__HBClusters_Hit2_ID = HitBasedTrkg__HBClusters->getShort("Hit2_ID",row);
	HitBasedTrkg__HBClusters_Hit3_ID = HitBasedTrkg__HBClusters->getShort("Hit3_ID",row);
	HitBasedTrkg__HBClusters_Hit4_ID = HitBasedTrkg__HBClusters->getShort("Hit4_ID",row);
	HitBasedTrkg__HBClusters_Hit5_ID = HitBasedTrkg__HBClusters->getShort("Hit5_ID",row);
	HitBasedTrkg__HBClusters_Hit6_ID = HitBasedTrkg__HBClusters->getShort("Hit6_ID",row);
	HitBasedTrkg__HBClusters_Hit7_ID = HitBasedTrkg__HBClusters->getShort("Hit7_ID",row);
	HitBasedTrkg__HBClusters_Hit8_ID = HitBasedTrkg__HBClusters->getShort("Hit8_ID",row);
	HitBasedTrkg__HBClusters_Hit9_ID = HitBasedTrkg__HBClusters->getShort("Hit9_ID",row);
	HitBasedTrkg__HBClusters_Hit10_ID = HitBasedTrkg__HBClusters->getShort("Hit10_ID",row);
	HitBasedTrkg__HBClusters_Hit11_ID = HitBasedTrkg__HBClusters->getShort("Hit11_ID",row);
	HitBasedTrkg__HBClusters_Hit12_ID = HitBasedTrkg__HBClusters->getShort("Hit12_ID",row);
	HitBasedTrkg__HBClusters_avgWire = HitBasedTrkg__HBClusters->getFloat("avgWire",row);
	HitBasedTrkg__HBClusters_fitChisqProb = HitBasedTrkg__HBClusters->getFloat("fitChisqProb",row);
	HitBasedTrkg__HBClusters_fitSlope = HitBasedTrkg__HBClusters->getFloat("fitSlope",row);
	HitBasedTrkg__HBClusters_fitSlopeErr = HitBasedTrkg__HBClusters->getFloat("fitSlopeErr",row);
	HitBasedTrkg__HBClusters_fitInterc = HitBasedTrkg__HBClusters->getFloat("fitInterc",row);
	HitBasedTrkg__HBClusters_fitIntercErr = HitBasedTrkg__HBClusters->getFloat("fitIntercErr",row);
	HitBasedTrkg__HBClusters_size = HitBasedTrkg__HBClusters->getByte("size",row);
	return 0;
} 

int TIdentificatorCLAS12::get_REC__RICH(int row){
	REC__RICH_index = REC__RICH->getShort("index",row);
	REC__RICH_pindex = REC__RICH->getShort("pindex",row);
	REC__RICH_detector = REC__RICH->getByte("detector",row);
	REC__RICH_sector = REC__RICH->getByte("sector",row);
	REC__RICH_layer = REC__RICH->getByte("layer",row);
	REC__RICH_energy = REC__RICH->getFloat("energy",row);
	REC__RICH_time = REC__RICH->getFloat("time",row);
	REC__RICH_path = REC__RICH->getFloat("path",row);
	REC__RICH_chi2 = REC__RICH->getFloat("chi2",row);
	REC__RICH_x = REC__RICH->getFloat("x",row);
	REC__RICH_y = REC__RICH->getFloat("y",row);
	REC__RICH_z = REC__RICH->getFloat("z",row);
	REC__RICH_hx = REC__RICH->getFloat("hx",row);
	REC__RICH_hy = REC__RICH->getFloat("hy",row);
	REC__RICH_hz = REC__RICH->getFloat("hz",row);
	REC__RICH_status = REC__RICH->getShort("status",row);
	return 0;
} 

int TIdentificatorCLAS12::get_REC__Track(int row){
	REC__Track_index = REC__Track->getShort("index",row);
	REC__Track_pindex = REC__Track->getShort("pindex",row);
	REC__Track_detector = REC__Track->getByte("detector",row);
	REC__Track_sector = REC__Track->getByte("sector",row);
	REC__Track_status = REC__Track->getShort("status",row);
	REC__Track_q = REC__Track->getByte("q",row);
	REC__Track_chi2 = REC__Track->getFloat("chi2",row);
	REC__Track_NDF = REC__Track->getShort("NDF",row);
	return 0;
} 

int TIdentificatorCLAS12::get_RECHB__Event(int row){
	RECHB__Event_category = RECHB__Event->getLong("category",row);
	RECHB__Event_topology = RECHB__Event->getLong("topology",row);
	RECHB__Event_beamCharge = RECHB__Event->getFloat("beamCharge",row);
	RECHB__Event_liveTime = RECHB__Event->getDouble("liveTime",row);
	RECHB__Event_startTime = RECHB__Event->getFloat("startTime",row);
	RECHB__Event_RFTime = RECHB__Event->getFloat("RFTime",row);
	RECHB__Event_helicity = RECHB__Event->getByte("helicity",row);
	RECHB__Event_helicityRaw = RECHB__Event->getByte("helicityRaw",row);
	RECHB__Event_procTime = RECHB__Event->getFloat("procTime",row);
	return 0;
} 

int TIdentificatorCLAS12::get_FTOF__hbhits(int row){
	FTOF__hbhits_id = FTOF__hbhits->getShort("id",row);
	FTOF__hbhits_status = FTOF__hbhits->getShort("status",row);
	FTOF__hbhits_trackid = FTOF__hbhits->getShort("trackid",row);
	FTOF__hbhits_sector = FTOF__hbhits->getByte("sector",row);
	FTOF__hbhits_layer = FTOF__hbhits->getByte("layer",row);
	FTOF__hbhits_component = FTOF__hbhits->getShort("component",row);
	FTOF__hbhits_energy = FTOF__hbhits->getFloat("energy",row);
	FTOF__hbhits_time = FTOF__hbhits->getFloat("time",row);
	FTOF__hbhits_energy_unc = FTOF__hbhits->getFloat("energy_unc",row);
	FTOF__hbhits_time_unc = FTOF__hbhits->getFloat("time_unc",row);
	FTOF__hbhits_x = FTOF__hbhits->getFloat("x",row);
	FTOF__hbhits_y = FTOF__hbhits->getFloat("y",row);
	FTOF__hbhits_z = FTOF__hbhits->getFloat("z",row);
	FTOF__hbhits_x_unc = FTOF__hbhits->getFloat("x_unc",row);
	FTOF__hbhits_y_unc = FTOF__hbhits->getFloat("y_unc",row);
	FTOF__hbhits_z_unc = FTOF__hbhits->getFloat("z_unc",row);
	FTOF__hbhits_tx = FTOF__hbhits->getFloat("tx",row);
	FTOF__hbhits_ty = FTOF__hbhits->getFloat("ty",row);
	FTOF__hbhits_tz = FTOF__hbhits->getFloat("tz",row);
	FTOF__hbhits_adc_idx1 = FTOF__hbhits->getShort("adc_idx1",row);
	FTOF__hbhits_adc_idx2 = FTOF__hbhits->getShort("adc_idx2",row);
	FTOF__hbhits_tdc_idx1 = FTOF__hbhits->getShort("tdc_idx1",row);
	FTOF__hbhits_tdc_idx2 = FTOF__hbhits->getShort("tdc_idx2",row);
	FTOF__hbhits_pathLength = FTOF__hbhits->getFloat("pathLength",row);
	FTOF__hbhits_pathLengthThruBar = FTOF__hbhits->getFloat("pathLengthThruBar",row);
	return 0;
} 

int TIdentificatorCLAS12::get_RTPC__adc(int row){
	RTPC__adc_sector = RTPC__adc->getByte("sector",row);
	RTPC__adc_layer = RTPC__adc->getByte("layer",row);
	RTPC__adc_component = RTPC__adc->getShort("component",row);
	RTPC__adc_order = RTPC__adc->getByte("order",row);
	RTPC__adc_ADC = RTPC__adc->getInt("ADC",row);
	RTPC__adc_time = RTPC__adc->getFloat("time",row);
	RTPC__adc_ped = RTPC__adc->getShort("ped",row);
	return 0;
} 

int TIdentificatorCLAS12::get_FTCAL__adc(int row){
	FTCAL__adc_sector = FTCAL__adc->getByte("sector",row);
	FTCAL__adc_layer = FTCAL__adc->getByte("layer",row);
	FTCAL__adc_component = FTCAL__adc->getShort("component",row);
	FTCAL__adc_order = FTCAL__adc->getByte("order",row);
	FTCAL__adc_ADC = FTCAL__adc->getInt("ADC",row);
	FTCAL__adc_time = FTCAL__adc->getFloat("time",row);
	FTCAL__adc_ped = FTCAL__adc->getShort("ped",row);
	return 0;
} 

int TIdentificatorCLAS12::get_CTOF__rawhits(int row){
	CTOF__rawhits_id = CTOF__rawhits->getShort("id",row);
	CTOF__rawhits_status = CTOF__rawhits->getShort("status",row);
	CTOF__rawhits_component = CTOF__rawhits->getShort("component",row);
	CTOF__rawhits_energy_up = CTOF__rawhits->getFloat("energy_up",row);
	CTOF__rawhits_energy_down = CTOF__rawhits->getFloat("energy_down",row);
	CTOF__rawhits_time_up = CTOF__rawhits->getFloat("time_up",row);
	CTOF__rawhits_time_down = CTOF__rawhits->getFloat("time_down",row);
	CTOF__rawhits_energy_up_unc = CTOF__rawhits->getFloat("energy_up_unc",row);
	CTOF__rawhits_energy_down_unc = CTOF__rawhits->getFloat("energy_down_unc",row);
	CTOF__rawhits_time_up_unc = CTOF__rawhits->getFloat("time_up_unc",row);
	CTOF__rawhits_time_down_unc = CTOF__rawhits->getFloat("time_down_unc",row);
	return 0;
} 

int TIdentificatorCLAS12::get_CND__adc(int row){
	CND__adc_sector = CND__adc->getByte("sector",row);
	CND__adc_layer = CND__adc->getByte("layer",row);
	CND__adc_component = CND__adc->getShort("component",row);
	CND__adc_order = CND__adc->getByte("order",row);
	CND__adc_ADC = CND__adc->getInt("ADC",row);
	CND__adc_time = CND__adc->getFloat("time",row);
	CND__adc_ped = CND__adc->getShort("ped",row);
	return 0;
} 

int TIdentificatorCLAS12::get_FTOF__hits(int row){
	FTOF__hits_id = FTOF__hits->getShort("id",row);
	FTOF__hits_status = FTOF__hits->getShort("status",row);
	FTOF__hits_trackid = FTOF__hits->getShort("trackid",row);
	FTOF__hits_sector = FTOF__hits->getByte("sector",row);
	FTOF__hits_layer = FTOF__hits->getByte("layer",row);
	FTOF__hits_component = FTOF__hits->getShort("component",row);
	FTOF__hits_energy = FTOF__hits->getFloat("energy",row);
	FTOF__hits_time = FTOF__hits->getFloat("time",row);
	FTOF__hits_energy_unc = FTOF__hits->getFloat("energy_unc",row);
	FTOF__hits_time_unc = FTOF__hits->getFloat("time_unc",row);
	FTOF__hits_x = FTOF__hits->getFloat("x",row);
	FTOF__hits_y = FTOF__hits->getFloat("y",row);
	FTOF__hits_z = FTOF__hits->getFloat("z",row);
	FTOF__hits_x_unc = FTOF__hits->getFloat("x_unc",row);
	FTOF__hits_y_unc = FTOF__hits->getFloat("y_unc",row);
	FTOF__hits_z_unc = FTOF__hits->getFloat("z_unc",row);
	FTOF__hits_tx = FTOF__hits->getFloat("tx",row);
	FTOF__hits_ty = FTOF__hits->getFloat("ty",row);
	FTOF__hits_tz = FTOF__hits->getFloat("tz",row);
	FTOF__hits_adc_idx1 = FTOF__hits->getShort("adc_idx1",row);
	FTOF__hits_adc_idx2 = FTOF__hits->getShort("adc_idx2",row);
	FTOF__hits_tdc_idx1 = FTOF__hits->getShort("tdc_idx1",row);
	FTOF__hits_tdc_idx2 = FTOF__hits->getShort("tdc_idx2",row);
	FTOF__hits_pathLength = FTOF__hits->getFloat("pathLength",row);
	FTOF__hits_pathLengthThruBar = FTOF__hits->getFloat("pathLengthThruBar",row);
	return 0;
} 

int TIdentificatorCLAS12::get_FMTRec__Clusters(int row){
	FMTRec__Clusters_ID = FMTRec__Clusters->getShort("ID",row);
	FMTRec__Clusters_sector = FMTRec__Clusters->getByte("sector",row);
	FMTRec__Clusters_layer = FMTRec__Clusters->getByte("layer",row);
	FMTRec__Clusters_size = FMTRec__Clusters->getShort("size",row);
	FMTRec__Clusters_ETot = FMTRec__Clusters->getFloat("ETot",row);
	FMTRec__Clusters_seedE = FMTRec__Clusters->getFloat("seedE",row);
	FMTRec__Clusters_seedStrip = FMTRec__Clusters->getInt("seedStrip",row);
	FMTRec__Clusters_centroid = FMTRec__Clusters->getFloat("centroid",row);
	FMTRec__Clusters_centroidResidual = FMTRec__Clusters->getFloat("centroidResidual",row);
	FMTRec__Clusters_seedResidual = FMTRec__Clusters->getFloat("seedResidual",row);
	FMTRec__Clusters_Hit1_ID = FMTRec__Clusters->getShort("Hit1_ID",row);
	FMTRec__Clusters_Hit2_ID = FMTRec__Clusters->getShort("Hit2_ID",row);
	FMTRec__Clusters_Hit3_ID = FMTRec__Clusters->getShort("Hit3_ID",row);
	FMTRec__Clusters_Hit4_ID = FMTRec__Clusters->getShort("Hit4_ID",row);
	FMTRec__Clusters_Hit5_ID = FMTRec__Clusters->getShort("Hit5_ID",row);
	FMTRec__Clusters_trkID = FMTRec__Clusters->getShort("trkID",row);
	return 0;
} 

int TIdentificatorCLAS12::get_HEL__adc(int row){
	HEL__adc_sector = HEL__adc->getByte("sector",row);
	HEL__adc_layer = HEL__adc->getByte("layer",row);
	HEL__adc_component = HEL__adc->getShort("component",row);
	HEL__adc_order = HEL__adc->getByte("order",row);
	HEL__adc_ADC = HEL__adc->getInt("ADC",row);
	HEL__adc_time = HEL__adc->getFloat("time",row);
	HEL__adc_ped = HEL__adc->getShort("ped",row);
	return 0;
} 

int TIdentificatorCLAS12::get_DETECTOR__ccpb(int row){
	DETECTOR__ccpb_sector = DETECTOR__ccpb->getByte("sector",row);
	DETECTOR__ccpb_nphe = DETECTOR__ccpb->getFloat("nphe",row);
	DETECTOR__ccpb_time = DETECTOR__ccpb->getFloat("time",row);
	DETECTOR__ccpb_path = DETECTOR__ccpb->getFloat("path",row);
	return 0;
} 

int TIdentificatorCLAS12::get_REC__ForwardTagger(int row){
	REC__ForwardTagger_index = REC__ForwardTagger->getShort("index",row);
	REC__ForwardTagger_pindex = REC__ForwardTagger->getShort("pindex",row);
	REC__ForwardTagger_detector = REC__ForwardTagger->getByte("detector",row);
	REC__ForwardTagger_layer = REC__ForwardTagger->getByte("layer",row);
	REC__ForwardTagger_energy = REC__ForwardTagger->getFloat("energy",row);
	REC__ForwardTagger_time = REC__ForwardTagger->getFloat("time",row);
	REC__ForwardTagger_path = REC__ForwardTagger->getFloat("path",row);
	REC__ForwardTagger_chi2 = REC__ForwardTagger->getFloat("chi2",row);
	REC__ForwardTagger_x = REC__ForwardTagger->getFloat("x",row);
	REC__ForwardTagger_y = REC__ForwardTagger->getFloat("y",row);
	REC__ForwardTagger_z = REC__ForwardTagger->getFloat("z",row);
	REC__ForwardTagger_dx = REC__ForwardTagger->getFloat("dx",row);
	REC__ForwardTagger_dy = REC__ForwardTagger->getFloat("dy",row);
	REC__ForwardTagger_radius = REC__ForwardTagger->getFloat("radius",row);
	REC__ForwardTagger_size = REC__ForwardTagger->getShort("size",row);
	REC__ForwardTagger_status = REC__ForwardTagger->getShort("status",row);
	return 0;
} 

int TIdentificatorCLAS12::get_FTTRK__adc(int row){
	FTTRK__adc_sector = FTTRK__adc->getByte("sector",row);
	FTTRK__adc_layer = FTTRK__adc->getByte("layer",row);
	FTTRK__adc_component = FTTRK__adc->getShort("component",row);
	FTTRK__adc_order = FTTRK__adc->getByte("order",row);
	FTTRK__adc_ADC = FTTRK__adc->getInt("ADC",row);
	FTTRK__adc_time = FTTRK__adc->getFloat("time",row);
	FTTRK__adc_ped = FTTRK__adc->getShort("ped",row);
	FTTRK__adc_integral = FTTRK__adc->getInt("integral",row);
	FTTRK__adc_timestamp = FTTRK__adc->getLong("timestamp",row);
	return 0;
} 

int TIdentificatorCLAS12::get_FTOF__rawhits(int row){
	FTOF__rawhits_id = FTOF__rawhits->getShort("id",row);
	FTOF__rawhits_status = FTOF__rawhits->getShort("status",row);
	FTOF__rawhits_sector = FTOF__rawhits->getByte("sector",row);
	FTOF__rawhits_layer = FTOF__rawhits->getByte("layer",row);
	FTOF__rawhits_component = FTOF__rawhits->getShort("component",row);
	FTOF__rawhits_energy_left = FTOF__rawhits->getFloat("energy_left",row);
	FTOF__rawhits_energy_right = FTOF__rawhits->getFloat("energy_right",row);
	FTOF__rawhits_time_left = FTOF__rawhits->getFloat("time_left",row);
	FTOF__rawhits_time_right = FTOF__rawhits->getFloat("time_right",row);
	FTOF__rawhits_energy_left_unc = FTOF__rawhits->getFloat("energy_left_unc",row);
	FTOF__rawhits_energy_right_unc = FTOF__rawhits->getFloat("energy_right_unc",row);
	FTOF__rawhits_time_left_unc = FTOF__rawhits->getFloat("time_left_unc",row);
	FTOF__rawhits_time_right_unc = FTOF__rawhits->getFloat("time_right_unc",row);
	return 0;
} 

int TIdentificatorCLAS12::get_RECHB__TrackCross(int row){
	RECHB__TrackCross_index = RECHB__TrackCross->getShort("index",row);
	RECHB__TrackCross_pindex = RECHB__TrackCross->getShort("pindex",row);
	RECHB__TrackCross_detector = RECHB__TrackCross->getByte("detector",row);
	RECHB__TrackCross_sector = RECHB__TrackCross->getByte("sector",row);
	RECHB__TrackCross_layer = RECHB__TrackCross->getByte("layer",row);
	RECHB__TrackCross_c_x = RECHB__TrackCross->getFloat("c_x",row);
	RECHB__TrackCross_c_y = RECHB__TrackCross->getFloat("c_y",row);
	RECHB__TrackCross_c_z = RECHB__TrackCross->getFloat("c_z",row);
	RECHB__TrackCross_c_ux = RECHB__TrackCross->getFloat("c_ux",row);
	RECHB__TrackCross_c_uy = RECHB__TrackCross->getFloat("c_uy",row);
	RECHB__TrackCross_c_uz = RECHB__TrackCross->getFloat("c_uz",row);
	RECHB__TrackCross_status = RECHB__TrackCross->getShort("status",row);
	return 0;
} 

int TIdentificatorCLAS12::get_DETECTOR__scpb(int row){
	DETECTOR__scpb_sector = DETECTOR__scpb->getByte("sector",row);
	DETECTOR__scpb_paddle = DETECTOR__scpb->getByte("paddle",row);
	DETECTOR__scpb_edep = DETECTOR__scpb->getFloat("edep",row);
	DETECTOR__scpb_time = DETECTOR__scpb->getFloat("time",row);
	DETECTOR__scpb_path = DETECTOR__scpb->getFloat("path",row);
	return 0;
} 

int TIdentificatorCLAS12::get_FTHODO__hits(int row){
	FTHODO__hits_sector = FTHODO__hits->getByte("sector",row);
	FTHODO__hits_layer = FTHODO__hits->getByte("layer",row);
	FTHODO__hits_component = FTHODO__hits->getShort("component",row);
	FTHODO__hits_x = FTHODO__hits->getFloat("x",row);
	FTHODO__hits_y = FTHODO__hits->getFloat("y",row);
	FTHODO__hits_z = FTHODO__hits->getFloat("z",row);
	FTHODO__hits_energy = FTHODO__hits->getFloat("energy",row);
	FTHODO__hits_time = FTHODO__hits->getFloat("time",row);
	FTHODO__hits_hitID = FTHODO__hits->getShort("hitID",row);
	FTHODO__hits_clusterID = FTHODO__hits->getShort("clusterID",row);
	return 0;
} 

int TIdentificatorCLAS12::get_RECFT__Event(int row){
	RECFT__Event_EvCAT = RECFT__Event->getShort("EvCAT",row);
	RECFT__Event_startTime = RECFT__Event->getFloat("startTime",row);
	return 0;
} 

int TIdentificatorCLAS12::get_CTOF__adc(int row){
	CTOF__adc_sector = CTOF__adc->getByte("sector",row);
	CTOF__adc_layer = CTOF__adc->getByte("layer",row);
	CTOF__adc_component = CTOF__adc->getShort("component",row);
	CTOF__adc_order = CTOF__adc->getByte("order",row);
	CTOF__adc_ADC = CTOF__adc->getInt("ADC",row);
	CTOF__adc_time = CTOF__adc->getFloat("time",row);
	CTOF__adc_ped = CTOF__adc->getShort("ped",row);
	return 0;
} 

int TIdentificatorCLAS12::get_EVENT__detector(int row){
	EVENT__detector_scsector = EVENT__detector->getInt("scsector",row);
	EVENT__detector_scpaddle = EVENT__detector->getInt("scpaddle",row);
	EVENT__detector_ecsector = EVENT__detector->getInt("ecsector",row);
	EVENT__detector_ccnphe = EVENT__detector->getFloat("ccnphe",row);
	EVENT__detector_sctime = EVENT__detector->getFloat("sctime",row);
	EVENT__detector_scpath = EVENT__detector->getFloat("scpath",row);
	EVENT__detector_ectime = EVENT__detector->getFloat("ectime",row);
	EVENT__detector_ecpath = EVENT__detector->getFloat("ecpath",row);
	EVENT__detector_ecin = EVENT__detector->getFloat("ecin",row);
	EVENT__detector_ecout = EVENT__detector->getFloat("ecout",row);
	EVENT__detector_ectot = EVENT__detector->getFloat("ectot",row);
	EVENT__detector_ecu = EVENT__detector->getFloat("ecu",row);
	EVENT__detector_ecv = EVENT__detector->getFloat("ecv",row);
	EVENT__detector_ecw = EVENT__detector->getFloat("ecw",row);
	return 0;
} 

int TIdentificatorCLAS12::get_DC__doca(int row){
	DC__doca_LR = DC__doca->getByte("LR",row);
	DC__doca_doca = DC__doca->getFloat("doca",row);
	DC__doca_sdoca = DC__doca->getFloat("sdoca",row);
	DC__doca_time = DC__doca->getFloat("time",row);
	DC__doca_stime = DC__doca->getFloat("stime",row);
	return 0;
} 

int TIdentificatorCLAS12::get_HTCC__tdc(int row){
	HTCC__tdc_sector = HTCC__tdc->getByte("sector",row);
	HTCC__tdc_layer = HTCC__tdc->getByte("layer",row);
	HTCC__tdc_component = HTCC__tdc->getShort("component",row);
	HTCC__tdc_order = HTCC__tdc->getByte("order",row);
	HTCC__tdc_TDC = HTCC__tdc->getInt("TDC",row);
	return 0;
} 

int TIdentificatorCLAS12::get_HEL__online(int row){
	HEL__online_helicity = HEL__online->getByte("helicity",row);
	HEL__online_helicityRaw = HEL__online->getByte("helicityRaw",row);
	return 0;
} 

int TIdentificatorCLAS12::get_REC__Event(int row){
	REC__Event_category = REC__Event->getLong("category",row);
	REC__Event_topology = REC__Event->getLong("topology",row);
	REC__Event_beamCharge = REC__Event->getFloat("beamCharge",row);
	REC__Event_liveTime = REC__Event->getDouble("liveTime",row);
	REC__Event_startTime = REC__Event->getFloat("startTime",row);
	REC__Event_RFTime = REC__Event->getFloat("RFTime",row);
	REC__Event_helicity = REC__Event->getByte("helicity",row);
	REC__Event_helicityRaw = REC__Event->getByte("helicityRaw",row);
	REC__Event_procTime = REC__Event->getFloat("procTime",row);
	return 0;
} 

int TIdentificatorCLAS12::get_REC__Particle(int row){
	REC__Particle_pid = REC__Particle->getInt("pid",row);
	REC__Particle_px = REC__Particle->getFloat("px",row);
	REC__Particle_py = REC__Particle->getFloat("py",row);
	REC__Particle_pz = REC__Particle->getFloat("pz",row);
	REC__Particle_vx = REC__Particle->getFloat("vx",row);
	REC__Particle_vy = REC__Particle->getFloat("vy",row);
	REC__Particle_vz = REC__Particle->getFloat("vz",row);
	REC__Particle_charge = REC__Particle->getByte("charge",row);
	REC__Particle_beta = REC__Particle->getFloat("beta",row);
	REC__Particle_chi2pid = REC__Particle->getFloat("chi2pid",row);
	REC__Particle_status = REC__Particle->getShort("status",row);
	return 0;
} 

int TIdentificatorCLAS12::get_REC__Traj(int row){
	REC__Traj_pindex = REC__Traj->getShort("pindex",row);
	REC__Traj_index = REC__Traj->getShort("index",row);
	REC__Traj_detector = REC__Traj->getByte("detector",row);
	REC__Traj_layer = REC__Traj->getByte("layer",row);
	REC__Traj_x = REC__Traj->getFloat("x",row);
	REC__Traj_y = REC__Traj->getFloat("y",row);
	REC__Traj_z = REC__Traj->getFloat("z",row);
	REC__Traj_cx = REC__Traj->getFloat("cx",row);
	REC__Traj_cy = REC__Traj->getFloat("cy",row);
	REC__Traj_cz = REC__Traj->getFloat("cz",row);
	REC__Traj_path = REC__Traj->getFloat("path",row);
	return 0;
} 

int TIdentificatorCLAS12::get_FTOF__tdc(int row){
	FTOF__tdc_sector = FTOF__tdc->getByte("sector",row);
	FTOF__tdc_layer = FTOF__tdc->getByte("layer",row);
	FTOF__tdc_component = FTOF__tdc->getShort("component",row);
	FTOF__tdc_order = FTOF__tdc->getByte("order",row);
	FTOF__tdc_TDC = FTOF__tdc->getInt("TDC",row);
	return 0;
} 

int TIdentificatorCLAS12::get_HitBasedTrkg__HBSegments(int row){
	HitBasedTrkg__HBSegments_id = HitBasedTrkg__HBSegments->getShort("id",row);
	HitBasedTrkg__HBSegments_status = HitBasedTrkg__HBSegments->getShort("status",row);
	HitBasedTrkg__HBSegments_sector = HitBasedTrkg__HBSegments->getByte("sector",row);
	HitBasedTrkg__HBSegments_superlayer = HitBasedTrkg__HBSegments->getByte("superlayer",row);
	HitBasedTrkg__HBSegments_Cluster_ID = HitBasedTrkg__HBSegments->getShort("Cluster_ID",row);
	HitBasedTrkg__HBSegments_Hit1_ID = HitBasedTrkg__HBSegments->getShort("Hit1_ID",row);
	HitBasedTrkg__HBSegments_Hit2_ID = HitBasedTrkg__HBSegments->getShort("Hit2_ID",row);
	HitBasedTrkg__HBSegments_Hit3_ID = HitBasedTrkg__HBSegments->getShort("Hit3_ID",row);
	HitBasedTrkg__HBSegments_Hit4_ID = HitBasedTrkg__HBSegments->getShort("Hit4_ID",row);
	HitBasedTrkg__HBSegments_Hit5_ID = HitBasedTrkg__HBSegments->getShort("Hit5_ID",row);
	HitBasedTrkg__HBSegments_Hit6_ID = HitBasedTrkg__HBSegments->getShort("Hit6_ID",row);
	HitBasedTrkg__HBSegments_Hit7_ID = HitBasedTrkg__HBSegments->getShort("Hit7_ID",row);
	HitBasedTrkg__HBSegments_Hit8_ID = HitBasedTrkg__HBSegments->getShort("Hit8_ID",row);
	HitBasedTrkg__HBSegments_Hit9_ID = HitBasedTrkg__HBSegments->getShort("Hit9_ID",row);
	HitBasedTrkg__HBSegments_Hit10_ID = HitBasedTrkg__HBSegments->getShort("Hit10_ID",row);
	HitBasedTrkg__HBSegments_Hit11_ID = HitBasedTrkg__HBSegments->getShort("Hit11_ID",row);
	HitBasedTrkg__HBSegments_Hit12_ID = HitBasedTrkg__HBSegments->getShort("Hit12_ID",row);
	HitBasedTrkg__HBSegments_avgWire = HitBasedTrkg__HBSegments->getFloat("avgWire",row);
	HitBasedTrkg__HBSegments_fitChisqProb = HitBasedTrkg__HBSegments->getFloat("fitChisqProb",row);
	HitBasedTrkg__HBSegments_fitSlope = HitBasedTrkg__HBSegments->getFloat("fitSlope",row);
	HitBasedTrkg__HBSegments_fitSlopeErr = HitBasedTrkg__HBSegments->getFloat("fitSlopeErr",row);
	HitBasedTrkg__HBSegments_fitInterc = HitBasedTrkg__HBSegments->getFloat("fitInterc",row);
	HitBasedTrkg__HBSegments_fitIntercErr = HitBasedTrkg__HBSegments->getFloat("fitIntercErr",row);
	HitBasedTrkg__HBSegments_SegEndPoint1X = HitBasedTrkg__HBSegments->getFloat("SegEndPoint1X",row);
	HitBasedTrkg__HBSegments_SegEndPoint1Z = HitBasedTrkg__HBSegments->getFloat("SegEndPoint1Z",row);
	HitBasedTrkg__HBSegments_SegEndPoint2X = HitBasedTrkg__HBSegments->getFloat("SegEndPoint2X",row);
	HitBasedTrkg__HBSegments_SegEndPoint2Z = HitBasedTrkg__HBSegments->getFloat("SegEndPoint2Z",row);
	HitBasedTrkg__HBSegments_size = HitBasedTrkg__HBSegments->getByte("size",row);
	return 0;
} 

int TIdentificatorCLAS12::get_FTCAL__hits(int row){
	FTCAL__hits_idx = FTCAL__hits->getByte("idx",row);
	FTCAL__hits_idy = FTCAL__hits->getByte("idy",row);
	FTCAL__hits_x = FTCAL__hits->getFloat("x",row);
	FTCAL__hits_y = FTCAL__hits->getFloat("y",row);
	FTCAL__hits_z = FTCAL__hits->getFloat("z",row);
	FTCAL__hits_energy = FTCAL__hits->getFloat("energy",row);
	FTCAL__hits_time = FTCAL__hits->getFloat("time",row);
	FTCAL__hits_hitID = FTCAL__hits->getShort("hitID",row);
	FTCAL__hits_clusterID = FTCAL__hits->getShort("clusterID",row);
	return 0;
} 

int TIdentificatorCLAS12::get_FMTRec__Hits(int row){
	FMTRec__Hits_ID = FMTRec__Hits->getShort("ID",row);
	FMTRec__Hits_sector = FMTRec__Hits->getByte("sector",row);
	FMTRec__Hits_layer = FMTRec__Hits->getByte("layer",row);
	FMTRec__Hits_strip = FMTRec__Hits->getInt("strip",row);
	FMTRec__Hits_fitResidual = FMTRec__Hits->getFloat("fitResidual",row);
	FMTRec__Hits_trkingStat = FMTRec__Hits->getInt("trkingStat",row);
	FMTRec__Hits_clusterID = FMTRec__Hits->getShort("clusterID",row);
	FMTRec__Hits_trkID = FMTRec__Hits->getShort("trkID",row);
	return 0;
} 

int TIdentificatorCLAS12::get_BAND__tdc(int row){
	BAND__tdc_sector = BAND__tdc->getByte("sector",row);
	BAND__tdc_layer = BAND__tdc->getByte("layer",row);
	BAND__tdc_component = BAND__tdc->getShort("component",row);
	BAND__tdc_order = BAND__tdc->getByte("order",row);
	BAND__tdc_TDC = BAND__tdc->getInt("TDC",row);
	return 0;
} 

int TIdentificatorCLAS12::get_RF__adc(int row){
	RF__adc_sector = RF__adc->getByte("sector",row);
	RF__adc_layer = RF__adc->getByte("layer",row);
	RF__adc_component = RF__adc->getShort("component",row);
	RF__adc_order = RF__adc->getByte("order",row);
	RF__adc_ADC = RF__adc->getInt("ADC",row);
	RF__adc_time = RF__adc->getFloat("time",row);
	RF__adc_ped = RF__adc->getShort("ped",row);
	return 0;
} 

int TIdentificatorCLAS12::get_FMTRec__Crosses(int row){
	FMTRec__Crosses_ID = FMTRec__Crosses->getShort("ID",row);
	FMTRec__Crosses_sector = FMTRec__Crosses->getByte("sector",row);
	FMTRec__Crosses_region = FMTRec__Crosses->getByte("region",row);
	FMTRec__Crosses_x = FMTRec__Crosses->getFloat("x",row);
	FMTRec__Crosses_y = FMTRec__Crosses->getFloat("y",row);
	FMTRec__Crosses_z = FMTRec__Crosses->getFloat("z",row);
	FMTRec__Crosses_err_x = FMTRec__Crosses->getFloat("err_x",row);
	FMTRec__Crosses_err_y = FMTRec__Crosses->getFloat("err_y",row);
	FMTRec__Crosses_err_z = FMTRec__Crosses->getFloat("err_z",row);
	FMTRec__Crosses_ux = FMTRec__Crosses->getFloat("ux",row);
	FMTRec__Crosses_uy = FMTRec__Crosses->getFloat("uy",row);
	FMTRec__Crosses_uz = FMTRec__Crosses->getFloat("uz",row);
	FMTRec__Crosses_Cluster1_ID = FMTRec__Crosses->getShort("Cluster1_ID",row);
	FMTRec__Crosses_Cluster2_ID = FMTRec__Crosses->getShort("Cluster2_ID",row);
	FMTRec__Crosses_trkID = FMTRec__Crosses->getShort("trkID",row);
	return 0;
} 

int TIdentificatorCLAS12::get_MC__Particle(int row){
	MC__Particle_pid = MC__Particle->getInt("pid",row);
	MC__Particle_px = MC__Particle->getFloat("px",row);
	MC__Particle_py = MC__Particle->getFloat("py",row);
	MC__Particle_pz = MC__Particle->getFloat("pz",row);
	MC__Particle_vx = MC__Particle->getFloat("vx",row);
	MC__Particle_vy = MC__Particle->getFloat("vy",row);
	MC__Particle_vz = MC__Particle->getFloat("vz",row);
	MC__Particle_vt = MC__Particle->getFloat("vt",row);
	return 0;
} 

int TIdentificatorCLAS12::get_ECAL__tdc(int row){
	ECAL__tdc_sector = ECAL__tdc->getByte("sector",row);
	ECAL__tdc_layer = ECAL__tdc->getByte("layer",row);
	ECAL__tdc_component = ECAL__tdc->getShort("component",row);
	ECAL__tdc_order = ECAL__tdc->getByte("order",row);
	ECAL__tdc_TDC = ECAL__tdc->getInt("TDC",row);
	return 0;
} 

int TIdentificatorCLAS12::get_RICH__tdc(int row){
	RICH__tdc_sector = RICH__tdc->getByte("sector",row);
	RICH__tdc_layer = RICH__tdc->getByte("layer",row);
	RICH__tdc_component = RICH__tdc->getShort("component",row);
	RICH__tdc_order = RICH__tdc->getByte("order",row);
	RICH__tdc_TDC = RICH__tdc->getInt("TDC",row);
	return 0;
} 

int TIdentificatorCLAS12::get_LTCC__tdc(int row){
	LTCC__tdc_sector = LTCC__tdc->getByte("sector",row);
	LTCC__tdc_layer = LTCC__tdc->getByte("layer",row);
	LTCC__tdc_component = LTCC__tdc->getShort("component",row);
	LTCC__tdc_order = LTCC__tdc->getByte("order",row);
	LTCC__tdc_TDC = LTCC__tdc->getInt("TDC",row);
	return 0;
} 

int TIdentificatorCLAS12::get_ECAL__moments(int row){
	ECAL__moments_distU = ECAL__moments->getFloat("distU",row);
	ECAL__moments_distV = ECAL__moments->getFloat("distV",row);
	ECAL__moments_distW = ECAL__moments->getFloat("distW",row);
	ECAL__moments_m1u = ECAL__moments->getFloat("m1u",row);
	ECAL__moments_m1v = ECAL__moments->getFloat("m1v",row);
	ECAL__moments_m1w = ECAL__moments->getFloat("m1w",row);
	ECAL__moments_m2u = ECAL__moments->getFloat("m2u",row);
	ECAL__moments_m2v = ECAL__moments->getFloat("m2v",row);
	ECAL__moments_m2w = ECAL__moments->getFloat("m2w",row);
	ECAL__moments_m3u = ECAL__moments->getFloat("m3u",row);
	ECAL__moments_m3v = ECAL__moments->getFloat("m3v",row);
	ECAL__moments_m3w = ECAL__moments->getFloat("m3w",row);
	return 0;
} 

int TIdentificatorCLAS12::get_EVENT__particle(int row){
	EVENT__particle_status = EVENT__particle->getByte("status",row);
	EVENT__particle_charge = EVENT__particle->getByte("charge",row);
	EVENT__particle_pid = EVENT__particle->getInt("pid",row);
	EVENT__particle_mass = EVENT__particle->getFloat("mass",row);
	EVENT__particle_px = EVENT__particle->getFloat("px",row);
	EVENT__particle_py = EVENT__particle->getFloat("py",row);
	EVENT__particle_pz = EVENT__particle->getFloat("pz",row);
	EVENT__particle_vx = EVENT__particle->getFloat("vx",row);
	EVENT__particle_vy = EVENT__particle->getFloat("vy",row);
	EVENT__particle_vz = EVENT__particle->getFloat("vz",row);
	EVENT__particle_dcstat = EVENT__particle->getByte("dcstat",row);
	EVENT__particle_ecstat = EVENT__particle->getByte("ecstat",row);
	EVENT__particle_scstat = EVENT__particle->getByte("scstat",row);
	EVENT__particle_ccstat = EVENT__particle->getByte("ccstat",row);
	EVENT__particle_lcstat = EVENT__particle->getByte("lcstat",row);
	return 0;
} 

int TIdentificatorCLAS12::get_RUN__trigger(int row){
	RUN__trigger_id = RUN__trigger->getInt("id",row);
	RUN__trigger_trigger = RUN__trigger->getInt("trigger",row);
	return 0;
} 

int TIdentificatorCLAS12::get_HitBasedTrkg__HBHits(int row){
	HitBasedTrkg__HBHits_id = HitBasedTrkg__HBHits->getShort("id",row);
	HitBasedTrkg__HBHits_status = HitBasedTrkg__HBHits->getShort("status",row);
	HitBasedTrkg__HBHits_sector = HitBasedTrkg__HBHits->getByte("sector",row);
	HitBasedTrkg__HBHits_superlayer = HitBasedTrkg__HBHits->getByte("superlayer",row);
	HitBasedTrkg__HBHits_layer = HitBasedTrkg__HBHits->getByte("layer",row);
	HitBasedTrkg__HBHits_wire = HitBasedTrkg__HBHits->getShort("wire",row);
	HitBasedTrkg__HBHits_TDC = HitBasedTrkg__HBHits->getInt("TDC",row);
	HitBasedTrkg__HBHits_trkDoca = HitBasedTrkg__HBHits->getFloat("trkDoca",row);
	HitBasedTrkg__HBHits_docaError = HitBasedTrkg__HBHits->getFloat("docaError",row);
	HitBasedTrkg__HBHits_LR = HitBasedTrkg__HBHits->getByte("LR",row);
	HitBasedTrkg__HBHits_LocX = HitBasedTrkg__HBHits->getFloat("LocX",row);
	HitBasedTrkg__HBHits_LocY = HitBasedTrkg__HBHits->getFloat("LocY",row);
	HitBasedTrkg__HBHits_X = HitBasedTrkg__HBHits->getFloat("X",row);
	HitBasedTrkg__HBHits_Z = HitBasedTrkg__HBHits->getFloat("Z",row);
	HitBasedTrkg__HBHits_B = HitBasedTrkg__HBHits->getFloat("B",row);
	HitBasedTrkg__HBHits_TProp = HitBasedTrkg__HBHits->getFloat("TProp",row);
	HitBasedTrkg__HBHits_TFlight = HitBasedTrkg__HBHits->getFloat("TFlight",row);
	HitBasedTrkg__HBHits_clusterID = HitBasedTrkg__HBHits->getShort("clusterID",row);
	HitBasedTrkg__HBHits_trkID = HitBasedTrkg__HBHits->getByte("trkID",row);
	return 0;
} 

int TIdentificatorCLAS12::get_TimeBasedTrkg__TBHits(int row){
	TimeBasedTrkg__TBHits_id = TimeBasedTrkg__TBHits->getShort("id",row);
	TimeBasedTrkg__TBHits_status = TimeBasedTrkg__TBHits->getShort("status",row);
	TimeBasedTrkg__TBHits_sector = TimeBasedTrkg__TBHits->getByte("sector",row);
	TimeBasedTrkg__TBHits_superlayer = TimeBasedTrkg__TBHits->getByte("superlayer",row);
	TimeBasedTrkg__TBHits_layer = TimeBasedTrkg__TBHits->getByte("layer",row);
	TimeBasedTrkg__TBHits_wire = TimeBasedTrkg__TBHits->getShort("wire",row);
	TimeBasedTrkg__TBHits_TDC = TimeBasedTrkg__TBHits->getInt("TDC",row);
	TimeBasedTrkg__TBHits_doca = TimeBasedTrkg__TBHits->getFloat("doca",row);
	TimeBasedTrkg__TBHits_docaError = TimeBasedTrkg__TBHits->getFloat("docaError",row);
	TimeBasedTrkg__TBHits_trkDoca = TimeBasedTrkg__TBHits->getFloat("trkDoca",row);
	TimeBasedTrkg__TBHits_timeResidual = TimeBasedTrkg__TBHits->getFloat("timeResidual",row);
	TimeBasedTrkg__TBHits_fitResidual = TimeBasedTrkg__TBHits->getFloat("fitResidual",row);
	TimeBasedTrkg__TBHits_LR = TimeBasedTrkg__TBHits->getByte("LR",row);
	TimeBasedTrkg__TBHits_X = TimeBasedTrkg__TBHits->getFloat("X",row);
	TimeBasedTrkg__TBHits_Z = TimeBasedTrkg__TBHits->getFloat("Z",row);
	TimeBasedTrkg__TBHits_B = TimeBasedTrkg__TBHits->getFloat("B",row);
	TimeBasedTrkg__TBHits_Alpha = TimeBasedTrkg__TBHits->getFloat("Alpha",row);
	TimeBasedTrkg__TBHits_TProp = TimeBasedTrkg__TBHits->getFloat("TProp",row);
	TimeBasedTrkg__TBHits_TFlight = TimeBasedTrkg__TBHits->getFloat("TFlight",row);
	TimeBasedTrkg__TBHits_T0 = TimeBasedTrkg__TBHits->getFloat("T0",row);
	TimeBasedTrkg__TBHits_TStart = TimeBasedTrkg__TBHits->getFloat("TStart",row);
	TimeBasedTrkg__TBHits_clusterID = TimeBasedTrkg__TBHits->getShort("clusterID",row);
	TimeBasedTrkg__TBHits_trkID = TimeBasedTrkg__TBHits->getByte("trkID",row);
	TimeBasedTrkg__TBHits_time = TimeBasedTrkg__TBHits->getFloat("time",row);
	TimeBasedTrkg__TBHits_beta = TimeBasedTrkg__TBHits->getFloat("beta",row);
	TimeBasedTrkg__TBHits_tBeta = TimeBasedTrkg__TBHits->getFloat("tBeta",row);
	return 0;
} 

int TIdentificatorCLAS12::get_CTOF__hits(int row){
	CTOF__hits_id = CTOF__hits->getShort("id",row);
	CTOF__hits_status = CTOF__hits->getShort("status",row);
	CTOF__hits_trkID = CTOF__hits->getShort("trkID",row);
	CTOF__hits_sector = CTOF__hits->getByte("sector",row);
	CTOF__hits_layer = CTOF__hits->getByte("layer",row);
	CTOF__hits_component = CTOF__hits->getShort("component",row);
	CTOF__hits_energy = CTOF__hits->getFloat("energy",row);
	CTOF__hits_time = CTOF__hits->getFloat("time",row);
	CTOF__hits_energy_unc = CTOF__hits->getFloat("energy_unc",row);
	CTOF__hits_time_unc = CTOF__hits->getFloat("time_unc",row);
	CTOF__hits_x = CTOF__hits->getFloat("x",row);
	CTOF__hits_y = CTOF__hits->getFloat("y",row);
	CTOF__hits_z = CTOF__hits->getFloat("z",row);
	CTOF__hits_x_unc = CTOF__hits->getFloat("x_unc",row);
	CTOF__hits_y_unc = CTOF__hits->getFloat("y_unc",row);
	CTOF__hits_z_unc = CTOF__hits->getFloat("z_unc",row);
	CTOF__hits_tx = CTOF__hits->getFloat("tx",row);
	CTOF__hits_ty = CTOF__hits->getFloat("ty",row);
	CTOF__hits_tz = CTOF__hits->getFloat("tz",row);
	CTOF__hits_adc_idx1 = CTOF__hits->getShort("adc_idx1",row);
	CTOF__hits_adc_idx2 = CTOF__hits->getShort("adc_idx2",row);
	CTOF__hits_tdc_idx1 = CTOF__hits->getShort("tdc_idx1",row);
	CTOF__hits_tdc_idx2 = CTOF__hits->getShort("tdc_idx2",row);
	CTOF__hits_pathLength = CTOF__hits->getFloat("pathLength",row);
	CTOF__hits_pathLengthThruBar = CTOF__hits->getFloat("pathLengthThruBar",row);
	return 0;
} 

int TIdentificatorCLAS12::get_RECHB__Particle(int row){
	RECHB__Particle_pid = RECHB__Particle->getInt("pid",row);
	RECHB__Particle_px = RECHB__Particle->getFloat("px",row);
	RECHB__Particle_py = RECHB__Particle->getFloat("py",row);
	RECHB__Particle_pz = RECHB__Particle->getFloat("pz",row);
	RECHB__Particle_vx = RECHB__Particle->getFloat("vx",row);
	RECHB__Particle_vy = RECHB__Particle->getFloat("vy",row);
	RECHB__Particle_vz = RECHB__Particle->getFloat("vz",row);
	RECHB__Particle_charge = RECHB__Particle->getByte("charge",row);
	RECHB__Particle_beta = RECHB__Particle->getFloat("beta",row);
	RECHB__Particle_chi2pid = RECHB__Particle->getFloat("chi2pid",row);
	RECHB__Particle_status = RECHB__Particle->getShort("status",row);
	return 0;
} 

int TIdentificatorCLAS12::get_RAW__tdc(int row){
	RAW__tdc_crate = RAW__tdc->getByte("crate",row);
	RAW__tdc_slot = RAW__tdc->getByte("slot",row);
	RAW__tdc_channel = RAW__tdc->getShort("channel",row);
	RAW__tdc_order = RAW__tdc->getByte("order",row);
	RAW__tdc_TDC = RAW__tdc->getInt("TDC",row);
	return 0;
} 

int TIdentificatorCLAS12::get_RICH__photons(int row){
	RICH__photons_id = RICH__photons->getShort("id",row);
	RICH__photons_type = RICH__photons->getShort("type",row);
	RICH__photons_hit_index = RICH__photons->getShort("hit_index",row);
	RICH__photons_hadron_index = RICH__photons->getShort("hadron_index",row);
	RICH__photons_start_time = RICH__photons->getFloat("start_time",row);
	RICH__photons_analytic_the = RICH__photons->getFloat("analytic_the",row);
	RICH__photons_analytic_phi = RICH__photons->getFloat("analytic_phi",row);
	RICH__photons_analytic_path = RICH__photons->getFloat("analytic_path",row);
	RICH__photons_analytic_time = RICH__photons->getFloat("analytic_time",row);
	RICH__photons_analytic_EtaC = RICH__photons->getFloat("analytic_EtaC",row);
	RICH__photons_analytic_aeron = RICH__photons->getFloat("analytic_aeron",row);
	RICH__photons_analytic_elpr = RICH__photons->getFloat("analytic_elpr",row);
	RICH__photons_analytic_pipr = RICH__photons->getFloat("analytic_pipr",row);
	RICH__photons_analytic_kpr = RICH__photons->getFloat("analytic_kpr",row);
	RICH__photons_analytic_prpr = RICH__photons->getFloat("analytic_prpr",row);
	RICH__photons_analytic_bgpr = RICH__photons->getFloat("analytic_bgpr",row);
	RICH__photons_traced_the = RICH__photons->getFloat("traced_the",row);
	RICH__photons_traced_phi = RICH__photons->getFloat("traced_phi",row);
	RICH__photons_traced_hitx = RICH__photons->getFloat("traced_hitx",row);
	RICH__photons_traced_hity = RICH__photons->getFloat("traced_hity",row);
	RICH__photons_traced_hitz = RICH__photons->getFloat("traced_hitz",row);
	RICH__photons_traced_path = RICH__photons->getFloat("traced_path",row);
	RICH__photons_traced_time = RICH__photons->getFloat("traced_time",row);
	RICH__photons_traced_nrfl = RICH__photons->getShort("traced_nrfl",row);
	RICH__photons_traced_nrfr = RICH__photons->getShort("traced_nrfr",row);
	RICH__photons_traced_1rfl = RICH__photons->getShort("traced_1rfl",row);
	RICH__photons_traced_EtaC = RICH__photons->getFloat("traced_EtaC",row);
	RICH__photons_traced_aeron = RICH__photons->getFloat("traced_aeron",row);
	RICH__photons_traced_elpr = RICH__photons->getFloat("traced_elpr",row);
	RICH__photons_traced_pipr = RICH__photons->getFloat("traced_pipr",row);
	RICH__photons_traced_kpr = RICH__photons->getFloat("traced_kpr",row);
	RICH__photons_traced_prpr = RICH__photons->getFloat("traced_prpr",row);
	RICH__photons_traced_bgpr = RICH__photons->getFloat("traced_bgpr",row);
	return 0;
} 

int TIdentificatorCLAS12::get_HEL__flip(int row){
	HEL__flip_run = HEL__flip->getInt("run",row);
	HEL__flip_event = HEL__flip->getInt("event",row);
	HEL__flip_timestamp = HEL__flip->getLong("timestamp",row);
	HEL__flip_helicity = HEL__flip->getByte("helicity",row);
	HEL__flip_helicityRaw = HEL__flip->getByte("helicityRaw",row);
	HEL__flip_pair = HEL__flip->getByte("pair",row);
	HEL__flip_pattern = HEL__flip->getByte("pattern",row);
	HEL__flip_status = HEL__flip->getByte("status",row);
	return 0;
} 

int TIdentificatorCLAS12::get_TimeBasedTrkg__Trajectory(int row){
	TimeBasedTrkg__Trajectory_id = TimeBasedTrkg__Trajectory->getShort("id",row);
	TimeBasedTrkg__Trajectory_detector = TimeBasedTrkg__Trajectory->getShort("detector",row);
	TimeBasedTrkg__Trajectory_layer = TimeBasedTrkg__Trajectory->getByte("layer",row);
	TimeBasedTrkg__Trajectory_x = TimeBasedTrkg__Trajectory->getFloat("x",row);
	TimeBasedTrkg__Trajectory_y = TimeBasedTrkg__Trajectory->getFloat("y",row);
	TimeBasedTrkg__Trajectory_z = TimeBasedTrkg__Trajectory->getFloat("z",row);
	TimeBasedTrkg__Trajectory_tx = TimeBasedTrkg__Trajectory->getFloat("tx",row);
	TimeBasedTrkg__Trajectory_ty = TimeBasedTrkg__Trajectory->getFloat("ty",row);
	TimeBasedTrkg__Trajectory_tz = TimeBasedTrkg__Trajectory->getFloat("tz",row);
	TimeBasedTrkg__Trajectory_B = TimeBasedTrkg__Trajectory->getFloat("B",row);
	TimeBasedTrkg__Trajectory_path = TimeBasedTrkg__Trajectory->getFloat("path",row);
	return 0;
} 

int TIdentificatorCLAS12::get_DETECTOR__ecpb(int row){
	DETECTOR__ecpb_sector = DETECTOR__ecpb->getByte("sector",row);
	DETECTOR__ecpb_etot = DETECTOR__ecpb->getFloat("etot",row);
	DETECTOR__ecpb_ein = DETECTOR__ecpb->getFloat("ein",row);
	DETECTOR__ecpb_eout = DETECTOR__ecpb->getFloat("eout",row);
	DETECTOR__ecpb_time = DETECTOR__ecpb->getFloat("time",row);
	DETECTOR__ecpb_path = DETECTOR__ecpb->getFloat("path",row);
	DETECTOR__ecpb_x = DETECTOR__ecpb->getFloat("x",row);
	DETECTOR__ecpb_y = DETECTOR__ecpb->getFloat("y",row);
	DETECTOR__ecpb_z = DETECTOR__ecpb->getFloat("z",row);
	return 0;
} 

int TIdentificatorCLAS12::get_CND__hits(int row){
	CND__hits_id = CND__hits->getShort("id",row);
	CND__hits_status = CND__hits->getShort("status",row);
	CND__hits_trkID = CND__hits->getShort("trkID",row);
	CND__hits_sector = CND__hits->getByte("sector",row);
	CND__hits_layer = CND__hits->getByte("layer",row);
	CND__hits_component = CND__hits->getShort("component",row);
	CND__hits_energy = CND__hits->getFloat("energy",row);
	CND__hits_time = CND__hits->getFloat("time",row);
	CND__hits_energy_unc = CND__hits->getFloat("energy_unc",row);
	CND__hits_time_unc = CND__hits->getFloat("time_unc",row);
	CND__hits_x = CND__hits->getFloat("x",row);
	CND__hits_y = CND__hits->getFloat("y",row);
	CND__hits_z = CND__hits->getFloat("z",row);
	CND__hits_x_unc = CND__hits->getFloat("x_unc",row);
	CND__hits_y_unc = CND__hits->getFloat("y_unc",row);
	CND__hits_z_unc = CND__hits->getFloat("z_unc",row);
	CND__hits_tx = CND__hits->getFloat("tx",row);
	CND__hits_ty = CND__hits->getFloat("ty",row);
	CND__hits_tz = CND__hits->getFloat("tz",row);
	CND__hits_tlength = CND__hits->getFloat("tlength",row);
	CND__hits_pathlength = CND__hits->getFloat("pathlength",row);
	CND__hits_indexLadc = CND__hits->getShort("indexLadc",row);
	CND__hits_indexRadc = CND__hits->getShort("indexRadc",row);
	CND__hits_indexLtdc = CND__hits->getShort("indexLtdc",row);
	CND__hits_indexRtdc = CND__hits->getShort("indexRtdc",row);
	return 0;
} 

int TIdentificatorCLAS12::InitBanks(){
	if (fFactory->hasSchema("BMTRec::Hits"))
		BMTRec__Hits = new hipo::bank(fFactory->getSchema("BMTRec::Hits"));
	if (fFactory->hasSchema("RAW::adc"))
		RAW__adc = new hipo::bank(fFactory->getSchema("RAW::adc"));
	if (fFactory->hasSchema("BAND::adc"))
		BAND__adc = new hipo::bank(fFactory->getSchema("BAND::adc"));
	if (fFactory->hasSchema("RUN::config"))
		RUN__config = new hipo::bank(fFactory->getSchema("RUN::config"));
	if (fFactory->hasSchema("RICH::clusters"))
		RICH__clusters = new hipo::bank(fFactory->getSchema("RICH::clusters"));
	if (fFactory->hasSchema("RECHB::Scintillator"))
		RECHB__Scintillator = new hipo::bank(fFactory->getSchema("RECHB::Scintillator"));
	if (fFactory->hasSchema("REC::RingCher"))
		REC__RingCher = new hipo::bank(fFactory->getSchema("REC::RingCher"));
	if (fFactory->hasSchema("BSTRec::LayerEffs"))
		BSTRec__LayerEffs = new hipo::bank(fFactory->getSchema("BSTRec::LayerEffs"));
	if (fFactory->hasSchema("RTPC::pos"))
		RTPC__pos = new hipo::bank(fFactory->getSchema("RTPC::pos"));
	if (fFactory->hasSchema("TimeBasedTrkg::TBCrosses"))
		TimeBasedTrkg__TBCrosses = new hipo::bank(fFactory->getSchema("TimeBasedTrkg::TBCrosses"));
	if (fFactory->hasSchema("HitBasedTrkg::HBTracks"))
		HitBasedTrkg__HBTracks = new hipo::bank(fFactory->getSchema("HitBasedTrkg::HBTracks"));
	if (fFactory->hasSchema("CVTRec::Cosmics"))
		CVTRec__Cosmics = new hipo::bank(fFactory->getSchema("CVTRec::Cosmics"));
	if (fFactory->hasSchema("RECHB::Cherenkov"))
		RECHB__Cherenkov = new hipo::bank(fFactory->getSchema("RECHB::Cherenkov"));
	if (fFactory->hasSchema("BSTRec::Clusters"))
		BSTRec__Clusters = new hipo::bank(fFactory->getSchema("BSTRec::Clusters"));
	if (fFactory->hasSchema("CVTRec::Trajectory"))
		CVTRec__Trajectory = new hipo::bank(fFactory->getSchema("CVTRec::Trajectory"));
	if (fFactory->hasSchema("RECHB::Calorimeter"))
		RECHB__Calorimeter = new hipo::bank(fFactory->getSchema("RECHB::Calorimeter"));
	if (fFactory->hasSchema("TimeBasedTrkg::TBSegmentTrajectory"))
		TimeBasedTrkg__TBSegmentTrajectory = new hipo::bank(fFactory->getSchema("TimeBasedTrkg::TBSegmentTrajectory"));
	if (fFactory->hasSchema("CTOF::tdc"))
		CTOF__tdc = new hipo::bank(fFactory->getSchema("CTOF::tdc"));
	if (fFactory->hasSchema("REC::Cherenkov"))
		REC__Cherenkov = new hipo::bank(fFactory->getSchema("REC::Cherenkov"));
	if (fFactory->hasSchema("BMTRec::LayerEffs"))
		BMTRec__LayerEffs = new hipo::bank(fFactory->getSchema("BMTRec::LayerEffs"));
	if (fFactory->hasSchema("FTOF::adc"))
		FTOF__adc = new hipo::bank(fFactory->getSchema("FTOF::adc"));
	if (fFactory->hasSchema("MC::Lund"))
		MC__Lund = new hipo::bank(fFactory->getSchema("MC::Lund"));
	if (fFactory->hasSchema("DETECTOR::lcpb"))
		DETECTOR__lcpb = new hipo::bank(fFactory->getSchema("DETECTOR::lcpb"));
	if (fFactory->hasSchema("MC::Header"))
		MC__Header = new hipo::bank(fFactory->getSchema("MC::Header"));
	if (fFactory->hasSchema("CND::clusters"))
		CND__clusters = new hipo::bank(fFactory->getSchema("CND::clusters"));
	if (fFactory->hasSchema("TimeBasedTrkg::TBCovMat"))
		TimeBasedTrkg__TBCovMat = new hipo::bank(fFactory->getSchema("TimeBasedTrkg::TBCovMat"));
	if (fFactory->hasSchema("RICH::hits"))
		RICH__hits = new hipo::bank(fFactory->getSchema("RICH::hits"));
	if (fFactory->hasSchema("RECHB::Track"))
		RECHB__Track = new hipo::bank(fFactory->getSchema("RECHB::Track"));
	if (fFactory->hasSchema("MC::True"))
		MC__True = new hipo::bank(fFactory->getSchema("MC::True"));
	if (fFactory->hasSchema("BST::adc"))
		BST__adc = new hipo::bank(fFactory->getSchema("BST::adc"));
	if (fFactory->hasSchema("MC::Event"))
		MC__Event = new hipo::bank(fFactory->getSchema("MC::Event"));
	if (fFactory->hasSchema("HitBasedTrkg::HBCrosses"))
		HitBasedTrkg__HBCrosses = new hipo::bank(fFactory->getSchema("HitBasedTrkg::HBCrosses"));
	if (fFactory->hasSchema("FTOF::clusters"))
		FTOF__clusters = new hipo::bank(fFactory->getSchema("FTOF::clusters"));
	if (fFactory->hasSchema("REC::TrackCross"))
		REC__TrackCross = new hipo::bank(fFactory->getSchema("REC::TrackCross"));
	if (fFactory->hasSchema("REC::Scintillator"))
		REC__Scintillator = new hipo::bank(fFactory->getSchema("REC::Scintillator"));
	if (fFactory->hasSchema("ECAL::peaks"))
		ECAL__peaks = new hipo::bank(fFactory->getSchema("ECAL::peaks"));
	if (fFactory->hasSchema("TimeBasedTrkg::TBClusters"))
		TimeBasedTrkg__TBClusters = new hipo::bank(fFactory->getSchema("TimeBasedTrkg::TBClusters"));
	if (fFactory->hasSchema("TAGGER::tgpb"))
		TAGGER__tgpb = new hipo::bank(fFactory->getSchema("TAGGER::tgpb"));
	if (fFactory->hasSchema("LTCC::clusters"))
		LTCC__clusters = new hipo::bank(fFactory->getSchema("LTCC::clusters"));
	if (fFactory->hasSchema("RICH::hadrons"))
		RICH__hadrons = new hipo::bank(fFactory->getSchema("RICH::hadrons"));
	if (fFactory->hasSchema("DC::tdc"))
		DC__tdc = new hipo::bank(fFactory->getSchema("DC::tdc"));
	if (fFactory->hasSchema("BSTRec::Crosses"))
		BSTRec__Crosses = new hipo::bank(fFactory->getSchema("BSTRec::Crosses"));
	if (fFactory->hasSchema("HTCC::rec"))
		HTCC__rec = new hipo::bank(fFactory->getSchema("HTCC::rec"));
	if (fFactory->hasSchema("FTHODO::adc"))
		FTHODO__adc = new hipo::bank(fFactory->getSchema("FTHODO::adc"));
	if (fFactory->hasSchema("FTOF::matchedclusters"))
		FTOF__matchedclusters = new hipo::bank(fFactory->getSchema("FTOF::matchedclusters"));
	if (fFactory->hasSchema("ECAL::adc"))
		ECAL__adc = new hipo::bank(fFactory->getSchema("ECAL::adc"));
	if (fFactory->hasSchema("FMT::adc"))
		FMT__adc = new hipo::bank(fFactory->getSchema("FMT::adc"));
	if (fFactory->hasSchema("BMT::adc"))
		BMT__adc = new hipo::bank(fFactory->getSchema("BMT::adc"));
	if (fFactory->hasSchema("ECAL::calib"))
		ECAL__calib = new hipo::bank(fFactory->getSchema("ECAL::calib"));
	if (fFactory->hasSchema("BMTRec::Crosses"))
		BMTRec__Crosses = new hipo::bank(fFactory->getSchema("BMTRec::Crosses"));
	if (fFactory->hasSchema("RUN::rf"))
		RUN__rf = new hipo::bank(fFactory->getSchema("RUN::rf"));
	if (fFactory->hasSchema("BSTRec::Hits"))
		BSTRec__Hits = new hipo::bank(fFactory->getSchema("BSTRec::Hits"));
	if (fFactory->hasSchema("HitBasedTrkg::HBSegmentTrajectory"))
		HitBasedTrkg__HBSegmentTrajectory = new hipo::bank(fFactory->getSchema("HitBasedTrkg::HBSegmentTrajectory"));
	if (fFactory->hasSchema("LTCC::adc"))
		LTCC__adc = new hipo::bank(fFactory->getSchema("LTCC::adc"));
	if (fFactory->hasSchema("CND::tdc"))
		CND__tdc = new hipo::bank(fFactory->getSchema("CND::tdc"));
	if (fFactory->hasSchema("HTCC::adc"))
		HTCC__adc = new hipo::bank(fFactory->getSchema("HTCC::adc"));
	if (fFactory->hasSchema("RAW::vtp"))
		RAW__vtp = new hipo::bank(fFactory->getSchema("RAW::vtp"));
	if (fFactory->hasSchema("CVTRec::Tracks"))
		CVTRec__Tracks = new hipo::bank(fFactory->getSchema("CVTRec::Tracks"));
	if (fFactory->hasSchema("RECFT::Particle"))
		RECFT__Particle = new hipo::bank(fFactory->getSchema("RECFT::Particle"));
	if (fFactory->hasSchema("TimeBasedTrkg::TBTracks"))
		TimeBasedTrkg__TBTracks = new hipo::bank(fFactory->getSchema("TimeBasedTrkg::TBTracks"));
	if (fFactory->hasSchema("DETECTOR::icpb"))
		DETECTOR__icpb = new hipo::bank(fFactory->getSchema("DETECTOR::icpb"));
	if (fFactory->hasSchema("RF::tdc"))
		RF__tdc = new hipo::bank(fFactory->getSchema("RF::tdc"));
	if (fFactory->hasSchema("HEADER::info"))
		HEADER__info = new hipo::bank(fFactory->getSchema("HEADER::info"));
	if (fFactory->hasSchema("FTCAL::clusters"))
		FTCAL__clusters = new hipo::bank(fFactory->getSchema("FTCAL::clusters"));
	if (fFactory->hasSchema("RAW::scaler"))
		RAW__scaler = new hipo::bank(fFactory->getSchema("RAW::scaler"));
	if (fFactory->hasSchema("BMTRec::Clusters"))
		BMTRec__Clusters = new hipo::bank(fFactory->getSchema("BMTRec::Clusters"));
	if (fFactory->hasSchema("TimeBasedTrkg::TBSegments"))
		TimeBasedTrkg__TBSegments = new hipo::bank(fFactory->getSchema("TimeBasedTrkg::TBSegments"));
	if (fFactory->hasSchema("RECHB::ForwardTagger"))
		RECHB__ForwardTagger = new hipo::bank(fFactory->getSchema("RECHB::ForwardTagger"));
	if (fFactory->hasSchema("REC::Calorimeter"))
		REC__Calorimeter = new hipo::bank(fFactory->getSchema("REC::Calorimeter"));
	if (fFactory->hasSchema("REC::CovMat"))
		REC__CovMat = new hipo::bank(fFactory->getSchema("REC::CovMat"));
	if (fFactory->hasSchema("RAW::epics"))
		RAW__epics = new hipo::bank(fFactory->getSchema("RAW::epics"));
	if (fFactory->hasSchema("REC::VertDoca"))
		REC__VertDoca = new hipo::bank(fFactory->getSchema("REC::VertDoca"));
	if (fFactory->hasSchema("FTHODO::clusters"))
		FTHODO__clusters = new hipo::bank(fFactory->getSchema("FTHODO::clusters"));
	if (fFactory->hasSchema("RUN::scaler"))
		RUN__scaler = new hipo::bank(fFactory->getSchema("RUN::scaler"));
	if (fFactory->hasSchema("BAND::hits"))
		BAND__hits = new hipo::bank(fFactory->getSchema("BAND::hits"));
	if (fFactory->hasSchema("ECAL::hits"))
		ECAL__hits = new hipo::bank(fFactory->getSchema("ECAL::hits"));
	if (fFactory->hasSchema("ECAL::clusters"))
		ECAL__clusters = new hipo::bank(fFactory->getSchema("ECAL::clusters"));
	if (fFactory->hasSchema("FT::particles"))
		FT__particles = new hipo::bank(fFactory->getSchema("FT::particles"));
	if (fFactory->hasSchema("HitBasedTrkg::HBClusters"))
		HitBasedTrkg__HBClusters = new hipo::bank(fFactory->getSchema("HitBasedTrkg::HBClusters"));
	if (fFactory->hasSchema("REC::RICH"))
		REC__RICH = new hipo::bank(fFactory->getSchema("REC::RICH"));
	if (fFactory->hasSchema("REC::Track"))
		REC__Track = new hipo::bank(fFactory->getSchema("REC::Track"));
	if (fFactory->hasSchema("RECHB::Event"))
		RECHB__Event = new hipo::bank(fFactory->getSchema("RECHB::Event"));
	if (fFactory->hasSchema("FTOF::hbhits"))
		FTOF__hbhits = new hipo::bank(fFactory->getSchema("FTOF::hbhits"));
	if (fFactory->hasSchema("RTPC::adc"))
		RTPC__adc = new hipo::bank(fFactory->getSchema("RTPC::adc"));
	if (fFactory->hasSchema("FTCAL::adc"))
		FTCAL__adc = new hipo::bank(fFactory->getSchema("FTCAL::adc"));
	if (fFactory->hasSchema("CTOF::rawhits"))
		CTOF__rawhits = new hipo::bank(fFactory->getSchema("CTOF::rawhits"));
	if (fFactory->hasSchema("CND::adc"))
		CND__adc = new hipo::bank(fFactory->getSchema("CND::adc"));
	if (fFactory->hasSchema("FTOF::hits"))
		FTOF__hits = new hipo::bank(fFactory->getSchema("FTOF::hits"));
	if (fFactory->hasSchema("FMTRec::Clusters"))
		FMTRec__Clusters = new hipo::bank(fFactory->getSchema("FMTRec::Clusters"));
	if (fFactory->hasSchema("HEL::adc"))
		HEL__adc = new hipo::bank(fFactory->getSchema("HEL::adc"));
	if (fFactory->hasSchema("DETECTOR::ccpb"))
		DETECTOR__ccpb = new hipo::bank(fFactory->getSchema("DETECTOR::ccpb"));
	if (fFactory->hasSchema("REC::ForwardTagger"))
		REC__ForwardTagger = new hipo::bank(fFactory->getSchema("REC::ForwardTagger"));
	if (fFactory->hasSchema("FTTRK::adc"))
		FTTRK__adc = new hipo::bank(fFactory->getSchema("FTTRK::adc"));
	if (fFactory->hasSchema("FTOF::rawhits"))
		FTOF__rawhits = new hipo::bank(fFactory->getSchema("FTOF::rawhits"));
	if (fFactory->hasSchema("RECHB::TrackCross"))
		RECHB__TrackCross = new hipo::bank(fFactory->getSchema("RECHB::TrackCross"));
	if (fFactory->hasSchema("DETECTOR::scpb"))
		DETECTOR__scpb = new hipo::bank(fFactory->getSchema("DETECTOR::scpb"));
	if (fFactory->hasSchema("FTHODO::hits"))
		FTHODO__hits = new hipo::bank(fFactory->getSchema("FTHODO::hits"));
	if (fFactory->hasSchema("RECFT::Event"))
		RECFT__Event = new hipo::bank(fFactory->getSchema("RECFT::Event"));
	if (fFactory->hasSchema("CTOF::adc"))
		CTOF__adc = new hipo::bank(fFactory->getSchema("CTOF::adc"));
	if (fFactory->hasSchema("EVENT::detector"))
		EVENT__detector = new hipo::bank(fFactory->getSchema("EVENT::detector"));
	if (fFactory->hasSchema("DC::doca"))
		DC__doca = new hipo::bank(fFactory->getSchema("DC::doca"));
	if (fFactory->hasSchema("HTCC::tdc"))
		HTCC__tdc = new hipo::bank(fFactory->getSchema("HTCC::tdc"));
	if (fFactory->hasSchema("HEL::online"))
		HEL__online = new hipo::bank(fFactory->getSchema("HEL::online"));
	if (fFactory->hasSchema("REC::Event"))
		REC__Event = new hipo::bank(fFactory->getSchema("REC::Event"));
	if (fFactory->hasSchema("REC::Particle"))
		REC__Particle = new hipo::bank(fFactory->getSchema("REC::Particle"));
	if (fFactory->hasSchema("REC::Traj"))
		REC__Traj = new hipo::bank(fFactory->getSchema("REC::Traj"));
	if (fFactory->hasSchema("FTOF::tdc"))
		FTOF__tdc = new hipo::bank(fFactory->getSchema("FTOF::tdc"));
	if (fFactory->hasSchema("HitBasedTrkg::HBSegments"))
		HitBasedTrkg__HBSegments = new hipo::bank(fFactory->getSchema("HitBasedTrkg::HBSegments"));
	if (fFactory->hasSchema("FTCAL::hits"))
		FTCAL__hits = new hipo::bank(fFactory->getSchema("FTCAL::hits"));
	if (fFactory->hasSchema("FMTRec::Hits"))
		FMTRec__Hits = new hipo::bank(fFactory->getSchema("FMTRec::Hits"));
	if (fFactory->hasSchema("BAND::tdc"))
		BAND__tdc = new hipo::bank(fFactory->getSchema("BAND::tdc"));
	if (fFactory->hasSchema("RF::adc"))
		RF__adc = new hipo::bank(fFactory->getSchema("RF::adc"));
	if (fFactory->hasSchema("FMTRec::Crosses"))
		FMTRec__Crosses = new hipo::bank(fFactory->getSchema("FMTRec::Crosses"));
	if (fFactory->hasSchema("MC::Particle"))
		MC__Particle = new hipo::bank(fFactory->getSchema("MC::Particle"));
	if (fFactory->hasSchema("ECAL::tdc"))
		ECAL__tdc = new hipo::bank(fFactory->getSchema("ECAL::tdc"));
	if (fFactory->hasSchema("RICH::tdc"))
		RICH__tdc = new hipo::bank(fFactory->getSchema("RICH::tdc"));
	if (fFactory->hasSchema("LTCC::tdc"))
		LTCC__tdc = new hipo::bank(fFactory->getSchema("LTCC::tdc"));
	if (fFactory->hasSchema("ECAL::moments"))
		ECAL__moments = new hipo::bank(fFactory->getSchema("ECAL::moments"));
	if (fFactory->hasSchema("EVENT::particle"))
		EVENT__particle = new hipo::bank(fFactory->getSchema("EVENT::particle"));
	if (fFactory->hasSchema("RUN::trigger"))
		RUN__trigger = new hipo::bank(fFactory->getSchema("RUN::trigger"));
	if (fFactory->hasSchema("HitBasedTrkg::HBHits"))
		HitBasedTrkg__HBHits = new hipo::bank(fFactory->getSchema("HitBasedTrkg::HBHits"));
	if (fFactory->hasSchema("TimeBasedTrkg::TBHits"))
		TimeBasedTrkg__TBHits = new hipo::bank(fFactory->getSchema("TimeBasedTrkg::TBHits"));
	if (fFactory->hasSchema("CTOF::hits"))
		CTOF__hits = new hipo::bank(fFactory->getSchema("CTOF::hits"));
	if (fFactory->hasSchema("RECHB::Particle"))
		RECHB__Particle = new hipo::bank(fFactory->getSchema("RECHB::Particle"));
	if (fFactory->hasSchema("RAW::tdc"))
		RAW__tdc = new hipo::bank(fFactory->getSchema("RAW::tdc"));
	if (fFactory->hasSchema("RICH::photons"))
		RICH__photons = new hipo::bank(fFactory->getSchema("RICH::photons"));
	if (fFactory->hasSchema("HEL::flip"))
		HEL__flip = new hipo::bank(fFactory->getSchema("HEL::flip"));
	if (fFactory->hasSchema("TimeBasedTrkg::Trajectory"))
		TimeBasedTrkg__Trajectory = new hipo::bank(fFactory->getSchema("TimeBasedTrkg::Trajectory"));
	if (fFactory->hasSchema("DETECTOR::ecpb"))
		DETECTOR__ecpb = new hipo::bank(fFactory->getSchema("DETECTOR::ecpb"));
	if (fFactory->hasSchema("CND::hits"))
		CND__hits = new hipo::bank(fFactory->getSchema("CND::hits"));
}

int TIdentificatorCLAS12::FillBanks(){
	if (fFactory->hasSchema("BMTRec::Hits"))
		 fEvent->getStructure(*BMTRec__Hits);
	if (fFactory->hasSchema("RAW::adc"))
		 fEvent->getStructure(*RAW__adc);
	if (fFactory->hasSchema("BAND::adc"))
		 fEvent->getStructure(*BAND__adc);
	if (fFactory->hasSchema("RUN::config"))
		 fEvent->getStructure(*RUN__config);
	if (fFactory->hasSchema("RICH::clusters"))
		 fEvent->getStructure(*RICH__clusters);
	if (fFactory->hasSchema("RECHB::Scintillator"))
		 fEvent->getStructure(*RECHB__Scintillator);
	if (fFactory->hasSchema("REC::RingCher"))
		 fEvent->getStructure(*REC__RingCher);
	if (fFactory->hasSchema("BSTRec::LayerEffs"))
		 fEvent->getStructure(*BSTRec__LayerEffs);
	if (fFactory->hasSchema("RTPC::pos"))
		 fEvent->getStructure(*RTPC__pos);
	if (fFactory->hasSchema("TimeBasedTrkg::TBCrosses"))
		 fEvent->getStructure(*TimeBasedTrkg__TBCrosses);
	if (fFactory->hasSchema("HitBasedTrkg::HBTracks"))
		 fEvent->getStructure(*HitBasedTrkg__HBTracks);
	if (fFactory->hasSchema("CVTRec::Cosmics"))
		 fEvent->getStructure(*CVTRec__Cosmics);
	if (fFactory->hasSchema("RECHB::Cherenkov"))
		 fEvent->getStructure(*RECHB__Cherenkov);
	if (fFactory->hasSchema("BSTRec::Clusters"))
		 fEvent->getStructure(*BSTRec__Clusters);
	if (fFactory->hasSchema("CVTRec::Trajectory"))
		 fEvent->getStructure(*CVTRec__Trajectory);
	if (fFactory->hasSchema("RECHB::Calorimeter"))
		 fEvent->getStructure(*RECHB__Calorimeter);
	if (fFactory->hasSchema("TimeBasedTrkg::TBSegmentTrajectory"))
		 fEvent->getStructure(*TimeBasedTrkg__TBSegmentTrajectory);
	if (fFactory->hasSchema("CTOF::tdc"))
		 fEvent->getStructure(*CTOF__tdc);
	if (fFactory->hasSchema("REC::Cherenkov"))
		 fEvent->getStructure(*REC__Cherenkov);
	if (fFactory->hasSchema("BMTRec::LayerEffs"))
		 fEvent->getStructure(*BMTRec__LayerEffs);
	if (fFactory->hasSchema("FTOF::adc"))
		 fEvent->getStructure(*FTOF__adc);
	if (fFactory->hasSchema("MC::Lund"))
		 fEvent->getStructure(*MC__Lund);
	if (fFactory->hasSchema("DETECTOR::lcpb"))
		 fEvent->getStructure(*DETECTOR__lcpb);
	if (fFactory->hasSchema("MC::Header"))
		 fEvent->getStructure(*MC__Header);
	if (fFactory->hasSchema("CND::clusters"))
		 fEvent->getStructure(*CND__clusters);
	if (fFactory->hasSchema("TimeBasedTrkg::TBCovMat"))
		 fEvent->getStructure(*TimeBasedTrkg__TBCovMat);
	if (fFactory->hasSchema("RICH::hits"))
		 fEvent->getStructure(*RICH__hits);
	if (fFactory->hasSchema("RECHB::Track"))
		 fEvent->getStructure(*RECHB__Track);
	if (fFactory->hasSchema("MC::True"))
		 fEvent->getStructure(*MC__True);
	if (fFactory->hasSchema("BST::adc"))
		 fEvent->getStructure(*BST__adc);
	if (fFactory->hasSchema("MC::Event"))
		 fEvent->getStructure(*MC__Event);
	if (fFactory->hasSchema("HitBasedTrkg::HBCrosses"))
		 fEvent->getStructure(*HitBasedTrkg__HBCrosses);
	if (fFactory->hasSchema("FTOF::clusters"))
		 fEvent->getStructure(*FTOF__clusters);
	if (fFactory->hasSchema("REC::TrackCross"))
		 fEvent->getStructure(*REC__TrackCross);
	if (fFactory->hasSchema("REC::Scintillator"))
		 fEvent->getStructure(*REC__Scintillator);
	if (fFactory->hasSchema("ECAL::peaks"))
		 fEvent->getStructure(*ECAL__peaks);
	if (fFactory->hasSchema("TimeBasedTrkg::TBClusters"))
		 fEvent->getStructure(*TimeBasedTrkg__TBClusters);
	if (fFactory->hasSchema("TAGGER::tgpb"))
		 fEvent->getStructure(*TAGGER__tgpb);
	if (fFactory->hasSchema("LTCC::clusters"))
		 fEvent->getStructure(*LTCC__clusters);
	if (fFactory->hasSchema("RICH::hadrons"))
		 fEvent->getStructure(*RICH__hadrons);
	if (fFactory->hasSchema("DC::tdc"))
		 fEvent->getStructure(*DC__tdc);
	if (fFactory->hasSchema("BSTRec::Crosses"))
		 fEvent->getStructure(*BSTRec__Crosses);
	if (fFactory->hasSchema("HTCC::rec"))
		 fEvent->getStructure(*HTCC__rec);
	if (fFactory->hasSchema("FTHODO::adc"))
		 fEvent->getStructure(*FTHODO__adc);
	if (fFactory->hasSchema("FTOF::matchedclusters"))
		 fEvent->getStructure(*FTOF__matchedclusters);
	if (fFactory->hasSchema("ECAL::adc"))
		 fEvent->getStructure(*ECAL__adc);
	if (fFactory->hasSchema("FMT::adc"))
		 fEvent->getStructure(*FMT__adc);
	if (fFactory->hasSchema("BMT::adc"))
		 fEvent->getStructure(*BMT__adc);
	if (fFactory->hasSchema("ECAL::calib"))
		 fEvent->getStructure(*ECAL__calib);
	if (fFactory->hasSchema("BMTRec::Crosses"))
		 fEvent->getStructure(*BMTRec__Crosses);
	if (fFactory->hasSchema("RUN::rf"))
		 fEvent->getStructure(*RUN__rf);
	if (fFactory->hasSchema("BSTRec::Hits"))
		 fEvent->getStructure(*BSTRec__Hits);
	if (fFactory->hasSchema("HitBasedTrkg::HBSegmentTrajectory"))
		 fEvent->getStructure(*HitBasedTrkg__HBSegmentTrajectory);
	if (fFactory->hasSchema("LTCC::adc"))
		 fEvent->getStructure(*LTCC__adc);
	if (fFactory->hasSchema("CND::tdc"))
		 fEvent->getStructure(*CND__tdc);
	if (fFactory->hasSchema("HTCC::adc"))
		 fEvent->getStructure(*HTCC__adc);
	if (fFactory->hasSchema("RAW::vtp"))
		 fEvent->getStructure(*RAW__vtp);
	if (fFactory->hasSchema("CVTRec::Tracks"))
		 fEvent->getStructure(*CVTRec__Tracks);
	if (fFactory->hasSchema("RECFT::Particle"))
		 fEvent->getStructure(*RECFT__Particle);
	if (fFactory->hasSchema("TimeBasedTrkg::TBTracks"))
		 fEvent->getStructure(*TimeBasedTrkg__TBTracks);
	if (fFactory->hasSchema("DETECTOR::icpb"))
		 fEvent->getStructure(*DETECTOR__icpb);
	if (fFactory->hasSchema("RF::tdc"))
		 fEvent->getStructure(*RF__tdc);
	if (fFactory->hasSchema("HEADER::info"))
		 fEvent->getStructure(*HEADER__info);
	if (fFactory->hasSchema("FTCAL::clusters"))
		 fEvent->getStructure(*FTCAL__clusters);
	if (fFactory->hasSchema("RAW::scaler"))
		 fEvent->getStructure(*RAW__scaler);
	if (fFactory->hasSchema("BMTRec::Clusters"))
		 fEvent->getStructure(*BMTRec__Clusters);
	if (fFactory->hasSchema("TimeBasedTrkg::TBSegments"))
		 fEvent->getStructure(*TimeBasedTrkg__TBSegments);
	if (fFactory->hasSchema("RECHB::ForwardTagger"))
		 fEvent->getStructure(*RECHB__ForwardTagger);
	if (fFactory->hasSchema("REC::Calorimeter"))
		 fEvent->getStructure(*REC__Calorimeter);
	if (fFactory->hasSchema("REC::CovMat"))
		 fEvent->getStructure(*REC__CovMat);
	if (fFactory->hasSchema("RAW::epics"))
		 fEvent->getStructure(*RAW__epics);
	if (fFactory->hasSchema("REC::VertDoca"))
		 fEvent->getStructure(*REC__VertDoca);
	if (fFactory->hasSchema("FTHODO::clusters"))
		 fEvent->getStructure(*FTHODO__clusters);
	if (fFactory->hasSchema("RUN::scaler"))
		 fEvent->getStructure(*RUN__scaler);
	if (fFactory->hasSchema("BAND::hits"))
		 fEvent->getStructure(*BAND__hits);
	if (fFactory->hasSchema("ECAL::hits"))
		 fEvent->getStructure(*ECAL__hits);
	if (fFactory->hasSchema("ECAL::clusters"))
		 fEvent->getStructure(*ECAL__clusters);
	if (fFactory->hasSchema("FT::particles"))
		 fEvent->getStructure(*FT__particles);
	if (fFactory->hasSchema("HitBasedTrkg::HBClusters"))
		 fEvent->getStructure(*HitBasedTrkg__HBClusters);
	if (fFactory->hasSchema("REC::RICH"))
		 fEvent->getStructure(*REC__RICH);
	if (fFactory->hasSchema("REC::Track"))
		 fEvent->getStructure(*REC__Track);
	if (fFactory->hasSchema("RECHB::Event"))
		 fEvent->getStructure(*RECHB__Event);
	if (fFactory->hasSchema("FTOF::hbhits"))
		 fEvent->getStructure(*FTOF__hbhits);
	if (fFactory->hasSchema("RTPC::adc"))
		 fEvent->getStructure(*RTPC__adc);
	if (fFactory->hasSchema("FTCAL::adc"))
		 fEvent->getStructure(*FTCAL__adc);
	if (fFactory->hasSchema("CTOF::rawhits"))
		 fEvent->getStructure(*CTOF__rawhits);
	if (fFactory->hasSchema("CND::adc"))
		 fEvent->getStructure(*CND__adc);
	if (fFactory->hasSchema("FTOF::hits"))
		 fEvent->getStructure(*FTOF__hits);
	if (fFactory->hasSchema("FMTRec::Clusters"))
		 fEvent->getStructure(*FMTRec__Clusters);
	if (fFactory->hasSchema("HEL::adc"))
		 fEvent->getStructure(*HEL__adc);
	if (fFactory->hasSchema("DETECTOR::ccpb"))
		 fEvent->getStructure(*DETECTOR__ccpb);
	if (fFactory->hasSchema("REC::ForwardTagger"))
		 fEvent->getStructure(*REC__ForwardTagger);
	if (fFactory->hasSchema("FTTRK::adc"))
		 fEvent->getStructure(*FTTRK__adc);
	if (fFactory->hasSchema("FTOF::rawhits"))
		 fEvent->getStructure(*FTOF__rawhits);
	if (fFactory->hasSchema("RECHB::TrackCross"))
		 fEvent->getStructure(*RECHB__TrackCross);
	if (fFactory->hasSchema("DETECTOR::scpb"))
		 fEvent->getStructure(*DETECTOR__scpb);
	if (fFactory->hasSchema("FTHODO::hits"))
		 fEvent->getStructure(*FTHODO__hits);
	if (fFactory->hasSchema("RECFT::Event"))
		 fEvent->getStructure(*RECFT__Event);
	if (fFactory->hasSchema("CTOF::adc"))
		 fEvent->getStructure(*CTOF__adc);
	if (fFactory->hasSchema("EVENT::detector"))
		 fEvent->getStructure(*EVENT__detector);
	if (fFactory->hasSchema("DC::doca"))
		 fEvent->getStructure(*DC__doca);
	if (fFactory->hasSchema("HTCC::tdc"))
		 fEvent->getStructure(*HTCC__tdc);
	if (fFactory->hasSchema("HEL::online"))
		 fEvent->getStructure(*HEL__online);
	if (fFactory->hasSchema("REC::Event"))
		 fEvent->getStructure(*REC__Event);
	if (fFactory->hasSchema("REC::Particle"))
		 fEvent->getStructure(*REC__Particle);
	if (fFactory->hasSchema("REC::Traj"))
		 fEvent->getStructure(*REC__Traj);
	if (fFactory->hasSchema("FTOF::tdc"))
		 fEvent->getStructure(*FTOF__tdc);
	if (fFactory->hasSchema("HitBasedTrkg::HBSegments"))
		 fEvent->getStructure(*HitBasedTrkg__HBSegments);
	if (fFactory->hasSchema("FTCAL::hits"))
		 fEvent->getStructure(*FTCAL__hits);
	if (fFactory->hasSchema("FMTRec::Hits"))
		 fEvent->getStructure(*FMTRec__Hits);
	if (fFactory->hasSchema("BAND::tdc"))
		 fEvent->getStructure(*BAND__tdc);
	if (fFactory->hasSchema("RF::adc"))
		 fEvent->getStructure(*RF__adc);
	if (fFactory->hasSchema("FMTRec::Crosses"))
		 fEvent->getStructure(*FMTRec__Crosses);
	if (fFactory->hasSchema("MC::Particle"))
		 fEvent->getStructure(*MC__Particle);
	if (fFactory->hasSchema("ECAL::tdc"))
		 fEvent->getStructure(*ECAL__tdc);
	if (fFactory->hasSchema("RICH::tdc"))
		 fEvent->getStructure(*RICH__tdc);
	if (fFactory->hasSchema("LTCC::tdc"))
		 fEvent->getStructure(*LTCC__tdc);
	if (fFactory->hasSchema("ECAL::moments"))
		 fEvent->getStructure(*ECAL__moments);
	if (fFactory->hasSchema("EVENT::particle"))
		 fEvent->getStructure(*EVENT__particle);
	if (fFactory->hasSchema("RUN::trigger"))
		 fEvent->getStructure(*RUN__trigger);
	if (fFactory->hasSchema("HitBasedTrkg::HBHits"))
		 fEvent->getStructure(*HitBasedTrkg__HBHits);
	if (fFactory->hasSchema("TimeBasedTrkg::TBHits"))
		 fEvent->getStructure(*TimeBasedTrkg__TBHits);
	if (fFactory->hasSchema("CTOF::hits"))
		 fEvent->getStructure(*CTOF__hits);
	if (fFactory->hasSchema("RECHB::Particle"))
		 fEvent->getStructure(*RECHB__Particle);
	if (fFactory->hasSchema("RAW::tdc"))
		 fEvent->getStructure(*RAW__tdc);
	if (fFactory->hasSchema("RICH::photons"))
		 fEvent->getStructure(*RICH__photons);
	if (fFactory->hasSchema("HEL::flip"))
		 fEvent->getStructure(*HEL__flip);
	if (fFactory->hasSchema("TimeBasedTrkg::Trajectory"))
		 fEvent->getStructure(*TimeBasedTrkg__Trajectory);
	if (fFactory->hasSchema("DETECTOR::ecpb"))
		 fEvent->getStructure(*DETECTOR__ecpb);
	if (fFactory->hasSchema("CND::hits"))
		 fEvent->getStructure(*CND__hits);
}

