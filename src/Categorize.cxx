#include "TIdentificatorCLAS12.h"

TString TIdentificatorCLAS12::GetCategorization(Int_t k)
{
  TString partId;
  partId = "not recognized";
  if (Pid(k) == 11) // fiducial and status cuts to be added
    partId="electron";

  if (Pid(k) == 22
      // &&0.95<Beta(k)&&Beta(k)<1.05

      )// fiducial and status cuts to be added
    partId="gamma";
  
  if (Pid(k) == 211 
      )// fiducial and status cuts to be added
    partId="pi+";

  if (Pid(k) == -211
      )// fiducial and status cuts to be added
    partId="pi-";

  if (Pid(k) == 2212
      )// fiducial and status cuts to be added
    partId="proton";

  if (Pid(k) == 321
      )// fiducial and status cuts to be added
    partId="K+";

  if (Pid(k) == -321
      )// fiducial and status cuts to be added
    partId="K-";

  if (Pid(k) == 2112
      )// fiducial and status cuts to be added
    partId="neutron";

  if (Pid(k) == 45
      )// fiducial and status cuts to be added
    partId="deuteron";
  
  return partId;
}
/*
TString TIdentificatorCLAS12::GetCategorization(Int_t k, const char* tt)
{
  TString partId;
  partId = "not recognized";

  if (number_dc != 0) {
    if (k == 0 &&
	Status(0) > 0 && Status(0) < 100 &&
	Charge(0) == -1 &&
	number_cc != 0 && number_ec != 0 && number_sc != 0 &&
	StatCC(0) > 0 && StatSC(0) > 0 &&
	StatDC(0) > 0 && StatEC(0) > 0 &&
	DCStatus(0) > 0 && SCStatus(0) == 33 &&
	( !mflag || Nphe(0) > (Sector(0)==0||Sector(0)==1)*25 //Added sector dependent cut. osoto_mod.
	  +(Sector(0)==2)*26 
	  +(Sector(0)==3)*21
	  +(Sector(0)==4 || Sector(0)==5 )*28)&&
	(mflag || Nphe(0)>25)&& //default hayk cc nphe.
	(!mflag || Momentum(0)>0.75)&& // Momentum triger added. osoto_mod.
	(!mflag || Ein(0)>0.06)&& // Inner stack energy cut. osoto_mod.
	(!mflag || (TimeEC(0)- TimeSC(0) - (PathEC(0)-PathSC(0) )/30) < 5*0.35 )&&// elapsed time between sc and ec (??). osoto_mod.
	(!mflag || SampFracCheck(tt))&& //osoto mod
	Etot(0) / 0.27 / 1.15 + 0.4 > Momentum(0) &&
	Etot(0) / 0.27 / 1.15 - 0.2 < Momentum(0) &&
	Ein(0) + Eout(0) > 0.8 * 0.27 * Momentum(0) &&
	Ein(0) + Eout(0) < 1.2 * 0.27 * Momentum(0) &&
	Eout(0) != 0 && FidCheckCut() == 1)
      partId = "electron";
    
    if (k > 0)
    {
      if (Charge(k) == 0 )
      { 
	TVector3 *ECxyz = new TVector3(XEC(k),YEC(k),ZEC(k));
	TVector3 *ECuvw = XYZToUVW(ECxyz);
	if(ECuvw->X()>40 && ECuvw->X()<410 && // u coordinate in ]40,410[
	   ECuvw->Y()<370 && // v coordinate in [0,370[
	   ECuvw->Z()<410 && // w coordinate in [0,410[
	   PathEC(k)/(Betta(k)*30) - PathEC(k)/30>-2.2&&
	   PathEC(k)/(Betta(k)*30) - PathEC(k)/30<1.3&&
	   TMath::Max(Etot(k), Ein(k)+Eout(k))/0.273>0.1)
         {
	   partId = "gamma";
	 }
              
      }
      //negative particles osoto_mod
      if (Charge(k) == -1 &&
	  Status(k) > 0 && Status(k) < 100 &&
	  StatDC(k) > 0 && DCStatus(k) > 0)
      {
	Float_t P=Momentum(k);
	Float_t T4 =TimeCorr4(0.1396,k);
	if (Etot(k)<0.15 && Ein(k)<0.085-0.5*Eout(k) &&
	    ( ( (!(StatCC(k)>0 && Nphe(k)>25)) &&
		( 
		 (0<P && P<=0.5 && T4>-0.87 && T4<0.63)
		 ||(0.5<P && P<=1.0 && T4>-0.55 && T4<0.37)
		 ||(1.0<P && P<=1.5 && T4>-0.55 && T4<0.38)
		 ||(1.5<P && P<=2.0 && T4>-0.60 && T4<0.44)
		 ||(2.0<P && P<=2.5 && T4>-1.00 && T4<0.45)
		  )
		)
	      ||(2.5<P && P<=3.0 && T4>-1.00 && T4<0.40)
	      || (3.0<P && T4>-2.00 && T4<0.45)
	      )
	    )
	  partId = "pi-";
      }
      if (    Status(k) > 0 && Status(k) < 100 &&
	      Charge(k) == -1 &&
	      number_cc != 0 && number_ec != 0 && number_sc != 0 &&
	      StatCC(k) > 0 && StatSC(k) > 0 &&
	      StatDC(k) > 0 && StatEC(k) > 0 &&
	      DCStatus(k) > 0 && SCStatus(k) == 33 &&
	      Nphe(0) > (Sector(k)==0||Sector(k)==1)*25 //Added sector dependent cut. osoto_mod.
	      +(Sector(k)==2)*26 
	      +(Sector(k)==3)*21
	      +(Sector(k)==4 || Sector(k)==5 )*28 &&
	      Momentum(k)>0.75&& // Momentum triger added. osoto_mod.
	      Ein(k)>0.06&& // Inner stack energy cut. osoto_mod.
	      Etot(k) / 0.27 / 1.15 + 0.4 > Momentum(k) &&
	      Etot(k) / 0.27 / 1.15 - 0.2 < Momentum(k) &&
	      Ein(k) + Eout(k) > 0.8 * 0.27 * Momentum(k) &&
	      Ein(k) + Eout(k) < 1.2 * 0.27 * Momentum(k) &&
	      Eout(k) != 0)
	partId = "s_electron";

      //positive particles
      if (Charge(k) == 1 &&
	  Status(k) > 0 && Status(k) < 100 &&
	  StatDC(k) > 0 && DCStatus(k) > 0)
      {
	if (Momentum(k)>=2.7)
	{
	  if (number_cc != 0 && StatCC(k) > 0 &&
	      Nphe(k) > 25 && Chi2CC(k) < 5 / 57.3)
	    partId = "high energy pion +";
	}
	if (Momentum(k) < 2.7)
	{
	  if (number_sc != 0 && StatSC(k) > 0 &&
	      ((Momentum(k) > 0 && Momentum(k) <= 0.25 && TimeCorr4(0.139,k) >= -1.45 && TimeCorr4(0.139,k) <= 1.05) ||
	       (Momentum(k) > 0.25 && Momentum(k) <= 0.5 && TimeCorr4(0.139,k) >= -1.44 && TimeCorr4(0.139,k) <= 1.05) || 
	       (Momentum(k) > 0.5 && Momentum(k) <= 0.75 && TimeCorr4(0.139,k) >= -1.44 && TimeCorr4(0.139,k) <= 1.05) ||
	       (Momentum(k) > 0.75 && Momentum(k) <= 1 && TimeCorr4(0.139,k) >= -1.4 && TimeCorr4(0.139,k) <= 1.05) ||
	       (Momentum(k) > 1 && Momentum(k) <= 1.25 && TimeCorr4(0.139,k) >= -1.35 && TimeCorr4(0.139,k) <= 1.03) ||
	       (Momentum(k) > 1.25 && Momentum(k) <= 1.5 && TimeCorr4(0.139,k) >= -1.35 && TimeCorr4(0.139,k) <= 0.95) ||
	       (Momentum(k) > 1.5 && Momentum(k) <= 1.75 && TimeCorr4(0.139,k) >= -1.35 && TimeCorr4(0.139,k) <= 0.87) ||
	       (Momentum(k) > 1.75 && Momentum(k) <= 2 && TimeCorr4(0.139,k) >= -1.25 && TimeCorr4(0.139,k) <= 0.68) ||
	       (Momentum(k) > 2 && Momentum(k) <= 2.25 && TimeCorr4(0.139,k) >= -0.95 && TimeCorr4(0.139,k) <= 0.65) ||
	       (Momentum(k) > 2.25 && Momentum(k) <= 2.5 && TimeCorr4(0.139,k) >= -1.05 && TimeCorr4(0.139,k) <= 0.61
		&& Mass2(k) < 0.5) ||
	       (Momentum(k) > 2.5 && Momentum(k) < 2.7 && TimeCorr4(0.139,k) >= -1.05 && TimeCorr4(0.139,k) <= 0.61
		&& Mass2(k) < 0.4)))
	    partId = "low energy pion +";
	}

	if (Momentum(k) < 2.)
	{
	  if (number_sc != 0 && StatSC(k) > 0 &&
	      ((Momentum(k) >= 1. &&
		TimeCorr4(0.938,k) >= -0.69 &&
		TimeCorr4(0.938,k) <= 1.38) ||
	       (Momentum(k) < 1. &&
		TimeCorr4(0.938,k) >= -3.78 &&
		TimeCorr4(0.938,k) <= 6.75)))
	    partId = "low energy proton";
	}
	
	if (Charge(k) == 1 && number_cc != 0 &&
	    number_ec != 0 && number_sc != 0 &&
	    StatCC(k) > 0 && StatSC(k) > 0 &&
	    StatDC(k) > 0 && StatEC(k) > 0 &&
	    DCStatus(k) > 0 && Nphe(k) > 25 &&
	    Etot(k) / 0.27 / 1.15 + 0.4 > Momentum(k) &&
	    Etot(k) / 0.27 / 1.15 - 0.2 < Momentum(k) &&
	    Ein(k) + Eout(k) > 0.8 * 0.27 * Momentum(k) &&
	    Ein(k) + Eout(k) < 1.2 * 0.27 * Momentum(k))
	  partId = "positron";
      }

    }
  }

  return partId;
}

TString TIdentificatorCLAS12::GetCategorizationOld(Int_t k){
    Int_t number_dc = fCT->GetNRows("DCPB");
    Int_t number_cc = fCT->GetNRows("CCPB");
    Int_t number_sc = fCT->GetNRows("SCPB");
    Int_t number_ec = fCT->GetNRows("ECPB");

    TString partId;

    partId = "not recognized";

    if (number_dc != 0) {
        if (k == 0 &&
                    Status(0) > 0 && Status(0) < 100 &&
                    Charge(0) == -1 &&
                    number_cc != 0 && number_ec != 0 && number_sc != 0 &&
                    StatCC(0) > 0 && StatSC(0) > 0 &&
                    StatDC(0) > 0 && StatEC(0) > 0 &&
                    DCStatus(0) > 0 && SCStatus(0) == 33 &&
                    Nphe(0) > 25 &&
                    Etot(0) / 0.27 + 0.4 > Momentum(0) &&
                    Etot(0) / 0.27 - 0.4 < Momentum(0) &&
                    Ein(0) + Eout(0) > 0.8 * 0.27 * Momentum(0) &&
                    Ein(0) + Eout(0) < 1.2 * 0.27 * Momentum(0) &&
                    FidCheckCut() == 1)
            partId = "electron";


        if (k > 0) {
            if (Charge(k)==0 && Betta(k)>0.95 && ECStatus(k)>0)
                partId = "photon";

            //positive particles
            if (Charge(k) == 1 &&
                        Status(k) > 0 && Status(k) < 100 &&
                        StatDC(k) > 0 && DCStatus(k) > 0) {
                if (Momentum(k)>=2.7) {
                    if (number_cc != 0 && StatCC(k) > 0 &&
                                Nphe(k) > 25 && Chi2CC(k) < 5 / 57.3)
                        partId = "high energy pion +";
                }

                if (Momentum(k) < 2.7) {
                    if (number_sc != 0 && StatSC(k) > 0 &&
                                ((Momentum(k) < 1 &&
                                    TimeCorr4(0.139,k) >= -1.46 &&
                                    TimeCorr4(0.139,k) <= 0.15) ||
                                (Momentum(k) >=1 &&
                                    TimeCorr4(0.139,k) >= -1.38 &&
                                    TimeCorr4(0.139,k) <= 0.53)))
                        partId = "low energy pion +";
                }

                if (Momentum(k) < 2.) {
                    if (number_sc != 0 && StatSC(k) > 0 &&
                                ((Momentum(k) >= 1. &&
                                    TimeCorr4(0.938,k) >= -0.69 &&
                                    TimeCorr4(0.938,k) <= 1.38) ||
                                (Momentum(k) < 1. &&
                                    TimeCorr4(0.938,k) >= -3.78 &&
                                    TimeCorr4(0.938,k) <= 6.75)))
                        partId = "low energy proton";
                }

                if (Charge(k) == 1 && number_cc != 0 &&
                            number_ec != 0 && number_sc != 0 &&
                            StatCC(k) > 0 && StatSC(k) > 0 &&
                            StatDC(k) > 0 && StatEC(k) > 0 &&
                            DCStatus(k) > 0 && Nphe(k) > 25 &&
                            Etot(k) / 0.27 + 0.4 > Momentum(k) &&
                            Etot(k) / 0.27 - 0.4 < Momentum(k) &&
                            Ein(k) + Eout(k) > 0.8 * 0.27 * Momentum(k) &&
                            Ein(k) + Eout(k) < 1.2 * 0.27 * Momentum(k))
                    partId = "positron";
            }
        }
    }
    return partId;
}



TString* TIdentificatorCLAS12::GetCategorization()
{
    Int_t number = fCT->GetNRows("EVNT");

    if (fPartIds != 0) delete [] fPartIds;
    fPartIds = new TString[number];

    if (number != 0) {
        for (Int_t i = 0; i < number; i++)
            fPartIds[i] = GetCategorization(i);
    }

    return fPartIds;
}



void TIdentificatorCLAS12::PrintCategorization()
{
    Int_t number = fCT->GetNRows("EVNT");
    fPartIds = GetCategorization();

    if (fPartIds->CompareTo("electron") == 0) {
        for (Int_t i = 0; i < number; i++)
            cout << *(fPartIds+i) << endl;
        cout << endl;
    }
}



void TIdentificatorCLAS12::PrintCategorization(TString* partIds)
{
    Int_t number = fCT->GetNRows("EVNT");

    if (partIds->CompareTo("electron") == 0) {
        for (Int_t i = 0; i < number; i++)
            cout << *(partIds+i) << endl;
        cout << endl;
    }
}

TString TIdentificatorCLAS12::GetCategorizationMin(Int_t k)
{
    Int_t number_dc = fCT->GetNRows("DCPB");
    Int_t number_cc = fCT->GetNRows("CCPB");
    Int_t number_sc = fCT->GetNRows("SCPB");
    Int_t number_ec = fCT->GetNRows("ECPB");

    TString partId;

    partId = "not recognized";

    if (number_dc != 0) 
    {
      if (k == 0 &&
          Status(0) > 0 &&
          Charge(0) == -1 &&
	  //number_cc != 0 &&
	  number_ec != 0 && number_sc != 0 &&
          StatSC(0) > 0 &&
          StatDC(0) > 0 && 
	  StatEC(0) > 0 &&
          DCStatus(0) > 0 &&
	  ECStatus(0) > 0 && SCStatus(0)>0
          )
          {
            partId = "electron";
          }

        //positive particles
      if (k > 0) 
      {
        if (Charge(k) == 1 && Status(k) > 0 && Status(k) < 100 &&
                        StatDC(k) > 0 && DCStatus(k) > 0 && FidCheckCutPiPlus(k) == 1) 
        {
          if (Momentum(k)>=2.7) 
          {
            if (number_cc != 0 && StatCC(k) > 0 && Nphe(k) > 25 && Chi2CC(k) < 5 / 57.3)
                        partId = "high energy pion +";
          }

          if (Momentum(k) < 2.7) 
          {
            if (number_sc != 0 && StatSC(k) > 0 && TimeCorr4(0.139,k) <= 0.55)
              partId = "low energy pion +";
          }

        }
        if (Status(k) > 0 && Status(k) < 100&&
            Charge(k) == 1 && number_cc != 0 && number_ec != 0 && number_sc != 0 &&
            StatCC(k) > 0 && StatSC(k) > 0 &&
            StatDC(k) > 0 && StatEC(k) > 0 &&
            DCStatus(k) > 0)
            partId = "positron";
        // Gamma ID
        if (Charge(k) == 0 )
          { 
            TVector3 *ECxyz = new TVector3(XEC(k),YEC(k),ZEC(k));
            TVector3 *ECuvw = XYZToUVW(ECxyz);
            if(ECuvw->X()>40 && ECuvw->X()<410 && // u coordinate in ]40,410[
              ECuvw->Y()<370 && // v coordinate in [0,370[
              ECuvw->Z()<410 && // w coordinate in [0,410[
              PathEC(k)/(Betta(k)*30) - PathEC(k)/30>-2.2&&
              PathEC(k)/(Betta(k)*30) - PathEC(k)/30<1.3&&
              TMath::Max(Etot(k), Ein(k)+Eout(k))/0.273>0.1
              )
              {
                partId = "gamma";
              }
          
        }
      }
    }

    return partId;
}
*/
