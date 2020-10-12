#include <getopt.h>
extern bool GSIM, EFLAG, QUIET;
extern long Ne, START;
extern Float_t HELIC;
extern TString INDIR, INFILE, REACTION, OFILE;
extern TRandom3 *rndm;

inline void printhelp()
{
  std::cout<<"####### Help #########\n"
    "\t[-n | --max-events] <N>          : Max number of events to be considered. Default: all\n"
    "\t[-s | --start-event] <N>         : Start event number. Default: 0\n"
    "\t[-d | --in-dir] <d>              : Path containing the root input files. All files will be read <d>/*.root\n"
    "\t                                   It has priority over -f option\n"
    "\t[-f | --in-file] <f>             : Input file. Default value data.root\n"
    "\t[-r | --reaction] <r>            : Select the mixing particle reaction.\n"
    "\t                                   e.g.: rho:pi+_pi-, :pi+_pi-. Default :pi+_pi-\n"
    "\t-e                               : exact match. Default false.\n"
    "\t-l                               : set helicity internally. Default read from data.\n"
    "\t-o                               : output file name. Default outfiles/pippim_all.root.\n"
    "\t-1                               : quiet mode, minimum printing for farm execution. default false.\n"
    "\t-h                               : Print help.\n"
    "#########################"	   <<std::endl;
  exit(0);
}

inline int parseopt(int argc, char* argv[])
{
  int c;
  int option_index = 0;
  static struct option long_options[] =
  {
    {"max-events",     required_argument,       0, 'n'},
    {"start-event",     required_argument,       0, 's'},
    {"help",     no_argument,       0, 'h'},
    {"in-dir", required_argument,0,'d'},
    {"in-file", required_argument,0,'f'},
    {"reaction", required_argument,0,'r'},
    {"output", required_argument,0,'o'},
    {"exact-match", no_argument,0,'e'},
    {"manual-helicity", no_argument,0,'l'},
    {"quiet", no_argument,0,'q'},
    {0, 0, 0, 0}
  };

  if(argc==1)
    printhelp();
  while ( (c = getopt_long(argc, argv, "hlen:d:f:r:o:s:q", long_options, &option_index))  != -1)
    switch (c)
      {
      case 'e':
        EFLAG = true;
        break;
      case 'l':
        HELIC = 111;
        break;
      case 'o':
        OFILE = optarg;
        break;
      case 'r':
        REACTION = optarg;
        break;
      case 'f':
        INFILE = optarg;
        break;
      case 'd':
        INDIR = optarg;
        break;
      case 's':
        START = atol(optarg);
        break;
      case 'n':
        Ne = atol(optarg);
        break;
      case 'q':
        QUIET = true;
        break;
      case 'h':
	printhelp();
	break;
      case '?':
        if (isprint (optopt))
          fprintf (stderr, "Unknown option `-%c'.\n", optopt);
        else
          fprintf (stderr,
                   "Unknown option character `\\x%x'.\n",
                   optopt);
        return 1;
      default:
        abort ();
      }
  return 0;
} 

inline Bool_t check_reaction()
{
  TDatabasePDG db;
  Ssiz_t i_c=0,i_p=0,sl=0;
  i_c = REACTION.Index(":");
  TString primary="", secondary="";
  primary = REACTION(i_p,i_c);
  if (primary != "" && !db.GetParticle(primary)){
    std::cout<<primary<<" --- primary particle doesn't exist on TDatabasePDG"<<std::endl;
    return false;
  }
  i_p=i_c+1;
  while ((i_c=REACTION.Index("_",i_p)) != -1){
    sl=i_c-i_p;
    secondary=REACTION(i_p,sl);
    if (secondary != "" && !db.GetParticle(secondary)){
      std::cout<<secondary<<" --- secondary particle doesn't exist on TDatabasePDG"<<std::endl;
      return false;
    }
    i_p=i_c+1;
  }
  secondary=REACTION(i_p,REACTION.Length());
  if (secondary != "" && !db.GetParticle(secondary)){
    std::cout<<secondary<<" --- secondary particle doesn't exist on TDatabasePDG"<<std::endl;
    return false;
  }

  return true;
    
}

inline TString get_primary(){
  Ssiz_t i_c=0,i_p=0;
  i_c = REACTION.Index(":");
  TString primary="";
  primary = REACTION(i_p,i_c);
  return primary;
}

inline TString pop_secondary(Ssiz_t &start){
  TString secondary = "";
  REACTION.Tokenize(secondary,start,"_");
  return secondary;
}


//############### Asymmetry estimation ###########
Float_t Ax0 = 0.0606961;
Float_t mx = 0.139071;

//Float_t Az0 = 0.0128893;
//Float_t z0 = 0.63336;
//Float_t cz = 1.47663;
Float_t Az0 = 0.0128893;
Float_t z0 = 0.63336;
Float_t cz = 0.2;

Float_t Am0 = 0.0420307;
Float_t m0 = 0.811456;
Float_t cm = 0.300714;

Float_t m_Az = 0.0264659;
Float_t m_Am = 0.0271934;
Float_t m_Ax = 0.0251313;
//Float_t max_Az = 0.3;
Float_t max_Az = 0.051;
Float_t pi = TMath::Pi();

Float_t ALU(Float_t x,Float_t z,Float_t m)
{
  Float_t Ax = Ax0 - mx*x;
  Float_t Az = cz*(z-z0)*(z-z0) + Az0;
  Float_t Am = -cm*(m-m0)*(m-m0) + Am0;
  return Ax*Az*Am/(m_Az*m_Am);
}

Float_t Ax(Float_t x)
{
  return Ax0 - mx*x;
}

Float_t  Az(Float_t z)
{
  return cz*(z-z0)*(z-z0) + Az0;
}

Float_t Am(Float_t m)
{
  return -cm*(m-m0)*(m-m0) + Am0;
}
///////// asymetry function /////////
Float_t asym_n(Float_t z,Float_t phiR,Float_t h)
{
  h = h==0?-1:h;
  return (1. + h*Az(z)*sin(phiR*pi/180.) )/(1+max_Az);
}

/*
Float_t get_helicity(Float_t phi,Float_t A = 0.1)
{
  Float_t w = rndm->Uniform()>=0.5?1:-1; 
  Float_t uval = rndm->Uniform(); 

  Float_t aval =  (1 + w*A*sin(phi*TMath::DegToRad()))/2.; // max efficiency.
  
  if (w>0)
    return (uval<=aval?1:-1);
  else 
    return (uval<=aval?-1:1);
}
*/
Float_t get_helicity(Float_t phi, Float_t A = 0.1, Float_t sinth = 1.0)
{
  Float_t w = rndm->Uniform()>=0.5?1:-1; 
  Float_t uval = rndm->Uniform(); 

  Float_t aval =  ( 1 + w*A*sin(phi*TMath::DegToRad())*sinth ) / 2.; // max efficiency.
  if (w>0)
    return (uval<=aval?1:-1);
  else 
    return (uval<=aval?-1:1);
}


Float_t get_UxS(Float_t phi, Float_t B = 0.2, Float_t C = 0.15){

  Float_t w = rndm->Uniform()>=0.5?1:-1; 
  Float_t uval = rndm->Uniform(); 

  Float_t aval =  ( 1 + w*(B*cos(phi*TMath::DegToRad()) + C*cos(2*phi*TMath::DegToRad())) ) / 2.; // max efficiency.
  if (w>0)
    return (uval<=aval?1:-1);
  else 
    return (uval<=aval?-1:1);
}

/// #########################
