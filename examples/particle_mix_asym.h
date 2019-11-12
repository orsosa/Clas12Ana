#include <getopt.h>
#include "TMath.h"
extern bool GSIM;
extern int data_type;
extern long Ne;
extern Float_t HELIC;
extern Int_t MULT;

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

/// #########################


inline void printhelp()
{
  std::cout<<"####### Help #########\n"
    "\t[-n | --max-events] <N>          : max number of events to be considered. Default: all\n"
    "\t[-l | --helicity] <l>            : fix helicity value. Default: read from file\n"
    "\t[-m | --multiplicity] <m>        : Multiply statistic\n" 
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
    {"help",     no_argument,       0, 'h'},
    {"helicity", required_argument,0,'l'},
    {"multiplicity", required_argument,0,'m'},
    {0, 0, 0, 0}
  };

  //  if(argc==1)
  //  printhelp();
  while ( (c = getopt_long(argc, argv, "hn:l:m:", long_options, &option_index))  != -1)
    switch (c)
      {
      case 'n':
        Ne = atoi(optarg);
        break;
      case 'm':
        MULT = atoi(optarg);
        break;

      case 'l':
	HELIC=atof(optarg);
	std::cout<<"Setting helicity to: " <<HELIC<<std::endl;
	break;
      case 'h':
	printhelp();
	break;
      case '?':
        if (optopt == 't')
          fprintf (stderr, "Option -%c data type [data | simrec | gsim ]. Default data.\n", optopt);
        else if (isprint (optopt))
          fprintf (stderr, "Unknown option `-%c'.\n", optopt);
        else
          fprintf (stderr,
                   "Unknown option character `\\x%x'.\n",
                   optopt);
        return 1;
      default:
        abort ();
      }
}
