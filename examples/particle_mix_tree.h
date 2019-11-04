#include <getopt.h>
extern bool GSIM, EFLAG;
extern int data_type;
extern long Ne;
extern Float_t HELIC;
extern TString INDIR,INFILE,REACTION,OFILE;

inline void printhelp()
{
  std::cout<<"####### Help #########\n"
    "\t[-t | --data-type] <target_type> : Data type [data | simrec | mc ]. Default: data\n"
    "\t[-n | --max-events] <N>          : Max number of events to be considered. Default: all\n"
    "\t[-l | --helicity] <l>            : Fix helicity value. Default: read from file\n"
    "\t[-d | --in-dir] <d>              : Path containing the root input files. All files will be read <d>/*.root\n"
    "\t                                   It has priority over -f option\n"
    "\t[-f | --in-file] <f>             : Input file. Default value data.root\n"
    "\t[-r | --reaction] <r>            : Select the mixing particle reaction.\n"
    "\t                                   e.g.: rho:pi+_pi-, :pi+_pi-. Default :pi+_pi-\n"
    "\t-e                               : exact match. Default false.\n"
    "\t-o                               : output file name. Default outfiles/pippim_all.root.\n"
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
    {"data-type",  required_argument, 0, 't'},
    {"helicity", required_argument,0,'l'},
    {"in-dir", required_argument,0,'d'},
    {"in-file", required_argument,0,'f'},
    {"reaction", required_argument,0,'r'},
    {"output", required_argument,0,'o'},
    {"exact-match", no_argument,0,'e'},
    {0, 0, 0, 0}
  };

  if(argc==1)
    printhelp();
  while ( (c = getopt_long(argc, argv, "het:n:l:d:f:r:o:", long_options, &option_index))  != -1)
    switch (c)
      {
      case 'e':
        EFLAG = true;
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
      case 'n':
        Ne = atol(optarg);
        break;
      case 't':
	if(!strcmp(optarg,"data"))
	  data_type = 0;
	else if (!strcmp(optarg,"simrec")) 
	  data_type = 1;
	else if (!strcmp(optarg,"mc"))
	  data_type = 2; 
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
          fprintf (stderr, "Option -%c data type [data | simrec | mc ]. Default data.\n", optopt);
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
