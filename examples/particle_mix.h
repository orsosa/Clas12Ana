#include <getopt.h>
extern bool GSIM;
extern int data_type;
extern long Ne;
extern Float_t HELIC;
extern TString INDIR,INFILE;

inline void printhelp()
{
  std::cout<<"####### Help #########\n"
    "\t[-t | --data-type] <target_type> : Data type [data | simrec | mc ]. Default: data\n"
    "\t[-n | --max-events] <N>          : Max number of events to be considered. Default: all\n"
    "\t[-l | --helicity] <l>            : Fix helicity value. Default: read from file\n"
    "\t[-d | --in-dir] <d>              : Path containing the root input files. All files will be read <d>/*.root\n"
    "\t                                   It has priority over -f option\n"
    "\t[-f | --in-file] <f>             : Input file. Default value data.root\n" 
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
    {0, 0, 0, 0}
  };

  if(argc==1)
    printhelp();
  while ( (c = getopt_long(argc, argv, "ht:n:l:d:f:", long_options, &option_index))  != -1)
    switch (c)
      {
      case 'f':
        INFILE = optarg;
        break;
      case 'd':
        INDIR = optarg;
        break;
      case 'n':
        Ne = atoi(optarg);
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
}
