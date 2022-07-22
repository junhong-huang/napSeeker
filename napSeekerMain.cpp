/***********************************************************************
* napSeeker: a software for discovering novel napRNAs from NAP-seq data
*$2021/9/09/$ @Jian-Hua Yang yangjh7@mail.sysu.edu.cn
************************************************************************/
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<ctype.h>
#include <getopt.h>
#include "BamReader.h"
#include "BamAux.h"
using namespace BamTools;
#include <map>
#include <algorithm>
#include <ios>
#include <iostream>
#include <string>
#include <ostream>
#include <fstream>
#include <iomanip>
#include <locale>
#include <sstream>
#include <vector>
using namespace std;

#include "bioUtils.h"
#include "faiFile.h"
#include "bedFile.h"
#include "samFile.h"
#include "homer_statistics.h"
#include "statistic.h"

#include "napSeeker.h"

char version[] = "napSeeker version 0.1";
void usage(void);

int main(int argc, char *argv[])
{
  char *outFile   = NULL;
  char *inputFile = NULL;
  char *faFile    = NULL;
  char *faiFile   = NULL;
  FILE *genomefp  = NULL;
  FILE *faifp     = NULL;
  FILE *outfp     = NULL;
  int showVersion = 0;
  int showHelp    = 0;
  int c           = 0;
  int i           = 0;
  struct parameterInfo paraInfo;
  /* parse commmand line parameters */

  if (argc == 1)
  {
    usage();
  }

  const char *shortOptions = "vhVkPfsSnNo:p:l:L:t:r:a:A:i:m:M:c:F:e:";

  const struct option longOptions[] =
  {
    { "verbose" , no_argument , NULL, 'v' },
    { "help" , no_argument , NULL, 'h' },
    { "version" , no_argument , NULL, 'V' },
    { "keep-dup" , no_argument , NULL, 'k' },
    { "pair" , no_argument , NULL, 'P' },
    { "fold" , no_argument , NULL, 'f' },
    { "strand" , no_argument , NULL, 's' },
    { "split" , no_argument , NULL, 'S' },
    { "norm" , no_argument , NULL, 'n' },
    { "nano" , no_argument , NULL, 'N' },
    { "outdir" , required_argument , NULL, 'o' },
    { "pvalue" , required_argument, NULL, 'p' },
    { "min-len" , required_argument, NULL, 'l' },
    { "max-len" , required_argument, NULL, 'L' },
    { "min" , required_argument, NULL, 'm' },
    { "max" , required_argument, NULL, 'M' },
    { "cov" , required_argument, NULL, 'c' },
    { "mfold" , required_argument, NULL, 'F' },
    { "sefold" , required_argument, NULL, 'e' },
    { "min-tag" , required_argument, NULL, 't' },
    { "rpm" , required_argument , NULL, 'r' },
    { "fa" , required_argument , NULL, 'a' },
    { "fai" , required_argument , NULL, 'A' },
    { "input" , required_argument , NULL, 'i' },
    {NULL, 0, NULL, 0} ,  /* Required at end of array. */
  };

  paraInfo.verbose    = 0;
  paraInfo.keepDup    = 0;
  paraInfo.genomeSize = 0;
  paraInfo.pvalue     = 1e-5;
  paraInfo.minTag      = 5;
  paraInfo.minLen      = 15;
  paraInfo.maxLen      = 5000;
  paraInfo.strand      = 0;
  paraInfo.rpm         = 0.05;
  paraInfo.norm        = 0;
  paraInfo.insertSize  = 1000000;
  paraInfo.minClusterLen = 20;
  paraInfo.maxClusterLen = 10000;
  paraInfo.fold          = 1;
  paraInfo.minCov        = 2.0;
  paraInfo.mfold         = 2.0;
  paraInfo.sefold        = 5.0;
  paraInfo.split         = 0;
  paraInfo.pairEnd       = 0;
  paraInfo.nanopore      = 0;

  while ((c = getopt_long(argc, argv, shortOptions, longOptions, NULL)) >= 0)
  {
    switch (c)
    {
    case 'v':
      paraInfo.verbose = 1;
      break;
    case 'h':
      showHelp = 1;
      break;
    case 'V':
      showVersion = 1;
      break;
    case 'k':
      paraInfo.keepDup = 1;
      break;
    case 'P':
      paraInfo.pairEnd = 1;
      break;
    case 's':
      paraInfo.strand = 1;
      break;
    case 'S':
      paraInfo.split = 1;
      break;
    case 'n':
      paraInfo.norm = 1;
      break;
    case 'N':
      paraInfo.nanopore = 1;
      break;      
    case 'f':
      paraInfo.fold = 1;
      break;
    case 'o':
      outFile  = optarg;
      break;
    case 'a':
      faFile  = optarg;
      break;
    case 'A':
      faiFile = optarg;
      break;
    case 'i':
      inputFile  = optarg;
      break;
    case 'p':
      paraInfo.pvalue = atof(optarg);
      break;
    case 't':
      paraInfo.minTag = atof(optarg);
      break;
    case 'c':
      paraInfo.minCov = atof(optarg);
      break;
    case 'F':
      paraInfo.mfold = atof(optarg);
      break;
    case 'e':
      paraInfo.sefold = atof(optarg);
      break;
    case 'l':
      paraInfo.minLen = atoi(optarg);
      break;
    case 'L':
      paraInfo.maxLen = atoi(optarg);
      break;
    case 'm':
      paraInfo.minClusterLen = atoi(optarg);
      break;
    case 'M':
      paraInfo.maxClusterLen = atoi(optarg);
      break;
    case '?':
      showHelp = 1;
      break;
    default:
      usage();
    }
  }

  // help for version
  if (showVersion)
  {
    fprintf(stderr, "%s", version);
    exit(1);
  }

  if (showHelp)
  {
    usage();
    exit(1);
  }

  if (inputFile == NULL)
  {
    fprintf(stderr, "ERROR: please set the option: --input <bam alignments>\n");
    usage();
  }

  if (faFile != NULL)
  {
    genomefp = (FILE *) fopen(faFile, "r");
    if (genomefp == NULL)
    {
      fprintf(stderr, "ERROR: Can't open genome file: %s\n", faFile);
      usage();
    }
  }
  else
  {
    fprintf(stderr, "ERROR: please set the option: --fa <genome file>\n");
    usage();
  }

  if (faiFile != NULL)
  {
    faifp = (FILE *) fopen(faiFile, "r");
    if (faifp == NULL)
    {
      fprintf(stderr, "ERROR: Can't open genome file: %s\n", faiFile);
      usage();
    }
  }
  else
  {
    fprintf(stderr, "ERROR: please set the option: --fai <fai file>\n");
    usage();
  }

  if (outFile == NULL)
  {
    outfp = stdout;
  }
  else
  {
    outfp = (FILE *) fopen(outFile, "w");
    if (outfp == NULL)
    {
      fprintf(stderr, "ERROR: Can't open %s\n", outFile);
      usage();
    }
  }

  runnapSeeker(&paraInfo, outfp, genomefp, faifp, inputFile);

  fclose(genomefp);
  fclose(faifp);
  fclose(outfp);
  return 0;
}

void usage(void)
{
  fprintf(stderr, "%s", "Usage:  napSeeker [options] --fa <genome seq> --fai <fai file> --input <BAM alignments>\n\
napSeeker: for discovering novel napRNAs from NAP-seq data\n\
[options]\n\
-v/--verbose                   : verbose information\n\
-V/--version                   : napSeeker version\n\
-h/--help                      : help informations\n\
-n/--norm                      : normalized reads to the locus number\n\
-N/--nano                      : nanopore sequencing\n\
-f/--fold                      : rnafold the sequence with length<500nt\n\
-k/--keep-dup                  : keep duplication, deault is false\n\
-P/--pair                      : input is paired-end format\n\
-s/--strand                    : strand-specific\n\
-i/--input <string>            : input file<BAM format>\n\
-o/--outfile <string>          : output file\n\
-t/--min-tag <double>          : minimum tag number for each start or end pos, default>=5.0 reads\n\
-c/--cov <double>              : minimum coverage in each contig, default>=2 reads\n\
-r/--rpm <double>              : minimum rpm value for each psi, default>=0.05\n\
-l/--min-len <int>             : minimum length of reads, default=15\n\
-L/--max-len <int>             : maximum length of reads, default=1000\n\
-m <int>                       : minimum cluster length default=20\n\
-M <int>                       : maximum cluster length, default=10000\n\
-F <double>                    : minimum fold-change between start/end and mean reads[default=2]\n\
-e <double>                    : maximum fold-change between start and end reads[default=5]\n\
");
  exit(1);
}
