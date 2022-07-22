/***********************************************************************
* napSeeker: a software for discovering novel napRNAs from NAP-seq data
*$2021/9/09/$ @Jian-Hua Yang yangjh7@mail.sysu.edu.cn
************************************************************************/
#ifndef napSeeker_HEAD_H
#define napSeeker_HEAD_H

#define pseudCount 1.0

typedef std::pair<string, double> strDbl;
typedef std::pair<int, double> intDbl;

typedef struct parameterInfo
{
  int    verbose;
  int    keepDup;
  int    pairEnd;
  int    minLen;
  int    maxLen;
  int    strand;
  int    split;
  int    norm;
  int    fold;
  int    insertSize;
  int    minClusterLen;
  int    maxClusterLen;
  int    nanopore;
  double minTag;
  double pvalue;
  double fdr;
  double totalReadNum;
  double rpm;
  double minCov;
  double mfold;
  double sefold;
  long   genomeSize;
} parameterInfo;

typedef struct profileInfo
{
  double *startProfile;
  double *endProfile;
  double *heightProfile;
} profileInfo;

struct peakInfo
{
  char *chrom;  // chromosome name
  int chromStart; // chromosome start
  int chromEnd; // chromosome end
  char strand; // strand
  double readNum;
  double startReadNum;
  double endReadNum;

  double startPvalue;
  double endPvalue;

  int peakStart;
  int peakEnd;

  double startFold;
  double endFold;
  double upFold;
  double downFold;
  double up20ntFold;
  double down20ntFold;
};

typedef struct peakInfo Cluster;

typedef vector<Cluster *> clusterVector;

void runnapSeeker(struct parameterInfo *paraInfo, FILE *outfp, FILE *genomefp,
                   FILE *faifp, char *inputBamFile);

double removeMisprimingReads(struct parameterInfo *paraInfo,
                             FILE *genomefp,
                             faidxMap &faiHash,
                             chromBed12Map &bedHash,
                             chromBed12Map &newBedHash);

void findNcapRNAs(struct parameterInfo *paraInfo, FILE *outfp,
                  map<string, int> &mapSize, FILE *genomefp,
                  faidxMap &faiHash, chromBed12Map &bedHash);

int getReadVals(struct parameterInfo *paraInfo, bed12Vector &bedList,
                profileInfo *profile, int geneLen, char strand);

int getBedReadVals(struct parameterInfo * paraInfo, FILE * outfp, FILE * genomefp,
                   faidx * fai, char *chrom, char strand,
                   int geneLen, profileInfo * readProfile, int ncapNum);

int findContig(struct parameterInfo *paraInfo, FILE *outfp,
               FILE *genomefp, faidx *fai, char *chrom,
               char strand, int geneLen, profileInfo *readProfile,
               int start, int end);

int outputClusterInfo(struct parameterInfo *paraInfo, FILE *outfp, FILE *gfp, faidx *fai, clusterVector &contigVector);

void freeClusterVector(clusterVector &contigVector);

void freeCluster(Cluster *pCluster);

double getEndsFold(double *posProfile, int arrayNum, int bestPos);

double getUpCovFold(double *covProfile, int arrayNum, int start, int end, int extendLen, double regMeanCov);

double getDownCovFold(double *covProfile, int arrayNum, int start, int end, int extendLen, double regMeanCov);

double getMeanCoverage(double *covProfile, int start, int end, int posNum);

int getMinCoverage(double *covProfile, int start, int end, double minCov, double minStartCov, double minEndCov);

char *drawShapes(char *structure);

double scoreCDbox(char *seq, struct parameterInfo *paraInfo);

int getOverlap(int start, int end, double peReads, clusterVector &candidateVector);

void freeProfiles(profileInfo *profile);

bool compIntDbl(intDbl a, intDbl b);

double round(double r);

int encodeInt(char ch);


#endif /* End napSeeker_HEAD_H */
