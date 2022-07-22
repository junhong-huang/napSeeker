/***********************************************************************
* napSeeker: a software for discovering novel napRNAs from NAP-seq data
*$2021/9/09/$ @Jian-Hua Yang yangjh7@mail.sysu.edu.cn
************************************************************************/
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<ctype.h>
#include<assert.h>
#include<math.h>
extern "C" {
#include "fold.h"
#include "fold_vars.h"
#include "utils.h"
#include "pair_mat.h"
}
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

#define CBOXLENGTH 7
#define DBOXLENGTH 4
#define BASES 4

/* suffix tree to improving speed */
double CboxSuffixMatrix[CBOXLENGTH - 1] = {10.916, 9.245, 6.988, 5.256, 3.544, 1.307};

/* Cbox matrix model is PSSM */
/* Note: -6.000 for escaping N or -(gap) in the pairwise sequence */
double CboxScores[CBOXLENGTH][BASES + 1] =
{ /*     A       C       G      T/U      N/- */
  /* 1 */{ 1.122, -2.116, 0.536, -4.448, -6.000},
  /* 2 */{ -5.448, -4.923, -2.116, 1.671, -6.000},
  /* 3 */{ -5.448, -4.923, 2.257, -5.448, -6.000},
  /* 4 */{ 1.732, -4.923, -4.923, -5.448, -6.000},
  /* 5 */{ -4.448, -4.923, -3.923, 1.712, -6.000},
  /* 6 */{ -5.448, -4.923, 2.237, -3.863, -6.000},
  /* 7 */{ 1.307, -3.338, -1.338, -0.804, -6.000}
};

void runnapSeeker(struct parameterInfo *paraInfo, FILE *outfp, FILE *genomefp,
                   FILE *faifp, char *inputBamFile)
{
  double totalReadNum  = 0;
  long   genomeSize    = 1;
  int skipDeletion     = 0;
  int skipSplice       = 2;
  int collapser        = 0;
  int skipBigClip      = 0;

  if (paraInfo->nanopore) skipBigClip = 1;

  chromBed12Map inputBedHash;
  chromBed12Map newInputBedHash;
  map<string, int> mapSize;

  char delims[]   = ",";
  char **bamFiles = NULL;
  int infoNum     = 0;

  faidxMap faiHash;
  fprintf(stderr, "#read genome fai file\n");
  readFai(faifp, faiHash);

  bamFiles = splitString(inputBamFile, delims, &infoNum);

  for (int i = 0; i < infoNum; i++)
  {
    char *bamFile = bamFiles[i];
    BamReader inputReader;
    if (!inputReader.Open(bamFile))
    {
      cerr << "Failed to open BAM file " << bamFile << endl;
      exit(1);
    }
    paraInfo->genomeSize = getGenomeSize(inputReader, mapSize);
    fprintf(stderr, "#read %s file to bed list\n", bamFile);
    if (paraInfo->pairEnd)
    {
      totalReadNum += readPEBamToBed12Map(inputReader, inputBedHash, skipDeletion, skipSplice, collapser, skipBigClip);
    }
    else
    {
      totalReadNum += readSEBamToBed12Map(inputReader, inputBedHash, skipDeletion, skipSplice, collapser, skipBigClip);
    }
    inputReader.Close();
  }

  fprintf(stderr, "#readNum=%.f\n", totalReadNum);
  fprintf(stderr, "#remove mispriming reads from bed list\n");
  totalReadNum = removeMisprimingReads(paraInfo, genomefp, faiHash, inputBedHash, newInputBedHash);
  fprintf(stderr, "#after removing the mispriming reads and remain reads=%.f\n", totalReadNum);

  if (paraInfo->norm)
  {
    //fprintf(stderr, "#normalized reads to their locus number\n");
    //totalReadNum = normalizedReadsByNH(newInputBedHash);
    fprintf(stderr, "#normalize reads\n");
    totalReadNum = normalizedBed12Reads(newInputBedHash);
  }

// output the head line
  fprintf(outfp, "#chrom\tstart\tend\tname\tscore\tstrand\t");
  fprintf(outfp, "contigLen\tclusterLen\treadNum\t");
  fprintf(outfp, "startReadNum\tendReadNum\t");
  fprintf(outfp, "startPvalue\tendPvalue\t");
  fprintf(outfp, "startFold\tendFold\t");
  fprintf(outfp, "upFold\tdownFold\t");
  fprintf(outfp, "up20ntFold\tdown20ntFold\t");
  fprintf(outfp, "clusterSeq\t");
  fprintf(outfp, "cdScore\t");
  fprintf(outfp, "structure\tshape\tmfe\n");
  paraInfo->totalReadNum = totalReadNum;

  fprintf(stderr, "#find non-cap RNAs\n");
  findNcapRNAs(paraInfo, outfp, mapSize, genomefp, faiHash, newInputBedHash);

  if (paraInfo->verbose)  fprintf(stderr, "#free bedHash...\n");
  freeChromBed12Map(inputBedHash);
}

double removeMisprimingReads(struct parameterInfo *paraInfo,
                             FILE *genomefp,
                             faidxMap &faiHash,
                             chromBed12Map &bedHash,
                             chromBed12Map &newBedHash)
{
  chromBed12Map::iterator it;
  int i = 0;
  int primerNum = 0;
  double totalNum = 0;
  char delims[]   = ":";
  char **infos    = NULL;
  int infoNum     = 0;

  for (it = bedHash.begin(); it != bedHash.end(); ++it)
  {
    bed12Vector bedList = it->second;
    char *chrom = (char *)it->first.c_str();
    string chromStr(chrom);
    if (faiHash.find(chromStr) == faiHash.end())
    {
      fprintf(stderr, "can't not find the chromosome %s, skip it.\n", chrom);
      continue;
    }
    faidx *fai = faiHash[chromStr];
    if (bedList.size() >= 1)
    {
      for (bed12Vector::iterator vecItr = bedList.begin(); vecItr != bedList.end(); vecItr++)
      {
        CBed12 *bed = *vecItr;

        infos = splitString(bed->name, delims, &infoNum);
        char *barcodeSeq = strClone(infos[infoNum - 2]);
        int adapterFlag = atoi(infos[infoNum - 1]);
        freeWords(infos, infoNum);
        if (adapterFlag != 3 && adapterFlag != 8) // not with 3'-adapter
        {
          newBedHash[chromStr].push_back(bed);
          totalNum += bed->score;
        }
        else //  with 3'-adapter
        {
          int barcodeLen = strlen(barcodeSeq);
          int start = bed->chromStart;
          int end = bed->chromEnd;
          int primerStart = end;
          int primerEnd = end + barcodeLen;
          if (bed->strand == '-')
          {
            primerEnd = start;
            primerStart = start - barcodeLen;
          }
          if (primerStart < 0) primerStart = 0;
          int matchNum = 0;
          char *primerSeq = faidxFetchSeq(genomefp, fai, primerStart, primerEnd, bed->strand);
          for (i = 0; i < barcodeLen; i++)
          {
            if (primerSeq[i] == barcodeSeq[i])
            {
              matchNum++;
            }
          }
          if (matchNum >= barcodeLen - 1)
          {
            primerNum++;
            if (paraInfo->verbose) fprintf(stderr, "mispriming read%d: %s %c\n", primerNum, bed->name, bed->strand);
          }
          else {
            newBedHash[chromStr].push_back(bed);
            totalNum += bed->score;
          }
          safeFree(primerSeq);
        }// else end
        safeFree(barcodeSeq);
      } // for end
    }// if sizes
  } // for bed hash
  fprintf(stderr, "#remove the mispriming reads=%d\n", primerNum);
  return totalNum;
}

void findNcapRNAs(struct parameterInfo *paraInfo, FILE *outfp,
                  map<string, int> &mapSize, FILE *genomefp,
                  faidxMap &faiHash, chromBed12Map &bedHash)
{
  chromBed12Map::iterator it;
  char strand = '.';
  int readNum = 0;
  int ncapNum = 1;
  for (it = bedHash.begin(); it != bedHash.end(); ++it)
  {
    bed12Vector bedList = it->second;
    char *chrom = (char *)it->first.c_str();
    string chromStr(chrom);
    if (faiHash.find(chromStr) == faiHash.end())
    {
      fprintf(stderr, "can't not find the chromosome %s, skip it.\n", chrom);
      continue;
    }
    faidx *fai = faiHash[chromStr];
    int geneLen = mapSize[it->first];
    if (bedList.size() >= 1)
    {
      if (paraInfo->strand)
      {
        strand = '+';
        if (strand == '+')
        {
          profileInfo *readProfile = (profileInfo *)safeMalloc(sizeof(struct profileInfo));
          fprintf(stderr, "#get read values from chrom %s: %c\n", chrom, strand);
          readNum = getReadVals(paraInfo, bedList, readProfile, geneLen, strand);
          ncapNum = getBedReadVals(paraInfo, outfp, genomefp, fai, chrom, strand, geneLen, readProfile, ncapNum);
          freeProfiles(readProfile);
        }
        strand = '-';
        if (strand == '-')
        {
          profileInfo *readProfile = (profileInfo *)safeMalloc(sizeof(struct profileInfo));
          fprintf(stderr, "#get read values from chrom %s: %c\n", chrom, strand);
          readNum = getReadVals(paraInfo, bedList, readProfile, geneLen, strand);
          ncapNum = getBedReadVals(paraInfo, outfp, genomefp, fai, chrom, strand, geneLen, readProfile, ncapNum);
          freeProfiles(readProfile);
        }
      } // if strand
    }// if sizes
  } // for bed hash
}

int getReadVals(struct parameterInfo *paraInfo, bed12Vector &bedList,
                profileInfo *profile, int geneLen, char strand)
{
  int i = 0;
  int tagNum = 0;
  char delims[]   = ":";
  char **infos    = NULL;
  int infoNum     = 0;
  profile->startProfile  = (double *)safeMalloc(sizeof(double) * geneLen);
  profile->endProfile    = (double *)safeMalloc(sizeof(double) * geneLen);
  profile->heightProfile = (double *)safeMalloc(sizeof(double) * geneLen);
  for (i = 0; i < geneLen; i++)
  {
    profile->startProfile[i]  = 0;
    profile->endProfile[i]    = 0;
    profile->heightProfile[i] = 0;
  }
  fprintf(stderr, "# safeZeroedMalloc for chrom length %d\n", geneLen);
  for (bed12Vector::iterator vecItr = bedList.begin(); vecItr != bedList.end(); vecItr++)
  {
    CBed12 *bed = *vecItr;
    int start = bed->chromStart;
    int end = bed->chromEnd;
    int readLen = end - start;
    if (readLen < paraInfo->minLen) continue;
    int startIdx = start;
    int endIdx = end - 1;
    if (bed->strand == '-')
    {
      startIdx = end - 1;
      endIdx = start;
    }
    if (bed->strand == strand)
    {
      infos = splitString(bed->name, delims, &infoNum);
      int adapterFlag = atoi(infos[infoNum - 1]);
      if (paraInfo->nanopore)
      {
        if (adapterFlag == 8)
        {
          profile->startProfile[startIdx] += bed->score;
          profile->endProfile[endIdx] += bed->score;
        }
      }
      else
      {
        if (adapterFlag == 5 || adapterFlag == 8)
          profile->startProfile[startIdx] += bed->score;
        if (adapterFlag == 3 || adapterFlag == 8)
          profile->endProfile[endIdx] += bed->score;
      }
      freeWords(infos, infoNum);
    }
    for (i = 0; i < bed->blockCount; i++) // modified and remove the 1
    {
      int bStart = bed->chromStart + bed->chromStarts[i];
      int bEnd   = bStart + bed->blockSizes[i];
      int j = 0;
      for (j = bStart; j < bEnd; j++)
      {
        profile->heightProfile[j] += bed->score;
      }
    }
    tagNum++;
  } // for end
  return tagNum;
}

int getBedReadVals(struct parameterInfo * paraInfo, FILE * outfp, FILE * genomefp,
                   faidx * fai, char *chrom, char strand,
                   int geneLen, profileInfo * readProfile, int ncapNum)
{
  int i = 0;
  int contigStart = 0;
  int contigEnd = 0;
  int startFlag = 0;
  int contigNum = 0;
  double *hProfile  = readProfile->heightProfile;
  double *sProfile  = readProfile->startProfile;
  double *eProfile  = readProfile->endProfile;
  for (i = 0; i < geneLen; i++)
  {
    double hNum = hProfile[i];
    if (hNum != 0 && startFlag == 0)
    {
      contigStart = i;
      startFlag = 1;
      continue;
    }
    if (hNum  == 0 && startFlag == 1)
    {
      contigEnd = i;
      contigNum = findContig(paraInfo, outfp, genomefp, fai, chrom, strand, geneLen, readProfile, contigStart, contigEnd);
      startFlag = 0;
    }
  }
  return contigNum;
}

int findContig(struct parameterInfo *paraInfo, FILE *outfp,
               FILE *genomefp, faidx *fai, char *chrom,
               char strand, int geneLen, profileInfo *readProfile,
               int start, int end)
{
  int i = 0;
  vector<intDbl> startIntDblVec;
  vector<intDbl> endIntDblVec;
  clusterVector candidateVector;
  int contigLen = end - start + 1;
  int contigNum = 0;
  double *hProfile  = readProfile->heightProfile;
  double *sProfile  = readProfile->startProfile;
  double *eProfile  = readProfile->endProfile;
  double startTotalNum = 0;
  double endTotalNum = 0;
  double maxSEFold = paraInfo->sefold; // start num/end number
  for (i = start; i < end; i++)
  {
    double sNum = sProfile[i];
    double eNum = eProfile[i];
    startTotalNum += sNum;
    endTotalNum   += eNum;
    if (sNum >= paraInfo->minTag)
    {
      startIntDblVec.push_back(make_pair(i, sNum));
    }
    if (eNum >= paraInfo->minTag)
    {
      endIntDblVec.push_back(make_pair(i, eNum));
    }
  }
  sort(startIntDblVec.begin(), startIntDblVec.end(), compIntDbl);
  sort(endIntDblVec.begin(), endIntDblVec.end(), compIntDbl);
  for (vector<intDbl>::iterator sit = startIntDblVec.begin(); sit != startIntDblVec.end(); ++sit )
  {
    int startPos = sit->first;
    double startNum = sit->second;
    for (vector<intDbl>::iterator eit = endIntDblVec.begin(); eit != endIntDblVec.end(); ++eit )
    {
      int chromStart = 0;
      int chromEnd = 0;
      int endPos = eit->first;
      double endNum = eit->second;
      if (startNum / endNum > maxSEFold || startNum / endNum < 1 / maxSEFold) continue;
      if (strand == '+')
      {
        chromStart = startPos;
        chromEnd = endPos + 1;
      }
      else if (strand == '-')
      {
        chromStart = endPos;
        chromEnd = startPos + 1;
      }
      int readLen = chromEnd - chromStart;
      if (readLen < paraInfo->minClusterLen || readLen > paraInfo->maxClusterLen) continue;
      double minStartCov = hProfile[chromStart] / maxSEFold;
      double minEndCov   = hProfile[chromEnd] / maxSEFold;
      if (getMinCoverage(hProfile, chromStart, chromEnd, paraInfo->minCov, minStartCov, minEndCov)) continue;
      double peReads = startNum + endNum;
      int overlapTag = getOverlap(chromStart, chromEnd, peReads, candidateVector);
      if (overlapTag == 0)
      {
        double startFold      = getEndsFold(sProfile, geneLen, startPos);
        double endFold        = getEndsFold(eProfile, geneLen, endPos);
        double meanCov        = getMeanCoverage(hProfile, chromStart, chromEnd, readLen);
        double upFold         = getUpCovFold(hProfile, geneLen, chromStart, chromEnd, 0, meanCov);
        double downFold       = getDownCovFold(hProfile, geneLen, chromStart, chromEnd, 0, meanCov);
        double up20ntFold     = getUpCovFold(hProfile, geneLen, chromStart, chromEnd, 20, meanCov);
        double down20ntFold   = getDownCovFold(hProfile, geneLen, chromStart, chromEnd, 20, meanCov);
        double p = 1 / (double)contigLen;
        double startPval = logbinomial((int)round(startTotalNum), (int)round(startNum), p, contigLen);
        double endPval = logbinomial((int)round(endTotalNum), (int)round(endNum), p, contigLen);
        if (startFold >= paraInfo->mfold
            && endFold >= paraInfo->mfold
            && upFold >= paraInfo->mfold
            && downFold >= paraInfo->mfold
            && up20ntFold >= paraInfo->mfold
            && down20ntFold >= paraInfo->mfold)
        {
          Cluster *pCluster = (Cluster *)safeMalloc(sizeof(Cluster));
          pCluster->chrom = strClone(chrom);
          pCluster->chromStart = chromStart;
          pCluster->chromEnd   = chromEnd;
          pCluster->strand     = strand;
          pCluster->readNum    = peReads;
          pCluster->startReadNum  = startNum;
          pCluster->endReadNum    = endNum;
          pCluster->startPvalue   = startPval;
          pCluster->endPvalue     = endPval;
          pCluster->peakStart  = start;
          pCluster->peakEnd    = end;
          pCluster->startFold  = startFold;
          pCluster->endFold    = endFold;
          pCluster->upFold     = upFold;
          pCluster->downFold   = downFold;
          pCluster->up20ntFold    = up20ntFold;
          pCluster->down20ntFold  = down20ntFold;
          candidateVector.push_back(pCluster);
        }
      }
    }
  }
  if (candidateVector.size() > 0)
  {
    contigNum = outputClusterInfo(paraInfo, outfp, genomefp, fai, candidateVector);
    freeClusterVector(candidateVector);
  }
  return contigNum;
}

int outputClusterInfo(struct parameterInfo *paraInfo, FILE *outfp, FILE *gfp, faidx *fai, clusterVector &contigVector)
{
  int clusterNum = 1;
  for (clusterVector::iterator vecItr = contigVector.begin(); vecItr != contigVector.end(); vecItr++)
  {
    Cluster *pCluster = *vecItr;
    if (pCluster->readNum >= paraInfo->minTag)
    {
      char *clusterSeq = faidxFetchSeq(gfp, fai, pCluster->chromStart, pCluster->chromEnd, pCluster->strand);
      double cdScore = scoreCDbox(clusterSeq, paraInfo);
      int contigLen = pCluster->chromEnd - pCluster->chromStart;
      int peakLen = pCluster->peakEnd - pCluster->peakStart;
      if (paraInfo->fold)
      {
        fprintf(outfp, "%s\t%d\t%d\tnapSeeker_%s_%d_%d\t%.0f\t%c\t", pCluster->chrom, pCluster->chromStart, pCluster->chromEnd, pCluster->chrom, pCluster->chromStart, pCluster->chromEnd, pCluster->readNum, pCluster->strand);
        fprintf(outfp, "%d\t%d\t%.5f\t", contigLen, peakLen, pCluster->readNum);
        fprintf(outfp, "%.5f\t%.5f\t", pCluster->startReadNum, pCluster->endReadNum);
        fprintf(outfp, "%g\t%g\t", pCluster->startPvalue, pCluster->endPvalue);
        fprintf(outfp, "%.5f\t%.5f\t", pCluster->startFold, pCluster->endFold);
        fprintf(outfp, "%.5f\t%.5f\t", pCluster->upFold, pCluster->downFold);
        fprintf(outfp, "%.5f\t%.5f\t", pCluster->up20ntFold, pCluster->down20ntFold);
        fprintf(outfp, "%s\t%.5f\t", clusterSeq, cdScore);
        if (contigLen < 500)
        {
          char *structure = (char *)safeMalloc(strlen(clusterSeq) + 1);
          double mfe = fold(clusterSeq, structure);
          char *shape = drawShapes(structure);
          fprintf(outfp, "%s\t%s\t%.3f\n", structure, shape, mfe);
          free_arrays;
          safeFree(shape);
          safeFree(structure);
        }
        else
        {
          fprintf(outfp, ".\t.\t0\n");
        }
      }
      fflush(outfp);
      safeFree(clusterSeq);
      clusterNum++;
    }
  }
  return (clusterNum - 1);
}

int encodeInt(char ch)
{
  /* translate character to number */

  ch = toupper(ch);

  if (ch == 'A') return 0;
  else if (ch == 'C') return 1;
  else if (ch == 'G') return 2;
  else if (ch == 'U' || ch == 'T') return 3;
  else return 4;
}

double scoreCDbox(char *seq, struct parameterInfo *paraInfo)
{
  int i, j;
  int cboxStart = 4;
  int cboxEnd   = 7;
  int cboxLen   = 7;
  double tmpScore  = 0;
  double bestCboxScore  = 0;
  double minCboxScore = 3.0; /* minimum score 3.825 for box C (U28 snoRNA) */
  int seqLen = strlen(seq);
  int boxDtag = 0;
  double DboxScore = 0;
  if (seqLen < 18) return 0;
  i = seqLen - 4; // for CTGA
  if (boxDtag != 1 && seq[i] == 'C' && seq[i + 1] == 'T' && seq[i + 2] == 'G' && seq[i + 3] == 'A')
  {
    boxDtag = 1;
    DboxScore = 7.967;
  }
  i = seqLen - 6; // for CTGANN
  if (boxDtag != 1 && seq[i] == 'C' && seq[i + 1] == 'T' && seq[i + 2] == 'G' && seq[i + 3] == 'A')
  {
    boxDtag = 1;
    DboxScore = 7.967;
  }
  i = seqLen - 9; // for CTGANNNNN
  if (boxDtag != 1 && seq[i] == 'C' && seq[i + 1] == 'T' && seq[i + 2] == 'G' && seq[i + 3] == 'A')
  {
    boxDtag = 1;
    DboxScore = 7.967;
  }
  i = seqLen - 7; // for CTGANNN
  if (boxDtag != 1 && seq[i] == 'C' && seq[i + 1] == 'T' && seq[i + 2] == 'G' && seq[i + 3] == 'A')
  {
    boxDtag = 1;
    DboxScore = 7.967;
  }
  i = seqLen - 8; // for CTGANNNN
  if (boxDtag != 1 && seq[i] == 'C' && seq[i + 1] == 'T' && seq[i + 2] == 'G' && seq[i + 3] == 'A')
  {
    boxDtag = 1;
    DboxScore = 7.967;
  }
  // for ATGA D-Box
  if (boxDtag != 1)
  {
    i = seqLen - 6; // for ATGANN
    if (seq[i] == 'A' && seq[i + 1] == 'T' && seq[i + 2] == 'G' && seq[i + 3] == 'A')
    {
      boxDtag = 2;
      DboxScore = 1.272;
    }
    i = seqLen - 9; // for ATGANNNNN
    if (seq[i] == 'A' && seq[i + 1] == 'T' && seq[i + 2] == 'G' && seq[i + 3] == 'A')
    {
      boxDtag = 2;
      DboxScore = 1.272;
    }
  }

  if (boxDtag != 1)
  {
    return 0;
  }

  for (i = cboxStart; i < cboxEnd; i++)
  {
    tmpScore = 0;
    char *tmpSeq = seq + i;
    for (j = 0; j < cboxLen; j++)
    {
      tmpScore += CboxScores[j][(int)(encodeInt(tmpSeq[j]))];
      /* suffix tree for box C  */
      if (j < (cboxLen - 1))
      {
        if ((tmpScore + CboxSuffixMatrix[j]) <= minCboxScore)
        {
          break;
        } // if score
      }// if boxlen
    }// for j
    if (tmpScore > bestCboxScore)
    {
      bestCboxScore = tmpScore;
    }
  }// for i

  if (bestCboxScore > minCboxScore)
  {
    return (bestCboxScore + DboxScore);
  }
  else
  {
    return 0;
  }
}

void freeClusterVector(clusterVector &contigVector)
{
  for (clusterVector::iterator vecItr = contigVector.begin(); vecItr != contigVector.end(); vecItr++)
  {
    Cluster *pCluster = *vecItr;
    freeCluster(pCluster);
  }
  contigVector.clear();
}

void freeCluster(Cluster *pCluster)
{
  safeFree(pCluster->chrom);
  safeFree(pCluster);
}

double getEndsFold(double *posProfile, int arrayNum, int bestPos)
{
  int i, j;
  double readNum = 0;
  int posNum = 0;
  double std = 0;
  double meanVal = 0;
  int extendLen = 100;
  int start = bestPos - extendLen;
  int end = bestPos + extendLen;
  if (start < 0) start = 0;
  if (end > arrayNum) end = arrayNum;
  for (i = start; i < end; i++)
  {
    if (posProfile[i] != 0 && i != bestPos)
    {
      readNum += posProfile[i];
      posNum  += 1;
    }
  }
  if (readNum == 0)
  {
    readNum = 0.5;
    posNum = 1;
  }
  meanVal = readNum / double(posNum);
  std = posProfile[bestPos] / meanVal;
  return std;
}

double getUpCovFold(double *covProfile, int arrayNum, int start, int end, int extendLen, double regMeanCov)
{
  int i;
  double meanVal = 0;
  if (extendLen == 0)
  {
    extendLen = end - start;
  }
  int upStart = start - extendLen;
  int upEnd = start;
  if (upStart < 0) upStart = 0;
  meanVal = getMeanCoverage(covProfile, upStart, upEnd, extendLen);
  double mfold = regMeanCov / meanVal;
  return mfold;
}

double getDownCovFold(double *covProfile, int arrayNum, int start, int end, int extendLen, double regMeanCov)
{
  int i;
  double meanVal = 0;
  if (extendLen == 0)
  {
    extendLen = end - start;
  }
  int downStart = end;
  int downEnd = end + extendLen;
  if (downEnd > arrayNum) downEnd = arrayNum;
  meanVal = getMeanCoverage(covProfile, downStart, downEnd, extendLen);
  double mfold = regMeanCov / (double)meanVal;
  return mfold;
}

int getMinCoverage(double *covProfile, int start, int end, double minCov, double minStartCov, double minEndCov)
{
  int i;
  double retVal = 0;
  for (i = start; i < end; i++)
  {
    if (covProfile[i] < minCov || covProfile[i] < minStartCov || covProfile[i] < minEndCov)
    {
      retVal = 1;
      break;
    }
  }
  return retVal;
}

double getMeanCoverage(double *covProfile, int start, int end, int posNum)
{
  int i;
  double mean = 0;
  double readNum = 0;
  for (i = start; i < end; i++)
  {
    readNum += covProfile[i];
  }
  if (readNum == 0)
  {
    readNum = 0.5;
  }
  if (posNum == 0)
  {
    posNum = 1;
  }
  mean = (double)readNum / (double)posNum;
  return mean;
}

int getOverlap(int start, int end, double peReads, clusterVector &candidateVector)
{
  int overlapTag = 0;
  for (clusterVector::iterator vecItr = candidateVector.begin(); vecItr != candidateVector.end();)
  {
    Cluster *pCluster = *vecItr;
    int oLen = overlapLength(start, end, pCluster->chromStart, pCluster->chromEnd);
    double seReads = pCluster->startReadNum + pCluster->endReadNum;
    if (oLen > 0)
    {
      if (peReads > seReads)
      {
        freeCluster(pCluster);
        vecItr = candidateVector.erase(vecItr);
        overlapTag = 0;
      }
      else {
        overlapTag = 1;
        break;
      }
    }
    else {
      vecItr++;
    }
  }
  return overlapTag;
}

void freeProfiles(profileInfo * profile)
{
  safeFree(profile->startProfile);
  safeFree(profile->endProfile);
  safeFree(profile->heightProfile);
  safeFree(profile);
}

double round(double r)
{
  return (r > 0.0) ? floor(r + 0.5) : ceil(r - 0.5);
}

bool compIntDbl(intDbl a, intDbl b)
{
  return (a.second > b.second);
}

char *drawShapes(char *structure)
{
  int i = 0;
  int j = 0;
  int hx = 0;
  int pp = 0;
  int cp = 0; // current postion
  char *shapes = (char *)safeMalloc(sizeof(char) * strlen(structure) + 1);
  short int *pairTable = make_pair_table(structure);
  int *blacketTable = (int *)safeMalloc(sizeof(int) * strlen(structure) + 1);
  int structLen = strlen(structure);
  blacketTable[0] = 0;
  for (i = 0; i < structLen; i++)
  {
    if (structure[i] == '(') hx += 1;
    blacketTable[i + 1] = hx;
  }
  for (i = 1; i <= structLen; i++)
  {
    cp = i - 1;
    if (pairTable[i] > 0)
    {
      if (structure[pp] == '(' && structure[cp] == ')')
      {
        shapes[j] = '[';
        shapes[++j] = ']';
        ++j;
      }
      if (abs(blacketTable[pairTable[i]] - blacketTable[pairTable[pp + 1]]) > 1
          && (structure[pp] == structure[cp]))
      {
        if (structure[cp] == '(')
        {
          shapes[j] = '[';
          ++j;
        }
        else
        {
          shapes[j] = ']';
          ++j;
        }
      } // if abs
      pp = cp; // previous postion
    }
  }// for structLen
  shapes[j] = '\0';
  if (strlen(shapes) < 1)
  {
    shapes[0] = '.';
    shapes[1] = '\0';
  }
  safeFree(pairTable);
  safeFree(blacketTable);
  return shapes;
}
