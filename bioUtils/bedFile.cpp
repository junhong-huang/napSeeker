#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<ctype.h>
#include<assert.h>
#include<math.h>
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

bool compareBed( const CBedInfo *oBed, const CBedInfo *tBed)
{
  return (oBed->chromStart < tBed->chromStart);
}

bool compareBed6( const CBed6Info *oBed, const CBed6Info *tBed)
{
  return (oBed->chromStart < tBed->chromStart);
}

int compareBedScore(const void *a, const void *b)
{
  CBed *oBed = *(struct CBedInfo **)a;
  CBed *tBed = *(struct CBedInfo **)b;
  if ((tBed->score - oBed->score) > 0)
    return 1;
  else if ((tBed->score - oBed->score) < 0)
    return -1;
  else
    return 0;
}

void freeBedList(CBed *bedList)
/*free bed list */
{
  CBed *m = NULL;
  CBed *p = NULL;
  m = bedList;
  while (m != NULL)
  {
    p = m;
    m = p->next;
    freeBedItem(p);
  }
}

void freeBedItem(CBed *bed)
{
  safeFree(bed->chrom);
  safeFree(bed);
}

void freeBed6Item(CBed6 *bed6)
{
  safeFree(bed6->chrom);
  safeFree(bed6->name);
  safeFree(bed6);
}

void freeBed12Item(CBed12 *bed12)
{
  safeFree(bed12->chrom);
  safeFree(bed12->name);
  safeFree(bed12->blockSizes);
  safeFree(bed12->chromStarts);
  safeFree(bed12);
}


void freeBedMapList(map<string, CBed *> &bedHash) /*free bed list */
{
  for (map<string, CBed *>::iterator curr = bedHash.begin(); curr != bedHash.end();
       curr++)
  {
    CBed *p = curr->second;
    freeBedList(p);
  }
  bedHash.clear();
}

void freeNameBed12Map(nameBed12Map &bedHash) /*free bed name hash */
{
  for (nameBed12Map::iterator mapItr = bedHash.begin(); mapItr != bedHash.end(); mapItr++) {
    CBed12 *bed = mapItr->second;
    freeBed12Item(bed);
  }
  bedHash.clear();
}

void freeChromStrandMap(chromStrandBed12Map &bedHash) /*free bed list */
{
  for (chromStrandBed12Map::iterator mapItr = bedHash.begin(); mapItr != bedHash.end(); mapItr++) {
    bed12Vector bedList = mapItr->second;
    freeBed12Vector(bedList);
  }
  bedHash.clear();
}

void freeBed12Vector(bed12Vector &bedList)
{
  for (bed12Vector::iterator vecItr = bedList.begin(); vecItr != bedList.end(); vecItr++)
  {
    CBed12 *bed = *vecItr;
    freeBed12Item(bed);
  }
  bedList.clear();
}

void freeChromStrandBed6Map(chromStrandBed6Map &bedHash) /*free bed list */
{
  for (chromStrandBed6Map::iterator mapItr = bedHash.begin(); mapItr != bedHash.end(); mapItr++) {
    bed6Vector bedList = mapItr->second;
    freeBed6Vector(bedList);
  }
  bedHash.clear();
}

void freeBed6Vector(bed6Vector &bedList)
{
  for (bed6Vector::iterator vecItr = bedList.begin(); vecItr != bedList.end(); vecItr++)
  {
    CBed6 *bed = *vecItr;
    freeBed6Item(bed);
  }
  bedList.clear();
}

void copyBed6(CBed6 *tBed, CBed6 *oBed)
// copy oBed to tBed
{
  tBed->chrom       = strClone(oBed->chrom);
  tBed->name        = strClone(oBed->name);
  tBed->chromStart  = oBed->chromStart;
  tBed->chromEnd    = oBed->chromEnd;
  tBed->score       = oBed->score;
  tBed->strand      = oBed->strand;
}

void copyBed(CBed *tBed, CBed *oBed)
// copy oBed to tBed
{
  tBed->chrom      = strClone(oBed->chrom);
  tBed->chromStart = oBed->chromStart;
  tBed->chromEnd   = oBed->chromEnd;
  tBed->score    = oBed->score;
  tBed->strand     = oBed->strand;
}

void sortBed(CBed **pList, int count, int (*compare )(const void *item1,  const void *item2))
// sort bed, I have done major modified from jksrc common.c slSort, thanks Jim
{
  //fprintf(stderr, "sort bed array...\n");
  CBed *list = *pList;
  if (count > 1)
  {
    CBed *el;
    CBed **array;
    int i;
    array = (CBed **)safeMalloc(count * sizeof(CBed *));
    for (el = list, i = 0; el != NULL; el = el->next, i++)
      array[i] = el;
    qsort(array, count, sizeof(array[0]), compare);
    list = array[0];
    for (i = 1; i < count; ++i)
    {
      list->next = array[i];
      list = array[i];
    }
    list->next = NULL;
    *pList = array[0];
    safeFree(array);
  }
  //fprintf(stderr, "sort end...\n");
}

CBed12 *parseBed12Line(char *line)
{
  int i = 0;
  int fieldNum = 0;
  char **fields = NULL;
  int tmpFieldNum = 0;
  char **tmpFields = NULL;
  char delim[] = ",";
  CBed12 *bed = NULL;
  fields = splitWhitespace(line, &fieldNum);
  if (fieldNum != 12)
  {
    freeWords(fields, fieldNum);
    fprintf(stderr, "the format of annotation file must be bed12\n");
    exit(1);
  }
  bed = (CBed12 *)safeMalloc(sizeof(CBed12));
  bed->chrom = strClone(fields[0]);
  bed->chromStart = atoi(fields[1]);
  bed->chromEnd = atoi(fields[2]);
  bed->name = strClone(fields[3]);
  bed->score = atof(fields[4]);
  bed->strand = fields[5][0];
  bed->thickStart = atoi(fields[6]);
  bed->thickEnd = atoi(fields[7]);
  bed->itemRgb = 0;
  bed->blockCount = atoi(fields[9]);
  bed->blockSizes = (int *)safeMalloc(sizeof(int) * bed->blockCount);
  bed->chromStarts = (int *)safeMalloc(sizeof(int) * bed->blockCount);
  tmpFields = splitString(fields[10], delim, &tmpFieldNum);
  for (i = 0; i < bed->blockCount; i++)
  {
    bed->blockSizes[i] = atoi(tmpFields[i]);
  }
  freeWords(tmpFields, tmpFieldNum);
  tmpFields = splitString(fields[11], delim, &tmpFieldNum);
  for (i = 0; i < bed->blockCount; i++)
  {
    bed->chromStarts[i] = atoi(tmpFields[i]);
  }
  freeWords(tmpFields, tmpFieldNum);
  freeWords(fields, fieldNum);

  return bed;
}

double readBedToBed6Map(FILE *bedfp, chromBed6Map &bed6Hash)
{
  double totalNum = 0;
  char *line = NULL;
  CBed6 *bedPtr = NULL;
  while (line = getLine(bedfp))
  {
    if (feof(bedfp) || line == NULL)
    {
      safeFree(line);
      break;
    }
    bedPtr = parseBed6Line(line); // get bed information
    string chrom = bedPtr->chrom;
    bed6Hash[chrom].push_back(bedPtr);
    safeFree(line); // free memory
    totalNum += 1;
  }
  return totalNum;
}

int getBed12ToNameMap(FILE *bedfp, nameBed12Map &bed12NameHash)
{
  int bedNum = 0;
  char *line = NULL;
  CBed12 *bed12 = NULL;
  while (line = getLine(bedfp))
  {
    if (feof(bedfp) || line == NULL)
    {
      safeFree(line);
      break;
    }
    bed12 = parseBed12Line(line); // get bed information
    string name = bed12->name;
    bed12NameHash[name] = bed12;
    safeFree(line); // free memory
    bedNum++;
  }
  return bedNum;
}

int getBedToMap(FILE *bedfp, chromStrandBed12Map &bed12Hash)
{
  int bedNum = 0;
  char *line = NULL;
  CBed12 *bed12 = NULL;
  while (line = getLine(bedfp))
  {
    if (feof(bedfp) || line == NULL)
    {
      safeFree(line);
      break;
    }
    bed12 = parseBed12Line(line); // get bed information
    string chromStrand = bed12->chrom;
    chromStrand = chromStrand + bed12->strand;
    bed12Hash[chromStrand].push_back(bed12);
    safeFree(line); // free memory
    bedNum++;
  }
  return bedNum;
}

int getBed6ToMap(FILE *bedfp, chromStrandBed6Map &bed6Hash)
{
  int bedNum = 0;
  char *line = NULL;
  CBed6 *bed6 = NULL;
  while (line = getLine(bedfp))
  {
    if (feof(bedfp) || line == NULL)
    {
      safeFree(line);
      break;
    }
    bed6 = parseBed6Line(line); // get bed information
    string chromStrand = bed6->chrom;
    chromStrand = chromStrand + bed6->strand;
    bed6Hash[chromStrand].push_back(bed6);
    safeFree(line); // free memory
    bedNum++;
  }
  return bedNum;
}

long getChromSize(FILE *gfp, chromSizeMap &mapSize)
{
  long genomeSize = 0;
  char *line = NULL;
  int fieldNum = 0;
  char **fields = NULL;
  while (line = getLine(gfp))
  {
    if (feof(gfp) || line == NULL)
    {
      safeFree(line);
      break;
    }
    fields = splitWhitespace(line, &fieldNum);
    if (fieldNum < 2)
    {
      freeWords(fields, fieldNum);
      fprintf(stderr, "the format of chromosome size is<chromName><TAB><chromSize>\n");
      exit(1);
    }
    string refName = fields[0];
    int refLen = atoi(fields[1]);
    mapSize[refName] = refLen;
    genomeSize += refLen;
    safeFree(line); // free memory
    freeWords(fields, fieldNum);
  }
  return genomeSize;
}


CBed6 *parseBed6Line(char *line)
{
  int i = 0;
  int fieldNum = 0;
  char **fields = NULL;
  char delim[] = ",";
  CBed6 *bed = NULL;
  fields = splitWhitespace(line, &fieldNum);
  if (fieldNum < 6)
  {
    freeWords(fields, fieldNum);
    fprintf(stderr, "the format of annotation file must be be6\n");
    exit(1);
  }
  bed = (CBed6 *)safeMalloc(sizeof(CBed6));
  bed->chrom = strClone(fields[0]);
  bed->chromStart = atoi(fields[1]);
  bed->chromEnd = atoi(fields[2]);
  bed->name = strClone(fields[3]);
  bed->score = atof(fields[4]);
  bed->strand = fields[5][0];
  freeWords(fields, fieldNum);
  return bed;
}

// warn bed12 information
void warnBed12(CBed12 *bed12)
{
  int blockCount, i;
  blockCount = bed12->blockCount;
  fprintf(stderr, "%s\t%d\t%d\t%s\t%.0f\t%c\t%d\t%d\t0\t%d\t",
          bed12->chrom, bed12->chromStart, bed12->chromEnd, bed12->name, bed12->score, bed12->strand, bed12->thickStart, bed12->thickEnd, bed12->blockCount);
  for (i = 0; i < blockCount; i++)
  {
    fprintf(stderr, "%d,", bed12->blockSizes[i]);
  }
  fprintf(stderr, "\t");
  for (i = 0; i < blockCount; i++)
  {
    fprintf(stderr, "%d,", bed12->chromStarts[i]);
  }
  fprintf(stderr, "\n");
}

CBed12 *twoBed6ToBed12(CBed6 *bed1, CBed6 *bed2)
// copy oBed to tBed
{
  CBed12 *bed = (CBed12 *)safeMalloc(sizeof(CBed12));
  bed->chrom = strClone(bed1->chrom);
  bed->chromStart = MIN(bed1->chromStart, bed2->chromStart);
  bed->chromEnd = MAX(bed1->chromEnd, bed2->chromEnd);
  bed->name = strClone(bed1->name);
  bed->score = (bed1->score + bed2->score) / 2;
  bed->strand = bed1->strand;
  bed->thickStart = bed->chromStart;
  bed->thickEnd = bed->chromStart;
  bed->itemRgb = 0;
  if (overlapLength(bed1->chromStart, bed1->chromEnd, bed2->chromStart, bed2->chromEnd) > 0)
  {
    bed->blockCount = 1;
    bed->blockSizes = (int *)safeMalloc(sizeof(int) * bed->blockCount);
    bed->chromStarts = (int *)safeMalloc(sizeof(int) * bed->blockCount);
    bed->blockSizes[0] = bed->chromEnd - bed->chromStart;
    bed->chromStarts[0] = bed->chromStart - bed->chromStart;
  }
  else
  {
    bed->blockCount = 2;
    bed->blockSizes = (int *)safeMalloc(sizeof(int) * bed->blockCount);
    bed->chromStarts = (int *)safeMalloc(sizeof(int) * bed->blockCount);
    if (bed1->chromStart < bed2->chromStart)
    {
      bed->blockSizes[0]  = bed1->chromEnd - bed1->chromStart;
      bed->chromStarts[0] = bed1->chromStart - bed->chromStart;
      bed->blockSizes[1]  = bed2->chromEnd - bed2->chromStart;
      bed->chromStarts[1] = bed2->chromStart - bed->chromStart;
    }
    else
    {
      bed->blockSizes[0]  = bed2->chromEnd - bed2->chromStart;
      bed->chromStarts[0] = bed2->chromStart - bed->chromStart;
      bed->blockSizes[1]  = bed1->chromEnd - bed1->chromStart;
      bed->chromStarts[1] = bed1->chromStart - bed->chromStart;
    }
  }
  return bed;
}

CBed12 *mergeTwoBed12ToBed12(CBed12 *bed1, CBed12 *bed2)
// copy oBed to tBed
{
  int i = 0;
  int j = 0;
  int firstTag = 0;
  int blockCount = 1;
  CBed12 *bed = (CBed12 *)safeMalloc(sizeof(CBed12));
  bed->chrom = strClone(bed1->chrom);
  bed->chromStart = MIN(bed1->chromStart, bed2->chromStart);
  bed->chromEnd = MAX(bed1->chromEnd, bed2->chromEnd);
  bed->name = strClone(bed1->name);
  bed->score = (bed1->score + bed2->score) / 2;
  bed->strand = bed1->strand;
  bed->thickStart = bed->chromStart;
  bed->thickEnd = bed->chromStart;
  bed->itemRgb = 0;
  int geneLen = bed->chromEnd - bed->chromStart;
  int *geneMatrix = (int *)safeMalloc(sizeof(int) * geneLen);
  for (i = 0; i < geneLen; i++)
  {
    geneMatrix[i] = 0;
  }
  for (i = 0; i < bed1->blockCount; i++)
  {
    int start = bed1->chromStart + bed1->chromStarts[i] - bed->chromStart;
    int end = start + bed1->blockSizes[i];
    for (j = start; j < end; j++)
    {
      geneMatrix[j] = 1;
    }
  }
  for (i = 0; i < bed2->blockCount; i++)
  {
    int start = bed2->chromStart + bed2->chromStarts[i] - bed->chromStart;
    int end = start + bed2->blockSizes[i];
    for (j = start; j < end; j++)
    {
      geneMatrix[j] = 1;
    }
  }
  for (i = 0; i < geneLen; i++)
  {
    if (i < geneLen - 1)
    {
      if (geneMatrix[i] == 0 && geneMatrix[i + 1] >= 1) blockCount += 1;
    }
  }

  bed->blockCount = blockCount;
  bed->blockSizes = (int *)safeMalloc(sizeof(int) * blockCount);
  bed->chromStarts = (int *)safeMalloc(sizeof(int) * blockCount);

  int blockStart = 0;
  int blockEnd = 0;
  int blockSize = 0;
  int bCount = 0;
  firstTag = 1;
  for (i = 0; i < geneLen; i++)
  {
    if (geneMatrix[i] != 0)
    {
      blockEnd = i;
      if (firstTag)
      {
        blockStart = i;
        firstTag = 0;
      }
    }
    if ((i > 0 && geneMatrix[i] == 0 && geneMatrix[i - 1] > 0) || i == geneLen - 1)
    {
      blockSize = blockEnd - blockStart + 1;
      bed->chromStarts[bCount] = blockStart;
      bed->blockSizes[bCount] = blockSize;
      firstTag = 1;
      bCount++;
    }
  }

  safeFree(geneMatrix);
  if (bCount != blockCount)
  {
    fprintf(stderr, "the blockCount is wrong: blockCount:%d != bCount:%d\n", blockCount, bCount);
    exit(1);
  }

  return bed;
}

CBed6 *mergeTwoBed6(CBed6 *bed1, CBed6 *bed2)
// copy oBed to tBed
{
  CBed6 *bed = (CBed6 *)safeMalloc(sizeof(CBed6));
  bed->chrom = strClone(bed1->chrom);
  bed->chromStart = MIN(bed1->chromStart, bed2->chromStart);
  bed->chromEnd = MAX(bed1->chromEnd, bed2->chromEnd);
  bed->name = strClone(bed1->name);
  bed->score = (bed1->score + bed2->score) / 2;
  bed->strand = bed1->strand;
  return bed;
}


CBed6 *createBed6(const string &chrom, int chromStart, int chromEnd,  const string &name, double score, char strand)
{
  CBed6 *bed = NULL;
  bed = (CBed6 *)safeMalloc(sizeof(CBed6));
  bed->chrom = strClone(const_cast<char *>(chrom.c_str()));
  bed->chromStart = chromStart;
  bed->chromEnd = chromEnd;
  bed->name = strClone(const_cast<char*>(name.c_str()));
  bed->score = score;
  bed->strand = strand;
  return bed;
}

double normalizedReadsByNH(chromBed6Map &bed6Hash)
{
  double totalNum = 0;
  chromBed6Map::iterator it;
  for (it = bed6Hash.begin(); it != bed6Hash.end(); ++it)
  {
    bed6Vector bedList = it->second;
    if (bedList.size() >= 1)
    {
      for (bed6Vector::iterator vecItr = bedList.begin(); vecItr != bedList.end(); vecItr++)
      {
        CBed6 *bed6 = *vecItr;
        bed6->score = bed6->score / bed6->locusNum;
        totalNum += bed6->score;
      } // for end
    }
  } // for bed hash
  return totalNum;
}

double normalizedBed6Reads(chromBed6Map &bed6Hash)
{
  double totalNum = 0;
  chromBed6Map::iterator it;
  map<string, int> mapReads;
  for (it = bed6Hash.begin(); it != bed6Hash.end(); ++it)
  {
    bed6Vector bedList = it->second;
    if (bedList.size() >= 1)
    {
      for (bed6Vector::iterator vecItr = bedList.begin(); vecItr != bedList.end(); vecItr++)
      {
        CBed6 *bed6 = *vecItr;
        mapReads[bed6->name] = 0;
      } // for end
    }
  } // for bed hash
  for (it = bed6Hash.begin(); it != bed6Hash.end(); ++it)
  {
    bed6Vector bedList = it->second;
    if (bedList.size() >= 1)
    {
      for (bed6Vector::iterator vecItr = bedList.begin(); vecItr != bedList.end(); vecItr++)
      {
        CBed6 *bed6 = *vecItr;
        mapReads[bed6->name]++;
      } // for end
    }
  } // for bed hash
  for (it = bed6Hash.begin(); it != bed6Hash.end(); ++it)
  {
    bed6Vector bedList = it->second;
    if (bedList.size() >= 1)
    {
      for (bed6Vector::iterator vecItr = bedList.begin(); vecItr != bedList.end(); vecItr++)
      {
        CBed6 *bed6 = *vecItr;
        bed6->score = bed6->score / mapReads[bed6->name];
        totalNum += bed6->score;
      } // for end
    }
  } // for bed hash
  return totalNum;
}

double normalizedBed12Reads(chromBed12Map &bed12Hash)
{
  double totalNum = 0;
  chromBed12Map::iterator it;
  map<string, int> mapReads;
  for (it = bed12Hash.begin(); it != bed12Hash.end(); ++it)
  {
    bed12Vector bedList = it->second;
    if (bedList.size() >= 1)
    {
      for (bed12Vector::iterator vecItr = bedList.begin(); vecItr != bedList.end(); vecItr++)
      {
        CBed12 *bed12 = *vecItr;
        mapReads[bed12->name] = 0;
      } // for end
    }
  } // for bed hash
  for (it = bed12Hash.begin(); it != bed12Hash.end(); ++it)
  {
    bed12Vector bedList = it->second;
    if (bedList.size() >= 1)
    {
      for (bed12Vector::iterator vecItr = bedList.begin(); vecItr != bedList.end(); vecItr++)
      {
        CBed12 *bed12 = *vecItr;
        mapReads[bed12->name]++;
      } // for end
    }
  } // for bed hash
  for (it = bed12Hash.begin(); it != bed12Hash.end(); ++it)
  {
    bed12Vector bedList = it->second;
    if (bedList.size() >= 1)
    {
      for (bed12Vector::iterator vecItr = bedList.begin(); vecItr != bedList.end(); vecItr++)
      {
        CBed12 *bed12 = *vecItr;
        bed12->score = bed12->score / mapReads[bed12->name];
        totalNum += bed12->score;
      } // for end
    }
  } // for bed hash
  return totalNum;
}

double removeMisprimReadBed12(char *primerSeq, int barcodeLen,
                             FILE *genomefp,
                             faidxMap &faiHash,
                             chromBed12Map &bedHash,
                             chromBed12Map &newBedHash)
{
  fprintf(stderr, "primer:%s\tbarcodeLen=%d\n", primerSeq, barcodeLen);
  chromBed12Map::iterator it;
  int i = 0;
  int j = 0;
  int primerNum = 0;
  double totalNum = 0;
  int primerSeqLen = strlen(primerSeq);
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
        int start = bed->chromStart;
        int end = bed->chromEnd;
        int primerStart = end + barcodeLen;
        int primerEnd = primerStart + primerSeqLen;
        if (bed->strand == '-')
        {
          primerEnd = start - barcodeLen;
          primerStart = primerEnd - primerSeqLen;
        }
        if (primerStart < 0) primerStart = 0;
        if (primerEnd < 0) primerEnd = 0;
        int matchNum = 0;
        char *fetchSeq = faidxFetchSeq(genomefp, fai, primerStart, primerEnd, bed->strand);
        int fetchLen = strlen(fetchSeq);
        for (i = 0; i < primerSeqLen && i < fetchLen; i++)
        {
          if (fetchSeq[i] == primerSeq[i])
          {
            matchNum++;
          }
        }
        if (matchNum >= primerSeqLen - 1)
        {
          primerNum++;
        }
        else {
          newBedHash[chromStr].push_back(bed);
          totalNum += bed->score;
        }
        safeFree(fetchSeq);
      } // for end
    }// if sizes
  } // for bed hash
  fprintf(stderr, "#remove the mispriming reads=%d\n", primerNum);
  return totalNum;
}
