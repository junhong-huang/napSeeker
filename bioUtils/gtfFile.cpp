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
#include "gtfFile.h"

gtfLine *parseGtfLine(char *line)
{
  int i = 0;
  int fieldNum = 0;
  char **fields = NULL;
  int tmpFieldNum = 0;
  char **tmpFields = NULL;
  char delim[] = ",";
  char chopper = '\t';
  gtfLine *gtf = NULL;
  fields = splitByChar(line, chopper, &fieldNum);
  if (fieldNum != 9)
  {
    freeWords(fields, fieldNum);
    fprintf(stderr, "the format of annotation file must be gtf/gff\n");
    exit(1);
  }
  gtf = (gtfLine *)safeMalloc(sizeof(gtfLine));
  gtf->seqName = strClone(fields[0]);
  gtf->source = strClone(fields[1]);
  gtf->feature = strClone(fields[2]);
  gtf->start = atoi(fields[3]);
  gtf->end = atoi(fields[4]);
  gtf->score = atof(fields[5]);
  gtf->strand = fields[6][0];
  gtf->frame = atoi(fields[7]);

  gtf->geneId = NULL;
  gtf->transcriptId = NULL;
  gtf->geneName = NULL;
  gtf->transcriptName = NULL;
  gtf->exonId = NULL;
  gtf->exonNumber = 0;

  tmpFields = splitByChar(fields[8], ';', &tmpFieldNum);

  for (i = 0; i < tmpFieldNum; i++)
  {
    char *attributes = skipStartWhitespace(tmpFields[i]);
    if (startStr(attributes, "gene_id"))
    {
      gtf->geneId = getQuotedString(attributes, '"');
    }
    else if (startStr(attributes, "transcript_id"))
    {
      gtf->transcriptId = getQuotedString(attributes, '"');
    }
    else if (startStr(attributes, "gene_name"))
    {
      gtf->geneName = getQuotedString(attributes, '"');
    }
    else if (startStr(attributes, "transcript_name"))
    {
      gtf->transcriptName = getQuotedString(attributes, '"');
    }
    else if (startStr(attributes, "exon_id"))
    {
      gtf->exonId = getQuotedString(attributes, '"');
    }
    else if (startStr(attributes, "exon_number"))
    {
      int exonFieldNum = 0;
      char **exonField = splitWhitespace(attributes, &exonFieldNum);
      if (exonFieldNum != 2)
      {
        freeWords(exonField, exonFieldNum);
        fprintf(stderr, "the format of exon_number is ERROR%s\n", attributes);
        exit(1);
      }
      gtf->exonNumber = atoi(exonField[1]);
    }
  }
  freeWords(tmpFields, tmpFieldNum);
  freeWords(fields, fieldNum);
  return gtf;
}

int getGtfToBed12Map(FILE *gtffp, chromStrandBed12Map &bed12Hash, const char *type)
{
  int bedNum = 0;
  char *line = NULL;
  char *preName = NULL;
  gtfVector gtfList;
  while (line = getLine(gtffp))
  {
    if (feof(gtffp) || line == NULL)
    {
      safeFree(line);
      break;
    }
    if (line[0] == '#')
    {
      safeFree(line);
      continue;
    }
    gtfLine *gtf = parseGtfLine(line); // get bed information
    if (sameString(gtf->feature, "CDS") || sameString(gtf->feature, "exon")
        || sameString(gtf->feature, "stop_codon") || sameString(gtf->feature, "start_codon") )
    {
      if (sameString(type, "transcript") && gtf->transcriptId != NULL)
      {
        if (preName != NULL && sameString(gtf->transcriptId, preName))
        {
          gtfList.push_back(gtf);
        }
        else
        {
          if (gtfList.size() >= 1)
          {
            CBed12 *bed12  = mergeGtfToBed12(gtfList, type);
            if (bed12 != NULL)
            {
              //warnBed12(bed12);
              string chromStrand = bed12->chrom;
              chromStrand = chromStrand + bed12->strand;
              bed12Hash[chromStrand].push_back(bed12);
            }
            freeGtfVector(gtfList);
            bedNum++;
          }
          preName = gtf->transcriptId;
          gtfList.push_back(gtf);
        }
      }
      if (sameString(type, "gene") && gtf->geneId != NULL)
      {
        if (preName != NULL && sameString(gtf->geneId, preName))
        {
          gtfList.push_back(gtf);
        }
        else
        {
          if (gtfList.size() >= 1)
          {
            CBed12 *bed12  = mergeGtfToBed12(gtfList, type);
            if (bed12 != NULL)
            {
              string chromStrand = bed12->chrom;
              chromStrand = chromStrand + bed12->strand;
              bed12Hash[chromStrand].push_back(bed12);
            }
            freeGtfVector(gtfList);
            bedNum++;
          }
          preName = gtf->geneId;
          gtfList.push_back(gtf);
        }
      }
    }
    else
    {
      freeGtf(gtf);
    }
    safeFree(line); // free memory
  }
  // last time
  if (gtfList.size() >= 1)
  {
    CBed12 *bed12  = mergeGtfToBed12(gtfList, type);
    if (bed12 != NULL)
    {
      //warnBed12(bed12);
      string chromStrand = bed12->chrom;
      chromStrand = chromStrand + bed12->strand;
      bed12Hash[chromStrand].push_back(bed12);
    }
    freeGtfVector(gtfList);
    bedNum++;
  }
  return bedNum;
}

CBed12 *mergeGtfToBed12(gtfVector &gtfList, const char *type)
{
  int i = 0;
  int minStart = 0;
  int maxEnd = 0;
  int firstTag = 1;
  char *chrom = gtfList[0]->seqName;
  char *tid = gtfList[0]->transcriptId;
  char *gid = gtfList[0]->geneId;
  char strand = gtfList[0]->strand;
  CBed12 *bed = NULL;
  // get min-start and max-end
  for (gtfVector::iterator vecItr = gtfList.begin(); vecItr != gtfList.end(); vecItr++)
  {
    gtfLine *gtf = *vecItr;
    if (firstTag)
    {
      minStart = gtf->start;
      maxEnd   = gtf->end;
      firstTag = 0;
      continue;
    }
    minStart = MIN(minStart, gtf->start);
    maxEnd   = MAX(maxEnd, gtf->end);
  }

  int geneLen = maxEnd - minStart + 1;
  if (geneLen <= 1) return bed;

  int *geneMatrix = (int *)safeMalloc(sizeof(int) * geneLen);
  for (i = 0; i < geneLen; i++) geneMatrix[i] = 0;
  // put position to matrix
  for (gtfVector::iterator vecItr = gtfList.begin(); vecItr != gtfList.end(); vecItr++)
  {
    gtfLine *gtf = *vecItr;
    int start = gtf->start - minStart;
    int end = gtf->end - minStart + 1;
    int val = 1;
    if (sameString(gtf->feature, "CDS") || sameString(gtf->feature, "stop_codon") || sameString(gtf->feature, "start_codon")) val = 2;
    for (i = start; i < end; i++)
    {
      geneMatrix[i] = MAX(geneMatrix[i], val);
    }
  }

  bed = (CBed12 *)safeMalloc(sizeof(CBed12));
  bed->chrom = strClone(chrom);
  bed->chromStart = minStart - 1;
  bed->chromEnd = maxEnd;
  if (sameString(type, "gene")) bed->name = strClone(gid);
  if (sameString(type, "transcript")) bed->name = strClone(tid);
  bed->score = 900;
  bed->strand = strand;
  bed->itemRgb = 0;
  int blockCount = 1;
  int thickStart = 0;
  int thickEnd = 0;
  firstTag = 1;
  for (i = 0; i < geneLen; i++)
  {
    if (i < geneLen - 1)
    {
      if (geneMatrix[i] == 0 && geneMatrix[i + 1] >= 1) blockCount += 1;
    }
    if (geneMatrix[i] == 2)
    {
      thickEnd = i+1;
      if (firstTag)
      {
        firstTag = 0;
        thickStart = i;
      }
    }
  }
  bed->thickStart = bed->chromStart + thickStart;
  bed->thickEnd = bed->chromStart + thickEnd;
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

  if (bCount != blockCount)
  {
    fprintf(stderr, "the blockCount is wrong: %d\n", blockCount);
    exit(1);
  }
  return bed;
}

void freeGtfVector(gtfVector &gtfList)
{
  for (gtfVector::iterator vecItr = gtfList.begin(); vecItr != gtfList.end(); vecItr++)
  {
    gtfLine *gtf = *vecItr;
    freeGtf(gtf);
  }
  gtfList.clear();
}

void freeGtf(gtfLine *gtf)
{
  safeFree(gtf->seqName);
  safeFree(gtf->source);
  safeFree(gtf->feature);
  if (gtf->geneId != NULL) safeFree(gtf->geneId);
  if (gtf->transcriptId != NULL) safeFree(gtf->transcriptId);
  if (gtf->exonId != NULL) safeFree(gtf->exonId);
  if (gtf->transcriptName != NULL) safeFree(gtf->transcriptName);
  if (gtf->geneName != NULL) safeFree(gtf->geneName);
  safeFree(gtf);
}
