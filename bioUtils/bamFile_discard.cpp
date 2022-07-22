#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<ctype.h>
#include<assert.h>
#include<math.h>
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
#include "bedFile.h"
#include "faiFile.h"
#include "bamFile.h"

template<typename T>
string NumberToString(T Number)
{
  ostringstream ss;
  ss << Number;
  return ss.str();
}

string BuildCigarString(const vector<CigarOp> &cigar)
{

  stringstream cigarString;
  for (size_t i = 0; i < cigar.size(); ++i)
  {
    switch (cigar[i].Type)
    {
    case ('M') :
    case ('I') :
    case ('D') :
    case ('N') :
    case ('S') :
    case ('H') :
    case ('P') :
      cigarString << cigar[i].Length << cigar[i].Type;
    }
  }
  return cigarString.str();
}

string PrintTag(const BamAlignment &bam, const string &tag)
{
  uint32_t uTagValue;
  int32_t sTagValue;
  ostringstream value;
  if (bam.GetTag(tag, uTagValue))
    value << uTagValue;
  else if (bam.GetTag(tag, sTagValue))
    value << sTagValue;
  else
  {
    cerr << "The requested tag ("
         << tag
         << ") was not found in the BAM file.  Exiting\n";
    exit(1);
  }
  return value.str();
}

void openBamFile(char *bamFile, BamReader &reader)
{
  if (!reader.Open(bamFile))
  {
    cerr << "Failed to open BAM file " << bamFile << endl;
    exit(1);
  }
  else
  {
    reader.LocateIndex();
  }
  if ( !reader.HasIndex())
  {
    cerr << "Failed to open BAI file " << bamFile << endl;
    //exit(1);
  }
}

void openBamNoIdxFile(char *bamFile, BamReader &reader)
{
  if (!reader.Open(bamFile))
  {
    cerr << "Failed to open BAM file " << bamFile << endl;
    exit(1);
  }
}

int readBamToLocusNum(BamReader &reader, map<string, int> &readLocus)
{
  int i = 0;
  reader.Rewind();

  BamAlignment bam;
  while (reader.GetNextAlignment(bam))
  {
    if (bam.IsMapped() == true)
    {
      string readName = bam.Name;
      readLocus[readName] += 1;
      i++;
    }
  }
  return i;
}

long getGenomeSize(BamReader & reader, chromSizeMap &mapSize)
{
  // retrieve names and lengths from reference data
  long genomeSize = 0;
  RefVector m_references = reader.GetReferenceData();
  RefVector::const_iterator refIter = m_references.begin();
  RefVector::const_iterator refEnd  = m_references.end();
  for ( ; refIter != refEnd; ++refIter)
  {
    string refName = (*refIter).RefName;
    int refLen = (*refIter).RefLength;
    mapSize[refName] = refLen;
    genomeSize += refLen;
  }
  return genomeSize;
}


int readSamToLocusNum(FILE *fp, map<string, int> &readLocus)
{
  char *strLine;
  char *tmpLine;
  char delims[] = "-";
  int headTag = 0;
  int fieldNum = 0;
  char **fields   = NULL;
  char **infos = NULL;
  int infoNum = 0;
  int i = 0;
  int flag = 0;
  double totalNum = 1;
  if (feof(fp))
  {
    return 0;
  }
  //chr1  147646  147749  U6-related  1000  -
  while (strLine = getLine(fp))
  {
    tmpLine  = strLine;
    tmpLine  = skipStartWhitespace(tmpLine);
    fieldNum = 0;
    fields   = NULL;
    if (tmpLine[0] == '@')
    {
      safeFree(strLine);
      continue;
    }

    fields = splitWhitespace(tmpLine, &fieldNum);
    if (fieldNum >= 10 && strcmp(fields[2], "*") != 0)
    {
      string readName(fields[0]);
      if (readLocus.find(readName) == readLocus.end())
      {
        readLocus[readName] = 1;
      }
      else
      {
        readLocus[readName] += 1;
      }
      i++;
    }
    freeWords(fields, fieldNum);
    safeFree(strLine);
  }
  return i;
}

CBed6 *bamToBed6(BamAlignment &bam, string &chrom)
{
  CBed6 *bedPtr = NULL;
  string readName  = bam.Name;
  int  chromStart = bam.Position;
  int  chromEnd   = bam.GetEndPosition(false, false);
  bedPtr             = (CBed6 *)safeMalloc(sizeof(CBed6));
  bedPtr->chrom      = strClone(const_cast<char *>(chrom.c_str()));
  bedPtr->name       = strClone(const_cast<char*>(readName.c_str()));
  bedPtr->chromStart = chromStart;
  bedPtr->chromEnd   = chromEnd;
  bedPtr->strand     = '+';
  if (bam.IsReverseStrand() == true) bedPtr->strand = '-';
  if (bam.IsSecondMate() == true)
  {
    if (bedPtr->strand == '-')
    {
      bedPtr->strand = '+';
    }
    else
    {
      bedPtr->strand = '-';
    }
  }
  bedPtr->score = 1.0;
  bedPtr->next  = NULL;
  return bedPtr;
}

void bamToBlocks(const BamAlignment &bam, const string &chrom, bed6Vector &bedBlocks, bool skipDeletion, bool skipSplice)
/*modified from bedTools, Thanks!*/
{
  int chromStart = bam.Position;
  char strand = '+';
  if (bam.IsReverseStrand() == true) strand = '-';
  int blockLength  = 0;

  //  Rip through the CIGAR ops and figure out if there is more
  //  than one block for this alignment
  vector<CigarOp>::const_iterator cigItr = bam.CigarData.begin();
  vector<CigarOp>::const_iterator cigEnd = bam.CigarData.end();
  for (; cigItr != cigEnd; ++cigItr) {
    switch (cigItr->Type) {
    case ('M') :
      blockLength += cigItr->Length;
    case ('I') : break;
    case ('S') : break;
    case ('D') :
      if (skipDeletion)
        blockLength += cigItr->Length;
      else {
        bedBlocks.push_back( createBed6(chrom, chromStart, chromStart + blockLength,
                                        bam.Name, 1.0, strand) );
        chromStart += cigItr->Length + blockLength;
        blockLength = 0;
      }
    case ('P') : break;
    case ('N') :
      if (skipSplice)
        blockLength += cigItr->Length;
      else {
        bedBlocks.push_back( createBed6(chrom, chromStart, chromStart + blockLength,
                                        bam.Name, 1.0, strand) );
        chromStart += cigItr->Length + blockLength;
        blockLength = 0;
      }
    case ('H') : break;                             // for 'H' - do nothing, move to next op
    default    :
      printf("ERROR: Invalid Cigar op type\n");   // shouldn't get here
      exit(1);
    }
  }
  bedBlocks.push_back( createBed6(chrom, chromStart, chromStart + blockLength,
                                  bam.Name, 1.0, strand) );
}

CBed12 *bamToBed12(const BamAlignment &bam, const string &chrom, bool skipDeletion, bool skipSplice)
{
  int i;
  CBed12 *bedPtr = NULL;
  string readName  = bam.Name;
  if (bam.IsFirstMate()) readName += "/1";
  if (bam.IsSecondMate()) readName += "/2";
  int  chromStart = bam.Position;
  int  chromEnd   = bam.GetEndPosition(false, false);
  bedPtr             = (CBed12 *)safeMalloc(sizeof(CBed12));
  bedPtr->chrom      = strClone(const_cast<char *>(chrom.c_str()));
  bedPtr->name       = strClone(const_cast<char*>(readName.c_str()));
  bedPtr->chromStart = chromStart;
  bedPtr->chromEnd   = chromEnd;
  bedPtr->strand     = '+';
  if (bam.IsReverseStrand() == true) bedPtr->strand = '-';
  bedPtr->score = 1.0;
  bedPtr->next  = NULL;
  bedPtr->thickStart = bedPtr->chromStart;
  bedPtr->thickEnd = bedPtr->chromStart;
  bedPtr->itemRgb = 0;
  bed6Vector bedBlocks;
  bamToBlocks(bam, chrom, bedBlocks, skipDeletion, skipSplice);
  int blockCount = bedBlocks.size();
  bedPtr->blockCount = blockCount;
  bedPtr->blockSizes = (int *)safeMalloc(sizeof(int) * bedPtr->blockCount);
  bedPtr->chromStarts = (int *)safeMalloc(sizeof(int) * bedPtr->blockCount);
  for (i = 0; i < blockCount - 1; i++)
  {
    bedPtr->blockSizes[i] = bedBlocks[i]->chromEnd - bedBlocks[i]->chromStart;
    bedPtr->chromStarts[i] = bedBlocks[i]->chromStart - bedPtr->chromStart;
  }
  return bedPtr;
}

double readSEBamToBed12Map(BamReader &reader, chromBed12Map &bed12Hash, bool skipDeletion, bool skipSplice)
// must at same chromosome
{
  int i            = 0;
  double totalNum  = 0;
  int insertSize   = 0;
  // get header & reference information
  reader.Rewind();
  string header  = reader.GetHeaderText();
  RefVector refs = reader.GetReferenceData();
  // rip through the BAM file and convert each mapped entry to BED
  BamAlignment bam;
  while (reader.GetNextAlignment(bam))
  {
    if (bam.IsMapped() == true)
    {
      string chrom = refs.at(bam.RefID).RefName;
      CBed12 *bed12Ptr = bamToBed12(bam, chrom, skipDeletion, skipSplice);
      bed12Hash[chrom].push_back(bed12Ptr);
      totalNum++;
    }
  }
  return totalNum;
}

double readPEBamToBed6MapNew(BamReader &reader, chromBed6Map &bed6Hash, int maxInsertLen)
// must at same chromosome
{
  int i            = 0;
  double totalNum  = 0;
  int insertSize   = 0;
  // get header & reference information
  reader.Rewind();
  string header  = reader.GetHeaderText();
  RefVector refs = reader.GetReferenceData();
  // rip through the BAM file and convert each mapped entry to BED
  BamAlignment bam1;
  BamAlignment bam2;
  while (reader.GetNextAlignment(bam1))
  {
    if (bam1.IsMapped() && bam1.IsPaired() && bam1.RefID == bam1.MateRefID) // same chrom
    {
      // get paired bam items
      reader.GetNextAlignment(bam2);
      if (bam2.IsMapped())
      {
        if (bam1.Name != bam2.Name)
        {
          fprintf(stderr, "*****WARNING: Query is marked as paired, \n but its mate does not occur next to it in your BAM file. Skipping.\n");
          continue;
        }
        else
        {
          string chrom1 = refs.at(bam1.RefID).RefName;
          string chrom2 = refs.at(bam2.RefID).RefName;
          CBed6 *bed1 = bamToBed6(bam1, chrom1);
          CBed6 *bed2 = bamToBed6(bam2, chrom2);
          CBed6 *bed6Ptr = mergeTwoBed6(bed1, bed2);
          freeBed6Item(bed1);
          freeBed6Item(bed2);
          if (abs(bed6Ptr->chromEnd - bed6Ptr->chromStart) > maxInsertLen)
          {
            freeBed6Item(bed6Ptr);
            continue;
          }
          else {
            bed6Hash[chrom1].push_back(bed6Ptr);
            totalNum++;
          }
        }
      }
    }
  }
  return totalNum;
}

double readPEBamToBed12Map(BamReader &reader, chromBed12Map &bed12Hash, bool skipDeletion, bool skipSplice)
// must at same chromosome
{
  int i            = 0;
  double totalNum  = 0;
  int insertSize   = 0;
  // get header & reference information
  reader.Rewind();
  string header  = reader.GetHeaderText();
  RefVector refs = reader.GetReferenceData();
  // rip through the BAM file and convert each mapped entry to BED
  BamAlignment bam1;
  BamAlignment bam2;
  while (reader.GetNextAlignment(bam1))
  {
    if (bam1.IsMapped() && bam1.IsPaired() && bam1.RefID == bam1.MateRefID) // same chrom
    {
      // get paired bam items
      reader.GetNextAlignment(bam2);
      if (bam2.IsMapped())
      {
        if (bam1.Name != bam2.Name)
        {
          fprintf(stderr, "*****WARNING: Query is marked as paired, \n but its mate does not occur next to it in your BAM file. Skipping.\n");
          continue;
        }
        else
        {
          string chrom1 = refs.at(bam1.RefID).RefName;
          string chrom2 = refs.at(bam2.RefID).RefName;
          CBed12 *bed1 = bamToBed12(bam1, chrom1, skipDeletion, skipSplice);
          CBed12 *bed2 = bamToBed12(bam2, chrom2, skipDeletion, skipSplice);
          CBed12 *bed12Ptr = mergeTwoBed12ToBed12(bed1, bed2);
          bed12Hash[chrom1].push_back(bed12Ptr);
          freeBed12Item(bed1);
          freeBed12Item(bed2);
          totalNum++;
        }
      }
    }
  }
  return totalNum;
}

double readBamToSignals(BamReader &reader, chromSignalMap &sigHash, int keepDup, int normalization)
{
  int i           = 0;
  char delims[]   = "-";
  int fieldNum    = 0;
  char **infos    = NULL;
  int infoNum     = 0;
  double totalNum = 0;
  string NHtag    = "NH";
  int    nhTagVal = 1;
  // get header & reference information
  reader.Rewind();
  BamAlignment bam;
  string header = reader.GetHeaderText();
  RefVector refs = reader.GetReferenceData();
  // rip through the BAM file and convert each mapped entry to BED
  for (RefVector::iterator vecItr = refs.begin(); vecItr != refs.end(); vecItr++)
  {
    RefData ref = *vecItr;
    string chrom = ref.RefName;
    int chromSize = ref.RefLength;
    double *sigValPlus = (double *)safeMalloc(sizeof(double) * chromSize);
    memset(sigValPlus, 0, sizeof(double)*chromSize);
    char strand = '+';
    string chromStr = chrom + strand;
    sigHash[chromStr] = sigValPlus;
    strand = '-';
    chromStr = chrom + strand;
    double *sigValMinus = (double *)safeMalloc(sizeof(double) * chromSize);
    memset(sigValMinus, 0, sizeof(double)*chromSize);
    sigHash[chromStr] = sigValMinus;
  }

  while (reader.GetNextAlignment(bam))
  {
    if (bam.IsMapped() == true)
    {
      string readName    = bam.Name;
      string chrom  = refs.at(bam.RefID).RefName;
      char strand = '+';
      if (bam.IsReverseStrand() == true) strand = '-';
      int mateFlag = 0;
      if (bam.AlignmentFlag & 0x1)
      {
        if (bam.AlignmentFlag & 0x40)
        {
          mateFlag = 1;
          if (strand == '+')
          {
            strand = '-';
          }
          else {
            strand = '+';
          }
        }
        if (bam.AlignmentFlag & 0x80) mateFlag = 2;
      }
      int alignPos = bam.Position;
      if (bam. GetTag(NHtag, nhTagVal))
      {
        if (nhTagVal < 1) nhTagVal = 1;
      }
      else
      {
        fprintf(stderr, "Error: there is not NH tag\n");
      }
      double score    = 1 / (double)nhTagVal;
      infos = splitString(const_cast<char *>(readName.c_str()), delims, &infoNum);
      if (infoNum == 2 && keepDup)
      {
        if (isdigit(infos[1][0]))
          score = atof(infos[1]) / (double)nhTagVal;
      }
      freeWords(infos, infoNum); // free memory
      totalNum += score;
      string chromStr = chrom + strand;
      double *sigChr = sigHash[chromStr];
      // iterate over cigar operations
      vector<CigarOp>::const_iterator cigarIter = bam.CigarData.begin();
      vector<CigarOp>::const_iterator cigarEnd  = bam.CigarData.end();
      for ( ; cigarIter != cigarEnd; ++cigarIter) {
        const CigarOp& op = (*cigarIter);
        int opLen = op.Length;
        switch ( op.Type ) {
        // increase end position on CIGAR chars [DMXN=]
        case Constants::BAM_CIGAR_DEL_CHAR      :
        case Constants::BAM_CIGAR_REFSKIP_CHAR  :
          alignPos += opLen;
          break;
        case Constants::BAM_CIGAR_MATCH_CHAR    :
          for (i = 0; i < opLen; i++) {
            sigChr[alignPos + i] += score;
            alignPos++;
          }
        // all other CIGAR chars do not affect end position
        default :
          break;
        } // swith option
      } // for cigar
    } // if mapped
  } // while
  return totalNum;
}


double readBamToBedSignalMap(BamReader &reader, chromBedMap &bedHash, int keepDup, int maxLocusNum)
{
  int i           = 0;
  char delims[]   = "-";
  int fieldNum    = 0;
  char **infos    = NULL;
  int infoNum     = 0;
  double totalNum = 0;
  CBed *bedPtr    = NULL;
  map<string, int> dupMap;
  // get header & reference information
  reader.Rewind();
  string header = reader.GetHeaderText();
  RefVector refs = reader.GetReferenceData();
  // rip through the BAM file and convert each mapped entry to BED
  BamAlignment bam;
  while (reader.GetNextAlignment(bam))
  {
    if (bam.IsMapped() == true)
    {
      string nhStr;
      uint8_t nhVal = 1;
      if (!bam.GetTag("NH", nhVal)) {
        fprintf(stderr, "The requested tag NH is not exist in bam file\n");
        exit(1);
      }
      else {
        //char *nh = const_cast<char*>(nhStr.c_str());
        //fprintf(stderr, "NH value is %d\n", nhVal);
        //nhVal = atoi(nh);
      }
      if (nhVal > maxLocusNum) continue;
      string readName    = bam.Name;
      string chrom  = refs.at(bam.RefID).RefName;
      char strand = '+';
      if (bam.IsReverseStrand() == true) strand = '-';
      int mateFlag = 0;
      if (bam.AlignmentFlag & 0x1)
      {
        if (bam.AlignmentFlag & 0x40) mateFlag = 1;
        if (bam.AlignmentFlag & 0x80) mateFlag = 2;
      }
      int chromStart = bam.Position;
      int chromEnd = bam.GetEndPosition(false, false);
      // get score
      double score    = 1 / (double)nhVal;
      infos = splitString(const_cast<char *>(readName.c_str()), delims, &infoNum);
      if (infoNum == 2 && keepDup)
      {
        if (isdigit(infos[1][0]))
          score = atof(infos[1]) / (double)nhVal;
      }
      freeWords(infos, infoNum); // free memory

      vector<CigarOp>::const_iterator cigarIter = bam.CigarData.begin();
      vector<CigarOp>::const_iterator cigarEnd  = bam.CigarData.end();
      for ( ; cigarIter != cigarEnd; ++cigarIter) {
        const CigarOp& op = (*cigarIter);
        int opLen = op.Length;
        switch ( op.Type ) {
        // increase end position on CIGAR chars [DMXN=]
        case Constants::BAM_CIGAR_DEL_CHAR      :
        case Constants::BAM_CIGAR_REFSKIP_CHAR  :
          chromStart += opLen;
          break;
        case Constants::BAM_CIGAR_MATCH_CHAR    :
          chromEnd = chromStart + opLen;
          bedPtr             = (CBed *)safeMalloc(sizeof(CBed));
          bedPtr->chrom      = strClone(const_cast<char *>(chrom.c_str()));
          bedPtr->chromStart = chromStart;
          bedPtr->chromEnd   = chromEnd;
          bedPtr->strand     = strand;
          bedPtr->mateFlag   = mateFlag;
          bedPtr->score      = score;
          bedPtr->next       = NULL;
          bedHash[chrom].push_back(bedPtr);
          chromStart += opLen;
        } // swith option
      } // for cigar
      totalNum += score;
    }
  }
  return totalNum;
}

double readBamToBed6Map(BamReader &reader, chromBed6Map &bed6Hash, int keepDup)
{
  int i           = 0;
  char delims[]   = "-";
  int fieldNum    = 0;
  char **infos    = NULL;
  int infoNum     = 0;
  double totalNum = 0;
  CBed6 *bedPtr    = NULL;
  map<string, int> dupMap;
  // get header & reference information
  reader.Rewind();
  string header = reader.GetHeaderText();
  RefVector refs = reader.GetReferenceData();
  // rip through the BAM file and convert each mapped entry to BED
  BamAlignment bam;
  while (reader.GetNextAlignment(bam))
  {
    if (bam.IsMapped() == true)
    {
      string readName    = bam.Name;
      string chrom  = refs.at(bam.RefID).RefName;
      char strand = '+';
      if (bam.IsReverseStrand() == true) strand = '-';
      int mateFlag = 0;
      if (bam.AlignmentFlag & 0x1)
      {
        if (bam.AlignmentFlag & 0x40) mateFlag = 1;
        if (bam.AlignmentFlag & 0x80) mateFlag = 2;
      }
      uint8_t nhVal = 1;
      if (!bam.GetTag("NH", nhVal)) {
        //if (warnning)
        //  fprintf(stderr, "The requested tag NH is not exist in bam file\n");
        //exit(1);
        //warnning = 0;
        nhVal = 1;
      }
      bedPtr             = (CBed6 *)safeMalloc(sizeof(CBed6));
      bedPtr->chrom      = strClone(const_cast<char *>(chrom.c_str()));
      bedPtr->chromStart = bam.Position;
      bedPtr->chromEnd   = bam.GetEndPosition(false, false);
      bedPtr->strand     = strand;
      bedPtr->name       = strClone(const_cast<char*>(readName.c_str()));
      bedPtr->mateFlag   = mateFlag;
      bedPtr->locusNum   = nhVal;
      //if (bam.IsReverseStrand() == true) bedPtr->strand = '-';
      bedPtr->score      = 1;
      infos = splitString(const_cast<char *>(readName.c_str()), delims, &infoNum);
      if (infoNum == 2 && keepDup)
      {
        if (isdigit(infos[1][0]))
          bedPtr->score = atof(infos[1]);
      }
      totalNum += bedPtr->score;
      freeWords(infos, infoNum); // free memory
      bedPtr->next = NULL;
      bed6Hash[chrom].push_back(bedPtr);
      i++;
    }
  }
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

double readBamToBedMap(BamReader &reader, chromBedMap &bedHash, map<string, int> &readLocus, int keepDup, int maxLocusNum)
{
  int i           = 0;
  char delims[]   = "-";
  int fieldNum    = 0;
  char **infos    = NULL;
  int infoNum     = 0;
  double totalNum = 0;
  CBed *bedPtr    = NULL;
  //map<string, int> dupMap;
  // get header & reference information
  reader.Rewind();
  string header = reader.GetHeaderText();
  RefVector refs = reader.GetReferenceData();
  // rip through the BAM file and convert each mapped entry to BED
  BamAlignment bam;
  while (reader.GetNextAlignment(bam))
  {
    if (bam.IsMapped() == true)
    {
      string readName    = bam.Name;
      int lociNum = 1;
      if (readLocus.find(readName) != readLocus.end())
      {
        lociNum = readLocus[readName];
      }
      if (lociNum > maxLocusNum) continue; // remove multiple locus
      string chrom  = refs.at(bam.RefID).RefName;
      char strand = '+';
      if (bam.IsReverseStrand() == true) strand = '-';
      /*string dupPos = chrom + NumberToString(bam.Position) + NumberToString(bam.GetEndPosition(false, false)) + strand;
      if (!keepDup)
      {
        if (dupMap.find(dupPos) != dupMap.end())
        {
          dupMap[dupPos] += 1;
          continue;
        }
        else {
          dupMap[dupPos]  = 1;
        }
      }*/
      int mateFlag = 0;
      if (bam.AlignmentFlag & 0x1)
      {
        if (bam.AlignmentFlag & 0x40) mateFlag = 1;
        if (bam.AlignmentFlag & 0x80) mateFlag = 2;
      }
      bedPtr             = (CBed *)safeMalloc(sizeof(CBed));
      bedPtr->chrom      = strClone(const_cast<char *>(chrom.c_str()));
      bedPtr->chromStart = bam.Position;
      bedPtr->chromEnd   = bam.GetEndPosition(false, false);
      bedPtr->strand     = strand;
      bedPtr->mateFlag   = mateFlag;
      //if (bam.IsReverseStrand() == true) bedPtr->strand = '-';
      bedPtr->score    = 1 / (double)lociNum;
      infos = splitString(const_cast<char *>(readName.c_str()), delims, &infoNum);
      if (infoNum == 2 && keepDup)
      {
        if (isdigit(infos[1][0]))
          bedPtr->score = atof(infos[1]) / (double)lociNum;
      }
      totalNum += bedPtr->score;
      freeWords(infos, infoNum); // free memory
      bedPtr->next = NULL;
      bedHash[chrom].push_back(bedPtr);
      i++;
    }
  }
  return totalNum;
}

int filterLowQualities(int qualityBase, int minScore, int minPercent, const char *qualities)
{
  int i = 0;
  int count = 0;
  int seqLen = strlen(qualities);
  for (i = 0; i < seqLen; i++) {
    int val = qualities[i] - qualityBase - minScore;
    //fprintf(stderr, "%c %d %d\n", qualities[i], qualities[i], val);
    if (val > 0) count++;
  }
  if (count / (double)seqLen * 100 >= minPercent) return 1;
  return 0;
}

double readPEBamToBed6Map(BamReader &reader, chromBed6Map &bed6Hash, int maxInsertLen)
{
  int i           = 0;
  char delims[]   = "-";
  int fieldNum    = 0;
  char **infos    = NULL;
  int infoNum     = 0;
  double totalNum = 0;
  int matePosition = 0;
  int insertSize   = 0;
  string readSeqName;
  map<string, int> pairStartHash;
  map<string, int> pairEndHash;
  CBed6 *bedPtr     = NULL;
  // get header & reference information
  reader.Rewind();
  string header = reader.GetHeaderText();
  RefVector refs = reader.GetReferenceData();
  // rip through the BAM file and convert each mapped entry to BED
  BamAlignment bam;
  while (reader.GetNextAlignment(bam))
  {
    if (bam.IsMapped() == true)
    {
      string readName  = bam.Name;
      if (readSeqName.compare(readName) != 0)
      {
        pairStartHash.clear();
        pairEndHash.clear();
      }
      readSeqName     = readName;
      int  chromStart = bam.Position;
      int  chromEnd   = bam.GetEndPosition(false, false);
      if (bam.IsSecondMate() != true) // skip mate1, keep mate2
      {
        pairStartHash[readName] = chromStart;
        pairEndHash[readName]   = chromEnd;
        continue;
      }
      /*if (pairStartHash.find(readName) == pairStartHash.end())
      {
        pairStartHash.clear();
        pairEndHash.clear();
        continue;
      }*/
      insertSize = bam.InsertSize;
      if (abs(insertSize) > maxInsertLen)
      {
        pairStartHash.clear();
        pairEndHash.clear();
        continue;
      }
      matePosition = bam.MatePosition;
      if (insertSize != 0 && bam.RefID == bam.MateRefID)
      {
        if (insertSize > 0)
        {
          chromEnd = chromStart + insertSize;
        }
        else
        {
          chromStart = matePosition;
        }
      }// if insert size
      else {
        pairStartHash.clear();
        pairEndHash.clear();
        continue;
      }
      bedPtr             = (CBed6 *)safeMalloc(sizeof(CBed6));
      string chrom       = refs.at(bam.RefID).RefName;
      bedPtr->chrom      = strClone(const_cast<char *>(chrom.c_str()));
      bedPtr->name       = strClone(const_cast<char*>(readName.c_str()));
      bedPtr->chromStart = chromStart;
      bedPtr->chromEnd   = chromEnd;
      bedPtr->strand     = '+';
      if (bam.IsReverseStrand() == true) bedPtr->strand = '-';
      if (bam.IsSecondMate() == true)
      {
        if (bedPtr->strand == '-')
        {
          bedPtr->strand = '+';
        }
        else
        {
          bedPtr->strand = '-';
        }
      }
      bedPtr->score = 1.0;
      totalNum += bedPtr->score;
      bedPtr->next = NULL;
      bed6Hash[chrom].push_back(bedPtr);
      pairStartHash.clear();
      pairEndHash.clear();
      i++;
    }
  }
  return totalNum;
}


double readPEBamSortedNameToBedMap(BamReader &reader, chromBedMap &bedHash, map<string, int> &readLocus,
                                   int keepDup, int maxLocusNum, int maxInsertLen, int readLen)
{
  int i           = 0;
  char delims[]   = "-";
  int fieldNum    = 0;
  char **infos    = NULL;
  int infoNum     = 0;
  double totalNum = 0;
  int matePosition = 0;
  int insertSize   = 0;
  int read1Len     = 0;
  int read2Len     = 0;
  int maxDiffLen   = 5;
  string readSeqName;
  map<string, int> pairHash;
  CBed *bedPtr     = NULL;
  // get header & reference information
  reader.Rewind();
  string header = reader.GetHeaderText();
  RefVector refs = reader.GetReferenceData();
  // rip through the BAM file and convert each mapped entry to BED
  BamAlignment bam;
  while (reader.GetNextAlignment(bam))
  {
    if (bam.IsMapped() == true)
    {
      string readName    = bam.Name;
      if (readSeqName.compare(readName) != 0)
      {
        pairHash.clear();
      }
      readSeqName = readName;
      int lociNum = 1;
      if (readLocus.find(readName) != readLocus.end())
      {
        lociNum = readLocus[readName];
      }
      if (lociNum > maxLocusNum * 2) continue; // remove multiple locus
      bedPtr             = (CBed *)safeMalloc(sizeof(CBed));
      string chrom       = refs.at(bam.RefID).RefName;
      bedPtr->chrom      = strClone(const_cast<char *>(chrom.c_str()));
      bedPtr->chromStart = bam.Position;
      bedPtr->chromEnd   = bam.GetEndPosition(false, false);
      read1Len           = bedPtr->chromEnd - bedPtr->chromStart;
      bedPtr->strand     = '+';
      if (bam.IsReverseStrand() == true) bedPtr->strand = '-';
      if (bam.IsSecondMate() == true)
      {
        if (bedPtr->strand == '-')
        {
          bedPtr->strand = '+';
        }
        else
        {
          bedPtr->strand = '-';
        }
      }
      matePosition = bam.MatePosition;
      insertSize = bam.InsertSize;
      if (insertSize != 0 && bam.RefID == bam.MateRefID && abs(insertSize) < maxInsertLen)
      {
        if (read1Len < readLen && abs(bedPtr->chromStart - matePosition) > maxDiffLen) {
          freeBedItem(bedPtr);
          continue;
        }
        string pairString = readName + NumberToString(bedPtr->chromStart);
        if (pairHash.find(pairString) == pairHash.end())
        {
          string readString = readName + NumberToString(matePosition);
          pairHash[readString] = read1Len;
          freeBedItem(bedPtr);
          continue;
        }
        else {
          read2Len = pairHash[pairString];
        }
        if (insertSize > 0)
        {
          bedPtr->chromEnd = bedPtr->chromStart + insertSize;
        }
        else
        {
          bedPtr->chromStart = matePosition;
        }
      }// if insert size
      else
      {
        freeBedItem(bedPtr);
        continue;
      }
      if ((bedPtr->chromEnd - bedPtr->chromStart) > readLen
          && (read1Len < readLen - maxDiffLen || read2Len < readLen - maxDiffLen)) // discard the degradome reads
      {
        freeBedItem(bedPtr);
        continue;
      }
      bedPtr->score = 1.0 / (double)lociNum;
      if (lociNum > 1) bedPtr->score *= 2;
      totalNum += bedPtr->score;
      bedPtr->next = NULL;
      bedHash[chrom].push_back(bedPtr);
      i++;
    }
  }
  return totalNum;
}

double readPEBamToBedMap(BamReader &reader, chromBedMap &bedHash, map<string, int> &readLocus,
                         int keepDup, int maxLocusNum, int maxInsertLen, int readLen)
{
  int i           = 0;
  char delims[]   = "-";
  int fieldNum    = 0;
  char **infos    = NULL;
  int infoNum     = 0;
  double totalNum = 0;
  int matePosition = 0;
  int insertSize   = 0;
  int read1Len     = 0;
  int read2Len     = 0;
  int maxDiffLen   = 5;
  map<string, int> pairHash;
  CBed *bedPtr     = NULL;
  // get header & reference information
  reader.Rewind();
  string header = reader.GetHeaderText();
  RefVector refs = reader.GetReferenceData();
  // rip through the BAM file and convert each mapped entry to BED
  BamAlignment bam;
  while (reader.GetNextAlignment(bam))
  {
    if (bam.IsMapped() == true)
    {
      string readName    = bam.Name;
      int lociNum = 1;
      if (readLocus.find(readName) != readLocus.end())
      {
        lociNum = readLocus[readName];
      }
      if (lociNum > maxLocusNum * 2) continue; // remove multiple locus
      /*const char* qualities = bam.Qualities.c_str();
      if (filterLowQualities(33, 20, 90, qualities) == 0){
        continue;
      }*/
      bedPtr             = (CBed *)safeMalloc(sizeof(CBed));
      string chrom       = refs.at(bam.RefID).RefName;
      bedPtr->chrom      = strClone(const_cast<char *>(chrom.c_str()));
      bedPtr->chromStart = bam.Position;
      bedPtr->chromEnd   = bam.GetEndPosition(false, false);
      read1Len           = bedPtr->chromEnd - bedPtr->chromStart;
      bedPtr->strand     = '+';
      if (bam.IsReverseStrand() == true) bedPtr->strand = '-';
      if (bam.IsSecondMate() == true)
      {
        if (bedPtr->strand == '-')
        {
          bedPtr->strand = '+';
        }
        else
        {
          bedPtr->strand = '-';
        }
      }
      matePosition = bam.MatePosition;
      insertSize = bam.InsertSize;
      if (insertSize != 0 && bam.RefID == bam.MateRefID && abs(insertSize) < maxInsertLen)
      {
        if (read1Len < readLen && abs(bedPtr->chromStart - matePosition) > maxDiffLen) {
          freeBedItem(bedPtr);
          continue;
        }
        string pairString = readName + NumberToString(bedPtr->chromStart);
        if (pairHash.find(pairString) == pairHash.end())
        {
          string readString = readName + NumberToString(matePosition);
          pairHash[readString] = read1Len;
          freeBedItem(bedPtr);
          continue;
        }
        else {
          read2Len = pairHash[pairString];
        }
        if (insertSize > 0)
        {
          bedPtr->chromEnd = bedPtr->chromStart + insertSize;
        }
        else
        {
          bedPtr->chromStart = matePosition;
        }
      }// if insert size
      else
      {
        freeBedItem(bedPtr);
        continue;
      }
      if ((bedPtr->chromEnd - bedPtr->chromStart) > readLen
          && (read1Len < readLen - maxDiffLen || read2Len < readLen - maxDiffLen)) // discard the degradome reads
      {
        freeBedItem(bedPtr);
        continue;
      }
      bedPtr->score = 1.0 / (double)lociNum;
      if (lociNum > 1) bedPtr->score *= 2;
      totalNum += bedPtr->score;
      bedPtr->next = NULL;
      bedHash[chrom].push_back(bedPtr);
      i++;
    }
  }
  return totalNum;
}

double readBamToSamMapNew(BamReader &reader, chromSamMap &samHash, int keepDup, int maxLocusNum)
{
  int i           = 0;
  char delims[]   = "-";
  int fieldNum    = 0;
  char **infos    = NULL;
  int infoNum     = 0;
  double totalNum = 0;
  int warnning    = 1;
  CSam *samPtr    = NULL;
  // get header & reference information
  reader.Rewind();
  string header = reader.GetHeaderText();
  RefVector refs = reader.GetReferenceData();
  // rip through the BAM file and convert each mapped entry to BED
  BamAlignment bam;
  while (reader.GetNextAlignment(bam))
  {
    if (bam.IsMapped() == true)
    {
      uint8_t nhVal = 1;
      if (!bam.GetTag("NH", nhVal)) {
        if (warnning)
          fprintf(stderr, "The requested tag NH is not exist in bam file\n");
        //exit(1);
        warnning = 0;
        nhVal = 1;
      }
      if (nhVal > maxLocusNum) continue;
      string readName    = bam.Name;
      samPtr             = (CSam *)safeMalloc(sizeof(CSam));
      string chrom       = refs.at(bam.RefID).RefName;
      samPtr->chrom      = strClone(const_cast<char *>(chrom.c_str()));
      samPtr->chromStart = bam.Position;
      samPtr->chromEnd   = bam.GetEndPosition(false, false);
      samPtr->readName   = strClone(const_cast<char*>(readName.c_str()));
      samPtr->strand     = '+';
      if (bam.IsReverseStrand() == true) samPtr->strand = '-';
      string cigar    = BuildCigarString(bam.CigarData);
      samPtr->cigar   = strClone(const_cast<char*>(cigar.c_str()));
      string readSeq  = bam.QueryBases;
      samPtr->readSeq = strClone(const_cast<char*>(readSeq.c_str()));
      string mdz;
      if (!bam.GetTag("MD", mdz)) {
        fprintf(stderr, "The requested tag MD is not exist in bam file\n");
        exit(1);
      }
      samPtr->mdz = strClone(const_cast<char*>(mdz.c_str()));
      samPtr->readNum    = 1 / (double)nhVal;
      infos = splitString(const_cast<char *>(readName.c_str()), delims, &infoNum);
      if (infoNum == 2 && keepDup)
      {
        if (isdigit(infos[1][0]))
          samPtr->readNum = atof(infos[1]) / (double)nhVal;
      }
      totalNum += samPtr->readNum;
      freeWords(infos, infoNum); // free memory
      samPtr->next = NULL;
      samHash[chrom].push_back(samPtr);
      i++;
    }// if mapped
  }// while read
  return totalNum;
}

double normalizedSamReads(chromSamMap &samHash)
{
  double totalNum = 0;
  chromSamMap::iterator it;
  map<string, int> mapReads;
  for (it = samHash.begin(); it != samHash.end(); ++it)
  {
    samVector samList = it->second;
    if (samList.size() >= 1)
    {
      for (samVector::iterator vecItr = samList.begin(); vecItr != samList.end(); vecItr++)
      {
        CSam *sam = *vecItr;
        mapReads[sam->readName] = 0;
      } // for end
    }
  } // for bed hash
  for (it = samHash.begin(); it != samHash.end(); ++it)
  {
    samVector samList = it->second;
    if (samList.size() >= 1)
    {
      for (samVector::iterator vecItr = samList.begin(); vecItr != samList.end(); vecItr++)
      {
        CSam *sam = *vecItr;
        mapReads[sam->readName]++;
      } // for end
    }
  } // for bed hash
  for (it = samHash.begin(); it != samHash.end(); ++it)
  {
    samVector samList = it->second;
    if (samList.size() >= 1)
    {
      for (samVector::iterator vecItr = samList.begin(); vecItr != samList.end(); vecItr++)
      {
        CSam *sam = *vecItr;
        sam->readNum = sam->readNum / mapReads[sam->readName];
        totalNum += sam->readNum;
      } // for end
    }
  } // for bed hash
  return totalNum;
}

double readBamToSamMap(BamReader &reader, chromSamMap &samHash, map<string, int> &readLocus, int keepDup, int maxLocusNum)
{
  int i           = 0;
  char delims[]   = "-";
  int fieldNum    = 0;
  char **infos    = NULL;
  int infoNum     = 0;
  double totalNum = 0;
  CSam *samPtr    = NULL;
  // get header & reference information
  reader.Rewind();
  string header = reader.GetHeaderText();
  RefVector refs = reader.GetReferenceData();
  // rip through the BAM file and convert each mapped entry to BED
  BamAlignment bam;
  while (reader.GetNextAlignment(bam))
  {
    if (bam.IsMapped() == true)
    {
      string readName    = bam.Name;
      int lociNum        = 1;
      if (readLocus.find(readName) != readLocus.end())
      {
        lociNum = readLocus[readName];
      }
      if (lociNum > maxLocusNum) continue; // remove multiple locus
      samPtr             = (CSam *)safeMalloc(sizeof(CSam));
      string chrom       = refs.at(bam.RefID).RefName;
      samPtr->chrom      = strClone(const_cast<char *>(chrom.c_str()));
      samPtr->chromStart = bam.Position;
      samPtr->chromEnd   = bam.GetEndPosition(false, false);
      samPtr->readName   = strClone(const_cast<char*>(readName.c_str()));
      samPtr->strand     = '+';
      if (bam.IsReverseStrand() == true) samPtr->strand = '-';
      string cigar    = BuildCigarString(bam.CigarData);
      samPtr->cigar   = strClone(const_cast<char*>(cigar.c_str()));
      string readSeq  = bam.QueryBases;
      samPtr->readSeq = strClone(const_cast<char*>(readSeq.c_str()));
      string mdz;
      if (!bam.GetTag("MD", mdz)) {
        fprintf(stderr, "The requested tag MD is not exist in bam file\n");
        exit(1);
      }
      samPtr->mdz = strClone(const_cast<char*>(mdz.c_str()));
      samPtr->readNum    = 1 / (double)lociNum;
      infos = splitString(const_cast<char *>(readName.c_str()), delims, &infoNum);
      if (infoNum == 2 && keepDup)
      {
        if (isdigit(infos[1][0]))
          samPtr->readNum = atof(infos[1]) / (double)lociNum;
      }
      totalNum += samPtr->readNum;
      freeWords(infos, infoNum); // free memory
      samPtr->next = NULL;
      samHash[chrom].push_back(samPtr);
      i++;
    }
  }
  return totalNum;
}

int getRegionCounts(BamReader &reader, char *chrom, int chromStart, int chromEnd, char strand)
{
  int regionReadNum = 0;
  reader.Rewind();
  int id = reader.GetReferenceID(chrom);
  BamRegion region(id, chromStart, id, chromEnd);
  // rip through the BAM file
  if ( (id != -1) && (reader.SetRegion(region)) )
  {
    BamAlignment bam;
    while (reader.GetNextAlignment(bam))
    {
      if (bam.IsMapped() == true)
      {
        char bamStrand   = '+';
        if (bam.IsReverseStrand() == true) bamStrand = '-';
        if (strand != bamStrand) continue;
        int regStart = bam.Position;
        int regEnd = bam.GetEndPosition(false, false);
        if (overlapLength(chromStart, chromEnd, regStart, regEnd) > 0)
        {
          regionReadNum += 1;
        }
      } // is mapped
    } //while bam
  } // if id AND
  return regionReadNum;
}


int getMDZ(char **fields, int num)
{
  int i = 0;
  int mdzIdx = 0;
  const char *mdz = "MD:Z:";
  for (i = 0; i < num; i++)
  {
    if (strncmp(fields[i], mdz, 5) == 0)
    {
      mdzIdx = i;
      break;
    }
  }
  return mdzIdx;
}

int getEndPos(char *cigar)
{
  const char *samcigar = "MIDNSHP=X";
  int start = 0;
  int end   = 0;
  int readLen   = 0;
  int i = 0;
  char c;
  for (i = 0; i < strlen(cigar); i++)
  {
    c = cigar[i];
    if (!isdigit(c))
    {
      cigar[i] = '\0';
      int len = atoi(cigar + start);
      switch ( c )
      {
      // increase end position on CIGAR chars [DMXN=]
      case 'D' :
      case 'M' :
      case 'X' :
      case 'N' :
      case '=' :
        readLen += len;
        break;
      }
      cigar[i] = c;
      start = i + 1;
    }
  }
  return readLen;
}

void freeBedVector(bedVector &bedList)
{
  for (bedVector::iterator vecItr = bedList.begin(); vecItr != bedList.end(); vecItr++)
  {
    CBed *bed = *vecItr;
    freeBedItem(bed);
  }
  bedList.clear();
}

void freeSamVector(samVector &samList)
{
  for (samVector::iterator vecItr = samList.begin(); vecItr != samList.end(); vecItr++)
  {
    CSam *sam = *vecItr;
    freeSamItem(sam);
  }
  samList.clear();
}

void freeChromBed12Map(chromBed12Map &bed12Hash) /*free bed map */
{
  for (chromBed12Map::iterator mapItr = bed12Hash.begin(); mapItr != bed12Hash.end(); mapItr++) {
    bed12Vector bed12List = mapItr->second;
    freeBed12Vector(bed12List);
  }
  bed12Hash.clear();
}

void freeChromBed6Map(chromBed6Map &bed6Hash) /*free bed map */
{
  for (chromBed6Map::iterator mapItr = bed6Hash.begin(); mapItr != bed6Hash.end(); mapItr++) {
    bed6Vector bed6List = mapItr->second;
    freeBed6Vector(bed6List);
  }
  bed6Hash.clear();
}

void freeChromBedMap(chromBedMap &bedHash) /*free bed map */
{
  for (chromBedMap::iterator mapItr = bedHash.begin(); mapItr != bedHash.end(); mapItr++) {
    bedVector bedList = mapItr->second;
    freeBedVector(bedList);
  }
  bedHash.clear();
}

void freeChromSignalMap(chromSignalMap &sigHash)
{
  for (chromSignalMap::iterator mapItr = sigHash.begin(); mapItr != sigHash.end(); mapItr++) {
    double *sigArray = mapItr->second;
    safeFree(sigArray);
  }
  sigHash.clear();
}

void freeChromSamMap(chromSamMap &samHash) /*free sam map */
{
  for (chromSamMap::iterator mapItr = samHash.begin(); mapItr != samHash.end(); mapItr++) {
    samVector samList = mapItr->second;
    freeSamVector(samList);
  }
  samHash.clear();
}

void freeSamItem(CSam *sam)
{
  safeFree(sam->chrom);
  safeFree(sam->readName);
  safeFree(sam->readSeq);
  safeFree(sam->mdz);
  safeFree(sam->cigar);
  safeFree(sam);
}

void copySam(CSam *tSam, CSam *oSam)
// copy oBed to tBed
{
  tSam->chrom      = strClone(oSam->chrom);
  tSam->chromStart = oSam->chromStart;
  tSam->chromEnd   = oSam->chromEnd;
  tSam->readLen    = oSam->readLen;
  tSam->readLen    = oSam->readLen;
  tSam->readNum    = oSam->readNum;
  tSam->readSeq    = strClone(oSam->readSeq);
  tSam->readName   = strClone(oSam->readName);
  tSam->mdz        = strClone(oSam->mdz);
  tSam->cigar      = strClone(oSam->cigar);
}

bool compareSam(const CSam *oSam, const CSam *tSam)
{
  return (oSam->chromStart < tSam->chromStart);
}

int compareSamReadNum(const void *a, const void *b)
{
  CSam *oSam = *(CSam **)a;
  CSam *tSam = *(CSam **)b;
  if ((tSam->readNum - oSam->readNum) > 0)
    return 1;
  else if ((tSam->readNum - oSam->readNum) < 0)
    return -1;
  else
    return 0;
}
