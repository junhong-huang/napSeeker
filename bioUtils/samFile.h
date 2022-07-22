#ifndef BAM_HEAD_H
#define BAM_HEAD_H

struct CSamInfo
{
	struct CSamInfo *next; // for list
	char *chrom;	// chromosome name
	int chromStart; // chromosome start
	int chromEnd; // chromosome end
	int readLen; // read length
	int misMatches; // mismatches
	int mateFlag; // mate flag
	char *readName; // read name
	char *readSeq; // read sequence
	char *mdz; // aligned mdz
	char *cigar; // aligned cigar
	double readNum; // normalized read number
	char strand; // strand
};

typedef struct CSamInfo CSam;

typedef vector<CBed *> bedVector;
typedef map<string, bedVector> chromBedMap;

typedef vector<CBed6 *> bed6Vector;
typedef map<string, bed6Vector> chromBed6Map;

typedef vector<CSam *> samVector;
typedef map<string, samVector> chromSamMap;

typedef map<string, double *> chromSignalMap;
typedef map<string, int> chromSizeMap;

string BuildCigarString(const vector<CigarOp> &cigar);

string PrintTag(const BamAlignment &bam, const string &tag);

void openBamFile(char *bamFile, BamReader &reader);

void openBamNoIdxFile(char *bamFile, BamReader &reader);

int readBamToLocusNum(BamReader &reader, map<string, int> &readLocus);

int readSamToLocusNum(FILE *fp, map<string, int> &readLocus);

int getMDZ(char **fields, int num);

int getEndPos(char *cigar);

long getGenomeSize(BamReader & reader, chromSizeMap &mapSize);

double readBamToSignals(BamReader &reader, chromSignalMap &sigHash, int keepDup, int normalization);

double readBamToBedSignalMap(BamReader &reader, chromBedMap &bedHash, int keepDup, int maxLocusNum);

double readBamToBed6Map(BamReader &reader, chromBed6Map &bed6Hash, int keepDup=0, int maxInsertLen=1000000, int skipSplice=0, int ccaTail=0);

double readBamToBedMap(BamReader &reader, chromBedMap &bedHash, map<string, int> &readLocus, int keepDup, int maxLocusNum);

double readBamToSamMap(BamReader &reader, chromSamMap &samHash, map<string, int> &readLocus, int keepDup, int maxLocusNum);

double readBamToSamMapNew(BamReader &reader, chromSamMap &samHash, int keepDup, int maxLocusNum, int skipSplice);

double readPEBamToBed6Map(BamReader &reader, chromBed6Map &bed6Hash, int maxInsertLen);

double readPEBamToBedMap(BamReader &reader, chromBedMap &bedHash, map<string, int> &readLocus, int keepDup, int maxLocusNum, int maxInsertLen, int readLen);

double readPEBamSortedNameToBedMap(BamReader &reader, chromBedMap &bedHash, map<string, int> &readLocus,
                                   int keepDup, int maxLocusNum, int maxInsertLen, int readLen);

double normalizedSamReads(chromSamMap &samHash);

int getRegionCounts(BamReader &reader, char *chrom, int chromStart, int chromEnd, char strand);

int filterLowQualities(int qualityBase, int minScore, int minPercent, const char *qualities);

void freeChromBed6Map(chromBed6Map &bed6Hash) /*free bed map */;

void freeChromBedMap(chromBedMap &bedHash);

void freeChromSignalMap(chromSignalMap &sigHash);

void freeBedVector(bedVector &bedList);

void freeSamVector(samVector &samList);

void freeSamItem(CSam *sam);

void freeChromSamMap(chromSamMap &samHash);

bool compareSam(const CSam *oSam, const CSam *tSam);

int compareSamReadNum(const void *a, const void *b);

void copySam(CSam *tSam, CSam *oSam);

CBed6 *bamToBed6(BamAlignment &bam, string &chrom);

void freeChromBed12Map(chromBed12Map &bed12Hash) ;

void bamToBlocks(const BamAlignment &bam, const string &chrom, bed6Vector &bedBlocks, bool skipDeletion, bool skipSplice);

CBed12 *bamToBed12(const BamAlignment &bam, const string &chrom, int skipDeletion=0, int skipSplice=0, int collapser=0);

double readSEBamToBed12Map(BamReader &reader, chromBed12Map &bed12Hash, int skipDeletion=0, int skipSplice=0, int collapser=0, int skipBigClip=0, int skipSoft=0);

double readPEBamToBed6MapNew(BamReader &reader, chromBed6Map &bed6Hash, int maxInsertLen=1000000, int skipSplice=0, int ccaTail=0);

double readPEBamToBed12Map(BamReader &reader, chromBed12Map &bed12Hash, int skipDeletion=0, int skipSplice=0, int collapser=0, int skipBigClip=0, int skipSoft=0);

int getAllClipLen(string &cigarData, int *leftClipLen, int *rightClipLen);

#endif /* End BAM_HEAD_H */
