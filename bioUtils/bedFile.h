#ifndef BED_HEAD_H
#define BED_HEAD_H

struct CBedInfo
{
        struct CBedInfo *next; // for list
        char *chrom;    // chromosome name
        int chromStart; // chromosome start
        int chromEnd; // chromosome end
        double score; // normalized read number score
        char strand; // strand
        int mateFlag; // mate flag
        int locusNum; // locus number from NH flag in STAR program
};

typedef struct CBedInfo CBed;

struct CBed6Info
{
        struct CBed6Info *next; // for list
        char *chrom;    // chromosome name
        int chromStart; // chromosome start
        int chromEnd; // chromosome end
        char *name;   /* Name of item */
        double score; // normalized read number score
        char strand; // strand
        int mateFlag; // mate flag
        int locusNum; // locus number from NH flag in STAR program
};

typedef struct CBed6Info CBed6;
typedef vector<CBed6 *> bed6Vector;
typedef map<string, bed6Vector> chromStrandBed6Map;
typedef map<string, bed6Vector> chromBed6Map;

struct Bed12Info
{
        struct CBed12Info *next; // for list
        char *chrom;    // chromosome name
        int chromStart; // chromosome start
        int chromEnd; // chromosome end
        char *name;     /* Name of item */
        double score; // normalized read number score
        char strand; // strand
        int thickStart; /* Start of where display should be thick (start codon for genes) */
        int thickEnd;   /* End of where display should be thick (stop codon for genes) */
        int itemRgb;    /* RGB 8 bits each */
        int blockCount; /* Number of blocks. */
        int *blockSizes;     /* Comma separated list of block sizes.  */
        int *chromStarts;    /* Start positions inside chromosome.  Relative to chromStart*/
};

typedef struct Bed12Info CBed12;

typedef vector<CBed12 *> bed12Vector;
typedef map<string, bed12Vector> chromStrandBed12Map;
typedef map<string, bed12Vector> chromBed12Map;
typedef map<string, CBed12 *> nameBed12Map;

typedef map<string, int> chromSizeMap;

struct MateBedInfo
{
        struct MateBedInfo *next; // for list
        char *chrom1;    // chromosome name
        int chromStart1; // chromosome start
        int chromEnd1; // chromosome end
        double score1; // normalized read number
        char strand1; // strand
        char *chrom2;    // chromosome name
        int chromStart2; // chromosome start
        int chromEnd2; // chromosome end
        double score2; // normalized read number
        char strand2; // strand
};

typedef struct MateBedInfo MateBed;

void freeBedItem(CBed *bed);

void freeBed6Item(CBed6 *bed6);

void freeBed12Item(CBed12 *bed12);

void freeBedMapList(map<string, CBed *> &bedHash);
/*free bed list */

void freeBedList(CBed *bedList);
/*free bed list */

void sortBed(CBed **pList, int count, int (*compare )(const void *item1,  const void *item2));

void copyBed6(CBed6 *tBed, CBed6 *oBed);

void copyBed(CBed *tBed, CBed *oBed);

bool compareBed( const CBedInfo *oBed, const CBedInfo *tBed);

bool compareBed6( const CBed6Info *oBed, const CBed6Info *tBed);

int compareBedScore(const void *a, const void *b);

CBed12 *parseBed12Line(char *line);

CBed6 *parseBed6Line(char *line);

int getBed6ToMap(FILE *bedfp, chromStrandBed6Map &bed6Hash);

void freeChromStrandBed6Map(chromStrandBed6Map &bedHash);

void freeBed6Vector(bed6Vector &bedList);

void freeNameBed12Map(nameBed12Map &bedHash);

int getBed12ToNameMap(FILE *bedfp, nameBed12Map &bed12NameHash);

int getBedToMap(FILE *bedfp, chromStrandBed12Map &bed12Hash);

void freeChromStrandMap(chromStrandBed12Map &bedHash);

void freeBed12Vector(bed12Vector &bedList);

double readBedToBed6Map(FILE *bedfp, chromBed6Map &bed6Hash);

long getChromSize(FILE *gfp, chromSizeMap &mapSize);

void warnBed12(CBed12 *bed12);

CBed12 *twoBed6ToBed12(CBed6 *bed1, CBed6 *bed2);

CBed6 *createBed6(const string &chrom, int chromStart, int chromEnd,  const string &name, double score, char strand);

CBed6 *mergeTwoBed6(CBed6 *bed1, CBed6 *bed2);

CBed12 *mergeTwoBed12ToBed12(CBed12 *bed1, CBed12 *bed2);

double normalizedReadsByNH(chromBed6Map &bed6Hash);

double normalizedBed6Reads(chromBed6Map &bed6Hash);

double normalizedBed12Reads(chromBed12Map &bed12Hash);

double removeMisprimReadBed12(char *primerSeq, int barcodeLen,
                             FILE *genomefp,
                             faidxMap &faiHash,
                             chromBed12Map &bedHash,
                             chromBed12Map &newBedHash);

#endif /* End BED_HEAD_H */
