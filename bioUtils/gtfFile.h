#ifndef GTF_HEAD_H
#define GTF_HEAD_H

/* modified from kent ucsc codes */
struct gtfInfo
/* A parsed line in a GTF file. */
{
        struct gtfInfo *next;  /* Next line in file */
        char *seqName;      /* Name of sequence. */
        char *source;   /* Program that made this line.  Not allocated here. */
        char *feature;  /* Type field. (Intron, CDS, etc). Not allocated here. */
        int start;      /* Start of feature in sequence. Starts with 0, not 1 */
        int end;        /* End of feature in sequence. End is not included. */
        double score;   /* Score. */
        char strand;    /* Strand of sequence feature is on. + or - or .*/
        char frame;     /* Frame feature is in. 1, 2, 3, or . */
        char *geneId;    /* gene_id in GTF, NULL in GFF.  Not allocated here. */
        char *transcriptId;       /* transcript_id in GTF, NULL in GFF. Not allocated here. */
        char *exonId;       /* exon_id in GTF, NULL in GFF. Not allocated here. */
        int exonNumber; /* O in GFF or if missing in GTF.  Otherwise exon number. */
        char *geneName;       /* gene_name or NULL in GTF, NULL in GFF. Not allocated here. */
        char *transcriptName; /* transcript_name or NULL in GTF, NULL in GFF. Not allocated here. */
};

typedef struct gtfInfo gtfLine;

typedef vector<gtfLine *> gtfVector;

gtfLine *parseGtfLine(char *line);

int getGtfToBed12Map(FILE *gtffp, chromStrandBed12Map &bed12Hash, const char *type);

CBed12 *mergeGtfToBed12(gtfVector &gtfList, const char *type);

void freeGtfVector(gtfVector &gtfList);

void freeGtf(gtfLine *gtf);

#endif /* End GTF_HEAD_H */
