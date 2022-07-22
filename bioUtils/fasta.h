/*
 fasta.h
 Header for module fasta 
 $Id: fasta.h,v 1.0 2007/03/29/ 11:59 jianhua Exp $
*/

#ifndef FASTA_H
#define FASTA_H

struct FastaInfo{
	char *seqName;
	long seqLen;
	char *seq;
	};
typedef struct FastaInfo FASTA;

extern int readSeq(FILE *fp,FASTA *fastaSeq);
/*read Fasta format Sequence*/
extern void freeFasta(FASTA *fastaSeq);
/*free Fasta sequence */

#endif /* End FASTA_H */
