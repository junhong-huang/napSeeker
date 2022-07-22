/*
 fasta.c
 Functions for module fasta 
 $Id: fasta.c,v 1.0 2007/04/03/ 09:46 jianhua Exp $
*/
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
#include "fasta.h"

int readSeq(FILE *fp,FASTA *fastaSeq)
/* read fasta sequence */
{
	char *name;
  char delims[] = " \r\t\n";	
  char *seq;
  char *pSeq = NULL;
  int  seqLen;
  char c;  
  long size = 512;
  long firstAllocBufSize = 128*1024;
  char lineBuf[1024];
   
	if (feof(fp)) {
		return 0;
		}
	
	while(fgets(lineBuf, sizeof(lineBuf), fp) != NULL) {
		if (lineBuf[0]=='>') {
		name = strtok(lineBuf+1, delims);	    
fastaSeq->seqName = (char *)safeMalloc(strlen(name)+1);
strcpy(fastaSeq->seqName, name); 
			break;
			}
		}
				
	seq = (char *)safeMalloc(firstAllocBufSize+1);
	seqLen = 0;
	size = firstAllocBufSize;

/* read sequence */	
	for(;;){
		c = fgetc(fp);		
		if (c==EOF || c=='>')
			break;			
		if (seqLen>=size) {
			size = size + size;
			pSeq = seq;
			seq = (char *) safeMalloc(size);
			memcpy(seq, pSeq, seqLen);
			safeFree(pSeq);
			pSeq = NULL;
			}
	  if (isalpha(c)){
		seq[seqLen] = c;	
		seqLen++;
	   }			
	  }

seq[seqLen] = '\0';
fastaSeq->seqLen = seqLen;
fastaSeq->seq = seq;

if (c=='>') {
	ungetc(c,fp);
	}	
return 1;
	}

void freeFasta(FASTA *fasta) {
	if (fasta != NULL) {
		safeFree(fasta->seqName);
		safeFree(fasta->seq);
		}
	fasta = NULL;
	}



