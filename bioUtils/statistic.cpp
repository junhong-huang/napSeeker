#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<ctype.h>
#include <getopt.h>
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
#include "cdflib.hpp"
#include "statistic.h"

double log10Val(double pval)
{
	double returnVal = 0;
	if (pval == 0)
	{
		returnVal = -325.00;
	}
	else {
		returnVal = log(pval) / log(10);
	}
	return returnVal;
}

double log10Lnp(double lnP)
{
	double returnVal = 0;
	returnVal = lnP / log(10);
	if (returnVal > 0) returnVal = 0;
	return returnVal;
}

double log10Pval(double pval)
{
	double returnVal = 0;
	if (pval == 0)
	{
		returnVal = -325.00;
	}
	else {
		returnVal = log(pval) / log(10);
	}
	if (returnVal > 0) returnVal = 0;
	return returnVal;
}

double simpson(int overlapRegionLen, int queryRegionLen, int sampleRegionLen) {
	double 	simpson = 0;
	int minVal = queryRegionLen;
	if (minVal > sampleRegionLen) minVal = sampleRegionLen;
	simpson = overlapRegionLen / (double)minVal;
	return simpson;
}

double jaccard(int overlapRegionLen, int queryRegionLen, int sampleRegionLen)
{
	double jaccardVal = 0;
	jaccardVal = overlapRegionLen / (double)(queryRegionLen + sampleRegionLen - overlapRegionLen);
	return jaccardVal;
}

double binomialPval(double s, double xn, double pr) {
	int    which = 1;
	double p = 0.5;
	double q = 0.5;
	int status = 2;
	double bound = 0;
	double ompr = 1 - pr;
	cdfbin(&which, &p, &q, &s, &xn, &pr, &ompr, &status, &bound);
	return q;
}

double poissonPval(double s, double lambda) {
	int    which = 1;
	double p = 0.5;
	double q = 0.5;
	int status = 2;
	double bound = 0;
	cdfpoi(&which, &p, &q, &s, &lambda, &status, &bound);
	return q;
}

double poissonSval(double q, double lambda) {
	int    which = 2;
	double p = 1.0 - q;
	int status = 2;
	double bound = 0;
	double s = 0;
	cdfpoi(&which, &p, &q, &s, &lambda, &status, &bound);
	return s;
}

long double hypergeometric(int n, long double p, int k, int r)
{
	/*
	 n = refgene number (background number)
	 p = gene number in goterm/refgene-number
	 k = gene number picked by methods (gene number identified from your study)
	 r = gene number picked in goterm(overlapping gene number between k and gene-number from go term)
	 */
	long double q;
	int np;
	int nq;
	int top;
	int i;
	long double logNchooseK;
	long double lfoo;
	long double sum;
	q = (1 - p);
	np = floor((long double) n * p + 0.5);
	nq = floor((long double) n * q + 0.5);
	logNchooseK = lNchooseK(n, k);
	top = k;
	if (np < k) {
		top = np;
	}
	lfoo = lNchooseK(np, top) + lNchooseK(nq, k - top);
//fprintf(stderr, "%d %.5f %d %d\n", n, lfoo, k, r);
	sum = 0;

	for (i = top; i >= r; i--) {
		sum = sum + exp(lfoo - logNchooseK);
		if (i > r) {
			lfoo = lfoo + log((long double) i / (long double) (np - i + 1))
			       + log(
			           (long double) (nq - k + i)
			           / (long double) (k - i + 1));
		}
	}
	return sum;
}

// ln factorial subroutine
long double lFactorial(int number) {
	long double returnValue = 0;
	for (int i = 2; i <= number; i++) {
		returnValue = returnValue + log(i);
	}
	return returnValue;
}

// ln N choose K subroutine
long double lNchooseK(int n, int k) {
	long double answer = 0;
	int i = 0;
	if (k > (n - k)) {
		k = (n - k);
	}
	for (i = n; i > (n - k); i--) {
		answer = answer + log(i);
	}
	answer = answer - lFactorial(k);
	return answer;
}
