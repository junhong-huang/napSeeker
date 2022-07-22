/* statistic head file */

#ifndef STATISTIC_H
#define STATISTIC_H

double log10Val(double pval);

double log10Lnp(double lnP);

double log10Pval(double pval);

double simpson(int overlapRegionLen, int queryRegionLen, int sampleRegionLen);

double jaccard(int overlapRegionLen, int queryRegionLen, int sampleRegionLen);

double binomialPval(double s, double xn, double pr);

double poissonPval(double s, double lambda);

double poissonSval(double q, double lambda);

long double hypergeometric(int n, long double p, int k, int r);

long double lFactorial(int number) ;

long double lNchooseK(int n, int k);

#endif /* End STATISTIC_H */