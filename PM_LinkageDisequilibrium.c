#include "PM.h"

// Binary no gaps Reference
float LD_ref (PMPatternPool * pp, int index1, int index2)
{
	int i, K = pp->PMPatternPoolCompactSize;
	float sample_size = 1.0/((float)pp->PMPatternPoolPloidy*pp->PMPatternPoolSampleSize);

	int c1 = pp->PMPatternPoolDerivedAlleleCount[index1];
	int c2 = pp->PMPatternPoolDerivedAlleleCount[index2];
	int c3 = 0;
	int c4 = 0;
	for(i=0;i<K;i=i+PADDING)
	{
		c3 += pm_popcnt_u64(pp->PMPatternPoolDataCompact[index1][i]&pp->PMPatternPoolDataCompact[index2][i]);
		c4 += pm_popcnt_u64(pp->PMPatternPoolDataCompact[index1][i+1]&pp->PMPatternPoolDataCompact[index2][i+1]);
		c3 += pm_popcnt_u64(pp->PMPatternPoolDataCompact[index1][i+2]&pp->PMPatternPoolDataCompact[index2][i+2]);
		c4 += pm_popcnt_u64(pp->PMPatternPoolDataCompact[index1][i+3]&pp->PMPatternPoolDataCompact[index2][i+3]);
	}	
	float p1 = c1*sample_size;
	float p2 = c2*sample_size;
	float x = (c3+c4)*sample_size;

	float Num = (x-p1*p2)*(x-p1*p2);
	float Den = p1*p2*(1.0-p1)*(1.0-p2);
	return Num/Den;
}

// Binary no gaps Version 2
float LD_test1 (PMPatternPool * pp, int index1, int index2)
{
	int i, K = pp->PMPatternPoolCompactSize;
	float sample_size = (float)pp->PMPatternPoolPloidy*pp->PMPatternPoolSampleSize;

	int c1 = pp->PMPatternPoolDerivedAlleleCount[index1];
	int c2 = pp->PMPatternPoolDerivedAlleleCount[index2];
	int c3 = 0;
	int c4 = 0;
	for(i=0;i<K;i=i+PADDING)
	{
		c3 += pm_popcnt_u64(pp->PMPatternPoolDataCompact[index1][i]&pp->PMPatternPoolDataCompact[index2][i]);
		c4 += pm_popcnt_u64(pp->PMPatternPoolDataCompact[index1][i+1]&pp->PMPatternPoolDataCompact[index2][i+1]);
		c3 += pm_popcnt_u64(pp->PMPatternPoolDataCompact[index1][i+2]&pp->PMPatternPoolDataCompact[index2][i+2]);
		c4 += pm_popcnt_u64(pp->PMPatternPoolDataCompact[index1][i+3]&pp->PMPatternPoolDataCompact[index2][i+3]);
	}	
	float H1 = (float)(c3+c4);
	float Den = (c1*c2*(sample_size-c1)*(sample_size-c2));
	float Num = (sample_size*H1-c1*c2)*(sample_size*H1-c1*c2);
	return Num/Den;	
}

