/*
	typedef struct 
	{
		// All
		unsigned long long int 	PMExecutionDetailsExecutionStart;
		unsigned long long int 	PMExecutionDetailsExecutionEnd;
		float			PMExecutionDetailsExecutionTime;

		// Unique Patterns
		unsigned long long int	PMExecutionDetailsPatternsStart;
		unsigned long long int	PMExecutionDetailsPatternsEnd;
		float			PMExeuctionDetailsPatternsTime;

		// Unique Patterns
		unsigned long long int	PMExecutionDetailsLDMatrixStart;
		unsigned long long int	PMExecutionDetailsLDMatrixEnd;
		float			PMExeuctionDetailsLDMatrixTime;



	} PMExecutionDetails;

*/

#include "PM.h"

PMExecutionDetails * PMExecutionDetails_new (void)
{
	PMExecutionDetails * ed = NULL;
	
	ed = (PMExecutionDetails*)malloc(sizeof(PMExecutionDetails));
	assert(ed!=NULL);

	ed->PMExecutionDetailsExecutionStart = 0llu;
	ed->PMExecutionDetailsExecutionEnd = 0llu;
	ed->PMExecutionDetailsExecutionTime = 0.0f;

	ed->PMExecutionDetailsPatternsStart = 0llu;
	ed->PMExecutionDetailsPatternsEnd = 0llu;
	ed->PMExecutionDetailsPatternsTime = 0.0f;

	ed->PMExecutionDetailsLDMatrixStart = 0llu;
	ed->PMExecutionDetailsLDMatrixEnd = 0llu;
	ed->PMExecutionDetailsLDMatrixTime = 0.0f;

	return ed;
}

void PMExecutionDetails_free (PMExecutionDetails * ed)
{
	assert(ed!=NULL);

	free(ed);
}

void PMExecutionDetails_reset (PMExecutionDetails * ed)
{
	assert(ed!=NULL);

	ed->PMExecutionDetailsExecutionStart = 0llu;
	ed->PMExecutionDetailsExecutionEnd = 0llu;
	ed->PMExecutionDetailsExecutionTime = 0.0f;

	ed->PMExecutionDetailsPatternsStart = 0llu;
	ed->PMExecutionDetailsPatternsEnd = 0llu;
	ed->PMExecutionDetailsPatternsTime = 0.0f;

	ed->PMExecutionDetailsLDMatrixStart = 0llu;
	ed->PMExecutionDetailsLDMatrixEnd = 0llu;
	ed->PMExecutionDetailsLDMatrixTime = 0.0f;

}

void PMExecutionDetails_startExecutionTimer (PMExecutionDetails * ed)
{
	assert(ed!=NULL);

	ed->PMExecutionDetailsExecutionStart = rdtsc();
}

void PMExecutionDetails_pauseExecutionTimer (PMExecutionDetails * ed)
{
	assert(ed!=NULL);

	ed->PMExecutionDetailsExecutionEnd = rdtsc();
	
	ed->PMExecutionDetailsExecutionTime+= (ed->PMExecutionDetailsExecutionEnd-ed->PMExecutionDetailsExecutionStart)/1.3e9;
}

float PMExecutionDetails_readExecutionTimer (PMExecutionDetails * ed)
{
	assert(ed!=NULL);

	return ed->PMExecutionDetailsExecutionTime;
}

void PMExecutionDetails_startPatternsTimer (PMExecutionDetails * ed)
{
	assert(ed!=NULL);

	ed->PMExecutionDetailsPatternsStart = rdtsc();
}

void PMExecutionDetails_pausePatternsTimer (PMExecutionDetails * ed)
{
	assert(ed!=NULL);

	ed->PMExecutionDetailsPatternsEnd = rdtsc();
	
	ed->PMExecutionDetailsPatternsTime+= (ed->PMExecutionDetailsPatternsEnd-ed->PMExecutionDetailsPatternsStart)/1.3e9;
}

float PMExecutionDetails_readPatternsTimer (PMExecutionDetails * ed)
{
	assert(ed!=NULL);

	return ed->PMExecutionDetailsPatternsTime;
}

void PMExecutionDetails_startLDMatrixTimer (PMExecutionDetails * ed)
{
	assert(ed!=NULL);

	ed->PMExecutionDetailsLDMatrixStart = rdtsc();
}

void PMExecutionDetails_pauseLDMatrixTimer (PMExecutionDetails * ed)
{
	assert(ed!=NULL);

	ed->PMExecutionDetailsLDMatrixEnd = rdtsc();
	
	ed->PMExecutionDetailsLDMatrixTime+= (ed->PMExecutionDetailsLDMatrixEnd-ed->PMExecutionDetailsLDMatrixStart)/1.3e9;
}

float PMExecutionDetails_readLDMatrixTimer (PMExecutionDetails * ed)
{
	assert(ed!=NULL);

	return ed->PMExecutionDetailsLDMatrixTime;
}

void PMExecutionDetails_print (PMExecutionDetails * ed, FILE * fp)
{
	assert(ed!=NULL);
	assert(fp!=NULL);

	fprintf(fp, "\nEXECUTION TIME BREAKDOWN\n");

	float totalTime = -1.0;
	totalTime = PMExecutionDetails_readExecutionTimer(ed);

	float patternsTime = -1.0;
	patternsTime = PMExecutionDetails_readPatternsTimer(ed);
	fprintf(fp, "[patterns] %f (%.2f%%)\n", patternsTime, (patternsTime/totalTime)*100.0);

	float ldTime = -1.0;
	ldTime = PMExecutionDetails_readLDMatrixTimer(ed);
	fprintf(fp, "[ld] %f (%.2f%%)\n", ldTime, (ldTime/totalTime)*100.0);


	fprintf(fp, "[total] %f\n", totalTime);
}
