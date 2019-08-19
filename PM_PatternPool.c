/*

	typedef struct
	{
		int 		PMPatternPoolSize; // Number of unique patterns
		int		PMPatternPoolSampleSize; 
		int		PMPatternPoolPloidy;
		float		PMPatternPoolLDThreshold;

		int		PMPatternPoolCompactSize; // number of uint64_t per SNP 
		char	*	PMPatternPoolPendingSNP; // temp buffer to be checked against pattern pool
		uint64_t *  	PMPatternPoolPendingSNPCompact;
		uint64_t ** 	PMPatternPoolDataCompact;

		int	*	PMPatternPoolCount; // number of SNPs per unique pattern
		int	**	PMPatternPoolPositions; // the SNP positions
		char	***	PMPatternPoolIDs; // the SNP IDs
		int	**	PMPatternPoolSNPIndeces; // this is the actual SNP index in the file/region
		int	*	PMPatternPoolDerivedAlleleCount; // number of 1s per unique pattern

		PMWindow *	PMPatternPoolWindow; // corresponding genomic window for loaded pattern pool

		float **	PMPatternPoolLDMatrix;
		char   *	PMPatternPoolInLDClass;
		int		PMPatternPoolLDClassListSize;
		PMLDClass **	PMPatternPoolLDClassList;
	

	} PMPatternPool;

*/

#include "PM.h"

PMPatternPool * PMPatternPool_new (void)
{
	PMPatternPool * pp = NULL;
	pp = (PMPatternPool *) malloc(sizeof(PMPatternPool));
	assert(pp!=NULL);

	pp->PMPatternPoolSize = 0;
	pp->PMPatternPoolSampleSize = 0;
	pp->PMPatternPoolPloidy = -1;
	pp->PMPatternPoolLDThreshold = 0.0f;
	
	pp->PMPatternPoolCompactSize = 0;
	pp->PMPatternPoolPendingSNP = NULL;
	pp->PMPatternPoolPendingSNPCompact = NULL;
	pp->PMPatternPoolDataCompact = NULL;

	pp->PMPatternPoolCount = NULL;
	pp->PMPatternPoolPositions = NULL;
	pp->PMPatternPoolIDs = NULL;
	pp->PMPatternPoolSNPIndeces = NULL;
	pp->PMPatternPoolDerivedAlleleCount = NULL;

	pp->PMPatternPoolWindow = NULL;

	pp->PMPatternPoolLDMatrix = NULL;
	pp->PMPatternPoolInLDClass = NULL;
	pp->PMPatternPoolLDClassListSize = 0;
	pp->PMPatternPoolLDClassList = NULL;

	return pp;
}

void PMPatternPool_free (PMPatternPool * pp)
{
	assert(pp!=NULL);

	// Pending SNP
	if(pp->PMPatternPoolPendingSNP!=NULL)
		free(pp->PMPatternPoolPendingSNP);

	// Pending SNP Compact
	if(pp->PMPatternPoolPendingSNPCompact!=NULL)
		free(pp->PMPatternPoolPendingSNPCompact);

	// Data Compact
	if(pp->PMPatternPoolDataCompact!=NULL)
	{
		int i;
		for(i=0;i<PATTERNS_MAX;i++)
			free(pp->PMPatternPoolDataCompact[i]);

		free(pp->PMPatternPoolDataCompact);	
	}

	// Positions
	if(pp->PMPatternPoolPositions!=NULL)
	{
		int i;
		for(i=0;i<PATTERNS_MAX;i++)
		{
			if(pp->PMPatternPoolPositions[i]!=NULL)
				free(pp->PMPatternPoolPositions[i]);
		}
		free(pp->PMPatternPoolPositions);
	}

	// IDs
	if(pp->PMPatternPoolIDs!=NULL)
	{
		int i;
		for(i=0;i<PATTERNS_MAX;i++)
		{
			int j;
			for(j=0;j<pp->PMPatternPoolCount[i];j++)
			{
				if(pp->PMPatternPoolIDs[i][j]!=NULL)
					free(pp->PMPatternPoolIDs[i][j]);
			}
			free(pp->PMPatternPoolIDs[i]);
			
		}
		free(pp->PMPatternPoolIDs);
	}

	// Indeces
	if(pp->PMPatternPoolSNPIndeces!=NULL)
	{
		int i;
		for(i=0;i<PATTERNS_MAX;i++)
		{
			if(pp->PMPatternPoolSNPIndeces[i]!=NULL)
				free(pp->PMPatternPoolSNPIndeces[i]);
		}
		free(pp->PMPatternPoolSNPIndeces);
	}

	// Count
	if(pp->PMPatternPoolCount!=NULL)
		free(pp->PMPatternPoolCount);


	// Derived Allele Count
	if(pp->PMPatternPoolDerivedAlleleCount!=NULL)
		free(pp->PMPatternPoolDerivedAlleleCount);


	// Window
	if(pp->PMPatternPoolWindow!=NULL)
		PMWindow_free(pp->PMPatternPoolWindow);

	// LD Matrix
	if(pp->PMPatternPoolLDMatrix!=NULL)
	{
		int i;
		for(i=0;i<PATTERNS_MAX;i++)
		{
			if(pp->PMPatternPoolLDMatrix[i]!=NULL)
				free(pp->PMPatternPoolLDMatrix[i]);	
		}
		free(pp->PMPatternPoolLDMatrix);
	}

	// In LD Class
	if(pp->PMPatternPoolInLDClass!=NULL)
		free(pp->PMPatternPoolInLDClass);

	// LD Class List
	if(pp->PMPatternPoolLDClassList!=NULL)
	{
		int i;
		for(i=0;i<pp->PMPatternPoolLDClassListSize;i++)
		{
			if(pp->PMPatternPoolLDClassList[i]!=NULL)
				PMLDClass_free(pp->PMPatternPoolLDClassList[i]);
		}
		free(pp->PMPatternPoolLDClassList);
		pp->PMPatternPoolLDClassListSize = 0;		
	}

	free(pp);
}

int PMPatternPool_allocate (PMPatternPool * pp, PMInput * i, PMCommandLine * cl)
{
	assert(pp!=NULL);
	assert(i!=NULL);

	int SampleSize = PMInput_getSampleSize(i);
	PMPatternPool_setSampleSize(pp, SampleSize);

	PMPatternPool_setPloidy (pp, cl);
	PMPatternPool_setLDThreshold (pp, cl);

	// Pending SNP
	pp->PMPatternPoolPendingSNP = realloc(pp->PMPatternPoolPendingSNP, sizeof(char)*(pp->PMPatternPoolPloidy*SampleSize+1));
	assert(pp->PMPatternPoolPendingSNP!=NULL);
	pp->PMPatternPoolPendingSNP[SampleSize]='\0';
		
	int j;
	for(j=0;j<SampleSize;j++)
		pp->PMPatternPoolPendingSNP[j] = 'X';

	// Pending SNP Compact
	int entriesPerSNP = pp->PMPatternPoolPloidy*SampleSize;
	int compactLength = (entriesPerSNP%WORD_WIDTH!=0) ? entriesPerSNP/WORD_WIDTH+1 : entriesPerSNP/WORD_WIDTH;
	pp->PMPatternPoolCompactSize = compactLength;

	int incr = 0;
	if(compactLength%PADDING!=0)
		incr = PADDING-(compactLength%PADDING);

	pp->PMPatternPoolPendingSNPCompact = realloc(pp->PMPatternPoolPendingSNPCompact, sizeof(uint64_t)*(compactLength+incr));
	assert(pp->PMPatternPoolPendingSNPCompact!=NULL);

	for(j=0;j<compactLength+incr;j++)
		pp->PMPatternPoolPendingSNPCompact[j] = 0llu;

	// Data Compact
	pp->PMPatternPoolDataCompact = realloc(pp->PMPatternPoolDataCompact, sizeof(uint64_t*)*PATTERNS_MAX);
	assert(pp->PMPatternPoolDataCompact!=NULL);

	for(j=0;j<PATTERNS_MAX;j++)
	{
		pp->PMPatternPoolDataCompact[j] = (uint64_t*)malloc(sizeof(uint64_t)*(compactLength+incr));
		assert(pp->PMPatternPoolDataCompact[j]!=NULL);

		int k;
		for(k=0;k<compactLength+incr;k++)
			pp->PMPatternPoolDataCompact[j][k] = 0llu;
	}

	// Positions
	pp->PMPatternPoolPositions = realloc(pp->PMPatternPoolPositions, sizeof(int*)*PATTERNS_MAX);
	assert(pp->PMPatternPoolPositions!=NULL);

	for(j=0;j<PATTERNS_MAX;j++)
		pp->PMPatternPoolPositions[j] = NULL;

	// IDs
	pp->PMPatternPoolIDs = realloc(pp->PMPatternPoolIDs, sizeof(char**)*PATTERNS_MAX);
	assert(pp->PMPatternPoolIDs!=NULL);

	for(j=0;j<PATTERNS_MAX;j++)
		pp->PMPatternPoolIDs[j] = NULL;	

	// Indeces 
 	pp->PMPatternPoolSNPIndeces = realloc(pp->PMPatternPoolSNPIndeces, sizeof(int*)*PATTERNS_MAX);
	assert(pp->PMPatternPoolSNPIndeces!=NULL);

	for(j=0;j<PATTERNS_MAX;j++)
		pp->PMPatternPoolSNPIndeces[j] = NULL;

	// Count
	pp->PMPatternPoolCount = realloc(pp->PMPatternPoolCount, sizeof(int)*PATTERNS_MAX);
	assert(pp->PMPatternPoolCount!=NULL);

	for(j=0;j<PATTERNS_MAX;j++)
		pp->PMPatternPoolCount[j] = 0;

	// Derived Allele Count
	pp->PMPatternPoolDerivedAlleleCount = realloc(pp->PMPatternPoolDerivedAlleleCount, sizeof(int)*PATTERNS_MAX);
	assert(pp->PMPatternPoolDerivedAlleleCount!=NULL);

	for(j=0;j<PATTERNS_MAX;j++)
		pp->PMPatternPoolDerivedAlleleCount[j] = 0;

	// Window
	pp->PMPatternPoolWindow = PMWindow_new();
	assert(pp->PMPatternPoolWindow!=NULL);

	// LD Matrix
	pp->PMPatternPoolLDMatrix = realloc(pp->PMPatternPoolLDMatrix, sizeof(float*)*PATTERNS_MAX);
	assert(pp->PMPatternPoolLDMatrix!=NULL);
	for(j=0;j<PATTERNS_MAX;j++)
	{
		pp->PMPatternPoolLDMatrix[j] = (float*)malloc(sizeof(float)*PATTERNS_MAX);
		assert(pp->PMPatternPoolLDMatrix[j]!=NULL);
		int k;
		for(k=0;k<PATTERNS_MAX;k++)
			pp->PMPatternPoolLDMatrix[j][k] = 0.0f;	
	}
	
	// In LD class
	pp->PMPatternPoolInLDClass = realloc(pp->PMPatternPoolInLDClass, sizeof(char)*PATTERNS_MAX);
	assert(pp->PMPatternPoolInLDClass!=NULL);
	for(j=0;j<PATTERNS_MAX;j++)
		pp->PMPatternPoolInLDClass[j] = FALSE;

	return SUCCESS;	
}

void PMPatternPool_reset (PMPatternPool * pp)
{
	assert(pp!=NULL);

	int SampleSize = PMPatternPool_getSampleSize(pp);
	pp->PMPatternPoolPendingSNP[SampleSize]='\0';

	PMPatternPool_setSize (pp, 0);
		
	int j;
	// Pending SNP
	for(j=0;j<SampleSize;j++)
		pp->PMPatternPoolPendingSNP[j] = 'X';	

	// Pending SNP Compact
	for(j=0;j<pp->PMPatternPoolCompactSize;j++)
		pp->PMPatternPoolPendingSNPCompact[j] = 0llu;	

	// Data
	for(j=0;j<PATTERNS_MAX;j++)
	{
		int k;
		for(k=0;k<pp->PMPatternPoolCompactSize;k++)
			pp->PMPatternPoolDataCompact[j][k] = 0llu;

	}

	// Positions
	for(j=0;j<PATTERNS_MAX;j++)
	{
		if(pp->PMPatternPoolPositions[j]!=NULL)
			free(pp->PMPatternPoolPositions[j]);

		pp->PMPatternPoolPositions[j] = NULL;
	}

	// IDs
	for(j=0;j<PATTERNS_MAX;j++)
		pp->PMPatternPoolIDs[j] = NULL;	

	// Indeces PMPatternPoolSNPIndeces
	for(j=0;j<PATTERNS_MAX;j++)
	{
		if(pp->PMPatternPoolSNPIndeces[j]!=NULL)
			free(pp->PMPatternPoolSNPIndeces[j]);

		pp->PMPatternPoolSNPIndeces[j] = NULL;
	}

	// Count
	for(j=0;j<PATTERNS_MAX;j++)
		pp->PMPatternPoolCount[j] = 0;

	// Derived Allele Count
	for(j=0;j<PATTERNS_MAX;j++)
		pp->PMPatternPoolDerivedAlleleCount[j] = 0;

	// Window
	if(pp->PMPatternPoolWindow!=NULL)
		PMWindow_free(pp->PMPatternPoolWindow);

	pp->PMPatternPoolWindow = PMWindow_new();
	assert(pp->PMPatternPoolWindow!=NULL);

	// LD Matrix
	for(j=0;j<PATTERNS_MAX;j++)
	{
		int k;
		for(k=0;k<PATTERNS_MAX;k++)
			pp->PMPatternPoolLDMatrix[j][k] = 0.0f;	
	}

	// In LD Class
	for(j=0;j<PATTERNS_MAX;j++)
		pp->PMPatternPoolInLDClass[j] = FALSE;

	// LD Class List
	if(pp->PMPatternPoolLDClassList!=NULL)
	{
		int i;
		for(i=0;i<pp->PMPatternPoolLDClassListSize;i++)
		{
			if(pp->PMPatternPoolLDClassList[i]!=NULL)
				PMLDClass_free(pp->PMPatternPoolLDClassList[i]);
		}
		free(pp->PMPatternPoolLDClassList);
		pp->PMPatternPoolLDClassListSize = 0;
		pp->PMPatternPoolLDClassList = NULL;		
	}
}

void PMPatternPool_setPloidy (PMPatternPool * pp, PMCommandLine * cl)
{
	assert(pp!=NULL);
	assert(cl!=NULL);

	char Ploidy[STRING_SIZE];
	PMCommandLine_getUserFlagInput (cl, PLOIDY_INDEX, Ploidy);
	int PloidyN = -1;
	if(strlen(Ploidy)==1)
		PloidyN = Ploidy[0]-48;
	else
		PloidyN = atoi(Ploidy);

	pp->PMPatternPoolPloidy = PloidyN;
	assert(PloidyN>=1);
}

int PMPatternPool_getPloidy (PMPatternPool * pp)
{
	assert(pp!=NULL);
	
	return pp->PMPatternPoolPloidy;
}

void PMPatternPool_setLDThreshold (PMPatternPool * pp, PMCommandLine * cl)
{
	assert(pp!=NULL);
	assert(cl!=NULL);

	char LDThreshold[STRING_SIZE];
	PMCommandLine_getUserFlagInput (cl, LDTHRESHOLD_INDEX, LDThreshold);
	float LDThresholdN = -1.0;
	if(strlen(LDThreshold)==1)
		LDThresholdN = (float)(LDThreshold[0]-48);
	else
		LDThresholdN = atof(LDThreshold);

	pp->PMPatternPoolLDThreshold = LDThresholdN;
	
	assert(LDThresholdN>=0.0);
}

float PMPatternPool_getLDThreshold (PMPatternPool * pp)
{
	assert(pp!=NULL);
	
	return pp->PMPatternPoolLDThreshold;
}

uint64_t * PMPatternPool_getPatternCompact (PMPatternPool * pp, int Index)
{	
	assert(pp!=NULL);
	assert(Index>=0);
	int size = PMPatternPool_getSize(pp);
	assert(Index<=size);

	return pp->PMPatternPoolDataCompact[Index];
}

char * PMPatternPool_getPendingSNPPtr (PMPatternPool * pp)
{
	assert(pp!=NULL);
	
	return pp->PMPatternPoolPendingSNP;
}

uint64_t * PMPatternPool_getPendingSNPCompactPtr (PMPatternPool * pp)
{
	assert(pp!=NULL);
	
	return pp->PMPatternPoolPendingSNPCompact;
}

int * PMPatternPool_getCount (PMPatternPool * pp)
{
	assert(pp!=NULL);

	return pp->PMPatternPoolCount;
}

void PMPatternPool_setSampleSize (PMPatternPool * pp, int SampleSize)
{
	assert(pp!=NULL);
	assert(SampleSize>=0);

	pp->PMPatternPoolSampleSize = SampleSize;
}

int PMPatternPool_getSampleSize (PMPatternPool * pp)
{
	assert(pp!=NULL);

	return pp->PMPatternPoolSampleSize;
}

int PMPatternPool_getSize (PMPatternPool * pp)
{
	assert(pp!=NULL);

	return pp->PMPatternPoolSize;
}

void PMPatternPool_setSize (PMPatternPool * pp, int Size)
{
	assert(pp!=NULL);
	assert(Size>=0);

	pp->PMPatternPoolSize = Size;
}

int PMPatternPool_getCountByIndex (PMPatternPool * pp, int Index)
{
	assert(pp!=NULL);
	assert(Index>=0);

	return pp->PMPatternPoolCount[Index];
}

void PMPatternPool_appendPosition (PMPatternPool * pp, int index, int Position)
{
	assert(pp!=NULL);

	int size = PMPatternPool_getCountByIndex(pp, index);
	
	pp->PMPatternPoolPositions[index] = realloc(pp->PMPatternPoolPositions[index], sizeof(int*)*size);
	assert(pp->PMPatternPoolPositions[index]);

	pp->PMPatternPoolPositions[index][size-1] = Position;	
}

void PMPatternPool_appendSNPIndex (PMPatternPool * pp, int index, int SNPIndex)
{
	assert(pp!=NULL);

	int size = PMPatternPool_getCountByIndex(pp, index);
	
	pp->PMPatternPoolSNPIndeces[index] = realloc(pp->PMPatternPoolSNPIndeces[index], sizeof(int*)*size);
	assert(pp->PMPatternPoolSNPIndeces[index]);

	pp->PMPatternPoolSNPIndeces[index][size-1] = SNPIndex;	
}

void PMPatternPool_appendID (PMPatternPool * pp, int index, char * ID)
{
	assert(pp!=NULL);

	int size = PMPatternPool_getCountByIndex(pp, index);

	assert(size>=1);
	
	pp->PMPatternPoolIDs[index] = realloc(pp->PMPatternPoolIDs[index], sizeof(char**)*size);
	assert(pp->PMPatternPoolIDs[index]);

	pp->PMPatternPoolIDs[index][size-1] = (char*)malloc(sizeof(char)*STRING_SIZE);
	assert(pp->PMPatternPoolIDs[index][size-1]);

	strcpy(pp->PMPatternPoolIDs[index][size-1], ID);	
}

void PMPatternPool_setDerivedAlleleCount (PMPatternPool * pp, int index, int * DerivedAlleleCount)
{
	assert(pp!=NULL);

	if(pp->PMPatternPoolDerivedAlleleCount[index]==0)
		pp->PMPatternPoolDerivedAlleleCount[index] = DerivedAlleleCount[0];
	else
	{
		//assert((pp->PMPatternPoolDerivedAlleleCount[index]==DerivedAlleleCount[0]) || (pp->PMPatternPoolSampleSize*pp->PMPatternPoolPloidy-pp->PMPatternPoolDerivedAlleleCount[index]==DerivedAlleleCount[0]));
	}
}

int PMPatternPool_getCompactSize (PMPatternPool * pp)
{
	assert(pp!=NULL);

	return pp->PMPatternPoolCompactSize;
}

int PMPatternPool_matchSNP (PMPatternPool * pp, int Position, char * ID, int * DerivedAlleleCount, PMWindow * w, int SNP_AAindex)
{
	assert(pp!=NULL);
	assert(Position>=0);
	assert(ID!=NULL);

	//printf("Processing SNP %d\n", SNP_AAindex);

//	char * pendingSNP = PMPatternPool_getPendingSNPPtr(pp);
	uint64_t * pendingSNPCompact = PMPatternPool_getPendingSNPCompactPtr(pp); 
	int patterns = PMPatternPool_getSize(pp);
	int sampleSize = PMPatternPool_getSampleSize(pp)*PMPatternPool_getPloidy(pp)+1;
	int * patternCount = PMPatternPool_getCount (pp);

//	char * patternSNP = NULL;
	int i;

	uint64_t * patternSNPCompact = NULL;
	int CompactSize = PMPatternPool_getCompactSize(pp);
	int PatternMatch = FALSE;

	//printf("\n\t\tPPS: %d %d\n", patterns, CompactSize);
	for(i=patterns-1;i!=-1;--i)
	{
		patternSNPCompact = PMPatternPool_getPatternCompact (pp, i);
		if((dataComparison (patternSNPCompact, pendingSNPCompact, CompactSize)==0))
		{
			//printf("Match");
			PatternMatch = TRUE;
			break;		
		}

		if(dataComparisonInversion (patternSNPCompact, pendingSNPCompact, CompactSize, sampleSize-1)==0)
		{
			//printf("RMatch");
			PatternMatch = TRUE;
			break;		
		}
	}

	if(PatternMatch==TRUE)
	{
		patternCount[i]++;
		PMPatternPool_appendPosition(pp, i, Position);
		PMPatternPool_appendID(pp, i, ID);
		PMPatternPool_appendSNPIndex (pp, i, SNP_AAindex);
		PMPatternPool_setDerivedAlleleCount (pp, i, DerivedAlleleCount);
		PMWindow_appendPattern (w, i, Position);
		//printf("appending to %d\n", i);
		return FAILED;
	}
	else
	{
		patternSNPCompact = PMPatternPool_getPatternCompact (pp, patterns);
		memcpy(patternSNPCompact, pendingSNPCompact, CompactSize*8);

		patternCount[patterns]++;
		PMPatternPool_appendPosition(pp, patterns, Position);
		PMPatternPool_appendID(pp, patterns, ID);
		PMPatternPool_appendSNPIndex (pp, patterns, SNP_AAindex);
		PMPatternPool_setDerivedAlleleCount (pp, patterns, DerivedAlleleCount);
		PMWindow_appendPattern (w, patterns, Position);
		patterns++;
		PMPatternPool_setSize(pp, patterns);
		//printf("appending to %d\n", patterns-1);
		return SUCCESS;
	}
	return SUCCESS;	
}

PMWindow * PMPatternPool_getWindow (PMPatternPool * pp)
{
	assert(pp!=NULL);

	return pp->PMPatternPoolWindow;
}

int PMPatternPool_loadRegionVCF (PMPatternPool * pp, PMInput * input, FILE * fp)
{
	assert(pp!=NULL);
	assert(input!=NULL);
	assert(fp!=NULL);

	char Chromosome [STRING_SIZE];
	PMInput_getChromosome (input, Chromosome);

	int proceed = TRUE;
	int patterns = 0;
	int line = 0;

	int AlleleSize=-1;
	char AlleleVector[STRING_SIZE];
	int GenotypeLocation = -1;
	resetString(AlleleVector, STRING_SIZE);
	int SNPs = 0;
	int Position;
	char ID[STRING_SIZE];
	int DerivedAlleleCount=0;
	int SNP_AAindex = -1;
	PMWindow * w = PMPatternPool_getWindow(pp);
	assert(w!=NULL);
	PMWindow_resetLeftmostRightmostSitePosition(w);
	while(proceed==TRUE)
	{
		PMWindow_incrementSites(w);
		//printf("Line %d\n", PMWindow_getSites(w));
		int status = PMVCFParser_parseFields (fp, Chromosome, &AlleleSize, AlleleVector, &GenotypeLocation, &Position, ID);
		PMWindow_setLeftmostSitePosition(w, Position);
		PMWindow_setRightmostSitePosition(w, Position);
		if(status==SUCCESS)
		{
			PMExecutionDetails_startPatternsTimer(PMExecutionDetails_g);
			status = PMVCFParser_parseSamples (fp, pp, AlleleSize, AlleleVector, GenotypeLocation, Position, &DerivedAlleleCount);
			if(status==SUCCESS)
			{
				SNPs++;	
				//printf("");
				PMWindow_incrementSNPs(w);
				SNP_AAindex++;
				status = PMPatternPool_matchSNP(pp, Position, ID, &DerivedAlleleCount, w, SNP_AAindex); // SUCCESS returned if new pattern added
				if(status==SUCCESS)
				{
					patterns++;
					assert(patterns==PMPatternPool_getSize(pp));
					if(patterns==PATTERNS_MAX)
						proceed = FALSE;
				}	
			}
			else
			{
				PMWindow_incrementDiscarded(w);
				PMWindow_appendNoSNP (w, -1, Position);

				// TODO: error unsuccessful line, report line number
				//assert(0);
			}
			PMExecutionDetails_pausePatternsTimer(PMExecutionDetails_g);

		}
		else
		{
			PMWindow_appendNoSNP (w, -1, Position);
			PMWindow_incrementDiscarded(w);
		}
		PMVCFParser_skipLine (fp);
		
		line++;
		char c = fgetc(fp);
		if(c==EOF)
			proceed = FALSE;
		else
			ungetc(c, fp);
	}


	return SUCCESS;
}

int PMPatternPool_loadRegion (PMPatternPool * pp, PMInput * input, FILE * fp)
{
	assert(pp!=NULL);
	assert(input!=NULL);
	assert(fp!=NULL);

	int Type = PMInput_getType(input);
	PMPatternPool_reset (pp);

	if(Type==INPUT_TYPE_VCF)
	{
		int ret = PMPatternPool_loadRegionVCF (pp, input, fp);
		int SNPchecksum = 0, i, patterns = PMPatternPool_getSize(pp);
		for(i=0;i<patterns;i++)
			SNPchecksum += PMPatternPool_getCountByIndex (pp, i);
		
		assert(SNPchecksum==PMWindow_getSNPs(PMPatternPool_getWindow(pp)));


		/*{
			PMWindow * w = PMPatternPool_getWindow(pp);
			int i;
			printf("\n");
			for(i=0;i<PMWindow_getSNPs(w);i++)
			{
				int pID = w->PMWindowPatternMap[i];
				int pos = w->PMWindowPositionMap[i];
				printf("%d: [%d] %d\n", i, pID, pos);

			}
			printf("\n");
		}*/

		return ret;
	}
	else
	{
		assert(Type==INPUT_TYPE_VCF);
	}

	

	return 0;
}

int PMPatternPool_getMaxCount (PMPatternPool * pp)
{
	assert(pp!=NULL);

	int max = 0;
	int i, size = PMPatternPool_getSize(pp);
	for(i=0;i<size;i++)
		if(PMPatternPool_getCountByIndex (pp, i)>max)
			max = PMPatternPool_getCountByIndex (pp, i);

	return max;
}

int PMPatternPool_getDominantPattern (PMPatternPool * pp)
{
	assert(pp!=NULL);

	int max = 0;
	int index = 0;
	int i, size = PMPatternPool_getSize(pp);
	for(i=0;i<size;i++)
		if(PMPatternPool_getCountByIndex (pp, i)>max)
		{
			max = PMPatternPool_getCountByIndex (pp, i);
			index = i;
		}

	return index;
}

void PMPatternPool_print (PMPatternPool * pp, FILE * fp)
{
	assert(pp!=NULL);
	assert(fp!=NULL);

	fprintf(fp, "WINDOW INFORMATION\n");

	PMWindow * w = PMPatternPool_getWindow (pp);
	assert(w!=NULL);

	int LeftmostPosition = -1;
	LeftmostPosition = PMWindow_getLeftmostSitePosition(w);

	int RightmostPosition = -1;
	RightmostPosition = PMWindow_getRightmostSitePosition(w);
	fprintf(fp, "[region] %d - %d\n", LeftmostPosition, RightmostPosition);

	int Sites = -1;
	Sites = PMWindow_getSites(w);	
	fprintf(fp, "[sites] %d\n", Sites);

	int SNPs = -1;
	SNPs = PMWindow_getSNPs(w);	
	fprintf(fp, "[snps] %d\n", SNPs);

	int Discarded = -1;
	Discarded = PMWindow_getDiscarded(w);	
	fprintf(fp, "[discarded] %d\n", Discarded);

	int Patterns = -1;
	Patterns = PMPatternPool_getSize (pp);
	fprintf(fp, "[patterns] %d\n", Patterns);

	int Dominant = -1;
	Dominant = PMPatternPool_getDominantPattern(pp);
	fprintf(fp, "[dominant] %d\n", Dominant);

	int MaxCount = -1;
	MaxCount = PMPatternPool_getMaxCount (pp);
	fprintf(fp, "[maxcount] %d\n", MaxCount);

	fprintf(fp, "\nPATTERN REPORT\n");
	int i;
	for(i=0;i<Patterns;i++)
	{
		fprintf(fp, "[pattern %d][derivedcount %d][snps %d] ", i, pp->PMPatternPoolDerivedAlleleCount[i], PMPatternPool_getCountByIndex(pp, i));
		int j;
		for(j=0;j<PMPatternPool_getCountByIndex(pp, i);j++)
		{
			//fprintf(fp, "%d (%s, %d)", pp->PMPatternPoolPositions[i][j], pp->PMPatternPoolIDs[i][j], pp->PMPatternPoolSNPIndeces[i][j]);
			fprintf(fp, "%d ", pp->PMPatternPoolSNPIndeces[i][j]);
			
			if(j<PMPatternPool_getCountByIndex(pp, i)-1)
				fprintf(fp, ", ");
		}
		fprintf(fp, "\n");
	}
	
	fprintf(fp, "\n");
	

}

void PMPatternPool_computeLD (PMPatternPool * pp)
{
	assert(pp!=NULL);
	PMExecutionDetails_startLDMatrixTimer(PMExecutionDetails_g);
	int Patterns = PMPatternPool_getSize (pp);
	//double LDchecksum = 0.0;
	//uint64_t st, et;
	int i, j;
/*	// REFERENCE LD IMPLEMENTATION
	LDchecksum = 0.0;
	st = rdtsc();
	for(i=0;i<Patterns;i++)
	{
		for(j=0;j<i;j++)
		{
			pp->PMPatternPoolLDMatrix[i][j] = LD_ref(pp, i, j); 
			LDchecksum += (double)pp->PMPatternPoolLDMatrix[i][j];
		}
	}
	et = rdtsc();
	printf("Reference LD time (binary-no gaps):        %f\n", (et-st)/1.3e9);
	printf("Reference LD checksum (binary-no gaps):    %f\n", LDchecksum);
	printf("Reference LD value 0-74 (binary-no gaps):  %f\n", LD_ref(pp, 0, 74));
*/	// END OF REFERENCE LD IMPLEMENTATION

	// TEST1 LD IMPLEMENTATION
	//LDchecksum = 0.0;
	//st = rdtsc();
	printf("Patterns %d\n", Patterns);
	for(i=0;i<Patterns;i++)
	{
		for(j=0;j<Patterns;j++)
		{
			pp->PMPatternPoolLDMatrix[i][j] = LD_test1(pp, i, j); 
			//LDchecksum += (double)pp->PMPatternPoolLDMatrix[i][j];		
		}
		
	}
	//et = rdtsc();
	//printf("Test1 LD time (binary-no gaps):        %f\n", (et-st)/1.3e9);
	//printf("Test1 LD checksum (binary-no gaps):    %f\n", LDchecksum);
	//printf("Test1 LD value 0-74 (binary-no gaps):  %f\n", LD_test1(pp, 0, 74));
	// END OF REFERENCE LD IMPLEMENTATION
	PMExecutionDetails_pauseLDMatrixTimer(PMExecutionDetails_g);
}

void PMPatternPool_appendLDClass(PMPatternPool * pp, PMLDClass * ldc)
{
	assert(pp!=NULL);
	assert(ldc!=NULL);

	pp->PMPatternPoolLDClassListSize++;
	pp->PMPatternPoolLDClassList = realloc(pp->PMPatternPoolLDClassList, sizeof(PMLDClass*)*pp->PMPatternPoolLDClassListSize);
	pp->PMPatternPoolLDClassList[pp->PMPatternPoolLDClassListSize-1] = ldc;
}

double getRegionLDNoNoise (int length, int * region, PMPatternPool * pp)
{
	PMWindow * w = PMPatternPool_getWindow(pp);
	double LD = 0.0;
	int t1, t2;
	for (t1=0;t1<=length-2;t1++)
	{
		for(t2=t1+1;t2<=length-1;t2++)
		{
			int pr1 = region[t1];//w->PMWindowPatternMap[t1];
			int pr2 = region[t2];//w->PMWindowPatternMap[t2];
			if(pr1!=-1 && pr2!=-1)
				LD += (double)pp->PMPatternPoolLDMatrix[pr1][pr2];
		}
	}
	return LD;
}

double getCrossLDNoNoise (int length, int * region, int length2, int * region2, PMPatternPool * pp)
{
	PMWindow * w = PMPatternPool_getWindow(pp);
	double LD = 0.0;
	int t1, t2;
	for (t1=0;t1<=length-2;t1++)
	{
		for(t2=0;t2<=length2-2;t2++)
		{
			int pr1 = region[t1];//w->PMWindowPatternMap[t1];
			int pr2 = region2[t2];//w->PMWindowPatternMap[t2];
			if(pr1!=-1 && pr2!=-1)
				LD += (double)pp->PMPatternPoolLDMatrix[pr1][pr2];
		}
	}
	return LD;
}


void PMPatternPool_groupPatternsInLD (PMPatternPool * pp, FILE * fpclassoutput)
{
	assert(pp!=NULL);
	int i, j, Patterns = PMPatternPool_getSize (pp);
	float LDThreshold = PMPatternPool_getLDThreshold(pp);
	uint64_t st, et;
	//printf("Match LD operation for %d patterns and LD threshold %f\n", Patterns, LDThreshold);
	int counter = 0;
	st = rdtsc();
	for(i=0;i<Patterns;i++)
	{
		if(pp->PMPatternPoolInLDClass[i]==FALSE)
		{	
			PMLDClass * ldc = PMLDClass_new ();
			PMLDClass_increaseSNPSize (ldc, PMPatternPool_getCountByIndex(pp, i));
			PMLDClass_appendPattern (ldc, i);
			PMLDClass_setFirstSNPIndex (ldc, pp, i);
			PMLDClass_setLastSNPIndex (ldc, pp, i);
			counter++;
			pp->PMPatternPoolInLDClass[i] = TRUE; // TRUE
			for(j=i+1;j<Patterns;j++)
			{
				if(pp->PMPatternPoolInLDClass[j]==FALSE)
				{
					if(pp->PMPatternPoolLDMatrix[i][j]>=PMPatternPool_getLDThreshold(pp))
					{
						pp->PMPatternPoolInLDClass[j] = TRUE; // TRUE
						PMLDClass_increaseSNPSize (ldc, PMPatternPool_getCountByIndex(pp, j));
						PMLDClass_appendPattern (ldc, j);
						PMLDClass_setFirstSNPIndex (ldc, pp, j);
						PMLDClass_setLastSNPIndex (ldc, pp, j);

					}
				}
			}
			if(PMLDClass_getSNPSize(ldc)>=MININUM_SNPS_IN_LD_CLASS)
			{
				PMLDClass_createQueryStringAndRegionView(ldc, pp);
				PMPatternPool_appendLDClass(pp, ldc);
				//PMLDClass_matchQuery2Region (ldc, PMPatternPool_getWindow(pp), (void*)pp);
			}
			else
				PMLDClass_free(ldc);

		}		
	}
	et = rdtsc();
	
	for(i=0;i<pp->PMPatternPoolLDClassListSize;i++)
	{
		//printf("Class %d\n", i);
		PMLDClass_print(pp->PMPatternPoolLDClassList[i], fpclassoutput, i, pp);
		//printf("\n");
	}

	printf("\nLD clustering time:        %f (counter %d classes %d)\n", (et-st)/1.3e9, counter, pp->PMPatternPoolLDClassListSize);





/* 

	//assert(0);
	PMWindow * w = PMPatternPool_getWindow(pp);
	int windowWidth = PMWindow_getSites(w);

	
	// Init LD
	double ** wLDmat = (double**)malloc(sizeof(double*)*windowWidth);
	for(i=0;i<windowWidth;i++)
	{
		wLDmat[i] = (double*)malloc(sizeof(double)*windowWidth);
		for(j=0;j<windowWidth;j++)
		{
			/*double LD1 = 0.0;
			int t1, t2;
			for (t1=i;t1<=j-1;t1++)
			{
				for(t2=t1+1;t2<=j;t2++)
				{
					int pr1 = w->PMWindowPatternMap[t1];
					int pr2 = w->PMWindowPatternMap[t2];
					LD1 += (double)pp->PMPatternPoolLDMatrix[pr1][pr2];
				}
			}*/
/*			//printf("%d %d\n", i, j);
			int pr1 = w->PMWindowPatternMap[i];
			int pr2 = w->PMWindowPatternMap[j];

			//if(pr1>=PMPatternPool_getSize(pp) || pr2>=PMPatternPool_getSize(pp))
			//	printf("%d %d %d %d %d\n", i, j, pr1, pr2, PMPatternPool_getSize(pp));
			if(pr1!=-1&&pr2!=-1)
				wLDmat[i][j] = (double)pp->PMPatternPoolLDMatrix[pr1][pr2];
			else
				wLDmat[i][j] = 0;
		}
	
	}
	

	// Apply DP
	int lastIndex = windowWidth-1;
	for(i=0;i<=lastIndex;i++)
	{	assert(wLDmat[i][i]<=1.000001);
		wLDmat[i][i] = 0.0;
	}
	for(i=2;i<=lastIndex;i++)
	{

		for(j=i-2;j>=0;j--)
		{


			wLDmat[i][j] += wLDmat[i][j+1] + wLDmat[i-1][j] - wLDmat[i-1][j+1];
			wLDmat[j][i] = wLDmat[i][j];

		}
	}




	int samples = 1092;


	FILE * fpTest = fopen("PatternMatch.test", "w");


	int r1startprev, r1endprev, r2startprev, r2endprev;

	printf("\n");
	int R1START, R1END, R2START, R2END;
	//PMWindow * w = PMPatternPool_getWindow(pp);
	int pcounter = 0;
	{	// One more subregion pairwise attempt
		
		int classes = pp->PMPatternPoolLDClassListSize;
		//printf("\n Processing number of classes %d\n", classes);
		int mi, mj;
		double maxomega = 0.0;
		r1startprev=-1; r2endprev=-1; r2startprev=-1; r2endprev=-1;

		for(mi=0;mi<classes-1;mi++)
		{
			//ldclass1->PMLDClassUsed = 10;
			//maxomega = 0.0;	
			//if(ldclass1->PMLDClassUsed>=0)
			for(mj=mi+1;mj<classes;mj++)
			{
				PMLDClass * ldclass1 = pp->PMPatternPoolLDClassList[mi];
				PMLDClass * ldclass2 = pp->PMPatternPoolLDClassList[mj];

				int r1, regions1 = ldclass1->PMLDClassValidSubregions;
				int r2, regions2 = ldclass2->PMLDClassValidSubregions;

				//ldclass1->PMLDClassUsed--;
			
				//printf("doing pairwise classes %d(regions %d) and %d(regions %d)\n", mi,regions1, mj, regions2);
				//if(ldclass1->PMLDClassUsed>=0)
				for(r1=0;r1<regions1;r1++)
				{	//OMEGA_ACCUM = 0.0;
					
					
					for(r2=0;r2<regions2;r2++)
					{
						//printf("doing subregions %d and %d:", r1, r2);

						int r1start = ldclass1->PMLDClassValidSubregionFirstSNP[r1]-1;// + ldclass1->PMLDClassFirstSNPIndex-1;
						int r1end = ldclass1->PMLDClassValidSubregionLastSNP[r1]-1;//+ ldclass1->PMLDClassFirstSNPIndex-1;

						int r2start = ldclass2->PMLDClassValidSubregionFirstSNP[r2]-1;//+ ldclass2->PMLDClassFirstSNPIndex-1;
						int r2end = ldclass2->PMLDClassValidSubregionLastSNP[r2]-1;//+ ldclass2->PMLDClassFirstSNPIndex-1;


						if(((r1end<=r2start) || (r2end<=r1start)))
						{


						int skip=0;
						int match = 0;
						if(r1end<=r2start)
						{
							if(r1start==r1startprev)
								match++;
							
							if(r1end==r1endprev)
								match++;

							if(r2start==r2startprev)
								match++;
							
							if(r2end==r2endprev)
								match++;

							//if(match>=2)						
							//	skip=1;


							if(r1start==r1startprev && r1end==r1endprev)
								skip=1;

							if(r2start==r2startprev && r2end==r2endprev)
								skip=1;

							if(r1startprev==r1start)
								skip = 1;
				
							//if(r1endprev>=r1start)
							//	skip=1;

	
						}
						else
						{
							if(r2start==r1startprev)
								match++;
							
							if(r2end==r1endprev)
								match++;

							if(r1start==r2startprev)
								match++;
							
							if(r1end==r2endprev)
								match++;

							//if(match>=2)						
							//	skip=1;

							if(r2start==r1startprev && r2end==r1endprev)
								skip=1;

							if(r1start==r2startprev && r1end==r2endprev)
								skip=1;

							if(r1start==r2startprev)
								skip=1;
							
							//if(r1endprev>=r2start)
							//	skip=1;

			
						}		
			
						//if(skip==1)
						//	assert(0);
					

							/*	if((r1end<=r2start)) //|| (r2end<=r1start))
								{	R1START = r1start;  R1END = r1end; R2START = r2start;  R2END = r2end;
																	}
								else
								{
									R1START = r2start;  R1END = r2end; R2START = r1start;  R2END = r1end;
								}	
*/
							//printf("Pair: %d %d %d %d (%d %d)\n", R1START, R1END, R2START, R2END, mi, mj);
//							skip = 0;
							//if(r1start<=1500 || r2start<=1500 )
							//	skip =1;

/*							if((((r1end<=r2start)&&(r2end-r1start<=1000)) || ((r2end<=r1start)&&(r1end-r2start<=1000)))&&skip==0)
							{
							pcounter++;
							float preLD1 = ldclass1->PMLDClassValidSubregionTotalLD[r1];
							float preLD2 = ldclass2->PMLDClassValidSubregionTotalLD[r2];

							int t1, t2;
							/**/ // No noise
/*							double LD1 = getRegionLDNoNoise (ldclass1->PMLDClassSubregionList[r1]->PMLDSubregionLength, ldclass1->PMLDClassSubregionList[r1]->PMLDSubregionPatternIDWindow, pp);//wLDmat[r1start][r1end];
//printf("LD1: %f previous %f\n", LD1, wLDmat[r1start][r1end]);
							double LD2 = getRegionLDNoNoise (ldclass2->PMLDClassSubregionList[r2]->PMLDSubregionLength, ldclass2->PMLDClassSubregionList[r2]->PMLDSubregionPatternIDWindow, pp); //wLDmat[r2start][r2end];
							double LD3 = 0.0;
							//if(r1end<=r2start)
							//	LD3 = wLDmat[r1start][r2end] - wLDmat[r1start][r2start-1] - wLDmat[r1end+1][r2end] + wLDmat[r1end+1][r2start-1];

							//if(r2end<=r1start)
							//	LD3 = wLDmat[r2start][r1end] - wLDmat[r2start][r1start-1] - wLDmat[r2end+1][r1end] + wLDmat[r2end+1][r1start-1];

							LD3 = getCrossLDNoNoise (ldclass1->PMLDClassSubregionList[r1]->PMLDSubregionLength, ldclass1->PMLDClassSubregionList[r1]->PMLDSubregionPatternIDWindow, ldclass2->PMLDClassSubregionList[r2]->PMLDSubregionLength, ldclass2->PMLDClassSubregionList[r2]->PMLDSubregionPatternIDWindow, pp);
							/**/



							/** // With noise
							double LD1 =wLDmat[r1start][r1end];
							double LD2 = wLDmat[r2start][r2end];
							double LD3 = 0.0;
							if(r1end<=r2start)
								LD3 = wLDmat[r1start][r2end] - wLDmat[r1start][r2start-1] - wLDmat[r1end+1][r2end] + wLDmat[r1end+1][r2start-1];

							if(r2end<=r1start)
								LD3 = wLDmat[r2start][r1end] - wLDmat[r2start][r1start-1] - wLDmat[r2end+1][r1end] + wLDmat[r2end+1][r1start-1];
							/**/



/*							int snps1 = r1end - r1start + 1;
							double snps1sel2 = ((double)snps1*(snps1-1))/2.0;
							int snps2 = r2end - r2start + 1;
							double snps2sel2 = ((double)snps2*(snps2-1))/2.0;
							double snps1snps2 = ((double)snps1)*snps2;
							//printf("jj0: %d %d\n", r1startprev, R1START);
							double omega = ((LD1+LD2)/(snps1sel2+snps2sel2)) / (((LD3+(1.0/(samples*2.0)))/snps1snps2));
							if(omega>=maxomega)
							{
								maxomega = omega;
								if((r1end<=r2start)) //|| (r2end<=r1start))
								{	R1START = r1start;  R1END = r1end; R2START = r2start;  R2END = r2end;
									r1startprev = R1START;
									r1endprev = R1END;
									//r2startprev = R2START;
									//r2endprev = R2END;
									printf(" Best match region %d %d - %d %d Max Score <%f> %f %f (%f) %d %d %d %d [%d %d]:", r1start, r1end, r2start, r2end, omega, LD1, LD2, LD3, snps1, snps2, r1startprev, r1endprev,mi, mj);									
								}
								else
								{
									R1START = r2start;  R1END = r2end; R2START = r1start;  R2END = r1end;
									r1startprev = R1START;
									r1endprev = R1END;
									//r2startprev = R2START;
									//r2endprev = R2END;
									printf(" Best match region %d %d - %d %d Max Score <%f> %f %f (%f) %d %d %d %d [%d %d]:", r2start, r2end, r1start, r1end, omega, LD1, LD2, LD3, snps1, snps2, r1startprev, r1endprev, mi, mj);									
								}

								//ldclass1->PMLDClassUsed--;
								//ldclass2->PMLDClassUsed--;

							//		printf("jj: %d %d\n", r1startprev, R1START);

									r1startprev = R1START;
									r1endprev = R1END;
									r2startprev = R2START;
									r2endprev = R2END;

							//		printf("jj2: %d %d\n", r1startprev, R1START);

								


		
		{// Final refinement process 2
			int MINSNPS = 5;
			//printf("\nFINAL REFINEMENT  2 Doing region %d %d - %d %d \n", R1START, R1END, R2START, R2END);
			double maxomega = 0.0;
			int maxiloc;int maxileft; int maxiright;
			int ri;
			for(ri=R1END+MINSNPS;ri<R2START-MINSNPS;ri++)
			//ri = (R1END+MINSNPS) + ((R2START-MINSNPS) - (R1END+MINSNPS))/2.0;
			{
				int li;
				for(li=R1START;li<=ri-MINSNPS;li++)
				{
					double lLD = wLDmat[li][ri];
					
					int rri;
					for(rri=ri+1+MINSNPS;rri<=R2END;rri++)
					{
						double rLD= wLDmat[ri+1][rri];
						double cLD =  wLDmat[li][rri] - wLDmat[li][ri+1+MINSNPS-1] - wLDmat[ri-MINSNPS+1][rri] + wLDmat[ri-MINSNPS+1][ri+1+MINSNPS-1];
						int snps1 = ri - li + 1;
						double snps1sel2 = ((double)snps1*(snps1-1))/2.0;
						int snps2 = rri - (ri+1) + 1;
						double snps2sel2 = ((double)snps2*(snps2-1))/2.0;
						double snps1snps2 = ((double)snps1)*snps2;
						
						double omegaF = ((lLD+rLD)/(snps1sel2+snps2sel2)) / ((cLD+(1.0/(samples*2.0)))/snps1snps2);
						//fprintf(fpTest, "%f\t%f\n", (float)ri, omegaF);
						if(omegaF>=maxomega)
						{	maxomega = omegaF;
							//printf("%d omega %f %f %f <%f> MAX\n", ri, lLD, rLD, cLD, omegaF);
							//printf("%d:<%f>, ", ri, omegaF);
							maxiloc = ri;
							maxileft = li;
							maxiright = rri;
						}
					
					}
				}
				//if(ri>1500)
				//	fprintf(fpTest, "%f\t%f\n", (float)maxiloc, maxomega);
			}
			printf(" MAX %f at %d (%d %d)", maxomega, maxiloc, maxileft, maxiright);
			fprintf(fpTest, "%f\t%f\n", (float)maxiloc, maxomega);
		}

							printf("\n");

	
							}
							}
						}
					}
				}
			}
		}
	}
	
		printf("Pairs: %d\n", pcounter);		
		
		int windowSize = PMWindow_getSNPs(w);
		for(i=0;i<windowSize;i++)
		{
			if(w->PMWindowPositionMap[i]>3345)//9468)//71345376)//345376.00)
			{
				REFERENCE_OP_INDEX = i;
				break;
			}	
		}
		fprintf(stdout, "OmegaPlus indicated index: %d\n", REFERENCE_OP_INDEX);
	
		int REFERENCE_SD_INDEX = 0;
		for(i=0;i<windowSize;i++)
		{
			if(w->PMWindowPositionMap[i]>3313)
			{
				REFERENCE_SD_INDEX = i;
				break;
			}	
		}
		fprintf(stdout, "Left index: %d\n", REFERENCE_SD_INDEX);

		int REFERENCE_SD_INDEX2 = 0;
		for(i=0;i<windowSize;i++)
		{
			if(w->PMWindowPositionMap[i]>3403)
			{
				REFERENCE_SD_INDEX2 = i;
				break;
			}	
		}
		fprintf(stdout, "Right index: %d\n", REFERENCE_SD_INDEX2);
		



	for(i=0;i<windowWidth;i++)
	{
		free(wLDmat[i]);
	}
	free(wLDmat);

	fclose(fpTest);
*/

}





void PMPatternPool_groupPatternsInLDWithLDNoise (PMPatternPool * pp)
{
	assert(pp!=NULL);
	int i, j, Patterns = PMPatternPool_getSize (pp);
	float LDThreshold = PMPatternPool_getLDThreshold(pp);
	uint64_t st, et;
	//printf("Match LD operation for %d patterns and LD threshold %f\n", Patterns, LDThreshold);
	int counter = 0;
	st = rdtsc();
	for(i=0;i<Patterns;i++)
	{
		if(pp->PMPatternPoolInLDClass[i]==FALSE)
		{	
			PMLDClass * ldc = PMLDClass_new ();
			PMLDClass_increaseSNPSize (ldc, PMPatternPool_getCountByIndex(pp, i));
			PMLDClass_appendPattern (ldc, i);
			PMLDClass_setFirstSNPIndex (ldc, pp, i);
			PMLDClass_setLastSNPIndex (ldc, pp, i);
			counter++;
			pp->PMPatternPoolInLDClass[i] = FALSE; // TRUE
			for(j=i+1;j<Patterns;j++)
			{
				if(pp->PMPatternPoolInLDClass[j]==FALSE)
				{
					if(pp->PMPatternPoolLDMatrix[i][j]>=PMPatternPool_getLDThreshold(pp))
					{
						pp->PMPatternPoolInLDClass[j] = FALSE; // TRUE
						PMLDClass_increaseSNPSize (ldc, PMPatternPool_getCountByIndex(pp, j));
						PMLDClass_appendPattern (ldc, j);
						PMLDClass_setFirstSNPIndex (ldc, pp, j);
						PMLDClass_setLastSNPIndex (ldc, pp, j);

					}
				}
			}
			if(PMLDClass_getSNPSize(ldc)>=MININUM_SNPS_IN_LD_CLASS)
			{
				PMLDClass_createQueryStringAndRegionView(ldc, pp);
				PMPatternPool_appendLDClass(pp, ldc);
				PMLDClass_matchQuery2Region (ldc, PMPatternPool_getWindow(pp), (void*)pp);
			}
			else
				PMLDClass_free(ldc);

		}		
	}
	et = rdtsc();
	
	for(i=0;i<pp->PMPatternPoolLDClassListSize;i++)
	{
		//printf("Class %d\n", i);
		PMLDClass_print(pp->PMPatternPoolLDClassList[i], stdout, i, pp);
		//printf("\n");
	}

	printf("\nLD clustering time:        %f (counter %d classes %d)\n", (et-st)/1.3e9, counter, pp->PMPatternPoolLDClassListSize);








	PMWindow * w = PMPatternPool_getWindow(pp);
	int windowWidth = PMWindow_getSites(w);

	
	// Init LD
	double ** wLDmat = (double**)malloc(sizeof(double*)*windowWidth);
	for(i=0;i<windowWidth;i++)
	{
		wLDmat[i] = (double*)malloc(sizeof(double)*windowWidth);
		for(j=0;j<windowWidth;j++)
		{
			/*double LD1 = 0.0;
			int t1, t2;
			for (t1=i;t1<=j-1;t1++)
			{
				for(t2=t1+1;t2<=j;t2++)
				{
					int pr1 = w->PMWindowPatternMap[t1];
					int pr2 = w->PMWindowPatternMap[t2];
					LD1 += (double)pp->PMPatternPoolLDMatrix[pr1][pr2];
				}
			}*/
			//printf("%d %d\n", i, j);
			int pr1 = w->PMWindowPatternMap[i];
			int pr2 = w->PMWindowPatternMap[j];

			//if(pr1>=PMPatternPool_getSize(pp) || pr2>=PMPatternPool_getSize(pp))
			//	printf("%d %d %d %d %d\n", i, j, pr1, pr2, PMPatternPool_getSize(pp));
			if(pr1!=-1&&pr2!=-1)
				wLDmat[i][j] = (double)pp->PMPatternPoolLDMatrix[pr1][pr2];
			else
				wLDmat[i][j] = 0;
		}
	
	}
	

	// Apply DP
	int lastIndex = windowWidth-1;
	for(i=0;i<=lastIndex;i++)
	{	assert(wLDmat[i][i]<=1.000001);
		wLDmat[i][i] = 0.0;
	}
	for(i=2;i<=lastIndex;i++)
	{

		for(j=i-2;j>=0;j--)
		{


			wLDmat[i][j] += wLDmat[i][j+1] + wLDmat[i-1][j] - wLDmat[i-1][j+1];
			wLDmat[j][i] = wLDmat[i][j];

		}
	}












	printf("\n");
	int R1START, R1END, R2START, R2END;
	//PMWindow * w = PMPatternPool_getWindow(pp);
	int pcounter = 0;
	{	// One more subregion pairwise attempt
		
		int classes = pp->PMPatternPoolLDClassListSize;
		//printf("\n Processing number of classes %d\n", classes);
		int mi, mj;
		double maxomega = 0.0;
		for(mi=0;mi<classes-1;mi++)
		{
			for(mj=mi+1;mj<classes;mj++)
			{
				
				PMLDClass * ldclass1 = pp->PMPatternPoolLDClassList[mi];
				PMLDClass * ldclass2 = pp->PMPatternPoolLDClassList[mj];

				int r1, regions1 = ldclass1->PMLDClassValidSubregions;
				int r2, regions2 = ldclass2->PMLDClassValidSubregions;
			
				//printf("doing pairwise classes %d(regions %d) and %d(regions %d)\n", mi,regions1, mj, regions2);
				for(r1=0;r1<regions1;r1++)
				{	//OMEGA_ACCUM = 0.0;
					for(r2=0;r2<regions2;r2++)
					{
						//printf("doing subregions %d and %d:", r1, r2);

						int r1start = ldclass1->PMLDClassValidSubregionFirstSNP[r1]-1;// + ldclass1->PMLDClassFirstSNPIndex-1;
						int r1end = ldclass1->PMLDClassValidSubregionLastSNP[r1]-1;//+ ldclass1->PMLDClassFirstSNPIndex-1;

						int r2start = ldclass2->PMLDClassValidSubregionFirstSNP[r2]-1;//+ ldclass2->PMLDClassFirstSNPIndex-1;
						int r2end = ldclass2->PMLDClassValidSubregionLastSNP[r2]-1;//+ ldclass2->PMLDClassFirstSNPIndex-1;						
						

					
						if((r1end<=r2start) || (r2end<=r1start))
						{

							if(((r1end<=r2start)&&(r2end-r1start<=500)) || ((r2end<=r1start)&&(r1end-r2start<=500)))
							{
							pcounter++;
							float preLD1 = ldclass1->PMLDClassValidSubregionTotalLD[r1];
							float preLD2 = ldclass2->PMLDClassValidSubregionTotalLD[r2];

							int t1, t2;
							double LD1 = wLDmat[r1start][r1end];
							double LD2 = wLDmat[r2start][r2end];
							double LD3 = 0.0;
							if(r1end<=r2start)
								LD3 = wLDmat[r1start][r2end] - wLDmat[r1start][r2start-1] - wLDmat[r1end+1][r2end] + wLDmat[r1end+1][r2start-1];

							if(r2end<=r1start)
								LD3 = wLDmat[r2start][r1end] - wLDmat[r2start][r1start-1] - wLDmat[r2end+1][r1end] + wLDmat[r2end+1][r1start-1];

							int snps1 = r1end - r1start + 1;
							double snps1sel2 = ((double)snps1*(snps1-1))/2.0;
							int snps2 = r2end - r2start + 1;
							double snps2sel2 = ((double)snps2*(snps2-1))/2.0;
							double snps1snps2 = ((double)snps1)*snps2;
							double omega = ((LD1+LD2)/(snps1sel2+snps2sel2)) / (LD3/snps1snps2);
							if(omega>=maxomega)
							{
								maxomega = omega;
								if((r1end<=r2start)) //|| (r2end<=r1start))
								{
									printf(" Best match region %d %d - %d %d Max Score <%f> %f %f %d %d:", r1start, r1end, r2start, r2end, omega, LD1, LD2, snps1, snps2);									R1START = r1start;  R1END = r1end; R2START = r2start;  R2END = r2end;
								}
								else
								{
									printf(" Best match region %d %d - %d %d Max Score <%f> %f %f %d %d:", r2start, r2end, r1start, r1end, omega, LD1, LD2, snps1, snps2);									R1START = r2start;  R1END = r2end; R2START = r1start;  R2END = r1end;
								}



		{// Final refinement process 2
			int MINSNPS = 3;
			//printf("\nFINAL REFINEMENT  2 Doing region %d %d - %d %d \n", R1START, R1END, R2START, R2END);
			double maxomega = 0.0;
			int maxiloc;
			int ri;
			for(ri=R1END+MINSNPS;ri<R2START-MINSNPS;ri++)
			{
				int li;
				for(li=R1START;li<=ri-MINSNPS;li++)
				{
					double lLD = wLDmat[li][ri];
					
					int rri;
					for(rri=ri+1+MINSNPS;rri<=R2END;rri++)
					{
						double rLD= wLDmat[ri+1][rri];
						double cLD =  wLDmat[li][rri] - wLDmat[li][ri+1+MINSNPS-1] - wLDmat[ri-MINSNPS+1][rri] + wLDmat[ri-MINSNPS+1][ri+1+MINSNPS-1];
						int snps1 = ri - li + 1;
						double snps1sel2 = ((double)snps1*(snps1-1))/2.0;
						int snps2 = rri - (ri+1) + 1;
						double snps2sel2 = ((double)snps2*(snps2-1))/2.0;
						double snps1snps2 = ((double)snps1)*snps2;
						double omegaF = ((lLD+rLD)/(snps1sel2+snps2sel2)) / ((cLD+(1.0/(50.0*2.0)))/snps1snps2);
						if(omegaF>=maxomega)
						{	maxomega = omegaF;
							//printf("%d omega %f %f %f <%f> MAX\n", ri, lLD, rLD, cLD, omegaF);
							//printf("%d:<%f>, ", ri, omegaF);
							maxiloc = ri;
						}
					
					}
				}
			}
			printf(" MAX %f at %d\n", maxomega, maxiloc);
		}

							printf("\n");

	
							}
							}
						}
					}
				}
			}
		}
	}
	/*
		printf("Pairs: %d\n", pcounter);		
		
		int windowSize = PMWindow_getSNPs(w);
		for(i=0;i<windowSize;i++)
		{
			if(w->PMWindowPositionMap[i]>1377)//9468)//71345376)//345376.00)
			{
				REFERENCE_OP_INDEX = i;
				break;
			}	
		}
		fprintf(stdout, "OmegaPlus indicated index: %d\n", REFERENCE_OP_INDEX);
	
		int REFERENCE_SD_INDEX = 0;
		for(i=0;i<windowSize;i++)
		{
			if(w->PMWindowPositionMap[i]>71344426)
			{
				REFERENCE_SD_INDEX = i;
				break;
			}	
		}
		fprintf(stdout, "Left index: %d\n", REFERENCE_SD_INDEX);

		int REFERENCE_SD_INDEX2 = 0;
		for(i=0;i<windowSize;i++)
		{
			if(w->PMWindowPositionMap[i]>71346336)
			{
				REFERENCE_SD_INDEX2 = i;
				break;
			}	
		}
		fprintf(stdout, "Right index: %d\n", REFERENCE_SD_INDEX2);
	*/	



	for(i=0;i<windowWidth;i++)
	{
		free(wLDmat[i]);
	}
	free(wLDmat);



}





void PMPatternPool_groupPatternsInLD2 (PMPatternPool * pp)
{
	assert(pp!=NULL);
	int i, j, Patterns = PMPatternPool_getSize (pp);
	float LDThreshold = PMPatternPool_getLDThreshold(pp);
	uint64_t st, et;
	//printf("Match LD operation for %d patterns and LD threshold %f\n", Patterns, LDThreshold);
	int counter = 0;
	st = rdtsc();
	for(i=0;i<Patterns;i++)
	{
		if(pp->PMPatternPoolInLDClass[i]==FALSE)
		{	
			PMLDClass * ldc = PMLDClass_new ();
			PMLDClass_increaseSNPSize (ldc, PMPatternPool_getCountByIndex(pp, i));
			PMLDClass_appendPattern (ldc, i);
			PMLDClass_setFirstSNPIndex (ldc, pp, i);
			PMLDClass_setLastSNPIndex (ldc, pp, i);
			//printf("%d SNP %d in LD with: ", counter, i);
			counter++;
			pp->PMPatternPoolInLDClass[i] = TRUE; // TRUE
			for(j=i+1;j<Patterns;j++)
			{
				if(pp->PMPatternPoolInLDClass[j]==FALSE)
				{
					//pp->PMPatternPoolLDMatrix[i][j] = LD_test1(pp, i, j); 
					//LDchecksum += (double)pp->PMPatternPoolLDMatrix[i][j];
					//LDchecksum += (double)LD_test1(pp, i, j);

					if(pp->PMPatternPoolLDMatrix[i][j]>=PMPatternPool_getLDThreshold(pp))
					{
						//printf(" %d, ", j);
						pp->PMPatternPoolInLDClass[j] = TRUE; // TRUE
						PMLDClass_increaseSNPSize (ldc, PMPatternPool_getCountByIndex(pp, j));
						PMLDClass_appendPattern (ldc, j);
						PMLDClass_setFirstSNPIndex (ldc, pp, j);
						PMLDClass_setLastSNPIndex (ldc, pp, j);

					}
				}
			}
			//printf(" (total snps %d)\n", PMLDClass_getSNPSize (ldc));
			if(PMLDClass_getSNPSize(ldc)>=MININUM_SNPS_IN_LD_CLASS)
			{
				PMLDClass_createQueryStringAndRegionView(ldc, pp);
				PMPatternPool_appendLDClass(pp, ldc);
				//printf("test: %s %d\nmatch: %s %d\n", ldc->PMLDClassRegionView, strlen(ldc->PMLDClassRegionView), ldc->PMLDClassQueryString, strlen(ldc->PMLDClassQueryString));
				/*int ll, myindex = 0;
				printf("\nRegion\n");
				for(ll=0;ll<strlen(ldc->PMLDClassRegionView);ll++)
				{
					if(ll%10==0)
						printf("\n%d:", myindex++);

					printf("%c ",ldc->PMLDClassRegionView[ll]);

				}*/
				// Leaks here
				PMLDClass_matchQuery2Region (ldc, PMPatternPool_getWindow(pp), (void*)pp);
				//PMLDClass_processSubregions (ldc, PMPatternPool_getWindow(pp), (void*)pp);
				//printf("Matches:\n");
				//for(ll=0;ll<ldc->PMLDClassMatchListSize;ll++)
				{
				//	printf("MM%d: %d(%d) %d %d\n", ll, ldc->extra[ll],ldc->PMLDClassMatchListScore[ll],  ldc->PMLDClassMatchListCoordinateI[ll], ldc->PMLDClassMatchListCoordinateJ[ll]);
				}
				//stringMatch (ldc->PMLDClassQueryString, ldc->PMLDClassRegionView);
			}
			else
				PMLDClass_free(ldc);

		}		
	}
	et = rdtsc();
	//printf("LD classification time:        %f (counter %d classes %d)\n", (et-st)/1.3e9, counter, pp->PMPatternPoolLDClassListSize);


	PMWindow * w = PMPatternPool_getWindow(pp);
	int windowWidth = PMWindow_getSites(w);

	
	// Init LD
	double ** wLDmat = (double**)malloc(sizeof(double*)*windowWidth);
	for(i=0;i<windowWidth;i++)
	{
		wLDmat[i] = (double*)malloc(sizeof(double)*windowWidth);
		for(j=0;j<windowWidth;j++)
		{
/*			/*double LD1 = 0.0;
			int t1, t2;
			for (t1=i;t1<=j-1;t1++)
			{
				for(t2=t1+1;t2<=j;t2++)
				{
					int pr1 = w->PMWindowPatternMap[t1];
					int pr2 = w->PMWindowPatternMap[t2];
					LD1 += (double)pp->PMPatternPoolLDMatrix[pr1][pr2];
				}
			}*/
			//printf("%d %d\n", i, j);
			int pr1 = w->PMWindowPatternMap[i];
			int pr2 = w->PMWindowPatternMap[j];
			
			wLDmat[i][j] = (double)pp->PMPatternPoolLDMatrix[pr1][pr2];
		}
	}


	// Apply DP
	int lastIndex = windowWidth-1;
	for(i=0;i<=lastIndex;i++)
	{	assert(wLDmat[i][i]<=1.000001);
		wLDmat[i][i] = 0.0;
	}
	for(i=2;i<=lastIndex;i++)
	{

		for(j=i-2;j>=0;j--)
		{


			wLDmat[i][j] += wLDmat[i][j+1] + wLDmat[i-1][j] - wLDmat[i-1][j+1];
			wLDmat[j][i] = wLDmat[i][j];

		}
	}


	
	
/*	
	for(i=0;i<pp->PMPatternPoolLDClassListSize;i++)
	{
		//printf("Class %d\n", i);
		//PMLDClass_print(pp->PMPatternPoolLDClassList[i], stdout, i);
		//printf("\n");
	}

*/
	printf("\n");
	int R1START, R1END, R2START, R2END;
	//PMWindow * w = PMPatternPool_getWindow(pp);
	int pcounter = 0;
	{	// One more subregion pairwise attempt
		
		int classes = pp->PMPatternPoolLDClassListSize;
		//printf("\n Processing number of classes %d\n", classes);
		int mi, mj;
		double maxomega = 0.0;
		for(mi=0;mi<classes-1;mi++)
		{
			for(mj=mi+1;mj<classes;mj++)
			{
				
				PMLDClass * ldclass1 = pp->PMPatternPoolLDClassList[mi];
				PMLDClass * ldclass2 = pp->PMPatternPoolLDClassList[mj];

				int r1, regions1 = ldclass1->PMLDClassValidSubregions;
				int r2, regions2 = ldclass2->PMLDClassValidSubregions;
			
				//printf("doing pairwise classes %d(regions %d) and %d(regions %d)\n", mi,regions1, mj, regions2);
				for(r1=0;r1<regions1;r1++)
				{	//OMEGA_ACCUM = 0.0;
					for(r2=0;r2<regions2;r2++)
					{
						//printf("doing subregions %d and %d:", r1, r2);

						int r1start = ldclass1->PMLDClassValidSubregionFirstSNP[r1]-1;// + ldclass1->PMLDClassFirstSNPIndex-1;
						int r1end = ldclass1->PMLDClassValidSubregionLastSNP[r1]-1;//+ ldclass1->PMLDClassFirstSNPIndex-1;

						int r2start = ldclass2->PMLDClassValidSubregionFirstSNP[r2]-1;//+ ldclass2->PMLDClassFirstSNPIndex-1;
						int r2end = ldclass2->PMLDClassValidSubregionLastSNP[r2]-1;//+ ldclass2->PMLDClassFirstSNPIndex-1;						
						

					
						if((r1end<=r2start) || (r2end<=r1start))
						{

							if(((r1end<=r2start)&&(r2end-r1start<=500)) || ((r2end<=r1start)&&(r1end-r2start<=500)))
							{

						

							pcounter++;
							//printf("\n%d %d <-> %d %d", r1start, r1end, r2start, r2end);

							// precomputed LD checks
							float preLD1 = ldclass1->PMLDClassValidSubregionTotalLD[r1];
							float preLD2 = ldclass2->PMLDClassValidSubregionTotalLD[r2];
							//printf("%f %f :", preLD1, preLD2);

							int t1, t2;
									/*double LD1 = 0.0;
									for (t1=r1start;t1<=r1end-1;t1++)
									{
										for(t2=t1+1;t2<=r1end;t2++)
										{
											int pr1 = w->PMWindowPatternMap[t1];
											int pr2 = w->PMWindowPatternMap[t2];
											LD1 += (double)pp->PMPatternPoolLDMatrix[pr1][pr2];
										}
									}*/
									//assert(LD1==wLDmat[r1start][r1end]);
							double LD1 = wLDmat[r1start][r1end];
									//printf(" %f ", LD1);

									/*double LD2 = 0.0;
									for (t1=r2start;t1<=r2end-1;t1++)
									{
										for(t2=t1+1;t2<=r2end;t2++)
										{
											int pr1 = w->PMWindowPatternMap[t1];
											int pr2 = w->PMWindowPatternMap[t2];
											LD2 += (double)pp->PMPatternPoolLDMatrix[pr1][pr2];
										}
									}*/
									//printf(" %f ", LD2);
							double LD2 = wLDmat[r2start][r2end];

									/*double LD3 = 0.0;
									for (t1=r1start;t1<=r1end;t1++)
									{
										for(t2=r2start;t2<=r2end;t2++)
										{
											int pr1 = w->PMWindowPatternMap[t1];
											int pr2 = w->PMWindowPatternMap[t2];
											LD3 += (double)pp->PMPatternPoolLDMatrix[pr1][pr2];
										}
									}*/
							double LD3 = 0.0;
							if(r1end<=r2start)
								LD3 = wLDmat[r1start][r2end] - wLDmat[r1start][r2start-1] - wLDmat[r1end+1][r2end] + wLDmat[r1end+1][r2start-1];

							if(r2end<=r1start)
								LD3 = wLDmat[r2start][r1end] - wLDmat[r2start][r1start-1] - wLDmat[r2end+1][r1end] + wLDmat[r2end+1][r1start-1];

									//printf(" %f : ", LD3);
							int snps1 = r1end - r1start + 1;
							double snps1sel2 = ((double)snps1*(snps1-1))/2.0;
									//printf("%f ", snps1sel2);
							int snps2 = r2end - r2start + 1;
							double snps2sel2 = ((double)snps2*(snps2-1))/2.0;
									//printf("%f ", snps2sel2);
							double snps1snps2 = ((double)snps1)*snps2;
									//printf("%f ", snps1snps2);
							double omega = ((LD1+LD2)/(snps1sel2+snps2sel2)) / (LD3/snps1snps2);
							if(omega>=maxomega)
							{
								//printf("\n\ndoing pairwise classes %d(regions %d) and %d(regions %d)\n", mi,regions1, mj, regions2);
								maxomega = omega;
								//printf( "\t\tomega: %f MAX\n", omega);
								if((r1end<=r2start)) //|| (r2end<=r1start))
								{
									printf(" Best match region %d %d - %d %d Max Score <%f> %f %f %d %d:", r1start, r1end, r2start, r2end, omega, LD1, LD2, snps1, snps2);
									R1START = r1start;  R1END = r1end; R2START = r2start;  R2END = r2end;
								}
								else
								{
									printf(" Best match region %d %d - %d %d Max Score <%f> %f %f %d %d:", r2start, r2end, r1start, r1end, omega, LD1, LD2, snps1, snps2);
									R1START = r2start;  R1END = r2end; R2START = r1start;  R2END = r1end;
								}



{ // Final refinement process 2
			int MINSNPS = 3;
			//printf("\nFINAL REFINEMENT  2 Doing region %d %d - %d %d \n", R1START, R1END, R2START, R2END);
			double maxomega = 0.0;
			int maxiloc;
			int ri;
			for(ri=R1END+MINSNPS;ri<R2START-MINSNPS;ri++)
			{
				int li;
				for(li=R1START;li<=ri-MINSNPS;li++)
				{
					double lLD = wLDmat[li][ri];
					
					int rri;
					for(rri=ri+1+MINSNPS;rri<=R2END;rri++)
					{
						double rLD= wLDmat[ri+1][rri];
						double cLD =  wLDmat[li][rri] - wLDmat[li][ri+1+MINSNPS-1] - wLDmat[ri-MINSNPS+1][rri] + wLDmat[ri-MINSNPS+1][ri+1+MINSNPS-1];
						int snps1 = ri - li + 1;
						double snps1sel2 = ((double)snps1*(snps1-1))/2.0;
						int snps2 = rri - (ri+1) + 1;
						double snps2sel2 = ((double)snps2*(snps2-1))/2.0;
						double snps1snps2 = ((double)snps1)*snps2;
						double omegaF = ((lLD+rLD)/(snps1sel2+snps2sel2)) / ((cLD+(1.0/(50.0*2.0)))/snps1snps2);
						if(omegaF>=maxomega)
						{	maxomega = omegaF;
							//printf("%d omega %f %f %f <%f> MAX\n", ri, lLD, rLD, cLD, omegaF);
							//printf("%d:<%f>, ", ri, omegaF);
							maxiloc = ri;
						}
					
					}
				}
			}
			printf(" MAX %f at %d\n", maxomega, maxiloc);
		}


	
							}
							else
							{
								//printf( "\t\tomega: %f \n", omega);
				
							}
							}
						}
						//else {//printf("\n");
						//	}
					}
				}

			}
		}

		printf("Pairs: %d\n", pcounter);		
		
		int windowSize = PMWindow_getSNPs(w);
		for(i=0;i<windowSize;i++)
		{
			if(w->PMWindowPositionMap[i]>1377)//9468)//71345376)//345376.00)
			{
				REFERENCE_OP_INDEX = i;
				break;
			}	
		}
		fprintf(stdout, "OmegaPlus indicated index: %d\n", REFERENCE_OP_INDEX);
	
		int REFERENCE_SD_INDEX = 0;
		for(i=0;i<windowSize;i++)
		{
			if(w->PMWindowPositionMap[i]>71344426)
			{
				REFERENCE_SD_INDEX = i;
				break;
			}	
		}
		fprintf(stdout, "Left index: %d\n", REFERENCE_SD_INDEX);

		int REFERENCE_SD_INDEX2 = 0;
		for(i=0;i<windowSize;i++)
		{
			if(w->PMWindowPositionMap[i]>71346336)
			{
				REFERENCE_SD_INDEX2 = i;
				break;
			}	
		}
		fprintf(stdout, "Right index: %d\n", REFERENCE_SD_INDEX2);
		

	}

	/*	{ // Final refinement process 1
		
			printf(" FINAL REFINEMENT Doing region %d %d - %d %d \n", R1START, R1END, R2START, R2END);
			double maxomega = 0.0;
			int ri;
			for(ri=R1END;ri<R2START;ri++)
			{
				//printf("loc %d  ", ri);

				// left
				double lLD = 0.0;
				int li1, li2;
				for(li1=R1START;li1<=ri-1;li1++)
				{
					for(li2=li1+1;li2<=ri;li2++)
					{
						int pr1 = w->PMWindowPatternMap[li1];
						int pr2 = w->PMWindowPatternMap[li2];
						lLD += (double)pp->PMPatternPoolLDMatrix[pr1][pr2];
					}
				}
				//printf("lLD %f", lLD);

				// right
				double rLD = 0.0;
				int ri1, ri2;
				for(ri1=ri+1;ri1<=R2END-1;ri1++)
				{
					for(ri2=ri1+1;ri2<=R2END;ri2++)
					{
						int pr1 = w->PMWindowPatternMap[ri1];
						int pr2 = w->PMWindowPatternMap[ri2];
						rLD += (double)pp->PMPatternPoolLDMatrix[pr1][pr2];
					}
				}
				//printf("rLD %f", rLD);


				// cross
				double cLD = 0.0;
				int ci1, ci2;
				for(ci1=R1START;ci1<=ri;ci1++)
				{
					for(ci2=ri+1;ci2<=R2END;ci2++)
					{
						int pr1 = w->PMWindowPatternMap[ci1];
						int pr2 = w->PMWindowPatternMap[ci2];
						cLD += (double)pp->PMPatternPoolLDMatrix[pr1][pr2];
					}
				}
				//printf("cLD %f", cLD);


				int snps1 = ri - R1START + 1;
				double snps1sel2 = ((double)snps1*(snps1-1))/2.0;
				//printf("%f ", snps1sel2);
				int snps2 = R2END - (ri+1) + 1;
				double snps2sel2 = ((double)snps2*(snps2-1))/2.0;
				//printf("%f ", snps2sel2);
				double snps1snps2 = ((double)snps1)*snps2;
				//printf("%f ", snps1snps2);

				
				double omegaF = ((lLD+rLD)/(snps1sel2+snps2sel2)) / (cLD/snps1snps2);
				//printf("omega %f", omegaF);
				if(omegaF>=maxomega)
				{	maxomega = omegaF;
					printf("loc %d %f MAX\n ", ri, omegaF);
					//printf(" MAX ");
				}

				//printf("\n");
			}
		

		}

	*/


		//R1START = 1135;  R1END = 1144; R2START = 1458;  R2END = 1463;
	/*	{ // Final refinement process 2
			int MINSNPS = 5;
			printf("\n\n FINAL REFINEMENT  2 Doing region %d %d - %d %d \n", R1START, R1END, R2START, R2END);
			double maxomega = 0.0;
			int ri;
			for(ri=R1END;ri<R2START;ri++)
			{
				int li;
				for(li=R1START;li<=ri-MINSNPS;li++)
				{
					double lLD = wLDmat[li][ri];
					
					int rri;
					for(rri=ri+1+MINSNPS;rri<=R2END;rri++)
					{
						double rLD= wLDmat[ri+1][rri];
						double cLD =  wLDmat[li][rri] - wLDmat[li][ri+1+MINSNPS-1] - wLDmat[ri-MINSNPS+1][rri] + wLDmat[ri-MINSNPS+1][ri+1+MINSNPS-1];
						int snps1 = ri - li + 1;
						double snps1sel2 = ((double)snps1*(snps1-1))/2.0;
						int snps2 = rri - (ri+1) + 1;
						double snps2sel2 = ((double)snps2*(snps2-1))/2.0;
						double snps1snps2 = ((double)snps1)*snps2;
						double omegaF = ((lLD+rLD)/(snps1sel2+snps2sel2)) / ((cLD+(1.0/(50.0*2.0)))/snps1snps2);
						if(omegaF>=maxomega)
						{	maxomega = omegaF;
							printf("%d omega %f %f %f <%f> MAX\n", ri, lLD, rLD, cLD, omegaF);
						}
					
					}
				}
			}
		}
*/

/*
		R1START = 1135;  R1END = 1144; R2START = 1458;  R2END = 1463;
		//1135 1144 - 1458 1463
		{ // Final refinement process 2
			int MINSNPS = 20;
			printf("\n\n FINAL REFINEMENT  2 Doing region %d %d - %d %d \n", R1START, R1END, R2START, R2END);
			double maxomega = 0.0;
			int ri;
			for(ri=R1END;ri<R2START;ri++)
			{
				//printf("loc %d  \n", ri);


				int li;
				for(li=R1START;li<=ri-MINSNPS;li++)
				{
					//printf("%d", li);
					// left
					double lLD = 0.0;
					/*double lLD = 0.0;
					int li1, li2;
					int snps1check = 1;
					for(li1=li;li1<=ri-1;li1++)
					{
						for(li2=li1+1;li2<=ri;li2++)
						{
							int pr1 = w->PMWindowPatternMap[li1];
							int pr2 = w->PMWindowPatternMap[li2];
							lLD += (double)pp->PMPatternPoolLDMatrix[pr1][pr2];
						}
						snps1check++;
					}
					//printf(": %f ", lLD);
					assert(wLDmat[li][ri]==lLD);
					*/
/*						lLD = wLDmat[li][ri];
					
					int rri;
					for(rri=ri+1+MINSNPS;rri<=R2END;rri++)
					{
						//if(li==R1START)
						// 	printf("%d", rri);

						double rLD=0.0;
						/*int ri1, ri2;
						int snps2check = 1;
						for(ri1=ri+1;ri1<=rri-1;ri1++)
						{
							for(ri2=ri1+1;ri2<=rri;ri2++)
							{
								int pr1 = w->PMWindowPatternMap[ri1];
								int pr2 = w->PMWindowPatternMap[ri2];
								rLD += (double)pp->PMPatternPoolLDMatrix[pr1][pr2];
							}
							snps2check++;
						}*/
/*						rLD = wLDmat[ri+1][rri];
						//if(li==R1START)
						//	printf(": %f ", rLD);
						
						// cross
						double cLD = 0.0;
						/*int ci1, ci2;
						for(ci1=li;ci1<=ri-MINSNPS;ci1++)
						{
							for(ci2=ri+1+MINSNPS;ci2<=rri;ci2++)
							{
								int pr1 = w->PMWindowPatternMap[ci1];
								int pr2 = w->PMWindowPatternMap[ci2];
								cLD += (double)pp->PMPatternPoolLDMatrix[pr1][pr2];
							}
						}
						*/
/*						cLD = wLDmat[li][rri] - wLDmat[li][ri+1+MINSNPS-1] - wLDmat[ri-MINSNPS+1][rri] + wLDmat[ri-MINSNPS+1][ri+1+MINSNPS-1];
						//assert(cLDtest==cLD);
						//if(li==R1START)
						//	printf(": %f ", cLD);

						int snps1 = ri - li + 1;
						double snps1sel2 = ((double)snps1*(snps1-1))/2.0;
						//printf("%f ", snps1sel2);
						int snps2 = rri - (ri+1) + 1;
						//printf("%d %d %d %d %d\n", snps1, snps1check, snps2, snps2check, snps1+snps2);
						double snps2sel2 = ((double)snps2*(snps2-1))/2.0;
						//printf("%f ", snps2sel2);
						double snps1snps2 = ((double)snps1)*snps2;

						

						double omegaF = ((lLD+rLD)/(snps1sel2+snps2sel2)) / ((cLD+(1.0/(50.0*2.0)))/snps1snps2);
						
						if(omegaF>=maxomega)
						{	maxomega = omegaF;
		
							printf("%d omega %f %f %f <%f> MAX\n", ri, lLD, rLD, cLD, omegaF);

							
						}
						
						
					
						
					}
					
					
				}
				//printf("\n");
				/*
				// left
				double lLD = 0.0;
				int li1, li2;
				for(li1=R1START;li1<=ri-1;li1++)
				{
					for(li2=li1+1;li2<=ri;li2++)
					{
						int pr1 = w->PMWindowPatternMap[li1];
						int pr2 = w->PMWindowPatternMap[li2];
						lLD += (double)pp->PMPatternPoolLDMatrix[pr1][pr2];
					}
				}
				//printf("lLD %f", lLD);

				// right
				double rLD = 0.0;
				int ri1, ri2;
				for(ri1=ri+1;ri1<=R2END-1;ri1++)
				{
					for(ri2=ri1+1;ri2<=R2END;ri2++)
					{
						int pr1 = w->PMWindowPatternMap[ri1];
						int pr2 = w->PMWindowPatternMap[ri2];
						rLD += (double)pp->PMPatternPoolLDMatrix[pr1][pr2];
					}
				}
				//printf("rLD %f", rLD);


				// cross
				double cLD = 0.0;
				int ci1, ci2;
				for(ci1=R1START;ci1<=ri;ci1++)
				{
					for(ci2=ri+1;ci2<=R2END;ci2++)
					{
						int pr1 = w->PMWindowPatternMap[ci1];
						int pr2 = w->PMWindowPatternMap[ci2];
						cLD += (double)pp->PMPatternPoolLDMatrix[pr1][pr2];
					}
				}
				//printf("cLD %f", cLD);


				int snps1 = ri - R1START + 1;
				double snps1sel2 = ((double)snps1*(snps1-1))/2.0;
				//printf("%f ", snps1sel2);
				int snps2 = R2END - (ri+1) + 1;
				double snps2sel2 = ((double)snps2*(snps2-1))/2.0;
				//printf("%f ", snps2sel2);
				double snps1snps2 = ((double)snps1)*snps2;
				//printf("%f ", snps1snps2);

				
				double omegaF = ((lLD+rLD)/(snps1sel2+snps2sel2)) / (cLD/snps1snps2);
				//printf("omega %f", omegaF);
				if(omegaF>=maxomega)
				{	maxomega = omegaF;
					//printf(" MAX ");
				}
				*/
				
/*			}
		

		}
*/

	//}
	fprintf(stdout, "OmegaPlus indicated index: %d\n", REFERENCE_OP_INDEX);





/*	// Do Pairwise subregions all to all
	{

		slsize = 0;
		sl1 = NULL;
		sl2 = NULL;
		sl3 = NULL;
		sl4 = NULL;
		cole = NULL;
		slomega = NULL;		

		PMWindow * w = PMPatternPool_getWindow(pp);
		int classes = pp->PMPatternPoolLDClassListSize;
		//printf("LD Classes to Process %d\n", classes);
		int mi, mj;
		for(mi=0;mi<classes-1;mi++)
		{
			for(mj=mi+1;mj<classes;mj++)
			{
				//printf("matching classes %d and %d\n", mi, mj);
				
				PMLDClass * ldclass1 = pp->PMPatternPoolLDClassList[mi];
				PMLDClass * ldclass2 = pp->PMPatternPoolLDClassList[mj];

				int r1, regions1 = ldclass1->PMLDClassValidSubregions;
				int r2, regions2 = ldclass2->PMLDClassValidSubregions;

				//printf("class %d with %d regions and class %d with %d regions\n", mi, regions1, mj, regions2);
				double OMEGA_ACCUM = 0.0;
				for(r1=0;r1<regions1;r1++)
				{	OMEGA_ACCUM = 0.0;
					for(r2=0;r2<regions2;r2++)
					{
						int r1start = ldclass1->PMLDClassValidSubregionFirstSNP[r1]-1;// + ldclass1->PMLDClassFirstSNPIndex-1;
						int r1end = ldclass1->PMLDClassValidSubregionLastSNP[r1]-1;//+ ldclass1->PMLDClassFirstSNPIndex-1;

						int r2start = ldclass2->PMLDClassValidSubregionFirstSNP[r2]-1;//+ ldclass2->PMLDClassFirstSNPIndex-1;
						int r2end = ldclass2->PMLDClassValidSubregionLastSNP[r2]-1;//+ ldclass2->PMLDClassFirstSNPIndex-1;						
						
						if((r1end<=r2start) || (r2end<=r1start))
						{	///printf("processing regions [%d %d] and [%d %d]", r1start, r1end, r2start, r2end);
							
							float LD1 = ldclass1->PMLDClassValidSubregionTotalLD[r1];
							float LD2 = ldclass2->PMLDClassValidSubregionTotalLD[r2];

							// Cross LD
							int kr1, kr2;
							double LD3 = 0.0;
							for (kr1=r1start;kr1<=r1end;kr1++)
							{
								for(kr2=r2start;kr2<=r2end;kr2++)
								{
									int pr1 = w->PMWindowPatternMap[kr1];
									int pr2 = w->PMWindowPatternMap[kr2];
									LD3 += (double) pp->PMPatternPoolLDMatrix[pr1][pr2];
								}
							}

							//LD3 /= (r1end-r1start+1+r2end-r2start+1);
	

							if(LD3==0.0)
								LD3 = 1.0/(1093.0*1093.0);

					
							printf("LD1 %f LD2 %f LD3 %f OMEGA: %e\n", LD1, LD2, LD3, (LD1+LD2)/LD3);
							experimentalRegionCount (r1start, r1end, r2start, r2end, ((double)(LD1+LD2))/LD3);
							OMEGA_ACCUM += (LD1+LD2)/LD3;
							////if(((r1end<=REFERENCE_OP_INDEX)&&(REFERENCE_OP_INDEX<=r2start))||((r2end<=REFERENCE_OP_INDEX)&&(REFERENCE_OP_INDEX<=r1start)))
							//	printf("Valid (%d)", REFERENCE_OP_INDEX);

							



							//printf("\n");
						}					
					}
					//printf("             OMEGAACCUM: %e\n", OMEGA_ACCUM);
				}	
			}
		}


	}

	//printf("and now lets check");
	int myc;
	double maxmyscore = 0;
	
	for(myc=0;myc<slsize;myc++)
	{	PMWindow * w = PMPatternPool_getWindow(pp);	
	//	printf("Region [%d %d][%d %d]: %d \t\t%e\t\t%e", sl1[myc], sl2[myc], sl3[myc], sl4[myc], cole[myc], slomega[myc], slomega[myc]*exp(cole[myc]));
		double myscore = slomega[myc]*exp(cole[myc]);
		if(myscore>=maxmyscore)
		{
			maxmyscore = myscore;
			 printf("region: [%d - %d][%d -  %d]", w->PMWindowPositionMap[sl1[myc]], w->PMWindowPositionMap[sl2[myc]], w->PMWindowPositionMap[sl3[myc]], w->PMWindowPositionMap[sl4[myc]]);
		}
		printf("\n");

	}
*/
}


