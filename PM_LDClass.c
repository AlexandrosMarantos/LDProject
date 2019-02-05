/*
	typedef struct
	{
		int 	PMLDClassSNPSize; // Total number of SNPs
		int 	PMLDClassPatternSize; // Total number of Patterns
		int * 	PMLDClassPatternList; // List of Pattern IDs for the different patterns in the class
		int	PMLDClassFirstSNPIndex;
		int	PMLDClassLastSNPIndex;
		char *	PMLDClassRegionView; // This is how the region of SNPs looks, represented as a string, when processing the current LD class

	} PMLDClass;

*/

#include "PM.h"

PMLDClass * PMLDClass_new (void)
{
	PMLDClass * ldc = NULL;
	ldc = (PMLDClass *)malloc(sizeof(PMLDClass));
	assert(ldc!=NULL);

	ldc->PMLDClassSNPSize = 0;
	ldc->PMLDClassPatternSize = 0;
	ldc->PMLDClassPatternList = NULL;
	ldc->PMLDClassFirstSNPIndex = -1;
	ldc->PMLDClassLastSNPIndex = -1;
	ldc->PMLDClassFirstSNPPosition = -1;
	ldc->PMLDClassLastSNPPosition = -1;
	ldc->PMLDClassQueryString = NULL;
	ldc->PMLDClassRegionView = NULL;

	ldc->PMLDClassMatchListSize = 0;
	ldc->PMLDClassMatchListScore = NULL;
	ldc->extra = NULL;
	ldc->PMLDClassMatchListCoordinateI = NULL;
	ldc->PMLDClassMatchListCoordinateJ = NULL;
	ldc->PMLDClassMatchListCoordinateJ2 = NULL;

	ldc->PMLDClassValidSubregions = 0;
	ldc->PMLDClassValidSubregionFirstSNP = NULL;
	ldc->PMLDClassValidSubregionLastSNP = NULL;
	ldc->PMLDClassValidSubregionTotalLD = NULL;
	ldc->PMLDClassValidSubregionAvgLD = NULL;
	ldc->PMLDClassValidSubregionUsedTimes = NULL;

	ldc->PMLDClassSubregionList = NULL;

	ldc->PMLDClassUsed = 5;

	return ldc;
}

void PMLDClass_increaseSNPSize (PMLDClass * ldc, int size)
{
	assert(ldc!=NULL);
	ldc->PMLDClassSNPSize += size;
}

int PMLDClass_getSNPSize (PMLDClass * ldc)
{
	assert(ldc!=NULL);
	return ldc->PMLDClassSNPSize;
}

void PMLDClass_free (PMLDClass * ldc)
{
	assert(ldc!=NULL);

	if(ldc->PMLDClassPatternList!=NULL)
		free(ldc->PMLDClassPatternList);

	if(ldc->PMLDClassQueryString!=NULL)
		free(ldc->PMLDClassQueryString);

	if(ldc->PMLDClassRegionView!=NULL)
		free(ldc->PMLDClassRegionView);


	if(ldc->PMLDClassMatchListScore!=NULL)
		free(ldc->PMLDClassMatchListScore);
	
	if(ldc->extra!= NULL)
		free(ldc->extra);

	if(ldc->PMLDClassMatchListCoordinateI!=NULL)
		free(ldc->PMLDClassMatchListCoordinateI);

	if(ldc->PMLDClassMatchListCoordinateJ!=NULL)
		free(ldc->PMLDClassMatchListCoordinateJ);
	
	if(ldc->PMLDClassMatchListCoordinateJ2!=NULL)
		free(ldc->PMLDClassMatchListCoordinateJ2);

	if(ldc->PMLDClassValidSubregionFirstSNP!=NULL)
		free(ldc->PMLDClassValidSubregionFirstSNP);

	if(ldc->PMLDClassValidSubregionLastSNP!=NULL)
		free(ldc->PMLDClassValidSubregionLastSNP);
	
	if(ldc->PMLDClassValidSubregionTotalLD!=NULL)
		free(ldc->PMLDClassValidSubregionTotalLD);

	if(ldc->PMLDClassValidSubregionAvgLD!=NULL)
		free(ldc->PMLDClassValidSubregionAvgLD);

	if(ldc->PMLDClassSubregionList!=NULL)
	{
		int i;
		for(i=0;i<ldc->PMLDClassValidSubregions;i++)
		{
			PMLDSubregion * t = ldc->PMLDClassSubregionList[i];
			if(t->PMLDSubregionPatternIDWindow!=NULL)
				free(t->PMLDSubregionPatternIDWindow);

			free(t);
		}
		free(ldc->PMLDClassSubregionList);		
	}

	free(ldc);
}

void PMLDClass_print(PMLDClass * ldc, FILE * fp, int index, void * pp)
{
	assert(ldc!=NULL);
	assert(fp!=NULL);








/* */	PMPatternPool * ppt = (PMPatternPool *)pp;
	PMWindow * w = PMPatternPool_getWindow(ppt);

	fprintf(fp, "[class %d] ", index);
	int ll, myindex = 0;
	for(ll=0;ll<strlen(ldc->PMLDClassRegionView);ll++)
	{
		if(ldc->PMLDClassRegionView[ll]=='A')
		{
			fprintf(fp, "%d ", w->PMWindowPositionMap[ll]);		
		}
	}
	fprintf(fp, "\n");






/* *

	//fprintf(fp, "REFERENCE OP POINT: %d\n", REFERENCE_OP_INDEX);
	fprintf(fp, "\nCLASS INFORMATION [ID: %d]\n[SNPs: %d]\n[PATTERNs: %d]\n", index, ldc->PMLDClassSNPSize, ldc->PMLDClassPatternSize);
	//fprintf(fp, "Total SNPs:            %d\n", ldc->PMLDClassSNPSize);
	//fprintf(fp, "Total Patterns:        %d\n", ldc->PMLDClassPatternSize);
	int i;
	fprintf(fp, "[PATTERN INDECES: ");
	for(i=0;i<ldc->PMLDClassPatternSize;i++)
	{
		fprintf(fp, "%d", ldc->PMLDClassPatternList[i]);
		if(i<ldc->PMLDClassPatternSize-1)
			fprintf(fp, ", ");
		else
			fprintf(fp, "]\n");
	}
	fprintf(fp, "[SNP RANGE: %d-%d]\n",ldc->PMLDClassFirstSNPIndex, ldc->PMLDClassLastSNPIndex);
	//fprintf(fp, "First SNP index:       %d\n", ldc->PMLDClassFirstSNPIndex);
	//fprintf(fp, "Last SNP index:        %d\n", ldc->PMLDClassLastSNPIndex);
	fprintf(fp, "[REGION: %d-%d]\n", ldc->PMLDClassFirstSNPPosition, ldc->PMLDClassLastSNPPosition);
	//fprintf(fp, "First SNP position:    %d\n", ldc->PMLDClassFirstSNPPosition);
	//fprintf(fp, "Last SNP position:     %d\n", ldc->PMLDClassLastSNPPosition);
	fprintf(fp, "VIEW:");
	int ll, myindex = 0;
	int realSNPIndex = ldc->PMLDClassFirstSNPIndex;
	for(ll=0;ll<strlen(ldc->PMLDClassRegionView);ll++)
	{
		if(ll%100==0)
		{
			printf("\n>[%d:%d]:\t\t", (myindex)*100, (myindex+1)*100-1);
			myindex++;
		}
		printf("%c",ldc->PMLDClassRegionView[ll]);
		realSNPIndex++;
		//if(realSNPIndex==REFERENCE_OP_INDEX)
		//	printf(" ");
	}
	fprintf(fp, "\n");
	fprintf(fp, "SUBREGIONs: (%d)\n", ldc->PMLDClassValidSubregions);
	for(i=0;i<ldc->PMLDClassValidSubregions;i++)
	{
	fprintf(fp, "\nSubregion%d:            %d - %d\n", i, ldc->PMLDClassValidSubregionFirstSNP[i], ldc->PMLDClassValidSubregionLastSNP[i]);
	fprintf(fp, "Total_pw_ld:           %f\n", ldc->PMLDClassValidSubregionTotalLD[i]);
	fprintf(fp, "Avg_pw_ld:             %f\n", ldc->PMLDClassValidSubregionAvgLD[i]);
	fprintf(fp, "Subregion check:");
		myindex = 0;
		for(ll=0;ll<strlen(ldc->PMLDClassRegionView);ll++)
		{
			if(ll%100==0)
			{
				printf("\n>[%d:%d]:\t\t", (myindex)*100, (myindex+1)*100-1);
				myindex++;
			}
			if(ll>=ldc->PMLDClassValidSubregionFirstSNP[i]-1 && ll<=ldc->PMLDClassValidSubregionLastSNP[i]-1)
				printf("X");
			else
				printf("%c",ldc->PMLDClassRegionView[ll]);
		}
		fprintf(fp, "\n");

		fprintf(fp, "Subregion Pattern check:\n");
		myindex = 0;
		PMLDSubregion * tsubregion = NULL;
		
			tsubregion = ldc->PMLDClassSubregionList[i];
			printf("L %d\t\t\t", tsubregion->PMLDSubregionLength);
			for(myindex=0;myindex<tsubregion->PMLDSubregionLength;myindex++)
				printf("%d ", tsubregion->PMLDSubregionPatternIDWindow[myindex]);
		
		printf("\n");
	}
/* */	
}

void PMLDClass_appendPattern (PMLDClass * ldc, int PatternID)
{
	assert(ldc!=NULL);
	assert(PatternID>=0);

	ldc->PMLDClassPatternSize++;
	ldc->PMLDClassPatternList = realloc(ldc->PMLDClassPatternList, sizeof(int)*ldc->PMLDClassPatternSize);
	assert(ldc->PMLDClassPatternList!=NULL);
	ldc->PMLDClassPatternList[ldc->PMLDClassPatternSize-1] = PatternID;	
}

void PMLDClass_setFirstSNPIndex (PMLDClass * ldc, void * pp, int patternIndex)
{
	assert(ldc!=NULL);
	assert(pp!=NULL);
	assert(patternIndex>=0);

	PMPatternPool * ppp = (PMPatternPool *)pp; 
	int firstPatternIndex = ppp->PMPatternPoolSNPIndeces[patternIndex][0];
	if(ldc->PMLDClassFirstSNPIndex==-1)
		ldc->PMLDClassFirstSNPIndex = firstPatternIndex;
	else
		if(firstPatternIndex<ldc->PMLDClassFirstSNPIndex)
			ldc->PMLDClassFirstSNPIndex = firstPatternIndex;

	PMWindow * w = PMPatternPool_getWindow (ppp);
	ldc->PMLDClassFirstSNPPosition = w->PMWindowPositionMap[ldc->PMLDClassFirstSNPIndex];
}

void PMLDClass_setLastSNPIndex (PMLDClass * ldc, void * pp, int patternIndex)
{
	assert(ldc!=NULL);
	assert(pp!=NULL);
	assert(patternIndex>=0);

	PMPatternPool * ppp = (PMPatternPool *)pp; 
	int lastPatternIndex = ppp->PMPatternPoolSNPIndeces[patternIndex][PMPatternPool_getCountByIndex(ppp, patternIndex)-1];
	if(ldc->PMLDClassLastSNPIndex==-1)
		ldc->PMLDClassLastSNPIndex = lastPatternIndex;
	else
		if(lastPatternIndex>ldc->PMLDClassLastSNPIndex)
			ldc->PMLDClassLastSNPIndex = lastPatternIndex;

	PMWindow * w = PMPatternPool_getWindow (ppp);
	ldc->PMLDClassLastSNPPosition = w->PMWindowPositionMap[ldc->PMLDClassLastSNPIndex];
}



void PMLDClass_createQueryStringAndRegionView(PMLDClass * ldc, void * pp)
{
	assert(ldc!=NULL);
	assert(pp!=NULL);
	PMPatternPool * ppp = (PMPatternPool *)pp;
	PMWindow * w = PMPatternPool_getWindow (ppp);

	ldc->PMLDClassQueryString = realloc(ldc->PMLDClassQueryString, sizeof(char)*(ldc->PMLDClassSNPSize+1));
	assert(ldc->PMLDClassQueryString!=NULL);
	ldc->PMLDClassQueryString[ldc->PMLDClassSNPSize] = '\0';
	int i;
	for(i=0;i<ldc->PMLDClassSNPSize;i++)
		ldc->PMLDClassQueryString[i] = 'A';
	
	assert(i==ldc->PMLDClassSNPSize);

	int j;//, length =  ldc->PMLDClassLastSNPIndex - ldc->PMLDClassFirstSNPIndex + 2;
	//ldc->PMLDClassRegionView = realloc(ldc->PMLDClassRegionView, sizeof(char)*length);
	ldc->PMLDClassRegionView = realloc(ldc->PMLDClassRegionView, sizeof(char)*(PMWindow_getSites(w)+1));
	assert(ldc->PMLDClassRegionView!=NULL);
	//ldc->PMLDClassRegionView[length-1] = '\0';
	ldc->PMLDClassRegionView[PMWindow_getSites(w)] = '\0';
	int index = 0;
	//for(i=ldc->PMLDClassFirstSNPIndex;i<=ldc->PMLDClassLastSNPIndex;i++)
	for(i=0;i<=PMWindow_getSites(w)-1;i++)
	{
		int patternIndex = w->PMWindowPatternMap[i];
		//printf("%d patternIndex: %d\n", i, patternIndex);

		for(j=0;j<ldc->PMLDClassPatternSize;j++)
		{
			if(patternIndex==ldc->PMLDClassPatternList[j])
			{
				ldc->PMLDClassRegionView[index++] = 'A';//j+48;//'A';
				j=ldc->PMLDClassPatternSize+1;
			}
		}
		if(j==ldc->PMLDClassPatternSize)
		{	//if(ppp->PMPatternPoolLDMatrix[0][patternIndex]>=0.6 && ppp->PMPatternPoolLDMatrix[0][patternIndex]<0.9)
			//	ldc->PMLDClassRegionView[index++] = 'C';
			//else
			//	if(ppp->PMPatternPoolLDMatrix[0][patternIndex]<0.6 && ppp->PMPatternPoolLDMatrix[0][patternIndex]>=0.3)
			//		ldc->PMLDClassRegionView[index++] = 'G';
			//	else
					ldc->PMLDClassRegionView[index++] = '_';				
		}

		/*if((index>=10 && index<=20)||(index>=40 && index<=60))
			ldc->PMLDClassRegionView[index++] = 'A';
		else
			ldc->PMLDClassRegionView[index++] = 'C';
		*/					
	}
	//assert(index==length-1);

	//printf("the view is: %s\n", ldc->PMLDClassRegionView);	
	//assert(0);	
}

void PMLDClass_keepMatchScore (PMLDClass * ldc, int max_val, int max_i, int max_j)
{
	assert(ldc!=NULL);	

	if(ldc->PMLDClassMatchListSize==0)
	{
		ldc->PMLDClassMatchListSize++;
		ldc->PMLDClassMatchListScore = realloc(ldc->PMLDClassMatchListScore, sizeof(int)*ldc->PMLDClassMatchListSize);
		ldc->extra = realloc(ldc->extra, sizeof(int)*ldc->PMLDClassMatchListSize);
		ldc->PMLDClassMatchListCoordinateI = realloc(ldc->PMLDClassMatchListCoordinateI, sizeof(int)*ldc->PMLDClassMatchListSize);
		ldc->PMLDClassMatchListCoordinateJ = realloc(ldc->PMLDClassMatchListCoordinateJ, sizeof(int)*ldc->PMLDClassMatchListSize);
		ldc->PMLDClassMatchListCoordinateJ2 = realloc(ldc->PMLDClassMatchListCoordinateJ2, sizeof(int)*ldc->PMLDClassMatchListSize);
		
		ldc->PMLDClassMatchListScore [ldc->PMLDClassMatchListSize-1] = max_val;
		ldc->extra [ldc->PMLDClassMatchListSize-1] = max_val;
		ldc->PMLDClassMatchListCoordinateI [ldc->PMLDClassMatchListSize-1] = max_i;
		ldc->PMLDClassMatchListCoordinateJ [ldc->PMLDClassMatchListSize-1] = max_j;
		ldc->PMLDClassMatchListCoordinateJ2 [ldc->PMLDClassMatchListSize-1] = max_j;
	}
	else
	{
		int i;
		for(i=0;i<ldc->PMLDClassMatchListSize;i++)
		{
			if(ldc->PMLDClassMatchListCoordinateJ[i]==max_j)
			{
				if(max_i>=ldc->PMLDClassMatchListCoordinateI[i] || max_val>ldc->PMLDClassMatchListScore[i])
				{
					ldc->PMLDClassMatchListScore[i] = max_val;
					ldc->extra[i] = max_val;
					ldc->PMLDClassMatchListCoordinateI[i] = max_i;
					i=ldc->PMLDClassMatchListSize+1;
				}
			}
		}
		if(i==ldc->PMLDClassMatchListSize)
		{
			ldc->PMLDClassMatchListSize++;
			ldc->PMLDClassMatchListScore = realloc(ldc->PMLDClassMatchListScore, sizeof(int)*ldc->PMLDClassMatchListSize);
			ldc->extra = realloc(ldc->extra, sizeof(int)*ldc->PMLDClassMatchListSize);
			ldc->PMLDClassMatchListCoordinateI = realloc(ldc->PMLDClassMatchListCoordinateI, sizeof(int)*ldc->PMLDClassMatchListSize);
			ldc->PMLDClassMatchListCoordinateJ = realloc(ldc->PMLDClassMatchListCoordinateJ, sizeof(int)*ldc->PMLDClassMatchListSize);
			ldc->PMLDClassMatchListCoordinateJ2 = realloc(ldc->PMLDClassMatchListCoordinateJ2, sizeof(int)*ldc->PMLDClassMatchListSize);
			ldc->PMLDClassMatchListScore [ldc->PMLDClassMatchListSize-1] = max_val;
			ldc->extra [ldc->PMLDClassMatchListSize-1] = max_val;
			ldc->PMLDClassMatchListCoordinateI [ldc->PMLDClassMatchListSize-1] = max_i;
			ldc->PMLDClassMatchListCoordinateJ [ldc->PMLDClassMatchListSize-1] = max_j;
			ldc->PMLDClassMatchListCoordinateJ2 [ldc->PMLDClassMatchListSize-1] = max_j;
		}
	}
} 

int isInLDClass (PMLDClass * ldc, int pID)
{
	int i;
	for(i=0;i<ldc->PMLDClassPatternSize;i++)
		if(ldc->PMLDClassPatternList[i]==pID)
			return pID;

	return -1;
}

void PMLDClass_matchQuery2Region (PMLDClass * ldc, PMWindow * w, void * pp)
{
	assert(ldc!=NULL);

	PMPatternPool * ppp = (PMPatternPool*)pp;
	
	//strcpy(ldc->PMLDClassRegionView, "BBBBAAAAAABBBBBBBBBBB");
	//strcpy(ldc->PMLDClassQueryString, "AAAA");
	
	char * string1 = ldc->PMLDClassQueryString;
	char * string2 = ldc->PMLDClassRegionView;

	int i, length1 = strlen(string1);
	int j, length2 = strlen(string2);

	//fprintf(stdout, "string matching: \n%s\n%s\n%d %d\n", string1, string2, length1, length2);
	//assert(0);

	int val_up, val_le, val_diag;

	//assert(length1+1<PATTERNS_MAX);
	//assert(length2+1<PATTERNS_MAX);
	/*int ** mat = (int**)malloc(sizeof(int*)*(length1+1));
	for(i=0;i<length1+1;i++)
	{	mat[i] = (int*)malloc(sizeof(int)*(length2+1));
		for(j=0;j<length2+1;j++)
			mat[i][j] = 0;
	}*/
	int max_val=-1;
	int max_i=-1;
	int max_j=-1;

	// Matrix
	/*for(i=1;i<length1+1;i++)
	{
		//printf("\n%d:\t", i);
		for(j=1;j<length2+1;j++)
		{
			val_up = mat[i-1][j] + INDEL_PEN;
			val_le = mat[i][j-1] + INDEL_PEN;
			
			if(string1[i-1]==string2[j-1])
				val_diag = mat[i-1][j-1] + MATCH_SCR;
			else
				val_diag = mat[i-1][j-1] + MISMATCH_PEN;

			mat[i][j] = MAX3(val_diag, val_up, val_le);
			//printf("%d ", mat[i][j]);
			if(mat[i][j]>=max_val)
			{	
				max_val = mat[i][j]; 
				max_i = i;
				max_j = j;
				
				//PMLDClass_keepMatchScore (ldc, max_val, max_i, max_j);
				//printf("NEW: max val %d max i %d max j %d\n", max_val, max_i, max_j);
			}
		}
	}*/
	//assert(0);
	PMLDClass_keepMatchScore (ldc, max_val, max_i, max_j);

	//printf("\nEND: max val %d max i %d max j %d\n", max_val, max_i, max_j);
	int ll;
	int listsize = 1;//ldc->PMLDClassMatchListSize;
	for(ll=0;ll<listsize;ll++)
	{
		//printf("\nM'%d: %d %d %d ->", ll, ldc->PMLDClassMatchListScore[ll], ldc->PMLDClassMatchListCoordinateI[ll], ldc->PMLDClassMatchListCoordinateJ[ll]);
		max_i = ldc->PMLDClassMatchListCoordinateI[ll];
		max_j = ldc->PMLDClassMatchListCoordinateJ[ll];

		/*// Trace back
		int cur_score = mat[max_i][max_j];
		int up_score = mat[max_i-1][max_j];
		int diag_score = mat[max_i-1][max_j-1];
		int left_score = mat[max_i][max_j-1];

		while(cur_score>0)
		{
			if(MAX3(up_score, diag_score, left_score)==up_score)
			{	
				max_i--;
			}
			else
			{
				if(MAX3(up_score, diag_score, left_score)==diag_score)
				{
					max_i--;
					max_j--;	
				}
				else
				{
					max_j--;
				}	
			}

		if(max_i==0)
			break;
	
		if(max_j==0)
			break;

		 cur_score = mat[max_i][max_j];
		 up_score = mat[max_i-1][max_j];
		 diag_score = mat[max_i-1][max_j-1];
		 left_score = mat[max_i][max_j-1];		
		
		}
		*/
		//printf("\n\n %d %d (length in SNPs: %d)\n", max_j, ldc->PMLDClassMatchListCoordinateJ[ll], ldc->PMLDClassMatchListCoordinateJ[ll]-max_j+1);
		//fflush(stdout);

		/**/int ex=0;
		float score = 0.0;
		float maxscore2process = 0.0;
		int maxpoint = -1;
		//printf("scores: ");
		float * scorelist = (float *) malloc(sizeof(float)*length2);
		for(ex=0;ex<length2;ex++)
		{
			if(string2[ex]=='A')
				score += MATCH_SCR;
			else
			{
				if(score==0)
					score += 0.0;//MISMATCH_PEN;
				else
					score += MISMATCH_PEN;
			}
			//printf("%f ", score);
			scorelist[ex] = score;
			if(score>=maxscore2process)
			{
				maxscore2process = score;
				maxpoint = ex;
			}
		}
		//printf("end point %d (%d)\n", maxpoint+1, ldc->PMLDClassMatchListCoordinateJ[ll]);
		//assert(0);
		// scan back
		ex = maxpoint;
		while(scorelist[ex]>0)
			ex--;

		
		//printf("start point %d (%d)\n", ex, max_j);
		 ldc->PMLDClassMatchListCoordinateJ[ll] = maxpoint+1;
		max_j = ex+1;

		/**/

		ldc->PMLDClassMatchListCoordinateJ2[ll] = max_j;

		assert(ldc->PMLDClassMatchListCoordinateJ2[ll]<=ldc->PMLDClassMatchListCoordinateJ[ll]);

		double totalLD = 0.0;
		if(ldc->PMLDClassMatchListCoordinateJ[ll]-max_j+1>=MINIMUM_SNPS_IN_LD_REGION)
		{
			int lli, llj;
			for(lli=max_j;lli<=ldc->PMLDClassMatchListCoordinateJ[ll]-1;lli++)
			{
			 	int pID1 = w->PMWindowPatternMap[lli-1];//+ldc->PMLDClassFirstSNPIndex-1];
				if(pID1!=-1)
					for(llj=lli+1;llj<=ldc->PMLDClassMatchListCoordinateJ[ll];llj++)
					{
						int pID2 = w->PMWindowPatternMap[llj-1];//+ldc->PMLDClassFirstSNPIndex-1];
						
						if(pID2!=-1)
							totalLD += (double)ppp->PMPatternPoolLDMatrix[pID1][pID2];					
					}					
			}

		//assert(0);
		//	printf("TotalLD: %f ", totalLD);
			ldc->extra[ll] = ldc->PMLDClassMatchListScore[ll];
			ldc->PMLDClassMatchListScore[ll] = 1;

			ldc->PMLDClassValidSubregions++;
			ldc->PMLDClassValidSubregionFirstSNP = realloc(ldc->PMLDClassValidSubregionFirstSNP, sizeof(int)*ldc->PMLDClassValidSubregions);
			ldc->PMLDClassValidSubregionLastSNP = realloc(ldc->PMLDClassValidSubregionLastSNP, sizeof(int)*ldc->PMLDClassValidSubregions);
			ldc->PMLDClassValidSubregionTotalLD = realloc(ldc->PMLDClassValidSubregionTotalLD, sizeof(int)*ldc->PMLDClassValidSubregions);
			ldc->PMLDClassValidSubregionAvgLD = realloc(ldc->PMLDClassValidSubregionAvgLD, sizeof(int)*ldc->PMLDClassValidSubregions);
			ldc->PMLDClassValidSubregionFirstSNP[ldc->PMLDClassValidSubregions-1] = max_j;
			ldc->PMLDClassValidSubregionLastSNP[ldc->PMLDClassValidSubregions-1] =  ldc->PMLDClassMatchListCoordinateJ[ll];
			ldc->PMLDClassValidSubregionTotalLD[ldc->PMLDClassValidSubregions-1] = totalLD;
			ldc->PMLDClassValidSubregionAvgLD[ldc->PMLDClassValidSubregions-1] = totalLD/(ldc->PMLDClassMatchListCoordinateJ[ll]-max_j+1);

			ldc->PMLDClassSubregionList = realloc(ldc->PMLDClassSubregionList, sizeof(PMLDSubregion*)*ldc->PMLDClassValidSubregions);
			ldc->PMLDClassSubregionList[ldc->PMLDClassValidSubregions-1] = (PMLDSubregion*)malloc(sizeof(PMLDSubregion));
			
			PMLDSubregion * subregionTemp = ldc->PMLDClassSubregionList[ldc->PMLDClassValidSubregions-1];
			subregionTemp->PMLDSubregionLength = ldc->PMLDClassMatchListCoordinateJ[ll]-max_j+1;
			subregionTemp->PMLDSubregionPatternIDWindow = (int*)malloc(sizeof(int)*subregionTemp->PMLDSubregionLength);

			int pi, pID;
			for(pi=max_j;pi<=ldc->PMLDClassMatchListCoordinateJ[ll];pi++)
			{
				int pID = w->PMWindowPatternMap[pi-1];
				subregionTemp->PMLDSubregionPatternIDWindow[pi-max_j] = isInLDClass (ldc, pID);
				//printf("%d: %d\n", pi, pID);	
			}		
 
			//for(pi=0;pi<subregionTemp->PMLDSubregionLength;pi++)
			//	printf("%d: %d\n", pi, subregionTemp->PMLDSubregionPatternIDWindow[pi]);

			//assert(0);
		}
		else
		{
		//	printf("SKIPPED");
			ldc->extra[ll] = ldc->PMLDClassMatchListScore[ll];
			ldc->PMLDClassMatchListScore[ll] = 0;
		}
		
	}
	
	
	//printf("START: max val %d max i %d max j %d\n", max_val, max_i, max_j);	


	//for(i=0;i<length1+1;i++)
	//{	free(mat[i]);
	//}
	//free(mat);
}

void PMLDClass_processSubregions (PMLDClass * ldc, PMWindow * w, void * pp)
{
	assert(ldc!=NULL);

	//PMPatternPool * ppp = (PMPatternPool*)pp;

	int windowSize = PMWindow_getSNPs(w);
	int i;
	for(i=0;i<windowSize;i++)
	{
		if(w->PMWindowPositionMap[i]>71345376.00)
		{
			REFERENCE_OP_INDEX = i;
			break;
		}	
	}

	int totalSubregions = ldc->PMLDClassValidSubregions;
	printf("total subregions %d\n", totalSubregions);

	

}

