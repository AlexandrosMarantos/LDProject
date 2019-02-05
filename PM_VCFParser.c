
#include "PM.h"

int PMVCFParser_getGenotypeLocation (char * string)
{
	assert(string!=NULL);
	
	char tstring [STRING_SIZE];
	int Location = 0;
	int i, length = strlen(string);
	int startIndex = 0;
	int endIndex = 0;
	for(i=0;i<length;i++)
	{
		if((string[i]==':'))
		{
			endIndex =  i;
			memcpy(tstring, &string[startIndex], endIndex-startIndex);
			tstring[endIndex]='\0';
			if(!strcmp(tstring, "GT"))
				return Location;
			
			Location++;
			startIndex = i+1;
		}
	}
	endIndex =  i;
	memcpy(tstring, &string[startIndex], endIndex);
	tstring[endIndex]='\0';
	if(!strcmp(tstring, "GT"))
		return Location;

	return -1;	
}

static inline void PMVCFParser_getGenotypeSampleData (char * string, int Location, char * data)
{
	//assert(string!=NULL);
	//assert(Location>=0);
	//assert(data!=NULL);

	int tLocation = 0;
	int i, length = strlen(string);
	int startIndex = 0;
	int endIndex = 0;
	for(i=0;i<length;i++)
	{
		if((string[i]==':'))
		{
			endIndex =  i;
			memcpy(data, &string[startIndex], endIndex-startIndex);
			data[endIndex]='\0';

			if(tLocation==Location)
				return;
			
			tLocation++;
			startIndex = i+1;
		}
	}
	endIndex =  i;
	memcpy(data, &string[startIndex], endIndex);
	data[endIndex]='\0';
	if(tLocation==Location)
		return;	

	assert(0);
	return;
}

int PMVCFParser_parseFields (FILE * fp, char * Chromosome, int * AlleleSize, char * AlleleVector, int * GenotypeLocation, int * Position, char * ID)
{
	assert(fp!=NULL);
	assert(Chromosome!=NULL);
	assert(AlleleSize!=NULL);
	assert(AlleleVector!=NULL);
	assert(GenotypeLocation!=NULL);
	assert(Position!=NULL);
	assert(ID!=NULL);

	char tstring[STRING_SIZE];
	int ret = FAILED;

	// CHROM
	ret = fscanf(fp, "%s", tstring);
	assert(ret==SUCCESS);

	//printf("%s\n", tstring);

	if(strcmp(Chromosome, "all"))
		if(strcmp(tstring, Chromosome))
			return FAILED;

	// POSITION
	ret = fscanf(fp, "%s", tstring);
	assert(ret==SUCCESS);

	if(strlen(tstring)==1)
		*Position = tstring[0]-48;
	else
		*Position = atoi(tstring);

	//printf("%s\n", tstring);



	// ID
	ret = fscanf(fp, "%s", tstring);
	assert(ret==SUCCESS);

	strcpy(ID, tstring);
	

	//printf("%s\n", tstring);



	// REF
	ret = fscanf(fp, "%s", tstring);
	assert(ret==SUCCESS);

	if(strlen(tstring)!=1) // TODO: Fix this properly
		return FAILED;

	AlleleVector[0] = tstring[0];
	*AlleleSize = 1;
	//printf("REF: %s\n", tstring);



	// ALT
	ret = fscanf(fp, "%s", tstring);
	assert(ret==SUCCESS);

	if(strlen(tstring)!=1) // TODO: Fix this properly
		return FAILED;

	AlleleVector[1] = tstring[0];
	*AlleleSize = 2;
	//printf("%s\n", tstring);



	// QUAL
	ret = fscanf(fp, "%s", tstring);
	assert(ret==SUCCESS);

	//printf("%s\n", tstring);



	// FILTER
	ret = fscanf(fp, "%s", tstring);
	assert(ret==SUCCESS);

	//printf("%s\n", tstring);

	//if(strcmp(tstring, "PASS"))
	//	return FAILED;



	// INFO
	ret = fscanf(fp, "%s", tstring);
	assert(ret==SUCCESS);

	//printf("%s\n", tstring);



	// FORMAT
	ret = fscanf(fp, "%s", tstring);
	assert(ret==SUCCESS);

	*GenotypeLocation = PMVCFParser_getGenotypeLocation(tstring);

	//printf("GTLOC: %s\n", tstring);



	return SUCCESS;
}

void PMVCFParser_skipLine (FILE * fp)
{
	assert(fp!=NULL);

	while(fgetc(fp)!='\n');
}

void dataShuffleKnuth(char * data, int startIndex, int endIndex)
{
	if(startIndex == endIndex)
		return;

	int i, index;
	char tmp;

	for (i = endIndex; i > startIndex; i--)
	{
		index = startIndex + (rand() % (i - startIndex + 1));

		tmp = data[index];
		data[index] = data[i];
		data[i] = tmp;
	}
}

static inline void getGTdata (char * string, char * stateVector, int statesTotal, char * sampleData)
{	
	int i, j=0, index=0, start=0, end=0, len = strlen(string);

	for(i=0;i<len;i++)
	{	
		if(string[i]>=48 && string[i]<=57)
		{
			index = string[i]-48;
			
			assert(index<statesTotal);
			
			sampleData[j++] = stateVector[index];
		}
		else
		{
  			if(string[i]=='.')
			{
				sampleData[j++] = 'N';
			}

			if(string[i]=='/')
			{
				end++;
			}

			if(string[i]=='|')
			{
				dataShuffleKnuth(sampleData, start, end);
				start = j;
				end = j;
			}			
		}
	}
	dataShuffleKnuth(sampleData, start, end);
}

static inline void PMVCFParser_getSampleAlleleVector(char * string, int AlleleSize, char * AlleleVector, int GenotypeLocation, int * sampleAlleleSize, char * sampleAlleleVector, int Ploidy)
{
	//assert(string!=NULL);
	//assert(AlleleSize>=2);
	//assert(AlleleVector!=NULL);
	//assert(GenotypeLocation>=0);
	//assert(sampleAlleleSize!=NULL);
	//assert(sampleAlleleVector!=NULL);

	char data [STRING_SIZE];
	PMVCFParser_getGenotypeSampleData (string, GenotypeLocation, data);

	//int PloidyLengthCheck = 2*Ploidy-1;
	//assert(strlen(data)==PloidyLengthCheck); // TODO: Fix this properly

	getGTdata (data, AlleleVector, strlen(AlleleVector), sampleAlleleVector);
	
	*sampleAlleleSize = Ploidy; // TODO: fix this properly
	return;

	int i;
	for(i=0;i<Ploidy;i++)
	{
		sampleAlleleVector[i] = AlleleVector[data[2*i]-48];	
	}
	//sampleAlleleVector[0] = AlleleVector[data[0]-48];
	//sampleAlleleVector[1] = AlleleVector[data[2]-48];
	//sampleAlleleVector[2] = '\0';
	sampleAlleleVector[Ploidy] = '\0';

	///////char temp [STRING_SIZE];
	//getGTdata (sampleAlleleVector, AlleleVector, 2, temp);
	//strcpy(temp, sampleAlleleVector);

	//assert(0);
}

int PMVCFParser_parseSamples (FILE * fp, PMPatternPool * pp, int AlleleSize, char * AlleleVector, int GenotypeLocation, int Position, int * DerivedAlleleCount)
{
	//assert(fp!=NULL);
	//assert(pp!=NULL);
	//assert(AlleleSize>=2);
	//assert(AlleleVector!=NULL);
	//assert(GenotypeLocation>=0);

	int Ploidy = PMPatternPool_getPloidy(pp);

	char tstring[STRING_SIZE];

	int sampleAlleleSize=0;
	char sampleAlleleVector [STRING_SIZE];
	resetString(sampleAlleleVector, STRING_SIZE);

	//char * newSNP = PMPatternPool_getPendingSNPPtr(pp);
	uint64_t * newSNPCompact = PMPatternPool_getPendingSNPCompactPtr(pp);
	int newSNPCompactIndex = 0;
	int newSNPCompactEntryLoad = 0;
	//newSNP[0] = '\0';
	newSNPCompact[0] = 0llu;
	int tAlleleCount = 0;
	int i, j, SampleSize = PMPatternPool_getSampleSize(pp);
	for(i=0;i<SampleSize;i++)
	{
		int ret = fscanf(fp, "%s", tstring);
		assert(ret==SUCCESS); // TODO: Add here error with VCF line before asserting

		PMVCFParser_getSampleAlleleVector(tstring, AlleleSize, AlleleVector, GenotypeLocation, &sampleAlleleSize, sampleAlleleVector, Ploidy);
		//assert(sampleAlleleSize==Ploidy);
		//strcat(newSNP, sampleAlleleVector);

		for(j=0;j!=sampleAlleleSize;++j)
		{
			uint64_t matchIndex = 0;//(uint64_t) getMatchIndex(AlleleVector, sampleAlleleVector[j]) & 1ull;
			//matchIndex = (AlleleVector[0]==sampleAlleleVector[j]);
			//uint64_t matchIndex = (AlleleVector[0]==sampleAlleleVector[j])? 0ull : 1ull;
			int k;
			for(k=1;k!=AlleleSize;++k)
			{
				if(AlleleVector[k]==sampleAlleleVector[j])
				{
					matchIndex = k;
					break;
				}
			}	

			newSNPCompactEntryLoad++;
			newSNPCompact[newSNPCompactIndex] = (newSNPCompact[newSNPCompactIndex]<<1) | matchIndex;
			if(newSNPCompactEntryLoad==WORD_WIDTH)
			{
				tAlleleCount += pm_popcnt_u64(newSNPCompact[newSNPCompactIndex]);
				newSNPCompactIndex++;
				newSNPCompactEntryLoad=0;
				newSNPCompact[newSNPCompactIndex] = 0llu;
			}			
		}
	}
	tAlleleCount += pm_popcnt_u64(newSNPCompact[newSNPCompactIndex]);
	*DerivedAlleleCount = tAlleleCount;

	//newSNP[Ploidy*SampleSize]='\0';

	if(tAlleleCount==Ploidy*SampleSize)
		return FAILED;

	if(tAlleleCount==0)
		return FAILED;
		
	return SUCCESS;
}
