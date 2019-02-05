/*

	typedef struct
	{
		int 	PMInputType;
		int	PMInputSampleSize;
		char ** PMInputSampleName;
		char	PMInputChromosome[STRING_SIZE];
		int	PMInputFirstDataLineNumber;
	

	} PMInput;

*/

#include "PM.h"

PMInput *PMInput_new(void)
{
	PMInput *i = NULL;
	i = (PMInput *)malloc(sizeof(PMInput));
	assert(i != NULL);

	i->PMInputType = INVALID;
	i->PMInputSampleSize = 0;
	i->PMInputSampleName = NULL;
	strcpy(i->PMInputChromosome, "all");

	return i;
}

void PMInput_free(PMInput *i)
{
	assert(i != NULL);

	if (i->PMInputSampleName != NULL)
	{
		int j, size = i->PMInputSampleSize;
		for (j = 0; j < size; j++)
			free(i->PMInputSampleName[j]);

		free(i->PMInputSampleName);
	}

	free(i);
}

void PMInput_setType(PMInput *i, int Type)
{
	assert(i != NULL);
	assert((Type == INPUT_TYPE_VCF));

	i->PMInputType = Type;
}

int PMInput_getType(PMInput *i)
{
	assert(i != NULL);

	return i->PMInputType;
}

int PMInput_getSampleSize(PMInput *i)
{
	assert(i != NULL);

	return i->PMInputSampleSize;
}

void PMInput_setSampleSize(PMInput *i, int SampleSize)
{
	assert(i != NULL);
	assert(SampleSize >= 0);

	i->PMInputSampleSize = SampleSize;
}

void PMInput_appendSampleName(PMInput *i, char *SampleName)
{
	assert(i != NULL);
	assert(SampleName != NULL);

	int size = PMInput_getSampleSize(i);
	size++;
	PMInput_setSampleSize(i, size);

	i->PMInputSampleName = realloc(i->PMInputSampleName, sizeof(char **) * size);
	assert(i->PMInputSampleName);

	i->PMInputSampleName[size - 1] = (char *)malloc(sizeof(char) * STRING_SIZE);
	assert(i->PMInputSampleName[size - 1]);
	strcpy(i->PMInputSampleName[size - 1], SampleName);	//SampleName = sample864 etc
}

void PMInput_loadSampleSize(PMInput *i, FILE *fp)
{
	assert(i != NULL);
	assert(fp != NULL);

	int Type = PMInput_getType(i);
	char tstring[STRING_SIZE];
	char fchar = '\0';
	if (Type == INPUT_TYPE_VCF)
	{
		int skip = VCF_FIELDS;	//8
		while (skip != 0)
		{
			int ret = fscanf(fp, "%s", tstring);
			assert(ret == SUCCESS);
			--skip;
		}

		while (fchar != '\n')
		{
			int ret = fscanf(fp, "%s", tstring);
			assert(ret == SUCCESS);
			PMInput_appendSampleName(i, tstring);
			fchar = fgetc(fp);
		}
	}
	else
	{
		assert(Type == INPUT_TYPE_VCF);
	}
}

void PMInput_setFirstDataLineNumber(PMInput *input, int LineNumber)
{
	assert(input != NULL);

	input->PMInputFirstDataLineNumber = LineNumber;
}

int PMInput_getFirstDataLineNumber(PMInput *input)
{
	assert(input != NULL);

	return input->PMInputFirstDataLineNumber;
}

FILE *PMInput_check(PMInput *i, PMCommandLine *cl)
{
	assert(i != NULL);
	assert(cl != NULL);

	char myFlagName[STRING_SIZE];
	strcpy(myFlagName, "input");

	int myFlagIndex = INVALID;
	myFlagIndex = PMCommandLine_getValidFlagNameIndex(cl, myFlagName);
	assert(myFlagIndex != INVALID);

	char InputPath[STRING_SIZE];
	PMCommandLine_getUserFlagInput(cl, myFlagIndex, InputPath);

	FILE *fp = fopen(InputPath, "r");
	if (fp != NULL)
	{
		PMInput_setType(i, INPUT_TYPE_VCF);
		int FileFlagLine = getFileFlagLineNumber(fp, "#CHROM");
		assert(FileFlagLine != INVALID);
		PMInput_setFirstDataLineNumber(i, FileFlagLine + 1);
		PMInput_loadSampleSize(i, fp);
		return fp;
	}
	else
	{
		ReportTerminalErrorMessage(stderr, "Missing input file");
		return fp;
	}
	return fp;
}

void PMInput_getChromosome(PMInput *i, char *Chromosome)
{
	assert(i != NULL);
	assert(Chromosome != NULL);

	strcpy(Chromosome, i->PMInputChromosome);
}

void PMInput_print(PMInput *i, FILE *fp)
{
	assert(i != NULL);
	assert(fp != NULL);

	fprintf(fp, "INPUT INFORMATION\n");

	int Type = PMInput_getType(i);
	char TypeS[STRING_SIZE];
	getTypeString(Type, TypeS);
	fprintf(fp, "[format] %s\n", TypeS);

	int SampleSize = PMInput_getSampleSize(i);
	fprintf(fp, "[samples] %d\n", SampleSize);

	char Chromosome[STRING_SIZE];
	PMInput_getChromosome(i, Chromosome);
	fprintf(fp, "[chromosome] %s\n", Chromosome);

	fprintf(fp, "\n");
}
