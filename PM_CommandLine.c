/*

	typedef struct
	{
		int	PMCommnadLineValidFlagSize;
		char **	PMCommandLineValidFlagName;
		int   * PMCommandLineValidFlagType;
		char ** PMCommandLineUserFlagInput;

	} PMCommandLine;

*/

#include "PM.h"

int PMCommandLine_getValidFlagSize (PMCommandLine * cl)
{
	assert(cl!=NULL);

	return cl->PMCommnadLineValidFlagSize;
}

void PMCommandLine_setValidFlagSize (PMCommandLine * cl, int Size)
{
	assert(cl!=NULL);
	assert(Size>=0);

	cl->PMCommnadLineValidFlagSize = Size;
}

void PMCommandLine_appendValidFlagName (PMCommandLine * cl, char * FlagName, int Type)
{
	assert(cl!=NULL);
	assert(FlagName!=NULL);
	assert((Type==REQUIRED)||(Type==OPTIONAL));

	int size = PMCommandLine_getValidFlagSize (cl);
	size++;
	PMCommandLine_setValidFlagSize (cl, size);

	cl->PMCommandLineValidFlagName = realloc(cl->PMCommandLineValidFlagName, sizeof(char**)*size);
	assert(cl->PMCommandLineValidFlagName!=NULL);

	cl->PMCommandLineValidFlagName[size-1] = (char*)malloc(sizeof(char)*STRING_SIZE);
	assert(cl->PMCommandLineValidFlagName[size-1]!=NULL);

	cl->PMCommandLineValidFlagType = realloc(cl->PMCommandLineValidFlagType, sizeof(int*)*size);
	assert(cl->PMCommandLineValidFlagType!=NULL);

	cl->PMCommandLineUserFlagInput = realloc(cl->PMCommandLineUserFlagInput, sizeof(char**)*size);
	assert(cl->PMCommandLineUserFlagInput!=NULL);

	cl->PMCommandLineUserFlagInput[size-1] = (char*)malloc(sizeof(char)*STRING_SIZE);
	assert(cl->PMCommandLineUserFlagInput[size-1]!=NULL);

	strcpy(cl->PMCommandLineValidFlagName[size-1], FlagName);
	cl->PMCommandLineValidFlagType[size-1] = Type;
	strcpy(cl->PMCommandLineUserFlagInput[size-1], "\0");
}

void PMCommandLine_create (PMCommandLine * cl)
{
	assert(cl!=NULL);

	PMCommandLine_appendValidFlagName (cl, "name", REQUIRED); // Run name
	PMCommandLine_appendValidFlagName (cl, "input", REQUIRED); // Input file
	PMCommandLine_appendValidFlagName (cl, "ploidy", REQUIRED); // Ploidy
	PMCommandLine_appendValidFlagName (cl, "ldthreshold", REQUIRED); // LDthreshold

	// TODO add more input arguments here
	
}

PMCommandLine * PMCommandLine_new (void)
{
	PMCommandLine * cl = NULL;

	cl = (PMCommandLine *)malloc(sizeof(PMCommandLine));
	assert(cl!=NULL);

	cl->PMCommnadLineValidFlagSize = 0;
	cl->PMCommandLineValidFlagName = NULL;
	cl->PMCommandLineValidFlagType = NULL;
	cl->PMCommandLineUserFlagInput = NULL;

	PMCommandLine_create (cl);

	return cl;
}

void PMCommandLine_free (PMCommandLine * cl)
{
	int i;
	int size = cl->PMCommnadLineValidFlagSize;
	if(size>0)
	{
		assert(cl->PMCommandLineValidFlagName!=NULL);

		for(i=0;i<size;i++)
			free(cl->PMCommandLineValidFlagName[i]);

		free(cl->PMCommandLineValidFlagName);

		if(cl->PMCommandLineUserFlagInput!=NULL)
		{
			for(i=0;i<size;i++)
			{	
				if(cl->PMCommandLineUserFlagInput[i]!=NULL)
					free(cl->PMCommandLineUserFlagInput[i]);		
			}

			free(cl->PMCommandLineUserFlagInput);
		}

		if(cl->PMCommandLineValidFlagType!=NULL)
			free(cl->PMCommandLineValidFlagType);		
	}

	free(cl);
}

void PMCommandLine_getValidFlagName (PMCommandLine * cl, int index, char * FlagName)
{
	assert(cl!=NULL);
	int size = PMCommandLine_getValidFlagSize(cl);
	assert(index>=0);
	assert(index<size);
	assert(FlagName!=NULL);

	strcpy(FlagName, cl->PMCommandLineValidFlagName[index]);
}

void PMCommandLine_getUserFlagInput (PMCommandLine * cl, int index, char * FlagInput)
{
	assert(cl!=NULL);
	int size = PMCommandLine_getValidFlagSize(cl);
	assert(index>=0);
	assert(index<size);
	assert(FlagInput!=NULL);

	strcpy(FlagInput, cl->PMCommandLineUserFlagInput[index]);
}

void PMCommandLine_setUserFlagInput (PMCommandLine * cl, int index, char * FlagInput)
{
	assert(cl!=NULL);
	int size = PMCommandLine_getValidFlagSize(cl);
	assert(index>=0);
	assert(index<size);
	assert(FlagInput!=NULL);

	strcpy(cl->PMCommandLineUserFlagInput[index], FlagInput);
}

void PMCommandLine_print (PMCommandLine * cl, FILE * fp)
{
	assert(cl!=NULL);
	int i;
	int size = PMCommandLine_getValidFlagSize(cl);
	fprintf(fp, "\nEXECUTION INFORMATION\n");
	for(i=0;i<size;i++)
	{
		char FlagName[STRING_SIZE];
		PMCommandLine_getValidFlagName (cl, i, FlagName);
		char UserInput[STRING_SIZE];
		PMCommandLine_getUserFlagInput (cl, i, UserInput);

		fprintf(fp, "[%s] %s\n", FlagName, UserInput);
	}
	fprintf(fp, "\n");
}

int  PMCommandLine_getValidFlagNameIndex (PMCommandLine * cl, char * currentFlagName)
{
	assert(cl!=NULL);
	assert(currentFlagName!=NULL);
	int i;
	int size = PMCommandLine_getValidFlagSize(cl);
	for(i=0;i<size;i++)
	{
		char FlagName[STRING_SIZE];
		PMCommandLine_getValidFlagName (cl, i, FlagName);
	
		if(!strcmp(currentFlagName, FlagName))
			return i;
	}
	return INVALID;
}

int PMCommandLine_getValidFlagType (PMCommandLine * cl, int index)
{
	assert(cl!=NULL);
	int size = PMCommandLine_getValidFlagSize(cl);
	assert(index>=0);
	assert(index<size);

	return cl->PMCommandLineValidFlagType[index];
}

int PMCommandLine_check (PMCommandLine * cl)
{
	assert(cl!=NULL);
	int i, size = PMCommandLine_getValidFlagSize(cl);
	for(i=0;i<size;i++)
	{
		int Type = PMCommandLine_getValidFlagType (cl, i);
		if(Type==REQUIRED)
		{
			char UserInput[STRING_SIZE];
			PMCommandLine_getUserFlagInput(cl, i, UserInput);

			if(!strcmp(UserInput, "\0"))
				return 0;
		}
	}

	return 1;	
}

int PMCommandLine_load (PMCommandLine * cl, int argc, char ** argv)
{
	int success = 0;
	int i;
	for(i=1;i<argc;i++)
	{
		char curFlagName[STRING_SIZE];
		strcpy(curFlagName, argv[i]);		
		removeFirstCharacter (curFlagName);
		
		int index = PMCommandLine_getValidFlagNameIndex (cl, curFlagName);
		if(index!=INVALID)
		{
			i++;
			if(i==argc)
			{
				ReportTerminalErrorMessage (stderr, "Incomplete command line input (missing argument after flag)");
				assert(i!=argc);
			}
			char curUserInput[STRING_SIZE];
			strcpy(curUserInput, argv[i]);
			PMCommandLine_setUserFlagInput(cl, index, curUserInput);
		}
		else
		{	
			ReportTerminalErrorMessage (stderr, "Invalid command line input flag"); 
			assert(index!=INVALID);
		}
	}
	success = PMCommandLine_check (cl);
	return success;
}

