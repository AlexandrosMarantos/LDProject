#include "PM.h"

void PMInit (void)
{
	PMExecutionDetails_g = PMExecutionDetails_new();
	
	PMCommandLine_g = PMCommandLine_new();
	PMInput_g = PMInput_new();
	PMPatternPool_g = PMPatternPool_new();

	srand(0); // TODO: User input seed

#ifndef _IPOPCNT	
	init_POPCNT_LUT16 ();
#endif

}

void PMClose (void)
{
	PMExecutionDetails_free(PMExecutionDetails_g);

	if(fpInput!=NULL)
		fclose(fpInput);

	PMCommandLine_free(PMCommandLine_g);
	PMInput_free(PMInput_g);
	PMPatternPool_free(PMPatternPool_g);
}

int main (int argc, char** argv)
{
	PMInit();
	
	PMExecutionDetails * ed = PMExecutionDetails_g;
	PMExecutionDetails_startExecutionTimer(ed);

	PMCommandLine * cl = PMCommandLine_g;
	int success =  PMCommandLine_load(cl, argc, argv);
	if(success!=SUCCESS)
	{
		ReportTerminalErrorMessage (stderr, "Incomplete command line input (missing required argument)");
		return 0;
	}

	char outputFilePath [STRING_SIZE];
	PMCommandLine_getUserFlagInput (cl, CMD_NAME_INDEX, outputFilePath);
	strcat(outputFilePath, "_report.txt");
	FILE * fpOutput = fopen(outputFilePath, "w");
	assert(fpOutput!=NULL);



	char outputFilePath2 [STRING_SIZE];
	PMCommandLine_getUserFlagInput (cl, CMD_NAME_INDEX, outputFilePath2);
	strcat(outputFilePath2, "_class_report.txt");
	FILE * fpOutput2 = fopen(outputFilePath2, "w");
	assert(fpOutput2!=NULL);

	PMCommandLine_print(cl, stdout);

	if(success)
	{
		PMInput * i = PMInput_g;
		fpInput = PMInput_check (i, cl);
		PMInput_print (i, stdout);

		if(fpInput!=NULL)
		{
			PMPatternPool * pp = PMPatternPool_g;

			success = PMPatternPool_allocate (pp, i, cl);
			assert(success==SUCCESS);

			// TODO: Repeat for all windows here
			{
				bpvec = calloc(sizeof(int), 100);
				PMPatternPool_loadRegion (pp, i, fpInput);
				PMPatternPool_print (pp, fpOutput);
				PMPatternPool_computeLD (pp);
				PMPatternPool_groupPatternsInLD (pp, fpOutput2);
				//printf("93-93: %f 93-98 %f 98-93 %f\n", pp->PMPatternPoolLDMatrix[93][93],pp->PMPatternPoolLDMatrix[93][98],pp->PMPatternPoolLDMatrix[98][93]);
				//int bpi;
				//for(bpi=0;bpi<100;bpi++)
				//	printf("%d: %d\n", bpi, bpvec[bpi]);
			}				
		}
	}

	PMExecutionDetails_pauseExecutionTimer(ed);
	PMExecutionDetails_print(ed, stdout);

	fclose(fpOutput);
	fclose(fpOutput2);
	PMClose();

	return 0;
}
