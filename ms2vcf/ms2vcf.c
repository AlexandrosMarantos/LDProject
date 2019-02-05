#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

int main()
{
	FILE * fp = fopen("ms.out", "r");
	char tstring [1000];



	fscanf(fp, "%s", tstring);
	fscanf(fp, "%s", tstring);
	int samples = atoi(tstring);
	printf("samples: %d\n", samples);
	while(strcmp(tstring, "segsites:"))
		fscanf(fp, "%s", tstring);

	fscanf(fp, "%s", tstring);
	
	int sites = atoi(tstring);
	printf("sites: %d\n", sites);	

	fscanf(fp, "%s", tstring);

	float * positions = (float*)malloc(sizeof(float)*sites);

	int i, j;
	for(i=0;i<sites;i++)
	{
		fscanf(fp, "%s", tstring);
		positions[i] = atof(tstring);
		if(positions[i]*1000000>=510700)
		{
		//	printf("the index is %d\n", i);
			//assert(0);
		}
		if(i==3342)
		{
			//printf("index: %f\n", positions[i]);
			//fflush(stdout);
			//assert(0);
		}
	}

	char ** matrix = (char**)malloc(sizeof(char*)*samples);
	
	char * tline = (char*)malloc(sizeof(char)*(sites+1));
	for(i=0;i<samples;i++)
	{
		matrix[i] = (char*)malloc(sizeof(char)*(sites+1));
		
		fscanf(fp, "%s", matrix[i]);
	}
	
	fclose(fp);



	fp = fopen ("vcf_from_ms.vcf", "w");

	fprintf(fp, "##fileformat=VCFv4.1\n");
	fprintf(fp, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t");

	for(i=0;i<samples;i=i+2)
	{
		fprintf(fp, "sample%d", i);

		if(i<samples-2)
			fprintf(fp, "\t");
		else
			fprintf(fp, "\n");
	}


	for(j=0;j<sites;j++)
	{
		fprintf(fp, "23\t%d\t.\tA\tC\t100\tPASS\t.\tGT\t", (int)j+1);//positions[j]);
		for(i=0;i<samples;i=i+2)
		{

			fprintf(fp, "%c|%c", matrix[i][j], matrix[i+1][j]);
			//fprintf(fp, "%d", rand()%2);

			if(i<samples-2)
				fprintf(fp, "\t");
			else
				fprintf(fp, "\n");

		}
	}

	fclose (fp);	


	return 0;
}
