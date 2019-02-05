
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <time.h>
#include <sys/time.h>
#include <limits.h>
#include <ctype.h>
#include <inttypes.h>
#include <math.h>

#ifdef _IPOPCNT
#include <immintrin.h>
#else
int bitcount(unsigned int n);
char POPCNT_LUT16 [0x1u << 16];
void init_POPCNT_LUT16 (void);
int pm_popcnt_u64 (unsigned long long val);
#endif

#define STRING_SIZE 1000
#define CMD_NAME_INDEX 0
#define PLOIDY_INDEX 2
#define LDTHRESHOLD_INDEX 3
#define PADDING 4
#define VALID 1
#define INVALID -1
#define SUCCESS 1
#define FAILED -1
#define TRUE 1
#define FALSE 0
#define OPTIONAL 0
#define REQUIRED 1
#define INPUT_TYPE_MST 0 // Transposed MS
#define INPUT_TYPE_FASTAT 1 // Transposed FASTA
#define INPUT_TYPE_VCF 2
#define PATTERNS_MAX 10000
#define VCF_FIELDS 8
#define WORD_WIDTH 64 // default type uint64_t
#define MININUM_SNPS_IN_LD_CLASS 5
#define MINIMUM_SNPS_IN_LD_REGION 7
#define MATCH_SCR 5
#define MISMATCH_PEN -1
#define INDEL_PEN -1

#define MAX3(a, b, c) ((a) > (b) ? ((a) > (c) ? (a) : (c)) : ((b) > (c) ? (b) : (c)))

int * bpvec;

int REFERENCE_OP_INDEX;
		int slsize ;
		int * sl1 ;
		int * sl2 ;
		int * sl3 ;
		int * sl4 ;
		int * cole ;
		double * slomega ;
void experimentalRegionCount (int s1, int s2, int s3, int s4, double omega);

FILE * fpInput; 

// PM_SupportFunctions
unsigned long long rdtsc(void);
void ReportTerminalErrorMessage (FILE * fp, char * text);
void removeFirstCharacter (char * string);
int getFileFlagLineNumber (FILE * fp, char * FileFlag);
void getTypeString (int Type, char * String);
void resetString (char * String, int Length);
int getMatchIndex (char * vector, char character);
int dataComparison (uint64_t * A, uint64_t * B, int length);
int dataComparisonInversion (uint64_t * A, uint64_t * B, int length, int sampleSize);
void stringMatch (char * string1, char * string2);

// PM_CommandLine
typedef struct
{
	int	PMCommnadLineValidFlagSize;
	char **	PMCommandLineValidFlagName;
	int   * PMCommandLineValidFlagType;
	char ** PMCommandLineUserFlagInput;

} PMCommandLine;

PMCommandLine * PMCommandLine_g;

int PMCommandLine_getValidFlagSize (PMCommandLine * cl);
void PMCommandLine_setValidFlagSize (PMCommandLine * cl, int Size);
void PMCommandLine_appendValidFlagName (PMCommandLine * cl, char * FlagName, int Type);
void PMCommandLine_create (PMCommandLine * cl);
PMCommandLine * PMCommandLine_new (void);
void PMCommandLine_free (PMCommandLine * cl);
int PMCommandLine_load (PMCommandLine * cl, int argc, char ** argv);
void PMCommandLine_print (PMCommandLine * cl, FILE * fp);
void PMCommandLine_getValidFlagName (PMCommandLine * cl, int index, char * FlagName);
void PMCommandLine_getUserFlagInput (PMCommandLine * cl, int index, char * FlagInput);
int  PMCommandLine_getValidFlagNameIndex (PMCommandLine * cl, char * currentFlagName);
void PMCommandLine_setUserFlagInput (PMCommandLine * cl, int index, char * FlagInput);
int PMCommandLine_getValidFlagType (PMCommandLine * cl, int index);
int PMCommandLine_check (PMCommandLine * cl);

// PM_Input
typedef struct
{
	int 	PMInputType;
	int	PMInputSampleSize;
	char ** PMInputSampleName;
	char	PMInputChromosome[STRING_SIZE];
	int	PMInputFirstDataLineNumber;	

} PMInput;

PMInput * PMInput_g;

PMInput * PMInput_new (void);
void PMInput_free (PMInput * i);
void PMInput_setType (PMInput * i, int Type);
int PMInput_getType (PMInput * i);
int PMInput_getSampleSize (PMInput * i);
void PMInput_setSampleSize (PMInput * i, int SampleSize);
void PMInput_appendSampleName (PMInput * i, char * SampleName);
void PMInput_loadSampleSize (PMInput * i, FILE * fp);
FILE * PMInput_check (PMInput * i, PMCommandLine * cl);
void PMInput_print (PMInput * i, FILE * fp);
void PMInput_getChromosome (PMInput * i, char * Chromosome);
void PMInput_setFirstDataLineNumber (PMInput * input, int LineNumber);
int PMInput_getFirstDataLineNumber (PMInput * input);

// PM_Window
typedef struct
{
	int	PMWindowSites;
	int	PMWindowSNPs;
	int	PMWindowDiscarded;
	int	PMWindowLeftmostSitePosition;
	int	PMWindowRightmostSitePosition;
	int *	PMWindowPatternMap;
	int *   PMWindowPositionMap;

} PMWindow;

PMWindow * PMWindow_new (void);
void PMWindow_free (PMWindow * w);
void PMWindow_setSites (PMWindow * w, int Sites);
int PMWindow_getSites (PMWindow * w);
void PMWindow_incrementSites (PMWindow * w);
void PMWindow_setSNPs (PMWindow * w, int SNPs);
int PMWindow_getSNPs (PMWindow * w);
void PMWindow_incrementSNPs (PMWindow * w);
void PMWindow_setDiscarded (PMWindow * w, int Discarded);
int PMWindow_getDiscarded (PMWindow * w);
void PMWindow_incrementDiscarded (PMWindow * w);
void PMWindow_setLeftmostSitePosition (PMWindow * w, int Position);
void PMWindow_setRightmostSitePosition (PMWindow * w, int Position);
int PMWindow_getLeftmostSitePosition (PMWindow * w);
int PMWindow_getRightmostSitePosition (PMWindow * w);
void PMWindow_resetLeftmostRightmostSitePosition (PMWindow * w);
void PMWindow_appendPattern (PMWindow * w, int PatternID, int SNPPosition);
void PMWindow_appendNoSNP (PMWindow * w, int PatternID, int SNPPosition);

typedef struct
{
	int PMLDSubregionLength;
	int * PMLDSubregionPatternIDWindow;

} PMLDSubregion;

// PM_LDCLass
typedef struct
{
	int 	PMLDClassSNPSize; // Total number of SNPs
	int 	PMLDClassPatternSize; // Total number of Patterns
	int * 	PMLDClassPatternList; // List of Pattern IDs for the different patterns in the class
	int	PMLDClassFirstSNPIndex;
	int	PMLDClassLastSNPIndex;
	int	PMLDClassFirstSNPPosition;
	int	PMLDClassLastSNPPosition;
	char *	PMLDClassQueryString;
	char *	PMLDClassRegionView; // This is how the region of SNPs looks, represented as a string, when processing the current LD class

	int	PMLDClassMatchListSize; // Total number of max vals detected during dp
	int *	PMLDClassMatchListScore;
	int *   extra;
	int *	PMLDClassMatchListCoordinateI;
	int *	PMLDClassMatchListCoordinateJ;
	int *   PMLDClassMatchListCoordinateJ2;


	int	PMLDClassValidSubregions;
	int *	PMLDClassValidSubregionFirstSNP;
	int *   PMLDClassValidSubregionLastSNP;
	float * PMLDClassValidSubregionTotalLD;
	float * PMLDClassValidSubregionAvgLD;
	int *   PMLDClassValidSubregionUsedTimes;

	PMLDSubregion ** PMLDClassSubregionList;

	int 	PMLDClassUsed;

} PMLDClass;

PMLDClass * PMLDClass_new (void);
void PMLDClass_increaseSNPSize (PMLDClass * ldc, int size);
int PMLDClass_getSNPSize (PMLDClass * ldc);
void PMLDClass_free (PMLDClass * ldc);
void PMLDClass_print(PMLDClass * ldc, FILE * fp, int index, void * pp);
void PMLDClass_appendPattern (PMLDClass * ldc, int PatternID);
void PMLDClass_setFirstSNPIndex (PMLDClass * ldc, void * pp, int patternIndex);
void PMLDClass_setLastSNPIndex (PMLDClass * ldc, void * pp, int patternIndex);
void PMLDClass_createQueryStringAndRegionView(PMLDClass * ldc, void * pp);
void PMLDClass_matchQuery2Region (PMLDClass * ldc, PMWindow * w, void * pp);
void PMLDClass_processSubregions (PMLDClass * ldc, PMWindow * w, void * pp);

// PM_PatternPool
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

PMPatternPool * PMPatternPool_g;

PMPatternPool * PMPatternPool_new (void);
void PMPatternPool_free (PMPatternPool * pp);
int PMPatternPool_allocate (PMPatternPool * pp, PMInput * i, PMCommandLine * cl);
void PMPatternPool_reset (PMPatternPool * pp);
void PMPatternPool_setPloidy (PMPatternPool * pp, PMCommandLine * cl);
int PMPatternPool_getPloidy (PMPatternPool * pp);
void PMPatternPool_setLDThreshold (PMPatternPool * pp, PMCommandLine * cl);
float PMPatternPool_getLDThreshold (PMPatternPool * pp);
uint64_t * PMPatternPool_getPatternCompact (PMPatternPool * pp, int Index);
char * PMPatternPool_getPendingSNPPtr (PMPatternPool * pp);
uint64_t * PMPatternPool_getPendingSNPCompactPtr (PMPatternPool * pp);
int * PMPatternPool_getCount (PMPatternPool * pp);
void PMPatternPool_setSampleSize (PMPatternPool * pp, int SampleSize);
int PMPatternPool_getSampleSize (PMPatternPool * pp);
void PMPatternPool_setSize (PMPatternPool * pp, int Size);
int PMPatternPool_getSize (PMPatternPool * pp);
int PMPatternPool_getCountByIndex (PMPatternPool * pp, int Index);
void PMPatternPool_appendPosition (PMPatternPool * pp, int index, int Position);
void PMPatternPool_appendID (PMPatternPool * pp, int index, char * ID);
void PMPatternPool_appendSNPIndex (PMPatternPool * pp, int index, int SNPIndex);
void PMPatternPool_setDerivedAlleleCount (PMPatternPool * pp, int index, int * DerivedAlleleCount);
int PMPatternPool_getCompactSize (PMPatternPool * pp);
int PMPatternPool_matchSNP (PMPatternPool * pp, int Position, char * ID, int * DerivedAlleleCount, PMWindow * w, int SNP_AAindex);
PMWindow * PMPatternPool_getWindow (PMPatternPool * pp);
int PMPatternPool_loadRegionVCF (PMPatternPool * pp, PMInput * input, FILE * fp);
int PMPatternPool_loadRegion (PMPatternPool * pp, PMInput * i, FILE * fp);
int PMPatternPool_getMaxCount (PMPatternPool * pp);
int PMPatternPool_getDominantPattern (PMPatternPool * pp);
void PMPatternPool_print (PMPatternPool * pp, FILE * fp);
void PMPatternPool_computeLD (PMPatternPool * pp); 
void PMPatternPool_groupPatternsInLD (PMPatternPool * pp, FILE * fpclassoutput);
void PMPatternPool_appendLDClass(PMPatternPool * pp, PMLDClass * ldc);

// PM_ExecutionDetails
typedef struct 
{
	// All
	unsigned long long int 	PMExecutionDetailsExecutionStart;
	unsigned long long int 	PMExecutionDetailsExecutionEnd;
	float			PMExecutionDetailsExecutionTime;

	// Unique Patterns
	unsigned long long int	PMExecutionDetailsPatternsStart;
	unsigned long long int	PMExecutionDetailsPatternsEnd;
	float			PMExecutionDetailsPatternsTime;

	// LD Matrix
	unsigned long long int	PMExecutionDetailsLDMatrixStart;
	unsigned long long int	PMExecutionDetailsLDMatrixEnd;
	float			PMExecutionDetailsLDMatrixTime;

} PMExecutionDetails;

PMExecutionDetails * PMExecutionDetails_g;

PMExecutionDetails * PMExecutionDetails_new (void);
void PMExecutionDetails_free (PMExecutionDetails * ed);
void PMExecutionDetails_reset (PMExecutionDetails * ed);
void PMExecutionDetails_startExecutionTimer (PMExecutionDetails * ed);
void PMExecutionDetails_pauseExecutionTimer (PMExecutionDetails * ed);
float PMExecutionDetails_readExecutionTimer (PMExecutionDetails * ed);
void PMExecutionDetails_startPatternsTimer (PMExecutionDetails * ed);
void PMExecutionDetails_pausePatternsTimer (PMExecutionDetails * ed);
float PMExecutionDetails_readPatternsTimer (PMExecutionDetails * ed);
void PMExecutionDetails_startLDMatrixTimer (PMExecutionDetails * ed);
void PMExecutionDetails_pauseLDMatrixTimer (PMExecutionDetails * ed);
float PMExecutionDetails_readLDMatrixTimer (PMExecutionDetails * ed);
void PMExecutionDetails_print (PMExecutionDetails * ed, FILE * fp);

// PM_VCFParser
int PMVCFParser_parseFields (FILE * fp, char * Chromosome, int * AlleleSize, char * AlleleVector, int * GenotypeLocation, int * Position, char * ID);
void PMVCFParser_skipLine (FILE * fp);
int PMVCFParser_getGenotypeLocation (char * string);
int PMVCFParser_parseSamples (FILE * fp, PMPatternPool * pp, int AlleleSize, char * AlleleVector, int GenotypeLocation, int Position, int * DerivedAlleleCount);

// PM_LinkageDisequilibrium
float LD_ref (PMPatternPool * pp, int index1, int index2);
float LD_test1 (PMPatternPool * pp, int index1, int index2);
