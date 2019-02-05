#include "PM.h"

#ifndef _IPOPCNT
int bitcount(unsigned int n)
{
	int count = 0;

	while (n)
	{
		count += n & 0x1u;
		n >>= 1;
	}

	return count;
}

void init_POPCNT_LUT16(void)
{
	int i;
	for (i = 0; i < (0x1u << 16); i++)
		POPCNT_LUT16[i] = bitcount(i);
}
int pm_popcnt_u64(unsigned long long val)
{
	//assert(0);
	return POPCNT_LUT16[val & 0xffffu] + POPCNT_LUT16[(val >> 16) & 0xffffu] + POPCNT_LUT16[(val >> 32) & 0xffffu] + POPCNT_LUT16[(val >> 48) & 0xffffu];
}
#endif

unsigned long long rdtsc(void)
{
	unsigned a, d;

	__asm__ volatile("rdtsc"
					 : "=a"(a), "=d"(d));

	return ((unsigned long long)a) | (((unsigned long long)d) << 32);
}

void ReportTerminalErrorMessage(FILE *fp, char *text)
{
	assert(fp != NULL);
	assert(text != NULL);
	fprintf(fp, "\nERROR: %s\n\n", text);
	exit(0);
}

void removeFirstCharacter(char *string)
{
	assert(string != NULL);
	int length = strlen(string);
	char tstring[STRING_SIZE];
	strcpy(tstring, string);
	memcpy(string, &tstring[1], length);
}

int getFileFlagLineNumber(FILE *fp, char *FileFlag)
{
	assert(fp != NULL);
	assert(FileFlag != NULL);

	char tstring[STRING_SIZE];
	int Lines = INVALID;
	while (fscanf(fp, "%s", tstring) != EOF)
	{
		if (fgetc(fp) == '\n')
			Lines++;

		if (!strcmp(tstring, FileFlag))
			return Lines + 2;
	}

	return Lines;
}

int getMatchIndex(char *vector, char character)
{
	//assert(vector!=NULL);

	int i;
	for (i = 0; i < strlen(vector); i++)
		if (vector[i] == character)
			return i;

	assert(0);
	return -1;
}

void getTypeString(int Type, char *String)
{
	assert(Type == INPUT_TYPE_VCF);
	assert(String != NULL);

	if (Type == INPUT_TYPE_VCF)
		strcpy(String, "VCF");
}

void resetString(char *String, int Length)
{
	assert(String != NULL);
	assert(Length >= 1);
	int i;
	for (i = 0; i < Length; i++)
		String[i] = '\0';
}

int dataComparison(uint64_t *A, uint64_t *B, int length)
{
	assert(A != NULL);
	assert(B != NULL);
	assert(length >= 1);

	int i;
	for (i = 0; i != length; ++i)
		if (A[i] != B[i])
		{ //printf("BPS: %d\n", i);
			bpvec[i]++;
			return 1;
		}
	return 0;
}

int dataComparisonInversion(uint64_t *A, uint64_t *B, int length, int sampleSize)
{
	assert(A != NULL);
	assert(B != NULL);
	assert(length >= 1);

	int i;
	for (i = 0; i != length - 1; ++i)
	{
		if ((~A[i]) != B[i])
		{ //printf("BPR: %d\n", i);
			//bpvec[i]++;
			return 1;
		}
	}
	char temp = ~B[length - 1];
	int shiftLast = WORD_WIDTH - (sampleSize - (length - 1) * WORD_WIDTH);
	temp = temp << shiftLast;
	temp = temp >> shiftLast;
	if (temp != A[length - 1])
		return 1;

	return 0;
}

void stringMatch(char *string1, char *string2)
{
	fprintf(stdout, "string matching: ");
	int i, length1 = strlen(string1);
	int j, length2 = strlen(string2);

	int val_up, val_le, val_diag;

	assert(length1 + 1 < PATTERNS_MAX);
	assert(length2 + 1 < PATTERNS_MAX);
	int **mat = (int **)malloc(sizeof(int *) * (length1 + 1));
	for (i = 0; i < length1 + 1; i++)
	{
		mat[i] = (int *)malloc(sizeof(int) * (length2 + 1));
		for (j = 0; j < length2 + 1; j++)
			mat[i][j] = 0;
	}
	int max_val = -1;
	int max_i = -1;
	int max_j = -1;

	// Matrix
	for (i = 1; i < length1 + 1; i++)
	{
		//printf("\n%d:\t", i);
		for (j = 1; j < length2 + 1; j++)
		{
			val_up = mat[i - 1][j] + INDEL_PEN;
			val_le = mat[i][j - 1] + INDEL_PEN;

			if (string1[i - 1] == string2[j - 1])
				val_diag = mat[i - 1][j - 1] + MATCH_SCR;
			else
				val_diag = mat[i - 1][j - 1] + MISMATCH_PEN;

			mat[i][j] = MAX3(val_diag, val_up, val_le);
			//printf("%d ", mat[i][j]);
			if (mat[i][j] >= max_val)
			{
				max_val = mat[i][j];
				max_i = i;
				max_j = j;
				printf("NEW: max val %d max i %d max j %d\n", max_val, max_i, max_j);
			}
		}
	}

	printf("END: max val %d max i %d max j %d\n", max_val, max_i, max_j);
	// Trace back
	int cur_score = mat[max_i][max_j];
	int up_score = mat[max_i - 1][max_j];
	int diag_score = mat[max_i - 1][max_j - 1];
	int left_score = mat[max_i][max_j - 1];

	while (cur_score > 0)
	{
		if (MAX3(up_score, diag_score, left_score) == up_score)
		{
			max_i--;
		}
		else
		{
			if (MAX3(up_score, diag_score, left_score) == diag_score)
			{
				max_i--;
				max_j--;
			}
			else
			{
				max_j--;
			}
		}

		if (max_i == 0)
			break;

		if (max_j == 0)
			break;

		cur_score = mat[max_i][max_j];
		up_score = mat[max_i - 1][max_j];
		diag_score = mat[max_i - 1][max_j - 1];
		left_score = mat[max_i][max_j - 1];
	}

	printf("START: max val %d max i %d max j %d\n", max_val, max_i, max_j);

	for (i = 0; i < length1 + 1; i++)
	{
		free(mat[i]);
	}
	free(mat);
}

void experimentalRegionCount(int s1, int s2, int s3, int s4, double omega)
{
	int i;
	for (i = 0; i < slsize; i++)
	{
		if (sl1[i] == s1 && sl2[i] == s2 && sl3[i] == s3 && sl4[i] == s4)
		{
			cole[i]++;
			slomega[i] += omega;
			return;
		}
	}

	slsize++;

	sl1 = realloc(sl1, sizeof(int) * slsize);
	sl2 = realloc(sl2, sizeof(int) * slsize);
	sl3 = realloc(sl3, sizeof(int) * slsize);
	sl4 = realloc(sl4, sizeof(int) * slsize);
	cole = realloc(cole, sizeof(int) * slsize);
	slomega = realloc(slomega, sizeof(double) * slsize);

	sl1[slsize - 1] = s1;
	sl2[slsize - 1] = s2;
	sl3[slsize - 1] = s3;
	sl4[slsize - 1] = s4;
	cole[slsize - 1] = 1;
	slomega[slsize - 1] = omega;
}
