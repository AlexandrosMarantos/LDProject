/*
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

*/

#include "PM.h"

PMWindow * PMWindow_new (void)
{
	PMWindow * w = NULL;
	w = (PMWindow*)malloc(sizeof(PMWindow));
	assert(w!=NULL);

	w->PMWindowSites = 0;
	w->PMWindowSNPs = 0;
	w->PMWindowDiscarded = 0;
	w->PMWindowLeftmostSitePosition = -1;
	w->PMWindowRightmostSitePosition =-1;
	w->PMWindowPatternMap = NULL;
	w->PMWindowPositionMap = NULL;

	return w;
}

void PMWindow_appendPattern (PMWindow * w, int PatternID, int SNPPosition)
{
	assert(w!=NULL);
	assert(PatternID>=0);
	assert(SNPPosition>=0);
	//printf("size: %d\n", w->PMWindowSNPs);
	w->PMWindowPatternMap = realloc(w->PMWindowPatternMap, sizeof(int)*w->PMWindowSites); // w->PMWindowSNPs
	assert(w->PMWindowPatternMap!=NULL);
	w->PMWindowPatternMap[w->PMWindowSites-1] = PatternID;

	w->PMWindowPositionMap = realloc(w->PMWindowPositionMap, sizeof(int)*w->PMWindowSites);
	assert(w->PMWindowPositionMap!=NULL);
	w->PMWindowPositionMap[w->PMWindowSites-1] = SNPPosition;
}

void PMWindow_appendNoSNP (PMWindow * w, int PatternID, int SNPPosition)
{
	assert(w!=NULL);
	assert(PatternID==-1);
	assert(SNPPosition>=0);
	//printf("size: %d\n", w->PMWindowSNPs);
	w->PMWindowPatternMap = realloc(w->PMWindowPatternMap, sizeof(int)*w->PMWindowSites);
	assert(w->PMWindowPatternMap!=NULL);
	w->PMWindowPatternMap[w->PMWindowSites-1] = PatternID;

	w->PMWindowPositionMap = realloc(w->PMWindowPositionMap, sizeof(int)*w->PMWindowSites);
	assert(w->PMWindowPositionMap!=NULL);
	w->PMWindowPositionMap[w->PMWindowSites-1] = SNPPosition;
}


void PMWindow_free (PMWindow * w)
{
	assert(w!=NULL);

	if(w->PMWindowPatternMap!=NULL)
		free(w->PMWindowPatternMap);

	if(w->PMWindowPositionMap!=NULL)
		free(w->PMWindowPositionMap);

	free(w);
}

void PMWindow_setLeftmostSitePosition (PMWindow * w, int Position)
{
	assert(w!=NULL);
	assert(Position>=1);

	if(w->PMWindowLeftmostSitePosition==INVALID)
		w->PMWindowLeftmostSitePosition = Position;
}

void PMWindow_setRightmostSitePosition (PMWindow * w, int Position)
{
	assert(w!=NULL);
	assert(Position>=1);
	assert(Position>=w->PMWindowLeftmostSitePosition);
	assert(w->PMWindowLeftmostSitePosition!=INVALID);

	w->PMWindowRightmostSitePosition = Position;
}

int PMWindow_getLeftmostSitePosition (PMWindow * w)
{
	assert(w!=NULL);

	return w->PMWindowLeftmostSitePosition;
}

int PMWindow_getRightmostSitePosition (PMWindow * w)
{
	assert(w!=NULL);

	return w->PMWindowRightmostSitePosition;
}

void PMWindow_resetLeftmostRightmostSitePosition (PMWindow * w)
{
	assert(w!=NULL);
	w->PMWindowLeftmostSitePosition = INVALID;
	w->PMWindowRightmostSitePosition = INVALID;
}

void PMWindow_setSites (PMWindow * w, int Sites)
{
	assert(w!=NULL);
	assert(Sites>=0);

	w->PMWindowSites = Sites;
}

int PMWindow_getSites (PMWindow * w)
{
	assert(w!=NULL);

	return w->PMWindowSites;
}

void PMWindow_incrementSites (PMWindow * w)
{
	assert(w!=NULL);

	w->PMWindowSites++;
}

void PMWindow_setSNPs (PMWindow * w, int SNPs)
{
	assert(w!=NULL);
	assert(SNPs>=0);

	w->PMWindowSNPs = SNPs;
}

int PMWindow_getSNPs (PMWindow * w)
{
	assert(w!=NULL);

	return w->PMWindowSNPs;
}

void PMWindow_incrementSNPs (PMWindow * w)
{
	assert(w!=NULL);

	w->PMWindowSNPs++;
}

void PMWindow_setDiscarded (PMWindow * w, int Discarded)
{
	assert(w!=NULL);
	assert(Discarded>=0);

	w->PMWindowDiscarded = Discarded;
}

int PMWindow_getDiscarded (PMWindow * w)
{
	assert(w!=NULL);

	return w->PMWindowDiscarded;
}

void PMWindow_incrementDiscarded (PMWindow * w)
{
	assert(w!=NULL);

	w->PMWindowDiscarded++;
}
