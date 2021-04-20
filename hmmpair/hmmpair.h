#pragma GCC diagnostic ignored "-Wwrite-strings"
#pragma GCC diagnostic ignored "-Wformat"
#pragma GCC diagnostic ignored "-Wchar-subscripts"

extern "C" {
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#ifndef __APPLE__
#include <malloc.h>
#endif
#include <float.h>

#include "squid.h"		/* general sequence analysis library    */
#include "msa.h"                /* squid's multiple alignment i/o       */
#include "stopwatch.h"          /* squid's process timing module        */

#include "structs.h"		/* data structures, macros, #define's   */
#include "funcs.h"		/* external functions                   */
//#include "version.h"

// my added includes
#include "assert.h"
}

#ifndef _MSC_VER
#include <sys/types.h>
#include <unistd.h>
#endif

#include <vector>
#include <string>
#include <list>
#include <algorithm>
#include <exception>
#include <set>
#include <map>

#include <MiscExceptions.h>
#include <vectorPlus.h>
#include <multiDimVector.h>
#include <CommaSepFileReader.h>
#include <NoUnderflowDouble.h>

extern char *nucs;

#undef log2 // gcc 3.4 defines this
inline double log2 (double x) {
	return log(x)/log(2.0);
}
inline double pow2 (double x) {
	return pow(2.0,x);
}
#define LOGPOW2DEFINED
class HitList : public std::list<std::pair<int,int> > {
public:
	void Init (int length); // one interval from [0,length)
	void Init (int first,int second); // one interval
	void InitEmpty (void);
	void Dump (FILE *file) const;
	int GetOverallLast (void) const;

	int TotalSize (void) const;
	int64_t SizeIn2D (int windowLen) const; // what part of the dynamic programming table must we look at; it's basically TotalSize() * windowLen, except that at the beginning of each interval in the hit list, we only have to worry about a triangular part of the dynamic programming table.
};
struct TopLevelMatch { // matches corresponding to starting at state=0 (start state)
	int windowLast,windowLen;
	float score;

	// redundant
	int windowFirst;

	bool operator < (const TopLevelMatch& t) const {
		// sort first by windowFirst
		if (windowFirst!=t.windowFirst) {
			return windowFirst<t.windowFirst;
		}
		// then by windowLen
		return windowLen<t.windowLen;
	}
};
struct CykscanStats {
	MultiplyArray3d<float> scores;
	vector<float> scoresPerWindowLast;
	std::string programParams;

	vector<float> hmmScoresPerWindowLast; // for things with 2nd struct whose scores shouldn't ever be higher than a pure HMM
	vector<vector<float> > fullHmmDynProgTable; // first dimension is windowLast, next is hmm state

	bool isValid;
	bool collectScores;
	bool collectAboveThreshold;
	bool collectScoresPerWindowLast; // for comparing with HMM, just collects the max score from Start node over all possible windowLen
	bool collectHmmScoresPerWindowLast;
	bool collectHmmFullDynProgTable;

	float bestScore;

	typedef std::list<TopLevelMatch> TopLevelMatchList;
	TopLevelMatchList topLevelMatchesAboveThreshold;
	double threshold;

	CykscanStats (void) {
		isValid=false;
		collectScores=false;
		collectAboveThreshold=false;
		collectScoresPerWindowLast=false;
		collectHmmScoresPerWindowLast=false;
		collectHmmFullDynProgTable=false;
	}
};
class ScanParams {
public:
/*
	bool nullgcmix;
	int nullgcmix_gc_count[GC_SEGMENTS];
	ScanParams () {
		nullgcmix=false;
	}
	 */
};
int FindRightPartner (const std::string& ss,int first);
float ToProb (float count,double totalProb);
float ToProb (double count,double totalProb);
float ToProb (float count,NoUnderflowDouble totalProb);
float ToProb (double count,NoUnderflowDouble totalProb);
float ToProb (NoUnderflowDouble count,NoUnderflowDouble totalProb);
float ToProb (float count,double totalProb,int numEvents);
float ToProb (double count,double totalProb,int numEvents);
void AddWindowToList(HitList& hitList,const std::pair<int,int>& thisWindow);

#include "PositionsToIgnore.h"
#include "CovarianceModel.h"
#include "InfernalHmm.h"
#include "HmmType1.h"
#include "ScanHmm.h"
#include "SymbolicMath.h"
#include "Optimize.h"
#include "Cm2HMM.h"
#include "stl_extra.h"
#include "NaryCounter.h"
#include "MarkovModelStats.h"
#include "Cm2HmmOptimize.h"
