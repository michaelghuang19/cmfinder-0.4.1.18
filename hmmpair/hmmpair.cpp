#include "hmmpair.h"

#define CMFINDER_PACKAGE_VERSION "0.4.something"

void HmmForwardScore_SS (char *stoFileName,double nonCanonPairThreshold,int numFlankingPolyN,bool uniformDistributionOfProfileHmmStartsAndEnds,const char *positionsToIgnoreFileName,bool fragmentaryForCovaryingPairs,bool fragmentaryForHmmTraining);

int HitList::GetOverallLast (void) const
{
	if (empty()) {
		return 0;
	}
	else {
		return back().second;
	}
}
void HitList::Init (int first,int second)
{
	clear();
	std::pair<int,int> fullWindow;
	fullWindow.first=first;
	fullWindow.second=second;
	push_back(fullWindow);
}
void HitList::Init (int length)
{
	Init(0,length);
}
void HitList::InitEmpty (void)
{
	clear();
}
void AddWindowToList(HitList& hitList,const std::pair<int,int>& thisWindow)
{
	assert(thisWindow.second>=thisWindow.first);

	if (hitList.empty()) {
		hitList.push_back(thisWindow);
	}
	else {
		std::pair<int,int>& prev=hitList.back();
		if (prev.second >= thisWindow.first) {
			// merge them
			assert(prev.first <= thisWindow.first); // else we weren't calculating the left extent correctly, or it's not sliding properly
			assert(thisWindow.second >= prev.second);
			prev.second=thisWindow.second;
		}
		else {
			// this is a new one
			hitList.push_back(thisWindow);
		}
	}
}
float ToProb (float count,double totalProb,int numEvents)
{
	float prob=(float)(count/totalProb);
	if (count==0 && totalProb==0) {
		prob=(float)(1)/(float)(numEvents);
	}
	return prob;
}
float ToProb (NoUnderflowDouble count,NoUnderflowDouble totalProb)
{
	if (count==NoUnderflowDouble(0) && totalProb==NoUnderflowDouble(0)) {
		return 1.0;
	}
	else {
		NoUnderflowDouble prob=count/totalProb;
		return (float)(prob.ToDouble_ZeroOnUnderflow());
	}
}
float ToProb (double count,NoUnderflowDouble totalProb)
{
	return ToProb(NoUnderflowDouble(count),totalProb);
}
float ToProb (float count,NoUnderflowDouble totalProb)
{
	return ToProb(NoUnderflowDouble(count),totalProb);
}
float ToProb (float count,double totalProb)
{
	return ToProb(count,totalProb,1); // dunno how many events
}
float ToProb (double count,double totalProb)
{
	return ToProb((float)count,totalProb);
}
float ToProb (double count,double totalProb,int numEvents)
{
	return ToProb((float)count,totalProb,numEvents);
}

int try_main (int argc,char *argv[])
{
	if (argc<7) {
          if (argc>=1) {
            if (strcmp(argv[1],"--version")==0) {
              printf("CMFINDER_PACKAGE_VERSION=%s.\n",CMFINDER_PACKAGE_VERSION);
              exit(0);
            }
          }

          fprintf(stderr,"insufficient number of params (argc=%d)\n",argc);
          fprintf(stderr,"usage: hmmpair <sto file> <max non-canonical pair freq> <fragmentary-policy> <size of poly-N flanking seqs in forward prob calculations> <uniformDistributionOfProfileHmmStartsAndEnds> <partition function data file (from Vienna output), or the string \"NULL\" to ignore> <positions to ignore file name or the string \"NULL\">\n");
          fprintf(stderr,"positionsToIgnoreFileName is in emblcsv format (each line is comma-separated, and the fields are seqId,start,end).  it's for masking transcription terminators, and it's assumed to be double-sided, i.e. we should ignore a nucleotide if there's a terminator overlapping it, or a rev-comp of a terminator overlapping it.  So it doesn't matter if start<end or start>end.  Each line should have a terminator stem (the lowest coordinate of a paired nucleotide and the highest coordinate).  Also pairs within this range are ignored.\n");
          fprintf(stderr,"<fragmentary-policy>=d means disable, <fragmentary-policy>=f means force fragmentary for calculation of covarying base pairs only, <fragmentary-policy>=F means force fragmentary for covarying pairs and do the simple hack for training the HMM\n");
          fprintf(stderr,"reasonable defaults are:\nhmmpair <sto file> 0.05 f 200 0 <partition-func-stuff> <NULL or terminator stem predictions if you have them>\n");
          fprintf(stderr,"note: hmmpair only considers the secondary structure in the #=GC SS_cons line.  Motif predictions that include pseudoknots might benefit from an extended version of hmmpair that considers these additional base pairs.\n");
          return 1;
	}
	int a=1;
	char *stoFileName=argv[a++];
	double nonCanonPairThreshold=atof(argv[a++]);
	std::string fragmentaryPolicy=argv[a++];
	int numFlankingPolyN=atoi(argv[a++]);
	bool uniformDistributionOfProfileHmmStartsAndEnds=atoi(argv[a++])!=0;
	const char *positionsToIgnoreFileName=argv[a++];
	bool fragmentaryForCovaryingPairs=false;
	bool fragmentaryForHmmTraining=false;
	bool fragmentaryPolicyOkay=false;
	if (fragmentaryPolicy=="d") {
		fragmentaryPolicyOkay=true;
	}
	if (fragmentaryPolicy=="f") {
		fragmentaryPolicyOkay=true;
		fragmentaryForCovaryingPairs=true;
	}
	if (fragmentaryPolicy=="F") {
		fragmentaryPolicyOkay=true;
		fragmentaryForCovaryingPairs=true;
		fragmentaryForHmmTraining=true;
	}
	if (!fragmentaryPolicyOkay) {
		throw SimpleStringException("bad fragmentary-policy string (\"%s\")",fragmentaryPolicy.c_str());
	}
	HmmForwardScore_SS (stoFileName,nonCanonPairThreshold,numFlankingPolyN,uniformDistributionOfProfileHmmStartsAndEnds,positionsToIgnoreFileName,fragmentaryForCovaryingPairs,fragmentaryForHmmTraining);
	return 0;
}

int main (int argc,char *argv[])
{
	try {
		return try_main(argc,argv);
	}
	catch (const std::exception& e) {
		printf("ERROR: %s\n",e.what());
		return 1;
	}
}
