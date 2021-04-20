#include "hmmpair.h"
#include <ctype.h>

extern "C" {
#include "prior.h"		/* mixture Dirichlet prior */
#include "cm_eweight.h"

#include <fold.h>
#include <part_func.h>
#include <fold_vars.h>
#include <PS_dot.h>
#include <utils.h>

struct plist *make_plist(int length, double pmin) {
  /* convert matrix of pair probs to plist */
  struct plist *pl;
  int i,j,k=0,maxl;
  maxl = 2*length;
  pl = (struct plist *)space(maxl*sizeof(struct plist));
  k=0;
  for (i=1; i<length; i++)
    for (j=i+1; j<=length; j++) {
      if (pr[iindx[i]-j]<pmin) continue;
      if (k>=maxl-1) {
	maxl *= 2;
	pl = (struct plist *)xrealloc(pl,maxl*sizeof(struct plist));
      }
      pl[k].i = i;
      pl[k].j = j;
      pl[k++].p = pr[iindx[i]-j];
    }
  pl[k].i=0;
  pl[k].j=0;
  pl[k++].p=0.;
  return pl;
}
}

// GSCConsensus.cpp
struct ValidSeqRegion { // for fragmentary sequences
	int first,last;
};
typedef vector<ValidSeqRegion> ValidSeqRegionVector;
typedef std::set<std::pair<int,int> > IntPairSet;
extern void GetNucPairForMsaPosition(int& leftNuc,int& rightNuc,bool& hasDegen,MSA *msa,int seqNumWithinMsa,int ssFirst,int ssLast);
extern std::string NormalizeSs (const char *ssRaw,int len);
extern void GSCWeightedConsensus_CountPairFreqs_NoAdjustFreq(IntPairSet& pairSet,double& doubleGapWeight,double& gapWeight,double& nonCanonWeight,double& canonWeight,MSA *msa,int ssFirst,int ssLast,const ValidSeqRegionVector& validSeqRegionVector);
extern bool isCanonPair[MAXABET][MAXABET];
typedef double NucCount[MAXABET+1];
extern void GSCWeightedFreq_OneColumn(bool& hasData,NucCount& count,MSA *msa,const ValidSeqRegionVector& validSeqRegionVector,int pos);

// ScanHMM.cpp
extern NoUnderflowDouble ForwardOfLiteralString (const HmmType1& hmm,const char *rnaSequence,int rnaSequenceLen,HmmType1::State state1,int windowLast1,HmmType1::State state2,int windowLast2);
extern NoUnderflowDouble ForwardOfLiteralString (const HmmType1& hmm,const char *rnaSequence,int rnaSequenceLen,HmmType1::State state1,int windowLast1,HmmType1::State state2,int windowLast2,DynProgTableScoreCollector_OwnMemory<NoUnderflowDouble>& dynProgTable);

// MLHeuristic.cpp
extern void BuildHmm_MaxLikeliPathCorrespondence(const CovarianceModel& cm_,bool doLocalAlignment,InfernalHmm& createdInfernalHmm,const InfernalHmm& startInfernalHmm,const ExtraCm2HmmInfo& extraCm2HmmInfo,bool applyInsertPropogationBug);

// NOTE: these calculate the probabilities global wrt both the database seq, and the query HMM.
// we adjust the first and last IL state probabilities to allow for very long flanking seqs, while still remaining probabilistic
// we do this adjustment by putting lots of flankingpolyN on left & right of the actual alignment

typedef NoUnderflowDouble HmmSsReal;
HmmSsReal CalcProbThatSeqIsEmittedByHmm(const HmmType1& hmm,const vector<char>& dsq)
{
	const char *rnaSequence=&(dsq[0]); // convert to const char *; note that 'dsq' is 0-based, since I made the sequence myself
	HmmSsReal p=ForwardOfLiteralString(hmm,rnaSequence,(int)(dsq.size()),-1,0,-1,0);
	return p;
}
HmmSsReal CalcProbThatPairIsAlignedByHmm(const HmmType1& hmm,
	const vector<char>& dsq,
	int pos1,int pos2,
	int state1,int state2)
{
	const char *rnaSequence=&(dsq[0]); // convert to const char *; note as above that 'dsq' is 0-based since I made it
	const int windowLast1=pos1+1;
	const int windowLast2=pos2+1;
	return ForwardOfLiteralString(hmm,rnaSequence,(int)(dsq.size()),state1,windowLast1,state2,windowLast2);
}

struct DegappedSeqFromMsa {
	vector<char> dsq;
	vector<int> msaPosToDegappedPos;
	vector<int> degappedPosToMsaPos;
	vector<int> msaPosToHmmState;
};
typedef vector<DegappedSeqFromMsa> DegappedSeqFromMsaVector;
void DegapSeqs_FindMap_AddPolyN(DegappedSeqFromMsaVector& degappedSeqFromMsaVector,MSA *msa,char **msa_dsq,int numFlankingPolyN)
{
	degappedSeqFromMsaVector.resize(msa->nseq);
	for (int i=0; i<msa->nseq; i++) {
		DegappedSeqFromMsa& s=degappedSeqFromMsaVector[i];
		s.dsq.reserve(msa->alen + 2*numFlankingPolyN);
		s.msaPosToDegappedPos.resize(msa->alen);
		s.degappedPosToMsaPos.resize(msa->alen+2*numFlankingPolyN); // bigger than it needs to be, but whatever
		s.dsq.insert(s.dsq.end(),numFlankingPolyN,SymbolIndex('N'));
		for (int p=0; p<msa->alen; p++) {
			char nuc=msa_dsq[i][p+1];  // +1 since sentinel
			if (nuc==DIGITAL_GAP) {
				s.msaPosToDegappedPos[p]=-1;
			}
			else {
				s.msaPosToDegappedPos[p]=(int)(s.dsq.size());
				s.degappedPosToMsaPos[(int)(s.dsq.size())]=p;
				s.dsq.push_back(nuc);
			}
		}
		s.dsq.insert(s.dsq.end(),numFlankingPolyN,SymbolIndex('N'));
	}
}
// adapted from ParsetreeCount in parsetree.c
void FragmentaryParsetreeCount(CM_t *cm, Parsetree_t *tr, char *dsq, float wgt,const ValidSeqRegion& validSeqRegion,vector<double>& nodeGotDataWeight)
{
	int tidx;			/* counter through positions in the parsetree        */
	int v,z;			/* parent, child state index in CM                   */

	/* trivial preorder traverse, since we're already numbered that way */
	if (tr->n==0) {
		return;
	}
	int first_emitr=tr->emitr[0];
	for (tidx = 0; tidx < tr->n; tidx++) {
		v = tr->state[tidx];        	/* index of parent state in CM */
		if (v != cm->M && cm->sttype[v] != E_st && cm->sttype[v] != B_st) 
		{
			z = tr->state[tr->nxtl[tidx]];      /* index of child state in CM  */

			int emitl=tr->emitl[tidx]-1; // make 0-based
			// remember since it's an HMM, emitr is fixed
			assertr(tr->emitr[tidx]==first_emitr);

			if (emitl>=validSeqRegion.first && emitl<validSeqRegion.last) {

				int ndidx=cm->ndidx[v];
				nodeGotDataWeight[ndidx]++;

				if (z == cm->M)                
					cm->end[v] += wgt;
				else if (v == 0 && z - cm->cfirst[v] >= cm->cnum[v])
					cm->begin[z] += wgt;
				else
					cm->t[v][z - cm->cfirst[v]] += wgt; 

				if (cm->sttype[v] == MP_st) 
					PairCount(cm->e[v], dsq[tr->emitl[tidx]], dsq[tr->emitr[tidx]], wgt);
				else if (cm->sttype[v] == ML_st || cm->sttype[v] == IL_st) 
					SingletCount(cm->e[v], dsq[tr->emitl[tidx]], wgt);
				else if (cm->sttype[v] == MR_st || cm->sttype[v] == IR_st) 
					SingletCount(cm->e[v], dsq[tr->emitr[tidx]], wgt);
			}
		}
	}
}    
void HmmForwardScore_SS_MakeHmm(CovarianceModel& get_cm,InfernalHmm& hmm,vector<int>& msaPosToHmmState,DegappedSeqFromMsaVector& degappedSeqFromMsaVector,int numFlankingPolyN,MSA *msa,bool addFlankingPolyNWhenTrainingCM,bool uniformDistributionOfProfileHmmStartsAndEnds,const ValidSeqRegionVector& validSeqRegionVector,bool fragmentaryForHmmTraining,const std::string& ss)
{
	// always use GSC weighting, but assume that they're already set

	// temporarily remove the SS from the MSA
	char *saved_ss=msa->ss_cons;
	std::string linearSs(msa->alen,'.');
	msa->ss_cons=(char *)(linearSs.c_str());

	// following code is from Infernal-0.7, cmbuild.c (and adapted from code in CmalignZasha.cpp)
	float gapthresh=(float)1.1; // so nothing from the user is an insert
	int use_rf=FALSE;
	char **dsq = DigitizeAlignment(msa->aseq, msa->nseq, msa->alen);
	CM_t *cm;
	Parsetree_t *mtr;
	HandModelmaker(msa, dsq, use_rf, gapthresh, &cm, &mtr);
	// don't re-balance tree: we're not really doing anything with this CM anyway, and I don't want to mess around with the data in 'mtr'

	Parsetree_t *fullSeqParsetree;
	vector<char> polyA;
	polyA.assign(msa->alen,'A');
	char *polyAdsq=DigitizeSequence(&(polyA.front()),(int)polyA.size());
	fullSeqParsetree=Transmogrify(cm, mtr, polyAdsq, &(polyA.front()), msa->alen);
	free(polyAdsq);

	vector<double> nodeGotDataWeight;
	nodeGotDataWeight.assign(cm->nodes,0.0);
	Parsetree_t **forest=new Parsetree_t *[msa->nseq];
	for (int idx = 0; idx < msa->nseq; idx++)
	{
		forest[idx] = Transmogrify(cm, mtr, dsq[idx], msa->aseq[idx], msa->alen);
		if (fragmentaryForHmmTraining) {
			//throw SimpleStringException("sorry, I didn't implement training the HMM with fragmentary, since it got kinda complicated.  BUT, I've just realized that HMMs don't have base pairs, so I don't have to worry about that, and that also restricts what kind of parsetrees I can get.");
			FragmentaryParsetreeCount(cm, forest[idx], dsq[idx], msa->wgt[idx],validSeqRegionVector[idx],nodeGotDataWeight);
		}
		else {
			ParsetreeCount(cm, forest[idx], dsq[idx], msa->wgt[idx]);
		}

		if (addFlankingPolyNWhenTrainingCM) {
			// simulate this by simply adding the counts to ParsetreeCount (see the function 'ParsetreeCount'; I've tried to give variables the same name where easy)
			float wgt=msa->wgt[idx]; // the weight for the current sequence
			int v;

			// for N insertions on left & right at beginning of seq, we'll have
			// the sequence: 1*(S->IL), N*(IL->IL), 1*(IL->IR), N*(IR->IR), 1*(IR->next), where 'next' is the next state to transition to
			// first, let's do the transitions we know about (ignoring the transition to 'next')
			// transitions from root, which is always state 0
			v=0;
			int IL_child_num=0; // left state; the IL is always state 1, and child 0
			cm->t[v][IL_child_num] += wgt; // 1*(S->IL)
			// transitions from IL
			v=1;
			cm->t[v][0] += wgt*numFlankingPolyN; // N*(IL->IL)
			cm->t[v][1] += wgt; // 1*(IL->IR)
			// transitions from IR
			v=2;
			cm->t[v][0] += wgt*numFlankingPolyN; // N*(IR->IR)

			// we should subtract out the first transition to 'next' in the actual parsetree,
			// but I doubt it'll make much of a difference, so I'm going to be lazy (z1101010199)
#if 0
			// now subtract out the first transition to 'next' in the actual parsetree
			if (forest[idx]->n > 0) {
				v = tr->state[tidx];
				assertr(v==0 && cm->sttype[v]==S_st); // first should be root, since we're not using --local here.  This allows us to simplify some code I'm copying from the first iteration of ParsetreeCount
				int z = tr->state[tr->nxtl[tidx]];      // 'z' in ParsetreeCount corresponds to 'next' in my comment above
				assertr(z!=cm->M); // shouldn't be local end, since we're not using --local
				assertr (!(v == 0 && z - cm->cfirst[v] >= cm->cnum[v])); // nor can it be local begin
				if (z==1) {
					// this was a transition into the IL state
					// just add another IL->IL state to account for it
					// subtract out the addition that ParsetreeCount made
					cm->t[v][z - cm->cfirst[v]] -= wgt; 
				}
			}
#endif
		}
	}

	if (fragmentaryForHmmTraining) {
		// extrapolate counts

		// I didn't implement this
		// see 'papers-talks/motifs_2007/notes-pipeline-lab-notebook.doc', dated 2012-02-08
		assertr(false);
	}

	// always do entropy weighting
	float eloss;		        /* target entropy loss, entropy-weights  */
	float etarget;		/* target entropy (background - eloss)   */
	eloss             = (float)0.54;
	etarget = (float)(2.0) - eloss;
	Prior_t *pri=NULL;                 /* mixture Dirichlet prior structure */
	pri = Prior_Default();
	float eff_nseq = (float)(CM_Eweight(cm, pri, (float) msa->nseq, etarget));
	CMRescale(cm, eff_nseq / (float) msa->nseq);
	float randomseq[MAXABET];     /* null sequence model                     */
	CMDefaultNullModel(randomseq);
	CMSetNullModel(cm, randomseq);
	PriorifyCM(cm, pri, NULL);
	CMLogoddsify(cm);
	float el_selfprob;
	el_selfprob     = (float)0.94;
	cm->el_selfsc = (float)(sreLOG2(el_selfprob));

	// (don't allow --local)
	CMLogoddsify(cm);
	CMHackInsertScores(cm);	/* "TEMPORARY" fix for bad priors */

	CovarianceModel cmSansStruct;
	cmSansStruct.MirrorFrom(cm);

	// convert to HMM
	// since I'm lazy, use the ML-heuristic, which'll do a pretty good job of preserving the probabilities, and the code is already written
	// this code is adapted from MakeCmHmmForOptimization and BuildHmm_MaxLikeliPathCorrespondence(cmFileName,...)
	bool doLocalAlignment=false;
	bool applyInsertPropogationBug=false;
	std::string programParams="(internal)";
	ExtraCm2HmmInfo extraInfoTemp;
	extraInfoTemp.actuallySolveScores=false;
	InfernalHmm tempHmm;
	char *cmFileName=(char *)"(internal)";
	Cm2Hmm_HmmBuildType hmmType=HmmBuildType_Original;
	Cm2Hmm_WithWeighting_NoCaching(tempHmm,hmmType,cmSansStruct,cmFileName,NULL,&extraInfoTemp);
	ExtraCm2HmmInfo extraCm2HmmInfo;
	extraCm2HmmInfo.actuallySolveScores=false;
	InfernalHmm dummyInfernalHmm;
	Cm2Hmm_WithWeighting_NoCaching (dummyInfernalHmm,hmmType,cmSansStruct,cmFileName,NULL,&extraCm2HmmInfo);
	BuildHmm_MaxLikeliPathCorrespondence(cmSansStruct,doLocalAlignment,hmm,tempHmm,extraCm2HmmInfo,applyInsertPropogationBug);
	hmm.MultiplyEmitsBy(0.25);
	hmm.HackInsertScoresToStrictProbs();

	//cmSansStruct.Save("data/t/test-score.cm");

	// work out where states went for each sequence.
	for (int i=0; i<msa->nseq; i++) {
		degappedSeqFromMsaVector[i].msaPosToHmmState.assign(msa->alen,-1);
		for (int tpos=0; tpos<forest[i]->n; tpos++) {
			CovarianceModel::State cmState  = CovarianceModel::IntToState(forest[i]->state[tpos]);
			assertr(cmState!=cmSansStruct.GetLastState()); // else this'd imply a local end, and we're not working with local CMs
			assertr(cmSansStruct.GetStateType(cmState)!=MR_st && cmSansStruct.GetStateType(cmState)!=MP_st); // that'd be odd, given that we're trying to make a linear one

			int oneBasedMsaPos=forest[i]->emitl[tpos];
			int msaPos=oneBasedMsaPos-1; // -1 since 'emitl' is in 1-based coords, while we want 0-based
			if (cmSansStruct.GetStateType(cmState)==ML_st) { // only map these
				InfernalHmm::State hmmLeftState=hmm.GetHmmLeftStateOfCmState(cmState);

				degappedSeqFromMsaVector[i].msaPosToHmmState[msaPos]=InfernalHmm::StateToInt(hmmLeftState);
			}
		}
	}
	// get stuff from actual msaPos independent of any gaps in individual seqs
	msaPosToHmmState.assign(msa->alen,-1);
	for (int tpos=0; tpos<fullSeqParsetree->n; tpos++) {
		CovarianceModel::State cmState  = CovarianceModel::IntToState(fullSeqParsetree->state[tpos]);
		assertr(cmState!=cmSansStruct.GetLastState()); // else this'd imply a local end, and we're not working with local CMs
		assertr(cmSansStruct.GetStateType(cmState)!=MR_st && cmSansStruct.GetStateType(cmState)!=MP_st); // that'd be odd, given that we're trying to make a linear one

		int oneBasedMsaPos=fullSeqParsetree->emitl[tpos];
		int msaPos=oneBasedMsaPos-1; // -1 since 'emitl' is in 1-based coords, while we want 0-based
		if (cmSansStruct.GetStateType(cmState)==ML_st) { // only map these
			InfernalHmm::State hmmLeftState=hmm.GetHmmLeftStateOfCmState(cmState);

			msaPosToHmmState[msaPos]=InfernalHmm::StateToInt(hmmLeftState);
		}
	}

	if (uniformDistributionOfProfileHmmStartsAndEnds) {
		InfernalHmm::State state,state2;
		int nuc;

		// set the probabilities of transitions I don't want to zero, since keeping the topology of the HMM the same helps to avoid annoying bugs

		// S -> IL is the only possibility, with probability 1.0
		state=InfernalHmm::IntToState(0);
		state2=InfernalHmm::IntToState(1);
		assertr(hmm.GetNthChildState(state,0)==state2);
		assertr(hmm.GetNumChildren(state)==3);
		hmm.SetTransitionLogScore(state,0,0);
		hmm.SetTransitionLogScore(state,1,(float)IMPOSSIBLE);
		hmm.SetTransitionLogScore(state,2,(float)IMPOSSIBLE);

		// IL -> IL has probability 1.0, and IL -> ML also has probability 1.0.  This doesn't sum to 1.0, but in principle we can divide by the total number of starts and ends (although there's no point in doing this)
		state=InfernalHmm::IntToState(1);
		state2=InfernalHmm::IntToState(1);
		assertr(hmm.GetNthChildState(state,0)==state2);
		assertr(hmm.GetNumChildren(state)==3);
		hmm.SetTransitionLogScore(state,0,0);
		hmm.SetTransitionLogScore(state,1,0);
		hmm.SetTransitionLogScore(state,2,(float)IMPOSSIBLE);
		// and emits also all have probability 1.0.
		for (nuc=0; nuc<Alphabet_size; nuc++) {
			hmm.SetSingletEmissionLogScore(state,nuc,0);
		}

		// IR -> IR with prob 1.0, and IR -> E with prob 1.0
		state2=hmm.GetActualLastState();
		state=state2;
		state--;
		assertr(hmm.GetNumChildren(state)==2);
		hmm.SetTransitionLogScore(state,0,0);
		hmm.SetTransitionLogScore(state,1,0);
		for (nuc=0; nuc<Alphabet_size; nuc++) {
			hmm.SetSingletEmissionLogScore(state,nuc,0);
		}

		// PASSTHRU -> IR with prob 1.0  (should be PASSTHRU since we had a linear/structureless CM)
		state2=state;
		state--;
		assertr(hmm.GetStateType(state)==PASSTHRU_st);
		assertr(hmm.GetNthChildState(state,0)==state2);
		assertr(hmm.GetNumChildren(state)==2);
		hmm.SetTransitionLogScore(state,0,0);
		hmm.SetTransitionLogScore(state,1,(float)IMPOSSIBLE);
	}

	get_cm.CopyFrom(cmSansStruct);

	// restore SS to MSA
	msa->ss_cons=saved_ss;

	// free memory
	FreeCM(cm);
	Free2DArray((void**)dsq, msa->nseq);
	FreeParsetree(mtr);
	for (int i=0; i<msa->nseq; i++) {
		FreeParsetree(forest[i]);
	}
	FreeParsetree(fullSeqParsetree);
	delete [] forest;
	Prior_Destroy(pri);
}
struct ProbCacheForBasePairElement {
	int ssFirst,ssLast; // used for sanity checking
	bool isValid;
	HmmSsReal p;
};
typedef vector<ProbCacheForBasePairElement> ProbCacheForBasePair;
void Init_ProbCacheForBasePair (ProbCacheForBasePair& probCacheForBasePair,int numSeqs)
{
	ProbCacheForBasePairElement e;
	e.isValid=false;
	probCacheForBasePair.assign(numSeqs,e);
}
class BestScore {
	int ssFirst;
	int i1,i2;
	bool isValid;
	HmmSsReal score;
public:
	BestScore () { isValid=true; score=0; }
	void ProcessNewScore (HmmSsReal score_,int ssFirst_,int i1_,int i2_) {
		if (!isValid || score_ > score) {
			isValid=true;
			score=score_;
			ssFirst=ssFirst_;
			i1=i1_;
			i2=i2_;
		}
	}
	void ProcessNewScore (const BestScore& t) {
		if (!t.isValid) {
			return;
		}
		if (!isValid || t.score > score) {
			isValid=true;
			score=t.score;
			ssFirst=t.ssFirst;
			i1=t.i1;
			i2=t.i2;
		}
	}
	HmmSsReal GetScore () const {
		assertr(isValid || score==NoUnderflowDouble(0));
		return score;
	}
	int GetSsFirst () const { return ssFirst; }
	int GetI1 () const { return i1; }
	int GetI2 () const { return i2; }
	void PrintForUser (FILE *out,const char *name) {
		printf("%s score=%lg,%lg,%d,%d,%d\n",name,score.ToDouble_ZeroOnUnderflow(),score.Log2(),ssFirst,i1,i2);
	}
};
struct HmmPairScores {
	BestScore bestPair,bestPairPartitionFunc;
	HmmSsReal pairSum,pairSumWithPartitionFunc;
};
HmmSsReal CalcProbThatPairIsAlignedByHmm(ProbCacheForBasePair& probCacheForBasePair,const HmmType1& hmm,const DegappedSeqFromMsaVector& degappedSeqFromMsaVector,int i,int ssFirst,int ssLast)
{
	if (probCacheForBasePair[i].isValid) {
		assertr(ssFirst==probCacheForBasePair[i].ssFirst && ssLast==probCacheForBasePair[i].ssLast); // otherwise the cache is not being cleared when we move to a new base pair
	}
	else {
		probCacheForBasePair[i].ssFirst=ssFirst;
		probCacheForBasePair[i].ssLast=ssLast;
		probCacheForBasePair[i].isValid=true;

		int leftNucDegappedPos=degappedSeqFromMsaVector[i].msaPosToDegappedPos[ssFirst];
		int rightNucDegappedPos=degappedSeqFromMsaVector[i].msaPosToDegappedPos[ssLast-1];
		int leftHmmState=degappedSeqFromMsaVector[i].msaPosToHmmState[ssFirst];
		int rightHmmState=degappedSeqFromMsaVector[i].msaPosToHmmState[ssLast-1];
		assertr(leftHmmState>=0 && rightHmmState>=0);
		probCacheForBasePair[i].p=CalcProbThatPairIsAlignedByHmm(hmm,
			degappedSeqFromMsaVector[i].dsq,
			leftNucDegappedPos,
			rightNucDegappedPos,
			leftHmmState,rightHmmState);
	}
	return probCacheForBasePair[i].p;
}
struct BasePair {
	int ssFirst,ssLast;
};
typedef std::list<BasePair> BasePairList;
class PartitionFunctionLookup {
	vector2d<double> lookup; // first dimension is sequence#, second dimension is msaSsFirst (i.e. the ssFirst in MSA coordinates)
	vector<int> ssLastBySsFirst;
  void PopulateLookup(MSA *msa,const DegappedSeqFromMsaVector& degappedSeqFromMsaVector,const BasePairList& basePairList,int numFlankingPolyN);
public:
	PartitionFunctionLookup (MSA *msa,const DegappedSeqFromMsaVector& degappedSeqFromMsaVector,const BasePairList& basePairList,int numFlankingPolyN);
	~PartitionFunctionLookup ();
	double GetBasePairProb (int seq,int msaSsFirst);
};
PartitionFunctionLookup::PartitionFunctionLookup (MSA *msa,const DegappedSeqFromMsaVector& degappedSeqFromMsaVector,const BasePairList& basePairList,int numFlankingPolyN)
{
	ssLastBySsFirst.assign(msa->alen,-1);
	for (BasePairList::const_iterator bpi=basePairList.begin(); bpi!=basePairList.end(); bpi++) {
		assertr(ssLastBySsFirst[bpi->ssFirst]==-1); // otherwise there are conflicting base pairs
		ssLastBySsFirst[bpi->ssFirst]=bpi->ssLast;
	}

	if (false) {
          // disable partition function stuff
		lookup.assign(msa->nseq,msa->alen,1.0);
		return;
	}
	lookup.assign(msa->nseq,msa->alen,0); // all base pairs are 0 prob unless explicitly mentioned

        PopulateLookup(msa,degappedSeqFromMsaVector,basePairList,numFlankingPolyN);
}
void PartitionFunctionLookup::PopulateLookup(MSA *msa,const DegappedSeqFromMsaVector& degappedSeqFromMsaVector,const BasePairList& basePairList,int numFlankingPolyN)
{
  char *string, *line;
  char *structure=NULL, *cstruc=NULL;
  char  *c;
  int   i, length, l, sym, r;
  double energy, min_en;
  double kT, sfact=1.07;
  int   pf=0, noPS=0, istty;
  int noconv=0;
  double tmp;
  bool d=false;

  int seqIndex,nucIndex;
  for (seqIndex=0; seqIndex<msa->nseq; seqIndex++) {
    int degappedSeqIndex=0;
    char *degappedSeq=(char *)malloc(msa->alen+1);
    for (nucIndex=0; nucIndex<msa->alen; nucIndex++) {
      char nuc=msa->aseq[seqIndex][nucIndex];
      if (isalpha(nuc)) {
        degappedSeq[degappedSeqIndex++]=toupper(nuc);
      }
    }
    degappedSeq[degappedSeqIndex]=0;;

    pf=1; /* always do partition func  -- actually that's all we want to do */
    do_backtrack = 1;

    string=degappedSeq;
    length = (int) strlen(string);
    structure = (char *) space((unsigned) length+1);
    for (l = 0; l < length; l++) {
      string[l] = toupper(string[l]);
      if (!noconv && string[l] == 'T') string[l] = 'U';
    }
    min_en = fold(string, structure);
    if (d) printf("seq & struct (seqIndex=%d)\n%s\n%s\n", seqIndex, string, structure);
    if (pf) {
      char *pf_struc;
      pf_struc = (char *) space((unsigned) length+1);
      if (dangles==1) {
        dangles=2;   /* recompute with dangles as in pf_fold() */
        min_en = energy_of_struct(string, structure);
        dangles=1;
      }
    
      kT = (temperature+273.15)*1.98717/1000.; /* in Kcal */
      pf_scale = exp(-(sfact*min_en)/kT/length);
      //if (length>2000) fprintf(stderr, "scaling factor %f\n", pf_scale);
    
      init_pf_fold(length);

      if (cstruc!=NULL) {
        strncpy(pf_struc, cstruc, length+1);
      }
      energy = pf_fold(string, pf_struc);

      if (do_backtrack) {
        plist *pl1,*pl2,*plt;
        char *cent;
        pl1 = make_plist(length, 1e-5);
        /*(void) PS_dot_plot_list(string, ffname, pl1, pl2, "");*/
        for (plt=pl1; plt->i>0; plt++) {
          int i=plt->i;
          int j=plt->j;
          double p=plt->p;

          int ssFirst=i-1;
          int ssLast=j-1+1; // +1 to make it half-open interval
          int msaSsFirst=degappedSeqFromMsaVector[seqIndex].degappedPosToMsaPos[ssFirst+numFlankingPolyN];
          int msaSsLast=degappedSeqFromMsaVector[seqIndex].degappedPosToMsaPos[ssLast+numFlankingPolyN];
          if (ssLastBySsFirst[msaSsFirst]==msaSsLast) {
            assertr(lookup[seqIndex][msaSsFirst]==0); // otherwise we have conflicting base pairs for ssFirst
            lookup[seqIndex][msaSsFirst]=p;
          }
        }
        free(pl1);
        free(pf_struc);
      }
      free_pf_arrays();
    }
    if (cstruc!=NULL) free(cstruc);
    (void) fflush(stdout);
    free(structure);
    free(degappedSeq);
  }
}
PartitionFunctionLookup::~PartitionFunctionLookup ()
{
}
double PartitionFunctionLookup::GetBasePairProb (int seq,int msaSsFirst)
{
	double p=lookup[seq][msaSsFirst];
	if (p==0) {
		return 1e-6; // RNAfold doesn't like these low probabilities, so I shouldn't return zero
	}
	else {
		return p;
	}
}
void GenerateBasePairList(BasePairList& basePairList,std::string ss)
{
	basePairList.clear();
	for (int ssFirst=0; ssFirst<(int)(ss.size()); ssFirst++) {
		if (ss[ssFirst]=='<') {
			int ssLast=FindRightPartner(ss,ssFirst);
			BasePair bp;
			bp.ssFirst=ssFirst;
			bp.ssLast=ssLast;
			basePairList.push_back(bp);
		}
	}
}

void DetermineSeqBoundsForFragmentary(ValidSeqRegionVector& validSeqRegionVector,MSA *msa,bool fragmentary)
{
	// code copied from GSCConsensus.cpp, but I don't want to risk introducing bugs in the original context of the code

	int globalOuterFirst=0;
	int globalOuterLast=msa->alen;
	validSeqRegionVector.resize(msa->nseq);
	if (fragmentary) {

		const char FRAGMENTARY_SUB_REGION[]="FRAGMENTARY_SUB_REGION";
		for (int t=0; t<msa->ngc; t++) {
			//printf("#=GC %s\n",msa->gc_tag[t]);
			if (strcmp(FRAGMENTARY_SUB_REGION,msa->gc_tag[t])==0) {
				std::string ss=NormalizeSs(msa->gc[t],msa->alen);
				//printf("#=GC %s : %s\n",msa->gc_tag[t],ss.c_str());
				for (int i=0; i<msa->alen; i++) {
					if (ss[i]=='<') {
						globalOuterFirst=i;
					}
					if (ss[i]=='>') {
						globalOuterLast=i+1;
					}
				}
				break;
			}
		}
		int grt=-1;
		for (int t=0; t<msa->ngr; t++) {
			if (strcmp(FRAGMENTARY_SUB_REGION,msa->gr_tag[t])==0) {
				grt=t;
				break;
			}
		}

		for (int i=0; i<msa->nseq; i++) {

			int outerFirst=globalOuterFirst;
			int outerLast=globalOuterLast;
			if (grt!=-1) {
				printf("applying #=GR FRAGMENTARY_SUB_REGION for seq #%d\n",i);
				std::string ss=NormalizeSs(msa->gr[grt][i],msa->alen);
				for (int a=0; a<msa->alen; a++) {
					if (ss[a]=='<') {
						outerFirst=a;
					}
					if (ss[a]=='>') {
						outerLast=a+1;
					}
				}
			}

			validSeqRegionVector[i].first=outerFirst;
			while (validSeqRegionVector[i].first<msa->alen) {
				if (isalpha(msa->aseq[i][validSeqRegionVector[i].first])) {
					break;
				}
				validSeqRegionVector[i].first++;
			}
			validSeqRegionVector[i].last=outerLast;
			while (validSeqRegionVector[i].last-1>=0) {
				if (isalpha(msa->aseq[i][validSeqRegionVector[i].last-1])) {
					break;
				}
				validSeqRegionVector[i].last--;
			}
		}
	}
	else {
		for (int i=0; i<msa->nseq; i++) {
			validSeqRegionVector[i].first=0;
			validSeqRegionVector[i].last=msa->alen;
		}
	}
}

void HmmForwardScore_SS (HmmPairScores& hmmPairScores,std::string& extraInfo,MSA *msa,char *ssRaw,double nonCanonPairThreshold,int numFlankingPolyN,bool addFlankingPolyNWhenTrainingCM,bool uniformDistributionOfProfileHmmStartsAndEnds,const char *positionsToIgnoreFileName,bool fragmentaryForCovaryingPairs,bool fragmentaryForHmmTraining)
{
	assertr(ssRaw!=NULL); // actually legitimate; user could have forgotten the SS line, but I don't think that'll happen

	if (!fragmentaryForCovaryingPairs && fragmentaryForHmmTraining) {
		throw SimpleStringException("not allowed, since we need to calculate validSeqRegionVector if we want to use it for HMM training");
	}
	ValidSeqRegionVector validSeqRegionVector;
	DetermineSeqBoundsForFragmentary(validSeqRegionVector,msa,fragmentaryForCovaryingPairs);

	std::string ss=NormalizeSs(ssRaw,msa->alen);

	BasePairList basePairList;
	GenerateBasePairList(basePairList,ss);

	char **dsq = DigitizeAlignment(msa->aseq, msa->nseq, msa->alen);
	DegappedSeqFromMsaVector degappedSeqFromMsaVector;
	DegapSeqs_FindMap_AddPolyN(degappedSeqFromMsaVector,msa,dsq,numFlankingPolyN);

	// make hmm
	CovarianceModel cm;
	InfernalHmm infernalHmm;
	vector<int> msaPosToHmmState;
	HmmForwardScore_SS_MakeHmm(cm,infernalHmm,msaPosToHmmState,degappedSeqFromMsaVector,numFlankingPolyN,msa,addFlankingPolyNWhenTrainingCM,uniformDistributionOfProfileHmmStartsAndEnds,validSeqRegionVector,fragmentaryForHmmTraining,ss);
	//infernalHmm.SaveInBinary("data/t/test-score.hmm");
	// convert to HmmType1
	HmmType1 hmm;
	hmm.Init(infernalHmm);

	Stopwatch_t     *watchp2i;
	watchp2i = StopwatchCreate();
	StopwatchZero(watchp2i);
	StopwatchStart(watchp2i);
	PositionsToIgnore positionsToIgnore(positionsToIgnoreFileName,msa);
	StopwatchStop(watchp2i);
	StopwatchDisplay(stdout,(char *)"\nPositionsToIgnore CPU time: ", watchp2i);
	StopwatchFree(watchp2i);

	PartitionFunctionLookup partitionFunctionLookup(msa,degappedSeqFromMsaVector,basePairList,numFlankingPolyN);

	// we'll calculate   Pr(emits pair|emits seq,hmm) = Pr(emits pair^emits seq | hmm) / Pr(emits seq | hmm)
	
	// we now initialize the denominator for all seqs.  (We do this for convenience, even though ~50% of motifs don't even have covarying base pairs, according to our criteria)
	vector<HmmSsReal> probThatHmmEmitsSeqVector;
	probThatHmmEmitsSeqVector.resize(msa->nseq);
	for (int i=0; i<msa->nseq; i++) {
		probThatHmmEmitsSeqVector[i]=CalcProbThatSeqIsEmittedByHmm(hmm,degappedSeqFromMsaVector[i].dsq);
	}

#if 0
	// test code
	char tseq[]="AGUCAUCGAGUCAGUCAGUCACGUAGUCACGUACGUAUCGAGUCAGUCAGUCAGUCACGU";
	char *tdsqp=DigitizeSequence(tseq,(int)(sizeof(tseq)));
	vector<char> tdsq;
	tdsq.insert(tdsq.end(),tdsqp+1,tdsqp+strlen(tseq));
	NoUnderflowDouble pp=CalcProbThatSeqIsEmittedByHmm(hmm,tdsq);
#endif

	hmmPairScores.pairSum=0;
	hmmPairScores.pairSumWithPartitionFunc=0;
	for (BasePairList::const_iterator bpi=basePairList.begin(); bpi!=basePairList.end(); bpi++) {
		int ssFirst=bpi->ssFirst;
		int ssLast=bpi->ssLast;
		// pair is [ssFirst,ssLast)

		// check if this pair is Watson-Crick/G-U frequently enough
		double doubleGapWeight,gapWeight,nonCanonWeight,canonWeight;
		IntPairSet pairSet;
		GSCWeightedConsensus_CountPairFreqs_NoAdjustFreq(pairSet,doubleGapWeight,gapWeight,nonCanonWeight,canonWeight,msa,ssFirst,ssLast,validSeqRegionVector);
		//printf("GSCWeightedConsensus_CountPairFreqs_NoAdjustFreq: %lg,%lg,%lg\n",gapWeight,nonCanonWeight,canonWeight);
		double norm=gapWeight + canonWeight + nonCanonWeight;
		if (norm>0) {
			gapWeight /= norm;
			canonWeight /= norm;
			nonCanonWeight /= norm;

			// okay, just do the same test as the regular case -- we don't include the doubleGapWeight
			// in 'norm', in other words, we ignore it
			// probably, the HMM won't align these well anyway
			if (gapWeight + nonCanonWeight <= nonCanonPairThreshold) {
				// good pair

				ProbCacheForBasePair probCacheForBasePair;
				Init_ProbCacheForBasePair(probCacheForBasePair,msa->nseq);

				BestScore bestOfBasePair,baseOfBasePairPartitionFunc;

				// iterate over pairs of seqs, finding a seq-pair having base-pairs at this position that are both canonical
				// (i.e. Watson-Crick or G-U) and where the left and right nucleotides in the pair are different in each seq
				for (int i1=0; i1<msa->nseq; i1++) {
					for (int i2=0; i2<msa->nseq; i2++) {

						if (i1<i2) { // it's wrong to test a seq against itself, and it's redundant to test pairs in both orientations

							int leftNuc1,rightNuc1,leftNuc2,rightNuc2;
							bool hasDegen1,hasDegen2;
							GetNucPairForMsaPosition(leftNuc1,rightNuc1,hasDegen1,msa,i1,ssFirst,ssLast);
							GetNucPairForMsaPosition(leftNuc2,rightNuc2,hasDegen2,msa,i2,ssFirst,ssLast);

							// NOTE: we don't bother with fragmentary, since this only restricts pairs.
							// all that would happen is that we reject a pair because it has a gap in one nucleotide, instead of
							// rejecting the pair because it's beyond the fragmentary range in that seq -- same thing either way

							bool ignore=false;
							// we require that BOTH paired nucleotides are ignored before we ignore everything
							// the reason is we don't want to ignore stems that partially overlap a terminator, which could be a P1 stem of an aptamer
							if (positionsToIgnore.IgnorePosition(i1,ssFirst) && positionsToIgnore.IgnorePosition(i1,ssLast-1)) {
								ignore=true;
							}
							if (positionsToIgnore.IgnorePosition(i2,ssFirst) && positionsToIgnore.IgnorePosition(i2,ssLast-1)) {
								ignore=true;
							}

							if (// ignore positions that we're supposed to ignore
								ignore ||
								// ignore degenerate nucleotides; too much of a hassle
								hasDegen1 || hasDegen2 ||  
								// ignore gaps (we want pairs with both nucs present)
								leftNuc1==-1 || rightNuc1==-1 || leftNuc2==-1 || rightNuc2==-1
								) {
								// ignore this
							}
							else {
								if (isCanonPair[leftNuc1][rightNuc1] &&
								isCanonPair[leftNuc2][rightNuc2] &&
								leftNuc1!=leftNuc2 &&
								rightNuc1!=rightNuc2) {

									// okay, the base-pairs in both sequences are canonical, and the nucleotides
									// differ at both positions

									HmmSsReal prob1=CalcProbThatPairIsAlignedByHmm(probCacheForBasePair,hmm,degappedSeqFromMsaVector,i1,ssFirst,ssLast);
									HmmSsReal prob2=CalcProbThatPairIsAlignedByHmm(probCacheForBasePair,hmm,degappedSeqFromMsaVector,i2,ssFirst,ssLast);
									HmmSsReal thisPairEmitProb=(prob1/probThatHmmEmitsSeqVector[i1]) * (prob2/probThatHmmEmitsSeqVector[i2]);
									double partFunc1=partitionFunctionLookup.GetBasePairProb(i1,ssFirst);
									double partFunc2=partitionFunctionLookup.GetBasePairProb(i2,ssFirst);
									HmmSsReal thisPairEmitProbPartitionFunc = thisPairEmitProb * partFunc1 * partFunc2;
	//								extraInfo += stringprintf("#=GF pairEmitProb [%d,%d) %d %d %lg (%lg).  pf=%lg (%lg) [%lg,%lg]\n",ssFirst,ssLast,i1,i2,(double)(thisPairEmitProb),log2(thisPairEmitProb),(double)(thisPairEmitProbPartitionFunc),log2(thisPairEmitProbPartitionFunc),partFunc1,partFunc2);

									bestOfBasePair.ProcessNewScore(thisPairEmitProb,ssFirst,i1,i2);
									baseOfBasePairPartitionFunc.ProcessNewScore(thisPairEmitProbPartitionFunc,ssFirst,i1,i2);
								}
							}
						}
					}
				}

				hmmPairScores.bestPair.ProcessNewScore(bestOfBasePair);
				hmmPairScores.bestPairPartitionFunc.ProcessNewScore(baseOfBasePairPartitionFunc);
				hmmPairScores.pairSum += bestOfBasePair.GetScore();
				hmmPairScores.pairSumWithPartitionFunc += baseOfBasePairPartitionFunc.GetScore();
			}
		}
	}

	Free2DArray((void**)dsq, msa->nseq);
}
void HmmForwardScore_SS (char *stoFileName,double nonCanonPairThreshold,int numFlankingPolyN,bool uniformDistributionOfProfileHmmStartsAndEnds,const char *positionsToIgnoreFileName,bool fragmentaryForCovaryingPairs,bool fragmentaryForHmmTraining)
{
	MSA             *msa;         /* a multiple sequence alignment           */
	MSAFILE         *afp;         /* open alignment file                     */
	if ((afp = MSAFileOpen(stoFileName,MSAFILE_UNKNOWN, NULL)) == NULL) {
		Die((char *)"Alignment file %s could not be opened for reading", stoFileName);
	}
	bool addFlankingPolyNWhenTrainingCM=true;
	Stopwatch_t     *watch;
	watch = StopwatchCreate();
	StopwatchZero(watch);
	StopwatchStart(watch);
	while ((msa = MSAFileRead(afp)) != NULL) {
		GSCWeights(msa->aseq, msa->nseq, msa->alen, msa->wgt);
		HmmPairScores hmmPairScores;
		std::string extraInfo;

		// only bother with the main SS_cons; for motif predictions, we won't have pseudoknots anyway
		HmmForwardScore_SS(hmmPairScores,extraInfo,msa,msa->ss_cons,nonCanonPairThreshold,numFlankingPolyN,addFlankingPolyNWhenTrainingCM,uniformDistributionOfProfileHmmStartsAndEnds,positionsToIgnoreFileName,fragmentaryForCovaryingPairs,fragmentaryForHmmTraining);

		printf("\n%s\n",extraInfo.c_str());

		hmmPairScores.bestPair.PrintForUser(stdout,"SCORE bestPair");
		printf("SCORE pairSum score=%lg,%lg\n",hmmPairScores.pairSum.ToDouble_ZeroOnUnderflow(),hmmPairScores.pairSum.Log2());
		hmmPairScores.bestPairPartitionFunc.PrintForUser(stdout,"SCORE bestPairPf");
		printf("SCORE pairSumPf score=%lg,%lg\n",hmmPairScores.pairSumWithPartitionFunc.ToDouble_ZeroOnUnderflow(),hmmPairScores.pairSumWithPartitionFunc.Log2());
		printf("RESPECTS-IGNOREPOSITIONS\n");

		break; // only do one
	}
	MSAFree(msa);
	MSAFileClose(afp);
	StopwatchStop(watch);
	StopwatchDisplay(stdout, (char *)"\nOverall CPU time: ", watch);
	StopwatchFree(watch);
}

// yes, I know I have similar code in CmalignZasha.cpp, but what I need here is very simple, and it's easier to just implement that rather than re-purpose the CmalignZasha.cpp code
// this is only appropriate for profile HMMs
class MsaByColumn {
	int numColumns,numSeqs;
	bool addInReverseOrder;
	typedef vector2d<vector<char> > Msa;
	Msa msa;
	typedef vector<std::string> StringVector;
	StringVector names;
	void Pad (FILE *f,char ch,int n) const;
public:
	MsaByColumn (int numColumns_,int numSeqs_,bool addInReverseOrder_); // addInReverseOrder==true --> calling the .Add method is like a push_front (but this is simulated by reversing things later, for efficiency
	~MsaByColumn ();
	void SetSeqName (int seq,const std::string& name);
	void Add (int col,int seq,char nucAscii);  // first nuc in (col,seq) is assumed to be the consensus one
	void DumpStockholm (const char *fileName) const;
	void DumpStockholm (FILE *f) const;
};
MsaByColumn::MsaByColumn (int numColumns_,int numSeqs_,bool addInReverseOrder_)
{
	numColumns=numColumns_;
	numSeqs=numSeqs_;
	addInReverseOrder=addInReverseOrder_;
	msa.resize(numSeqs,numColumns);
	names.resize(numSeqs);
}
MsaByColumn::~MsaByColumn ()
{
}
void MsaByColumn::SetSeqName (int seq,const std::string& name)
{
	names[seq]=name;
}
void MsaByColumn::Add (int col,int seq,char nuc)
{
	msa[seq][col].push_back(nuc);
}
void MsaByColumn::DumpStockholm (const char *fileName) const
{
	FILE *f=ThrowingFopen(fileName,"wt");
	DumpStockholm(f);
	fclose(f);
}
void MsaByColumn::DumpStockholm (FILE *f) const
{
	fprintf(f,"# STOCKHOLM 1.0\n");

	// find column widths
	vector<int> colWidths;
	colWidths.assign(numColumns,0);
	for (int col=0; col<numColumns; col++) {
		for (int seq=0; seq<numSeqs; seq++) {
			colWidths[col]=std::max(colWidths[col],(int)(msa[seq][col].size()));
		}
	}
	// find name width
	int nameWidth=0;
	for (int seq=0; seq<numSeqs; seq++) {
		nameWidth=std::max(nameWidth,(int)(names[seq].size()));
	}

	// dump seqs
	for (int seq=0; seq<numSeqs; seq++) {
		fprintf(f,"%s ",names[seq].c_str());
		Pad(f,' ',nameWidth-(int)(names[seq].size()));
		for (int col=0; col<numColumns; col++) {
			if (addInReverseOrder) {
				for (int i=(int)(msa[seq][col].size())-1; i>=0; i--) {
					fprintf(f,"%c",msa[seq][col][i]);
				}
			}
			else {
				for (size_t i=0; i<msa[seq][col].size(); i++) {
					fprintf(f,"%c",msa[seq][col][i]);
				}
			}
			Pad(f,'.',colWidths[col]-(int)(msa[seq][col].size()));
		}
		fprintf(f,"\n");
	}

	fprintf(f,"//\n");
}
void MsaByColumn::Pad (FILE *f,char ch,int n) const
{
	for (int i=0; i<n; i++) {
		fprintf(f,"%c",ch);
	}
}

