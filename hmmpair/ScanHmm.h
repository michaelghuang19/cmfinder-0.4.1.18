
extern FILE *scoreDumpFile;

/*
class ScoreCollectorAbstract {
public:
	// topLevelState is the entry state on the right side of the HMM.  Normally it's the end state of the HMM, but if doLocal==true, it might not be
	virtual void GotScore (int windowLast,int topLevelRightwardState,float score) = 0;
	virtual void GotInfernalScore (int windowLast,float score) = 0;
	// there are 3 decisions that could be made: a transition, a leftwardBeginsc, an endsc (the latter 2 are only for doLocal==true)
	virtual void DecidedTransition(int windowLast,int state,int bestChildNum) = 0;
	virtual void UsedLeftwardBeginsc (int windowLast,int state) = 0;
	virtual void UsedEndsc (int windowLast,int state,int linkNum,int leftWindowLast) = 0;
	// plus special case start state, for !doLocal
	virtual void UsedStartState (int windowLast,int state) = 0;
	// nuc could be degenerate
	virtual void DidEmission (int windowLast,int nuc) = 0;
};
*/

// score collector that does nothing
class ScoreCollector_DeepNULL {
public:
	inline void StartInterval (int firstWindowLast,int lastWindowLast) {
	}
	inline void NotifyNumStates (int windowLast,int numStates) {
	}
	inline void GotScore (int windowLast,int topLevelRightwardState,float score) {
	}
	inline void GotInfernalScore (int windowLast,float score) {
	}
	inline void DecidedTransition(int windowLast,int state,int bestChildNum) {
	}
	inline void DidEmission (int windowLast,int nuc) {
	}
	inline void UsedLeftwardBeginsc (int windowLast,int state) {
	}
	inline void UsedEndsc (int windowLast,int state,int linkNum,int leftWindowLast) {
	}
	inline void UsedStartState (int windowLast,int state) {
	}
	inline void NotifyOfScoreInTable(int windowLast,int state,float score) {
	}
	inline void NotifyOfScoreInTable(int windowLast,int state,double score) {
	}
};
struct TwoTotals; // forward decl
class ScoreCollector_NULL : public ScoreCollector_DeepNULL {
public:
	/*
	inline void NotifyOfScoreInTable(int windowLast,int state,const SymbolicMath::Expression& score) {
	}
	 */
	inline void NotifyOfScoreInTable(int windowLast,int state,const TwoTotals& score) {
	}
	inline void NotifyOfScoreInTable(int windowLast,int state,float score) {
	}
	inline void NotifyOfScoreInTable(int windowLast,int state,double score) {
	}
	inline void NotifyOfScoreInTable(int windowLast,int state,NoUnderflowDouble score) {
	}
};

// ScoreCollector that keeps the table
template <class Real>
class DynProgTableScoreCollector : public ScoreCollector_NULL {
	vector<vector<Real> > *dynProgTable;
public:
	DynProgTableScoreCollector (int rnaSequenceLen,vector<vector<Real> >& dynProgTable_) { // dynProgTable[windowLast][state]
		dynProgTable=&dynProgTable_;
		dynProgTable->resize(rnaSequenceLen+1);
	}
	inline void NotifyNumStates (int windowLast,int numStates) {
		(*dynProgTable)[windowLast].resize(numStates);
	}
	inline void NotifyOfScoreInTable(int windowLast,int state,Real score) {
		(*dynProgTable)[windowLast][state]=score;
	}

	inline Real GetScore (int windowLast,int state) const {
		return (*dynProgTable)[windowLast][state];
	}
};
template <class Real>
class DynProgTableScoreCollector_OwnMemory : public ScoreCollector_NULL {
	vector<vector<Real> > concreteDynProgTable;
public:
	DynProgTableScoreCollector_OwnMemory () {
	}
	void Init (int rnaSequenceLen) {
		concreteDynProgTable.resize(rnaSequenceLen+1);
	}
	DynProgTableScoreCollector_OwnMemory(int rnaSequenceLen) {
		Init(rnaSequenceLen);
	}
	~DynProgTableScoreCollector_OwnMemory () {
	}
	inline void NotifyNumStates (int windowLast,int numStates) {
		concreteDynProgTable[windowLast].resize(numStates);
	}
	inline void NotifyOfScoreInTable(int windowLast,int state,Real score) {
		concreteDynProgTable[windowLast][state]=score;
	}

	inline Real GetScore (int windowLast,int state) const {
		return concreteDynProgTable[windowLast][state];
	}
};


class TracebackInfo : public ScoreCollector_DeepNULL {
public:
	enum Decision {
		Transition,Beginsc,Endsc,StartState
	};
protected:
	int windowLen;
	int slidingWindowSize;
	vector2d<Decision> decisionType; // indexed by [windowLast%slidingWindowSize][state]
	vector2d<int> transitionTaken; // indexed by [windowLast%slidingWindowSize][state].  If decisionType=Transition, this is the childNum, if decisionType=Endsc, it's the linkNum, if decisionType=Beginsc, it's unused 
	vector2d<int> leftWindowLastForEndsc; // indexed same way as transitionTaken
	vector<int> emissionDone; // indexed by [windowLast%slidingWindowSize]
	vector<int> visitsToState; // indexed by [state]
	vector<float> scores; // indexed by [windowLast%slidingWindowSize]
	const HmmType1& hmm;
public:
	TracebackInfo (const TracebackInfo& t);
	TracebackInfo (const HmmType1& hmm_,int windowLen,int slidingWindowSize_=-1); // set slidingWindowSize=length(seq) to store all info.  If slidingWindowSize==-1, it'll be set to 4*windowLen, a more aggressive size for TransitionCounter
	~TracebackInfo ();

	Decision GetDecision(int windowLast,HmmType1::State state);
	int GetTransition (int windowLast,HmmType1::State state);
	int GetEmission (int windowLast);

	void CopyFrom (const TracebackInfo& t); // but HMM must have been set to the same
	const HmmType1& GetHmm (void) const;
	float GetScore (int windowLast);
	int GetSlidingWindowSize (void) const;
};

class AbstractTraceback {
public:
	virtual ~AbstractTraceback ();

	virtual void StartInterval (int firstWindowLast,int lastWindowLast);
	virtual void AddSample (int windowLast,float score);
	virtual void DoneSample (void);
	virtual void AddRightBegin (HmmType1::State state);
	virtual void AddEmit (HmmType1::State state,int windowLast,int nuc);
	virtual void AddLeftBegin (HmmType1::State state);
	virtual void AddTransition (HmmType1::State state,int child);
	virtual void AddEnd (HmmType1::State state,int link,int nucsSkipped);
	virtual void AddNumTransitionsInParse(int numTransitionsInParse);
	virtual void VisitingState(int state);
	virtual void DoneVisitingState(void);
	virtual bool IgnoreWindowLastBelowWindowLen (void) const = 0;
};

extern void ScanHmm (HitList& hmmHitList,const InfernalHmm& hmm,float minLodScoreForHit,CykscanStats& scancykStats,const char *rnaSequence, int rnaSequenceLen,int windowLen);
extern void ScanHmm_HmmType1Float (HitList& hmmHitList,const HitList& inputHitList,const InfernalHmm& hmmInput,float float_minLodScoreForHit,CykscanStats& cykscanStats,const char *rnaSequence, int windowLen);
extern void ScanHmm_HmmType1Float (HitList& hmmHitList,const HitList& inputHitList,const HmmType1& hmmInput,float float_minLodScoreForHit,CykscanStats& cykscanStats,const char *rnaSequence, int windowLen);
extern void ScanHmm_HmmType1Float_VectorizeScores (HitList& hmmHitList,const HitList& inputHitList,const InfernalHmm& hmmInput,float minLodScoreForHit,CykscanStats& cykscanStats,const char *rnaSequence, int windowLen,vector<float>& scores);
extern void ScanHmm_HmmType1Float_VectorizeScores_InfernalToo (HitList& hmmHitList,const HitList& inputHitList,const InfernalHmm& hmmInput,float minLodScoreForHit,CykscanStats& cykscanStats,const char *rnaSequence, int windowLen,vector<float>& hmmScores,vector<float>& cmScores);
/*
extern void ScanHmm_HmmType1Float_AbstractScores (HitList& hmmHitList,const HitList& inputHitList,const InfernalHmm& hmmInput,float minLodScoreForHit,CykscanStats& cykscanStats,const char *rnaSequence, int windowLen,ScoreCollectorAbstract *scoreCollector);
*/
extern void ScanHmm_HmmType1Float_Traceback (HitList& hmmHitList,const HitList& inputHitList,const HmmType1& hmmInput,float minLodScoreForHit,CykscanStats& cykscanStats,const char *rnaSequence, int windowLen,TransitionCounter& transitionCounter,double& totalScore);
extern void ScanHmm_HmmType1Float_TracebackInfo (HitList& hmmHitList,const HitList& inputHitList,const HmmType1& hmmInput,float minLodScoreForHit,CykscanStats& cykscanStats,const char *rnaSequence, int windowLen,TransitionCounter& transitionCounter,TracebackInfo& tracebackInfo);
void ScanHmm_HmmType1Float_AbstractTraceback (HitList& hmmHitList,const HitList& inputHitList,const HmmType1& hmmInput,float minLodScoreForHit,CykscanStats& cykscanStats,const char *rnaSequence, int windowLen,AbstractTraceback& traceback,TracebackInfo& tracebackInfo);


// class to do scanning that avoids expensive re-initialization
/*
class ScanHmmWithSubCm_HmmType1Float {
protected:

	const InfernalHmm *infernalInputHmm; // should be constant within the life of the object
	HmmType1 hmm;
	HmmOptimizedForRememberNucScanning *rememberNucHmm;
	HMMWithSubCMData hmmWithSubCMData;

	void *pvoid_scanHmmAugmenter; // stored as void pointer so I don't have to declare the type in this header file & put a huge amount of templated code here
	template <class ScoreCollector>
	void InternalScanOrVectorize (HitList& hmmHitList,const HitList& inputHitList,const CovarianceModel& cm,const InfernalHmm& hmmInput,float minLodScoreForHit,CykscanStats& cykscanStats,const char *rnaSequence,int windowLen,const HMMWithSubCMData& hmmWithSubCMData,ScoreCollector& scoreCollector); // scoreCollector!=NULL --> set them
public:
	ScanHmmWithSubCm_HmmType1Float (const CovarianceModel& cm,const InfernalHmm& hmmInput,int windowLen,const HMMWithSubCMData& hmmWithSubCMData);
	~ScanHmmWithSubCm_HmmType1Float ();
	void Scan (HitList& hmmHitList,const HitList& inputHitList,const CovarianceModel& cm,const InfernalHmm& hmmInput,float minLodScoreForHit,CykscanStats& cykscanStats,const char *rnaSequence,int windowLen,const HMMWithSubCMData& hmmWithSubCMData);
	void VectorizeScores  (HitList& hmmHitList,const HitList& inputHitList,const CovarianceModel& cm,const InfernalHmm& hmmInput,float minLodScoreForHit,CykscanStats& cykscanStats,const char *rnaSequence,int windowLen,const HMMWithSubCMData& hmmWithSubCMData,vector<float>& scores);
};
 */

extern void ForwardAlgToGetExpectedScore_NoCaching (const InfernalHmm& infernalHmm,int windowLen,double& expectedScore);
extern void ForwardAlgToGetExpectedScore_NoCaching (const InfernalHmm& infernalHmm,int windowLen,double& expectedScore,const InfernalHmm::State startState,const InfernalHmm::State endState);

/*
extern double InfiniteLengthForwardAlg_LastStateDoesntEmit_HmmType1Double (const InfernalHmm& infernalHmm,const InfernalHmm::State infernalStartState,const InfernalHmm::State infernalEndState,MarkovModelStats& markovModelStats);
extern double InfiniteLengthForwardAlg_LastStateDoesntEmit_HmmType1Double (const InfernalHmm& infernalHmm,MarkovModelStats& markovModelStats);
extern double InfiniteLengthForwardAlg_LastStateDoesntEmit_HmmType1Double (const HmmType1& hmm,MarkovModelStats& markovModelStats);
extern double InfiniteLengthForwardAlg_LastStateDoesntEmit_HmmWithPenaltiesDouble (const HmmWithPenalties& hmm,MarkovModelStats& markovModelStats);
extern double CalcExpectedEmitProb_0order_InfernalHmmDouble (const InfernalHmm& infernalHmm,InfernalHmm::State state,MarkovModelStats& markovModelStats);
 */

double ForwardOfLiteralString (const HmmType1& hmm,const char *rnaSequence,int rnaSequenceLen,bool allowOnlyFullString,bool epsilonOnlyEnds);
double ForwardOfLiteralString (vector<vector<double> >& dynProgTable /*[windowLast][state]*/,const HmmType1& hmm,const char *rnaSequence,int rnaSequenceLen,bool allowOnlyFullString,bool epsilonOnlyEnds);


// helper function for InfiniteLengthForwardAlg_LastStateDoesntEmit (see below)
/*
template <class Hmm,class Real>
inline Real CalcExpectedEmitProb_0order (const Hmm& hmm,const typename Hmm::State state,MarkovModelStats& markovModelStats)
{
	if (hmm.IsEmittingState(state)) {
		Real emitProb=0.0;
		for (int nuc=0; nuc<Alphabet_size; nuc++) {
			emitProb += markovModelStats.GetProbOfNuc_0order(nuc) * hmm.GetSingletEmissionProb(state,nuc);
		}
		return emitProb;
	}
	else {
		return 1;
	}
}
*/