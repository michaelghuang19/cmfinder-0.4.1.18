#include "hmmpair.h"

#include "QueueWithFindMin.h"

//#define SCAN_DUMP // deprecated in favor of '--score-dump-file' command flag

#ifdef SCAN_DUMP
FILE *scanDumpFile=NULL;
#endif

FILE *scoreDumpFile=NULL;

typedef DynProgTableScoreCollector<double> DoubleDynProgTableScoreCollector;

// gets them in an array
class ScoreCollector_Vectorize : public ScoreCollector_NULL {
public:
	vector<float> *scores,*infernal_scores;
	inline void GotScore (int windowLast,int topLevelRightwardState,float score) {
		(*scores)[windowLast]=score;
	}
	inline void GotInfernalScore (int windowLast,float score) {
		if (infernal_scores!=NULL) {
			(*infernal_scores)[windowLast]=score;
		}
	}
};



//////////////////////////////////
//
//   ViterbiOrForwardClosedSemiring : template parameter to ClassScanHmm that codes
//   for the Viterbi alg (using addition and max on log score) or the 
//   Forward alg (using multiply and add on actual probabilities).
//   A "closed semiring" is a mathematical object that defines the properties
//   that the Viterbi & Forward algs have in common that allow us to use
//   the same dynamic program alg to solve them.  It's the same as for
//   shortest path problems on graphs (or I should have read the textbook
//   more closely...)
//
//   A closed semiring has a "summary" operator, which says how to combine the values of many in-edges.  For Viterbi, summary op is max, for Forward it's plus.  Closed semiring has "extension" operator, which says how to extend the value of two adjacent edges in the graph.  For Viterbi, extension op is plus; for Forward it's multiply.
//   In the case of this alg, we also need operations to access the transition/emission scores, since we must use either the log probs for Viterbi, or the actual probs for Forward

template <class Hmm>
class Viterbi_ClosedSemiring_GenericHmm {
public:
	// operations

	inline static void AssertIsViterbi (void) {} // you'll get a compiler error if you call this on Forward_ClosedSemiRing

	// Extension operator, e.g. Addition for Viterbi
	inline static float ExtensionOperator (float x1,float x2) {
		return x1+x2;
	}
	// Summary operator, e.g. max for Viterbi
	inline static float SummaryOperator (float x1,float x2) {
		return std::max(x1,x2);
	}
	// Summary operator, with allowance for keeping track of what happened (traceback for the Viterbi alg)
	inline static bool SummaryOperatorWithMaxNotification (float& currValue,float valueToCombine) {
		if (valueToCombine>=currValue) {
			currValue=valueToCombine;
			return true;
		}
		else {
			return false;
		}
	}
	// Summary operator, with allowance for keeping track of the max child
	inline static void SummaryOperatorWithTraceback (float& currValue,float valueToCombine,int& bestChildNum,int childNum) {
		if (SummaryOperatorWithMaxNotification(currValue,valueToCombine)) {
			bestChildNum=childNum;
		}
	}
	// tells code whether Summary operator requires traceback
	inline static bool SummaryOperatorLikesTraceback (void) {
		return true;
	}

	// constants
	// impossible is the identity for SummaryOperator and the annihilator for ExtensionOperator
	inline static float GetImpossible (void) {
		return (float)IMPOSSIBLE;
	}
	// this is for numerical reasons in the Viterbi alg.  For Forward, it's same as GetImpossible
	inline static float GetExtraImpossible (void) {
		return -FLT_MAX;
	}
	inline static bool IsImpossible (float score) {
		return score==GetImpossible();
	}
	// synonym
	inline static float GetSummaryOperatorIdentity (void) {
		return GetImpossible();
	}
	// identity for ExtensionOperator
	inline static float GetExtensionOperatorIdentity (void) {
		return 0;
	}

	// should we use log probs, or actual probs?
	inline static float GetNthChildTransition (const Hmm& hmm,const typename Hmm::State state,const int childNum) {
		return hmm.GetNthChildTsc(state,childNum);
	}
	inline static float GetSingletEmission (const Hmm& hmm,const typename Hmm::State state,const int nuc) {
		return hmm.GetSingletEmissionScore(state,nuc);
	}
	inline static float GetLeftwardBegin (const Hmm& hmm,const typename Hmm::State state) {
		return hmm.GetLeftwardBeginsc(state);
	}
	inline static float GetRightwardBegin (const Hmm& hmm,const typename Hmm::State state) {
		return hmm.GetRightwardBeginsc(state);
	}
	inline static float GetEnd (const Hmm& hmm,const typename Hmm::State state,int link) {
		return hmm.GetEndscLinksToLeft_Endsc(state,link);
	}
	inline static float GetLocalEndSelfLoop (const Hmm& hmm) {
		return hmm.GetLocalEndSelfLoopScore();
	}

	inline static bool IsViterbi (void) {
		return true;
	}
	
	template <class Real>
	static void ProcessScoreAtInit(Real initBestScore,float float_minLodScoreForHit,int startPosToScan,int endPosToScan,int windowLen,const CykscanStats& cykscanStats) {
	}

	// this function is in ClosedSemiring because the alternative, having it in ClassScanHmm, means we have all sorts of >= type function for the Viterbi case that are undefined for SymbolicMath::Expression.  So, it's easiest to separate them.
	template <class Real,class ScoreCollector>
	inline static void ProcessScoreAtPos(Real& totalSummaryScore,Real bestScoreAtPos,float float_minLodScoreForHit,int startPosToScan,int endPosToScan,int windowLast,int windowLen,HitList& hmmHitList,const CykscanStats& cykscanStats,ScoreCollector& scoreCollector)
	{
		// in this context, Viterbi means scanning to find hits -- the standard functionality of this program

		//int pos=windowLast-1;

		if (bestScoreAtPos >= float_minLodScoreForHit) {
			std::pair<int,int> thisWindow;
			thisWindow.first=std::max(startPosToScan,windowLast-windowLen);
			thisWindow.second=windowLast;

			AddWindowToList(hmmHitList,thisWindow);
		}
		if (cykscanStats.collectHmmScoresPerWindowLast) {
			double pureHmmScore=cykscanStats.hmmScoresPerWindowLast[windowLast];
			double float_bestScoreAtPos=bestScoreAtPos;
			double scaledDiff=(float_bestScoreAtPos - pureHmmScore)/std::max(0.1,fabs(pureHmmScore+float_bestScoreAtPos));
			if (scaledDiff >= 0.0005) {
				fflush(stdout);
				assert(false);
				throw SimpleStringException("bestScoreAtPos > pureHmmScore, which isn't supposed to happen.  %s:%d, windowLast=%d (in interval [%d,%d) ), myScore=%lf, hmmScore=%lf\n",__FILE__,__LINE__,windowLast,startPosToScan,endPosToScan,float_bestScoreAtPos,pureHmmScore);
			}
		}
	}
};
// typical Hmm type we'll be using, so I don't have to change so much code
typedef Viterbi_ClosedSemiring_GenericHmm<HmmType1> Viterbi_ClosedSemiring;


// for comments, see class Viterbi_ClosedSemiring (above)
class Forward_ClosedSemiring {
public:
	// with doubles
	inline static double ExtensionOperator (double x1,double x2) {
		return x1*x2;
	}
	inline static double SummaryOperator (double x1,double x2) {
		return x1+x2;
	}
	inline static bool SummaryOperatorWithMaxNotification (double& currValue,double valueToCombine) {
		// no traceback stuff is needed (or meaningful) for Forward alg
		currValue=SummaryOperator(currValue,valueToCombine);
		return false;
	}
	inline static void SummaryOperatorWithTraceback (double& currValue,double valueToCombine,int& bestChildNum,int childNum) {
		SummaryOperatorWithMaxNotification(currValue,valueToCombine);
	}
	inline static bool SummaryOperatorLikesTraceback (void) {
		return false;
	}
	inline static float GetImpossible (void) {
		return 0.0;
	}
	inline static float GetExtraImpossible (void) {
		return 0.0;
	}
	inline static bool IsImpossible (double score) {
		return score==GetImpossible();
	}
	inline static float GetSummaryOperatorIdentity (void) {
		return GetImpossible();
	}
	inline static float GetExtensionOperatorIdentity (void) {
		return 1.0;
	}
	inline static float GetNthChildTransition (const HmmType1& hmm,const HmmType1::State state,const int childNum) {
		return hmm.GetNthChildTransitionProb(state,childNum);
	}
	inline static float GetSingletEmission (const HmmType1& hmm,const HmmType1::State state,const int nuc) {
		return hmm.GetSingletEmissionProb(state,nuc);
	}
	inline static float GetLeftwardBegin (const HmmType1& hmm,const HmmType1::State state) {
		return hmm.GetLeftwardBeginProb(state);
	}
	inline static float GetRightwardBegin (const HmmType1& hmm,const HmmType1::State state) {
		return hmm.GetRightwardBeginProb(state);
	}
	inline static float GetEnd (const HmmType1& hmm,const HmmType1::State state,int link) {
		return hmm.GetEndProb(state,link);
	}
	inline static float GetLocalEndSelfLoop (const HmmType1& hmm) {
		NotImplemented(); // not implemented for forward case, since I don't think I need it
	}
	inline static bool IsViterbi (void) {
		return false;
	}

	static void ProcessScoreAtInit(double initBestScore,float float_minLodScoreForHit,int startPosToScan,int endPosToScan,int windowLen,const CykscanStats& cykscanStats) {
	}
	template <class ScoreCollector>
	inline static void ProcessScoreAtPos(double& totalSummaryScore,double bestScoreAtPos,float float_minLodScoreForHit,int startPosToScan,int endPosToScan,int windowLast,int windowLen,HitList& hmmHitList,const CykscanStats& ckyscanStats,ScoreCollector& scoreCollector) {
		totalSummaryScore=bestScoreAtPos; // see SymbolicForward_ClosedSemiring for comments
	}
};

class Forward_ClosedSemiring_NoUnderflowDouble {
public:
	inline static NoUnderflowDouble ExtensionOperator (NoUnderflowDouble x1,NoUnderflowDouble x2) {
		NoUnderflowDouble r(x1);
		r *= x2;
		return r;
	}
	inline static NoUnderflowDouble SummaryOperator (NoUnderflowDouble x1,NoUnderflowDouble x2) {
		NoUnderflowDouble r(x1);
		r += x2;
		return r;
	}
	inline static bool SummaryOperatorWithMaxNotification (NoUnderflowDouble& currValue,NoUnderflowDouble valueToCombine) {
		// no traceback stuff is needed (or meaningful) for Forward alg
		currValue=SummaryOperator(currValue,valueToCombine);
		return false;
	}
	inline static void SummaryOperatorWithTraceback (NoUnderflowDouble& currValue,NoUnderflowDouble valueToCombine,int& bestChildNum,int childNum) {
		SummaryOperatorWithMaxNotification(currValue,valueToCombine);
	}
	inline static bool SummaryOperatorLikesTraceback (void) {
		return false;
	}
	inline static float GetImpossible (void) {
		return 0.0;
	}
	inline static float GetExtraImpossible (void) {
		return 0.0;
	}
	inline static bool IsImpossible (NoUnderflowDouble score) {
		return score==NoUnderflowDouble(0.0);
	}
	inline static float GetSummaryOperatorIdentity (void) {
		return GetImpossible();
	}
	inline static float GetExtensionOperatorIdentity (void) {
		return 1.0;
	}
	inline static float GetNthChildTransition (const HmmType1& hmm,const HmmType1::State state,const int childNum) {
		return hmm.GetNthChildTransitionProb(state,childNum);
	}
	inline static float GetSingletEmission (const HmmType1& hmm,const HmmType1::State state,const int nuc) {
		return hmm.GetSingletEmissionProb(state,nuc);
	}
	inline static float GetLeftwardBegin (const HmmType1& hmm,const HmmType1::State state) {
		return hmm.GetLeftwardBeginProb(state);
	}
	inline static float GetRightwardBegin (const HmmType1& hmm,const HmmType1::State state) {
		return hmm.GetRightwardBeginProb(state);
	}
	inline static float GetEnd (const HmmType1& hmm,const HmmType1::State state,int link) {
		return hmm.GetEndProb(state,link);
	}
	inline static float GetLocalEndSelfLoop (const HmmType1& hmm) {
		NotImplemented(); // not implemented for forward case, since I don't think I need it
	}
	inline static bool IsViterbi (void) {
		return false;
	}

	static void ProcessScoreAtInit(double initBestScore,float float_minLodScoreForHit,int startPosToScan,int endPosToScan,int windowLen,const CykscanStats& cykscanStats) {
	}
	template <class ScoreCollector>
	inline static void ProcessScoreAtPos(NoUnderflowDouble& totalSummaryScore,NoUnderflowDouble bestScoreAtPos,float float_minLodScoreForHit,int startPosToScan,int endPosToScan,int windowLast,int windowLen,HitList& hmmHitList,const CykscanStats& ckyscanStats,ScoreCollector& scoreCollector) {
		totalSummaryScore=bestScoreAtPos; // see SymbolicForward_ClosedSemiring for comments
	}
};
class ForwardCulmulative_ClosedSemiring_NoUnderflowDouble : public Forward_ClosedSemiring_NoUnderflowDouble {
public:
	template <class ScoreCollector>
	inline static void ProcessScoreAtPos(NoUnderflowDouble& totalSummaryScore,NoUnderflowDouble bestScoreAtPos,float float_minLodScoreForHit,int startPosToScan,int endPosToScan,int windowLast,int windowLen,HitList& hmmHitList,const CykscanStats& ckyscanStats,ScoreCollector& scoreCollector) {

		// culmulative.  See SymbolicForward for comments
		totalSummaryScore=SummaryOperator(totalSummaryScore,bestScoreAtPos);
	}
};

//////////////////////////////
//  ScanHmmAugmenter classes, for use as template params of ClassScanHmm that augment the scanning algorithm, mainly for secondary structure

// these enum values can be combined with boolean ORs
enum ScanHmmAugmenter_OverrideDynProgScoreOperation {
	OverrideDynProgScore_DontChange, // don't change the dyn prog score at all -- completely normal
	OverrideDynProgScore_LeaveOutEmitScoreButMoveWithEmit, // completely normal, except set emitScore=0 no matter what ("MoveWithEmit") means that if there was an emit, we should still look at the previous column
	OverrideDynProgScore_SetScore, // set the score to get_score (output param of the OverrideDynProgScore method
	OverrideDynProgScore_AddToScore, // completely normal, but add get_score to whatever the normal score is
	OverrideDynProgScore_AddToScore_And_NoEmit, // add get_score to whatever the normal score is, set emitScore==0 && use currTable (i.e. do NOT move with emit)
	OverrideDynProgScore_AddToScore_And_LeaveOutEmitScoreButMoveWithEmit, // combination
};
inline bool LeaveOutEmitScoreButMoveWithEmit (ScanHmmAugmenter_OverrideDynProgScoreOperation op) {
	return op==OverrideDynProgScore_LeaveOutEmitScoreButMoveWithEmit || op==OverrideDynProgScore_AddToScore_And_LeaveOutEmitScoreButMoveWithEmit;
}
inline bool AddToScore (ScanHmmAugmenter_OverrideDynProgScoreOperation op) {
	return op==OverrideDynProgScore_AddToScore || op==OverrideDynProgScore_AddToScore_And_LeaveOutEmitScoreButMoveWithEmit || op==OverrideDynProgScore_AddToScore_And_NoEmit;
}
template <class ScoreHelperType,class Hmm,class ScoreCollector,class ClosedSemiring>
class ScanHmmAugmenter_NULL {
public:
	inline void StartSequence (const Hmm& hmm,const char *rnaSequence) {}
	inline void StartInterval (const Hmm& hmm,int firstWindowLast,int lastWindowLast) {}
	inline void EndInterval (void) {}
	inline void StartWindowLast (int windowLast) {}
	inline void EndWindowLast (int windowLast) {}

	// allows the same calculation of score, but not counting emission, in case we want to splice something in, like oh say a sub-CM.  for overrideDynProgScoreOperation, see definition of enum, above.  The interpretation of get_score depends on the result of overrideDynProgScoreOperation.  windowLast and state define the dyn-prog-table entry being computed.  scoreCollector is in case this class wants to support it (which I'm probably too lazy to do).  emitScore is the normal emission score, in case the function wants to override the actual score, but still use the emit score.
	inline void OverrideDynProgScore (ScanHmmAugmenter_OverrideDynProgScoreOperation& overrideDynProgScoreOperation,typename ScoreHelperType::ScoreType &get_score,int windowLast,typename Hmm::State state,ScoreCollector& scoreCollector,typename ScoreHelperType::ScoreType emitScore,const Hmm& hmm) {
		overrideDynProgScoreOperation=OverrideDynProgScore_DontChange;
	}

	inline void NotifyOfScore (int windowLast,typename Hmm::State state,typename ScoreHelperType::ScoreType get_score) {}
	inline void ReNotifyOfScore_EndingWindowLast (int windowLast,typename Hmm::State state,typename ScoreHelperType::ScoreType get_score) {}
	
	inline bool OverrideFinalWindowLastScore (typename ScoreHelperType::ScoreType& get_score,int windowLast,typename ScoreHelperType::ScoreType score) { return false; }

	typename ScoreHelperType::ScoreType GetBestScoreInWindow (typename Hmm::State state,int windowLast,int& get_leftWindowLast) {assertr(false); }
	static bool Is_ScanHmmAugmenter_GenericLocal (void) {return false;}
};

template <class ScoreHelperType,class Hmm,class OtherScanHmmAugmenter>
class ScoreHmmsAsCm_Augmenter {
protected:
	const CovarianceModel& cm;
	const InfernalHmm& infernalHmm;
	const Hmm& hmm;
	OtherScanHmmAugmenter& otherScanHmmAugmenter;
public:
	ScoreHmmsAsCm_Augmenter (const CovarianceModel& cm_,const InfernalHmm& infernalHmm_,const Hmm& hmm_,OtherScanHmmAugmenter& otherScanHmmAugmenter_)
		: cm(cm_),infernalHmm(infernalHmm_),hmm(hmm_),otherScanHmmAugmenter(otherScanHmmAugmenter_)
	{
	}
	inline void StartSequence (const Hmm& hmm,const char *rnaSequence) {
		otherScanHmmAugmenter.StartSequence(hmm,rnaSequence);
	}
	inline void StartInterval (const Hmm& hmm,int firstWindowLast,int lastWindowLast) {
		otherScanHmmAugmenter.StartInterval(hmm,firstWindowLast,lastWindowLast); 
	}
	inline void EndInterval (void) {
		otherScanHmmAugmenter.EndInterval();
	}
	inline void StartWindowLast (int windowLast) {
		otherScanHmmAugmenter.StartWindowLast(windowLast);
	}
	inline void EndWindowLast (int windowLast) {
		otherScanHmmAugmenter.EndWindowLast(windowLast);
	}
	/*
	inline void OverrideDynProgScore (ScanHmmAugmenter_OverrideDynProgScoreOperation& overrideDynProgScoreOperation,typename ScoreHelperType::ScoreType &get_score,int windowLast,typename Hmm::State state,ScoreCollector& scoreCollector,typename ScoreHelperType::ScoreType emitScore,const Hmm& hmm) {
		otherScanHmmAugmenter.OverrideDynProgScore (overrideDynProgScoreOperation,get_score,windowLast,state,scoreCollector,emitScore,hmm);
	}
	*/
	inline void NotifyOfScore (int windowLast,typename Hmm::State state,typename ScoreHelperType::ScoreType get_score) {
		otherScanHmmAugmenter.NotifyOfScore (windowLast,state,get_score);
	}
	inline void ReNotifyOfScore_EndingWindowLast (int windowLast,typename Hmm::State state,typename ScoreHelperType::ScoreType get_score) {
		otherScanHmmAugmenter.ReNotifyOfScore_EndingWindowLast (windowLast,state,get_score);
	}
	typename ScoreHelperType::ScoreType GetBestScoreInWindow (typename Hmm::State state,int windowLast,int& get_leftWindowLast) {
		return otherScanHmmAugmenter.GetBestScoreInWindow (state,windowLast,get_leftWindowLast);
	}
	static bool Is_ScanHmmAugmenter_GenericLocal (void) {
		return OtherScanHmmAugmenter::Is_ScanHmmAugmenter_GenericLocal();
	}
};

// will remember scores, and give the max-value score between minWindowLen and maxWindowLen
class ScoreRememberer_VeryDumb {
protected:
	struct Datum {
		float score;
		int windowLast; // I told you this was really dumb, so we just keep this info around
	};
	vector<Datum> data; // circular array -- keeping the whole array would be even dumber, but unfortunatley that seems like too much RAM

	int circularArrayIndex;
	int minWindowLen,maxWindowLen;
public:
	ScoreRememberer_VeryDumb (void); // must call Init, then
	~ScoreRememberer_VeryDumb ();
	void Init (int minWindowLen,int maxWindowLen);
	void ReInit (void);
	void AddScore (int windowLast,float score);
	float GetBestScoreWithinWindow(int forWindowLast,int& get_leftWindowLast) const;
	bool IsInit (void) const;
	float GetScoreAt (int windowLast) const;
};
float ScoreRememberer_VeryDumb::GetScoreAt (int windowLast) const
{
	int dataArraySize=(int)(data.size());
	int index=windowLast % dataArraySize;
	assert(data[index].windowLast==windowLast);
	return data[index].score;
}
ScoreRememberer_VeryDumb::ScoreRememberer_VeryDumb (void)
{
	maxWindowLen=-1;
}
ScoreRememberer_VeryDumb::~ScoreRememberer_VeryDumb ()
{
}
bool ScoreRememberer_VeryDumb::IsInit (void) const
{
	return maxWindowLen!=-1;
}
void ScoreRememberer_VeryDumb::Init (int minWindowLen_,int maxWindowLen_)
{
	minWindowLen=minWindowLen_;
	maxWindowLen=maxWindowLen_;

	ReInit();
}
void ScoreRememberer_VeryDumb::ReInit (void)
{
	Datum datum;
	datum.windowLast=-maxWindowLen-1; // impossible value, effectively marking all cells invalid
	data.assign(maxWindowLen+1,datum);
	circularArrayIndex=0;
}
void ScoreRememberer_VeryDumb::AddScore (int windowLast,float score)
{
	assert(windowLast-data[circularArrayIndex].windowLast > maxWindowLen); // otherwise we didn't reserve enough space in data to remember everything, because we can't throw this away yet

	data[circularArrayIndex].score=score;
	data[circularArrayIndex].windowLast=windowLast;

	circularArrayIndex++;
	if (circularArrayIndex==(int)(data.size())) {
		circularArrayIndex=0;
	}
}
float ScoreRememberer_VeryDumb::GetBestScoreWithinWindow(int forWindowLast,int& get_leftWindowLast) const
{
	float bestScore=(float)IMPOSSIBLE;
	for (size_t i=0; i<data.size(); i++) {
		int distance=forWindowLast - data[i].windowLast;
		if (distance>=minWindowLen && distance<=maxWindowLen) {
			if (data[i].score>=bestScore) {
				bestScore=data[i].score;
				get_leftWindowLast=data[i].windowLast;
			}
		}
	}
	return bestScore;
}

// ScoreRememberer_Log: like ScoreRememberer_VeryDumb, but uses something log-time, like priority heap; used for local ends
template <class Queue>
class ScoreRememberer_Log_Base {
protected:
	Queue queue;
	int maxWindowLen;
	int firstWindowLastInQueue;

	void Leave (int windowLast) {
		if (firstWindowLastInQueue<0) {
			firstWindowLastInQueue=windowLast;
		}
		else {
			if (windowLast-firstWindowLastInQueue>maxWindowLen) {
				queue.Leave();
			}
		}
	}
public:
	ScoreRememberer_Log_Base (void) : queue(0) {
		maxWindowLen=-1;
	}
	~ScoreRememberer_Log_Base () {
	}
	void Init (int minWindowLen,int maxWindowLen_) { // minWindowLen must be 0
		assert(minWindowLen==0); // a restriction of this implementation
		maxWindowLen=maxWindowLen_;

		queue.resize(maxWindowLen+1);
		ReInit();
	}
	void ReInit (void) {
		firstWindowLastInQueue=-1;
		queue.ClearAll();
	}
	bool IsInit (void) const {
		return maxWindowLen!=-1;
	}
};
template <class Queue>
class ScoreRememberer_Log : public ScoreRememberer_Log_Base<Queue> {
public:
	ScoreRememberer_Log ()
		: ScoreRememberer_Log_Base<Queue> (){
	}
	void AddScore (int windowLast,float score) {
		this->Leave(windowLast);
		this->queue.Enter(score);
	}
	float GetBestScoreWithinWindow(int forWindowLast,int& get_leftWindowLast) const {
		return this->queue.FindMin();
	}
};
	struct ScoreRememberer_Log_WithLeftWindowLast_Datum {
		int windowLast;
		float score;
		bool operator > (const ScoreRememberer_Log_WithLeftWindowLast_Datum& t) const {
			return score>t.score;
		}
	};
template <class Queue>
class ScoreRememberer_Log_WithLeftWindowLast : public ScoreRememberer_Log_Base<Queue> {
public:
	ScoreRememberer_Log_WithLeftWindowLast () 
	: ScoreRememberer_Log_Base<Queue> () {
	}
	typedef ScoreRememberer_Log_WithLeftWindowLast_Datum Datum;
	void AddScore (int windowLast,float score) {
		this->Leave(windowLast);
		Datum datum;
		datum.windowLast=windowLast;
		datum.score=score;
		this->queue.Enter(datum);
	}
	float GetBestScoreWithinWindow(int forWindowLast,int& get_leftWindowLast) const {
		const Datum& datum=this->queue.FindMin();
		get_leftWindowLast=datum.windowLast;
		return datum.score;
	}
};

typedef ScoreRememberer_Log<QueueWithFindMin_CircularArrayAndHeap<float,std::greater<float> > > ScoreRememberer_PriorityHeap;
typedef ScoreRememberer_Log<OpaqueQueueWithFindMin_CircularArrayStlMultiset<float,std::greater<float> > > ScoreRememberer_StlSet;
typedef ScoreRememberer_Log<OpaqueQueueWithFindMin_CleverBinaryCircularArray<float,std::greater<float> > > ScoreRememberer_CleverBinaryCircularArray;
typedef ScoreRememberer_Log_WithLeftWindowLast<OpaqueQueueWithFindMin_CleverBinaryCircularArray<ScoreRememberer_Log_WithLeftWindowLast_Datum,std::greater<ScoreRememberer_Log_WithLeftWindowLast_Datum> > > ScoreRememberer_CleverBinaryCircularArray_WithLeftWindowLast;


// compare results of _VeryDumb and _PriorityHeap
class ScoreRememberer_SanityChecker {
protected:
	ScoreRememberer_VeryDumb dumb;
	ScoreRememberer_PriorityHeap heap;
	ScoreRememberer_StlSet set;
	ScoreRememberer_CleverBinaryCircularArray cleverCircular;
	ScoreRememberer_CleverBinaryCircularArray_WithLeftWindowLast cleverCircularWithLWL;
	int prevWindowLast;
public:
	ScoreRememberer_SanityChecker ();
	~ScoreRememberer_SanityChecker ();
	void Init (int minWindowLen,int maxWindowLen); // minWindowLen must be 0
	void ReInit (void);
	void AddScore (int windowLast,float score);
	float GetBestScoreWithinWindow(int forWindowLast,int& get_leftWindowLast) const;
	bool IsInit (void) const;
};
ScoreRememberer_SanityChecker::ScoreRememberer_SanityChecker ()
{
}
ScoreRememberer_SanityChecker::~ScoreRememberer_SanityChecker ()
{
}
void ScoreRememberer_SanityChecker::Init (int minWindowLen,int maxWindowLen)
{
	dumb.Init(minWindowLen,maxWindowLen);
	heap.Init(minWindowLen,maxWindowLen);
	set.Init(minWindowLen,maxWindowLen);
	cleverCircular.Init(minWindowLen,maxWindowLen);
	cleverCircularWithLWL.Init(minWindowLen,maxWindowLen);
}
void ScoreRememberer_SanityChecker::ReInit (void)
{
	prevWindowLast=-1;
	dumb.ReInit();
	heap.ReInit();
	set.ReInit();
	cleverCircular.ReInit();
	cleverCircularWithLWL.ReInit();
}
void ScoreRememberer_SanityChecker::AddScore (int windowLast,float score)
{
	assert(windowLast>prevWindowLast);
	prevWindowLast=windowLast;

	dumb.AddScore(windowLast,score);
	heap.AddScore(windowLast,score);
	set.AddScore(windowLast,score);
	cleverCircular.AddScore(windowLast,score);
	cleverCircularWithLWL.AddScore(windowLast,score);
}
float ScoreRememberer_SanityChecker::GetBestScoreWithinWindow(int forWindowLast,int& get_leftWindowLast) const
{
	int leftWindowLast;
	int dumb_leftWindowLast,clever_leftWindowLast;
	float dumbScore=dumb.GetBestScoreWithinWindow(forWindowLast,dumb_leftWindowLast);
	float heapScore=heap.GetBestScoreWithinWindow(forWindowLast,leftWindowLast);
	float setScore=set.GetBestScoreWithinWindow(forWindowLast,leftWindowLast);
	float cleverScore=cleverCircular.GetBestScoreWithinWindow(forWindowLast,leftWindowLast);
	float cleverScoreLWL=cleverCircularWithLWL.GetBestScoreWithinWindow(forWindowLast,clever_leftWindowLast);
	assertr(dumbScore==heapScore);
	assertr(dumbScore==setScore);
	assertr(dumbScore==cleverScore);
	assertr(dumbScore==cleverScoreLWL);
	assertr(dumb_leftWindowLast==clever_leftWindowLast);
	return dumbScore;
}
bool ScoreRememberer_SanityChecker::IsInit (void) const
{
	assertr(dumb.IsInit()==heap.IsInit());
	assertr(dumb.IsInit()==set.IsInit());
	assertr(dumb.IsInit()==cleverCircular.IsInit());
	assertr(dumb.IsInit()==cleverCircularWithLWL.IsInit());
	return dumb.IsInit();
}


// the simplest possible way of handling local ends; allow them to be infinitely long, and store the best score
// so far.  This is probably a bad idea for scans, since it'd allow amazing scores with local ends of megabases
// long (_something's_ gonna match).  (Caveat: I haven't tested this assertion empirically.)
// However, it works well for my visions of optimizing local-HMMs using a FINITE (i.e. regular) Forward Alg
// In this case, the finite-length (probably the win len) prevents the problem with infinitely long local ends,
// and the math expressions are nicer with storing only the best score (depending on Summary operator)
template <class ScoreHelperType,class Hmm,class ClosedSemiring>
class ScanHmmAugmenter_LocalWithInfiniteEnds {
protected:

	struct InfoOnState {
		typename ScoreHelperType::ScoreType summaryScore;
	};
	typedef vector<InfoOnState> InfoOnStateVector;
	InfoOnStateVector infoOnStateVector;

	const Hmm& hmm;
public:
	ScanHmmAugmenter_LocalWithInfiniteEnds (const Hmm& hmm_,int maxWindowLen_,const char *rnaSequence,int rnaSequenceLen);
	~ScanHmmAugmenter_LocalWithInfiniteEnds ();
	void StartInterval (const Hmm& hmm,int firstWindowLast,int lastWindowLast);
	void NotifyOfScore (int windowLast,typename Hmm::State state,typename ScoreHelperType::ScoreType get_score);
	void ReNotifyOfScore_EndingWindowLast (int windowLast,typename Hmm::State state,typename ScoreHelperType::ScoreType get_score);
	typename ScoreHelperType::ScoreType GetBestScoreInWindow (typename Hmm::State state,int windowLast,int& get_leftWindowLast);
	inline static bool Is_ScanHmmAugmenter_GenericLocal (void) { return true; }

	inline void StartSequence (const Hmm& hmm,const char *rnaSequence) {}
	inline void EndInterval (void) {}
	inline void StartWindowLast (int windowLast) {}
	inline void EndWindowLast(int windowLast) {}
	inline void OverrideDynProgScore (ScanHmmAugmenter_OverrideDynProgScoreOperation& overrideDynProgScoreOperation,typename ScoreHelperType::ScoreType &get_score,int windowLast,typename Hmm::State state,ScoreCollector_NULL& scoreCollector,typename ScoreHelperType::ScoreType emitScore,const Hmm& hmm) {
		overrideDynProgScoreOperation=OverrideDynProgScore_DontChange;
	}
	inline bool OverrideFinalWindowLastScore (typename ScoreHelperType::ScoreType& get_score,int windowLast,typename ScoreHelperType::ScoreType score) { return false; }
};
template <class ScoreHelperType,class Hmm,class ClosedSemiring>
ScanHmmAugmenter_LocalWithInfiniteEnds<ScoreHelperType,Hmm,ClosedSemiring>::ScanHmmAugmenter_LocalWithInfiniteEnds (const Hmm& hmm_,int maxWindowLen_,const char *rnaSequence,int rnaSequenceLen)
: hmm(hmm_)
{
	//InfoOnState infoOnState;
	infoOnStateVector.resize(hmm.GetNumStates());
	for (typename Hmm::State state=hmm.GetFirstState(); state!=hmm.GetLastState(); state++) {
		for (int link=0; link<hmm.GetNumEndscLinksToLeft(state); link++) {
			const typename Hmm::State leftLinkState=hmm.GetEndscLinksToLeft_State(state,link);
			infoOnStateVector[leftLinkState].summaryScore=ClosedSemiring::GetSummaryOperatorIdentity();
		}
	}
}
template <class ScoreHelperType,class Hmm,class ClosedSemiring>
ScanHmmAugmenter_LocalWithInfiniteEnds<ScoreHelperType,Hmm,ClosedSemiring>::~ScanHmmAugmenter_LocalWithInfiniteEnds ()
{
}
template <class ScoreHelperType,class Hmm,class ClosedSemiring>
void ScanHmmAugmenter_LocalWithInfiniteEnds<ScoreHelperType,Hmm,ClosedSemiring>::StartInterval (const Hmm& hmm,int firstWindowLast,int lastWindowLast)
{
	for (int state=0; state<hmm.GetNumStates(); state++) {
		InfoOnState& infoOnState=infoOnStateVector[state];
		infoOnState.summaryScore=ClosedSemiring::GetSummaryOperatorIdentity();
	}
}
template <class ScoreHelperType,class Hmm,class ClosedSemiring>
void ScanHmmAugmenter_LocalWithInfiniteEnds<ScoreHelperType,Hmm,ClosedSemiring>::NotifyOfScore (int windowLast,typename Hmm::State state,typename ScoreHelperType::ScoreType score)
{
	if (hmm.IsLeftStateOfLocalEnd(state)) {
		if (!hmm.IsRightStateOfThisLocalEndEmitting(state)) {
			// only store immediately if the right state isn't emitting (for more explanation on the issue, see HmmType1::CheckThatLocalEndsHaveUniqueRightState in HmmType1.cpp)
			InfoOnState& infoOnState=infoOnStateVector[state];
			infoOnState.summaryScore=ClosedSemiring::SummaryOperator(ClosedSemiring::ExtensionOperator(infoOnState.summaryScore,ClosedSemiring::GetLocalEndSelfLoop(hmm)),score);
		}
	}
}
template <class ScoreHelperType,class Hmm,class ClosedSemiring>
void ScanHmmAugmenter_LocalWithInfiniteEnds<ScoreHelperType,Hmm,ClosedSemiring>::ReNotifyOfScore_EndingWindowLast (int windowLast,typename Hmm::State state,typename ScoreHelperType::ScoreType score)
{
	// this function is different from above; differs on right state emittingness

	if (hmm.IsLeftStateOfLocalEnd(state)) {
		if (hmm.IsRightStateOfThisLocalEndEmitting(state)) {
			// only store now if right state is emitting (for more explanation on the issue, see HmmType1::CheckThatLocalEndsHaveUniqueRightState in HmmType1.cpp)
			InfoOnState& infoOnState=infoOnStateVector[state];
			infoOnState.summaryScore=ClosedSemiring::SummaryOperator(ClosedSemiring::ExtensionOperator(infoOnState.summaryScore,ClosedSemiring::GetLocalEndSelfLoop(hmm)),score);
		}
	}
}
template <class ScoreHelperType,class Hmm,class ClosedSemiring>
typename ScoreHelperType::ScoreType ScanHmmAugmenter_LocalWithInfiniteEnds<ScoreHelperType,Hmm,ClosedSemiring>::GetBestScoreInWindow (typename Hmm::State state,int windowLast,int& get_leftWindowLast)
{
	assert(hmm.IsLeftStateOfLocalEnd(state)); // either (1) we're being called when we shouldn't, or (2) the HMM is lying about us not being the left state of a local end
	InfoOnState& infoOnState=infoOnStateVector[state];
	return infoOnState.summaryScore;
}

// handler for local ends that only allows \epsilon, i.e. local ends are immediate, and cannot match strings other than \epsilon.  This is wrong in general, but my nose for theory says it might perform better with --ml-sample-baum-welch
template <class ScoreHelperType,class Hmm,class ClosedSemiring>
class ScanHmmAugmenter_LocalWithEpsilonEnds {
protected:

	struct InfoOnState {
		typename ScoreHelperType::ScoreType summaryScore;
	};
	typedef vector<InfoOnState> InfoOnStateVector;
	InfoOnStateVector infoOnStateVector;

	const Hmm& hmm;
public:
	ScanHmmAugmenter_LocalWithEpsilonEnds (const Hmm& hmm_,int maxWindowLen_,const char *rnaSequence,int rnaSequenceLen);
	~ScanHmmAugmenter_LocalWithEpsilonEnds ();
	void StartInterval (const Hmm& hmm,int firstWindowLast,int lastWindowLast);
	void NotifyOfScore (int windowLast,typename Hmm::State state,typename ScoreHelperType::ScoreType get_score);
	void ReNotifyOfScore_EndingWindowLast (int windowLast,typename Hmm::State state,typename ScoreHelperType::ScoreType get_score);
	typename ScoreHelperType::ScoreType GetBestScoreInWindow (typename Hmm::State state,int windowLast,int& get_leftWindowLast);
	inline static bool Is_ScanHmmAugmenter_GenericLocal (void) { return true; }

	inline void StartSequence (const Hmm& hmm,const char *rnaSequence) {}
	inline void EndInterval (void) {}
	inline void StartWindowLast (int windowLast) {}
	inline void EndWindowLast(int windowLast) {}
	inline void OverrideDynProgScore (ScanHmmAugmenter_OverrideDynProgScoreOperation& overrideDynProgScoreOperation,typename ScoreHelperType::ScoreType &get_score,int windowLast,typename Hmm::State state,ScoreCollector_NULL& scoreCollector,typename ScoreHelperType::ScoreType emitScore,const Hmm& hmm) {
		overrideDynProgScoreOperation=OverrideDynProgScore_DontChange;
	}
	inline bool OverrideFinalWindowLastScore (typename ScoreHelperType::ScoreType& get_score,int windowLast,typename ScoreHelperType::ScoreType score) { return false; }
};
template <class ScoreHelperType,class Hmm,class ClosedSemiring>
ScanHmmAugmenter_LocalWithEpsilonEnds<ScoreHelperType,Hmm,ClosedSemiring>::ScanHmmAugmenter_LocalWithEpsilonEnds (const Hmm& hmm_,int maxWindowLen_,const char *rnaSequence,int rnaSequenceLen)
: hmm(hmm_)
{
	//InfoOnState infoOnState;
	infoOnStateVector.resize(hmm.GetNumStates());
	for (typename Hmm::State state=hmm.GetFirstState(); state!=hmm.GetLastState(); state++) {
		for (int link=0; link<hmm.GetNumEndscLinksToLeft(state); link++) {
			const typename Hmm::State leftLinkState=hmm.GetEndscLinksToLeft_State(state,link);
			infoOnStateVector[leftLinkState].summaryScore=ClosedSemiring::GetSummaryOperatorIdentity();
		}
	}
}
template <class ScoreHelperType,class Hmm,class ClosedSemiring>
ScanHmmAugmenter_LocalWithEpsilonEnds<ScoreHelperType,Hmm,ClosedSemiring>::~ScanHmmAugmenter_LocalWithEpsilonEnds ()
{
}
template <class ScoreHelperType,class Hmm,class ClosedSemiring>
void ScanHmmAugmenter_LocalWithEpsilonEnds<ScoreHelperType,Hmm,ClosedSemiring>::StartInterval (const Hmm& hmm,int firstWindowLast,int lastWindowLast)
{
	for (int state=0; state<hmm.GetNumStates(); state++) {
		InfoOnState& infoOnState=infoOnStateVector[state];
		infoOnState.summaryScore=ClosedSemiring::GetSummaryOperatorIdentity();
	}
}
template <class ScoreHelperType,class Hmm,class ClosedSemiring>
void ScanHmmAugmenter_LocalWithEpsilonEnds<ScoreHelperType,Hmm,ClosedSemiring>::NotifyOfScore (int windowLast,typename Hmm::State state,typename ScoreHelperType::ScoreType score)
{
	if (hmm.IsLeftStateOfLocalEnd(state)) {
		if (!hmm.IsRightStateOfThisLocalEndEmitting(state)) {
			// only store immediately if the right state isn't emitting (for more explanation on the issue, see HmmType1::CheckThatLocalEndsHaveUniqueRightState in HmmType1.cpp)
			InfoOnState& infoOnState=infoOnStateVector[state];
			infoOnState.summaryScore=score; // don't summarize -- only allow epsilon ends
		}
	}
}
template <class ScoreHelperType,class Hmm,class ClosedSemiring>
void ScanHmmAugmenter_LocalWithEpsilonEnds<ScoreHelperType,Hmm,ClosedSemiring>::ReNotifyOfScore_EndingWindowLast (int windowLast,typename Hmm::State state,typename ScoreHelperType::ScoreType score)
{
	if (hmm.IsLeftStateOfLocalEnd(state)) {
		if (hmm.IsRightStateOfThisLocalEndEmitting(state)) {
			// only store now if right state is emitting (for more explanation on the issue, see HmmType1::CheckThatLocalEndsHaveUniqueRightState in HmmType1.cpp)
			InfoOnState& infoOnState=infoOnStateVector[state];
			infoOnState.summaryScore=score; // don't summarize -- only allow epsilon ends
		}
	}
}
template <class ScoreHelperType,class Hmm,class ClosedSemiring>
typename ScoreHelperType::ScoreType ScanHmmAugmenter_LocalWithEpsilonEnds<ScoreHelperType,Hmm,ClosedSemiring>::GetBestScoreInWindow (typename Hmm::State state,int windowLast,int& get_leftWindowLast)
{
	assert(hmm.IsLeftStateOfLocalEnd(state)); // either (1) we're being called when we shouldn't, or (2) the HMM is lying about us not being the left state of a local end
	InfoOnState& infoOnState=infoOnStateVector[state];
	return infoOnState.summaryScore;
}

// implements aspects of the local aligning-HMM that are easiest to do outside.  So far, it remembers the windows of scores.
template <class ScoreHelperType,class Hmm,class ClosedSemiring,class ScoreRememberer=ScoreRememberer_CleverBinaryCircularArray> // ScoreRememberer_SanityChecker
class ScanHmmAugmenter_GenericLocal {
protected:

	struct InfoOnState {
		ScoreRememberer scoreRememberer;
	};
	typedef vector<InfoOnState> InfoOnStateVector;
	InfoOnStateVector infoOnStateVector;

	const Hmm& hmm;
	int maxWindowLen;

	// scanning state
	int currStartOfInterval;
	int currIntervalSize;
	const char *rnaSequence;
public:
	ScanHmmAugmenter_GenericLocal (const Hmm& hmm_,int maxWindowLen_,const char *rnaSequence,int rnaSequenceLen);
	~ScanHmmAugmenter_GenericLocal ();
	void StartSequence (const Hmm& hmm,const char *rnaSequence);
	void StartInterval (const Hmm& hmm,int firstWindowLast,int lastWindowLast);
	void EndInterval (void);
	void StartWindowLast (int windowLast);
	//void OverrideDynProgScore (ScanHmmAugmenter_OverrideDynProgScoreOperation& overrideDynProgScoreOperation,typename ScoreHelperType::ScoreType &get_score,int windowLast,typename Hmm::State state,ScoreCollector& scoreCollector,typename ScoreHelperType::ScoreType emitScore,const Hmm& hmm);
	inline void EndWindowLast(int windowLast) {}
	inline void OverrideDynProgScore (ScanHmmAugmenter_OverrideDynProgScoreOperation& overrideDynProgScoreOperation,typename ScoreHelperType::ScoreType &get_score,int windowLast,typename HmmType1::State state,ScoreCollector_DeepNULL& scoreCollector,typename ScoreHelperType::ScoreType emitScore,const HmmType1& hmmWithPenalties) {
		overrideDynProgScoreOperation=OverrideDynProgScore_DontChange;
	}
	inline void OverrideDynProgScore (ScanHmmAugmenter_OverrideDynProgScoreOperation& overrideDynProgScoreOperation,typename ScoreHelperType::ScoreType &get_score,int windowLast,typename HmmType1::State state,ScoreCollector_Vectorize& scoreCollector,typename ScoreHelperType::ScoreType emitScore,const HmmType1& hmmWithPenalties) {
		overrideDynProgScoreOperation=OverrideDynProgScore_DontChange;
	}
	void NotifyOfScore (int windowLast,typename Hmm::State state,typename ScoreHelperType::ScoreType get_score);
	void ReNotifyOfScore_EndingWindowLast (int windowLast,typename Hmm::State state,typename ScoreHelperType::ScoreType get_score);

	typename ScoreHelperType::ScoreType GetBestScoreInWindow (typename Hmm::State state,int windowLast,int& get_leftWindowLast);
	inline static bool Is_ScanHmmAugmenter_GenericLocal (void) { return true; }
	inline bool OverrideFinalWindowLastScore (typename ScoreHelperType::ScoreType& get_score,int windowLast,typename ScoreHelperType::ScoreType score) { return false; }
};
template <class ScoreHelperType,class Hmm,class ClosedSemiring,class ScoreRememberer>
ScanHmmAugmenter_GenericLocal<ScoreHelperType,Hmm,ClosedSemiring,ScoreRememberer>::ScanHmmAugmenter_GenericLocal (const Hmm& hmm_,int maxWindowLen_,const char *rnaSequence,int rnaSequenceLen)
: hmm(hmm_)
{
	maxWindowLen=maxWindowLen_;

	//InfoOnState infoOnState;
	infoOnStateVector.resize(hmm.GetNumStates());
	for (typename Hmm::State state=hmm.GetFirstState(); state!=hmm.GetLastState(); state++) {
		for (int link=0; link<hmm.GetNumEndscLinksToLeft(state); link++) {
			const typename Hmm::State leftLinkState=hmm.GetEndscLinksToLeft_State(state,link);
			if (!infoOnStateVector[leftLinkState].scoreRememberer.IsInit()) {
				infoOnStateVector[leftLinkState].scoreRememberer.Init(0,maxWindowLen);
			}
		}
	}
}
template <class ScoreHelperType,class Hmm,class ClosedSemiring,class ScoreRememberer>
ScanHmmAugmenter_GenericLocal<ScoreHelperType,Hmm,ClosedSemiring,ScoreRememberer>::~ScanHmmAugmenter_GenericLocal ()
{
}
template <class ScoreHelperType,class Hmm,class ClosedSemiring,class ScoreRememberer>
void ScanHmmAugmenter_GenericLocal<ScoreHelperType,Hmm,ClosedSemiring,ScoreRememberer>::StartSequence (const Hmm& hmm,const char *rnaSequence_)
{
	rnaSequence=rnaSequence_;
}
template <class ScoreHelperType,class Hmm,class ClosedSemiring,class ScoreRememberer>
void ScanHmmAugmenter_GenericLocal<ScoreHelperType,Hmm,ClosedSemiring,ScoreRememberer>::StartInterval (const Hmm& hmm,int firstWindowLast,int lastWindowLast)
{
	currStartOfInterval=firstWindowLast;
	currIntervalSize=lastWindowLast-firstWindowLast;

	for (int state=0; state<hmm.GetNumStates(); state++) {
		InfoOnState& infoOnState=infoOnStateVector[state];
		if (infoOnState.scoreRememberer.IsInit()) {
			infoOnState.scoreRememberer.ReInit();
		}
	}
}
template <class ScoreHelperType,class Hmm,class ClosedSemiring,class ScoreRememberer>
void ScanHmmAugmenter_GenericLocal<ScoreHelperType,Hmm,ClosedSemiring,ScoreRememberer>::EndInterval (void)
{
	// nothing to do
}
template <class ScoreHelperType,class Hmm,class ClosedSemiring,class ScoreRememberer>
void ScanHmmAugmenter_GenericLocal<ScoreHelperType,Hmm,ClosedSemiring,ScoreRememberer>::StartWindowLast (int windowLast)
{
	// nichts zu tun
}
template <class ScoreHelperType,class Hmm,class ClosedSemiring,class ScoreRememberer>
void ScanHmmAugmenter_GenericLocal<ScoreHelperType,Hmm,ClosedSemiring,ScoreRememberer>::NotifyOfScore (int windowLast,typename Hmm::State state,typename ScoreHelperType::ScoreType score)
{
	windowLast -= currStartOfInterval;

	if (hmm.IsLeftStateOfLocalEnd(state)) {
		if (!hmm.IsRightStateOfThisLocalEndEmitting(state)) {
			// only store immediately if right state is NOT emitting (for more explanation on the issue, see HmmType1::CheckThatLocalEndsHaveUniqueRightState in HmmType1.cpp)
			InfoOnState& infoOnState=infoOnStateVector[state];
			if (infoOnState.scoreRememberer.IsInit()) {
				infoOnState.scoreRememberer.AddScore(windowLast,score);
			}
		}
	}
}
template <class ScoreHelperType,class Hmm,class ClosedSemiring,class ScoreRememberer>
void ScanHmmAugmenter_GenericLocal<ScoreHelperType,Hmm,ClosedSemiring,ScoreRememberer>::ReNotifyOfScore_EndingWindowLast (int windowLast,typename Hmm::State state,typename ScoreHelperType::ScoreType score)
{
	windowLast -= currStartOfInterval;

	if (hmm.IsLeftStateOfLocalEnd(state)) {
		if (hmm.IsRightStateOfThisLocalEndEmitting(state)) {
			// only store now if right state is emitting (for more explanation on the issue, see HmmType1::CheckThatLocalEndsHaveUniqueRightState in HmmType1.cpp)
			InfoOnState& infoOnState=infoOnStateVector[state];
			if (infoOnState.scoreRememberer.IsInit()) {
				infoOnState.scoreRememberer.AddScore(windowLast,score);
			}
		}
	}
}
template <class ScoreHelperType,class Hmm,class ClosedSemiring,class ScoreRememberer>
typename ScoreHelperType::ScoreType ScanHmmAugmenter_GenericLocal<ScoreHelperType,Hmm,ClosedSemiring,ScoreRememberer>::GetBestScoreInWindow (typename Hmm::State state,int windowLast,int& get_leftWindowLast)
{
	windowLast -= currStartOfInterval;

	InfoOnState& infoOnState=infoOnStateVector[state];
	assert(infoOnState.scoreRememberer.IsInit());  // if we're being asked, we should have anticipated that; if we didn't, something's wrong.
	assert(hmm.IsLeftStateOfLocalEnd(state));  // and the HMM should know if we were supposed to be initialized

	return infoOnState.scoreRememberer.GetBestScoreWithinWindow(windowLast,get_leftWindowLast);
}

////////////////////////////////
//
//   Score collectors, for use as template params of ClassScanHmm
//

////////////////////////////
// TracebackScoreCollector

int TracebackInfo::GetSlidingWindowSize (void) const
{
	return slidingWindowSize;
}
TracebackInfo::TracebackInfo (const HmmType1& hmm_,int windowLen_,int slidingWindowSize_)
: hmm(hmm_)
{
	windowLen=windowLen_;
	slidingWindowSize=slidingWindowSize_+1;
	if (slidingWindowSize_==-1) {
		slidingWindowSize=windowLen*16; // probably overly generous, but I don't feel like implementing something to expand it on demand
	}
	scores.resize(slidingWindowSize);
	decisionType.resize(slidingWindowSize,hmm.GetNumStates());
	transitionTaken.resize(slidingWindowSize,hmm.GetNumStates());
	leftWindowLastForEndsc.resize(slidingWindowSize,hmm.GetNumStates());
	emissionDone.resize(slidingWindowSize);
	visitsToState.resize(hmm.GetNumStates());
}
TracebackInfo::~TracebackInfo ()
{
}
TracebackInfo::TracebackInfo (const TracebackInfo& t)
: hmm(t.hmm)
{
	CopyFrom(t);
}
TracebackInfo::Decision TracebackInfo::GetDecision(int windowLast,HmmType1::State state)
{
	Decision decisionHere=decisionType[windowLast%slidingWindowSize][state];
	return decisionHere;
}
int TracebackInfo::GetTransition (int windowLast,HmmType1::State state)
{
	int childNum=transitionTaken[windowLast%slidingWindowSize][state];
	return childNum;
}
int TracebackInfo::GetEmission (int windowLast)
{
	return emissionDone[windowLast%slidingWindowSize];
}
void TracebackInfo::CopyFrom (const TracebackInfo& t)
{
	windowLen=t.windowLen;
	slidingWindowSize=t.slidingWindowSize;
	decisionType=t.decisionType;
	transitionTaken=t.transitionTaken;
	emissionDone=t.emissionDone;
	visitsToState=t.visitsToState;
	scores=t.scores;
	leftWindowLastForEndsc=t.leftWindowLastForEndsc;
}
const HmmType1& TracebackInfo::GetHmm (void) const
{
	return hmm;
}
float TracebackInfo::GetScore (int windowLast)
{
	return scores[windowLast%slidingWindowSize];
}

template <class TransitionCounter>
class TracebackScoreCollector_templ : public TracebackInfo {
protected:
	TransitionCounter& transitionCounter;
	double totalScore;
	int numSamples;
	int scanningWindowLen;
public:
	TracebackScoreCollector_templ (const HmmType1& _hmm,TransitionCounter& _transitionCounter,int windowLen,int scanningWindowLen_);
	TracebackScoreCollector_templ (const TracebackInfo& t,TransitionCounter& transitionCounter_,int windowLen_,int scanningWindowLen_);
	~TracebackScoreCollector_templ ();

	void GotScore (int windowLast,int topLevelRightwardState,float score);
	void DecidedTransition(int windowLast,int state,int bestChildNum);
	void DidEmission (int windowLast,int nuc);
	void UsedLeftwardBeginsc (int windowLast,int state);
	void UsedEndsc (int windowLast,int state,int linkNum,int leftWindowLast);
	void UsedStartState (int windowLast,int state);

	double GetAvgScore (void);
	double GetTotalScore (void);
	void StartInterval (int firstWindowLast,int lastWindowLast);
};
template <class TransitionCounter>
TracebackScoreCollector_templ<TransitionCounter>::TracebackScoreCollector_templ (const HmmType1& hmm_,TransitionCounter& _transitionCounter,int windowLen_,int scanningWindowLen_)
: TracebackInfo(hmm_,windowLen_)
, transitionCounter(_transitionCounter)
{
	scanningWindowLen=scanningWindowLen_;
	totalScore=0;
	numSamples=0;
}
template <class TransitionCounter>
TracebackScoreCollector_templ<TransitionCounter>::TracebackScoreCollector_templ (const TracebackInfo& t,TransitionCounter& transitionCounter_,int windowLen_,int scanningWindowLen_)
: TracebackInfo(t)
, transitionCounter(transitionCounter_)
{
	windowLen=windowLen_;
	scanningWindowLen=scanningWindowLen_;
}
template <class TransitionCounter>
TracebackScoreCollector_templ<TransitionCounter>::~TracebackScoreCollector_templ ()
{
}
template <class TransitionCounter>
double TracebackScoreCollector_templ<TransitionCounter>::GetAvgScore (void)
{
	return totalScore/(double)(numSamples);
}
template <class TransitionCounter>
double TracebackScoreCollector_templ<TransitionCounter>::GetTotalScore (void)
{
	return totalScore;
}
template <class TransitionCounter>
void TracebackScoreCollector_templ<TransitionCounter>::StartInterval (int firstWindowLast,int lastWindowLast)
{
	transitionCounter.StartInterval(firstWindowLast,lastWindowLast);
}
template <class TransitionCounter>
void TracebackScoreCollector_templ<TransitionCounter>::GotScore (int windowLast,int topLevelRightwardState,float score)
{
	scores[windowLast%slidingWindowSize]=score;

	if (windowLast<scanningWindowLen && transitionCounter.IgnoreWindowLastBelowWindowLen()) {
		// don't take this seriously -- it doesn't happen often & likely represents non-typical alignments
		return;
	}


	// do traceback
	transitionCounter.AddSample(windowLast,score);
	totalScore += score;
	numSamples++;

	int numTransitionsInParse=0;
	if (score==(float)IMPOSSIBLE) {
		// in this case, the parses can be weird, so to avoid crashing, just return early
		transitionCounter.AddNumTransitionsInParse(numTransitionsInParse);
		transitionCounter.DoneSample();
		return;
	}

	HmmType1::State state;
	state=topLevelRightwardState;
	transitionCounter.AddRightBegin(state);
	int currWindowLast=windowLast;
	while (1) {

		transitionCounter.VisitingState(state);

		assert(currWindowLast<=windowLast);
		assert(windowLast-currWindowLast<slidingWindowSize);
		if (windowLast-currWindowLast>=slidingWindowSize) {
			assert(false);
			throw SimpleStringException("Internal error %s:%d: exceeded max slidingWindowSize in TracebackScoreCollector, i.e. while doing HMM alignment traceback.  Did you enter the correct windowLen parameter for the given Covariance Model?  I made that mistake once too.  It crashed my program.",__FILE__,__LINE__);
		}

		visitsToState[state]++;

		int nextWindowLast=currWindowLast;
		if (hmm.IsEmittingState(state)) {
			transitionCounter.AddEmit(state,currWindowLast,emissionDone[currWindowLast%slidingWindowSize]);
			nextWindowLast=currWindowLast-1;
			assert(currWindowLast>=0);
		}
		
		Decision decisionHere=decisionType[currWindowLast%slidingWindowSize][state];
		if (decisionHere==StartState || decisionHere==Beginsc) {
			if (decisionHere==Beginsc) {
				transitionCounter.AddLeftBegin(state);
			}
			break;
		}
		else {
			if (decisionHere==Transition) {
				int childNum=transitionTaken[currWindowLast%slidingWindowSize][state];
				assert(childNum>=0 && childNum<hmm.GetNumChildren(state));
				transitionCounter.AddTransition(state,childNum);
				numTransitionsInParse++;

				state=hmm.GetNthChildState(state,childNum);
			}
			else {
				assert(decisionHere==Endsc);
				assert(hmm.DoLocal()); // else Endsc is invalid
				int linkNum=transitionTaken[currWindowLast%slidingWindowSize][state];
				assert(linkNum>=0 && linkNum<hmm.GetNumEndscLinksToLeft(state));
				nextWindowLast=leftWindowLastForEndsc[currWindowLast%slidingWindowSize][state];
				transitionCounter.AddEnd(state,linkNum,currWindowLast-nextWindowLast);
				state=hmm.GetEndscLinksToLeft_State(state,linkNum);
			}
		}

		transitionCounter.DoneVisitingState();
		currWindowLast=nextWindowLast;
	}
	transitionCounter.DoneVisitingState();

	transitionCounter.AddNumTransitionsInParse(numTransitionsInParse);
	transitionCounter.DoneSample();
}
template <class TransitionCounter>
void TracebackScoreCollector_templ<TransitionCounter>::UsedStartState (int windowLast,int state)
{
	decisionType[windowLast%slidingWindowSize][state]=StartState;
}
template <class TransitionCounter>
void TracebackScoreCollector_templ<TransitionCounter>::UsedLeftwardBeginsc (int windowLast,int state)
{
	decisionType[windowLast%slidingWindowSize][state]=Beginsc;
}
template <class TransitionCounter>
void TracebackScoreCollector_templ<TransitionCounter>::UsedEndsc (int windowLast,int state,int linkNum,int leftWindowLast)
{
	decisionType[windowLast%slidingWindowSize][state]=Endsc;
	transitionTaken[windowLast%slidingWindowSize][state]=linkNum;
	leftWindowLastForEndsc[windowLast%slidingWindowSize][state]=leftWindowLast;
}
template <class TransitionCounter>
void TracebackScoreCollector_templ<TransitionCounter>::DecidedTransition(int windowLast,int state,int bestChildNum)
{
	decisionType[windowLast%slidingWindowSize][state]=Transition;
	transitionTaken[windowLast%slidingWindowSize][state]=bestChildNum;
}
template <class TransitionCounter>
void TracebackScoreCollector_templ<TransitionCounter>::DidEmission (int windowLast,int nuc)
{
	if (nuc>=0 && nuc<Alphabet_size) {
		emissionDone[windowLast%slidingWindowSize]=nuc;
	}
	else {
		emissionDone[windowLast%slidingWindowSize]=windowLast%4; // who cares -- we shouldn't be training on this kind of data anyway.  windowLast%4 should spread things out, but I don't really care how well.
	}
}
typedef TracebackScoreCollector_templ<TransitionCounter> TracebackScoreCollector;


// EmissionHandler customized for the finite-length Forward Alg, summed over sequences from a 0th-order model -- which means we ignore what the nuc was
template <class ScoreHelperType,class Hmm,class ClosedSemiring>
class EmissionHandlerForward {
public:
	MarkovModelStats& markovModelStats;
	EmissionHandlerForward (MarkovModelStats& markovModelStats_)
		: markovModelStats(markovModelStats_)
	{
	}
	class EmissionType_Simple {
	protected:
		MarkovModelStats& markovModelStats;
	public:
		inline EmissionType_Simple (int nuc,const EmissionHandlerForward& emissionHandler) 
			: markovModelStats(emissionHandler.markovModelStats)
		{
			// nuc doesn't matter
		}
		inline bool IsEmitOkay (void) const {
			return true;
		}
		inline typename ScoreHelperType::ScoreType GetEmissionScore (const Hmm& hmm,typename Hmm::State state,int windowLast) const {
			
			assert(hmm.IsEmittingState(state)); // else why were we called

			typename ScoreHelperType::ScoreType emitProb=0.0;
			for (int nuc=0; nuc<Alphabet_size; nuc++) {
				emitProb += markovModelStats.GetProbOfNuc_0order(nuc) * hmm.GetSingletEmissionProb(state,nuc);
			}
			return emitProb;
		}
	};

	typedef EmissionType_Simple EmissionType_Degenerate; // we don't care
};
class SymbolicScoreHelper {
public:
	typedef SymbolicMath::Expression ScoreType;

	inline static ScoreType FloatToScoreType (float t) {
		throw SimpleStringException("meaningless function for symbolic -- hope it's not called. %s:%d",__FILE__,__LINE__);
	}
	inline static float ScoreTypeToFloat (ScoreType t) {
		throw SimpleStringException("meaningless function for symbolic -- hope it's not called. %s:%d",__FILE__,__LINE__);
	}

	inline static bool IsSymbolic (void) {
		return true;
	}
};
// usual EmissionHandler for Viterbi scans
template <class ScoreHelperType,class Hmm,class ClosedSemiring>
class EmissionHandlerGeneric {
public:
	// handles emissions in cases where it's just a simple nuc
	class EmissionType_Simple {
	protected:
		int nuc;
	public:
		inline EmissionType_Simple (int _nuc,const EmissionHandlerGeneric& emissionHandler) {
			nuc=_nuc;
			assert(nuc>=0 && nuc<Alphabet_size); // else it's degenerate
		}
		inline bool IsEmitOkay (void) const {
			return true;
		}
		inline typename ScoreHelperType::ScoreType GetEmissionScore (const Hmm& hmm,typename Hmm::State state,int windowLast) const {
			assert(hmm.IsEmittingState(state)); // else why were we called
			return ClosedSemiring::GetSingletEmission(hmm,state,nuc);
		}
	};
	// handles emissions in cases where it's a degenerate nuc
	class EmissionType_Degenerate {
	protected:
		int degenerateNuc;
	public:
		inline EmissionType_Degenerate (int _degenerateNuc,const EmissionHandlerGeneric& emissionHandler) {
			degenerateNuc=_degenerateNuc;
			assert(degenerateNuc>=Alphabet_size); // else it's simple, and you should have used the more efficient class
		}
		inline bool IsEmitOkay (void) const {
			return true;
		}
		inline typename ScoreHelperType::ScoreType GetEmissionScore (const Hmm& hmm,typename Hmm::State state,int windowLast) const {
			assert(hmm.IsEmittingState(state)); // else why were we called

			//  NOTE: I believe this function still makes sense for the Forward alg, even though in terms of actual probabilities this becomes the geometric mean.  In terms of expected score, I think the arithmetic mean is still justified

			// with exploded hmms (store-pair technique, using class HmmWithPenalties), some emissions are considered IMPOSSIBLE.  With these states, it's not so easy to marginalize scores the way that Infernal does.  So, if I detect a state with an IMPOSSIBLE emit score, I revert to just taking the maximum score.  Otherwise, I marginalize like Infernal

			bool hasImpossibleScore=false;
			for (int concreteNuc=0; concreteNuc<Alphabet_size; concreteNuc++) {
				float score=ClosedSemiring::GetSingletEmission(hmm,state,concreteNuc);
				if (score == (float)IMPOSSIBLE) {
					hasImpossibleScore=true;
					break;
				}
				else {
					assert(score>-100); // if it's less than -100, but not IMPOSSIBLE, this is pretty suspicious, so I'd like to flag it.  My concern is that there's some other score being used as impossible that this code didn't expect.  Strictly speaking, though, a score <-100 is valid, but there'd have to be a vast number of training examples, or a weak prior.
				}
			}

			if (hasImpossibleScore) {
				// take max
				float emitScore=-FLT_MAX;
				for (int concreteNuc = 0; concreteNuc <Alphabet_size; concreteNuc++) {
					if (Degenerate[degenerateNuc][concreteNuc]) {
						emitScore=std::max(emitScore,ClosedSemiring::GetSingletEmission(hmm,state,concreteNuc));
					}
				}
				return ScoreHelperType::FloatToScoreType(emitScore);
			}
			else {
				// adapted from alphabet.c
				float emitScore=0;
				assert(Alphabet_size==4);
				for (int concreteNuc = 0; concreteNuc <Alphabet_size; concreteNuc++) {

					if (Degenerate[degenerateNuc][concreteNuc]) {
						emitScore += ClosedSemiring::GetSingletEmission(hmm,state,concreteNuc) / (float)(DegenCount[degenerateNuc]);
					}
				}

				return ScoreHelperType::FloatToScoreType(emitScore);
			}
		}
	};
};

// for HMMs that want to know if they emitted a particular state at a particular position, for my HmmPair scoring of CMfinder motifs
// (see AlignmentConsensusAndScoring.cpp)
template <class ScoreHelperType,class Hmm,class ClosedSemiring>
class EmissionHandler_RestrictPositions {
	struct StateAndPositionConstraint {
		typename Hmm::State state;
		int windowLast;
	};
	StateAndPositionConstraint stateMustEmitHere1,stateMustEmitHere2;
public:
	EmissionHandler_RestrictPositions () {
		stateMustEmitHere1.state=-1; // neuter the check
		stateMustEmitHere2.state=-1;
	}
	void SetStateMustEmitHere (typename Hmm::State state1,int windowLast1,typename Hmm::State state2,int windowLast2) {
		assertr(state1!=state2 || state1==-1); // impossible&redundant, also weird, also I can't think of a situation in which this is useful
		stateMustEmitHere1.state=state1;
		stateMustEmitHere1.windowLast=windowLast1;
		stateMustEmitHere2.state=state2;
		stateMustEmitHere2.windowLast=windowLast2;
	}

	EmissionHandlerGeneric<ScoreHelperType,Hmm,ClosedSemiring> subHandler;
	inline bool ViolatesConstraint (typename Hmm::State state,int windowLast) const {
		if (state==stateMustEmitHere1.state) {
			return windowLast!=stateMustEmitHere1.windowLast;
		}
		if (state==stateMustEmitHere2.state) {
			return windowLast!=stateMustEmitHere2.windowLast;
		}
		return false;
	}

	class EmissionType_Simple {
	protected:
		int nuc;
		const EmissionHandler_RestrictPositions& emissionHandler;
	public:
		inline EmissionType_Simple (int nuc_,const EmissionHandler_RestrictPositions& emissionHandler_)
			: emissionHandler(emissionHandler_)
		{
			nuc=nuc_;
			assert(nuc>=0 && nuc<Alphabet_size); // else it's degenerate
		}
		inline bool IsEmitOkay (void) const {
			return true;
		}
		inline typename ScoreHelperType::ScoreType GetEmissionScore (const Hmm& hmm,typename Hmm::State state,int windowLast) const {
			assert(hmm.IsEmittingState(state)); // else why were we called
			if (emissionHandler.ViolatesConstraint(state,windowLast)) {
				return ClosedSemiring::GetImpossible();
			}
			else {
				return ClosedSemiring::GetSingletEmission(hmm,state,nuc);
			}
		}
	};
	class EmissionType_Degenerate {
	protected:
		int degenerateNuc;
		const EmissionHandler_RestrictPositions& emissionHandler;
		typename EmissionHandlerGeneric<ScoreHelperType,Hmm,ClosedSemiring>::EmissionType_Degenerate subDegenerate;
	public:
		inline EmissionType_Degenerate (int _degenerateNuc,const EmissionHandler_RestrictPositions& emissionHandler_) 
			: emissionHandler(emissionHandler_),subDegenerate(_degenerateNuc,emissionHandler_.subHandler)
		{
			degenerateNuc=_degenerateNuc;
			assert(degenerateNuc>=Alphabet_size); // else it's simple, and you should have used the more efficient class
		}
		inline bool IsEmitOkay (void) const {
			return true;
		}
		inline typename ScoreHelperType::ScoreType GetEmissionScore (const Hmm& hmm,typename Hmm::State state,int windowLast) const {
			assert(hmm.IsEmittingState(state)); // else why were we called
			if (emissionHandler.ViolatesConstraint(state,windowLast)) {
				return ClosedSemiring::GetImpossible();
			}
			else {
				return subDegenerate.GetEmissionScore(hmm,state,windowLast);
			}
		}
	};
};


///////////////
// ClassScanHmm: various HMM algs
 // put doLocal as template param because I think it might make the compiler more aggressive.  Similarly about worryAboutStartStateInTransitions

template <class ScoreHelperType,class Hmm,class ScoreCollector,class ClosedSemiring,class ScanHmmAugmenter,bool doLocal,class EmissionHandler=EmissionHandlerGeneric<ScoreHelperType,Hmm,ClosedSemiring> >
class ClassScanHmm {
private:
	//typedef FloatScoreHelper ScoreHelperType;
protected:
	class CmHmmToFastHmm : public std::map<int,Hmm *> {
	public:
		~CmHmmToFastHmm () {
			typename std::map<int,Hmm *>::iterator i;
			for (i=this->begin(); i!=this->end(); i++) {
				delete i->second;
			}
		}
	};
	static CmHmmToFastHmm cmHmmToFastHmm;

	// handles emissions in case where no emissions are allowed (for initializing the score tables)
	class EmissionType_NoEmission {
	public:
		inline bool IsEmitOkay (void) const {
			return false;
		}
		inline typename ScoreHelperType::ScoreType GetEmissionScore (const Hmm& hmm,typename Hmm::State state,int windowLast) const {
			assert(hmm.IsEmittingState(state)); // else why were we called

			// when no emission is allowed, that's it
			return ClosedSemiring::GetImpossible();
		}
	};

	template <bool worryAboutStartStateInTransitions>
		inline static typename ScoreHelperType::ScoreType /* best emitless score*/ FindBestTransitionScoreEmitless (int& get_bestTransitionChildNum,const Hmm& hmm,vector<typename ScoreHelperType::ScoreType> *tableForChildState,const typename Hmm::State state,const typename Hmm::State startStateForScanning)
	{
		typename ScoreHelperType::ScoreType bestTransitionScoreEmitless=ClosedSemiring::GetImpossible();
		int bestTransitionChildNum=-1;
		if (doLocal && hmm.GetNumChildren(state)==0) { // this case should only happen in doLocal mode.  I put doLocal in the test s.t. the compiler could kill this code if doLocal==false
			// the first node's states don't have children, so we're forced to use leftwardBeginsc, which will happen later
		}
		else {
			if (hmm.GetNumChildren(state)==0) {
				// now I'm allowing this, since it allows me to avoid an explicit case for the start state, which was getting complicated
				bestTransitionScoreEmitless=ClosedSemiring::GetExtensionOperatorIdentity();
				bestTransitionChildNum=-1;
			}
			else {

				bestTransitionScoreEmitless=ClosedSemiring::GetExtraImpossible(); // make sure some transition is picked, even if they're all impossible
				for (int childNum=0; childNum<hmm.GetNumChildren(state); childNum++) {
					const typename Hmm::State childState=hmm.GetNthChildState(state,childNum);
					if (worryAboutStartStateInTransitions && childState<startStateForScanning) {
						// skip this state
					}
					else {
						assert(childState>=startStateForScanning);
						const typename ScoreHelperType::ScoreType tsc=ClosedSemiring::GetNthChildTransition(hmm,state,childNum);
						assert(childState<=state);
						const typename ScoreHelperType::ScoreType transitionScoreEmitless=ClosedSemiring::ExtensionOperator(tsc,(*tableForChildState)[childState]);
						ClosedSemiring::SummaryOperatorWithTraceback(bestTransitionScoreEmitless,transitionScoreEmitless,bestTransitionChildNum,childNum);
					}
				}
				if (ClosedSemiring::SummaryOperatorLikesTraceback()) {
					assert((bestTransitionChildNum>=0 && bestTransitionChildNum<hmm.GetNumChildren(state)) || startStateForScanning!=0); // if startStateForScanning!=0 (not the actual start state), then this condition doesn't necessarily hold, and I don't really care.
				}
			}
		}

		if (worryAboutStartStateInTransitions) {
			if (bestTransitionChildNum==-1 && state==startStateForScanning) {
				// in this case, it's 0, and we're allowed to start here (kind of)
				bestTransitionScoreEmitless=ClosedSemiring::GetExtensionOperatorIdentity();
			}
		}

		get_bestTransitionChildNum=bestTransitionChildNum;
		return bestTransitionScoreEmitless;
	}

	template <class EmissionType,bool worryAboutStartStateInTransitions>
		inline static typename ScoreHelperType::ScoreType ScanHmmOnePos_CommitDoLocal(const Hmm& hmm,ScoreCollector& scoreCollector,ScanHmmAugmenter& scanHmmAugmenter,vector<typename ScoreHelperType::ScoreType> *prevTable,vector<typename ScoreHelperType::ScoreType> *currTable,const int pos,const EmissionType& emissionType,const typename Hmm::State startState,const typename Hmm::State endState,CykscanStats& cykscanStats,bool startOnlyAtFirstPos)
	{
		const int windowLast=pos+1;

		assert(doLocal==hmm.DoLocal()); // make sure the caller got it right

		//typename Hmm::State topLevelState=Hmm::GetInvalidState();

		typename ScoreHelperType::ScoreType bestScoreAtPos=ClosedSemiring::GetImpossible();
		typename Hmm::State bestTopLevelState=Hmm::GetInvalidState();

		/*  There are 4 cases to worry about:
		- the score is the cost of the emission (if any) PLUS the cost of the maximal tsc+child.  For doLocal==false, this is the only possibility
		- (doLocal only) cost is leftwardBeginsc PLUS emit cost, as if we started here.  WARNING: we must include the emit cost since that is part of the cost of the state in alpha(state,windowLast,windowLen) in the CYK alg.
		- (doLocal only) cost is the maximum over all link-to-lefts of endsc plus the cost of the left state
		- (doLocal only) I'm the top-level root, and my cost is rightwardBeginsc PLUS my best score otherwise
		*/

		if (cykscanStats.collectHmmFullDynProgTable) {
			cykscanStats.fullHmmDynProgTable[windowLast].resize(hmm.GetNumStates());
		}
		scoreCollector.NotifyNumStates(windowLast,hmm.GetNumStates());

		typename Hmm::State state=startState;
		if (doLocal) {

			assertr(startState==hmm.GetFirstState()); // I'm not sure how doLocal interacts with the case where startState!=hmm.GetFirstState() or endState+1!=hmm.GetLastState(), or with scanHmmAugmenter, since the local begins are weird.  Currently I don't think this needs to be dealt with
		}

		// no need for a separate initialization of the start state -- it's just like the other states, except it has no children, and I'm not enthusiastic about maintaining this code
		scoreCollector.UsedStartState(windowLast,state); // do have to notify the ScoreCollector
		state--; // so the loop will make it to 0

		while (1) {

			state++;
			if (state==endState+1) {
				break;
			}

			typename ScoreHelperType::ScoreType emitScore;
			vector<typename ScoreHelperType::ScoreType> *tableForChildState;
			if (hmm.IsEmittingState(state)) {
				emitScore=emissionType.GetEmissionScore(hmm,state,windowLast);
				tableForChildState=prevTable;
			}
			else {
				emitScore=ClosedSemiring::GetExtensionOperatorIdentity();
				tableForChildState=currTable;
			}
			//don't feel like dealing with the errors this generates for the Symbolic case: assert (ClosedSemiring::IsViterbi() || emitScore>=0); // forward alg scores must be positive since they're probabilities, although I allow ones >1

			if (hmm.IsEmittingState(state) && !emissionType.IsEmitOkay()) {
				assert(ClosedSemiring::IsImpossible(emitScore) || ScoreHelperType::IsSymbolic());
			}

			ScanHmmAugmenter_OverrideDynProgScoreOperation overrideDynProgScoreOperation;
			typename ScoreHelperType::ScoreType scoreFromAugmenter;
			scanHmmAugmenter.OverrideDynProgScore(overrideDynProgScoreOperation,scoreFromAugmenter,windowLast,state,scoreCollector,emitScore,hmm);
			if (LeaveOutEmitScoreButMoveWithEmit(overrideDynProgScoreOperation)) {
				emitScore=ClosedSemiring::GetExtensionOperatorIdentity();
			}
			if (!doLocal && startOnlyAtFirstPos && state==hmm.GetFirstState() && windowLast!=0) {
				overrideDynProgScoreOperation=OverrideDynProgScore_SetScore;
				scoreFromAugmenter=ClosedSemiring::GetImpossible();
			}
			if (overrideDynProgScoreOperation==OverrideDynProgScore_AddToScore_And_NoEmit) {
				emitScore=ClosedSemiring::GetExtensionOperatorIdentity();
				tableForChildState=currTable;
			}

			if (overrideDynProgScoreOperation==OverrideDynProgScore_SetScore) {
				(*currTable)[state]=scoreFromAugmenter;
			}
			else {
				if ((hmm.IsEmittingState(state) && !emissionType.IsEmitOkay()) || ClosedSemiring::IsImpossible(emitScore)) { // important to explicitly check if emitScore is impossible, since it may have been modified by a ScanHmmAugmenter.  The first test is important for the symbolic case
					// don't bother with transition stuff
					(*currTable)[state]=ClosedSemiring::GetImpossible();
				}
				else {

					if (!doLocal) {
						typename ScoreHelperType::ScoreType bestTransitionScoreEmitless=ClosedSemiring::GetImpossible();
						int bestTransitionChildNum;
						bestTransitionScoreEmitless=FindBestTransitionScoreEmitless<worryAboutStartStateInTransitions> (bestTransitionChildNum,hmm,tableForChildState,state,startState);
						// easy -- it's just the transition score, so we might as well do this now
						typename ScoreHelperType::ScoreType bestScoreChildAndTransition=ClosedSemiring::ExtensionOperator(emitScore,bestTransitionScoreEmitless);
						(*currTable)[state]=bestScoreChildAndTransition;
						// not worth adapting to Symbolic case: assert (ClosedSemiring::IsViterbi() || bestScoreChildAndTransition>=0); // forward alg scores must be positive since they're probabilities, although I allow one's >1
						if (bestTransitionChildNum!=-1) {
							scoreCollector.DecidedTransition(windowLast,state,bestTransitionChildNum);
						}
					}
					else {
						assert(doLocal);

						// Looking at scancyk.c, there's the following disjoint possibilities: (1) transition & emit, if any (the normal case), (2) [if this is a left HMM state] leftwardBeginsc + emit, if any, (3) [if this is a right HMM state] endsc + score at up-to-window-length skip [at corresponding left state] + emit, if any.

						int bestTransitionChildNum=-1;
						typename ScoreHelperType::ScoreType bestScore=ClosedSemiring::GetSummaryOperatorIdentity();

						assert(!(hmm.IsEmittingState(state) && !emissionType.IsEmitOkay()));

						bestScore=ClosedSemiring::ExtensionOperator(emitScore,FindBestTransitionScoreEmitless<worryAboutStartStateInTransitions> (bestTransitionChildNum,hmm,tableForChildState,state,startState));
#ifdef _DEBUG
						if (scoreDumpFile!=NULL) {
							fprintf(scoreDumpFile,"\tstate = %d, emits = %s, emit=%g\n",state,hmm.IsEmittingState(state)?"true":"false",ScoreHelperType::ScoreTypeToFloat(emitScore));
							typename Hmm::State childState=hmm.GetInvalidState();
							if (bestTransitionChildNum>=0) {
								childState=hmm.GetNthChildState(state,bestTransitionChildNum);
							}
							fprintf(scoreDumpFile,"\t\twinLast=%d.  best transition + emit = %g (transition to %d)\n",windowLast,ScoreHelperType::ScoreTypeToFloat(bestScore),childState);
						}
#endif

						// leftward beginsc
						bool pickedBeginscOverTransition=false;
						bool tryLeft=true;
						if (startOnlyAtFirstPos) {
							if (hmm.IsEmittingState(state)) {
								if (windowLast!=1) { // only do it when there's a start
									tryLeft=false;
								}
							}
							else {
								if (windowLast!=0) {
									tryLeft=false;
								}
							}
						}
						if (tryLeft) {
							typename ScoreHelperType::ScoreType leftwardBeginsc=ClosedSemiring::ExtensionOperator(emitScore,ClosedSemiring::GetLeftwardBegin(hmm,state));
							if (ClosedSemiring::SummaryOperatorWithMaxNotification(bestScore,leftwardBeginsc)) {  // >= to make sure it makes a decision even if no emission was allowed (so bestScore is not set) and beginsc is IMPOSSIBLE.
								pickedBeginscOverTransition=true;
								bestScore=leftwardBeginsc;
							}
#ifdef _DEBUG
							if (scoreDumpFile!=NULL) {
								fprintf(scoreDumpFile,"\t\twinLast=%d.  leftBeginsc + emit = %g\n",windowLast,ScoreHelperType::ScoreTypeToFloat(leftwardBeginsc));
							}
#endif
						}

						// endsc
						bool pickedEndsc=false;
						int endscLinkNumPicked=-1;
						int leftWindowLast=-1;
						for (int linkNum=0; linkNum<hmm.GetNumEndscLinksToLeft(state); linkNum++) {

							const typename Hmm::State leftLinkedState=hmm.GetEndscLinksToLeft_State(state,linkNum);
							assert(leftLinkedState<state); // that's what left means in our representation, and I'd also like to ensure that there aren't any self-loops here
							const typename ScoreHelperType::ScoreType endsc=ClosedSemiring::GetEnd(hmm,state,linkNum);
							const typename ScoreHelperType::ScoreType leftLinkedStateScore=scanHmmAugmenter.GetBestScoreInWindow(leftLinkedState,windowLast,leftWindowLast); // setting leftWindowLast is optional
							const typename ScoreHelperType::ScoreType thisEndscScore = ClosedSemiring::ExtensionOperator(emitScore,ClosedSemiring::ExtensionOperator( endsc , leftLinkedStateScore));

#ifdef _DEBUG
							if (scoreDumpFile!=NULL) {
								fprintf(scoreDumpFile,"\t\twinLast=%d.  local end + emit = %g  (emit=%g, leftLinkedState=%d, leftLinkedStateScore=%g, endsc=%g)\n",windowLast,ScoreHelperType::ScoreTypeToFloat(thisEndscScore),ScoreHelperType::ScoreTypeToFloat(emitScore),leftLinkedState,ScoreHelperType::ScoreTypeToFloat(leftLinkedStateScore),ScoreHelperType::ScoreTypeToFloat(endsc));
							}
#endif

							if (ClosedSemiring::SummaryOperatorWithMaxNotification(bestScore,thisEndscScore)) {
								pickedEndsc=true;
								endscLinkNumPicked=linkNum;
								bestScore=thisEndscScore;
								//printf("HMM endsc (%d,%d)\n",state,windowLast);
							}
						}

						if (pickedEndsc) {
							// notify ScoreCollector
							assert(endscLinkNumPicked>=0);
							scoreCollector.UsedEndsc(windowLast,state,endscLinkNumPicked,leftWindowLast);
						}
						else {
							if (pickedBeginscOverTransition) {
								// notify ScoreCollector
								scoreCollector.UsedLeftwardBeginsc(windowLast,state);
							}
							else {
								// conventional transition to child (with optional emission)
								// notify ScoreCollector
								assert(bestTransitionChildNum!=-1 || !ClosedSemiring::IsViterbi());
								scoreCollector.DecidedTransition(windowLast,state,bestTransitionChildNum);
							}
						}
						(*currTable)[state]=bestScore;
					}
				}
			}

			if (AddToScore(overrideDynProgScoreOperation)) {
				(*currTable)[state] += scoreFromAugmenter;
			}

			scanHmmAugmenter.NotifyOfScore(windowLast,state,(*currTable)[state]);

			//printf("%d,%d,%f\n",pos,state,(*currTable)[state]);
			if (cykscanStats.collectHmmFullDynProgTable) {
				cykscanStats.fullHmmDynProgTable[windowLast][state]=ScoreHelperType::ScoreTypeToFloat((*currTable)[state]);
			}
			scoreCollector.NotifyOfScoreInTable(windowLast,state,(*currTable)[state]);
		}

		if (!doLocal) {
			bestScoreAtPos=(*currTable)[endState];
			bestTopLevelState=endState;
		}
		else {
			for (state=startState; state<=endState; state++) {
				// re-notify scores, for benefit of local ends -- see large text in exception thrown in HmmType1::CheckThatLocalEndsHaveUniqueRightState (HmmType1.cpp)
				scanHmmAugmenter.ReNotifyOfScore_EndingWindowLast(windowLast,state,(*currTable)[state]);
			}

			// the last state is included as a possible exit, using its rightwardBeginsc.  The code in scancyk.c does NOT give state 0 any special treatment.  This state's transition scores are clobbered to be impossible, so its scores are impossible, so when LOCAL_BEGIN is set, state 0's score gets clobbered always.
			bestTopLevelState=endState;
			bestScoreAtPos=ClosedSemiring::GetSummaryOperatorIdentity();
			for (state=startState; state<=endState; state++) {
				const typename ScoreHelperType::ScoreType thisScore=ClosedSemiring::ExtensionOperator((*currTable)[state],ClosedSemiring::GetRightwardBegin(hmm,state));
				if (ClosedSemiring::SummaryOperatorWithMaxNotification(bestScoreAtPos,thisScore)) {
					bestScoreAtPos=thisScore;
					bestTopLevelState=state;
				}
			}
		}
		assert(doLocal || bestTopLevelState==endState); // if we're not doing local, should always end at the end.
		if (!ScoreHelperType::IsSymbolic()) {
			scoreCollector.GotScore(windowLast,bestTopLevelState,ScoreHelperType::ScoreTypeToFloat(bestScoreAtPos));
			assertr(ClosedSemiring::IsViterbi() || ScoreHelperType::ScoreTypeToFloat(bestScoreAtPos)>=0); // forward alg scores must be positive since they're probabilities, although I allow one's >1.  When this assert has failed in the past, the problem was that things underflowed in floats, which led to various problems.
		}

		typename ScoreHelperType::ScoreType overrideScore;
		if (scanHmmAugmenter.OverrideFinalWindowLastScore (overrideScore,windowLast,bestScoreAtPos)) {
			bestScoreAtPos=overrideScore;
		}

		if (scoreDumpFile!=NULL) {
			for(state=startState; state<=endState; state++) {
				fprintf(scoreDumpFile,"%d,%d,%g\n",windowLast,state,ScoreHelperType::ScoreTypeToFloat((*currTable)[state]));
			}
			fprintf(scoreDumpFile,"best: %d,%g\n",bestTopLevelState,ScoreHelperType::ScoreTypeToFloat(bestScoreAtPos));
		}

		return bestScoreAtPos;
	}

	template <class EmissionType>
		inline static typename ScoreHelperType::ScoreType ScanHmmOnePos(const Hmm& hmm,ScoreCollector& scoreCollector,ScanHmmAugmenter& scanHmmAugmenter,vector<typename ScoreHelperType::ScoreType> *prevTable,vector<typename ScoreHelperType::ScoreType> *currTable,const int pos,const EmissionType& emissionType,const typename Hmm::State startState,const typename Hmm::State endState,CykscanStats& cykscanStats,bool startOnlyAtFirstPos)
	{
		if (hmm.DoLocal()) {
			if (hmm.GetFirstState()==startState) {
				return ScanHmmOnePos_CommitDoLocal<EmissionType,false>(hmm,scoreCollector,scanHmmAugmenter,prevTable,currTable,pos,emissionType,startState,endState,cykscanStats,startOnlyAtFirstPos);
			}
			else {
				return ScanHmmOnePos_CommitDoLocal<EmissionType,true>(hmm,scoreCollector,scanHmmAugmenter,prevTable,currTable,pos,emissionType,startState,endState,cykscanStats,startOnlyAtFirstPos);
			}
		}
		else {
			if (hmm.GetFirstState()==startState) {
				return ScanHmmOnePos_CommitDoLocal<EmissionType,false>(hmm,scoreCollector,scanHmmAugmenter,prevTable,currTable,pos,emissionType,startState,endState,cykscanStats,startOnlyAtFirstPos);
			}
			else {
				return ScanHmmOnePos_CommitDoLocal<EmissionType,true>(hmm,scoreCollector,scanHmmAugmenter,prevTable,currTable,pos,emissionType,startState,endState,cykscanStats,startOnlyAtFirstPos);
			}
		}
	}

	inline static typename ScoreHelperType::ScoreType InitTable(vector<typename ScoreHelperType::ScoreType> *prevTable,const Hmm& hmm,const int startPosToScan,ScoreCollector& scoreCollector,ScanHmmAugmenter& scanHmmAugmenter,const typename Hmm::State startState,const typename Hmm::State endState,CykscanStats& cykscanStats,bool startOnlyAtFirstPos)
	{
		// this is how we should initialize the tables -- non-emitting states can do transitions
		EmissionType_NoEmission emissionType;
		return ScanHmmOnePos(hmm,scoreCollector,scanHmmAugmenter,NULL,prevTable,startPosToScan-1,emissionType,startState,endState,cykscanStats,startOnlyAtFirstPos);
	}

	inline static void ScanHmmWindow(HitList& hmmHitList,const int startPosToScan,const int endPosToScan,const Hmm& hmm,float float_minLodScoreForHit,CykscanStats& cykscanStats,const char *rnaSequence,const int windowLen,ScoreCollector& scoreCollector,ScanHmmAugmenter& scanHmmAugmenter,vector<typename ScoreHelperType::ScoreType>& table1,vector<typename ScoreHelperType::ScoreType>& table2,typename ScoreHelperType::ScoreType& get_totalSummaryScore,const typename Hmm::State startState,const typename Hmm::State endState,const EmissionHandler& emissionHandler,bool startOnlyAtFirstPos) 
	{
		vector<typename ScoreHelperType::ScoreType> *currTable,*prevTable;
		currTable=&table1;
		prevTable=&table2;

		// for Forward alg, useless for Viterbi
		typename ScoreHelperType::ScoreType totalSummaryScore=ClosedSemiring::GetSummaryOperatorIdentity();

		// initialize prevTable
		scanHmmAugmenter.StartWindowLast(startPosToScan);
		const typename ScoreHelperType::ScoreType initBestScore=InitTable(prevTable,hmm,startPosToScan,scoreCollector,scanHmmAugmenter,startState,endState,cykscanStats,startOnlyAtFirstPos); // initBestScore is really just for debugging
		ClosedSemiring::ProcessScoreAtInit(initBestScore,float_minLodScoreForHit,startPosToScan,endPosToScan,windowLen,cykscanStats);

		int pos;
		for (pos=startPosToScan; pos<endPosToScan; pos++) {
			//printf("pos=%d\n",pos);
			const int windowLast=pos+1; // for concept-compatibility with my code in ScancykBranchAndBound.cpp

			const int actualNucLetter=rnaSequence[pos];
			scoreCollector.DidEmission(windowLast,actualNucLetter);

			if (scoreDumpFile!=NULL) {
				fprintf(scoreDumpFile,"starting scan windowLast=%d, nuc=%c\n",windowLast,nucs[actualNucLetter]);
			}

			scanHmmAugmenter.StartWindowLast(windowLast);
			typename ScoreHelperType::ScoreType bestScoreAtPos;
			if (actualNucLetter>=0 && actualNucLetter<Alphabet_size) {
				// normal, literal nucleotide
				typename EmissionHandler::EmissionType_Simple emissionType(actualNucLetter,emissionHandler);
				bestScoreAtPos=ScanHmmOnePos(hmm,scoreCollector,scanHmmAugmenter,prevTable,currTable,pos,emissionType,startState,endState,cykscanStats,startOnlyAtFirstPos);
			}
			else {
				// degenerate nucleotide
				typename EmissionHandler::EmissionType_Degenerate emissionType(actualNucLetter,emissionHandler);
				bestScoreAtPos=ScanHmmOnePos(hmm,scoreCollector,scanHmmAugmenter,prevTable,currTable,pos,emissionType,startState,endState,cykscanStats,startOnlyAtFirstPos);
			}
			scanHmmAugmenter.EndWindowLast(windowLast);

			//printf("%d,%f\n",windowLast,bestScoreAtPos);

			ClosedSemiring::ProcessScoreAtPos(totalSummaryScore,bestScoreAtPos,float_minLodScoreForHit,startPosToScan,endPosToScan,windowLast,windowLen,hmmHitList,cykscanStats,scoreCollector);

			std::swap(currTable,prevTable);
		}

		get_totalSummaryScore=totalSummaryScore;
	}

	static const typename Hmm::State FindLowestChildOfStartState(const Hmm& hmm,const typename Hmm::State startState,const typename Hmm::State endState)
	{
		typename Hmm::State lowestChildOfStartState=startState; // assume there aren't any states before
		for (typename Hmm::State state=startState; state<=endState; state++) {
			for (int childNum=0; childNum<hmm.GetNumChildren(state); childNum++) {
				typename Hmm::State childState=hmm.GetNthChildState(state,childNum);
				lowestChildOfStartState=std::min(lowestChildOfStartState,childState);
			}
		}
		return lowestChildOfStartState;
	}

public:
	typedef ScanHmmAugmenter ScanHmmAugmenterType;

	inline static void ScanHmm (HitList& hmmHitList,const HitList& inputHitList,const Hmm& hmm,float float_minLodScoreForHit,CykscanStats& cykscanStats,const char *rnaSequence, int windowLen,ScoreCollector& scoreCollector,ScanHmmAugmenter& scanHmmAugmenter,typename ScoreHelperType::ScoreType& get_totalSummaryScore,const typename Hmm::State startState,const typename Hmm::State endState,const EmissionHandler& emissionHandler=EmissionHandler(),bool startOnlyAtFirstPos=false)
	{
		assertr(doLocal==hmm.DoLocal()); // we must be called with the right templating.  Strictly enforce this.
		if (doLocal) {
			assertr(ScanHmmAugmenter::Is_ScanHmmAugmenter_GenericLocal());
		}

		// okay, some cut&pasting here
#ifdef SCAN_DUMP
		if (scanDumpFile==NULL) {
			scanDumpFile=fopen("scanDumpFile.txt","wt");
		}
#endif
        // this strategy didn't work; I'm actually going to insert the code to check it explicitly:		const typename Hmm::State lowestChildOfStartState=FindLowestChildOfStartState(hmm,startState,endState);

		// only used for Forward alg
		typename ScoreHelperType::ScoreType totalSummaryScore=ClosedSemiring::GetSummaryOperatorIdentity();

		hmmHitList.clear();

		vector<typename ScoreHelperType::ScoreType> table1,table2,emitlessScoreTable;
		table1.resize(hmm.GetNumStates());
		table2.resize(hmm.GetNumStates());

		scanHmmAugmenter.StartSequence(hmm,rnaSequence);

		if (cykscanStats.collectHmmFullDynProgTable) {
			cykscanStats.fullHmmDynProgTable.resize(inputHitList.back().second+1);
		}

		HitList::const_iterator inputHitListIterator;
		for (inputHitListIterator=inputHitList.begin(); inputHitListIterator!=inputHitList.end(); inputHitListIterator++) {

			int startPosToScan=inputHitListIterator->first;
			int endPosToScan=inputHitListIterator->second;

			scanHmmAugmenter.StartInterval(hmm,startPosToScan,endPosToScan);
			scoreCollector.StartInterval(startPosToScan,endPosToScan);
			typename ScoreHelperType::ScoreType subtotalSummaryScore;
			ScanHmmWindow(hmmHitList,startPosToScan,endPosToScan,hmm,float_minLodScoreForHit,cykscanStats,rnaSequence,windowLen,scoreCollector,scanHmmAugmenter,table1,table2,subtotalSummaryScore,startState,endState,emissionHandler,startOnlyAtFirstPos);
			scanHmmAugmenter.EndInterval();

			if (! ClosedSemiring::IsViterbi()) {
				totalSummaryScore=ClosedSemiring::SummaryOperator(totalSummaryScore,subtotalSummaryScore);
			}
		}

		get_totalSummaryScore=totalSummaryScore;
	}
	// convenience -- default startState,endState
	static void ScanHmm (HitList& hmmHitList,const HitList& inputHitList,const Hmm& hmm,float float_minLodScoreForHit,CykscanStats& cykscanStats,const char *rnaSequence, int windowLen,ScoreCollector& scoreCollector,ScanHmmAugmenter& scanHmmAugmenter,typename ScoreHelperType::ScoreType& get_totalSummaryScore,const EmissionHandler& emissionHandler=EmissionHandler(),bool startOnlyAtFirstPos=false)
	{
		const typename Hmm::State startState=hmm.GetFirstState();
		const typename Hmm::State endState=hmm.GetActualLastState();
		ScanHmm(hmmHitList,inputHitList,hmm,float_minLodScoreForHit,cykscanStats,rnaSequence,windowLen,scoreCollector,scanHmmAugmenter,get_totalSummaryScore,startState,endState,emissionHandler,startOnlyAtFirstPos);
	}

	static void ScanHmm (HitList& hmmHitList,const HitList& inputHitList,const InfernalHmm& hmmInput,float float_minLodScoreForHit,CykscanStats& cykscanStats,const char *rnaSequence, int windowLen,ScoreCollector& scoreCollector,ScanHmmAugmenter& scanHmmAugmenter,typename ScoreHelperType::ScoreType& get_totalSummaryScore,const InfernalHmm::State infernalStartState,const InfernalHmm::State infernalEndState,bool enableHmmTranslationCaching=true,const EmissionHandler& emissionHandler=EmissionHandler(),bool startOnlyAtFirstPos=false)
	{
		const typename Hmm::State startState=InfernalHmm::StateToInt(infernalStartState);
		const typename Hmm::State endState=InfernalHmm::StateToInt(infernalEndState);

		Hmm *hmmUse;
		if (enableHmmTranslationCaching) {
			typename CmHmmToFastHmm::iterator findHmmIter=cmHmmToFastHmm.find(hmmInput.GetUniqueId());
			if (findHmmIter==cmHmmToFastHmm.end()) {
				// not built yet
				hmmUse=new Hmm;
				hmmUse->Init(hmmInput);
				const typename CmHmmToFastHmm::value_type insertee(hmmInput.GetUniqueId(),hmmUse);
				cmHmmToFastHmm.insert(insertee);
			}
			else {
				hmmUse=findHmmIter->second;
			}
		}
		else {
			hmmUse=new Hmm;
			hmmUse->Init(hmmInput);
		}

#if 0
		Stopwatch_t     *watch;
		watch = StopwatchCreate();
		StopwatchZero(watch);
		StopwatchStart(watch);
#endif

		ScanHmm (hmmHitList,inputHitList,*hmmUse,float_minLodScoreForHit,cykscanStats,rnaSequence,windowLen,scoreCollector,scanHmmAugmenter,get_totalSummaryScore,startState,endState,emissionHandler,startOnlyAtFirstPos);

		if (!enableHmmTranslationCaching) {
			delete hmmUse;
		}

#if 0
		StopwatchStop(watch);
		StopwatchDisplay(stdout, "\nScanHMM CPU time: ", watch);
#endif
	}

	// convenience -- default startState,endState
	static void ScanHmm (HitList& hmmHitList,const HitList& inputHitList,const InfernalHmm& hmmInput,float float_minLodScoreForHit,CykscanStats& cykscanStats,const char *rnaSequence, int windowLen,ScoreCollector& scoreCollector,ScanHmmAugmenter& scanHmmAugmenter,typename ScoreHelperType::ScoreType& get_totalSummaryScore,bool enableHmmTranslationCaching=true,const EmissionHandler& emissionHandler=EmissionHandler(),bool startOnlyAtFirstPos=false)
	{
		const typename InfernalHmm::State startState=hmmInput.GetFirstState();
		const typename InfernalHmm::State endState=hmmInput.GetActualLastState();
		ScanHmm(hmmHitList,inputHitList,hmmInput,float_minLodScoreForHit,cykscanStats,rnaSequence,windowLen,scoreCollector,scanHmmAugmenter,get_totalSummaryScore,startState,endState,enableHmmTranslationCaching,emissionHandler,startOnlyAtFirstPos);
	}
};
template <class ScoreHelper,class Hmm,class ScoreCollector,class ClosedSemiring,class ScanHmmAugmenter,bool doLocal,class EmissionHandler>
typename ClassScanHmm<ScoreHelper,Hmm,ScoreCollector,ClosedSemiring,ScanHmmAugmenter,doLocal,EmissionHandler>::CmHmmToFastHmm ClassScanHmm<ScoreHelper,Hmm,ScoreCollector,ClosedSemiring,ScanHmmAugmenter,doLocal,EmissionHandler>::cmHmmToFastHmm;



// forward alg for AlignmentConsensusAndScoring.cpp
class NoUnderflowDoubleScoreHelper {
public:
	typedef NoUnderflowDouble ScoreType;

	inline static ScoreType GetZeroScore (void) {
		return 0;
	}
	inline static ScoreType GetImpossibleScore (void) {
		return IMPOSSIBLE;
	}
	inline static ScoreType GetNegativeInfinityScore (void) {
		return -DBL_MAX;
	}
	inline static ScoreType FloatToScoreType (float t) {
		return t;
	}
	inline static float ScoreTypeToFloat (ScoreType t) {
		return (float)(t.ToDouble_ZeroOnUnderflow());
	}
	inline static bool IsSymbolic (void) {
		return false;
	}
};
template <class ScoreCollector>
NoUnderflowDouble ForwardOfLiteralString_internal (const HmmType1& hmm,const char *rnaSequence,int rnaSequenceLen,HmmType1::State state1,int windowLast1,HmmType1::State state2,int windowLast2,ScoreCollector& scoreCollector)
{
	CykscanStats cykscanStats;
	NoUnderflowDouble totalSummaryScore;
	const int windowLen=rnaSequenceLen;
	HitList hmmHitList,inputHitList;
	inputHitList.Init(windowLen);
	float minLodScoreForHit=0;

	typedef Forward_ClosedSemiring_NoUnderflowDouble ClosedSemiring;
	EmissionHandler_RestrictPositions<NoUnderflowDoubleScoreHelper,HmmType1,ClosedSemiring> emissionHandler;
	emissionHandler.SetStateMustEmitHere(state1,windowLast1,state2,windowLast2);

	if (hmm.DoLocal()) {
		assertr(false); // I've been assuming we wouldn't use local HMMs for this, although it should work -- in theory
	}
	else {
		bool allowOnlyFullString=true;
		ScanHmmAugmenter_NULL<NoUnderflowDoubleScoreHelper,HmmType1,ScoreCollector,ClosedSemiring> nullAugmenter;
		ClassScanHmm<NoUnderflowDoubleScoreHelper,HmmType1,ScoreCollector,ClosedSemiring,ScanHmmAugmenter_NULL<NoUnderflowDoubleScoreHelper,HmmType1,ScoreCollector,ClosedSemiring>,false,EmissionHandler_RestrictPositions<NoUnderflowDoubleScoreHelper,HmmType1,ClosedSemiring> >::ScanHmm(hmmHitList,inputHitList,hmm,minLodScoreForHit,cykscanStats,rnaSequence,windowLen,scoreCollector,nullAugmenter,totalSummaryScore,emissionHandler,allowOnlyFullString);
	}

	return totalSummaryScore;
}
NoUnderflowDouble ForwardOfLiteralString (const HmmType1& hmm,const char *rnaSequence,int rnaSequenceLen,HmmType1::State state1,int windowLast1,HmmType1::State state2,int windowLast2)
{
	ScoreCollector_NULL scoreCollector;
	return ForwardOfLiteralString_internal (hmm,rnaSequence,rnaSequenceLen,state1,windowLast1,state2,windowLast2,scoreCollector);
}
NoUnderflowDouble ForwardOfLiteralString (const HmmType1& hmm,const char *rnaSequence,int rnaSequenceLen,HmmType1::State state1,int windowLast1,HmmType1::State state2,int windowLast2,DynProgTableScoreCollector_OwnMemory<NoUnderflowDouble>& dynProgTable)
{
	return ForwardOfLiteralString_internal (hmm,rnaSequence,rnaSequenceLen,state1,windowLast1,state2,windowLast2,dynProgTable);
}
