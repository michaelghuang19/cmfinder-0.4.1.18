#include "hmmpair.h"

////////////////////
// BackwardHmmType1

BackwardHmmType1::BackwardHmmType1 ()
{
}
BackwardHmmType1::~BackwardHmmType1 ()
{
}
void BackwardHmmType1::Init (const HmmType1& hmm)
{
	HmmType1::Init(hmm); // setup some things

	for (State state=GetFirstState(); state<GetLastState()/2; state++) {
		State otherState=ReverseState(state);
		std::swap(stateInfoVector[state],stateInfoVector[otherState]);
		if (DoLocal()) {
			std::swap(localStateInfoVector[state],localStateInfoVector[otherState]);
		}
		for (int nuc=0; nuc<Alphabet_size; nuc++) {
			std::swap(singletEmissionScores[nuc][state],singletEmissionScores[nuc][otherState]);
		}
		std::swap(forwardAlgInfoVector[state],forwardAlgInfoVector[otherState]);
	}

	// CM/HMM mapping is weird, so remove it for now (It could work, I just don't feel like thinking it thru)
	hmm2CmStateVector.clear();
	cm2HmmStateVector.clear();

	// transitions
	// first clear them
	for (State state=GetFirstState(); state<GetLastState(); state++) {
		stateInfoVector[state].numChildren=0;
	}
	// now cycle thru
	for (State sourceState=GetActualLastState(); sourceState>=GetFirstState(); sourceState--) {
		State thisState=ReverseState(sourceState);
		for (int child=0; child<hmm.GetNumChildren(sourceState); child++) {
			State sourceChildState=hmm.GetNthChildState(sourceState,child);
			State thisChildState=ReverseState(sourceChildState);
			if (stateInfoVector[thisChildState].numChildren==0) {
				stateInfoVector[thisChildState].firstChild=thisState;
				stateInfoVector[thisChildState].tsc[0]=hmm.GetNthChildTsc(sourceState,child);
				forwardAlgInfoVector[thisChildState].transitionProbs[0]=hmm.GetNthChildTransitionProb(sourceState,child);
			}
			else {
				assertr(stateInfoVector[thisChildState].firstChild+stateInfoVector[thisChildState].numChildren==thisState); // or we did it in the wrong order
				stateInfoVector[thisChildState].tsc[stateInfoVector[thisChildState].numChildren]=hmm.GetNthChildTsc(sourceState,child);
				forwardAlgInfoVector[thisChildState].transitionProbs[stateInfoVector[thisChildState].numChildren]=hmm.GetNthChildTransitionProb(sourceState,child);
			}
			stateInfoVector[thisChildState].numChildren++;
		}
	}

	if (DoLocal()) {
		// swap left/right begins, and clear local ends for now
		for (State state=GetFirstState(); state!=GetLastState(); state++) {
			std::swap(localStateInfoVector[state].leftwardBeginsc,localStateInfoVector[state].rightwardBeginsc);
			std::swap(forwardAlgInfoVector[state].leftBeginProb,forwardAlgInfoVector[state].rightBeginProb);
			localStateInfoVector[state].endscLinksToLeft.resize(0);
		}
		// set up local ends
		for (State sourceState=GetFirstState(); sourceState!=GetLastState(); sourceState++) {
			State thisState=ReverseState(sourceState);
			for (int link=0; link<hmm.GetNumEndscLinksToLeft(sourceState); link++) {
				State sourceLeftState=hmm.GetEndscLinksToLeft_State(sourceState,link);
				State thisLeftState=ReverseState(sourceLeftState);

				assertr(thisLeftState>thisState); // left/rightness is reversed
				int newLink=(int)(localStateInfoVector[thisLeftState].endscLinksToLeft.size());
				localStateInfoVector[thisLeftState].endscLinksToLeft.resize(newLink+1);
				localStateInfoVector[thisLeftState].endscLinksToLeft[newLink].hmmLeftState=thisState;
				localStateInfoVector[thisLeftState].endscLinksToLeft[newLink].endsc=hmm.GetEndscLinksToLeft_Endsc(sourceState,link);
				forwardAlgInfoVector[thisLeftState].ends.resize(newLink+1);
				forwardAlgInfoVector[thisLeftState].ends[newLink].leftLink=thisState;
				forwardAlgInfoVector[thisLeftState].ends[newLink].endProb=hmm.GetEndProb(sourceState,link);
			}
		}
	}

	CheckThatLocalEndsHaveUniqueRightState();
}


///////////////////////
// HmmType1

void HmmType1::CheckThatLocalEndsHaveUniqueRightState (void)
{
	if (!DoLocal()) {
		return; // non-local HMMs are fine
	}

	State state;
	for (state=GetFirstState(); state!=GetLastState(); state++) {
		localStateInfoVector[state].isLeftStateOfLocalEnd=false;
	}

	for (state=GetFirstState(); state!=GetLastState(); state++) {
		for (int link=0; link<GetNumEndscLinksToLeft(state); link++) {
			State leftState=GetEndscLinksToLeft_State(state,link);
			if (localStateInfoVector[leftState].isLeftStateOfLocalEnd) {
				throw SimpleStringException("Ooops.  I assume that local HMMs have the following property: consider all local ends.  For any left state of a local end, there is only 1 right state, i.e. there are never two local ends that share the same left state.  (A stronger property, that they never share a right state either, is convenient, but only necessary for things like --ml-sample-baum-welch.)  Why is this important?  In the code in ScanHMM.cpp, if the right state of a local end does NOT emit, then matching \\epsilon to the local end corresponds to going to the left state in the _same_ windowLast.  However, if the right state DOES emit, matching \\epsilon means going to the left state in windowLast-1 -- this is the same issue in transitions where we need a different column in the dyn prog table depending on whether the curr state emits or not.  If the left state has multiple right states (via local ends), and some of these right states emit & some don't, then we need to be clever about whether or not to return the current number, depending on which right state we're coming from.  If there's a unique right state, it's easy.  So, I just assume there's a unique right state.  (A less demanding assumption would be just that either all right states emit, or none of them do.)");
			}
			localStateInfoVector[leftState].isLeftStateOfLocalEnd=true;
			localStateInfoVector[leftState].rightStateOfLocalEndEmits=IsEmittingState(state);
		}
	}
}
void HmmType1::Init (const HmmType1& hmm)
{
	numStates=hmm.numStates;
	stateInfoVector=hmm.stateInfoVector;
	singletEmissionScores=hmm.singletEmissionScores;
	doLocal=hmm.doLocal;
	localStateInfoVector=hmm.localStateInfoVector;
	hmm2CmStateVector=hmm.hmm2CmStateVector;
	cm2HmmStateVector=hmm.cm2HmmStateVector;
	forwardAlgInfoVector=hmm.forwardAlgInfoVector;
	localEndSelfLoopScore=hmm.localEndSelfLoopScore;
	localEndSelfLoopProb=hmm.localEndSelfLoopProb;
}
void HmmType1::Init (const InfernalHmm& unreversedSourceHmm)
{
	InfernalHmm sourceHmm;
	sourceHmm.CopyReverseOf(unreversedSourceHmm);

	numStates=sourceHmm.GetNumStates();

	stateInfoVector.resize(numStates);
	singletEmissionScores.resize(Alphabet_size);
	for (int nuc=0; nuc<Alphabet_size; nuc++) {
		singletEmissionScores[nuc].resize(numStates);
	}

	localEndSelfLoopScore=unreversedSourceHmm.GetLocalEndSelfLoopScore();
	localEndSelfLoopProb=unreversedSourceHmm.GetLocalEndSelfLoopProb();

	doLocal=false;
	if (sourceHmm.DoLocal()) {
		doLocal=true;
		localStateInfoVector.resize(numStates);
	}

	hmm2CmStateVector=sourceHmm.GetCmStateVector();
	cm2HmmStateVector=sourceHmm.GetHmmStateVector();

	CovarianceModel::State state;
	for (state=sourceHmm.GetFirstState(); state!=sourceHmm.GetLastState(); state++) {
		int i_state=CovarianceModel::StateToInt(state);

		int stateType=sourceHmm.GetStateType(state);
		assert(stateType==ML_st || stateType==IL_st || stateType==E_st || stateType==D_st || stateType==PASSTHRU_st || stateType==S_st); // these are the only types that HMM should have
		if (stateType==ML_st || stateType==IL_st) {
			stateInfoVector[i_state].isEmitting=true;
		}
		else {
			if (stateType==E_st || stateType==D_st || stateType==S_st || stateType==PASSTHRU_st) {
				stateInfoVector[i_state].isEmitting=false;
			}
			else {
				throw SimpleStringException("Internal error: Unexpected state type in HMM represented as a CovarianceModel");
			}
		}

		for (int nuc=0; nuc<Alphabet_size; nuc++) {
			float floatScore;
			if (stateInfoVector[i_state].isEmitting) {
				floatScore=sourceHmm.GetSingletEmissionScore(state,nuc);
			}
			else {
				floatScore=0;
			}
			singletEmissionScores[nuc][i_state]=ScoreHelperType::FloatToScoreType(floatScore);
		}

		stateInfoVector[i_state].numChildren=sourceHmm.GetNumChildren(state);
		assert(stateInfoVector[i_state].numChildren<=MAX_CHILDREN);
		stateInfoVector[i_state].firstChild=-1;
		stateInfoVector[i_state].isRightState=sourceHmm.IsRightState(state);

		if (stateInfoVector[i_state].numChildren>0) {
			stateInfoVector[i_state].firstChild=CovarianceModel::StateToInt(sourceHmm.GetNthChildState(state,0));

			for (int childNum=0; childNum<sourceHmm.GetNumChildren(state); childNum++) {
				if (GetNthChildState(i_state,childNum)!=CovarianceModel::StateToInt(sourceHmm.GetNthChildState(state,childNum))) {
					throw SimpleStringException("Internal error: CovarianceModel apparently doesn't have consecutively-numbered states any more.  I hate it when things change on me, so I'm going to die now.");
				}
				stateInfoVector[i_state].tsc[childNum]=ScoreHelperType::FloatToScoreType(sourceHmm.GetNthChildTsc(state,childNum));
			}
		}

		if (sourceHmm.DoLocal()) {

			localStateInfoVector[i_state].leftwardBeginsc=sourceHmm.GetLeftwardBeginsc(state);
			localStateInfoVector[i_state].rightwardBeginsc=sourceHmm.GetRightwardBeginsc(state);

			int numLinks=sourceHmm.GetNumEndscLinksToLeft(state);
			localStateInfoVector[i_state].endscLinksToLeft.resize(numLinks);
			for (int link=0; link<numLinks; link++) {
				localStateInfoVector[i_state].endscLinksToLeft[link].hmmLeftState=InfernalHmm::StateToInt(sourceHmm.GetEndscLinkToLeft_State(state,link));
				localStateInfoVector[i_state].endscLinksToLeft[link].endsc=sourceHmm.GetEndscLinkToLeft_Endsc(state,link);
			}
		}
	}

	// Forward Alg info
	forwardAlgInfoVector.resize(GetNumStates());
	for (State state=GetFirstState(); state!=GetLastState(); state++) {
		for (int i=0; i<GetNumChildren(state); i++) {
			forwardAlgInfoVector[state].transitionProbs[i]=(float)(pow2(GetNthChildTsc(state,i)));
		}
		assert(Alphabet_size==MAXABET); // I used MAXABET to define the bounds of singletEmissionProbs in ScanHMM.h
		for (int nuc=0; nuc<Alphabet_size; nuc++) {
			if (IsEmittingState(state)) {
				forwardAlgInfoVector[state].singletEmissionProbs[nuc]=(float)(pow2(GetSingletEmissionScore(state,nuc)));
			}
		}

		if (DoLocal()) {
			forwardAlgInfoVector[state].leftBeginProb=(float)(pow2(GetLeftwardBeginsc(state)));
			forwardAlgInfoVector[state].rightBeginProb=(float)(pow2(GetRightwardBeginsc(state)));
			forwardAlgInfoVector[state].ends.resize(GetNumEndscLinksToLeft(state));
			for (int link=0; link<GetNumEndscLinksToLeft(state); link++) {
				forwardAlgInfoVector[state].ends[link].leftLink=GetEndscLinksToLeft_State(state,link);
				forwardAlgInfoVector[state].ends[link].endProb=(float)(pow2(GetEndscLinksToLeft_Endsc(state,link)));
			}
		}
	}

	CheckThatLocalEndsHaveUniqueRightState();
}
void HmmType1::Dump (FILE *out) const
{
	fprintf(out,"# states: %d\n",GetNumStates());
	for (State state=GetFirstState(); state!=GetLastState(); state++) {
		fprintf(out,"state #%d : %s\n",state,IsEmittingState(state) ? "EMITTING" : "silent");
		if (IsEmittingState(state)) {
			fprintf(out,"\t");
			for (int nuc=0; nuc<MAXABET; nuc++) {
				fprintf(out,"%c: %f , ",nucs[nuc],GetSingletEmissionScore(state,nuc));
			}
			fprintf(out,"\n");
		}
		fprintf(out,"\t");
		for (int child=0; child<GetNumChildren(state); child++) {
			fprintf(out,"#%d: %f , ",GetNthChildState(state,child),GetNthChildTsc(state,child));
		}
		fprintf(out,"\n");
	}
}
void HmmType1::Dump (const char *fileName) const
{
	FILE *out=ThrowingFopen(fileName,"wt");
	Dump(out);
	fclose(out);
}



////////////////////
// TransitionCounter

TransitionCounter::TransitionCounter (const HmmType1& _hmm,int windowLen_)
: hmm(_hmm)
{
	windowLen=windowLen_;

	transitionsInParseHistogram.assign(1000,0);

	transitionCounts.resize(hmm.GetNumStates());
	emitCounts.resize(hmm.GetNumStates());
	if (hmm.DoLocal()) {
		leftBeginCounts.assign(hmm.GetNumStates(),0);
		rightBeginCounts.assign(hmm.GetNumStates(),0);
		endCounts.resize(hmm.GetNumStates());
	}

	HmmType1::State state;
	for (state=hmm.GetFirstState(); state!=hmm.GetLastState(); state++) {
		transitionCounts[state].assign(hmm.GetNumChildren(state),0);
		if (hmm.IsEmittingState(state)) {
			emitCounts[state].assign(Alphabet_size,0);
		}
		if (hmm.DoLocal()) {
			endCounts[state].assign(hmm.GetNumEndscLinksToLeft(state),0);
		}
	}

	numSamples=0;
	totalViterbiScore=0;
}
TransitionCounter::~TransitionCounter ()
{
}
int TransitionCounter::GetWindowLen (void) const
{
	return windowLen;
}
double TransitionCounter::GetAverageViterbiScore (void) const
{
	return totalViterbiScore/(double)(numSamples);
}
void TransitionCounter::AddNumTransitionsInParse(int numTransitionsInParse)
{
	int oldSize=(int)(transitionsInParseHistogram.size());
	if (numTransitionsInParse>=oldSize) {
		transitionsInParseHistogram.insert(transitionsInParseHistogram.end(),numTransitionsInParse*2-oldSize,0);
	}
	transitionsInParseHistogram[numTransitionsInParse]++;
}
void TransitionCounter::AddSample (int windowLast,float viterbiScore)
{
	numSamples++;
	totalViterbiScore += viterbiScore;
}
int TransitionCounter::GetNumSamples (void) const
{
	return numSamples;
}
int TransitionCounter::GetTransitionCount_Unreversed (int unreversed_fromState,int unreversed_toState) const
{
	int fromState=unreversed_toState;
	int toState=unreversed_fromState;
	return transitionCounts[fromState][hmm.GetChildNum_Slow(fromState,toState)];
}
int TransitionCounter::GetTransitionCount_Unreversed (CovarianceModel::State unreversed_fromState,CovarianceModel::State unreversed_toState) const
{
	return GetTransitionCount_Unreversed(CovarianceModel::StateToInt(unreversed_fromState),CovarianceModel::StateToInt(unreversed_toState));
}
int TransitionCounter::GetEmitCount (int state,int nuc) const
{
	return emitCounts[state][nuc];
}
int TransitionCounter::GetEmitCount (CovarianceModel::State state,int nuc) const
{
	return GetEmitCount(CovarianceModel::StateToInt(state),nuc);
}
double TransitionCounter::GetEntryProbability_Unreversed (int unreversed_toState,int unreversed_alternateToStateFirst,int unreversed_alternateToStateLast) const
{
	assert(unreversed_alternateToStateLast>=unreversed_alternateToStateFirst);
	assert(unreversed_toState>=unreversed_alternateToStateFirst && unreversed_toState<unreversed_alternateToStateLast);

	int fromState=unreversed_toState;
	int alternateFromStateFirst=unreversed_alternateToStateFirst;
	int alternateFromStateLast=unreversed_alternateToStateLast;

	int total=0;
	int numerator=-1;
	for (int currFromState=alternateFromStateFirst; currFromState<alternateFromStateLast; currFromState++) {

		int totalForState=0;
		for (int childNum=0; childNum<hmm.GetNumChildren(currFromState); childNum++) {
			if (currFromState!=hmm.GetNthChildState(currFromState,childNum)) {
				totalForState += transitionCounts[currFromState][childNum];
			}
		}
		if (currFromState==fromState) {
			numerator=totalForState;
		}

		total += totalForState;
	}
	assert(numerator>=0);

	if (total==0) {
		// special case for start node
		assert(numerator==0);
		return 0; // path is never taken at all, so prob is 0
	}
	else {
		double result=(double)(numerator)/(double)(total);
		assert(result>=0 && result<=1);
		return result;
	}
} 
double TransitionCounter::GetEntryProbability_Unreversed (CovarianceModel::State unreversed_toState,CovarianceModel::State unreversed_alternateToStateFirst,CovarianceModel::State unreversed_alternateToStateLast) const
{
	return GetEntryProbability_Unreversed(CovarianceModel::StateToInt(unreversed_toState),CovarianceModel::StateToInt(unreversed_alternateToStateFirst),CovarianceModel::StateToInt(unreversed_alternateToStateLast));
}
double TransitionCounter::GetTransitionFrequency_Unreversed (int unreversed_fromState,int unreversed_toState,double pseudocount) const
{
	int fromState=unreversed_toState;
	int toState=unreversed_fromState;

	if (fromState==toState) {
		assert(false); // we shouldn't need to explore self-loops
		return 1;
	}

	double prChild=-1; // we'll actually set this later
	double prState=0;
	for (size_t i=0; i<transitionCounts[fromState].size(); i++) {
		HmmType1::State childState=hmm.GetNthChildState(fromState,(int)i);
		if (childState!=fromState) { // avoid self-loops
			prState += (double)(transitionCounts[fromState][i]) + pseudocount;
			if (childState==toState) {
				assert(prChild==-1); // should only set this once
				prChild=(double)(transitionCounts[fromState][i]) + pseudocount;
			}
		}
	}
	assert(prChild!=-1); // should set this

	if (prState==0) {
		assert(prChild==0);
		return 0;
	}
	else {
		double result=prChild/prState;
		assert(result>=0 && result<=1);
		return result;
	}
}
double TransitionCounter::GetTransitionFrequency_Unreversed (CovarianceModel::State unreversed_fromState,CovarianceModel::State unreversed_toState,double pseudocount) const
{
	return GetTransitionFrequency_Unreversed (CovarianceModel::StateToInt(unreversed_fromState),CovarianceModel::StateToInt(unreversed_toState),pseudocount);
}
double TransitionCounter::GetEmitFrequency (int state,int nuc,double pseudocount) const
{
	double prNuc=(double)(emitCounts[state][nuc]) + pseudocount;

	double prState=0;
	for (size_t i=0; i<emitCounts[state].size(); i++) {
		prState += (double)(emitCounts[state][i]) + pseudocount;
	}

	if (prState==0) {
		assert(prNuc==0);
		return 0;
	}
	else {
		double result=prNuc/prState;
		assert(result>=0 && result<=1);
		return result;
	}
}
double TransitionCounter::GetEmitFrequency (CovarianceModel::State state,int nuc,double pseudocount) const
{
	return GetEmitFrequency(CovarianceModel::StateToInt(state),nuc,pseudocount);
}
double TransitionCounter::GetNumTransitions (int state,int child) const
{
	return transitionCounts[state][child];
}
double TransitionCounter::GetNumEmissions (int state,int nuc) const
{
	return emitCounts[state][nuc];
}
double TransitionCounter::GetNumLeftBegins (int state) const
{
	return leftBeginCounts[state];
}
double TransitionCounter::GetNumRightBegins (int state) const
{
	return rightBeginCounts[state];
}
double TransitionCounter::GetNumEnds (int state,int link) const
{
	return endCounts[state][link];
}
void TransitionCounter::AddTransition (int state,int childNum)
{
	transitionCounts[state][childNum]++;
}
void TransitionCounter::AddEmit (int state,int windowLast,int nuc)
{
	if (nuc<0 || nuc>=Alphabet_size) {
		throw SimpleStringException("Got degenerate nucleotide while collecting Viterbi traceback statistics; this is not cool.\n");
	}
	emitCounts[state][nuc]++;
}
void TransitionCounter::AddLeftBegin (int state)
{
	if (hmm.DoLocal()) {
		leftBeginCounts[state]++;
	}
	else {
		assertr(state==hmm.GetFirstState());
	}
}
void TransitionCounter::AddRightBegin (int state)
{
	if (hmm.DoLocal()) {
		rightBeginCounts[state]++;
	}
	else {
		assertr(state==hmm.GetActualLastState());
	}
}
void TransitionCounter::AddEnd (int state,int link,int nucsSkipped)
{
	assertr(hmm.DoLocal());
	endCounts[state][link]++;
}
double TransitionCounter::GetNumParsesWithNumTransitions (int numTransitions) const
{
	if (numTransitions>=(int)(transitionsInParseHistogram.size())) {
		return 0;
	}
	return transitionsInParseHistogram[numTransitions];
}
void TransitionCounter::DumpParseLenHistogram (FILE *out)
{
	fprintf(out,"Parse-length histogram (Viterbi):\n");
	{
		size_t highest=transitionsInParseHistogram.size();
		while (highest>1 && transitionsInParseHistogram[highest-1]==0) {
			highest--;
		}
		for (size_t i=0; i<highest; i++) {
			fprintf(out,"\t%u transitions: %d \t(%lg)\n",i,transitionsInParseHistogram[i],(double)(transitionsInParseHistogram[i])/(double)(numSamples));
		}
	}
}
bool TransitionCounter::IgnoreWindowLastBelowWindowLen (void) const
{
	return true;
}
void TransitionCounter::DoneVisitingState(void)
{
}
void TransitionCounter::VisitingState(int state)
{
}
void TransitionCounter::Dump (FILE *out)
{
	fprintf(out,"Total samples: %d\n",numSamples);

	DumpParseLenHistogram(out);

	HmmType1::State state,childState;
	for (state=hmm.GetFirstState(); state!=hmm.GetLastState(); state++) {
		fprintf(out,"state #%d:\n",state);
		fprintf(out,"\tTransitions:\n");
		for (int childNum=0; childNum<hmm.GetNumChildren(state); childNum++) {
			childState=hmm.GetNthChildState(state,childNum);
			fprintf(out,"\t\t--> state #%d: %d\n",childState,transitionCounts[state][childNum]);
		}
		if (hmm.IsEmittingState(state)) {
			fprintf(out,"\tEmissions:\n");
			for (int nuc=0; nuc<Alphabet_size; nuc++) {
				fprintf(out,"\t\temit %c: %d\n",nucs[nuc],emitCounts[state][nuc]);
			}
		}

		if (hmm.DoLocal()) {
			fprintf(out,"\tLocal begins:\n");
			if (hmm.IsLeftwardLocalBeginValid(state)) {
				fprintf(out,"\t\tleftward : %d\n",leftBeginCounts[state]);
			}
			if (hmm.IsRightwardLocalBeginValid(state)) {
				fprintf(out,"\t\trightward: %d\n",rightBeginCounts[state]);
			}
			fprintf(out,"\tLocal ends:\n");
			for (int link=0; link<hmm.GetNumEndscLinksToLeft(state); link++) {
				fprintf(out,"\t\tto state %d: %d\n",hmm.GetEndscLinksToLeft_State(state,link),endCounts[state][link]);
			}
		}
	}
}
void TransitionCounter::StartInterval (int firstWindowLast,int lastWindowLast)
{
}
void TransitionCounter::DoneSample (void)
{
}
