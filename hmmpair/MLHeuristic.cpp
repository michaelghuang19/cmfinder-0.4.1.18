#include "hmmpair.h"

void SetProbToHmmTransition(InfernalHmm& infernalHmm,InfernalHmm::State hmmState,InfernalHmm::State toHmmState,double setProb)
{
	//printf("SetProbToHmmTransition.  %d -> %d  =  %lf\n",InfernalHmm::StateToInt(hmmState),InfernalHmm::StateToInt(toHmmState),addProb);
	if (toHmmState<hmmState) {
		std::swap(toHmmState,hmmState);
	}
	int child=infernalHmm.GetChildNumEitherWay_Slow(hmmState,toHmmState);
	infernalHmm.SetTransitionProbDirectly(hmmState,child,(float)(setProb));
}
void AddProbToHmmTransition(InfernalHmm& infernalHmm,InfernalHmm::State hmmState,InfernalHmm::State toHmmState,double addProb)
{
	//printf("AddProbToHmmTransition.  %d -> %d  +=  %lf\n",InfernalHmm::StateToInt(hmmState),InfernalHmm::StateToInt(toHmmState),addProb);
	if (toHmmState<hmmState) {
		std::swap(toHmmState,hmmState);
	}
	int child=infernalHmm.GetChildNumEitherWay_Slow(hmmState,toHmmState);
	float p=infernalHmm.GetTransitionProbDirectly(hmmState,child);
	infernalHmm.SetTransitionProbDirectly(hmmState,child,(float)(p+addProb));
}
void AddProbToHmmEmission(InfernalHmm& infernalHmm,InfernalHmm::State hmmState,int nuc,double addProb)
{
	float p=infernalHmm.GetEmissionProbDirectly(hmmState,nuc);
	infernalHmm.SetEmissionProbDirectly(hmmState,nuc,(float)(p+addProb));
}
/* After visiting an IL or IR state from a split set state, go to the next states after the IL/IR state, and update the probabilities of what happens on the other side of the HMM.  e.g. if you go from MP to IL to IR, then on the right side of the HMM, there's an MP to IR transition.
*/
void FullJumpViaInsert(const CovarianceModel& cm,InfernalHmm& infernalHmm,CovarianceModel::State cmToState,double cmTransitionProb,double selfLoopInfiniteSeries,InfernalHmm::State otherSideHmmState,bool isLeft)
{
	for (int child=0; child<cm.GetNumChildren(cmToState); child++) {
		CovarianceModel::State cmNextState=cm.GetNthChildState(cmToState,child);
		double thisTransitionProb=cm.GetNthChildTransitionProb(cmToState,child);
		//if (!cm.IsInsertState(cmNextState)) {
		if (cmNextState!=cmToState) {
			InfernalHmm::State otherSideHmmNextState=isLeft ? infernalHmm.GetHmmLeftStateOfCmState(cmNextState) : infernalHmm.GetHmmRightStateOfCmState(cmNextState);
			assert(otherSideHmmNextState!=InfernalHmm::GetInvalidState());
			AddProbToHmmTransition(infernalHmm,otherSideHmmState,otherSideHmmNextState,cmTransitionProb*selfLoopInfiniteSeries*thisTransitionProb);
		}
	}
}
double GetSelfLoopInfiniteSeries(const CovarianceModel& cm,CovarianceModel::State cmInsertState)
{
	assert(cm.IsInsertState(cmInsertState));
	int selfLoopChild=0;
	assert(cm.GetNthChildState(cmInsertState,selfLoopChild)==cmInsertState); // I'm expecting it to be self-loop
	double selfLoopProb=cm.GetNthChildTransitionProb(cmInsertState,selfLoopChild);
	assert(selfLoopProb>=0 && selfLoopProb<1);
	double selfLoopInfiniteSeries=1/(1-selfLoopProb);
	return selfLoopInfiniteSeries;
}
void FindExpectedEntriesIntoState (VectorByCmState<double>& numExpectedEntriesIntoState,const CovarianceModel& cm,bool applyInsertPropogationBug)
{
	numExpectedEntriesIntoState.assign(cm.GetNumStates(),0.0);

	if (cm.DoLocal()) {
		for (CovarianceModel::State state=cm.GetFirstState(); state!=cm.GetLastState(); state++) {
			numExpectedEntriesIntoState[state]=cm.GetBeginProb(state);
		}
	}
	else {
		numExpectedEntriesIntoState[cm.GetFirstState()]=1.0;
	}

	for (CovarianceModel::State state=cm.GetFirstState(); state!=cm.GetLastState(); state++) {

		if (cm.IsBifurcation(state)) {
			// for BIF, both children _must_ be visited
			numExpectedEntriesIntoState[cm.GetLeftBifurcationChild(state)] += numExpectedEntriesIntoState[state];
			numExpectedEntriesIntoState[cm.GetRightBifurcationChild(state)] += numExpectedEntriesIntoState[state];
		}
		else {
			// we don't have to explicitly deal with local ends, since no more states get visited

			// first, detect any self loops & adjust our expected entries.  Note: our children should also get this boost.  e.g., suppose self-loop prob is 0.5 and transition prob is 0.5.  We're clearly eventually going to go thru the transition, so we should multiply it by 1/(1-selfLoop)
			double thisStateExpectedEntries=numExpectedEntriesIntoState[state];
			double selfLoopMultiple=1.0;
			for (int child=0; child<cm.GetNumChildren(state); child++) {
				CovarianceModel::State toState=cm.GetNthChildState(state,child);
				double transitionProb=cm.GetNthChildTransitionProb(state,child);
				if (state==toState) { // only self-loops for now
					assert(selfLoopMultiple==1.0); // shouldn't have >1 self-loops
					if (transitionProb<0 || transitionProb>=1) {
						throw SimpleStringException("Input CM has insert state self-loop with prob not in range 0..1.  %s:%d",__FILE__,__LINE__);
					}
					selfLoopMultiple=1.0/(1.0-transitionProb);
				}
			}
			numExpectedEntriesIntoState[state] *= selfLoopMultiple;
			if (!applyInsertPropogationBug) {
				thisStateExpectedEntries=numExpectedEntriesIntoState[state]; // this should propogate to children, but not doing this has improved results in one test, so I want to investigate more
			}

			for (int child=0; child<cm.GetNumChildren(state); child++) {
				CovarianceModel::State toState=cm.GetNthChildState(state,child);
				double transitionProb=cm.GetNthChildTransitionProb(state,child);
				if (state!=toState) { // ignore self-loops
					numExpectedEntriesIntoState[toState] += transitionProb * thisStateExpectedEntries;
				}
			}
		}
	}
}
void BuildHmm_MaxLikeliPathCorrespondence_EmitsOnly(InfernalHmm& infernalHmm,const CovarianceModel& cm,const VectorByCmState<double>& numExpectedEntriesIntoState)
{
	infernalHmm.ZeroAllEmitProbs();

	CovarianceModel::State cmState;
	for (cmState=cm.GetFirstState(); cmState!=cm.GetLastState(); cmState++) {
		if (cm.DoLocal() && cm.GetNode(cmState)==cm.GetFirstNode()) {
			// avoid this
		}
		else {
			double expectedEnteringIntoState=numExpectedEntriesIntoState[cmState]; // adjust by the probability of being in the original split set state
			switch (cm.GetNumSymbolsEmitted(cmState)) {
				case 0:
					// nothing to do
					break;
				case 1:
					// simply copy them to the appropriate state (left/right)
					{
						bool isLeftEmitting=cm.EmitsLeft(cmState);
						InfernalHmm::State hmmState=isLeftEmitting ? infernalHmm.GetHmmLeftStateOfCmState(cmState) : infernalHmm.GetHmmRightStateOfCmState(cmState);
						assert(hmmState!=InfernalHmm::GetInvalidState());
						for (int nuc=0; nuc<Alphabet_size; nuc++) {
							AddProbToHmmEmission(infernalHmm,hmmState,nuc,expectedEnteringIntoState*cm.GetSingletEmissionProb(cmState,nuc));
						}
					}
					break;
				case 2:
					// project the pair emission probabilities onto a left & right
					{
						assert(Alphabet_size==4);
						double leftProbs[4]={0,0,0,0},rightProbs[4]={0,0,0,0},total=0;
						for (int leftNuc=0; leftNuc<Alphabet_size; leftNuc++) {
							for (int rightNuc=0; rightNuc<Alphabet_size; rightNuc++) {
								double pairProb=pow2(cm.GetPairEmissionScore(cmState,leftNuc,rightNuc));
								total += pairProb;
								leftProbs[leftNuc] += pairProb;
								rightProbs[rightNuc] += pairProb;
							}
						}
						assertr(fabs(total-1.0)<1e-4); // should sum to 1

						// now simply set them in
						{
							InfernalHmm::State leftHmmState=infernalHmm.GetHmmLeftStateOfCmState(cmState);
							for (int leftNuc=0; leftNuc<Alphabet_size; leftNuc++) {
								AddProbToHmmEmission(infernalHmm,leftHmmState,leftNuc,expectedEnteringIntoState*leftProbs[leftNuc]);
							}
						}

						{
							InfernalHmm::State rightHmmState=infernalHmm.GetHmmRightStateOfCmState(cmState);
							for (int rightNuc=0; rightNuc<Alphabet_size; rightNuc++) {
								AddProbToHmmEmission(infernalHmm,rightHmmState,rightNuc,expectedEnteringIntoState*rightProbs[rightNuc]);
							}
						}
					}
					break;
			}
		}
	}

	infernalHmm.NormalizeEmissionsToStrictProbabilitiesViaProbabilities();
}
/* NEW STRATEGY: for emissions, same as before; seemed to work.
For transitions, for each node, consider all sub-paths from split set state to split set state.  For each sub-path, compute it's probability, & add the total probability to each transition in the corresponding HMM path.
Note: this code resembles the code in Cm2Hmm.cpp that enumerates all paths to figure out the inequalities, but I think it's sufficiently, that it'd require a lot of refactoring and consequently much risk to a very important part of my code.
*/
typedef std::list<CovarianceModel::State> CmStateList;
typedef std::list<InfernalHmm::State> HmmStateList;
struct CorrespondingCmAndHmmPath {
	double cmPathProb;
	CmStateList cmPath;
	HmmStateList leftHmmPath,rightHmmPath;
};
typedef std::list<CorrespondingCmAndHmmPath> CorrespondingCmAndHmmPathList;
void EnumerateCorrespondingCmAndHmmSubPaths_Recurse (CorrespondingCmAndHmmPathList& correspondingCmAndHmmPathList,const CovarianceModel& cm,const InfernalHmm& infernalHmm,const CovarianceModel::Node cmNode,const CorrespondingCmAndHmmPath& subPathSoFar)
{
	const CovarianceModel::State cmState=subPathSoFar.cmPath.back();

	if (cm.GetNode(cmState)!=cmNode) {
		// got to next node - add sub-path & terminate this line of search 
		correspondingCmAndHmmPathList.push_back(subPathSoFar);
		return;
	}

	for (int child=0; child<cm.GetNumChildren(cmState); child++) {

		const CovarianceModel::State cmToState=cm.GetNthChildState(cmState,child);
		if (cmState!=cmToState) { // avoid self-loops

			CorrespondingCmAndHmmPath nextSubPathSoFar=subPathSoFar;

			double transitionProb=cm.GetNthChildTransitionProb(cmState,child);
			nextSubPathSoFar.cmPathProb *= transitionProb;

			nextSubPathSoFar.cmPath.push_back(cmToState);
			InfernalHmm::State leftHmmState=infernalHmm.GetHmmLeftStateOfCmState(cmToState);
			InfernalHmm::State rightHmmState=infernalHmm.GetHmmRightStateOfCmState(cmToState);
			if (leftHmmState!=InfernalHmm::GetInvalidState()) {
				assert(nextSubPathSoFar.leftHmmPath.back()!=leftHmmState);
				nextSubPathSoFar.leftHmmPath.push_back(leftHmmState);
			}
			if (rightHmmState!=InfernalHmm::GetInvalidState()) {
				assert(nextSubPathSoFar.rightHmmPath.back()!=rightHmmState);
				nextSubPathSoFar.rightHmmPath.push_back(rightHmmState);
			}

			EnumerateCorrespondingCmAndHmmSubPaths_Recurse (correspondingCmAndHmmPathList,cm,infernalHmm,cmNode,nextSubPathSoFar);
		}
	}
}
/*
enumerate all CM sub-paths (and the corresponding HMM sub-path) starting at a split set state of one node and ending at the split set state of the next node, except
(1) no emissions
(2) never explore self-loops
(3) BIF nodes have no sub-paths
*/
void EnumerateCorrespondingCmAndHmmSubPaths (CorrespondingCmAndHmmPathList& correspondingCmAndHmmPathList,const CovarianceModel& cm,const InfernalHmm& infernalHmm,CovarianceModel::Node cmNode)
{
	if (cm.IsBifurcationNode(cmNode)) {
		// no sub-paths, really
		return;
	}

	for (CovarianceModel::State cmState=cm.GetFirstStateOfNode(cmNode); cmState!=cm.GetLastSplitSetStateOfNode(cmNode); cmState++) {
		if (cm.DoLocal() && cm.GetNode(cmState)==cm.GetFirstNode()) {
		}
		else {
			CorrespondingCmAndHmmPath subPathSoFar;
			subPathSoFar.cmPathProb=1;
			subPathSoFar.cmPath.push_back(cmState);
			InfernalHmm::State leftHmmState=infernalHmm.GetHmmLeftStateOfCmState(cmState);
			InfernalHmm::State rightHmmState=infernalHmm.GetHmmRightStateOfCmState(cmState);
			assert(leftHmmState!=InfernalHmm::GetInvalidState() && rightHmmState!=InfernalHmm::GetInvalidState());
			subPathSoFar.leftHmmPath.push_back(leftHmmState);
			subPathSoFar.rightHmmPath.push_back(rightHmmState);

			EnumerateCorrespondingCmAndHmmSubPaths_Recurse(correspondingCmAndHmmPathList,cm,infernalHmm,cmNode,subPathSoFar);
		}
	}
}
void AdjustSubPathProbsForStartStateAndSelfLoopProbs(CorrespondingCmAndHmmPathList& correspondingCmAndHmmPathList,const CovarianceModel& cm,const VectorByCmState<double>& numExpectedEntriesIntoState)
{
	for (CorrespondingCmAndHmmPathList::iterator i=correspondingCmAndHmmPathList.begin(); i!=correspondingCmAndHmmPathList.end(); i++) {
		CorrespondingCmAndHmmPath& subPath=*i;

		// adjust prob by start state
		CovarianceModel::State cmStateState=subPath.cmPath.front();
		subPath.cmPathProb *=numExpectedEntriesIntoState[cmStateState];

		// now adjust it be infinite series of insertion states
		// zlh: changed my mind again -- it's a hassle, but it seems right // zlh: NO, I changed my mind.  All paths with the insert state will have the same adjustment, so all transitions out of the insert state will have the same multiple, so it's just a multiple factor that will be divided out when we normalize the probabilites.
		for (CmStateList::iterator i=subPath.cmPath.begin(); i!=subPath.cmPath.end(); i++) {
			CovarianceModel::State cmState=*i;
			if (cm.IsInsertState(cmState)) {
				int cmSelfLoopChildNum=cm.GetChildNum_Slow(cmState,cmState);
				double selfLoopProb=cm.GetNthChildTransitionProb(cmState,cmSelfLoopChildNum);
				subPath.cmPathProb *= 1.0/(1.0-selfLoopProb);
			}
		}
	}
}
void AddTransitionProbsForHmmPath(InfernalHmm& infernalHmm,const HmmStateList& hmmPath,double addProb)
{
	InfernalHmm::State hmmPrevState=InfernalHmm::GetInvalidState();
	for (HmmStateList::const_iterator i=hmmPath.begin(); i!=hmmPath.end(); i++) {
		InfernalHmm::State hmmState=*i;
		if (hmmPrevState!=InfernalHmm::GetInvalidState()) {
			AddProbToHmmTransition(infernalHmm,hmmPrevState,hmmState,addProb);
		}
		hmmPrevState=hmmState;
	}
}
void MakeMLHeuristicFromStrictProbabilisticCM(InfernalHmm& infernalHmm,const CovarianceModel& cm,const ExtraCm2HmmInfo& extraCm2HmmInfo,bool applyInsertPropogationBug)
{
	VectorByCmState<double> numExpectedEntriesIntoState;
	FindExpectedEntriesIntoState(numExpectedEntriesIntoState,cm,applyInsertPropogationBug);

	// first do emissions, which are easier
	BuildHmm_MaxLikeliPathCorrespondence_EmitsOnly(infernalHmm,cm,numExpectedEntriesIntoState);

	// clear all transition probs
	infernalHmm.ZeroAllTransitionProbs();

	for (CovarianceModel::Node cmNode=cm.GetFirstNode(); cmNode!=cm.GetLastNode(); cmNode++) {

		// find paths
		CorrespondingCmAndHmmPathList correspondingCmAndHmmPathList;
		EnumerateCorrespondingCmAndHmmSubPaths (correspondingCmAndHmmPathList,cm,infernalHmm,cmNode);
		AdjustSubPathProbsForStartStateAndSelfLoopProbs(correspondingCmAndHmmPathList,cm,numExpectedEntriesIntoState);

		// now we can just go thru each side of the HMM
		for (CorrespondingCmAndHmmPathList::const_iterator i=correspondingCmAndHmmPathList.begin(); i!=correspondingCmAndHmmPathList.end(); i++) {
			const CorrespondingCmAndHmmPath& correspondingCmAndHmmPath=*i;
			AddTransitionProbsForHmmPath(infernalHmm,correspondingCmAndHmmPath.leftHmmPath,correspondingCmAndHmmPath.cmPathProb);
			AddTransitionProbsForHmmPath(infernalHmm,correspondingCmAndHmmPath.rightHmmPath,correspondingCmAndHmmPath.cmPathProb);
		}

	}

	// can handle local stuff separately
	for (CovarianceModel::State cmState=cm.GetFirstState(); cmState!=cm.GetLastState(); cmState++) {
		InfernalHmm::State hmmLeftState=infernalHmm.GetHmmLeftStateOfCmState(cmState);
		InfernalHmm::State hmmRightState=infernalHmm.GetHmmRightStateOfCmState(cmState);
		if (cm.GetBeginProb(cmState)!=0) {
			float halfProb=(float)(cm.GetBeginProb(cmState));
			assertr(hmmLeftState!=InfernalHmm::GetInvalidState() && hmmRightState!=InfernalHmm::GetInvalidState());
			infernalHmm.SetLeftwardBeginProbDirectly(hmmLeftState,halfProb);
			infernalHmm.SetRightwardBeginProbDirectly(hmmRightState,halfProb);
		}
		if (cm.GetEndProb(cmState)!=0 && cm.GetStateType(cmState)!=E_st) {
			assertr(!cm.IsInsertState(cmState)); // formula doesn't work with insert states (and their nefarious self-loops), and I didn't think anyone used them anyway
			bool found=false;
			for (int link=0; link<infernalHmm.GetNumEndscLinksToLeft(hmmRightState); link++) {
				if (hmmLeftState==infernalHmm.GetEndscLinkToLeft_State(hmmRightState,link)) {
					assertr(!found);
					float p=(float)(cm.GetEndProb(cmState) * numExpectedEntriesIntoState[cmState]);
					infernalHmm.SetEndProbDirectly(hmmRightState,link,p);
					found=true;
				}
			}
			assertr(found);
		}
	}

	// and add in prob for insert state self-loops
	for (CovarianceModel::State cmState=cm.GetFirstState(); cmState!=cm.GetLastState(); cmState++) {
		if (cm.DoLocal() && cm.GetNode(cmState)==cm.GetFirstNode()) {
		}
		else {
			if (cm.IsInsertState(cmState)) {
				int cmSelfLoopChildNum=cm.GetChildNum_Slow(cmState,cmState);
				double selfLoopProb=cm.GetNthChildTransitionProb(cmState,cmSelfLoopChildNum);
				double expectedInsertStateEntries=numExpectedEntriesIntoState[cmState];

				InfernalHmm::State hmmState=cm.EmitsLeft(cmState) ? infernalHmm.GetHmmLeftStateOfCmState(cmState) : infernalHmm.GetHmmRightStateOfCmState(cmState);
				//int hmmSelfLoopChildNum=infernalHmm.GetChildNum_Slow(hmmState,hmmState);

				SetProbToHmmTransition(infernalHmm,hmmState,hmmState,expectedInsertStateEntries*selfLoopProb);
			}
		}
	}

	infernalHmm.NormalizeTransitionsToStrictProbabilitiesViaProbabilities();
}
void FindEquivalentProbCM(CovarianceModel& probCM,const CovarianceModel& cm)
{
	probCM.CopyFrom(cm); // get structure
	VectorByCmState<NoUnderflowDouble_NoAutoCast> pathSumAtState,emitTotalAtState;
	pathSumAtState.resize(cm.GetNumStates());
	emitTotalAtState.resize(cm.GetNumStates());
	for (CovarianceModel::State state=cm.GetActualLastState(); state>=cm.GetFirstState(); state--) {

		if (state==CovarianceModel::IntToState(88)) {
			//int q=9;
		}

		if (cm.GetStateType(state)==B_st || cm.GetStateType(state)==E_st) {
			switch (cm.GetStateType(state)) {
				case B_st:
					pathSumAtState[state]=pathSumAtState[cm.GetLeftBifurcationChild(state)] * pathSumAtState[cm.GetRightBifurcationChild(state)];
					emitTotalAtState[state]=1;
					break;
				case E_st:
					assertr(cm.GetNumChildren(state)==0);
					pathSumAtState[state]=1;
					emitTotalAtState[state]=1;
					break;
			}
		}
		else {
			NoUnderflowDouble_NoAutoCast emit=0.0;
			if (cm.IsEmittingState(state)) {
				if (cm.EmitsLeftAndRight(state)) {
					for (int left=0; left<Alphabet_size; left++) {
						for (int right=0; right<Alphabet_size; right++) {
							emit += cm.GetPairEmissionProbDirectly(state,left,right);
						}
					}

					// normalize these emits, now that we've used the info
					for (int left=0; left<Alphabet_size; left++) {
						for (int right=0; right<Alphabet_size; right++) {
							float p=ToProb(cm.GetPairEmissionProbDirectly(state,left,right),emit);
							probCM.SetPairEmissionProbDirectly(state,left,right,p);
						}
					}
				}
				else {
					for (int nuc=0; nuc<Alphabet_size; nuc++) {
						emit += cm.GetEmissionProbDirectly(state,nuc);
					}

					// normalize these emits, now that we've used the info
					for (int nuc=0; nuc<Alphabet_size; nuc++) {
						float p=ToProb(cm.GetEmissionProbDirectly(state,nuc),emit);
						probCM.SetEmissionProbDirectly(state,nuc,p);
					}
				}
			}
			else {
				emit=1;
			}
			emitTotalAtState[state]=emit;

			NoUnderflowDouble_NoAutoCast transitionProb=0;
			NoUnderflowDouble_NoAutoCast selfLoopMultiplier=1;
			if (cm.DoLocal()) {
				assertr(cm.GetStateType(state)!=E_st);
				transitionProb += cm.GetEndProbDirectly(state);
			}
			if (state==CovarianceModel::IntToState(298)) {
				//int q=9;
			}
			for (int child=0; child<cm.GetNumChildren(state); child++) {

				CovarianceModel::State childState=cm.GetNthChildState(state,child);

				if (childState==state) {
					NoUnderflowDouble_NoAutoCast selfLoopProb=emit*(double)(cm.GetTransitionProbDirectly(state,child));
					if (selfLoopProb>=1) {
						throw SimpleStringException("Trying to probabilize a CM that has self-loop prob>1.  State=#%d, self loop prob=%lg",CovarianceModel::StateToInt(state),selfLoopProb.ToDouble_ZeroOnUnderflow());
					}
					selfLoopMultiplier=1.0/(1.0-selfLoopProb);
				}
				else {
					transitionProb += pathSumAtState[childState] * (double)(cm.GetTransitionProbDirectly(state,child));
				}
			}

			pathSumAtState[state]=selfLoopMultiplier*emit*transitionProb;
		}
	}

	// now adjust transitions & local ends
	for (CovarianceModel::State state=cm.GetFirstState(); state!=cm.GetLastState(); state++) {
		if (cm.GetStateType(state)==B_st || cm.GetStateType(state)==E_st) {
			switch (cm.GetStateType(state)) {
				case B_st:
					// you're excused
					assertr(cm.GetEndProb(state)==0); // this has been true, so we don't need to adjust anything
					break;
				case E_st:
					// you're also excused
					assertr(cm.GetNumChildren(state)==0);
					break;
			}
		}
		else {
			NoUnderflowDouble_NoAutoCast totalProb=0;
			NoUnderflowDouble_NoAutoCast selfLoopProb=0;
			if (cm.DoLocal()) {
				totalProb += cm.GetEndProbDirectly(state);
			}
			for (int child=0; child<cm.GetNumChildren(state); child++) {
				CovarianceModel::State childState=cm.GetNthChildState(state,child);
				if (childState==state) {
					assertr(selfLoopProb==0);
					selfLoopProb=emitTotalAtState[state]*(double)(cm.GetTransitionProbDirectly(state,child));
					probCM.SetTransitionProbDirectly(state,child,selfLoopProb.ToFloat());
				}
				else {
					totalProb += cm.GetTransitionProbDirectly(state,child)*pathSumAtState[childState];
				}
			}
			if (selfLoopProb>=1.0) {
				throw SimpleStringException("While probabilizing a CM: a self-loop path has probability >1, so the sum of products of paths will be infinite, so this algorithm cannot find an equivalent probabilistic model.");
			}
			totalProb /= (1.0-selfLoopProb); // adjust it s.t. everything else fits into 1-selfLoopProb

			if (totalProb==0) {
				// don't change anything; these should remain 0
			}
			else {
				if (cm.DoLocal()) {
					float endProb=cm.GetEndProbDirectly(state);
					float p=ToProb(endProb,totalProb);
					probCM.SetEndProbDirectly(state,p);
				}
				for (int child=0; child<cm.GetNumChildren(state); child++) {
					CovarianceModel::State childState=cm.GetNthChildState(state,child);
					if (childState==state) {
						// do nothing
					}
					else {
						float p=ToProb(pathSumAtState[childState]*cm.GetTransitionProbDirectly(state,child),totalProb);
						probCM.SetTransitionProbDirectly(state,child,p);
					}
				}
			}
		}
	}

	// and local begins
	if (cm.DoLocal()) {
		NoUnderflowDouble_NoAutoCast totalProb=0;
		for (CovarianceModel::State state=cm.GetFirstState(); state!=cm.GetLastState(); state++) {
			totalProb += cm.GetBeginProbDirectly(state)*pathSumAtState[state];
		}
		for (CovarianceModel::State state=cm.GetFirstState(); state!=cm.GetLastState(); state++) {
			float p=ToProb(pathSumAtState[state]*cm.GetBeginProbDirectly(state),totalProb);
			probCM.SetBeginProbDirectly(state,p);
		}
	}

	probCM.SetLODScores();
}
void BuildHmm_MaxLikeliPathCorrespondence(const CovarianceModel& cm_,bool doLocalAlignment,InfernalHmm& createdInfernalHmm,const InfernalHmm& startInfernalHmm,const ExtraCm2HmmInfo& extraCm2HmmInfo,bool applyInsertPropogationBug)
{
	CovarianceModel cm;
	cm.CopyFrom(cm_);
	InfernalHmm infernalHmm;
	infernalHmm.CopyFrom(startInfernalHmm);

	cm.MultiplyEmitsBy(0.25);
	cm.HackInsertScoresToStrictProbs();

	CovarianceModel probCM;
	FindEquivalentProbCM(probCM,cm); // ::CMRenormalize doesn't work with --local
	//probCM.DumpCsv("data/rsearch/prob-cm.csv");

	MakeMLHeuristicFromStrictProbabilisticCM(infernalHmm,probCM,extraCm2HmmInfo,applyInsertPropogationBug);

	InfernalHmm outputHmm;
	outputHmm.CopyFromWithEscHack(infernalHmm);
	//outputHmm.CopyFrom(infernalHmm);
	createdInfernalHmm.CopyFrom(outputHmm);
}
void BuildHmm_MaxLikeliPathCorrespondence(char *cmFileName,bool doLocalAlignment,InfernalHmm& createdInfernalHmm,const std::string& programParams_,Cm2Hmm_HmmBuildType hmmType,bool applyInsertPropogationBug)
{
	std::string programParams(programParams_);

	CovarianceModel cm;
	InfernalHmm infernalHmm;
	ExtraCm2HmmInfo extraCm2HmmInfo;
	MakeCmHmmForOptimization (cm,infernalHmm,extraCm2HmmInfo,cmFileName,doLocalAlignment,programParams,hmmType,false);

	BuildHmm_MaxLikeliPathCorrespondence(cm,doLocalAlignment,createdInfernalHmm,infernalHmm,extraCm2HmmInfo,applyInsertPropogationBug);
	createdInfernalHmm.AddBuildDescription(programParams);
}
void BuildHmm_MaxLikeliPathCorrespondence(char *cmFileName,bool doLocalAlignment,const char *hmmBinFileName,const std::string& programParams,Cm2Hmm_HmmBuildType hmmType,bool applyInsertPropogationBug)
{
	InfernalHmm outputHmm;
	BuildHmm_MaxLikeliPathCorrespondence(cmFileName,doLocalAlignment,outputHmm,programParams,hmmType,applyInsertPropogationBug);
	outputHmm.SaveInBinary(hmmBinFileName);
}
