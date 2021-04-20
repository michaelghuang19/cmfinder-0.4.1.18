#include "hmmpair.h"

namespace lpsolve {
	typedef int lprec; // dummy definition
}

// note: what I call "normal states" and "(post) insert states", Sean Eddy calls "split set" and "insert set" in Eddy, S.R. (2002), BMC Bioinformatics, 3:18 (see Table 3 and text on that page).  I guess I should ideally have read that paper before doing this...

// my LP problems were unbounded, so just in case the problem recurs (this was from my very early work, but the lower bound seems pretty conservative)
#define HACK_LOWERBOUND -10000.0

//#define ENABLE_CACHING // building the HMM can take kind of a while.  Caching is somewhat deprecated (it seems better to deal with loading/saving HMMs in specific files), and I don't expect it to make it into a production version of this code.

#ifdef CM2HMM_ONLY
// this caching thing doesn't make sense for distributed code -- it's just a complication to explain
#undef ENABLE_CACHING
#endif

bool enableHmmCaching=true;

//#define DEBUG_DUMP // enable dumping of HMM build work, into hmmdump.txt
FILE *dumpFile=NULL;

//#define DUMP_INTERMEDIATE_HMMS // dump intermediate HMMs (structure only) as it builds each piece, & splices them

// THE FOLLOWING EQUIVOCATION IS NO LONGER NECESSARY - THIS VARIABLE WILL ALWAYS BE SET TO 'TRUE'
// I'm not sure if this is necessary/good.  currently I think it is
// if true, then every CM node has at least 1 left & at least 1 right node in the HMM
// e.g. a MATL node will have a dummy PASSTHRU node on the right
// if false, MATL nodes wouldn't have anything on the right
// advantage of true: it's well-defined how to move thru the HMM based on the CM, whereas
// without the passthru nodes, it's weird.  Also, without the passthru node, we get more info
// in the HMM, but this isn't reflected in the CM, so it doesn't help with an upper bound
// of the CM (e.g. a CM has ROOT,MATL,MATL,MATL,MATL,MATR.  if everyCmNodeHasLeftAndRightHmmNodes==false,
// we can hook up the ROOT directly to MATR & make better decisions on what state to enter in the MATR,
// but the CM doesn't do this anyway).
#define everyCmNodeHasLeftAndRightHmmNodes true

// the HMM consists of Left-emitting stuff from the CM, and Right-emitting stuff from the CM
enum HMM_Direction {
	HMM_Left,HMM_Right
};

struct CmNodeTypeData {
	// total # of states in the Covariance Model
	int numCmStates;
	// all nodes start with something normal, like a S_st, MP, ML, MR, D, that has to do with their function.  Then they can have an IL and/or IR state.  the normal states of a node always connect to the insert states of that node and the normal states of the next node.  The insert states of a node connect only to the normal states of the next node.
	int numNormalStates,numInsertPostStates;
	// state types, for verifying that they're the way I expect
	const int *stateTypes;
	
	// starts/ends segment
	bool isStartingNode,isEndingNode;

	// stuff about HMM encoding
	// # of HMM states to the left, # to the right
	int numHmmLeftStates,numHmmRightStates;
	// how many of those states are IL/IR
	int numHmmLeftInsertStates,numHmmRightInsertStates; // these are either 0 or 1
};
static int BEGR_stateTypes[]={S_st,IL_st,INVALID_st},
	BEGL_stateTypes[]={S_st,INVALID_st},
	ROOT_stateTypes[]={S_st,IL_st,IR_st,INVALID_st},
	BIF_stateTypes[]={B_st,INVALID_st},
	END_stateTypes[]={E_st,INVALID_st},
	MATP_stateTypes[]={MP_st,ML_st,MR_st,D_st,IL_st,IR_st,INVALID_st},
	MATL_stateTypes[]={ML_st,D_st,IL_st,INVALID_st},
	MATR_stateTypes[]={MR_st,D_st,IR_st,INVALID_st};
static CmNodeTypeData BEGR_data,BEGL_data,ROOT_data,BIF_data,END_data,MATP_data,MATL_data,MATR_data;
static bool isInit_CmNodeTypeData=false;
void InferStuffFromStateTypes(CmNodeTypeData& data,const int *stateTypes)
{
	data.stateTypes=stateTypes; // that was easy :-)

	data.isStartingNode=stateTypes[0]==S_st;
	data.isEndingNode=stateTypes[0]==E_st || stateTypes[0]==EL_st || stateTypes[0]==B_st;

	data.numCmStates=0;
	data.numNormalStates=0;
	data.numInsertPostStates=0;
	data.numHmmLeftStates=0;
	data.numHmmRightStates=0;
	data.numHmmLeftInsertStates=0;
	data.numHmmRightInsertStates=0;
	for (int i=0; stateTypes[i]!=INVALID_st; i++) {
		bool isInsertState=stateTypes[i]==IL_st || stateTypes[i]==IR_st;
		bool isStartingOrEndingState=stateTypes[i]==S_st || stateTypes[i]==E_st || stateTypes[i]==EL_st || stateTypes[i]==B_st;
		bool isLeft=stateTypes[i]==IL_st || stateTypes[i]==ML_st || stateTypes[i]==MP_st || stateTypes[i]==D_st || isStartingOrEndingState;
		bool isRight=stateTypes[i]==IR_st || stateTypes[i]==MR_st || stateTypes[i]==MP_st || stateTypes[i]==D_st || isStartingOrEndingState;

		data.numCmStates++;
		if (isInsertState) {
			data.numInsertPostStates++;
		}
		else {
			data.numNormalStates++;
		}
		if (stateTypes[i]==MP_st) {
			// don't put anything extra for HMM for MP_st - ML_st and MR_st should put it in
		}
		else {
			if (isLeft) {
				data.numHmmLeftStates++;
				if (isInsertState) {
					data.numHmmLeftInsertStates++;
				}
			}
			if (isRight) {
				data.numHmmRightStates++;
				if (isInsertState) {
					data.numHmmRightInsertStates++;
				}
			}
		}
	}

	if (stateTypes[0]==ML_st) {
		data.numHmmRightStates--; // generic alg overcounts because it thinks the D_st is on both sides
	}
	if (stateTypes[0]==MR_st) {
		data.numHmmLeftStates--; // generic alg overcounts because it thinks the D_st is on both sides
	}

	if (everyCmNodeHasLeftAndRightHmmNodes) {
		if (data.numHmmLeftStates==0) {
			data.numHmmLeftStates++;
		}
		if (data.numHmmRightStates==0) {
			data.numHmmRightStates++;
		}
	}
}
void InitCmNodeTypeData(void)
{
	InferStuffFromStateTypes(BEGR_data,BEGR_stateTypes);
	InferStuffFromStateTypes(BEGL_data,BEGL_stateTypes);
	InferStuffFromStateTypes(ROOT_data,ROOT_stateTypes);
	InferStuffFromStateTypes(BIF_data,BIF_stateTypes);
	InferStuffFromStateTypes(END_data,END_stateTypes);
	InferStuffFromStateTypes(MATP_data,MATP_stateTypes);
	InferStuffFromStateTypes(MATL_data,MATL_stateTypes);
	InferStuffFromStateTypes(MATR_data,MATR_stateTypes);

	isInit_CmNodeTypeData=true;
}
CmNodeTypeData GetCmNodeTypeData (int nodeType)
{
	if (!isInit_CmNodeTypeData) {
		InitCmNodeTypeData();
	}

	switch (nodeType) {
case BIF_nd:
	return BIF_data;
case END_nd:
	return END_data;
case MATP_nd:
	return MATP_data;
case MATL_nd:
	return MATL_data;
case MATR_nd:
	return MATR_data;
case BEGL_nd:
	return BEGL_data;
case BEGR_nd:
	return BEGR_data;
case ROOT_nd:
	return ROOT_data;
case DUMMY_nd:
	assert(false); // I thought DUMMY nodes get removed once the CM is made
	break;
default:
	assert(false);
	throw SimpleStringException("unknown node type when converting CM to HMM\n");
	break;
	}

	// keep compiler happy
	CmNodeTypeData data;
	data.isStartingNode=false;
	return data;
}

// make sure we got everything about CovarianceModels right
// don't use asserts -- I want to run this at run time
bool ValidateCmNodeTypeData (const CovarianceModel& cm)
{
	CovarianceModel::Node node;
	for (node=cm.GetFirstNode(); node!=cm.GetLastNode(); node++) {
		int thisNodeType=cm.GetNodeType(node);
		CmNodeTypeData data=GetCmNodeTypeData(thisNodeType);
		CovarianceModel::State state;
		int numCmStates=0;
		int numInsertPostStates=0;
		for (state=cm.GetFirstStateOfNode(node); state!=cm.GetLastStateOfNode(node); state++) {
			if (cm.GetStateType(state)==IL_st || cm.GetStateType(state)==IR_st) {
				numInsertPostStates++;
			}

			CovarianceModel::Node nextNode(node);
			nextNode++;

			// the curr node is not an ending node, then the next node should be physically next
			if (nextNode!=cm.GetLastNode() && !cm.IsBifuricationNode(node) && !cm.IsEndNode(node)) {
				CovarianceModel::State state=cm.GetFirstStateOfNode(node); // I don't trust GetLastStateOfNode, since I might change the alg
				while (cm.GetNode(state)==node) {
					state++;
				}
				// now that we're past 'node', we should be at 'nextNode', since nextNode should be physically next, as well as numerically
				if (cm.GetNode(state)!=nextNode) {
					assert(false);
					return false;
				}
			}

			if (nextNode!=cm.GetLastNode() && !cm.IsBifuricationNode(node)) { // last node doesn't link to anything & I'll ignore bifurication nodes
				CmNodeTypeData nextNodeData=GetCmNodeTypeData(cm.GetNodeType(nextNode));
				CovarianceModel::State expectedFirstChild;
				int numTransitions=cm.GetNumChildren(state);
				int expectedNumTransitions;

				if (cm.GetNodeType(node)==END_nd) {
					// special case for end nodes
					expectedNumTransitions=0;
				}
				else {
					if (numCmStates<data.numNormalStates) {
						// should link to IL/IR states + all normal states in the next node
						expectedNumTransitions=data.numInsertPostStates + nextNodeData.numNormalStates;
						expectedFirstChild=cm.GetFirstStateOfNode(node).PlusInt(data.numNormalStates);
					}
					else {
						// isa IL/IR state.  should link to itself and any remaining IL/IR states & all normal states im next node
						expectedNumTransitions=(data.numCmStates-numCmStates) + nextNodeData.numNormalStates;
						expectedFirstChild=state; // I should be my first child if I'm an IL/IR
					}
				}

				// verify we have the exact child states we expect
				if (numTransitions!=expectedNumTransitions) {
					assert(false);
					return false;
				}
				for (int childNum=0; childNum<numTransitions; childNum++) {
					if (cm.GetNthChildState(state,childNum)!=expectedFirstChild.PlusInt(childNum)) {
						assert(false);
						return false;
					}
				}
			}

			numCmStates++;
		}
		if (cm.GetNodeType(node)!=BIF_nd) { // ignore bifurications
			if (data.numCmStates!=numCmStates || data.numInsertPostStates!=numInsertPostStates || data.numNormalStates+data.numInsertPostStates!=data.numCmStates) {
				assert(false);
				return false;
			}
		}
	}
	return true;
}

int GetNumExtraMATPStatesVersusOriginalHmmBuildType (Cm2Hmm_HmmBuildType hmmType)
{
	switch (hmmType) {
		case HmmBuildType_Original:
			return 0;
		case HmmBuildType_separateMPandMLMR:
			return 1;
		case HmmBuildType_separateMPMLMRD:
			return 2;
	}
	assert(false);
	throw SimpleStringException("internal error %s:%d",__FILE__,__LINE__);
}

int GetNumNormalHmmStatesInNextLeftNode(const CovarianceModel& cm,CovarianceModel::Node currNode,Cm2Hmm_HmmBuildType hmmBuildType)
{
	// check that it's not an ending node, in which case there is no next node
	{
		int currNodeType=cm.GetNodeType(currNode);
		CmNodeTypeData data=GetCmNodeTypeData(currNodeType);
		if (data.isEndingNode) {
			return 0;
		}
	}

	// else, see what's next
	currNode++;
	while (1) {
		int currNodeType=cm.GetNodeType(currNode);
		CmNodeTypeData data=GetCmNodeTypeData(currNodeType);
		if (data.isEndingNode) {
			return 1;
		}

		if (data.numHmmLeftStates>0) {
			int num=data.numHmmLeftStates - data.numHmmLeftInsertStates;
			if (currNodeType==MATP_nd) {
				num += GetNumExtraMATPStatesVersusOriginalHmmBuildType(hmmBuildType);
			}
			return num;
		}

		currNode++;
	}
}

CovarianceModel::Node FindPreviousRightNode_NotAtStart(const CovarianceModel& cm,CovarianceModel::Node currNode)
{
	// check that it's not a starting node, in which case there is no previous node, and this violates our precondition
	assert(!GetCmNodeTypeData(cm.GetNodeType(currNode)).isStartingNode);

	// else see what's next
	currNode--;
	while (1) {
		int currNodeType=cm.GetNodeType(currNode);
		CmNodeTypeData data=GetCmNodeTypeData(currNodeType);
		if (data.isStartingNode) {
			return currNode;
		}

		if (data.numHmmLeftStates>0) {
			return currNode;
		}

		currNode++;
	}
}

// find the first (lowest-numbered) state in the given CM node that has a state on the Right side of the HMM
CovarianceModel::State GetFirstHmmRightState(const HmmAndBuildInfo& newHMM,const CovarianceModel& cm,CovarianceModel::Node prevRightNodeInCM)
{
	CovarianceModel::State firstCmState=cm.GetFirstStateOfNode(prevRightNodeInCM);
	CovarianceModel::State lastCmState=cm.GetLastStateOfNode(prevRightNodeInCM);

	CovarianceModel::State firstState=CovarianceModel::GetInvalidState();
	CovarianceModel::State cmState;
	for (cmState=firstCmState; cmState!=lastCmState; cmState++) {
		CovarianceModel::State hmmRightState=newHMM.cm2HmmState[cmState].hmmRightState;
		if (hmmRightState!=CovarianceModel::GetInvalidState()) {
			if (firstState==CovarianceModel::GetInvalidState()) {
				firstState=hmmRightState;
			}
			else {
				if (CovarianceModel::StateToInt(hmmRightState) < CovarianceModel::StateToInt(firstState)) {
					firstState=hmmRightState;
				}
			}
		}
	}

	assert(firstState!=CovarianceModel::GetInvalidState()); // should have got something
	return firstState;
}

// does all the work to set up the children of an HMM state, based on a CM source state
// the tricky thing about this function is that the children are relative to a start-to-end traversal of the HMM.  Thus, on the HMM that corresponds to Right-side emissions in the CM, we have to set up the children backwards, since CM nodes appear in the reverse order on the right-emitting side.
void SetupChildren (HmmAndBuildInfo& newHMM,const CovarianceModel& cm,CovarianceModel::State cmState,CovarianceModel::Node cmNode,HMM_Direction hmmDirection,CovarianceModel::Node cmFirstNode,CovarianceModel::Node cmLastNode,Cm2Hmm_HmmBuildType hmmBuildType)
{
	assert(hmmDirection==HMM_Left || hmmDirection==HMM_Right);

	CovarianceModel::State firstCmStateOfNode=cm.GetFirstStateOfNode(cmNode);

	CovarianceModel::State hmmState;
	if (hmmDirection==HMM_Left) {
		hmmState=newHMM.cm2HmmState[cmState].hmmLeftState;
	}
	else {
		hmmState=newHMM.cm2HmmState[cmState].hmmRightState;
	}

	int thisNodeType=cm.GetNodeType(cmNode);
	CmNodeTypeData thisNodeData=GetCmNodeTypeData(thisNodeType);

	assert(CovarianceModel::IsStateInRange(cmState,firstCmStateOfNode,firstCmStateOfNode.PlusInt(thisNodeData.numCmStates)));

	int thisNode_separate_extraStates=0;
	if (thisNodeType==MATP_nd) {
		thisNode_separate_extraStates=GetNumExtraMATPStatesVersusOriginalHmmBuildType(hmmBuildType);
	}

	if (hmmDirection==HMM_Left) {

		int numNormalStatesInNextLeftNode=GetNumNormalHmmStatesInNextLeftNode(cm,cmNode,hmmBuildType);
		if (CovarianceModel::IsStateInRange(cmState,firstCmStateOfNode,firstCmStateOfNode.PlusInt(thisNodeData.numNormalStates))) {
			// normal state on left side is connected to insert states in this node on the left side & normal states in next node in CM

			newHMM.hmm.SetFirstChild(hmmState,newHMM.cm2HmmState[firstCmStateOfNode].hmmLeftState.PlusInt(thisNodeData.numHmmLeftStates+thisNode_separate_extraStates-thisNodeData.numHmmLeftInsertStates)); // first insert state in this node
			newHMM.hmm.SetNumChildren(hmmState,thisNodeData.numHmmLeftInsertStates + numNormalStatesInNextLeftNode);
		}
		else {
			assert(cm.GetStateType(cmState)!=MP_st); // this shouldn't happen for MP state, since separateMPandMLMR logic is not applied here

			// insert state on left side is connected to itself and normal states in the next node in the CM
			newHMM.hmm.SetFirstChild(hmmState,newHMM.cm2HmmState[cmState].hmmLeftState);
			assert(thisNodeData.numHmmLeftInsertStates==1); // I'm an insert state, and I can only have <=1 IL_st, and <=1 IR_st
			newHMM.hmm.SetNumChildren(hmmState,1 + numNormalStatesInNextLeftNode);
		}
	}
	else {
		assert(hmmDirection==HMM_Right);
		if (CovarianceModel::IsStateInRange(cmState,firstCmStateOfNode,firstCmStateOfNode.PlusInt(thisNodeData.numNormalStates))) {
			// normal state on right side is connected to all states (insert and normal) of the _previous_ CM node that emits on the right
			if (thisNodeData.isStartingNode || cmNode==cmFirstNode) { // if it's an actual start node, or if it's acting like a start node (presumably because it's in DoLocal mode)
				// start node on right is connected to nothing (it's really an end node);
				newHMM.hmm.SetNoChildren(hmmState);
			}
			else {
				CovarianceModel::Node prevRightNode=FindPreviousRightNode_NotAtStart(cm,cmNode);
				//CovarianceModel::State firstStateOfPrevRightNode=cm.GetFirstStateOfNode(prevRightNode);
				int prevRightNodeType=cm.GetNodeType(prevRightNode);
				CmNodeTypeData prevRightNodeData=GetCmNodeTypeData(prevRightNodeType);
				newHMM.hmm.SetFirstChild(hmmState,GetFirstHmmRightState(newHMM,cm,prevRightNode));
				int extraStateForPrevNode=0;
				if (prevRightNodeType==MATP_nd) {
					extraStateForPrevNode=GetNumExtraMATPStatesVersusOriginalHmmBuildType(hmmBuildType);
				}
				newHMM.hmm.SetNumChildren(hmmState,prevRightNodeData.numHmmRightStates + extraStateForPrevNode);
			}
		}
		else {
			// insert state on right side is connected to itself & all normal states of the _current_ CM node that emit on the right
			newHMM.hmm.SetFirstChild(hmmState,GetFirstHmmRightState(newHMM,cm,cmNode));
			newHMM.hmm.SetNumChildren(hmmState,thisNodeData.numHmmRightStates+thisNode_separate_extraStates);
		}
	}
}

void ReverseMapCm2HmmState(HmmAndBuildInfo& newHMM,const CovarianceModel& cm)
{
	// set up reverse mapping 'hmm2CmStateVector'
	// WARNING: might be ambiguous since some HMM states are used by multiple CM states, but for now it doesn't matter.  So, I'm using a list of states
	newHMM.hmm2CmStateVector.clear(); // forget anything that's here already
	newHMM.hmm2CmStateVector.resize(newHMM.hmm.GetNumStates());
	CovarianceModel::State cmState;
	for (cmState=cm.GetFirstState(); cmState!=cm.GetLastState(); cmState++) {

		CovarianceModel::State hmmLeftState=newHMM.cm2HmmState[cmState].hmmLeftState;
		if (hmmLeftState!=CovarianceModel::GetInvalidState()) {
			newHMM.hmm2CmStateVector[hmmLeftState].push_back(cmState);
		}

		CovarianceModel::State hmmRightState=newHMM.cm2HmmState[cmState].hmmRightState;
		if (hmmRightState!=CovarianceModel::GetInvalidState()) {
			newHMM.hmm2CmStateVector[hmmRightState].push_back(cmState);
		}
	}
}


bool FirstNodeLooksLocallyClobbered(const CovarianceModel& cm) // do an explicit test, to allow us to feed a local-CM in and pretend it's global.  (Maybe it sounds crazy, but it's all in the name of scientific research.)
{
	CovarianceModel::State cmRootState=CovarianceModel::IntToState(0);
	for (int child=0; child<cm.GetNumChildren(cmRootState); child++) {
		if (cm.GetTransitionProbDirectly(cmRootState,child)!=0) {
			return false;
		}
	}
	return true;
}

/*
Converts the block of the CM in the closed interval [firstNode,lastNode]
to an HMM structure (i.e. not filling in the transition/emission weights)

firstNode must be a ROOT or BEGIN node
lastNode must be END, local END or BIFURICATION NODE
there cannot be any bifurication node in the range firstNode..(lastNode-1)
*/
void Cm2Hmm_Structurally_Block (HmmAndBuildInfo& newHMM,bool overlapMergedSubHmms,const CovarianceModel& sourceCM,CovarianceModel::Node firstNode,CovarianceModel::Node lastNode,ExtraCm2HmmInfo *extraCm2HmmInfo,Cm2Hmm_HmmBuildType hmmBuildType)
{
	assert(sourceCM.IsRootOrBeginNode(firstNode) || sourceCM.DoLocal() || FirstNodeLooksLocallyClobbered(sourceCM)); // firstNode must be BEGL, BEGR or ROOT, unless we're doing local, in which case we skip the root
	assert(sourceCM.IsEndNode(lastNode) || sourceCM.IsBifuricationNode(lastNode)); // last state must be END or BIF

	CovarianceModel::Node lastNodePlus1(lastNode);
	lastNodePlus1++;

	// other sanity checks
	CovarianceModel::Node cmNode;
	for (cmNode=firstNode; cmNode!=lastNode; cmNode++) {
		assert(cmNode==lastNode || !sourceCM.IsBifuricationNode(cmNode)); // can't have bifurication node, except at end
		assert(!(sourceCM.GetNodeType(cmNode)==DUMMY_nd)); // I thought those were removed by the time the CM was built
	}

	// work out how many states we'll need in the HMM
	int numHmmStates=0;
	for (CovarianceModel::Node node=firstNode; node!=lastNodePlus1; node++) {

		const CmNodeTypeData nodeTypeData=GetCmNodeTypeData(sourceCM.GetNodeType(node));
		numHmmStates += nodeTypeData.numHmmLeftStates + nodeTypeData.numHmmRightStates;

		if (sourceCM.GetNodeType(node)==MATP_nd) {
			numHmmStates += 2*GetNumExtraMATPStatesVersusOriginalHmmBuildType(hmmBuildType);
		}
	}
	if (overlapMergedSubHmms) {
		numHmmStates--; // the end states are over counted since they're both left & right
	}

	newHMM.hmm.Init(numHmmStates);
	Cm2HmmState nullCm2HmmState;
	nullCm2HmmState.hmmLeftState=CovarianceModel::GetInvalidState();
	nullCm2HmmState.hmmRightState=CovarianceModel::GetInvalidState();
	newHMM.cm2HmmState.assign(sourceCM.GetNumStates(),nullCm2HmmState);
	CovarianceModel::State firstHmmState=newHMM.hmm.GetFirstState(),lastHmmState=newHMM.hmm.GetLastState();

	// first set up the new HMM states & which CM states map to them (this info will be used in setting up the children of those states)
	for (cmNode=firstNode; cmNode!=lastNodePlus1; cmNode++) {

		assert(firstHmmState!=lastHmmState); // else we miscalculated somewhere

		int thisNodeType=sourceCM.GetNodeType(cmNode);
		CmNodeTypeData thisNodeData=GetCmNodeTypeData(thisNodeType);
		CovarianceModel::State firstCmState=sourceCM.GetFirstStateOfNode(cmNode);
		CovarianceModel::State lastCmState=sourceCM.GetLastStateOfNode(cmNode);

		if (thisNodeData.isEndingNode) {
			if (overlapMergedSubHmms) {
				lastHmmState--;
				assert(firstHmmState==lastHmmState); // end in the middle
				newHMM.hmm.SetStateType(firstHmmState,PASSTHRU_st);
				newHMM.cm2HmmState[firstCmState].hmmLeftState=firstHmmState;
				newHMM.cm2HmmState[firstCmState].hmmRightState=firstHmmState;
				assert(thisNodeData.numCmStates==1); // should have checked this already, but whatever
				newHMM.leftToRightPassthruState=firstHmmState;
				firstHmmState++;
			}
			else {
				lastHmmState--;
				newHMM.hmm.SetStateType(firstHmmState,PASSTHRU_st);
				newHMM.hmm.SetStateType(lastHmmState,PASSTHRU_st);
				newHMM.cm2HmmState[firstCmState].hmmLeftState=firstHmmState;
				newHMM.cm2HmmState[firstCmState].hmmRightState=lastHmmState;
				assert(thisNodeData.numCmStates==1); // should have checked this already, but whatever
				newHMM.leftToRightPassthruState=firstHmmState;
				firstHmmState++;
				assert(firstHmmState==lastHmmState); // end in the middle
			}
		}
		else {
			// do normal states of node
			switch (sourceCM.GetNodeType(cmNode)) {
				case ROOT_nd:
				case BEGL_nd:
				case BEGR_nd:
					// left side start
					newHMM.hmm.SetStateType(firstHmmState,S_st);
					newHMM.cm2HmmState[firstCmState].hmmLeftState=firstHmmState;
					firstHmmState++;
					// right side start
					lastHmmState--;
					newHMM.hmm.SetStateType(lastHmmState,E_st);
					newHMM.cm2HmmState[firstCmState].hmmRightState=lastHmmState;
					break;
				case MATL_nd:
					// ML_st
					newHMM.hmm.SetStateType(firstHmmState,ML_st);
					newHMM.cm2HmmState[firstCmState].hmmLeftState=firstHmmState;
					newHMM.cm2HmmState[firstCmState].hmmRightState=CovarianceModel::GetInvalidState();
					firstHmmState++;
					// D_st
					newHMM.hmm.SetStateType(firstHmmState,D_st);
					newHMM.cm2HmmState[firstCmState.PlusInt(1)].hmmLeftState=firstHmmState;
					newHMM.cm2HmmState[firstCmState.PlusInt(1)].hmmRightState=CovarianceModel::GetInvalidState();
					firstHmmState++;
					if (everyCmNodeHasLeftAndRightHmmNodes) {
						lastHmmState--;
						newHMM.hmm.SetStateType(lastHmmState,PASSTHRU_st);
						newHMM.cm2HmmState[firstCmState].hmmRightState=lastHmmState;
						newHMM.cm2HmmState[firstCmState.PlusInt(1)].hmmRightState=lastHmmState;
					}
					break;
				case MATR_nd:
					// MR_st
					lastHmmState--;
					newHMM.hmm.SetStateType(lastHmmState,ML_st);
					newHMM.cm2HmmState[firstCmState].hmmLeftState=CovarianceModel::GetInvalidState();
					newHMM.cm2HmmState[firstCmState].hmmRightState=lastHmmState;
					// D_st
					lastHmmState--;
					newHMM.hmm.SetStateType(lastHmmState,D_st);
					newHMM.cm2HmmState[firstCmState.PlusInt(1)].hmmLeftState=CovarianceModel::GetInvalidState();
					newHMM.cm2HmmState[firstCmState.PlusInt(1)].hmmRightState=lastHmmState;
					if (everyCmNodeHasLeftAndRightHmmNodes) {
						newHMM.hmm.SetStateType(firstHmmState,PASSTHRU_st);
						newHMM.cm2HmmState[firstCmState].hmmLeftState=firstHmmState;
						newHMM.cm2HmmState[firstCmState.PlusInt(1)].hmmLeftState=firstHmmState;
						firstHmmState++;
					}
					break;
				case MATP_nd:
					switch (hmmBuildType) {
						case HmmBuildType_Original: // also called "Compacted" or "type 0"
							// normal case

							// left match
							newHMM.hmm.SetStateType(firstHmmState,ML_st);
							// left delete
							newHMM.hmm.SetStateType(firstHmmState.PlusInt(1),D_st);
							// right stuff
							lastHmmState--;
							lastHmmState--;
							// right match
							newHMM.hmm.SetStateType(lastHmmState,ML_st);
							// right delete
							newHMM.hmm.SetStateType(lastHmmState.PlusInt(1),D_st);

							// hook up CM to equivalent HMM states
							// MP_st
							newHMM.cm2HmmState[firstCmState.PlusInt(0)].hmmLeftState=firstHmmState;
							newHMM.cm2HmmState[firstCmState.PlusInt(0)].hmmRightState=lastHmmState;
							// ML_st
							newHMM.cm2HmmState[firstCmState.PlusInt(1)].hmmLeftState=firstHmmState;
							newHMM.cm2HmmState[firstCmState.PlusInt(1)].hmmRightState=lastHmmState.PlusInt(1);
							// MR_st
							newHMM.cm2HmmState[firstCmState.PlusInt(2)].hmmLeftState=firstHmmState.PlusInt(1);
							newHMM.cm2HmmState[firstCmState.PlusInt(2)].hmmRightState=lastHmmState;
							// D_st
							newHMM.cm2HmmState[firstCmState.PlusInt(3)].hmmLeftState=firstHmmState.PlusInt(1);
							newHMM.cm2HmmState[firstCmState.PlusInt(3)].hmmRightState=lastHmmState.PlusInt(1);

							firstHmmState++;
							firstHmmState++;
							break;
						case HmmBuildType_separateMPandMLMR: // also called "expanded" or "type 1"

							// with extra states

							// left match MP
							newHMM.hmm.SetStateType(firstHmmState,ML_st);
							// left match ML
							newHMM.hmm.SetStateType(firstHmmState.PlusInt(1),ML_st);
							// left delete
							newHMM.hmm.SetStateType(firstHmmState.PlusInt(2),D_st);
							// right stuff
							lastHmmState--;
							lastHmmState--;
							lastHmmState--;
							// right match MP
							newHMM.hmm.SetStateType(lastHmmState,ML_st);
							// right match MR
							newHMM.hmm.SetStateType(lastHmmState.PlusInt(1),ML_st);
							// right delete
							newHMM.hmm.SetStateType(lastHmmState.PlusInt(2),D_st);

							// hook up CM to equivalent HMM states
							// MP_st
							newHMM.cm2HmmState[firstCmState.PlusInt(0)].hmmLeftState=firstHmmState;
							newHMM.cm2HmmState[firstCmState.PlusInt(0)].hmmRightState=lastHmmState;
							// ML_st
							newHMM.cm2HmmState[firstCmState.PlusInt(1)].hmmLeftState=firstHmmState.PlusInt(1);
							newHMM.cm2HmmState[firstCmState.PlusInt(1)].hmmRightState=lastHmmState.PlusInt(2);
							// MR_st
							newHMM.cm2HmmState[firstCmState.PlusInt(2)].hmmLeftState=firstHmmState.PlusInt(2);
							newHMM.cm2HmmState[firstCmState.PlusInt(2)].hmmRightState=lastHmmState.PlusInt(1);
							// D_st
							newHMM.cm2HmmState[firstCmState.PlusInt(3)].hmmLeftState=firstHmmState.PlusInt(2);
							newHMM.cm2HmmState[firstCmState.PlusInt(3)].hmmRightState=lastHmmState.PlusInt(2);

							firstHmmState++;
							firstHmmState++;
							firstHmmState++;
							break;
						case HmmBuildType_separateMPMLMRD:

							// with more extra states

							// left match MP
							newHMM.hmm.SetStateType(firstHmmState,ML_st);
							// left match ML
							newHMM.hmm.SetStateType(firstHmmState.PlusInt(1),ML_st);
							// left delete MR
							newHMM.hmm.SetStateType(firstHmmState.PlusInt(2),D_st);
							// left delete D
							newHMM.hmm.SetStateType(firstHmmState.PlusInt(3),D_st);
							// right stuff
							lastHmmState--;
							lastHmmState--;
							lastHmmState--;
							lastHmmState--;
							// right match MP
							newHMM.hmm.SetStateType(lastHmmState,ML_st);
							// right match MR
							newHMM.hmm.SetStateType(lastHmmState.PlusInt(1),ML_st);
							// right delete ML
							newHMM.hmm.SetStateType(lastHmmState.PlusInt(2),D_st);
							// right delete D
							newHMM.hmm.SetStateType(lastHmmState.PlusInt(3),D_st);

							// hook up CM to equivalent HMM states
							// MP_st
							newHMM.cm2HmmState[firstCmState.PlusInt(0)].hmmLeftState=firstHmmState;
							newHMM.cm2HmmState[firstCmState.PlusInt(0)].hmmRightState=lastHmmState;
							// ML_st
							newHMM.cm2HmmState[firstCmState.PlusInt(1)].hmmLeftState=firstHmmState.PlusInt(1);
							newHMM.cm2HmmState[firstCmState.PlusInt(1)].hmmRightState=lastHmmState.PlusInt(2);
							// MR_st
							newHMM.cm2HmmState[firstCmState.PlusInt(2)].hmmLeftState=firstHmmState.PlusInt(2);
							newHMM.cm2HmmState[firstCmState.PlusInt(2)].hmmRightState=lastHmmState.PlusInt(1);
							// D_st
							newHMM.cm2HmmState[firstCmState.PlusInt(3)].hmmLeftState=firstHmmState.PlusInt(3);
							newHMM.cm2HmmState[firstCmState.PlusInt(3)].hmmRightState=lastHmmState.PlusInt(3);

							firstHmmState++;
							firstHmmState++;
							firstHmmState++;
							firstHmmState++;
							break;
					}
					break;
				default:
					assert(false);
					break;
			}

			// do insert states
			CovarianceModel::State firstCmInsertState=firstCmState.PlusInt(thisNodeData.numNormalStates);
			CovarianceModel::State lastCmInsertState=lastCmState;
			CovarianceModel::State cmState;
			for (cmState=firstCmInsertState; cmState!=lastCmInsertState; cmState++) {
				if (sourceCM.GetStateType(cmState)==IL_st) {
					newHMM.hmm.SetStateType(firstHmmState,IL_st);
					newHMM.cm2HmmState[cmState].hmmLeftState=firstHmmState;
					newHMM.cm2HmmState[cmState].hmmRightState=CovarianceModel::GetInvalidState();
					firstHmmState++;
				}
				else {
					assert(sourceCM.GetStateType(cmState)==IR_st);

					lastHmmState--;
					newHMM.hmm.SetStateType(lastHmmState,IL_st);
					newHMM.cm2HmmState[cmState].hmmLeftState=CovarianceModel::GetInvalidState();
					newHMM.cm2HmmState[cmState].hmmRightState=lastHmmState;
				}
			}
		}
	}

	// now figure out child structure
	for (cmNode=firstNode; cmNode!=lastNodePlus1; cmNode++) {

		CovarianceModel::State cmState;
		for (cmState=sourceCM.GetFirstStateOfNode(cmNode); cmState!=sourceCM.GetLastStateOfNode(cmNode); cmState++) {

			if (newHMM.cm2HmmState[cmState].hmmLeftState!=CovarianceModel::GetInvalidState()) {
				SetupChildren (newHMM,sourceCM,cmState,cmNode,HMM_Left,firstNode,lastNode,hmmBuildType);
			}
			if (newHMM.cm2HmmState[cmState].hmmRightState!=CovarianceModel::GetInvalidState()) {
				SetupChildren (newHMM,sourceCM,cmState,cmNode,HMM_Right,firstNode,lastNode,hmmBuildType);
			}
		}
	}

	if (!overlapMergedSubHmms) {
		// fix up left-to-right connection
		newHMM.hmm.SetFirstChild(newHMM.leftToRightPassthruState,newHMM.leftToRightPassthruState.PlusInt(1)); // this gets set right, but I want to be explicit.  (No way!  This is family code, buster.)
		newHMM.hmm.SetNumChildren(newHMM.leftToRightPassthruState,1); // this gets set wrong, so I have to be explicit.  (Look buster, you can get your point across without using naughty words.)
	}

	ReverseMapCm2HmmState(newHMM,sourceCM);
}

void DumpHmm (FILE *file,const HmmAndBuildInfo& newHMM,const CovarianceModel& cm)
{
	throw SimpleStringException("DumpHmm not implemented");
}

void DumpHmmAndBuildInfo (FILE *file,const HmmAndBuildInfo& newHMM,const CovarianceModel& cm)
{
	DumpHmm(file,newHMM,cm);

	CovarianceModel::State state;
	fprintf(file,"CM had %u states\n",newHMM.cm2HmmState.size());
	for (state=CovarianceModel::IntToState(0); state!=CovarianceModel::IntToState((int)newHMM.cm2HmmState.size()); state++) {
		fprintf(file,"CM state %d maps to HMM: left = %d, right = %d\n",CovarianceModel::StateToInt(state),CovarianceModel::StateToInt(newHMM.cm2HmmState[state].hmmLeftState),CovarianceModel::StateToInt(newHMM.cm2HmmState[state].hmmRightState));
	}
}

CovarianceModel::Node FindFirstEndingNode (const CovarianceModel& sourceCM,CovarianceModel::Node firstNode)
{
	CovarianceModel::Node lastNode=sourceCM.GetLastNode();
	CovarianceModel::Node currNode;
	for (currNode=firstNode; currNode!=lastNode; currNode++) {
		int thisNodeType=sourceCM.GetNodeType(currNode);
		CmNodeTypeData thisNodeData=GetCmNodeTypeData(thisNodeType);

		if (thisNodeData.isEndingNode) {
			return currNode;
		}
	}
	return lastNode;
}

void SetupReverseMapping(ScoreVariablesInfo& scoreVariablesInfo,const HmmAndBuildInfo& newHMM,const CovarianceModel& sourceCM)
{
	TransitionOrEmissionInfo dummy;
	dummy.isUsed=false;
	scoreVariablesInfo.globalVariableToTransitionOrEmissionVector.assign(scoreVariablesInfo.numVariables,dummy);

	CovarianceModel::State hmmState;
	for (hmmState=newHMM.hmm.GetFirstState(); hmmState!=newHMM.hmm.GetLastState(); hmmState++) {

		for (unsigned int i=0;	i<scoreVariablesInfo.transitionToVariableNumVector[hmmState].size(); i++) {
			int globalVar=scoreVariablesInfo.transitionToVariableNumVector[hmmState][i];
			if (globalVar!=-1) {
				CovarianceModel::State toState=newHMM.hmm.GetNthChildState(hmmState,i);

				assert(!scoreVariablesInfo.globalVariableToTransitionOrEmissionVector[globalVar].isUsed);
				scoreVariablesInfo.globalVariableToTransitionOrEmissionVector[globalVar].isUsed=true;
				scoreVariablesInfo.globalVariableToTransitionOrEmissionVector[globalVar].type=TransitionOrEmissionInfo::Transition;
				scoreVariablesInfo.globalVariableToTransitionOrEmissionVector[globalVar].edgeInfo.fromState=hmmState;
				scoreVariablesInfo.globalVariableToTransitionOrEmissionVector[globalVar].edgeInfo.toState=toState;
				scoreVariablesInfo.globalVariableToTransitionOrEmissionVector[globalVar].edgeInfo.childNum=(int)i;
			}
		}
		for (unsigned int i=0;	i<scoreVariablesInfo.emissionToVariableNumVector[hmmState].size(); i++) {
			int globalVar=scoreVariablesInfo.emissionToVariableNumVector[hmmState][i];
			if (globalVar!=-1) {

				assert(!scoreVariablesInfo.globalVariableToTransitionOrEmissionVector[globalVar].isUsed);
				scoreVariablesInfo.globalVariableToTransitionOrEmissionVector[globalVar].isUsed=true;
				scoreVariablesInfo.globalVariableToTransitionOrEmissionVector[globalVar].type=TransitionOrEmissionInfo::Emission;
				scoreVariablesInfo.globalVariableToTransitionOrEmissionVector[globalVar].emissionInfo.state=hmmState;
				scoreVariablesInfo.globalVariableToTransitionOrEmissionVector[globalVar].emissionInfo.nuc=(int)i;
			}
		}

		if (newHMM.hmm.DoLocal()) {

			int globalVar;

			globalVar=scoreVariablesInfo.leftwardBeginToVariableNumberVector[hmmState];
			if (globalVar!=-1) {
				assert(!scoreVariablesInfo.globalVariableToTransitionOrEmissionVector[globalVar].isUsed);
				scoreVariablesInfo.globalVariableToTransitionOrEmissionVector[globalVar].isUsed=true;
				scoreVariablesInfo.globalVariableToTransitionOrEmissionVector[globalVar].type=TransitionOrEmissionInfo::LeftwardBegin;
				scoreVariablesInfo.globalVariableToTransitionOrEmissionVector[globalVar].hmmState=hmmState;
			}

			globalVar=scoreVariablesInfo.rightwardBeginToVariableNumberVector[hmmState];
			if (globalVar!=-1) {
				assert(!scoreVariablesInfo.globalVariableToTransitionOrEmissionVector[globalVar].isUsed);
				scoreVariablesInfo.globalVariableToTransitionOrEmissionVector[globalVar].isUsed=true;
				scoreVariablesInfo.globalVariableToTransitionOrEmissionVector[globalVar].type=TransitionOrEmissionInfo::RightwardBegin;
				scoreVariablesInfo.globalVariableToTransitionOrEmissionVector[globalVar].hmmState=hmmState;
			}

			for (int link=0; link<newHMM.hmm.GetNumEndscLinksToLeft(hmmState); link++) {
				globalVar=scoreVariablesInfo.localEndToVariableNumberVector[hmmState][link];
				if (globalVar!=-1) {
					assert(!scoreVariablesInfo.globalVariableToTransitionOrEmissionVector[globalVar].isUsed);
					scoreVariablesInfo.globalVariableToTransitionOrEmissionVector[globalVar].isUsed=true;
					scoreVariablesInfo.globalVariableToTransitionOrEmissionVector[globalVar].type=TransitionOrEmissionInfo::LocalEnd;
					scoreVariablesInfo.globalVariableToTransitionOrEmissionVector[globalVar].hmmState=hmmState;
					scoreVariablesInfo.globalVariableToTransitionOrEmissionVector[globalVar].localEndLinkNum=link;
				}
			}
		}
	}
}

void SetupTransitionAndEmissionVariables(ScoreVariablesInfo& scoreVariablesInfo,HmmAndBuildInfo& newHMM,const CovarianceModel& sourceCM)
{
	scoreVariablesInfo.transitionToVariableNumVector.resize(newHMM.hmm.GetNumStates());
	scoreVariablesInfo.emissionToVariableNumVector.resize(newHMM.hmm.GetNumStates());

	if (newHMM.hmm.DoLocal()) {
		scoreVariablesInfo.leftwardBeginToVariableNumberVector.resize(newHMM.hmm.GetNumStates());
		scoreVariablesInfo.rightwardBeginToVariableNumberVector.resize(newHMM.hmm.GetNumStates());
		scoreVariablesInfo.localEndToVariableNumberVector.resize(newHMM.hmm.GetNumStates());
	}

	int nextVariableNum=0;
	CovarianceModel::State hmmState;
	for (hmmState=newHMM.hmm.GetFirstState(); hmmState!=newHMM.hmm.GetLastState(); hmmState++) {

		// transition variables
		switch (newHMM.hmm.GetNumChildren(hmmState)) {
			case 0:
				// nothing to do
				break;
			case 1:
				// with only 1 transition, it doesn't matter -- just set it to 0 (this case can arise with PASSTHRU nodes)
				newHMM.hmm.SetTransitionLogScore(hmmState,0,0.0);
				// resize to 1, for convenience
				scoreVariablesInfo.transitionToVariableNumVector[hmmState].resize(1);
				scoreVariablesInfo.transitionToVariableNumVector[hmmState][0]=-1;
				break;
			default:
				// the interesting case -- multiple transitions
				scoreVariablesInfo.transitionToVariableNumVector[hmmState].resize(newHMM.hmm.GetNumChildren(hmmState));
				for (int childNum=0; childNum<newHMM.hmm.GetNumChildren(hmmState); childNum++) {
					if (newHMM.hmm.GetNthChildState(hmmState,childNum)==hmmState) {

						// self-loops are easy -- they must have the same score as the corresponding loop in the CM

						assert(newHMM.hmm.IsInsertState(hmmState)); // this is the only self loops I'm expecting
						assert(newHMM.hmm2CmStateVector[hmmState].size()==1); // insert state should only map to one state
						CovarianceModel::State cmState=newHMM.hmm2CmStateVector[hmmState].front();
						assert(sourceCM.GetNthChildState(cmState,0)==cmState && sourceCM.IsInsertState(cmState)); // this is what I expect based on the stereotyped structure of CMs in the infernal code -- otherwise, I'd have to search for the self-loop

						// NOTE: the following 2 lines are a workaround to a scary (possible) g++ optimizer bug that makes this code crash in release mode, but do fine in debug.  By using the variable 'tsc' (rather than just putting 'sourceCM.GetNthChildTsc(cmState,0)' directly as a param), it doesn't crash, & appears to work.  This happens for RF00032 and RF00016, and I'd guess everything else.  The code appears to work properly with this fix.  I've checked hmm-dump.txt for RF00032 and it achieves the same inflations; the actual weights are different by around 1e-6, but I already know the exact solution is different on different platforms.  Pretty scary, because I'm really not sure what caused this bug, and I can't think of an obvious, benign reason for it to suddenly show up when I changed to code to do structure entirely first, and then solve the scores.  Oh well.
						float tsc=sourceCM.GetNthChildTsc(cmState,0);
						newHMM.hmm.SetTransitionLogScore(hmmState,childNum,tsc);
						scoreVariablesInfo.transitionToVariableNumVector[hmmState][childNum]=-1;

						// while I'm at it, might as well set the emission scores, which should also be the same
						scoreVariablesInfo.emissionToVariableNumVector[hmmState].resize(Alphabet_size);
						for (int nuc=0; nuc<Alphabet_size; nuc++) {
							scoreVariablesInfo.emissionToVariableNumVector[hmmState][nuc]=-1;
							newHMM.hmm.SetSingletEmissionLogScore(hmmState,nuc,sourceCM.GetSingletEmissionScore(cmState,nuc));
						}
					}
					else {
						// add a new variable
						scoreVariablesInfo.transitionToVariableNumVector[hmmState][childNum]=nextVariableNum;
						nextVariableNum++;
					}
				}
				break;
		}

		// emission variables
		// we already did IL_st, IR_st, since they're self-looping, and we can just copy the probs from the CM.  Now, we have to do ML_st and MR_st.
		assert(newHMM.hmm.GetNumSymbolsEmitted(hmmState)==0 || newHMM.hmm.GetNumSymbolsEmitted(hmmState)==1); // HMMs don't have MP_st
		if (newHMM.hmm.GetNumSymbolsEmitted(hmmState)==1 && !newHMM.hmm.IsInsertState(hmmState)) {
			scoreVariablesInfo.emissionToVariableNumVector[hmmState].resize(Alphabet_size);
			for (int nuc=0; nuc<Alphabet_size; nuc++) {
				scoreVariablesInfo.emissionToVariableNumVector[hmmState][nuc]=nextVariableNum;
				nextVariableNum++;
			}
		}

		// local begin/end variables
		if (newHMM.hmm.DoLocal()) {

			if (newHMM.hmm.GetLeftwardBeginsc(hmmState)==(float)IMPOSSIBLE) {
				scoreVariablesInfo.leftwardBeginToVariableNumberVector[hmmState]=-1;
			}
			else {
				scoreVariablesInfo.leftwardBeginToVariableNumberVector[hmmState]=nextVariableNum;
				nextVariableNum++;
			}

			if (newHMM.hmm.GetRightwardBeginsc(hmmState)==(float)IMPOSSIBLE) {
				scoreVariablesInfo.rightwardBeginToVariableNumberVector[hmmState]=-1;
			}
			else {
				scoreVariablesInfo.rightwardBeginToVariableNumberVector[hmmState]=nextVariableNum;
				nextVariableNum++;
			}

			scoreVariablesInfo.localEndToVariableNumberVector[hmmState].resize(newHMM.hmm.GetNumEndscLinksToLeft(hmmState));
			for (int link=0; link<newHMM.hmm.GetNumEndscLinksToLeft(hmmState); link++) {
				scoreVariablesInfo.localEndToVariableNumberVector[hmmState][link]=nextVariableNum;
				nextVariableNum++;
			}
		}
	}

	scoreVariablesInfo.numVariables=nextVariableNum;
}


int BindLocalVariable (int globalVariableNum,vector<int>& globalToLocalVariables,int& numLocalVariables)
{
	int localVariableNum=globalToLocalVariables[globalVariableNum];
	if (localVariableNum==-1) {
		localVariableNum=numLocalVariables;
		globalToLocalVariables[globalVariableNum]=localVariableNum;
		numLocalVariables++;
	}
	return localVariableNum;
}

TemporarilyModifyInequality::TemporarilyModifyInequality (Inequality& _inequalitySoFar)
: inequalitySoFar(_inequalitySoFar)
{
	numVariablesAdded=0;
	startingScore=inequalitySoFar.rhs;
	startingWeight=inequalitySoFar.weight;
	starting_sumOfConstantsInHmm=inequalitySoFar.sumOfConstantsInHmm;
	starting_hmmInsertStatesInPath_size=inequalitySoFar.hmmInsertStatesInPath.size();
	for (int i=0; i<MAXABET; i++) {
		startingNucEmitCount[i]=inequalitySoFar.nucEmitCount[i];
	}
}
TemporarilyModifyInequality::~TemporarilyModifyInequality ()
{
	for (int i=0; i<numVariablesAdded; i++) {
		assert(!inequalitySoFar.lhs.empty());
		inequalitySoFar.lhs.pop_back();
	}

	inequalitySoFar.rhs=startingScore;
	inequalitySoFar.weight=startingWeight;
	inequalitySoFar.sumOfConstantsInHmm=starting_sumOfConstantsInHmm;

	assert(inequalitySoFar.hmmInsertStatesInPath.size()>=starting_hmmInsertStatesInPath_size); // how could it get smaller as we do more of the path?!
	while (inequalitySoFar.hmmInsertStatesInPath.size()>starting_hmmInsertStatesInPath_size) {
		inequalitySoFar.hmmInsertStatesInPath.pop_back();
	}

	for (int i=0; i<MAXABET; i++) {
		inequalitySoFar.nucEmitCount[i]=startingNucEmitCount[i];
	}
}
void TemporarilyModifyInequality::PushInsertState (InfernalHmm::State insertState)
{
	inequalitySoFar.hmmInsertStatesInPath.push_back(insertState);
	assert(inequalitySoFar.hmmInsertStatesInPath.size()<=2);
}
void TemporarilyModifyInequality::AddScore (float addToScore)
{
	inequalitySoFar.rhs += addToScore;
}
void TemporarilyModifyInequality::MultiplyWeight (double mult)
{
	inequalitySoFar.weight *= mult;
}
void TemporarilyModifyInequality::AddVariable (int globalVariableNum,vector<int>& globalToLocalVariables,int& numLocalVariables)
{
	InequalityTerm newTerm;
	int localVariableNum=globalToLocalVariables[globalVariableNum];
	if (localVariableNum==-1) {
		localVariableNum=numLocalVariables;
		globalToLocalVariables[globalVariableNum]=localVariableNum;
		numLocalVariables++;
	}
	newTerm.variableNum=localVariableNum;
	inequalitySoFar.lhs.push_back(newTerm);
	numVariablesAdded++;
}

void Cm2Hmm_MakeInequalitiesForPath_HmmTransition(TemporarilyModifyInequality& temporarilyModifyInequality,const CovarianceModel::State hmmState,CovarianceModel::State& hmmNextState,const HmmAndBuildInfo& newHMM,const ScoreVariablesInfo& scoreVariablesInfo,vector<int>& globalToLocalVariables,int& numLocalVariables,WeightedInequalitiesInfo *weightedInequalitiesInfo)
{
	if (hmmNextState==CovarianceModel::GetInvalidState()) {
		// stays the same, and no transition cost
		hmmNextState=hmmState;
	}
	else {
		if (hmmNextState==hmmState) {
			// nothing to do -- this HMM state didn't change on that CM transition
		}
		else {
			// check if we're going to an insert state
			if (newHMM.hmm.IsInsertState(hmmNextState)) {
				temporarilyModifyInequality.PushInsertState(hmmNextState);
			}

			// update inequality to take into account this transition
			CovarianceModel::State fromState=hmmState,toState=hmmNextState;
			if (fromState>=toState) {
				// we must be on the right side of the HMM, so we're going backwards
				std::swap(fromState,toState);
			}
			int hmmChildNum=newHMM.hmm.GetChildNum_Slow(fromState,toState);
			int globalVariableNum=scoreVariablesInfo.transitionToVariableNumVector[fromState][hmmChildNum];
			if (globalVariableNum==-1) {
				// oh, don't bother
			}
			else {
				// okay, add the variable for this transition
				temporarilyModifyInequality.AddVariable(globalVariableNum,globalToLocalVariables,numLocalVariables);

				// and weight the inequality 
				if (weightedInequalitiesInfo!=NULL) {
					temporarilyModifyInequality.MultiplyWeight(weightedInequalitiesInfo->transitionCounter->GetTransitionFrequency_Unreversed(fromState,toState));
				}
			}
		}
	}
}


void DumblyWorkOutHmmStates(CovarianceModel::State& hmmFirstNormalState,CovarianceModel::State& hmmLastNormalState,const HmmAndBuildInfo& newHMM,CovarianceModel::Node cmNode,const CovarianceModel& sourceCM,const CovarianceModel::State Cm2HmmState::*leftOrRight)
{
	// do this the dumb way that works;

	CmNodeTypeData nodeData=GetCmNodeTypeData(sourceCM.GetNodeType(cmNode));
	CovarianceModel::State cmFirstNormalState=sourceCM.GetFirstStateOfNode(cmNode);
	CovarianceModel::State cmLastNormalState=cmFirstNormalState.PlusInt(nodeData.numNormalStates);

	hmmFirstNormalState=newHMM.cm2HmmState[cmFirstNormalState].*leftOrRight;
	hmmLastNormalState=newHMM.cm2HmmState[cmFirstNormalState].*leftOrRight;

	CovarianceModel::State cmRootState;
	for (cmRootState=cmFirstNormalState; sourceCM.IsStateInRange(cmRootState,cmFirstNormalState,cmLastNormalState); cmRootState++) {
		CovarianceModel::State thisHmmState=newHMM.cm2HmmState[cmRootState].*leftOrRight;
		hmmFirstNormalState=std::min(hmmFirstNormalState,thisHmmState);
		hmmLastNormalState=std::max(hmmLastNormalState,thisHmmState);
	}
	// make half-open
	hmmLastNormalState++;
	assert(hmmFirstNormalState!=CovarianceModel::GetInvalidState() && hmmLastNormalState!=CovarianceModel::GetInvalidState() && hmmFirstNormalState!=CovarianceModel::GetInvalidState() && hmmLastNormalState!=CovarianceModel::GetInvalidState());
	assert(hmmFirstNormalState<hmmLastNormalState); // every node should have at least 1 "normal" state
}

// this function is like 'Cm2Hmm_MakeInequalitiesForPath' (see below), but does things for the local case.  In the local case, there's really no need to walk paths recursively, because it's pretty specific to a state
// NOTE: it's actually an assumption that paths aren't necessary; this is true for the current version of CMs, but might become false if CMs are changed.  The key requirement is that local begins/ends can only be applied to split-set states.
void Cm2Hmm_MakeInequalitiesForLocals(HmmAndBuildInfo& newHMM,const CovarianceModel& sourceCM,const ScoreVariablesInfo& scoreVariablesInfo,const CovarianceModel::Node cmStartPathNode,const CovarianceModel::Node cmEndPathNode,InequalityList& inequalityList,vector<int>& globalToLocalVariables,int& numLocalVariables)
{
	if (!newHMM.hmm.DoLocal()) {
		return;
	}

	CovarianceModel::State cmState;
	cmState=sourceCM.GetFirstStateOfNode(cmStartPathNode);
	while (cmState!=sourceCM.GetFirstStateOfNode(cmEndPathNode)) {

		InfernalHmm::State hmmLeftState=newHMM.cm2HmmState[cmState].hmmLeftState;
		InfernalHmm::State hmmRightState=newHMM.cm2HmmState[cmState].hmmRightState;

		// local begin: [left] + [right] >= CM-beginsc
		float cmLocalBeginsc=sourceCM.GetBeginsc(cmState);
		if (cmLocalBeginsc!=(float)IMPOSSIBLE) {
			int leftGlobalVar=-1,rightGlobalVar=-1;
			if (hmmLeftState!=InfernalHmm::GetInvalidState()) {
				leftGlobalVar=scoreVariablesInfo.leftwardBeginToVariableNumberVector[hmmLeftState];
			}
			if (hmmRightState!=InfernalHmm::GetInvalidState()) {
				rightGlobalVar=scoreVariablesInfo.rightwardBeginToVariableNumberVector[hmmRightState];
			}
			if (leftGlobalVar!=-1 && rightGlobalVar!=-1) {
				int leftLocalVar=BindLocalVariable(leftGlobalVar,globalToLocalVariables,numLocalVariables);
				int rightLocalVar=BindLocalVariable(rightGlobalVar,globalToLocalVariables,numLocalVariables);

				Inequality ineq;
				ineq.inequalityType=IneqType_GE;
				ineq.weight=1.0;
				ineq.nucEmitCount[0]=ineq.nucEmitCount[1]=ineq.nucEmitCount[2]=ineq.nucEmitCount[3]=0;
				ineq.sumOfConstantsInHmm=0;
				ineq.rhs=cmLocalBeginsc;

				InequalityTerm leftTerm,rightTerm;
				leftTerm.variableNum=leftLocalVar;
				rightTerm.variableNum=rightLocalVar;
				ineq.lhs.push_back(leftTerm);
				ineq.lhs.push_back(rightTerm);

				assert(!ineq.lhs.empty());
				inequalityList.push_back(ineq);
			}
		}

		// local end:
		// It's easy to set the HMM local end score once we know the emit scores.  Unfortunately, in order
		// to set it up as a Linear Program (or Non-Linear Program) we need to add an endsc variable,
		// and several inequalities, because the actual equations uses min/max operations -- if we didn't add the
		// constraints, we may get impossible values
		// Anyway, if cmState is non emitting, then [emitsc]==CM endsc (this is an easy case, although I'm pretty sure it never arises, since delete states can never be the consensus CM position)
		// if cmState emits left or right only, then for all nucs X: [hmm emit X] + [hmm endsc] >= [CM emit X] + [CM endsc]
		// if cmState is MP, then for all nucs X,Y: [hmm emit X] + [hmm emit Y] + [hmm endsc] >= [CM emit X-Y pair] + [CM endsc]
		float cmEndsc=sourceCM.GetEndsc(cmState);
		if (cmEndsc!=(float)IMPOSSIBLE && scoreVariablesInfo.localEndToVariableNumberVector[hmmRightState].size()>0) { // otherwise not set

			bool foundIt=false;
			for (int link=0; link<newHMM.hmm.GetNumEndscLinksToLeft(hmmRightState); link++) {

				if (newHMM.hmm.GetEndscLinkToLeft_State(hmmRightState,link)==hmmLeftState) {
					// this is the link we want
					foundIt=true;

					int globalEndVar=scoreVariablesInfo.localEndToVariableNumberVector[hmmRightState][link];
					if (globalEndVar!=-1) {
						int localEndVar=BindLocalVariable(globalEndVar,globalToLocalVariables,numLocalVariables);
						float cmEndsc=sourceCM.GetEndsc(cmState);
						assertr(cmEndsc!=(float)IMPOSSIBLE); // else why did we create a variable?

						InequalityTerm endTerm;
						endTerm.variableNum=localEndVar;
						Inequality baseIneq;
						baseIneq.inequalityType=IneqType_GE;
						baseIneq.weight=1.0;
						baseIneq.nucEmitCount[0]=baseIneq.nucEmitCount[1]=baseIneq.nucEmitCount[2]=baseIneq.nucEmitCount[3]=0;
						baseIneq.sumOfConstantsInHmm=0;

						baseIneq.lhs.push_back(endTerm);
						baseIneq.rhs=cmEndsc;

						switch (sourceCM.GetNumSymbolsEmitted(cmState)) {
							case 0:
								// although I wouldn't expect this to happen for delete states, it can happen for start states that are children of a bifurcation
								{
									Inequality ineq=baseIneq;
									assert(!ineq.lhs.empty());
									inequalityList.push_back(ineq);
								}
								break;
							case 1:
								for (int nuc=0; nuc<Alphabet_size; nuc++) {
									InfernalHmm::State hmmState;
									if (sourceCM.EmitsLeft(cmState)) {
										hmmState=hmmLeftState;
									}
									else {
										hmmState=hmmRightState;
									}
									int globalEmitVar=scoreVariablesInfo.emissionToVariableNumVector[hmmState][nuc];
									assertr(globalEmitVar!=-1); // we shouldn't be in an insert state, so this should be a variable
									int localEmitVar=BindLocalVariable(globalEmitVar,globalToLocalVariables,numLocalVariables);
									InequalityTerm localEmitTerm;
									localEmitTerm.variableNum=localEmitVar;

									Inequality ineq=baseIneq;
									ineq.lhs.push_back(localEmitTerm);
									ineq.rhs += sourceCM.GetSingletEmissionScore(cmState,nuc);
									assert(!ineq.lhs.empty());
									inequalityList.push_back(ineq);
								}
								break;
							case 2:
								for (int leftNuc=0; leftNuc<Alphabet_size; leftNuc++) {
									for (int rightNuc=0; rightNuc<Alphabet_size; rightNuc++) {
										int globalLeftEmitVar=scoreVariablesInfo.emissionToVariableNumVector[hmmLeftState][leftNuc];
										assertr(globalLeftEmitVar!=-1); // we shouldn't be in an insert state, so this should be a variable
										int globalRightEmitVar=scoreVariablesInfo.emissionToVariableNumVector[hmmRightState][rightNuc];
										assertr(globalRightEmitVar!=-1); // we shouldn't be in an insert state, so this should be a variable

										int localLeftEmitVar=BindLocalVariable(globalLeftEmitVar,globalToLocalVariables,numLocalVariables);
										int localRightEmitVar=BindLocalVariable(globalRightEmitVar,globalToLocalVariables,numLocalVariables);
										InequalityTerm localLeftEmitTerm,localRightEmitTerm;
										localLeftEmitTerm.variableNum=localLeftEmitVar;
										localRightEmitTerm.variableNum=localRightEmitVar;

										Inequality ineq=baseIneq;
										ineq.lhs.push_back(localLeftEmitTerm);
										ineq.lhs.push_back(localRightEmitTerm);
										ineq.rhs += sourceCM.GetPairEmissionScore(cmState,leftNuc,rightNuc);
										assert(!ineq.lhs.empty());
										inequalityList.push_back(ineq);
									}
								}
								break;
						}
					}
				}
			}
			assertr(foundIt);
		}

		cmState++;
	}
}

// recursively make the inequalities for each node along a path of consecutive nodes.
void Cm2Hmm_MakeInequalitiesForPath(HmmAndBuildInfo& newHMM,const CovarianceModel& sourceCM,const ScoreVariablesInfo& scoreVariablesInfo,const CovarianceModel::Node cmStartPathNode,const CovarianceModel::Node cmEndPathNode,const CovarianceModel::State cmRootState,const CovarianceModel::State hmmLeftState,const CovarianceModel::State hmmRightState,const bool doneEmission,InequalityList& inequalityList,vector<int>& globalToLocalVariables,int& numLocalVariables,Inequality& inequalitySoFar,WeightedInequalitiesInfo *weightedInequalitiesInfo)
{
	//fprintf(stderr,"MakeInequalitiesForPath: root state=%d,done emit=%c\n",CovarianceModel::StateToInt(cmRootState),doneEmission?'T':'F');

	// check for base case
	if (sourceCM.GetNode(cmRootState)==cmEndPathNode) {

		// it's the base case, end of path
		inequalitySoFar.pathEndState=cmRootState;

		/* // nope, not for now at least
		// add to totalWeight (for later normalization)
		totalWeight += inequalitySoFar.weight;
		*/

		if (weightedInequalitiesInfo!=NULL) {
			// do the left side of the HMM in-edges at the end of the walk (search for comment with "JHKOOIUYOIUY" for details on why)
			CovarianceModel::State hmmLeftFirstNormalState,hmmLeftLastNormalState;
			DumblyWorkOutHmmStates(hmmLeftFirstNormalState,hmmLeftLastNormalState,newHMM,cmEndPathNode,sourceCM,&Cm2HmmState::hmmLeftState);
			double leftProb=weightedInequalitiesInfo->transitionCounter->GetEntryProbability_Unreversed(hmmLeftState,hmmLeftFirstNormalState,hmmLeftLastNormalState);
			inequalitySoFar.weight *= leftProb;
		}

		// add this equation to the list
		if (!inequalitySoFar.lhs.empty()) { // can get empty inequalities for Start states in right bifurcation, which have only 1 transition and no emits
			inequalityList.push_back(inequalitySoFar);
		}
	}
	else {

		// not base case

		// take another step
		if (sourceCM.IsEmitting(cmRootState) && !doneEmission) {

			if (sourceCM.GetNumSymbolsEmitted(cmRootState)==1) {

				// case 1: only 1 symbol is emitted

				// NOTE: we could be a position where cmRootState emits 1 symbol, but both hmmLeftState AND hmmRightState are emitting.  This happens, e.g., in ROOT node in a transition from the IL state to the IR state, where hmmLeftState is still on the IL state.
				CovarianceModel::State hmmEmitState;
				if (sourceCM.EmitsLeft(cmRootState)) {
					assert(newHMM.hmm.IsEmitting(hmmLeftState));
					hmmEmitState=hmmLeftState;
				}
				else {
					assert(sourceCM.EmitsRight(cmRootState));
					assert(newHMM.hmm.IsEmitting(hmmRightState));
					hmmEmitState=hmmRightState;
				}

				// special case: given the current implementation of infernal-0.54, and the fact that we don't need to solve much for IL/IR nodes, we'll have many IL/IR nodes where the emissions don't have variables assigned for any symbol.  In this case, we only need to do the recursion for one of the vars.  This saves time, and may allow us to build more sophisticated sets of equations
				if (sourceCM.IsInsertState(cmRootState)) {
					// is IL/IR: activate special case, but verify it's true

					// sanity checking code
					for (int nuc=0; nuc<Alphabet_size; nuc++) {
#ifdef _DEBUG
					  float emitScore=sourceCM.GetSingletEmissionScore(cmRootState,nuc);
						assert (scoreVariablesInfo.emissionToVariableNumVector[hmmEmitState][nuc]==-1); // IL/IR states should never need vars for emissions
						assert(newHMM.hmm.GetSingletEmissionScore(hmmEmitState,nuc)==emitScore); // I should have set it
#endif
					}

					// NOTE: don't add anything to inequalitySoFar.sumOfConstantsInHmm or inequalitySoFar.nucEmitCount, since this is an insert state, and insert states (since their emits are not vars, and are equal) get treated weirdly in other code

					// we got here -- it's safe, so recurse
					Cm2Hmm_MakeInequalitiesForPath(newHMM,sourceCM,scoreVariablesInfo,cmStartPathNode,cmEndPathNode,cmRootState,hmmLeftState,hmmRightState,true,inequalityList,globalToLocalVariables,numLocalVariables,inequalitySoFar,weightedInequalitiesInfo);
				}
				else {

					// normal case
					// do the emission
					for (int nuc=0; nuc<Alphabet_size; nuc++) {

						TemporarilyModifyInequality temporarilyModifyInequality(inequalitySoFar);

						// do weighting.  Technically, weighting shouldn't care about whether it's an insert or match (i.e. whether or not we need a variable for the emission cost).  However, we can put this code here: for the insert case, we're really exploring all nucs at once (since we don't use a variable), and the total probability for the possible emissions is 1
						if (weightedInequalitiesInfo!=NULL) {
							temporarilyModifyInequality.MultiplyWeight(weightedInequalitiesInfo->transitionCounter->GetEmitFrequency(hmmEmitState,nuc));
						}

						inequalitySoFar.nucEmitCount[nuc]++;

						float emitScore=sourceCM.GetSingletEmissionScore(cmRootState,nuc);
						if (scoreVariablesInfo.emissionToVariableNumVector[hmmEmitState][nuc]==-1) {
							assert(newHMM.hmm.GetSingletEmissionScore(hmmEmitState,nuc)==emitScore); // I should have set it this way, if no variable was assigned, in which case, they cancel out
							assert(false); // now that I'm doing inequalitySoFar.sumOfConstantsInHmm, all emits except insert states should have variables.  Otherwise, the code has to change to accomodate this.
						}
						else {
							// put it in
							int globalVariableNum=scoreVariablesInfo.emissionToVariableNumVector[hmmEmitState][nuc];
							temporarilyModifyInequality.AddVariable (globalVariableNum,globalToLocalVariables,numLocalVariables);
							temporarilyModifyInequality.AddScore(emitScore);
						}

						// and recurse
						Cm2Hmm_MakeInequalitiesForPath(newHMM,sourceCM,scoreVariablesInfo,cmStartPathNode,cmEndPathNode,cmRootState,hmmLeftState,hmmRightState,true,inequalityList,globalToLocalVariables,numLocalVariables,inequalitySoFar,weightedInequalitiesInfo);
					}
				}
			}
			else {
				assert(sourceCM.GetNumSymbolsEmitted(cmRootState)==2); // only remaining possibility

				for (int leftNuc=0; leftNuc<Alphabet_size; leftNuc++) {
					for (int rightNuc=0; rightNuc<Alphabet_size; rightNuc++) {

						TemporarilyModifyInequality temporarilyModifyInequality(inequalitySoFar);
						float emitScore=sourceCM.GetPairEmissionScore(cmRootState,leftNuc,rightNuc);

						inequalitySoFar.nucEmitCount[leftNuc]++;
						inequalitySoFar.nucEmitCount[rightNuc]++;

						// weight by join prob
						if (weightedInequalitiesInfo!=NULL) {
							temporarilyModifyInequality.MultiplyWeight(weightedInequalitiesInfo->transitionCounter->GetEmitFrequency(hmmLeftState,leftNuc));
							temporarilyModifyInequality.MultiplyWeight(weightedInequalitiesInfo->transitionCounter->GetEmitFrequency(hmmRightState,rightNuc));
						}

						int globalLeftVariableNum=scoreVariablesInfo.emissionToVariableNumVector[hmmLeftState][leftNuc];
						int globalRightVariableNum=scoreVariablesInfo.emissionToVariableNumVector[hmmRightState][rightNuc];
						assert(globalLeftVariableNum!=-1 && globalRightVariableNum!=-1); // MP_st should map to two ML_st in the HMM, and since they're not insert states, we should have variables for both of their emissions

						temporarilyModifyInequality.AddVariable (globalLeftVariableNum,globalToLocalVariables,numLocalVariables);
						temporarilyModifyInequality.AddVariable (globalRightVariableNum,globalToLocalVariables,numLocalVariables);
						temporarilyModifyInequality.AddScore (emitScore);

						// and recurse
						Cm2Hmm_MakeInequalitiesForPath(newHMM,sourceCM,scoreVariablesInfo,cmStartPathNode,cmEndPathNode,cmRootState,hmmLeftState,hmmRightState,true,inequalityList,globalToLocalVariables,numLocalVariables,inequalitySoFar,weightedInequalitiesInfo);
					}
				}
			}
		}
		else {

			// take a transition to another state
			for (int childNum=0; childNum<sourceCM.GetNumChildren(cmRootState); childNum++) {
				CovarianceModel::State cmNextState=sourceCM.GetNthChildState(cmRootState,childNum);
				if (cmNextState!=cmRootState) { // no need to try self loops

					TemporarilyModifyInequality temporarilyModifyInequality(inequalitySoFar);

					float cmTsc=sourceCM.GetNthChildTsc(cmRootState,childNum);
					temporarilyModifyInequality.AddScore(cmTsc);

					CovarianceModel::State hmmNextLeftState=newHMM.cm2HmmState[cmNextState].hmmLeftState;
					CovarianceModel::State hmmNextRightState=newHMM.cm2HmmState[cmNextState].hmmRightState;

					Cm2Hmm_MakeInequalitiesForPath_HmmTransition(temporarilyModifyInequality,hmmLeftState,hmmNextLeftState,newHMM,scoreVariablesInfo,globalToLocalVariables,numLocalVariables,weightedInequalitiesInfo);
					Cm2Hmm_MakeInequalitiesForPath_HmmTransition(temporarilyModifyInequality,hmmRightState,hmmNextRightState,newHMM,scoreVariablesInfo,globalToLocalVariables,numLocalVariables,weightedInequalitiesInfo);

					// and recurse
					Cm2Hmm_MakeInequalitiesForPath(newHMM,sourceCM,scoreVariablesInfo,cmStartPathNode,cmEndPathNode,cmNextState,hmmNextLeftState,hmmNextRightState,false,inequalityList,globalToLocalVariables,numLocalVariables,inequalitySoFar,weightedInequalitiesInfo);
				}
			}
		}
	}
}


void SetHmmScores(HmmAndBuildInfo& newHMM,const ScoreVariablesInfo& scoreVariablesInfo,const vector<int>& globalToLocalVariables,const vector<float>& localVariablesToValue)
{
	CovarianceModel::State hmmState;
	for (hmmState=newHMM.hmm.GetFirstState(); hmmState!=newHMM.hmm.GetLastState(); hmmState++) {
		for (unsigned int i=0;	i<scoreVariablesInfo.transitionToVariableNumVector[hmmState].size(); i++) {
			int globalVar=scoreVariablesInfo.transitionToVariableNumVector[hmmState][i];
			if (globalVar!=-1) {
				int localVar=globalToLocalVariables[globalVar];
				if (localVar!=-1) {
					int childNum=(int)i;
					newHMM.hmm.SetTransitionLogScore(hmmState,childNum,localVariablesToValue[localVar]);
				}
			}
		}
		for (unsigned int i=0;	i<scoreVariablesInfo.emissionToVariableNumVector[hmmState].size(); i++) {
			int globalVar=scoreVariablesInfo.emissionToVariableNumVector[hmmState][i];
			if (globalVar!=-1) {
				int localVar=globalToLocalVariables[globalVar];
				if (localVar!=-1) {
					int nuc=i;
					newHMM.hmm.SetSingletEmissionLogScore(hmmState,nuc,localVariablesToValue[localVar]);
				}
			}
		}

		// local stuff
		if (newHMM.hmm.DoLocal()) {
			int globalVar;

			globalVar=scoreVariablesInfo.leftwardBeginToVariableNumberVector[hmmState];
			if (globalVar!=-1) {
				int localVar=globalToLocalVariables[globalVar];
				if (localVar!=-1) {
					newHMM.hmm.SetLeftwardBeginsc(hmmState,localVariablesToValue[localVar]);
				}
			}

			globalVar=scoreVariablesInfo.rightwardBeginToVariableNumberVector[hmmState];
			if (globalVar!=-1) {
				int localVar=globalToLocalVariables[globalVar];
				if (localVar!=-1) {
					newHMM.hmm.SetRightwardBeginsc(hmmState,localVariablesToValue[localVar]);
				}
			}

			for (int link=0; link<(int)(scoreVariablesInfo.localEndToVariableNumberVector[hmmState].size()); link++) {

				globalVar=scoreVariablesInfo.localEndToVariableNumberVector[hmmState][link];
				if (globalVar!=-1) {
					int localVar=globalToLocalVariables[globalVar];
					if (localVar!=-1) {
						newHMM.hmm.SetEndscLinkToLeft_Endsc(hmmState,link,localVariablesToValue[localVar]);
					}
				}
			}
		}
	}
}

// warning: I think this function takes some short cuts that only work when walking delete state (well, obviously it doesn't do any emissions)
void WalkHmmDeleteStates(const HmmAndBuildInfo& newHMM,const CovarianceModel& cm)
{
	float cmCost=0,hmmCost=0;
	CovarianceModel::State cmState,hmmLeftState,hmmRightState;
	cmState=cm.GetFirstState();
	hmmLeftState=newHMM.cm2HmmState[cmState].hmmLeftState;
	hmmRightState=newHMM.cm2HmmState[cmState].hmmRightState;
	while (cm.GetStateType(cmState)!=E_st && cm.GetStateType(cmState)!=B_st) {

		fprintf(dumpFile,"CM=(%d) $%f   /   HMM=(%d,%d) $%f\n",CovarianceModel::StateToInt(cmState),cmCost,CovarianceModel::StateToInt(hmmLeftState),CovarianceModel::StateToInt(hmmRightState),hmmCost);

		CovarianceModel::State cmOkayChildState=CovarianceModel::GetInvalidState();
		float tscCmCost=0; // make compiler realize it is initialized
		for (int childNum=0; childNum<cm.GetNumChildren(cmState); childNum++) {
			CovarianceModel::State cmChildState=cm.GetNthChildState(cmState,childNum);
			int cmStateType=cm.GetStateType(cmChildState);
			if (cmStateType==D_st || cmStateType==PASSTHRU_st || cmStateType==E_st || cmStateType==B_st) {
				assert(cmOkayChildState==CovarianceModel::GetInvalidState()); // I think there should be only one of these
				cmOkayChildState=cmChildState;
				tscCmCost=cm.GetNthChildTsc(cmState,childNum);
			}
		}
		assert(cmOkayChildState!=CovarianceModel::GetInvalidState()); // should have got one

		cmState=cmOkayChildState;
		cmCost += tscCmCost;

		CovarianceModel::State hmmNextLeftState=newHMM.cm2HmmState[cmState].hmmLeftState;
		CovarianceModel::State hmmNextRightState=newHMM.cm2HmmState[cmState].hmmRightState;
		if (hmmNextLeftState==CovarianceModel::GetInvalidState()) {
			hmmNextLeftState=hmmLeftState;
		}
		if (hmmNextRightState==CovarianceModel::GetInvalidState()) {
			hmmNextRightState=hmmRightState;
		}
		if (hmmNextLeftState!=hmmLeftState) {
			CovarianceModel::State fromState=hmmLeftState, toState=hmmNextLeftState;
			if (fromState>=toState) {
				std::swap(fromState,toState);
			}
			int childNum=newHMM.hmm.GetChildNum_Slow(fromState,toState);
			hmmCost += newHMM.hmm.GetNthChildTsc(fromState,childNum);
		}
		if (hmmNextRightState!=hmmRightState) {
			CovarianceModel::State fromState=hmmRightState, toState=hmmNextRightState;
			if (fromState>=toState) {
				std::swap(fromState,toState);
			}
			int childNum=newHMM.hmm.GetChildNum_Slow(fromState,toState);
			hmmCost += newHMM.hmm.GetNthChildTsc(fromState,childNum);
		}
		hmmLeftState=hmmNextLeftState;
		hmmRightState=hmmNextRightState;
	}
}

/* // this isn't the way I want to do this
void NormalizeInequalityWeights(InequalityList& inequalityList,WeightedInequalitiesInfo *weightedInequalitiesInfo,double totalWeight)
{
	double pseudocountTotal=weightedInequalitiesInfo->pseudocount * (double)(inequalityList.size());

	InequalityList::iterator ineqIter;
	for (ineqIter=inequalityList.begin(); ineqIter!=inequalityList.end(); ineqIter++) {
		ineqIter->weight = (ineqIter->weight+weightedInequalitiesInfo->pseudocount)/(totalWeight + pseudocountTotal);
	}
}
*/

vector<int> ReverseGlobalToLocalVariables(const vector<int>& globalToLocalVariables,int numLocalVariables)
{
	vector<int> localToGlobalVariables;
	localToGlobalVariables.resize(numLocalVariables);

	int globalVar;
	for (globalVar=0; globalVar<(int)(globalToLocalVariables.size()); globalVar++) {
		int localVar=globalToLocalVariables[globalVar];
		if (localVar!=-1) {
			localToGlobalVariables[localVar]=globalVar;
		}
	}

	return localToGlobalVariables;
}

void Cm2Hmm_SolveScoresForPath(HmmAndBuildInfo& newHMM,const CovarianceModel& sourceCM,const ScoreVariablesInfo& scoreVariablesInfo,CovarianceModel::Node cmStartPathNode,CovarianceModel::Node cmEndPathNode,CovarianceModel::Node firstNode,CovarianceModel::Node lastNode,float& get_avgInflation,const int nodesToSpanWhileSolvingScores,HmmCommittee::CommitteeBuildInfo *committeeBuildInfo,WeightedInequalitiesInfo *weightedInequalitiesInfo,ExtraCm2HmmInfo *extraCm2HmmInfo)
{
	// start in all normal nodes in cmStartPathNode, and end in any normal node in cmEndPathNode
	assert(cmStartPathNode!=cmEndPathNode); // else path is 0 len
	CovarianceModel::State cmFirstStateOfStartPathNode=sourceCM.GetFirstStateOfNode(cmStartPathNode);

	// first initialize mapping of overall variables, to variables we'll use in this Linear Program
	vector<int> globalToLocalVariables;
	globalToLocalVariables.assign(scoreVariablesInfo.numVariables,-1); // they all don't have a local mapping
	int numLocalVariables=0;

	// we'll build up a set of inequalities
	InequalityList inequalityList;

	//double totalWeight=0;
	CmNodeTypeData startNodeData=GetCmNodeTypeData(sourceCM.GetNodeType(cmStartPathNode));
	CovarianceModel::State cmRootState; // the root of this path thru the CM - try all normal nodes
	CovarianceModel::State cmFirstNormalState=cmFirstStateOfStartPathNode;
	CovarianceModel::State cmLastNormalState=cmFirstStateOfStartPathNode.PlusInt(startNodeData.numNormalStates);
	//CovarianceModel::State cmActualLastNormalState=cmLastNormalState.PlusInt(-1);

	// NOTE (JHKOOIUYOIUY): the left side of the HMM goes up in state#, the right side goes down.  To handle this, we check the in-edges for the left side at cmEndPathNode (at the end of walking a path), and the in-edges for the right side at the beginning using cmStartPathNode (i.e. in this function).  Things may also be backwards because TransitionCounter (ScanHMM.h) uses the reversed HMM used for scanning, i.e. HmmType1 (ScanHMM.h), rather than the one that we build in this file using CovarianceModel.
	CovarianceModel::State hmmRightFirstNormalState,hmmRightLastNormalState;
	DumblyWorkOutHmmStates(hmmRightFirstNormalState,hmmRightLastNormalState,newHMM,cmStartPathNode,sourceCM,&Cm2HmmState::hmmRightState);

	for (cmRootState=sourceCM.GetFirstStateOfNode(cmStartPathNode); sourceCM.IsStateInRange(cmRootState,cmFirstNormalState,cmLastNormalState); cmRootState++) {

		// find the HMM states
		CovarianceModel::State hmmLeftState=newHMM.cm2HmmState[cmRootState].hmmLeftState;
		CovarianceModel::State hmmRightState=newHMM.cm2HmmState[cmRootState].hmmRightState;
		assert(hmmLeftState!=CovarianceModel::GetInvalidState() && hmmRightState!=CovarianceModel::GetInvalidState()); // my current code requires that we always get a state, at least for normal states (as opposed to post insert states)

		// prepare root state for search
		Inequality inequalitySoFar;
		inequalitySoFar.lhs.clear();
		inequalitySoFar.rhs=0.0;
		inequalitySoFar.weight=1.0;
		inequalitySoFar.pathStartState=cmRootState;
		inequalitySoFar.sumOfConstantsInHmm=0;
		for (int i=0; i<MAXABET; i++) {
			inequalitySoFar.nucEmitCount[i]=0;
		}
		if (weightedInequalitiesInfo==NULL) {
			// no need to do anything with the weight
		}
		else {
			// start off with in-edge probabilities, on right side only (search for "JHKOOIUYOIUY" for explanation)
			double rightProb=weightedInequalitiesInfo->transitionCounter->GetEntryProbability_Unreversed(hmmRightState,hmmRightFirstNormalState,hmmRightLastNormalState);
			inequalitySoFar.weight *= rightProb;
		}

		// and start exploring paths
		Cm2Hmm_MakeInequalitiesForPath(newHMM,sourceCM,scoreVariablesInfo,cmStartPathNode,cmEndPathNode,cmRootState,hmmLeftState,hmmRightState,false,inequalityList,globalToLocalVariables,numLocalVariables,inequalitySoFar,weightedInequalitiesInfo);
		assert(inequalitySoFar.lhs.empty());

	}
	if (newHMM.hmm.DoLocal()) {
		Cm2Hmm_MakeInequalitiesForLocals(newHMM,sourceCM,scoreVariablesInfo,cmStartPathNode,cmEndPathNode,inequalityList,globalToLocalVariables,numLocalVariables);
	}

	if (weightedInequalitiesInfo!=NULL) {
		// adjust weights by pseudocounts
		double numSamples=weightedInequalitiesInfo->transitionCounter->GetNumSamples();
		double pseudocount=weightedInequalitiesInfo->pseudocount;
		for (InequalityList::iterator ineqIter=inequalityList.begin(); ineqIter!=inequalityList.end(); ineqIter++) {
			assert(ineqIter->weight>=0 && ineqIter->weight<=1);  // should be a probability of the traversal, and not have any pseudocount-type adjustment
			ineqIter->weight += pseudocount/numSamples;
		}
	}

	if (extraCm2HmmInfo!=NULL) {
		extraCm2HmmInfo->inequalitiesAndLocalVariables[cmStartPathNode].numLocalVariables=numLocalVariables;
		extraCm2HmmInfo->inequalitiesAndLocalVariables[cmStartPathNode].inequalityList=inequalityList;
		extraCm2HmmInfo->inequalitiesAndLocalVariables[cmStartPathNode].globalToLocalVariables=globalToLocalVariables;
		extraCm2HmmInfo->inequalitiesAndLocalVariables[cmStartPathNode].localToGlobalVariables=ReverseGlobalToLocalVariables(globalToLocalVariables,numLocalVariables);
	}

	vector<float> localVariableToValue;
	bool actuallySolveScores=false;
	if (extraCm2HmmInfo!=NULL) {
		actuallySolveScores=extraCm2HmmInfo->actuallySolveScores;
	}

#ifdef DEBUG_DUMP
	if (actuallySolveScores) {
		if (IsFirstHmmBuilt(committeeBuildInfo)) {
			DumpInequalities(dumpFile,newHMM,cmStartPathNode,cmEndPathNode,inequalityList,globalToLocalVariables,numLocalVariables,scoreVariablesInfo,localVariableToValue,committeeBuildInfo);
			fprintf(dumpFile,"Avg inflation here = %f\n",get_avgInflation);
		}
	}
#endif

	if (actuallySolveScores) {
		SetHmmScores(newHMM,scoreVariablesInfo,globalToLocalVariables,localVariableToValue);
	}
}

void Cm2Hmm_FindScores (HmmAndBuildInfo& newHMM,const CovarianceModel& sourceCM,const ScoreVariablesInfo& scoreVariablesInfo,CovarianceModel::Node firstNode,CovarianceModel::Node lastNode,const int nodesToSpanWhileSolvingScores,HmmCommittee::CommitteeBuildInfo *committeeBuildInfo,WeightedInequalitiesInfo *weightedInequalitiesInfo,ExtraCm2HmmInfo *extraCm2HmmInfo)
{
	float totalAvgInflation=0;

	CovarianceModel::Node cmStartPathNode,cmEndPathNode;
	cmStartPathNode=firstNode;
	bool keepGoing=true;
	while (keepGoing && cmStartPathNode!=lastNode) {

		//fprintf(stderr,"cmStartPathNode=%d\n",CovarianceModel::NodeToInt(cmStartPathNode));

		cmEndPathNode=cmStartPathNode;
		for (int i=0; i<nodesToSpanWhileSolvingScores; i++) {
			if (cmEndPathNode==lastNode) {
				keepGoing=false; // we're at the end, just do this last segment, and don't re-loop
				break;
			}
			cmEndPathNode++;
		}

		float avgInflation;
		Cm2Hmm_SolveScoresForPath(newHMM,sourceCM,scoreVariablesInfo,cmStartPathNode,cmEndPathNode,firstNode,lastNode,avgInflation,nodesToSpanWhileSolvingScores,committeeBuildInfo,weightedInequalitiesInfo,extraCm2HmmInfo);
		totalAvgInflation += avgInflation;

		cmStartPathNode=cmEndPathNode;
	}

#ifdef DEBUG_DUMP
	fprintf(dumpFile,"total avg inflation for this sub-CM: %f\n",totalAvgInflation);
#endif

	//WalkHmmDeleteStates(newHMM,sourceCM);
}

/*
Adjusts Cm2Hmm mapping data in preparation for merging profile HMMs.
In particular, this implements the case where some states in the HMM must be given higher numbers
in order to make room for states that will be inserted.
cm: CM we're making the profile HMM for
cm2HmmState: Cm2Hmm mapping data
firstHmmStateToInc: the lowest-numbered HMM state that needs to be incremented to a higher number, i.e. the splice point.  Everything higher than this is increased.  If firstHmmStateToInc is the invalid state, then everything's incremented
incrementBy: how many states are added
*/
void IncCm2HmmState(const CovarianceModel& cm,Cm2HmmStateVector& cm2HmmState,InfernalHmm::State firstHmmStateToInc,int incrementBy)
{
	CovarianceModel::State cmState;
	for (cmState=cm.GetFirstState(); cmState!=cm.GetLastState(); cmState++) {
		if (cm2HmmState[cmState].hmmLeftState!=InfernalHmm::GetInvalidState()) {
			if (cm2HmmState[cmState].hmmLeftState>firstHmmStateToInc || firstHmmStateToInc==InfernalHmm::GetInvalidState()) { // sic -- left side only inc's if it's _greater_.  This is for the case of overlapping HMM states where firstHmmStateToInc will correspond to E_st, which is both left & right
				cm2HmmState[cmState].hmmLeftState += incrementBy;
			}
		}
		if (cm2HmmState[cmState].hmmRightState!=InfernalHmm::GetInvalidState()) {
			if (cm2HmmState[cmState].hmmRightState>=firstHmmStateToInc || firstHmmStateToInc==InfernalHmm::GetInvalidState()) {
				cm2HmmState[cmState].hmmRightState += incrementBy;
			}
		}
	}
}
// moves left state if >=
void IncCm2HmmState_Literal (const CovarianceModel& cm,Cm2HmmStateVector& cm2HmmState,InfernalHmm::State firstHmmStateToInc,int incrementBy)
{
	CovarianceModel::State cmState;
	for (cmState=cm.GetFirstState(); cmState!=cm.GetLastState(); cmState++) {
		if (cm2HmmState[cmState].hmmLeftState!=InfernalHmm::GetInvalidState()) {
			if (cm2HmmState[cmState].hmmLeftState>=firstHmmStateToInc || firstHmmStateToInc==InfernalHmm::GetInvalidState()) {
				cm2HmmState[cmState].hmmLeftState += incrementBy;
			}
		}
		if (cm2HmmState[cmState].hmmRightState!=InfernalHmm::GetInvalidState()) {
			if (cm2HmmState[cmState].hmmRightState>=firstHmmStateToInc || firstHmmStateToInc==InfernalHmm::GetInvalidState()) {
				cm2HmmState[cmState].hmmRightState += incrementBy;
			}
		}
	}
}

/*
Internal function, see other overload of MergeCm2HmmState
*/
void MergeCm2HmmState(Cm2HmmStateVector& cm2HmmState,const Cm2HmmStateVector& cm2HmmStateToAdd,CovarianceModel::State state,CovarianceModel::State Cm2HmmState::*hmmState)
{
	bool has=cm2HmmState[state].*hmmState!=CovarianceModel::GetInvalidState();
	bool hasToAdd=cm2HmmStateToAdd[state].*hmmState!=CovarianceModel::GetInvalidState();
	if (has && hasToAdd) {
		assert(false); // they both do the same state??
		throw SimpleStringException("Internal error %s:%d",__FILE__,__LINE__);
	}
	if (!has && hasToAdd) {
		// get it from ToAdd
		cm2HmmState[state].*hmmState=cm2HmmStateToAdd[state].*hmmState;
	}
}

/*
For merging sub-HMMs, must merge the Cm2Hmm mapping data.

cm: source CM
cm2HmmState: one part of the Cm2Hmm mapping data, and also the output
cm2HmmStateToAdd: the other part of the Cm2Hmm mapping data
*/
void MergeCm2HmmState(const CovarianceModel& cm,Cm2HmmStateVector& cm2HmmState,const Cm2HmmStateVector& cm2HmmStateToAdd)
{
	assert(cm2HmmState.size()==cm2HmmStateToAdd.size());

	CovarianceModel::State state;
	for (state=cm.GetFirstState(); state!=cm.GetLastState(); state++) {
		MergeCm2HmmState(cm2HmmState,cm2HmmStateToAdd,state,&Cm2HmmState::hmmLeftState);
		MergeCm2HmmState(cm2HmmState,cm2HmmStateToAdd,state,&Cm2HmmState::hmmRightState);
	}
}

void SetStatePassthruToNextState(InfernalHmm& hmm,InfernalHmm::State state)
{
#ifdef _DEBUG
	int stateType=hmm.GetStateType(state);
	assert(stateType==E_st || stateType==PASSTHRU_st); // Other states could work, but it'd be suspicious -- worth checking
#endif

	hmm.SetNumChildren(state,1);
	hmm.SetFirstChild(state,state.PlusInt(1));
	hmm.SetTransitionLogScore(state,0,0);
}

/*
Purpose: build an HMM structurally-only based on a CM, or a part of a CM.  The function calls itself recursively to build
parts of the CM to deal with bifurication nodes.

Params:
firstNode - the CM node of type BEGL, BEGR or ROOT that will correspond to the start state of the HMM
justBuildToBIFNode - if justBuildToBIFNode==false, and we get to a BIF node, the function should recurse and build each component, and
stitch them together.  if justBuildToBIFNode==true, and we get to a BIF node, the we're already recursing, so just build up to the BIF node as if it were an END node.  Regarless of the value of justBuildToBIFNode, if we get to an END node, that's just the end.
*/
void Cm2Hmm_Structurally (HmmAndBuildInfo& newHMM,bool overlapMergedSubHmms,const CovarianceModel& sourceCM,const CovarianceModel::Node firstNode,ExtraCm2HmmInfo *extraCm2HmmInfo,Cm2Hmm_HmmBuildType hmmBuildType,const bool justBuildToBIFNode=false)
{
	CovarianceModel::Node endingNode=FindFirstEndingNode(sourceCM,firstNode);

#ifdef DUMP_INTERMEDIATE_HMMS
	std::string baseFileName="debug_",block="block_",endFileName=".hmm.csv";
	char buf[256];
#endif

	int endingNodeType=sourceCM.GetNodeType(endingNode);
	assert(endingNodeType==END_nd || endingNodeType==BIF_nd);
	if (justBuildToBIFNode || endingNodeType==END_nd) {
		// base case
		Cm2Hmm_Structurally_Block(newHMM,overlapMergedSubHmms,sourceCM,firstNode,endingNode,extraCm2HmmInfo,hmmBuildType);

#ifdef DUMP_INTERMEDIATE_HMMS
		sprintf(buf,"%d",CovarianceModel::NodeToInt(firstNode));
		std::string fileName=baseFileName+block+buf+endFileName;
		newHMM.hmm.SetCm2HmmState(newHMM.cm2HmmState);
		newHMM.hmm.DumpCsv(fileName.c_str(),sourceCM);
#endif

		SolveScoresPath thisPath;
		thisPath.cmFirstNode=firstNode;
		thisPath.cmEndingNode=endingNode;
		newHMM.solveScoresPathList.push_back(thisPath);

#ifdef DEBUG_DUMP
		DumpHmmAndBuildInfo (dumpFile,newHMM,sourceCM);
#endif
	}
	else {
		// BIF case
		// Build 3 sub-HMMs: one rooted at the current firstNode, and one at each of the child states of the BIF
		// then hook these HMMs together.

		// work out where the subtrees are
		CovarianceModel::State bifState=sourceCM.GetFirstStateOfNode(endingNode);
		if (!sourceCM.IsBifurication(bifState)) {
			throw SimpleStringException("BIF nodes should have only a single B_st.  What kind of a world do we live in?");
		}
		CovarianceModel::State leftChildState=sourceCM.GetLeftBifurifactionChild(bifState);
		CovarianceModel::Node leftChildNode=sourceCM.GetNode(leftChildState);
		CovarianceModel::State rightChildState=sourceCM.GetRightBifurifactionChild(bifState);
		CovarianceModel::Node rightChildNode=sourceCM.GetNode(rightChildState);

		// and recurse
		HmmAndBuildInfo rootHmm,leftChildHmm,rightChildHmm;
		Cm2Hmm_Structurally(rootHmm,overlapMergedSubHmms,sourceCM,firstNode,extraCm2HmmInfo,hmmBuildType,true);
		assert(rootHmm.leftToRightPassthruState!=CovarianceModel::GetInvalidState());
		Cm2Hmm_Structurally(leftChildHmm,overlapMergedSubHmms,sourceCM,leftChildNode,extraCm2HmmInfo,hmmBuildType,false);
		Cm2Hmm_Structurally(rightChildHmm,overlapMergedSubHmms,sourceCM,rightChildNode,extraCm2HmmInfo,hmmBuildType,false);

#ifdef DEBUG_DUMP
		fprintf(dumpFile,"Recursed HMM construction.  root @%d,  left child @%d,  right child @%d\n",CovarianceModel::NodeToInt(firstNode),CovarianceModel::NodeToInt(leftChildNode),CovarianceModel::NodeToInt(rightChildNode));
		fprintf(dumpFile,"Starting HMM merge.\n");
		fprintf(dumpFile,"\nroot hmm: \n");
		DumpHmm(dumpFile,rootHmm,sourceCM);
		fprintf(dumpFile,"\nleft child hmm: \n");
		DumpHmm(dumpFile,leftChildHmm,sourceCM);
		fprintf(dumpFile,"\nright child hmm: \n");
		DumpHmm(dumpFile,rightChildHmm,sourceCM);
#endif

		assertr(rootHmm.leftToRightPassthruState!=CovarianceModel::GetInvalidState());

		if (overlapMergedSubHmms) {

			// original stuff, with overlapping

			// order is: left part of rootHmm, left child, right child, right part of rootHmm

			// first splice rightChildHmm onto leftChildHmm s.t. the start state of rightChildHmm takes over the end state of leftChildHmm.  (So we have to increment the states #s in rightChildHmm in order to make room for leftChildHmm's states)
			// in the cm2HmmState vectors
			IncCm2HmmState(sourceCM,rightChildHmm.cm2HmmState,CovarianceModel::GetInvalidState(),leftChildHmm.hmm.GetNumStates()-1);
			MergeCm2HmmState(sourceCM,leftChildHmm.cm2HmmState,rightChildHmm.cm2HmmState);
			// and in the HMMs
			CovarianceModel::State originalLeftChildEndState=leftChildHmm.hmm.GetActualLastState();
			assert(leftChildHmm.hmm.GetStateType(originalLeftChildEndState)==E_st); // this is what we're clobbering
			leftChildHmm.hmm.AddStates(rightChildHmm.hmm.GetNumStates()-1); // -1 for the E_st we're clobbering
			rightChildHmm.hmm.MoveStatesHigher(rightChildHmm.hmm.GetFirstState(),originalLeftChildEndState);
			leftChildHmm.hmm.CopyStatesVerbatimFrom(rightChildHmm.hmm,originalLeftChildEndState,leftChildHmm.hmm.GetLastState());

			// now splice the curr leftChildHmm (which is a concatenation of the original leftChildHmm and rightChildHmm) into the (middle-ish) part of rootHmm that corresponded to the original E_st in the CM that rootHmm was built from (i.e. rootHmm.leftToRightPassthruState).
			// move the right-side of rootHmm up to make room, and also move the left-to-right passthru state up, because it knows what its children are.
			int numLeftChildStatesThatAreMoving=leftChildHmm.hmm.GetNumStates()-1; // -1 since its E_st is going to be tragically lost in the move
			// in cm2HmmState
			IncCm2HmmState(sourceCM,rootHmm.cm2HmmState,rootHmm.leftToRightPassthruState,numLeftChildStatesThatAreMoving);
			// in the actual HMM
			rootHmm.hmm.MoveStatesHigher(rootHmm.leftToRightPassthruState,rootHmm.leftToRightPassthruState.PlusInt(numLeftChildStatesThatAreMoving)); 
			// prepare leftChildHmm to be spliced in
			// in cm2HmmState
			IncCm2HmmState(sourceCM,leftChildHmm.cm2HmmState,CovarianceModel::GetInvalidState(),CovarianceModel::StateToInt(rootHmm.leftToRightPassthruState));
			// in actual HMM
			leftChildHmm.hmm.MoveStatesHigher(leftChildHmm.hmm.GetFirstState(),rootHmm.leftToRightPassthruState);
			// and copy it in
			rootHmm.hmm.CopyStatesVerbatimFrom(leftChildHmm.hmm,rootHmm.leftToRightPassthruState,rootHmm.leftToRightPassthruState.PlusInt(numLeftChildStatesThatAreMoving));
			MergeCm2HmmState(sourceCM,rootHmm.cm2HmmState,leftChildHmm.cm2HmmState);
		}
		else {
			// alternate, sans clobbering

			IncCm2HmmState_Literal(sourceCM,rightChildHmm.cm2HmmState,CovarianceModel::GetInvalidState(),leftChildHmm.hmm.GetNumStates());
			MergeCm2HmmState(sourceCM,leftChildHmm.cm2HmmState,rightChildHmm.cm2HmmState);
			CovarianceModel::State originalLeftChildPastEnd=leftChildHmm.hmm.GetLastState();
			leftChildHmm.hmm.AddStates(rightChildHmm.hmm.GetNumStates());
			rightChildHmm.hmm.MoveStatesHigher(rightChildHmm.hmm.GetFirstState(),originalLeftChildPastEnd);
			leftChildHmm.hmm.CopyStatesVerbatimFrom(rightChildHmm.hmm,originalLeftChildPastEnd,leftChildHmm.hmm.GetLastState());
			SetStatePassthruToNextState(leftChildHmm.hmm,originalLeftChildPastEnd.PlusInt(-1)); // hook up the two parts trivially

			int numLeftChildStatesThatAreMoving=leftChildHmm.hmm.GetNumStates();
			IncCm2HmmState_Literal(sourceCM,rootHmm.cm2HmmState,rootHmm.leftToRightPassthruState.PlusInt(+1),numLeftChildStatesThatAreMoving);
			rootHmm.hmm.MoveStatesHigher(rootHmm.leftToRightPassthruState.PlusInt(+1),rootHmm.leftToRightPassthruState.PlusInt(numLeftChildStatesThatAreMoving+1)); 

			IncCm2HmmState_Literal(sourceCM,leftChildHmm.cm2HmmState,CovarianceModel::GetInvalidState(),CovarianceModel::StateToInt(rootHmm.leftToRightPassthruState)+1);
			leftChildHmm.hmm.MoveStatesHigher(leftChildHmm.hmm.GetFirstState(),rootHmm.leftToRightPassthruState.PlusInt(+1));
			rootHmm.hmm.CopyStatesVerbatimFrom(leftChildHmm.hmm,rootHmm.leftToRightPassthruState.PlusInt(+1),rootHmm.leftToRightPassthruState.PlusInt(numLeftChildStatesThatAreMoving+1));
			MergeCm2HmmState(sourceCM,rootHmm.cm2HmmState,leftChildHmm.cm2HmmState);

			// stitch these together
			SetStatePassthruToNextState(rootHmm.hmm,rootHmm.leftToRightPassthruState);
			SetStatePassthruToNextState(rootHmm.hmm,rootHmm.leftToRightPassthruState.PlusInt(numLeftChildStatesThatAreMoving));
		}

		// re-build HMM-to-CM state mapping
		ReverseMapCm2HmmState(rootHmm,sourceCM);

		// merge solveScoresPathList
		rootHmm.solveScoresPathList.insert(rootHmm.solveScoresPathList.end(),leftChildHmm.solveScoresPathList.begin(),leftChildHmm.solveScoresPathList.end());
		rootHmm.solveScoresPathList.insert(rootHmm.solveScoresPathList.end(),rightChildHmm.solveScoresPathList.begin(),rightChildHmm.solveScoresPathList.end());

		// et voila -- copy it back
		newHMM.hmm.CopyFrom(rootHmm.hmm);
		newHMM.cm2HmmState=rootHmm.cm2HmmState;
		newHMM.hmm2CmStateVector=rootHmm.hmm2CmStateVector;
		newHMM.solveScoresPathList=rootHmm.solveScoresPathList;
		newHMM.leftToRightPassthruState=CovarianceModel::GetInvalidState(); // we shouldn't need this any more, and don't have a meaningful answer anyway
#ifdef DEBUG_DUMP
		fprintf(dumpFile,"\nResulting hmm:\n");
		DumpHmm(dumpFile,newHMM,sourceCM);
#endif
#ifdef DUMP_INTERMEDIATE_HMMS
		sprintf(buf,"%d",CovarianceModel::NodeToInt(firstNode));
		std::string fileName=baseFileName+buf+endFileName;
		newHMM.hmm.SetCm2HmmState(newHMM.cm2HmmState);
		newHMM.hmm.DumpCsv(fileName.c_str(),sourceCM);
#endif
	}
}

void SetupLeftRightness(HmmAndBuildInfo& newHMM,const CovarianceModel& sourceCM)
{
	// set left/right-ness of HMM states
	InfernalHmm::State hmmState;
	for (hmmState=newHMM.hmm.GetFirstState(); hmmState!=newHMM.hmm.GetLastState(); hmmState++) {
		newHMM.hmm.SetLeftState(hmmState,false);
		newHMM.hmm.SetRightState(hmmState,false);
	}
	CovarianceModel::State cmState;
	for (cmState=sourceCM.GetFirstState(); cmState!=sourceCM.GetLastState(); cmState++) {

		InfernalHmm::State hmmLeftState=newHMM.cm2HmmState[cmState].hmmLeftState;
		if (hmmLeftState!=InfernalHmm::GetInvalidState()) {
			newHMM.hmm.SetLeftState(hmmLeftState,true);
		}

		InfernalHmm::State hmmRightState=newHMM.cm2HmmState[cmState].hmmRightState;
		if (hmmRightState!=InfernalHmm::GetInvalidState()) {
			newHMM.hmm.SetRightState(hmmRightState,true);
		}
	}
}

void Cm2Hmm_Structurally(HmmAndBuildInfo& newHMM,const CovarianceModel& sourceCM,ExtraCm2HmmInfo *extraCm2HmmInfo,Cm2Hmm_HmmBuildType hmmBuildType)
{
	CovarianceModel::Node cmFirstNodeForConversion=sourceCM.GetFirstNode();
	if (sourceCM.DoLocal() || FirstNodeLooksLocallyClobbered(sourceCM)) {
		cmFirstNodeForConversion++;  // don't do anything for the whole first node, since that gets clobbered by the ConfigLocal function in modelconfig.c.
	}

	bool overlapMergedSubHmms=!sourceCM.DoLocal();

	Cm2Hmm_Structurally(newHMM,overlapMergedSubHmms,sourceCM,cmFirstNodeForConversion,extraCm2HmmInfo,hmmBuildType);

	newHMM.hmm.AllocHmmData();

	SetupLeftRightness(newHMM,sourceCM);
}

void InitLocalData (InfernalHmm& hmm,const CovarianceModel& sourceCM,bool setupLocalBegin,bool setupLocalEnd)
{
	assert(sourceCM.DoLocal()); // otherwise there's no reason to call this function

	hmm.SetDoLocal(true);

	hmm.SetLocalEndSelfLoopScore(sourceCM.GetLocalEndSelfScore());

	InfernalHmm::State hmmState;
	for (hmmState=hmm.GetFirstState(); hmmState!=hmm.GetLastState(); hmmState++) {
		if (setupLocalBegin) {
			hmm.SetLeftwardBeginsc(hmmState,(float)IMPOSSIBLE);
			hmm.SetRightwardBeginsc(hmmState,(float)IMPOSSIBLE);
		}
		if (setupLocalEnd) {
			hmm.SetNumEndscLinksToLeft(hmmState,0);
		}
	}

	CovarianceModel::State cmState;
	for (cmState=sourceCM.GetFirstState(); cmState!=sourceCM.GetLastState(); cmState++) {

		if (sourceCM.GetNode(cmState)!=sourceCM.GetFirstNode()) { // don't do anything for the whole first node, since that gets clobbered by the ConfigLocal function in modelconfig.c. 
			InfernalHmm::State hmmLeftState=hmm.GetHmmLeftStateOfCmState(cmState);
			InfernalHmm::State hmmRightState=hmm.GetHmmRightStateOfCmState(cmState);

			if (setupLocalBegin) {
				float beginsc=sourceCM.GetBeginsc(cmState);
				if (beginsc!=(float)IMPOSSIBLE) {
					

					if (hmmLeftState==InfernalHmm::GetInvalidState() || hmmRightState==InfernalHmm::GetInvalidState()) {
						assert(false);
						throw SimpleStringException("internal error %s:%d",__FILE__,__LINE__);
					}

					float hmmOldLeftBeginsc=hmm.GetLeftwardBeginsc(hmmLeftState);
					float hmmOldRightBeginsc=hmm.GetRightwardBeginsc(hmmRightState);

					// scores must sum to original beginsc (or more), for now I'll just divide the score evenly into 2.
					float newBeginsc=beginsc/(float)2.0;

					// sometimes multiple CM states can use the same HMM left or right state; so, we must take the max of the beginsc
					hmm.SetLeftwardBeginsc(hmmLeftState,std::max(hmmOldLeftBeginsc,newBeginsc));
					hmm.SetRightwardBeginsc(hmmRightState,std::max(hmmOldRightBeginsc,newBeginsc));
				}
			}
			if (setupLocalEnd) {

				if (sourceCM.GetEndsc(cmState)!=(float)IMPOSSIBLE) {

					// link right state to this left state, for simulating endsc
					if (hmmLeftState==InfernalHmm::GetInvalidState() || hmmRightState==InfernalHmm::GetInvalidState()) {
						assert(false);
						throw SimpleStringException("For making an HMM local, it's pretty tricky if the insert states (IL,IR) in the CM have non-IMPOSSIBLE endsc scores, because in my mapping, IL,IR states don't have both a left&right hmmState, and it's not easy to see what it should be.  So, I just assumed it wouldn't happen, making an ass out of me and you, the anonymous user of this program.  Sorry.");
					}

					int numLinks=hmm.GetNumEndscLinksToLeft(hmmRightState);
					hmm.SetNumEndscLinksToLeft(hmmRightState,numLinks+1);
					hmm.SetEndscLinkToLeft_State(hmmRightState,numLinks,hmmLeftState);
				}
			}
		}
	}

	if (setupLocalEnd) {
		hmm.SetAdjustedEndscGlobally(sourceCM);
	}
}

void SetLocal(HmmAndBuildInfo& newHMM,const CovarianceModel& sourceCM)
{
	InitLocalData(newHMM.hmm,sourceCM,true,true);
}

void Cm2Hmm (HmmAndBuildInfo& newHMM,Cm2Hmm_HmmBuildType hmmBuildType,const CovarianceModel& sourceCM,const int nodesToSpanWhileSolvingScores,HmmCommittee::CommitteeBuildInfo *committeeBuildInfo=NULL,WeightedInequalitiesInfo *weightedInequalitiesInfo=NULL,ExtraCm2HmmInfo *extraCm2HmmInfo=NULL)
{
	Cm2Hmm_Structurally(newHMM,sourceCM,extraCm2HmmInfo,hmmBuildType);
	newHMM.hmm.SetHmm2CmState(newHMM.hmm2CmStateVector);
	newHMM.hmm.SetCm2HmmState(newHMM.cm2HmmState);

	newHMM.hmm.SetDoLocal(false);
	if (sourceCM.DoLocal()) {
		SetLocal(newHMM,sourceCM);
	}

	ScoreVariablesInfo scoreVariablesInfo;
	SetupTransitionAndEmissionVariables(scoreVariablesInfo,newHMM,sourceCM);
	SetupReverseMapping(scoreVariablesInfo,newHMM,sourceCM);
#ifdef DEBUG_DUMP
	DumpVariables(dumpFile,scoreVariablesInfo,newHMM);
#endif
	if (extraCm2HmmInfo!=NULL) {
		extraCm2HmmInfo->scoreVariablesInfo=scoreVariablesInfo;
		InequalitiesAndLocalVariables dummyInequalitiesAndLocalVariables;
		dummyInequalitiesAndLocalVariables.numLocalVariables=0;
		extraCm2HmmInfo->inequalitiesAndLocalVariables.assign(sourceCM.GetNumNodes(),dummyInequalitiesAndLocalVariables);
	}

	SolveScoresPathList::const_iterator pathIter;
	for (pathIter=newHMM.solveScoresPathList.begin(); pathIter!=newHMM.solveScoresPathList.end(); pathIter++) {
		Cm2Hmm_FindScores(newHMM,sourceCM,scoreVariablesInfo,pathIter->cmFirstNode,pathIter->cmEndingNode,nodesToSpanWhileSolvingScores,committeeBuildInfo,weightedInequalitiesInfo,extraCm2HmmInfo);
	}

	if (sourceCM.DoLocal()) {
		newHMM.hmm.SetAdjustedEndscGlobally(sourceCM);
	}
}

void Cm2Hmm_WithWeighting_NoCaching (InfernalHmm& hmm,Cm2Hmm_HmmBuildType hmmBuildType,const CovarianceModel& sourceCM,const char *cmFileName,WeightedInequalitiesInfo *weightedInequalitiesInfo,ExtraCm2HmmInfo *extraCm2HmmInfo,const int nodesToSpanWhileSolvingScores)
{
	// verify that the structure of the CM is what we expect, so we can just assume it for the rest of this function
	if (!ValidateCmNodeTypeData(sourceCM)) {
		throw SimpleStringException("(while converting Covariance Model to HMM) CM doesn't have the structure we expect");
	}

	HmmAndBuildInfo newHMM;

	Cm2Hmm(newHMM,hmmBuildType,sourceCM,nodesToSpanWhileSolvingScores,NULL,weightedInequalitiesInfo,extraCm2HmmInfo);
	newHMM.hmm.SetFromCmFileName(cmFileName);

#ifdef DEBUG_DUMP
	DumpHmmAndBuildInfo (dumpFile,newHMM,sourceCM);
#endif

	hmm.CopyFrom(newHMM.hmm);
}

void Cm2Hmm (InfernalHmm& hmm,Cm2Hmm_HmmBuildType hmmBuildType,const CovarianceModel& sourceCM,const char *cmFileName,const std::string& programParams,bool forceCreate)
{
	const int nodesToSpanWhileSolvingScores=1; // making it higher didn't seem to lead to any more optimal scores (using RF00032==Histone3, which is small), and it takes forever to build the HMM -- it spends most of its time entering the constraints to lp_solve, since lp_solve's sparse matrix format requires it to do a lot of copying & reallocation.  I could probably optimize this by eliminating redundant constraints myself, but it doesn't look like we can improve things this way anyway.

	if (dumpFile==NULL) {
		dumpFile=ThrowingFopen("hmm-dump.txt","wt");
	}

	InfernalHmm creatingHmm;
#ifdef ENABLE_CACHING
	bool loadHmmFromCache=false;
	std::string cacheFileName;
	if (enableHmmCaching) {
		if (overrideHmmCacheFileName!=NULL) {
			cacheFileName=overrideHmmCacheFileName;
		}
		else {
			cacheFileName=MakeCm2HmmCacheFileName(cmFileName,sourceCM.DoLocal(),nodesToSpanWhileSolvingScores,1,hmmBuildType);
		}
		fprintf(stderr,"Trying to load HMM from cached file '%s'\n",cacheFileName.c_str());
		printf("Trying to load HMM from cached file '%s'\n",cacheFileName.c_str());

		loadHmmFromCache=!forceCreate;
		if (loadHmmFromCache) {
			loadHmmFromCache=creatingHmm.LoadInBinary(cacheFileName.c_str());
		}
	}
	if (loadHmmFromCache) {
#ifdef DEBUG_DUMP
		fprintf(dumpFile,"Loaded HMM from cache.\n");
		fprintf(stderr,"Loaded HMM from cache.\n");
#ifdef DEBUG_DUMP
	creatingHmm.DumpInfernalHmm(dumpFile,sourceCM);
#endif
#endif
	}
	else {
#endif
	fprintf(stderr,"Building HMM based on CM...\n");

	Cm2Hmm_WithWeighting_NoCaching(creatingHmm,hmmBuildType,sourceCM,cmFileName,NULL,NULL,nodesToSpanWhileSolvingScores);
	std::string build("(default construction of HMM) ");
	build += programParams;
	creatingHmm.AddBuildDescription(build);
	creatingHmm.SetHmmBuildType(hmmBuildType);

	fprintf(stderr,"Built HMM.\n");

#ifdef ENABLE_CACHING

		if (enableHmmCaching) {
			creatingHmm.SaveInBinary(cacheFileName.c_str());
		}
	}
#endif

	hmm.CopyFrom(creatingHmm);

	fflush(dumpFile);
}
