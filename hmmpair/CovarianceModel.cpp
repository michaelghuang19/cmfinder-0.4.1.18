#include "hmmpair.h"

int CovarianceModelBase::nextUniqueId=0;
bool CovarianceModelBase::enableInsertHack=true; // by default, we should be using this for compatibility with Infernal, and also because as I write this, my code is dependent on it
void CovarianceModelBase::DisableInsertHack (void)
{
	enableInsertHack=false;
}
const char *CovarianceModelBase::GetCmFileName (void) const
{
	if (cmFileNameLoaded.empty()) {
		throw SimpleStringException("Internal error: GetCmFileName called, but we didn't load this CM from a file.");
	}
	return cmFileNameLoaded.c_str();
}
bool CovarianceModelBase::IsStructurallyTheSame (const CovarianceModelBase& t) const
{
	if (DoLocal()!=t.DoLocal()) {
		return false;
	}
	if (GetNumStates()!=t.GetNumStates()) {
		return false;
	}
	for (State state=GetFirstState(); state!=GetLastState(); state++) {
		if (GetNumChildren(state)!=t.GetNumChildren(state)) {
			return false;
		}
		for (int child=0; child<GetNumChildren(state); child++) {
			if (GetNthChildState(state,child)!=t.GetNthChildState(state,child)) {
				return false;
			}
		}
	}
	return true;
}
CovarianceModelBase::CovarianceModelBase (void)
{
	rsearchQueryLen=-1;
	cmMemoryOwned=true;
	cm=NULL;
	cmfp=NULL;
	isRsearchCM=false;
	isSmithWaterCM=false;

	uniqueId=nextUniqueId;
	nextUniqueId++;
}
void CovarianceModelBase::Destruct (void)
{
	if (cm!=NULL && cmMemoryOwned) {
		FreeCM(cm);
		cm=NULL;
	}
	if (cmfp!=NULL) {
		CMFileClose(cmfp);
		cmfp=NULL;
	}
	cmMemoryOwned=true;
}
CovarianceModelBase::~CovarianceModelBase ()
{
	Destruct();
}
int CovarianceModelBase::GetUniqueId (void) const
{
	return uniqueId;
}
void CovarianceModelBase::DumpCsv (const char *csvFileName) const
{
	FILE *csv=ThrowingFopen(csvFileName,"wt");
	DumpCsv(csv);
	fclose(csv);
}
void CovarianceModelBase::DumpCsv (FILE *csv) const
{
	fprintf(csv,"Local,%d\n",DoLocal()?1:0);
	fprintf(csv,"# states,%d\n",GetNumStates());
	fprintf(csv,"\n");
	fprintf(csv,"state#,state type,,local begin,local end,,");
	fprintf(csv,"consensus emit,emit A,emit C,emit G,emit U,");
	fprintf(csv,",first child,transition1,transition2...\n");
#if 0
	fprintf(csv,",,,,,,");
	for (int left=0; left<Alphabet_size; left++) {
		for (int right=0; right<Alphabet_size; right++) {
			fprintf(csv,"emit %c-%c,",nucs[left],nucs[right]);
		}
	}
#endif
	for (State state=GetFirstState(); state!=GetLastState(); state++) {
		fprintf(csv,"%d,%s,,%g,%g,,",StateToInt(state),GetStateTypeName(state),GetBeginsc(state),GetEndsc(state));
		if (IsEmittingState(state)) {
			if (GetNumSymbolsEmitted(state)==1) {
				if (state!=GetConsensusStateOfNode(GetNode((state)))) {
					fprintf(csv,","); // don't both with non-consensus emit
				}
				else {
					float bestScore=-FLT_MAX;
					std::string best;
					for (int nuc=0; nuc<Alphabet_size; nuc++) {
						float thisScore=GetSingletEmissionScore(state,nuc);
						std::string thisNuc=stringprintf("%c",nucs[nuc]);
						if (thisScore>bestScore) {
							best.clear();
						}
						if (thisScore>=bestScore) {
							if (!best.empty()) {
								best += ",";
							}
							best += thisNuc;
							bestScore=thisScore;
						}
					}
					fprintf(csv,"%s,",best.c_str());
				}
			}
			if (GetNumSymbolsEmitted(state)==2) {
				float bestScore=-FLT_MAX;
				std::string best;
				for (int left=0; left<Alphabet_size; left++) {
					for (int right=0; right<Alphabet_size; right++) {
						float thisScore=GetPairEmissionScore(state,left,right);
						std::string thisNuc=stringprintf("%c-%c",nucs[left],nucs[right]);
						if (thisScore>bestScore) {
							best.clear();
						}
						if (thisScore>=bestScore) {
							if (!best.empty()) {
								best += ",";
							}
							best += thisNuc;
							bestScore=thisScore;
						}
					}
				}
				fprintf(csv,"%s,",best.c_str());
			}
			if (EmitsLeft(state)) {
				if (EmitsRight(state)) {
					for (int right=0; right<Alphabet_size; right++) {
						fprintf(csv,"%g,",GetPairEmissionScore(state,0,right));
					}
				}
				else {
					for (int nuc=0; nuc<Alphabet_size; nuc++) {
						fprintf(csv,"%g,",GetSingletEmissionScore(state,nuc));
					}
				}
			}
			else {
				for (int nuc=0; nuc<Alphabet_size; nuc++) {
					fprintf(csv,"%g,",GetSingletEmissionScore(state,nuc));
				}
			}
		}
		else {
			fprintf(csv,",");
			for (int i=0; i<4; i++) {
				fprintf(csv,",");
			}
		}

		fprintf(csv,",");
		if (GetStateType(state)==B_st) {
			fprintf(csv,"%d",StateToInt(GetLeftBifurcationChild(state)));
		}
		else {
			if (GetNumChildren(state)==0) {
				//fprintf(csv,"");
			}
			else {
				fprintf(csv,"%d,",StateToInt(GetNthChildState(state,0)));
				for (int child=0; child<GetNumChildren(state); child++) {
					fprintf(csv,"%g,",GetNthChildTsc(state,child));
				}
			}
		}

		fprintf(csv,"\n");

		if (EmitsLeft(state) && EmitsRight(state)) {
			for (int left=1; left<Alphabet_size; left++) {
				fprintf(csv,",,,,,,left,");
				for (int right=0; right<Alphabet_size; right++) {
					fprintf(csv,"%g,",GetPairEmissionScore(state,left,right));
				}
				fprintf(csv,"\n");
			}
		}
	}
}
void CovarianceModelBase::Save (char *cmFileName)
{
	FILE *out=ThrowingFopen(cmFileName,"wt");
	CMFileWrite(out,cm,FALSE);
	fclose(out);
}
bool CovarianceModelBase::Load (char *cmFileName,bool doLocalAlignment)
{
	cmFileNameLoaded=cmFileName;

	if (strcmp(cmFileName,"dummy")==0) {
		// pretend it's a file
		return true;
	}

	// if it looks like it came from tRNAscan-SE, make sure enableInsertHack==true
	if (strlen(cmFileName)>=4) {

		const char *cursor=cmFileName+strlen(cmFileName);
		while (cursor!=cmFileName) {
			if (*cursor=='/' || *cursor=='\\') {
				cursor++;
				break;
			}
			cursor--;
		}

		if ((cursor[0]=='t' || cursor[0]=='T')
		&& (cursor[1]=='r' || cursor[1]=='R')
		&& (cursor[2]=='n' || cursor[2]=='N')
		&& (cursor[3]=='a' || cursor[3]=='A')) {
			if (enableInsertHack) {
				throw SimpleStringException("CM file name begins with 'trna', which suggests it might have come from tRNAscan-SE.  Are you sure you didn't mean to use the --no-insert-hack flag?  To get around this error (__IF__ it's a false alarm), change the file name.\n");
			}
		}
	}

	isRsearchCM=false;
	isSmithWaterCM=false;
	char RSEARCHCM[]=".rsearchcm";
	if (strlen(cmFileName)>strlen(RSEARCHCM)) {
		const char *cursor=cmFileName+strlen(cmFileName)-strlen(RSEARCHCM);
		if (strcmp(cursor,RSEARCHCM)==0) {
			isRsearchCM=true;
		}
	}

	if ((cmfp = CMFileOpen(cmFileName, NULL)) == NULL)
		throw SimpleStringException("Failed to open covariance model save file %s\n", cmFileName);

	{
		if (! CMFileRead(cmfp, &cm))
			throw SimpleStringException("Failed to read a CM from %s -- file corrupt?\n", cmFileName);
	}
	if (cm == NULL) 
		throw SimpleStringException("%s empty?\n", cmFileName);

	if (isRsearchCM) {
		CMLogoddsify(cm);
		if (doLocalAlignment) {
			cm->flags |= CM_LOCAL_BEGIN;
			cm->flags |= CM_LOCAL_END;
		}
		else {
			// I'm allowing global.  throw SimpleStringException("To use an RSEARCH CM, you must use the --local flag, since RSEARCH CMs are implicitly local.  (I could make it so there's an alternative, but I haven't.)");
		}
		// read extra data, but don't do other things to the model.  This data always appears in an RSEARCH CM, even if we're using it in global mode
		for (State state=GetFirstState(); state!=GetLastState(); state++) {
			float beginsc,endsc;
			fscanf(cmfp->f,"%g\t%g\n",&beginsc,&endsc);
			if (DoLocal()) {
				SetBeginsc(state,beginsc);
				SetBeginProbDirectly(state,(float)pow2(beginsc));
				SetEndsc(state,endsc);
				SetEndProbDirectly(state,(float)pow2(endsc));
			}
		}
		fscanf(cmfp->f,"%d\n",&rsearchQueryLen);
		if (rsearchQueryLen==-2) {
			// newer format
			int format;
			fscanf(cmfp->f,"%d\n",&format);
			fscanf(cmfp->f,"%d\n",&rsearchQueryLen);
			if (format>=2) {
				int sw;
				fscanf(cmfp->f,"%d\n",&sw);
				isSmithWaterCM=sw!=0;
			}
			std::string programParams;
			int ch;
			while (1) {
				ch=fgetc(cmfp->f);
				if (ch=='\n') {
					break;
				}
				programParams += (char)ch;
			}
			printf("RSEARCH CM created with params: %s\n",programParams.c_str());
		}
	}
	else {
		if (doLocalAlignment) {

			//throw SimpleStringException("Local alignments aren't properly handled, since local-ends turned out to be more tricky than I expected.");

			ConfigLocal(cm, 0.5, 0.5);
		}

		CMLogoddsify(cm);
		if (enableInsertHack) {
			CMHackInsertScores(cm);	/* make insert emissions score zero. "TEMPORARY" FIX. */
		}
	}

	return true;
}

namespace tRNAscanSE { // stuff related to reverse engineering of tRNAscan-SE by Todd Lowe & Sean Eddy

	// this is from tRNAscan-SE-1.23/structs.h
	enum NodeType {
		BIFURC_NODE=0,
		MATP_NODE=1,
		MATL_NODE=2,
		MATR_NODE=3,
		BEGINL_NODE=4,
		BEGINR_NODE=5,
		ROOT_NODE=6,
		END_NODE=BIFURC_NODE
	};
	enum StateType {
		DEL_ST=0,
		MATP_ST=1,
		MATL_ST=2,
		MATR_ST=3,
		INSL_ST=4,
		INSR_ST=5,
		STATETYPES=6, /* MATP nodes contain 6 states */
		BEGIN_ST=DEL_ST,
		BIFURC_ST=DEL_ST,
		END_ST=DEL_ST,
	};

	struct NodeTypeInfo {
		NodeType nodeType;
		int infernalNodeType;

		int numSplitSetStates;
		int numInsertStates;
		StateType *states;
		bool firstStateIsStart;
	};
	StateType matpStates[]={MATP_ST,MATL_ST,MATR_ST,DEL_ST,INSL_ST,INSR_ST};
	StateType matlStates[]={MATL_ST,DEL_ST,INSL_ST};
	StateType matrStates[]={MATR_ST,DEL_ST,INSR_ST};
	StateType rootStates[]={DEL_ST,INSL_ST,INSR_ST};
	StateType beglStates[]={DEL_ST};
	StateType begrStates[]={DEL_ST,INSL_ST};
	StateType bifStates[]={BIFURC_ST};
	NodeTypeInfo nodeTypeInfo[]={
		{BIFURC_NODE,BIF_nd,1,0,bifStates,false}, // BIF handled specially, though
		{MATP_NODE,MATP_nd,4,2,matpStates,false},
		{MATL_NODE,MATL_nd,2,1,matlStates,false},
		{MATR_NODE,MATR_nd,2,1,matrStates,false},
		{BEGINL_NODE,BEGL_nd,1,0,beglStates,true},
		{BEGINR_NODE,BEGR_nd,1,1,begrStates,true},
		{ROOT_NODE,ROOT_nd,1,2,rootStates,true}
	};
	NodeTypeInfo GetNodeTypeInfo (NodeType nodeType) {
		for (unsigned int i=0; i<sizeof(nodeTypeInfo)/sizeof(NodeTypeInfo); i++) {
			if (nodeTypeInfo[i].nodeType==nodeType) {
				return nodeTypeInfo[i];
			}
		}
		throw SimpleStringException("GetNodeTypeInfo failed for node type %d",(int)nodeType);
	}

	struct StateTypeInfo {
		StateType stateType;
		int infernalStateType;

		int numEmits;
	};
	StateTypeInfo stateTypeInfo[]={
		{DEL_ST,D_st,0},
		{MATP_ST,MP_st,2},
		{MATL_ST,ML_st,1},
		{MATR_ST,MR_st,1},
		{INSL_ST,IL_st,1},
		{INSR_ST,IR_st,1},
		{BIFURC_ST,B_st,0},
		{BEGIN_ST,S_st,0},
		{END_ST,E_st,0},
		{END_ST,EL_st,0}
	};
	StateTypeInfo GetStateTypeInfo (StateType stateType) {
		for (unsigned int i=0; i<sizeof(stateTypeInfo)/sizeof(StateTypeInfo); i++) {
			if (stateTypeInfo[i].stateType==stateType) {
				return stateTypeInfo[i];
			}
		}
		throw SimpleStringException("GetStateTypeInfo failed for state type %d",(int)stateType);
	}
	StateType InfernalToCoveStateType (int infernalStateType) {
		for (unsigned int i=0; i<sizeof(stateTypeInfo)/sizeof(StateTypeInfo); i++) {
			if (stateTypeInfo[i].infernalStateType==infernalStateType) {
				return stateTypeInfo[i].stateType;
			}
		}
		throw SimpleStringException("InfernalToCoveStateType failed for state type %d",infernalStateType);
	}

	typedef float qTransitionMatrix[STATETYPES][STATETYPES];
	struct TransitionMatrix {qTransitionMatrix array;};
	typedef float qSingletEmit[4];
	struct SingletEmit {qSingletEmit array;} ;
	typedef float qPairEmit[4][4];
	struct PairEmit {qPairEmit array;} ;
};
void CopySingletEmitFromCove(CM_t *cm,int state,const tRNAscanSE::SingletEmit& coveEmit)
{
	for (int nuc=0; nuc<4; nuc++) {
		cm->e[state][nuc]=coveEmit.array[nuc];
		cm->esc[state][nuc]=(float)(sreLOG2(coveEmit.array[nuc]));
	}
}
void CovarianceModelBase::LoadCove2 (const char *cmFileName)
{
	Destruct();

	cmFileNameLoaded=cmFileName;

	bool useLocalEndStates=false; // I wonder what EL_st states are for, looks like it's always E_st

	FILE *cmFile=ThrowingFopen(cmFileName,"rt");
	char version[256]="";
	fscanf(cmFile,"### cove %s\n",version);
	if (strcmp(version,"V2")!=0) {
		throw SimpleStringException("Cove-format file didn't begin with \"### cove V2\" --> it's not in Cove V2 format.");
	}
	int numNodes;
	fscanf(cmFile, "%d \tnodes\n", &numNodes);

	// can allocate Infernal structure now, possibly using more memory than we need (because we're using an upper bound for required # of states), but I don't think that matters
	int upperBoundOnNumStates=numNodes*6;
	int upperBoundOnNumNodes=numNodes*2;
	cm=CreateCM(upperBoundOnNumNodes,upperBoundOnNumStates);
	assert(cm->nodes==upperBoundOnNumNodes); // just verifying my expectation of the CreateCM function
	assert(cm->flags==0); // shouldn't have any set, right?

	cm->name=(char *)MallocOrDie(strlen(cmFileName)+1);
	strcpy(cm->name,cmFileName);

	cm->null[0]=cm->null[1]=cm->null[2]=cm->null[3]=0.25; // default NULL model

	vector<int> coveLeftChild(numNodes),coveRightChild(numNodes);
	vector<tRNAscanSE::NodeType> coveNodeType(numNodes);
	vector<tRNAscanSE::SingletEmit> coveEmitIL(numNodes),coveEmitIR(numNodes),coveEmitML(numNodes),coveEmitMR(numNodes);
	vector<tRNAscanSE::TransitionMatrix> coveTransitionMatrix(numNodes);
	vector<tRNAscanSE::PairEmit> coveEmitMP(numNodes);
	vector<int> infernalNodeToCoveNode(upperBoundOnNumNodes); // upper bound on # of nodes, if every cove node were an END node

	// first, read in the file
	for (int currNode=0; currNode<numNodes; currNode++) {
		int currNodeFromFile;
		fscanf(cmFile, "### node %d",&currNodeFromFile);
		if (currNode!=currNodeFromFile) {
			throw SimpleStringException("Cove V2 file's nodes are not in sequential order, but this dinky code assumes they are.\n");
		}

		int thisCoveNodeType;
		fscanf(cmFile, " type %d\n", &thisCoveNodeType);
		coveNodeType[currNode]=(tRNAscanSE::NodeType)thisCoveNodeType;

		int leftCoveChildNode,rightCoveChildNode;
		fscanf(cmFile, "%d  %d\n", &leftCoveChildNode,&rightCoveChildNode);
		if (leftCoveChildNode!=-1 && leftCoveChildNode!=currNode+1) {
			throw SimpleStringException("Cove V2 file's left child (of node #%d) is not the next node or -1, which I was assuming it would be",currNode);
		}
		coveLeftChild[currNode]=leftCoveChildNode;
		coveRightChild[currNode]=rightCoveChildNode;

		for (int from=0; from<tRNAscanSE::STATETYPES; from++) {
			for (int to=0; to<tRNAscanSE::STATETYPES; to++) {
				fscanf(cmFile,"%f ",&(coveTransitionMatrix[currNode].array[from][to]));
			}
			fscanf(cmFile, "\n");
		}

		for (int i = 0; i < 4; i++) {
			fscanf(cmFile, "%f ", &(coveEmitIL[currNode].array[i]));
		}
		fscanf(cmFile, "# INSL\n");
		for (int i = 0; i < 4; i++) {
			fscanf(cmFile, "%f ", &(coveEmitIR[currNode].array[i]));
		}
		fscanf(cmFile, "# INSR\n");
		for (int i = 0; i < 4; i++) {
			for (int j=0; j<4; j++) {
				fscanf(cmFile, "%f ", &(coveEmitMP[currNode].array[i][j]));
			}
			fscanf(cmFile, "# MATP\n");
		}
		for (int i = 0; i < 4; i++) {
			fscanf(cmFile, "%f ", &(coveEmitML[currNode].array[i]));
		}
		fscanf(cmFile, "# MATL\n");
		for (int i = 0; i < 4; i++) {
			fscanf(cmFile, "%f ", &(coveEmitMR[currNode].array[i]));
		}
		fscanf(cmFile, "# MATR\n");
	}

	fclose(cmFile);

	// now, linearize the states, as in Infernal
	int nextInfernalState=0;
	int nextInfernalNode=0;
	for (int currNode=0; currNode<numNodes; currNode++) {

		infernalNodeToCoveNode[nextInfernalNode]=currNode;

		cm->nodemap[nextInfernalNode]=nextInfernalState;
		tRNAscanSE::NodeTypeInfo nodeTypeInfo=tRNAscanSE::GetNodeTypeInfo (coveNodeType[currNode]);

		cm->ndtype[nextInfernalNode]=nodeTypeInfo.infernalNodeType;
		if (coveNodeType[currNode]==tRNAscanSE::BIFURC_NODE) {
			if (coveLeftChild[currNode]==-1 && coveRightChild[currNode]==-1) {
				// it's really an end node
				cm->ndtype[nextInfernalNode]=END_nd;
				cm->sttype[nextInfernalState]=EL_st;
				cm->ndidx[nextInfernalState]=nextInfernalNode;
				cm->cfirst[nextInfernalState]=-1;
				cm->cnum[nextInfernalState]=0;
				if (currNode+1==numNodes || !useLocalEndStates) {
					// final end state
					cm->sttype[nextInfernalState]=E_st;
				}
			}
			else {
				if (coveLeftChild[currNode]==-1 || coveRightChild[currNode]==-1) {
					throw SimpleStringException("Cove V2 had Bif/End node (#%d) that had exactly one child.  I agree that's senseless.",currNode);
				}

				// bif node
				cm->sttype[nextInfernalState]=B_st;
				cm->ndidx[nextInfernalState]=nextInfernalNode;
				// set up the children later
			}
			nextInfernalState++;
		}
		else {

			int numSplitSetStates=nodeTypeInfo.numSplitSetStates;
			int numInsertStates=nodeTypeInfo.numInsertStates;

			int nextSplitSetStates;
			if (coveLeftChild[currNode]==-1 && coveRightChild[currNode]==-1) {
				// next node is a phantom END node
				nextSplitSetStates=1;
			}
			else {
				assert(currNode+1<numNodes); // else this should have been an END node
				tRNAscanSE::NodeTypeInfo nextNodeTypeInfo=tRNAscanSE::GetNodeTypeInfo (coveNodeType[currNode+1]);
				nextSplitSetStates=nextNodeTypeInfo.numSplitSetStates;
			}

			int firstChildState=nextInfernalState+numSplitSetStates;
			for (int currState=0; currState<numSplitSetStates+numInsertStates; currState++) {
				tRNAscanSE::StateType coveStateType=nodeTypeInfo.states[currState];
				tRNAscanSE::StateTypeInfo stateTypeInfo=GetStateTypeInfo(coveStateType);
				cm->sttype[nextInfernalState]=stateTypeInfo.infernalStateType;
				if (nodeTypeInfo.firstStateIsStart && currState==0) {
					cm->sttype[nextInfernalState]=S_st;
				}
				cm->ndidx[nextInfernalState]=nextInfernalNode;
				cm->cfirst[nextInfernalState]=firstChildState;
				cm->cnum[nextInfernalState]=nextSplitSetStates+nodeTypeInfo.numInsertStates;
				if (cm->sttype[nextInfernalState]==IR_st && numInsertStates==2) {
					// special case -- IR doesn't go to IL
					cm->cfirst[nextInfernalState]++;
					cm->cnum[nextInfernalState]--;
				}
				
				// emits
				switch (stateTypeInfo.numEmits) {
					case 0:
						// nothing to do for emits
						break;
					case 1:
						if (coveStateType==tRNAscanSE::INSL_ST) { // MSVC++ screws up tabbing if I use switch/case
							CopySingletEmitFromCove(cm,nextInfernalState,coveEmitIL[currNode]);
						}
						if (coveStateType==tRNAscanSE::INSR_ST) { // MSVC++ screws up tabbing if I use switch/case
							CopySingletEmitFromCove(cm,nextInfernalState,coveEmitIR[currNode]);
						}
						if (coveStateType==tRNAscanSE::MATL_ST) { // MSVC++ screws up tabbing if I use switch/case
							CopySingletEmitFromCove(cm,nextInfernalState,coveEmitML[currNode]);
						}
						if (coveStateType==tRNAscanSE::MATR_ST) { // MSVC++ screws up tabbing if I use switch/case
							CopySingletEmitFromCove(cm,nextInfernalState,coveEmitMR[currNode]);
						}
						break;
					case 2:
						{
							for (int leftNuc=0; leftNuc<4; leftNuc++) {
								for (int rightNuc=0; rightNuc<4; rightNuc++) {
									float emitProb=coveEmitMP[currNode].array[leftNuc][rightNuc];
									cm->e[nextInfernalState][GetPairIndex(leftNuc,rightNuc)]=emitProb;
									cm->esc[nextInfernalState][GetPairIndex(leftNuc,rightNuc)]=(float)(sreLOG2(emitProb));
								}
							}
						}
						break;
					default:
						assert(false);
						break;
				}

				nextInfernalState++;
			}

			// (annoying) special case: a node of regular type (like MATL) can also be an END node.  In this case, we have to put in the INFERNAL END node afterwards.
			if (coveLeftChild[currNode]==-1 && coveRightChild[currNode]==-1) {
				nextInfernalNode++;
				infernalNodeToCoveNode[nextInfernalNode]=currNode;
				cm->ndtype[nextInfernalNode]=END_nd;
				cm->nodemap[nextInfernalNode]=nextInfernalState;
				cm->sttype[nextInfernalState]=EL_st;
				cm->ndidx[nextInfernalState]=nextInfernalNode;
				cm->cfirst[nextInfernalState]=-1;
				cm->cnum[nextInfernalState]=0;
				if (currNode+1==numNodes || !useLocalEndStates) {
					// final end state
					cm->sttype[nextInfernalState]=E_st;
				}
				nextInfernalState++;
			}
		}

		nextInfernalNode++;
	}
	cm->M=nextInfernalState;
	cm->nodes=nextInfernalNode;

	// do transitions, now that every state has its type
	for (int state=0; state<cm->M; state++) {
		int currInfernalNode=cm->ndidx[state];
		int currNode=infernalNodeToCoveNode[currInfernalNode];

		if (cm->sttype[state]!=B_st) {
			tRNAscanSE::StateType coveFromStateType=tRNAscanSE::InfernalToCoveStateType(cm->sttype[state]);
			for (int i=0; i<cm->cnum[state]; i++) {
				int childState=cm->cfirst[state]+i;
				tRNAscanSE::StateType coveToStateType=tRNAscanSE::InfernalToCoveStateType(cm->sttype[childState]);
				float p=coveTransitionMatrix[currNode].array[coveFromStateType][coveToStateType];
				cm->t[state][i]=p;
				cm->tsc[state][i]=(float)(sreLOG2(p));
			}
		}
	}

	// now set up right children of BIF nodes
	// first reverse-map nodes
	vector<int> coveNodeToInfernalNode(numNodes);
	for (int infernalNode=0; infernalNode<cm->nodes; infernalNode++) {
		coveNodeToInfernalNode[infernalNodeToCoveNode[infernalNode]]=infernalNode;
	}
	// now do it
	for (int currNode=0; currNode<numNodes; currNode++) {
		if (coveNodeType[currNode]==tRNAscanSE::BIFURC_NODE) {
			if (coveLeftChild[currNode]!=-1 && coveRightChild[currNode]!=-1) {
				int bifState=cm->nodemap[coveNodeToInfernalNode[currNode]];
				int leftChildState=cm->nodemap[coveNodeToInfernalNode[coveLeftChild[currNode]]];
				int rightChildState=cm->nodemap[coveNodeToInfernalNode[coveRightChild[currNode]]];
				cm->cfirst[bifState]=leftChildState;
				cm->cnum[bifState]=rightChildState;
			}
		}
	}

	// now set up parent pointers & full stid
	// first clear parents
	for (int state=0; state<cm->M; state++) {
		cm->plast[state]=-1;
		cm->pnum[state]=0;
	}
	// now set everything
	for (State state=GetFirstState(); state!=GetLastState(); state++) {
		cm->stid[StateToInt(state)]=DeriveUniqueStateCode(GetNodeType(GetNode(state)), GetStateType(state));
		if (IsBifurcation(state)) {
			State childState=GetLeftBifurcationChild(state);
			cm->plast[StateToInt(childState)]=StateToInt(state);
			cm->pnum[StateToInt(childState)]++;
			childState=GetRightBifurcationChild(state);
			cm->plast[StateToInt(childState)]=StateToInt(state);
			cm->pnum[StateToInt(childState)]++;
		}
		else {
			for (int childNum=0; childNum<GetNumChildren(state); childNum++) {
				State childState=GetNthChildState(state,childNum);
				cm->plast[StateToInt(childState)]=StateToInt(state);
				cm->pnum[StateToInt(childState)]++;
			}
		}
	}
}
const char *CovarianceModelBase::GetName () const {
	return "";
}
void CovarianceModelBase::CopyFrom (const CovarianceModelBase& t)
{
	Destruct();

	isRsearchCM=t.isRsearchCM;
	isSmithWaterCM=t.isSmithWaterCM;
	cmFileNameLoaded=t.cmFileNameLoaded;

	cm=CreateCM(t.GetNumNodes(),t.GetNumStates());

	if (t.GetName()==NULL) {
		cm->name=NULL;
	}
	else {
		cm->name=(char *)MallocOrDie(strlen(t.GetName())+1);
		strcpy(cm->name,t.GetName());
	}

	COPY_ARRAY(sttype,t.GetNumStates()+1);
	COPY_STATE_ARRAY(ndidx);
	COPY_ARRAY(stid,t.GetNumStates()+1);
	COPY_STATE_ARRAY(cfirst);
	COPY_STATE_ARRAY(cnum);
	COPY_STATE_ARRAY(plast);
	COPY_STATE_ARRAY(pnum);

	COPY_ARRAY(nodemap,t.GetNumNodes());
	COPY_ARRAY(ndtype,t.GetNumNodes());

	COPY_ARRAY_2D(t,t.GetNumStates(),MAXCONNECT);
	COPY_ARRAY_2D(e,t.GetNumStates(),Alphabet_size*Alphabet_size);
	COPY_STATE_ARRAY(begin);
	COPY_STATE_ARRAY(end);
	COPY_ARRAY_2D(tsc,t.GetNumStates(),MAXCONNECT);
	COPY_ARRAY_2D(esc,t.GetNumStates(),Alphabet_size*Alphabet_size);
	COPY_STATE_ARRAY(beginsc);
	COPY_STATE_ARRAY(endsc);

	cm->flags  = t.cm->flags;

	for (int nuc=0; nuc<Alphabet_size; nuc++) {
		cm->null[nuc]=t.cm->null[nuc];
	}
}
void CovarianceModelBase::TakeFrom (CM_t *t)
{
	Destruct();

	cmMemoryOwned=true;
	cm=t;
	cmfp=NULL;
}
void CovarianceModelBase::MirrorFrom (CM_t *t)
{
	Destruct();

	cmMemoryOwned=false;
	cm=t;
	cmfp=NULL;
}

void CovarianceModelBase::Realloc2d(float **(&array),int oldNumStates,int newNumStates,int sizeDim2)
{
	float **old=array;
	array=FMX2Alloc(newNumStates, sizeDim2);
	for (int i=0; i<oldNumStates; i++) {
		for (int j=0; j<sizeDim2; j++) {
			array[i][j]=old[i][j];
		}
	}
	FMX2Free(old);
}


void CovarianceModelBase::Init (int numStates)
{
	cm=CreateCM(0,numStates);

	char dummyName[]="unnamed";
	cm->name = sre_strdup(dummyName, (int)(strlen(dummyName)));
}
void CovarianceModelBase::DivideOutUniformNullModel(void)
{
	int v, x, y;
	for (v = 0; v < cm->M; v++)
	{
		if (cm->sttype[v] == MP_st)
			for (x = 0; x < Alphabet_size; x++)
				for (y = 0; y < Alphabet_size; y++) {
					cm->e[v][x*Alphabet_size+y] *= 16.0;
					cm->esc[v][x*Alphabet_size+y] = (float)sreLOG2(cm->e[v][x*Alphabet_size+y]);
				}

		if (cm->sttype[v] == ML_st || cm->sttype[v] == MR_st ||
			cm->sttype[v] == IL_st || cm->sttype[v] == IR_st)
			for (x = 0; x < Alphabet_size; x++) {
				cm->e[v][x] *= 4.0;
				cm->esc[v][x] = (float)sreLOG2(cm->e[v][x]);
			}
	}
}
void CovarianceModelBase::CMRenormalize (void)
{
	::CMRenormalize(cm);// doesn't handle local ends, so I wrote this code myself
	/*
	double totalProb;
	totalProb=0;
	if (DoLocal()) {
		for (v=0; v<cm->M; v++) {
			totalProb += cm->begin[v];
		}
		for (v=0; v<cm->M; v++) {
			cm->begin[v]=ToProb(cm->begin[v],totalProb);
		}
	}
	for (State state=GetFirstState(); state!=GetLastState(); state++) {
		v=StateToInt(state);
		if (IsEmittingState(state)) {
			int n=Alphabet_size;
			if (EmitsLeftAndRight(state)) {
				n=Alphabet_size*Alphabet_size;
			}
			totalProb=0;
			for (int nuc=0; nuc<n; nuc++) {
				totalProb += cm->e[v][nuc];
			}
			for (int nuc=0; nuc<n; nuc++) {
				cm->e[v][nuc]=ToProb(cm->e[v][nuc],totalProb);
			}
		}

		totalProb=0;
		if (DoLocal()) {
			totalProb += cm->end[v];
		}
		for (int child=0; child<GetNumChildren(state); child++) {
			totalProb += cm->t[v][child];
		}

		bool doIt=true;
		if (DoLocal()) {
			doIt=!(cm->end[v]==0 && GetNumChildren(state)==0);
		}
		if (doIt) {
			if (DoLocal()) {
				cm->end[v]=ToProb(cm->end[v],totalProb);
			}
		}
	}
	*/

	SetLODScores();

	if (enableInsertHack) {
		::CMHackInsertScores(cm);
	}
}
void CovarianceModelBase::SetLODScores (void)
{
	int v, x, y;
	// copied from ::CMLogoddsify, but in this case we just want to Log
	for (v = 0; v < cm->M; v++)
	{
		if (cm->sttype[v] != B_st && cm->sttype[v] != E_st)
			for (x = 0; x < cm->cnum[v]; x++)
				cm->tsc[v][x] = (float)sreLOG2(cm->t[v][x]);

		if (cm->sttype[v] == MP_st)
			for (x = 0; x < Alphabet_size; x++)
				for (y = 0; y < Alphabet_size; y++)
					cm->esc[v][x*Alphabet_size+y] = (float)sreLOG2(cm->e[v][x*Alphabet_size+y]);

		if (cm->sttype[v] == ML_st || cm->sttype[v] == MR_st ||
			cm->sttype[v] == IL_st || cm->sttype[v] == IR_st)
			for (x = 0; x < Alphabet_size; x++)
				cm->esc[v][x] = (float)sreLOG2(cm->e[v][x]);

		/* These work even if begin/end distributions are inactive 0's,
		* sreLOG2 will set beginsc, endsc to -infinity.
		*/
		cm->beginsc[v] = (float)sreLOG2(cm->begin[v]);
		cm->endsc[v]   = (float)sreLOG2(cm->end[v]);
	}
}
void CovarianceModelBase::MultiplyEmitsBy (float x)
{
	float logX=(float)(sreLOG2(x));
	for (State state=GetFirstState(); state!=GetLastState(); state++) {
		if (IsEmittingState(state)) {
			if (EmitsLeftAndRight(state)) {
				for (int l=0; l<Alphabet_size; l++) {
					for (int r=0; r<Alphabet_size; r++) {
						float p=GetPairEmissionScore(state,l,r);
						p += 2*logX;
						SetPairEmissionScore(state,l,r,p);
						SetPairEmissionProbDirectly(state,l,r,(float)(pow2(p)));
					}
				}
			}
			else {
				for (int nuc=0; nuc<Alphabet_size; nuc++) {
					float p=GetSingletEmissionScore(state,nuc);
					p += logX;
					SetSingletEmissionLogScore(state,nuc,p);
					SetEmissionProbDirectly(state,nuc,(float)(pow2(p)));
				}
			}
		}
	}
}
void CovarianceModelBase::HackInsertScoresToStrictProbs (void)
{
	if (enableInsertHack) {
		for (State state=GetFirstState(); state!=GetLastState(); state++) {
			if (IsInsertState(state)) {
				for (int nuc=0; nuc<Alphabet_size; nuc++) {
					SetEmissionProbDirectly(state,nuc,0.25);
					SetSingletEmissionLogScore(state,nuc,-2.0);
				}
			}
		}
	}
}
void CovarianceModelBase::SetTransitionProbDirectly (State state,int child,float prob)
{
	cm->t[StateToInt(state)][child]=prob;
}
float CovarianceModelBase::GetTransitionProbDirectly (State state,int child) const
{
	return cm->t[StateToInt(state)][child];
}
void CovarianceModelBase::SetEmissionProbDirectly (State state,int nuc,float prob)
{
	cm->e[StateToInt(state)][nuc]=prob;
}
float CovarianceModelBase::GetEmissionProbDirectly (State state,int nuc) const
{
	return cm->e[StateToInt(state)][nuc];
}
float CovarianceModelBase::GetPairEmissionProbDirectly (State state,int leftNuc,int rightNuc) const
{
	return cm->e[StateToInt(state)][leftNuc*Alphabet_size + rightNuc];
}
void CovarianceModelBase::SetPairEmissionProbDirectly(State state,int leftNuc,int rightNuc,float prob)
{
	cm->e[StateToInt(state)][leftNuc*Alphabet_size + rightNuc]=prob;
}
void CovarianceModelBase::NormalizeEmissionsToStrictProbabilitiesViaProbabilities (void)
{
	for (State state=GetFirstState(); state!=GetLastState(); state++) {

		switch (GetNumSymbolsEmitted(state)) {
			case 1:
				{
					double totalProb=0;
					for (int nuc=0; nuc<Alphabet_size; nuc++) {
						totalProb += GetEmissionProbDirectly(state,nuc);
					}
					for (int nuc=0; nuc<Alphabet_size; nuc++) {
						float t=GetEmissionProbDirectly(state,nuc);
						float prob=(float)(t/totalProb);
						if (t==0 && totalProb==0) {
							prob=(float)(1.0)/(float)(Alphabet_size);
						}
						SetEmissionProbDirectly(state,nuc,prob);
						SetSingletEmissionLogScore(state,nuc,(float)(sreLOG2(prob)));
					}
				}
				break;
			case 2:
				{
					double totalProb=0;
					for (int leftNuc=0; leftNuc<Alphabet_size; leftNuc++) {
						for (int rightNuc=0; rightNuc<Alphabet_size; rightNuc++) {
							totalProb += GetPairEmissionProbDirectly(state,leftNuc,rightNuc);
						}
					}
					for (int leftNuc=0; leftNuc<Alphabet_size; leftNuc++) {
						for (int rightNuc=0; rightNuc<Alphabet_size; rightNuc++) {
							float t=GetPairEmissionProbDirectly(state,leftNuc,rightNuc);
							float prob=(float)(t/totalProb);
							if (t==0 && totalProb==0) {
								prob=(float)(1.0)/(float)(Alphabet_size*Alphabet_size);
							}
							SetPairEmissionProbDirectly(state,leftNuc,rightNuc,prob);
							SetPairEmissionScore(state,leftNuc,rightNuc,(float)(sreLOG2(prob)));
						}
					}
				}
				break;
		}
	}
}
void CovarianceModelBase::NormalizeTransitionsToStrictProbabilitiesViaProbabilities (void)
{
	for (State state=GetFirstState(); state!=GetLastState(); state++) {

		double totalProb=0;
		for (int child=0; child<GetNumChildren(state); child++) {
			totalProb += GetTransitionProbDirectly(state,child);
		}
		for (int child=0; child<GetNumChildren(state); child++) {
			float t=GetTransitionProbDirectly(state,child);
			float prob=(float)(t/totalProb);
			if (t==0 && totalProb==0) {
				prob=1;
			}
			SetTransitionProbDirectly(state,child,prob);
			SetTransitionLogScore(state,child,(float)(sreLOG2(prob)));
		}
	}
}
void CovarianceModelBase::ZeroAllEmitProbs (void)
{
	for (State state=GetFirstState(); state!=GetLastState(); state++) {
		assert(GetNumSymbolsEmitted(state)!=2); // only meant for HMMs
		if (IsEmitting(state)) {
			for (int nuc=0; nuc<Alphabet_size; nuc++) {
				SetEmissionProbDirectly(state,nuc,(float)0);
			}
		}
	}
}
void CovarianceModelBase::ZeroAllTransitionProbsExceptSelfLoops (void)
{
	for (State state=GetFirstState(); state!=GetLastState(); state++) {
		for (int child=0; child<GetNumChildren(state); child++) {
			State toState=GetNthChildState(state,child);
			if (state!=toState) { // no a self-loop
				SetTransitionProbDirectly(state,child,(float)0);
			}
		}
	}
}
void CovarianceModelBase::ZeroAllTransitionProbs (void)
{
	for (State state=GetFirstState(); state!=GetLastState(); state++) {
		for (int child=0; child<GetNumChildren(state); child++) {
			SetTransitionProbDirectly(state,child,(float)0);
		}
	}
}
void CovarianceModelBase::SetFirstChild (State state,State firstChild)
{
	cm->cfirst[StateToInt(state)]=StateToInt(firstChild);
}
void CovarianceModelBase::SetNumChildren (State state,int numChildren)
{
	assert(numChildren>=0 && numChildren<MAXCONNECT);
	cm->cnum[StateToInt(state)]=numChildren;
}
void CovarianceModelBase::SetNoChildren (State state)
{
	cm->cfirst[StateToInt(state)]=-1;
	cm->cnum[StateToInt(state)]=0;
}
void CovarianceModelBase::SetStateType (State state,int stateType)
{
	cm->sttype[StateToInt(state)]=stateType;
}
void CovarianceModelBase::SetTransitionLogScore (State state,int childNum,float score)
{
	assert(childNum>=0 && childNum<GetNumChildren(state));
	cm->tsc[StateToInt(state)][childNum]=score;
}
void CovarianceModelBase::SetSingletEmissionLogScore (State state,int symbol,float score)
{
	assert(symbol>=0 && symbol<Alphabet_size);
	cm->esc[StateToInt(state)][symbol]=score;
}
void CovarianceModelBase::SetEndsc (State state,float score)
{
	cm->endsc[StateToInt(state)]=score;
}
void CovarianceModelBase::SetBeginsc (State state,float score)
{
	cm->beginsc[StateToInt(state)]=score;
}
CMConsensus_t *CovarianceModelBase::CreateCMConsensus (float x,float y) const
{
	CM_t *nonConstCM=(CM_t *)cm;
	return ::CreateCMConsensus(nonConstCM,x,y);
}
int CovarianceModelBase::GetNumTransitionsTotal () const
{
	int n=0;
	for (State state=GetFirstState(); state!=GetLastState(); state++) {
		n += GetNumChildren(state);
	}
	return n;
}
int CovarianceModelBase::GetNumBifTotal () const
{
	int n=0;
	for (State state=GetFirstState(); state!=GetLastState(); state++) {
		if (IsBifurcation(state)) {
			n++;
		}
	}
	return n;
}
int CovarianceModelBase::GetNumBifChildrenTotal () const
{
	return GetNumBifTotal();
}


CmNodeLength::CmNodeLength (void)
{
	numNodes=0;
}
CmNodeLength::CmNodeLength (int numNodes_)
{
	numNodes=numNodes_;
	lengthVector.resize(numNodes);
}
CmNodeLength::~CmNodeLength ()
{
}
int CmNodeLength::GetNumNodes () const
{
	return numNodes;
}
void CmNodeLength::Zero (void)
{
	Len l;
	l.left=0;
	l.right=0;
	lengthVector.assign(numNodes,l);
}
int CmNodeLength::fileFormat=0;
void CmNodeLength::Save (const char *fileName)
{
	FILE *f=ThrowingFopen(fileName,"wt");
	fprintf(f,"format\t%d\n",fileFormat);
	fprintf(f,"numNodes\t%d\n",numNodes);
	for (int i=0; i<numNodes; i++) {
		fprintf(f,"%d\t%d\t%d\n",i,lengthVector[i].left,lengthVector[i].right);
	}
	fclose(f);
}
void CmNodeLength::Load (const char *fileName)
{
	CommaSepFileReader f(fileName,'\t');

	f.ReadLineOrFail();
	int format=f.GetFieldAsInt(1);
	if (format>fileFormat) {
		throw SimpleStringException("CmNodeLength::Load: file format is unknown -- either this is an old executable, or the file is corrupt.");
	}

	f.ReadLineOrFail();
	numNodes=f.GetFieldAsInt(1);
	lengthVector.resize(numNodes);

	for (int i=0; i<numNodes; i++) {
		f.ReadLineOrFail();
		lengthVector[i].left=f.GetFieldAsInt(1);
		lengthVector[i].right=f.GetFieldAsInt(2);
	}
}
void CmNodeLength::Set (CovarianceModel::Node node,int leftLength,int rightLength)
{
	Set(CovarianceModel::NodeToInt(node),leftLength,rightLength);
}
void CmNodeLength::Set (int node,int leftLength,int rightLength)
{
	lengthVector[node].left=leftLength;
	lengthVector[node].right=rightLength;
}
void CmNodeLength::IncNodeByEmits (CovarianceModel::Node node,bool emitLeft,bool emitRight)
{
	int i=CovarianceModel::NodeToInt(node);
	if (emitLeft) {
		lengthVector[i].left++;
	}
	if (emitRight) {
		lengthVector[i].right++;
	}
}
void CmNodeLength::MaxWith (const CmNodeLength& nl)
{
	for (int i=0; i<numNodes; i++) {
		lengthVector[i].left=std::max(lengthVector[i].left,nl.lengthVector[i].left);
		lengthVector[i].right=std::max(lengthVector[i].right,nl.lengthVector[i].right);
	}
}
int CmNodeLength::GetLeft (CovarianceModel::Node node) const
{
	return GetLeft(CovarianceModel::NodeToInt(node));
}
int CmNodeLength::GetLeft (int node) const
{
	return lengthVector[node].left;
}
int CmNodeLength::GetRight (CovarianceModel::Node node) const
{
	return GetRight(CovarianceModel::NodeToInt(node));
}
int CmNodeLength::GetRight (int node) const
{
	return lengthVector[node].right;
}
