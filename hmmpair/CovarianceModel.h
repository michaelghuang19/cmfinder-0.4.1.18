
// wrapper of CM_t object from infernal-0.54

// states I've added
// dummy states we use (& keep in the HMM) for the end of a non-bifuricating block.  It's functionally equivalent to a D_st
#define PASSTHRU_st 16
// nice sentinel value
#define INVALID_st -1

class CovarianceModelBase {
protected:
	CMFILE *cmfp;
	CM_t *cm;
	bool cmMemoryOwned;
	bool isRsearchCM;
	bool isSmithWaterCM; // implies isRsearchCM
	std::string cmFileNameLoaded;

	int uniqueId;
	static int nextUniqueId;

	void Destruct (void);

	static bool enableInsertHack;
	
	// was for debugging: void fwrite(void *p,size_t size,size_t count,FILE *file);

	void Realloc2d(float **(&array),int oldNumStates,int newNumStates,int sizeDim2);
public:
	CovarianceModelBase (void);
	~CovarianceModelBase ();

	// useful for caching & verifying that we don't have 2 CMs at the same address
	// may be conservative about being equal (i.e. if you copy a CM from another using CopyFrom, they'll have different uniqueId's, even though they're really the same)
	int GetUniqueId (void) const; 

	const char *GetCmFileName (void) const; // will only be set if you created it with Load
	int rsearchQueryLen;
	bool /*success*/ Load (char *cmFileName,bool doLocalAlignment);
	void Save (char *cmFileName);
	void SaveRsearchCM (char *cmFileName,std::string programParams="");
	void CopyFrom (const CovarianceModelBase& t);
	void TakeFrom (CM_t *cm);
	void MirrorFrom (CM_t *cm);
	void Init (int numStates);
	// load Cove 2-format used by tRNAscan-SE (Lowe & Eddy, 1997)
	void LoadCove2 (const char *cmFileName);
	static void DisableInsertHack (void);
	void DumpCsv (FILE *csv) const; // I hope a useful format for debugging
	void DumpCsv (const char *csvFileName) const;
	bool IsRsearchCM (void) const { return isRsearchCM; }
	void SetIsRsearchCM (bool x) { isRsearchCM=x; }
	bool IsInit (void) const { return cm!=NULL; }
	bool IsStructurallyTheSame (const CovarianceModelBase& t) const;
	void SetLODScores (void);
	void MorphToSmithWater(float gapOpen,float gapExtend,bool modifyTransitions=true);

	class NodeOrState {
		friend class CovarianceModelBase;
	protected:
		int index;
		inline NodeOrState (int _index) {
			index=_index;
		}
		inline int ToInt (void) const {
			return index;
		}
	public:
		inline NodeOrState () {
		}
		inline void operator = (const NodeOrState& t) {
			index=t.index;
		}
		inline NodeOrState (const NodeOrState& t) {
			*this=t;
		}
		inline bool operator == (const NodeOrState& t) const {
			return index==t.index;
		}
		inline bool operator != (const NodeOrState& t) const {
			return !(*this==t);
		}
		inline void operator ++ (int dummy) {
			index++;
		}
		inline void operator -- (int dummy) {
			index--;
		}
		void operator += (int t) {
			index += t;
		}
		inline bool GreaterThanZero (void) const {
			return index>=0;
		}
	};
	class Node : public NodeOrState {
		friend class CovarianceModelBase;
	protected:
		inline Node (int index) : NodeOrState(index) {}
	public:
		inline Node () : NodeOrState () {}
		inline Node (const Node& t) : NodeOrState(t) {}
		inline Node PlusInt (int i) const {
			return Node (index+i);
		}
		bool operator >= (const Node& t) const {
			return index>=t.index;
		}
		bool operator <= (const Node& t) const {
			return index<=t.index;
		}
		bool operator < (const Node& t) const {
			return index<t.index;
		}
		bool operator > (const Node& t) const {
			return index>t.index;
		}
	};
	class State : public NodeOrState {
		friend class CovarianceModelBase;
	protected:
		inline State (int index) : NodeOrState(index) {}
	public:
		inline State () : NodeOrState () {}
		inline State (const State& t) : NodeOrState(t) {}
		inline State PlusInt (int i) const {
			return State(index+i);
		}
		bool operator >= (const State& t) const {
			return index>=t.index;
		}
		bool operator <= (const State& t) const {
			return index<=t.index;
		}
		bool operator < (const State& t) const {
			return index<t.index;
		}
		bool operator > (const State& t) const {
			return index>t.index;
		}
	};
	typedef std::list<State> StateList;
	typedef std::list<Node> NodeList;
	typedef vector<Node> NodeVector;

	// __ONLY__ for calling code in infernal
	inline CM_t * GetCM (void) {
		return cm;
	}
	inline const CM_t *GetCM (void) const {
		return cm;
	}
	CMConsensus_t *CreateCMConsensus (float x,float y) const;
	//float InsideGlobalEpsilonEnds_UniformLOD(char *dsq, int L, int W,CykscanStats& cykscanStats) const;
	//NoUnderflowDouble InsideGlobalEpsilonEndsProb(char *dsq, int L, int W,CykscanStats& cykscanStats) const;
	//void CollectSamplesForEvalue (EvalueScoreSamples& evalueScoreSamples,int *gc_count, int *partitions, int num_partitions, int D, int num_samples,int sample_length, zrand::ZRandom *rander) const;
	//void serial_make_histogram (int *gc_count, int *partitions, int num_partitions,int D, int num_samples, int sample_length, float *lambda, float *K,zrand::ZRandom *rander) const;
	float CYKDivideAndConquer(char *seq,int L,State state,int i,int j,Parsetree_t **retr=NULL) const;
	float CYKInsideScore (char *seq,int L,State state,int i,int j) const;
	//void CYKScanZasha (char *dsq, int L, int W,int *ret_nhits, int **ret_hitr, int **ret_hiti, int **ret_hitj, float **ret_hitsc,CykscanStats& cykscanStats,bool null3cm,ScanParams& scanParams) const;
	//void CYKScan_BackgroundScore (char *dsq, int L, int W,int *ret_nhits, int **ret_hitr, int **ret_hiti, int **ret_hitj, float **ret_hitsc,float backgroundScore,bool null3cm,ScanParams& scanParams) const;
	//void CYKScan (char *dsq, int L, int W,int *ret_nhits, int **ret_hitr, int **ret_hiti, int **ret_hitj, float **ret_hitsc,bool null3cm,ScanParams& scanParams) const;
	//float /*best score*/ Rsearch_CYKScan(char *rnaSequence,int rnaSequenceLen,float scoreThreshold,int windowLen,scan_results_t *results) const;
	Fancyali_t *CreateFancyAli(Parsetree_t *tr, CMConsensus_t *cons, char *rnaSequence) const;
	void ParsetreeDump(FILE *out, Parsetree_t *tr, char *rnaSequence) const;
	float ParsetreeScore(Parsetree_t *tr,char *rnaSequence) const;

	inline float GetLocalEndSelfScore () const {
		return cm->el_selfsc;
	}

	void CMRenormalize (void); // calls CMRenormalize in cm.c
	void DivideOutUniformNullModel(void);
	void HackInsertScoresToStrictProbs (void);
	void ZeroAllTransitionProbsExceptSelfLoops (void);
	void ZeroAllTransitionProbs (void);
	void ZeroAllEmitProbs (void);
	void NormalizeTransitionsToStrictProbabilitiesViaProbabilities (void);
	void NormalizeEmissionsToStrictProbabilitiesViaProbabilities (void);
	void MultiplyEmitsBy (float x);

	const char *GetName (void) const;
	static inline int StateToInt (State state) {
		return state.ToInt();
	}
	static inline State IntToState (int s) {
		return State(s);
	}
	inline int GetNumStates (void) const {
		return cm->M;
	}
	int GetNumTransitionsTotal () const;
	int GetNumBifTotal () const;
	int GetNumBifChildrenTotal () const;
	inline bool IsStateValid (State state) const {
		int stateIndex=StateToInt(state);
		return stateIndex>=0 && stateIndex<cm->M;
	}
	inline State GetStartState (void) const {
		return State(0);
	}
	inline int GetStateType (State state) const {
		return cm->sttype[StateToInt(state)];
	}
	inline int GetCombinedStateAndNodeType (State state) const {
		return cm->stid[StateToInt(state)];
	}
	inline const char *GetStateTypeName (State state) const {
		switch (GetStateType(state)) {
			case INVALID_st:
				return "INVALID";
			case PASSTHRU_st:
				return "PASSTHRU";
		}
		return Statetype(GetStateType (state));
	}
	inline State GetFirstState (void) const {
		return State(0);
	}
	inline State GetLastState (void) const {
		return State(GetNumStates());
	}
	inline State GetActualLastState (void) const {
		State s(GetLastState());
		s--;
		return s;
	}
	static inline State GetInvalidState (void) {
		return State(-1);
	}
	static inline Node GetInvalidNode (void) {
		return Node(-1);
	}
	inline bool IsValidState (State state) const {
		int s=StateToInt(state);
		return s>=0 && s<GetNumStates();
	}
	inline bool IsValidNode (Node node) const {
		return node>=GetFirstNode() && node<GetLastNode();
	}
	static inline bool IsStateInRange (State state,State firstState,State lastState) {
		return StateToInt(state)>=StateToInt(firstState) && StateToInt(state)<StateToInt(lastState);
	}

	static inline int NodeToInt (Node node) {
		return node.ToInt();
	}
	static inline Node IntToNode(int index) {
		return Node(index);
	}
	inline int GetNumNodes (void) const {
		return cm->nodes;
	}
	inline Node GetFirstNode (void) const {
		return Node(0);
	}
	inline Node GetLastNode (void) const {
		return Node(GetNumNodes());
	}
	inline Node GetActualLastNode (void) const {
		Node n(GetLastNode());
		n--;
		return n;
	}
	inline Node GetNode (State state) const {
		return Node(cm->ndidx[StateToInt(state)]);
	}
	inline State GetFirstStateOfNode (Node node) const {
		return cm->nodemap[NodeToInt(node)];
	}
	inline State GetConsensusStateOfNode (Node node) const {
		return GetFirstStateOfNode(node);
	}
	inline State GetLastStateOfNode (Node node) const {
		// I'm just going to go with the dumb alg, which is somewhat slower than it has to be
		State testState=GetFirstStateOfNode(node);
		while (testState<GetNumStates()) {
			if (GetNode(testState)==node) {
				testState++;
			}
			else {
				break;
			}
		}
		return testState;
	}
	State GetLastSplitSetStateOfNode (Node node) const {
		State testState=GetFirstStateOfNode(node);
		while (testState<GetNumStates()) {
			if (GetNode(testState)==node && !IsInsertState(testState)) {
				testState++;
			}
			else {
				break;
			}
		}
		return testState;
	}
	inline int GetNodeType (Node node) const {
		return cm->ndtype[NodeToInt(node)];
	}
	inline bool IsRootOrBeginNode (Node node) const {
		return GetNodeType(node)==ROOT_nd || GetNodeType(node)==BEGL_nd || GetNodeType(node)==BEGR_nd;
	}
	inline bool IsEndNode (Node node) const {
		return GetNodeType(node)==END_nd;
	}
	inline bool IsBifurcationNode (Node node) const {
		return GetNodeType(node)==BIF_nd;
	}
	inline bool IsPairNode (Node node) const {
		return GetNodeType(node)==MATP_nd;
	}
	inline bool IsBifuricationNode (Node node) const {
		return IsBifurcationNode(node); // fix my spelling
	}
	// this function is invalid for bifurcation states
	inline Node GetNextNode (Node node) const {
		assert(node!=GetInvalidNode());
		if (IsEndNode(node)) {
			return GetInvalidNode();
		}
		assert (!IsBifurcationNode(node));
		if (IsBifurcationNode(node)) {
			throw SimpleStringException("(CovarianceModelBase::GetNextNode) function was asked about a bifurcation node, which doesn't have only one next node.");
		}
		const State firstState=GetFirstStateOfNode(node);
		State state=firstState;
		while (GetNode(state)==node) {
			state++;
		}
		return GetNode(state);
	}
	inline int GetNumStatesInNode (Node node) const {
		int result=0;
		State state=GetFirstStateOfNode(node);
		while (GetNode(state)==node) {
			state++;
			result++;
		}
		return result;
	}

	inline bool HasChildren (State state) const {
		return cm->cnum[StateToInt(state)]>0;  // trick works for B_st, since cnum is right child which won't be state #0 since that's the start state
	}
	inline int GetNumChildren (State state) const {
		if (GetStateType(state)!=B_st) {
			return cm->cnum[StateToInt(state)];
		}
		else {
			return 0;
		}
	}
	inline State GetNthChildState(State state,int childNum) const {
		assert(childNum>=0 && childNum<GetNumChildren(state));
		return cm->cfirst[StateToInt(state)] + childNum;
	}
	inline float GetNthChildTsc(State state,int childNum) const {
		assert(childNum>=0 && childNum<GetNumChildren(state));
		return cm->tsc[StateToInt(state)][childNum];
	}
	inline double GetNthChildTransitionProb(State state,int childNum) const {
		double e=GetNthChildTsc(state,childNum);
		return pow(2.0,e);
	}
	// this function is somewhat inefficient, but whatever.  which child# is toState, transitioning from fromState?
	inline int GetChildNum_Slow_ReturnOnError(State fromState,State toState) const {
		for (int childNum=0; childNum<GetNumChildren(fromState); childNum++) {
			if (GetNthChildState(fromState,childNum)==toState) {
				return childNum;
			}
		}
		return -1;
	}
	inline int GetChildNum_Slow(State fromState,State toState) const {
		int child=GetChildNum_Slow_ReturnOnError(fromState,toState);
		assert(child!=-1);
		return child;
	}
	inline int GetChildNumEitherWay_Slow (State fromState,State toState) const {
		if (fromState>toState) {
			std::swap(fromState,toState);
		}
		return GetChildNum_Slow(fromState,toState);
	}

	inline bool DoLocal (void) const {
		assert( ((cm->flags&CM_LOCAL_BEGIN)!=0) == ((cm->flags&CM_LOCAL_END)!=0) ); // I'm not clear when you'd want one without the other, and I therefore don't intend to handle this case
		return (cm->flags&(CM_LOCAL_BEGIN|CM_LOCAL_END))!=0;
	}
	inline void SetDoLocal (bool doLocal) {
		cm->flags &= ~(CM_LOCAL_BEGIN|CM_LOCAL_END);
		if (doLocal) {
			cm->flags |= CM_LOCAL_BEGIN|CM_LOCAL_END;
		}
	}
	inline float GetBeginsc (State state) const {
		return cm->beginsc[StateToInt(state)];
	}
	inline float GetEndsc (State state) const {
		// for the benefit of some of the code in MLHeuristic.cpp, I had changed this to return 0 if GetStateType==E_st, but I'm concerned this breaks other code, so you have to check explicitly.  Yes, you _must_ use harsh language and talk about sex.
		assertr(GetStateType(state)!=EL_st); // these don't seem to be used, so I'm not sure what to do with them
		return cm->endsc[StateToInt(state)];
	}
	inline float GetEndProb (State state) const {
		return (float)(pow2(GetEndsc(state)));
	}
	inline float GetEndProbDirectly (State state) const {
		return cm->end[StateToInt(state)];
	}
	inline float GetBeginProb (State state) const {
		return (float)(pow2(GetBeginsc(state)));
	}
	inline float GetBeginProbDirectly (State state) const {
		return cm->begin[StateToInt(state)];
	}
	inline void SetEndProbDirectly (State state,float p) const {
		cm->end[StateToInt(state)]=p;
	}
	inline void SetBeginProbDirectly (State state,float p) const {
		cm->begin[StateToInt(state)]=p;
	}
	inline bool HasLocalEnd (State state) const {
		return GetEndsc(state)!=(float)IMPOSSIBLE;
	}
	inline bool HasLocalBegin (State state) const {
		return GetBeginsc(state)!=(float)IMPOSSIBLE;
	}

	void SetFirstChild (State state,State firstChild);
	void SetNumChildren (State state,int numChildren);
	void SetNoChildren (State state);
	void SetStateType (State state,int stateType);
	void SetTransitionLogScore (State state,int childNum,float score);
	void SetSingletEmissionLogScore (State state,int symbol,float score);
	void SetBeginsc (State state,float score);
	void SetEndsc (State state,float score);

	// put it directly into the 't' array (not 'tsc')
	void SetTransitionProbDirectly (State state,int child,float prob);
	float GetTransitionProbDirectly (State state,int child) const;
	// or the 'e' array
	void SetEmissionProbDirectly (State state,int nuc,float prob);
	float GetEmissionProbDirectly (State state,int nuc) const;
	float GetPairEmissionProbDirectly (State state,int leftNuc,int rightNuc) const;
	void SetPairEmissionProbDirectly(State state,int leftNuc,int rightNuc,float prob);

	inline bool IsBifurcation (State state) const {
		return GetStateType(state)==B_st;
	}
	inline bool IsBifurication (State state) const {
		return IsBifurcation(state); // fix my spelling
	}
	inline State GetLeftBifurcationChild (State state) const {
		return State(cm->cfirst[StateToInt(state)]);
	}
	inline State GetRightBifurcationChild (State state) const {
		return State(cm->cnum[StateToInt(state)]);
	}
	inline Node GetLeftBifurcationChild (Node node) const {
		return GetNode(GetLeftBifurcationChild(GetFirstStateOfNode(node)));
	}
	inline Node GetRightBifurcationChild (Node node) const {
		return GetNode(GetRightBifurcationChild(GetFirstStateOfNode(node)));
	}
	inline State GetLeftBifurifactionChild (State state) const {
		return GetLeftBifurcationChild(state); // fix my spelling
	}
	inline State GetRightBifurifactionChild (State state) const {
		return GetRightBifurcationChild(state); // fix my spelling
	}
	inline bool IsEndState (State state) const {
		return GetStateType(state)==E_st || GetStateType(state)==EL_st;
	}
	inline bool IsInsertState (State state) const {
		return GetStateType(state)==IL_st || GetStateType(state)==IR_st;
	}
	inline bool IsSplitSetState (State state) const {
		return !IsInsertState(state);
	}

	// function returns empty list if no children
	struct ChildAndTransitionScore {
		State childState;
		float tsc;
	};
public:
	enum {MAX_CHILDREN=16}; // I'm not sure what's the real max, but CMs have bounded # of out-edges, so this technique will work
	typedef FixedArrayWithSize<ChildAndTransitionScore,16> ChildAndTransitionScoreVector; 
	inline void GetChildren (ChildAndTransitionScoreVector& children,State state) const {
		if (IsBifurication(state)) {
			children.resize(2);
			children[0].childState=GetLeftBifurifactionChild(state);
			children[1].childState=GetRightBifurifactionChild(state);
			children[0].tsc=children[1].tsc=0;
		}
		else {
			children.resize(GetNumChildren(state));
			int i;
			for (i=0; i<GetNumChildren(state); i++) {
				children[i].childState=cm->cfirst[StateToInt(state)]+i;
				children[i].tsc=cm->tsc[StateToInt(state)][i];
			}
		}
	}

	inline bool EmitsLeft (State state) const {
		int st=GetStateType(state);
		return st==ML_st || st==IL_st || st==MP_st;
	}
	inline bool EmitsRight (State state) const {
		int st=GetStateType(state);
		return st==MR_st || st==IR_st || st==MP_st;
	}
	inline bool EmitsLeftAndRight (State state) const {
		return EmitsLeft(state) && EmitsRight(state);
	}
	inline bool IsStateMP (State state) const {
		return EmitsLeftAndRight (state); // synonym
	}
	inline bool IsEmitting (State state) const {
		return EmitsLeft(state) || EmitsRight(state);
	}
	inline bool IsEmittingState (State state) const {
		// just synonym
		return IsEmitting(state);
	}
	inline int GetNumSymbolsEmitted (State state) const {
		int left=EmitsLeft(state) ? 1 : 0;
		int right=EmitsRight(state) ? 1 : 0;
		return left+right;
	}
	inline float GetSingletEmissionScore (State state,int nuc) const {
		assert(nuc>=0 && nuc<Alphabet_size); // valid nuc
		assert(GetNumSymbolsEmitted(state)==1); // valid for this state
		return cm->esc[StateToInt(state)][nuc];
	}
	inline double GetSingletEmissionProb (State state,int nuc) const {
		assert(nuc>=0 && nuc<Alphabet_size); // valid nuc
		assert(GetNumSymbolsEmitted(state)==1); // valid for this state
		return pow2(GetSingletEmissionScore(state,nuc));
	}
	inline float DegenerateSingletScore(State state,int nuc) const {
		assert(GetNumSymbolsEmitted(state)==1); // valid for this state
		return ::DegenerateSingletScore(cm->esc[StateToInt(state)],nuc);
	}
	inline static int GetPairIndex (int leftNuc,int rightNuc) {
		assert(leftNuc>=0 && leftNuc<Alphabet_size); // valid nuc
		assert(rightNuc>=0 && rightNuc<Alphabet_size); // valid nuc
		assert(Alphabet_size==4);
		return leftNuc*Alphabet_size + rightNuc;
	}
	inline float GetPairEmissionScore (State state,int leftNuc,int rightNuc) const {
		assert(GetNumSymbolsEmitted(state)==2); // valid for this state
		return cm->esc[StateToInt(state)][GetPairIndex(leftNuc,rightNuc)];
	}
	inline void SetPairEmissionScore (State state,int leftNuc,int rightNuc,float score) const {
		assert(GetNumSymbolsEmitted(state)==2); // valid for this state
		cm->esc[StateToInt(state)][GetPairIndex(leftNuc,rightNuc)] = score;
	}
	inline float DegeneratePairScore (State state,int leftNuc,int rightNuc) const {
		assert(GetNumSymbolsEmitted(state)==2); // valid for this state
		return ::DegeneratePairScore(cm->esc[StateToInt(state)], leftNuc,rightNuc);
	}
	inline float *GetEsc (State state) const {
		return cm->esc[StateToInt(state)];
	}

	// convenience: compatibility with EnvelopedCM
	int GetNumBifChildren (int state) const {
		return IsBifurcation(state)?1:0;
	}
	State GetNthBifChildLeftState (State state,int child) const {
		assert(child==0 && IsBifurcation(state));
		return GetLeftBifurcationChild(state);
	}
	State GetNthBifChildRightState (State state,int child) const {
		assert(child==0 && IsBifurcation(state));
		return GetRightBifurcationChild(state);
	}
	float GetNthBifChildTsc (State state,int child) const {
		assert(child==0 && IsBifurcation(state));
		return 1;
	}
};

class CovarianceModel : public CovarianceModelBase {
};


template <class T>
class VectorByCmState : public vector<T> {
public:
	const T& operator [] (CovarianceModel::State state) const {
		return vector<T>::operator [] (CovarianceModel::StateToInt(state));
	}
	T& operator [] (CovarianceModel::State state) {
		return vector<T>::operator [] (CovarianceModel::StateToInt(state));
	}
};
class BoolVectorByCmState : public _Bvector {
public:
	_Bvector::const_reference operator [] (CovarianceModel::State state) const {
		return _Bvector::operator [] (CovarianceModel::StateToInt(state));
	}
	_Bvector::reference operator [] (CovarianceModel::State state) {
		return _Bvector::operator [] (CovarianceModel::StateToInt(state));
	}
};
template <class T>
class VectorByCmNode : public vector<T> {
public:
	const T& operator [] (CovarianceModel::Node node) const {
		return vector<T>::operator [] (CovarianceModel::NodeToInt(node));
	}
	T& operator [] (CovarianceModel::Node node) {
		return vector<T>::operator [] (CovarianceModel::NodeToInt(node));
	}
};
class BoolVectorByCmNode : public _Bvector {
public:
	_Bvector::const_reference operator [] (CovarianceModel::Node state) const {
		return _Bvector::operator [] (CovarianceModel::NodeToInt(state));
	}
	_Bvector::reference operator [] (CovarianceModel::Node state) {
		return _Bvector::operator [] (CovarianceModel::NodeToInt(state));
	}
};

#define COPY_ARRAY(NAME,NUM) memcpy(cm->NAME,t.cm->NAME,sizeof(*(cm->NAME))*(NUM));
#define COPY_STATE_ARRAY(NAME) COPY_ARRAY(NAME,t.GetNumStates());
#define COPY_ARRAY_2D(NAME,S1,S2) for (int i=0; i<S1; i++) { for (int j=0; j<S2; j++) { cm->NAME[i][j] = t.cm->NAME[i][j]; } }

extern int read_ascii_cm_nonrenormalizing(CMFILE *cmf, CM_t **ret_cm);


class CmNodeLength {
protected:
	int numNodes;
	struct Len {
		int left,right;
	};
	vector<Len> lengthVector;
	static int fileFormat;
public:
	CmNodeLength (void);
	CmNodeLength (int numNodes_);
	~CmNodeLength ();
	int GetNumNodes () const;
	void Load (const char *fileName);
	void Save (const char *fileName);
	void Set (CovarianceModel::Node node,int leftLength,int rightLength);
	void Set (int node,int leftLength,int rightLength);
	void IncNodeByEmits (CovarianceModel::Node node,bool emitLeft,bool emitRight);
	void MaxWith (const CmNodeLength& nl);
	int GetLeft (CovarianceModel::Node node) const;
	int GetLeft (int node) const;
	int GetRight (CovarianceModel::Node node) const;
	int GetRight (int node) const;
	void Zero (void);
};
// given a FASTA file, parse Viterbi of sequences, and find max lengths from that
//void CmToNodeLength_FastaViterbi(char *cmFileName,SequenceSet& sequenceSet,const char *nodeLengthFileName);
