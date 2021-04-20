
extern bool enableHmmCaching;

struct SolveScoresPath {
	CovarianceModel::Node cmFirstNode,cmEndingNode;
};
typedef std::list<SolveScoresPath> SolveScoresPathList;

typedef std::list<int> GlobalVarList;

// collects the info we use as we build an HMM model
struct HmmAndBuildInfo {
	InfernalHmm hmm;
	Cm2HmmStateVector cm2HmmState;
	Hmm2CmStateVector hmm2CmStateVector;
	CovarianceModel::State leftToRightPassthruState; // the old end of the CM block, which is now in the middle of the HMM

	SolveScoresPathList solveScoresPathList; // remember this info so we can solve scores easily at the end
};

typedef VectorByCmState<vector<int> > TransitionToVariableNumVector; // first dim is HMM state, then transition #, and we get the variable number for the Linear Program.  If variable number is -1, then we don't use LP to find it (and just set it 0)
struct EmissionInfo {
	InfernalHmm::State state;
	int nuc;
};
struct TransitionOrEmissionInfo {
	//bool isEmission; // else is transition
	enum Type {
		Emission,Transition,LeftwardBegin,RightwardBegin,LocalEnd
	};
	Type type;
	bool isUsed; // sanity check

	InfernalHmm::State hmmState; // for LeftwardBegin,RightwardBegin,LocalEnd
	int localEndLinkNum; // for LocalEnd
	InfernalHmm::EdgeInfo edgeInfo; // for Transition
	EmissionInfo emissionInfo; // for Emission
};
typedef vector<TransitionOrEmissionInfo> TransitionOrEmissionInfoVector;

extern void GetGlobalVarsFromInfernalHmm (vector<double>& globalVars,const InfernalHmm& infernalHmm,const TransitionOrEmissionInfoVector& transitionOrEmissionInfoVector);
extern void SetGlobalVarsIntoInfernalHmm (InfernalHmm& infernalHmm,const vector<double>& globalVars,const TransitionOrEmissionInfoVector& transitionOrEmissionInfoVector);

struct ScoreVariablesInfo {
	TransitionToVariableNumVector transitionToVariableNumVector,emissionToVariableNumVector,localEndToVariableNumberVector;
	VectorByCmState<int> leftwardBeginToVariableNumberVector,rightwardBeginToVariableNumberVector;
	TransitionOrEmissionInfoVector globalVariableToTransitionOrEmissionVector;
	int numVariables;
};

struct InequalitiesAndLocalVariables {
	vector<int> globalToLocalVariables;
	vector<int> localToGlobalVariables;
	int numLocalVariables;
	InequalityList inequalityList;
};
typedef VectorByCmNode<InequalitiesAndLocalVariables> InequalitiesAndLocalVariablesByCmNode;

struct ExtraCm2HmmInfo {
	// outputs of Cm2Hmm
	InequalitiesAndLocalVariablesByCmNode inequalitiesAndLocalVariables;
	ScoreVariablesInfo scoreVariablesInfo;

	// inputs
	bool actuallySolveScores; // set this to false to save time when we only want the inequalities, not the actual solutions

	ExtraCm2HmmInfo (void) {
		actuallySolveScores=true;
	}
};

extern void AddInequalitiesInTermsOfGlobalVars(InequalityList& inequalityList,const ExtraCm2HmmInfo& extraCm2HmmInfo,const CovarianceModel::Node cmNode);
extern void AddCrossProductOfInequalityLists(InequalityList& inequalityList,InequalityListList& inequalityListList);
extern void ConvertInequalityListVarNums (InequalityList& inequalityList,const vector<int>& varNumMapping);

extern void GetGlobalVariableNumsForNode (std::list<int>& globalVarNums,const CovarianceModel& cm,InfernalHmm& infernalHmm,const CovarianceModel::Node cmNode,const ExtraCm2HmmInfo& extraCm2HmmInfo);
extern void GetLocalVariablesValueForNode(vector<float>& localVariablesToValue,const CovarianceModel& cm,InfernalHmm& infernalHmm,const CovarianceModel::Node cmNode,const ExtraCm2HmmInfo& extraCm2HmmInfo);
extern void SetLocalVariablesValueForNode(vector<float>& localVariablesToValue,const CovarianceModel& cm,InfernalHmm& infernalHmm,const CovarianceModel::Node cmNode,const ExtraCm2HmmInfo& extraCm2HmmInfo);


class TemporarilyModifyInequality {
protected:
	Inequality& inequalitySoFar;
	int numVariablesAdded;
	float startingScore;
	double startingWeight;
	float starting_sumOfConstantsInHmm;
	size_t starting_hmmInsertStatesInPath_size;
	int startingNucEmitCount[MAXABET];
public:
	TemporarilyModifyInequality (Inequality& _inequalitySoFar);
	~TemporarilyModifyInequality ();

	void AddScore (float addToScore);
	void AddVariable (int globalVariableNum,vector<int>& globalToLocalVariables,int& numlocalVariables);
	void MultiplyWeight (double mult);
	void PushInsertState (InfernalHmm::State insertState);
};

// this looks like it's only relevant to the "EMish" algorithm for optimizing HMMs, which is deprecated because it's vastly inferior to the infinite-length forward algorithm
struct WeightedInequalitiesInfo {
	const TransitionCounter *transitionCounter;
	double pseudocount; // should be 1, but maybe 0.5 makes sense.  Setting this to 0 is very dangerous, since then you can get divide by 0 for states that were never visited.  While in the ensuing debugging, you might miss your bus, and be stranded at school forever.  (Not really; it's only a 30-minute walk.)
};

// HMM committees (a series of HMMs to apply) are deprecated, since it doesn't help much when the HMMs are of the same type (i.e. all compact or all expanded) and are optimized with the inf-len forward alg.  Anyway, the same functionality can be specified on the command line.
class HmmCommittee {
protected:
	std::string description;
	int recommendedCommitteeSize; // is computed by evaluate_everyone, otherwise just unknown (0)
public:
	int GetRecommendedCommitteeSize (void) const; // will 'Die' if unknown

	typedef std::list<InfernalHmm> HmmList;
	HmmList hmmList;

	HmmCommittee (void);
	~HmmCommittee ();

	bool /*success*/ LoadInBinary (const char *fileName,int effectiveCommitteeSize);
	bool SaveInBinary (const char *fileName,const char *additionalDescription);

	void DumpMemberBuildInfo (FILE *out);

	// information stored on a linear program, so when we want to solve it again, we can generate a new solution.  For now, I'm just planning on finding all solutions at once, and then picking from them.  Another simplification, is that I assume the only distinction between solutions is which constraints they satisfy perfectly (i.e. which slack variables are 0).  In practice, there seem to be multiple solutions in this sense, but clearly you should be able to generate intermediates.
	struct LinearProgramStoredInfo {

		unsigned int numTimesGeneratedSolution; // our pattern is to first go thru all solutions, and then to pick randomly
		vector<vector<float> > localVariablesToValuePerSolution; // first dimension is solution #, 2nd is the local variables' values in that solution
		vector<float> avgInflationPerSolution;

		/*
		vector<double> sumOfSlack; // each constraint has a slack variable.  This is the sum of the slack variables in the optimal solutions generated so far
		_Bvector slackHasBeenZero; // each constraint may have been solved perfectly in one of the solutions, if so slack has zero is set to true

		double optimalValueOfObjectiveFunc; // when we come up with alternate solutions, we want them to all be just as optimal, i.e. have the same objective func value

		bool isFirstSolutionFound;
		LinearProgramStoredInfo (void) {
			isFirstSolutionFound=false;
		}
		*/
	};
	// allows me to add ideas later.  Usually passed around as a pointer, where NULL implies no need to do anything about committee.  Stores input parameters, and work-in-progress
	struct CommitteeBuildInfo {
		VectorByCmNode<LinearProgramStoredInfo> cmStartNodeToLinearProgramInfo;
		int currCommitteeMemberNum; // first is 0, then 1, ...
	};
};

extern void Cm2Hmm (InfernalHmm& hmm,Cm2Hmm_HmmBuildType hmmBuildType,const CovarianceModel& sourceCM,const char *cmFileName,const std::string& programParams,bool forceCreate=false); // cmFileName just used for caching the translation to an HMM

extern void Cm2Hmm (HmmCommittee& hmmCommittee,const CovarianceModel& sourceCM,const char *cmFileName,int committeeSize,int effectiveCommitteeSize);

// slightly lower-level function
void Cm2Hmm_WithWeighting_NoCaching (InfernalHmm& hmm,Cm2Hmm_HmmBuildType hmmBuildType,const CovarianceModel& sourceCM,const char *cmFileName,WeightedInequalitiesInfo *weightedInequalitiesInfo,ExtraCm2HmmInfo *extraCm2HmmInfo=NULL,const int nodesToSpanWhileSolvingScores=1);

void InitLocalData (InfernalHmm& hmm,const CovarianceModel& sourceCM,bool setupLocalBegin,bool setupLocalEnd);


// relatively generic Linear Programming stuff, here for convenience
bool /* solved to completion */ SolveLP (vector<double>& optimalVars,const vector<double>& coefficients,const InequalityList& inequalityList,int maxSolverIters=INT_MAX);
