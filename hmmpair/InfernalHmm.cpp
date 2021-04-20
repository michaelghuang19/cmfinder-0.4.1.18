#include "hmmpair.h"

///////////////////
// InfernalHmm

InfernalHmm::InfernalHmm (void)
{
	loadedFileFormat=GetCurrFileFormatNum();
	hmmBuildType=HmmBuildType_Original; // assume this for now
}
std::string InfernalHmm::GetBuildDescription (void) const
{
	return "GetBuildDescription not implemented";
}
void InfernalHmm::SetHmmBuildType (Cm2Hmm_HmmBuildType hmmBuildType_)
{
	hmmBuildType=hmmBuildType_;
}
void InfernalHmm::SetFromCmFileName (const std::string& _fromCmFileName)
{
	fromCmFileName=_fromCmFileName;
}
std::string InfernalHmm::GetLoadedHmmFileName (void) const
{
	return loadedHmmFileName;
}
void InfernalHmm::AddBuildDescription (const std::string addedBuildDescription)
{
	fullBuildDescription += ";";
	fullBuildDescription += addedBuildDescription;
}
int InfernalHmm::GetFileFormat (void) const
{
	return loadedFileFormat;
}
void InfernalHmm::Init (int numStates)
{
	CovarianceModelBase::Init(numStates);

	EndscLinksToLeft emptyLinksToLeft;
	emptyLinksToLeft.resize(0);
	endscLinksToLeftVector.assign(numStates,emptyLinksToLeft);

	otherStateInfoVector.resize(numStates);
	hmm2CmStateVector.resize(numStates);
	cm2HmmStateVector.clear();
	nonSavedInfoVector.clear();
}
void InfernalHmm::BuildNonSavedInfoIfNecessary (void)
{
	if (nonSavedInfoVector.empty()) {

		// must build
		nonSavedInfoVector.resize(GetNumStates());

		State parentState;
		for (parentState=GetFirstState(); parentState!=GetLastState(); parentState++) {
			for (int childNum=0; childNum<GetNumChildren(parentState); childNum++) {

				ParentAndMyChildNum parentAndMyChildNum;
				parentAndMyChildNum.parentState=parentState;
				parentAndMyChildNum.myChildNum=childNum;
				State myState=GetNthChildState(parentState,childNum);

				nonSavedInfoVector[myState].parentAndMyChildNumVector.push_back(parentAndMyChildNum);
			}
		}

		State state;
		for (state=GetFirstState(); state!=GetLastState(); state++) {
			nonSavedInfoVector[state].childNumOfSelfLoop=-1; // assume nothing
			for (int childNum=0; childNum<GetNumChildren(state); childNum++) {
				State childState=GetNthChildState(state,childNum);
				if (state==childState) {
					nonSavedInfoVector[state].childNumOfSelfLoop=childNum;
				}
			}
		}
	}
}
int InfernalHmm::GetChildNumOfSelfLoop (State state) const
{
	assert(!nonSavedInfoVector.empty()); // didn't call BuildNonSavedInfoIfNecessary
	return nonSavedInfoVector[state].childNumOfSelfLoop;
}
void InfernalHmm::SetHmm2CmState (const Hmm2CmStateVector& _hmm2CmStateVector)
{
	hmm2CmStateVector=_hmm2CmStateVector;
}
void InfernalHmm::SetCm2HmmState (const Cm2HmmStateVector& _cm2HmmStateVector)
{
	cm2HmmStateVector=_cm2HmmStateVector;
}
void InfernalHmm::CopyFrom (const InfernalHmm& t)
{
	CovarianceModelBase::CopyFrom(t);
	endscLinksToLeftVector=t.endscLinksToLeftVector;
	otherStateInfoVector=t.otherStateInfoVector;
	hmm2CmStateVector=t.hmm2CmStateVector;
	cm2HmmStateVector=t.cm2HmmStateVector;
	nonSavedInfoVector=t.nonSavedInfoVector;
	loadedFileFormat=t.loadedFileFormat;
	fullBuildDescription=t.fullBuildDescription;
	fromCmFileName=t.fromCmFileName;
	hmmBuildType=t.hmmBuildType;
	loadedHmmFileName=t.loadedHmmFileName;
	localEndSelfLoopScore=t.localEndSelfLoopScore;
}
void InfernalHmm::MirrorFromWithHackedExtraInfo(CM_t *cm)
{
	MirrorFrom(cm);
	endscLinksToLeftVector.resize(GetNumStates());
	otherStateInfoVector.resize(GetNumStates());
	hmm2CmStateVector.resize(GetNumStates());
	cm2HmmStateVector.clear();
}
void InfernalHmm::ClobberIR (void)
{
	State state=GetFirstState();
	state++;
	state++;
	if (GetStateType(state)!=IR_st) {
		throw SimpleStringException("InfernalHmm::ClobberIR: expected IR to be 3rd state.");
	}

	// change IR to a D state, with the same transitions.  I think this is a reasonable hack, and whatever we do the IR state will have low probability and there's only 1 IR state, so it won't make that much of a difference.
	SetStateType(state,D_st);

	// now verify there's nothing else right-emitting.
	for (state=GetFirstState(); state!=GetLastState(); state++) {
		if (EmitsRight(state)) {
			throw SimpleStringException("InfernalHmm::ClobberIR: found unexpected right-emitting state.");
		}
	}
}
void InfernalHmm::CopyFromWithEscHack(const InfernalHmm& t)
{
	CopyFrom(t);
	
	for (State state=GetFirstState(); state!=GetLastState(); state++) {
		if (IsEmitting(state)) {
			for (int nuc=0; nuc<Alphabet_size; nuc++) {
				float esc=GetSingletEmissionScore(state,nuc);
				esc += (float)2.0;
				SetSingletEmissionLogScore (state,nuc,esc);
			}
		}
	}
}
void InfernalHmm::SetLeftState (State state,bool isLeftState)
{
	otherStateInfoVector[state].isLeftState=isLeftState;
}
void InfernalHmm::SetRightState (State state,bool isRightState)
{
	otherStateInfoVector[state].isRightState=isRightState;
}
void InfernalHmm::SetLeftwardBeginsc (State state,float sc)
{
	otherStateInfoVector[state].leftwardBeginsc=sc;
}
void InfernalHmm::SetRightwardBeginsc (State state,float sc)
{
	otherStateInfoVector[state].rightwardBeginsc=sc;
}
void InfernalHmm::DumpCsv (const char *fileName,const CovarianceModel& cm)
{
	FILE *csv=ThrowingFopen(fileName,"wt");
	DumpCsv(csv,cm);
	fclose(csv);
}
void InfernalHmm::DumpCsv (FILE *csv,const CovarianceModel& cm)
{
	fprintf(csv,"Local,%d\n",DoLocal()?1:0);
	fprintf(csv,"# states,%d\n",GetNumStates());
	fprintf(csv,"\n");
	fprintf(csv,"state#,state type,,");
	if (DoLocal()) {
		fprintf(csv,"left beginsc,rightbeginsc,end #1 left,end #1 score,end #2 left,end #2 score,,");
	}
	fprintf(csv,"emit A,emit C,emit G,emit U,,first child,transition1,transition2...\n");
	for (State state=GetFirstState(); state!=GetLastState(); state++) {
		fprintf(csv,"%d,%s,,",StateToInt(state),GetStateTypeName(state));
		if (DoLocal()) {
			fprintf(csv,"%g,%g,",GetLeftwardBeginsc(state),GetRightwardBeginsc(state));
			for (int link=0; link<2; link++) {
				if (GetNumEndscLinksToLeft(state)>link) {
					fprintf(csv,"%d,%g,",StateToInt(GetEndscLinkToLeft_State(state,link)),GetEndscLinkToLeft_Endsc(state,link));
				}
				else {
					fprintf(csv,",,");
				}
			}
			fprintf(csv,",");
		}
		for (int nuc=0; nuc<Alphabet_size; nuc++) {
			if (IsEmittingState(state)) {
				fprintf(csv,"%g,",GetSingletEmissionScore(state,nuc));
			}
			else {
				fprintf(csv,",");
			}
		}

		fprintf(csv,",");
		if (GetNumChildren(state)==0) {
		  // nothing to emit
		}
		else {
			fprintf(csv,"%d,",StateToInt(GetNthChildState(state,0)));
			for (int child=0; child<GetNumChildren(state); child++) {
				fprintf(csv,"%g,",GetNthChildTsc(state,child));
			}
		}

		fprintf(csv,"\n");
	}

	fprintf(csv,"\n\nCM to HMM mappings\n");
	if (cm2HmmStateVector.empty()) {
		fprintf(csv,"MAPPINGS HAVE NOT BE CALCULATED & RECORDED YET.\n");
	}
	else {
		fprintf(csv,"\nCM state , HMM left state , HMM right state\n");
		for (CovarianceModel::State cmState=cm.GetFirstState(); cmState!=cm.GetLastState(); cmState++) {
			State hmmLeftState=cm2HmmStateVector[cmState].hmmLeftState;
			State hmmRightState=cm2HmmStateVector[cmState].hmmRightState;
			char leftStr[32],rightStr[32];
			if (hmmLeftState==GetInvalidState()) {
			  // nothing to emit
			}
			else {
				sprintf(leftStr,"%d",StateToInt(hmmLeftState));
			}
			if (hmmRightState==GetInvalidState()) {
			  // nothing to emit
			}
			else {
				sprintf(rightStr,"%d",StateToInt(hmmRightState));
			}
			fprintf(csv,"%d,%s,%s\n",CovarianceModel::StateToInt(cmState),leftStr,rightStr);
		}
		fprintf(csv,"\nCM node,{HMM state set}\n");
		for (CovarianceModel::Node cmNode=cm.GetFirstNode(); cmNode!=cm.GetLastNode(); cmNode++) {

			fprintf(csv,"%d,",CovarianceModel::NodeToInt(cmNode));

			StateFromNodeList hmmStateList;
			GetHmmStatesByCmNode (hmmStateList,cm,cmNode);
			for (StateFromNodeList::iterator i=hmmStateList.begin(); i!=hmmStateList.end(); i++) {
				fprintf(csv,"%d,",StateToInt(i->state));
			}
			fprintf(csv,"\n");
		}
	}
}
void InfernalHmm::AllocHmmData (void)
{
	endscLinksToLeftVector.resize(GetNumStates());
	otherStateInfoVector.resize(GetNumStates());
}
int InfernalHmm::GetCurrFileFormatNum (void)
{
	return 5;
}
void InfernalHmm::LoadInDeprecatedBinaryFormat (FILE *file)
{
	char firstFourChars[5];
	fread(firstFourChars,4,1,file);
	firstFourChars[4]=0;
	// old binary format -- and BTW those 4 bytes we read we the # states in the original format
	int numStates=*(int *)(firstFourChars);

	if (numStates==0) {

		// new format
		int format;
		fread(&format,sizeof(format),1,file);
		loadedFileFormat=format;

		//fprintf(stderr,"Reading InfernalHmm file format #%d\n",format);
		if (format<0 || format>=GetCurrFileFormatNum()+1) {
			throw SimpleStringException("Unknown file format for InfernalHmm");
		}

		if (format>=2) {

			// read cmFileName and fullBuildDescription
			size_t s;
			fread(&s,sizeof(s),1,file);
			char *p;
			p=new char[s+1];
			fread(p,1,s,file);
			p[s]=0;
			fromCmFileName=p;
			delete [] p;

			fread(&s,sizeof(s),1,file);
			p=new char[s+1];
			fread(p,1,s,file);
			p[s]=0;
			fullBuildDescription=p;
			delete [] p;
		}

		if (format>=3) {
			fread(&hmmBuildType,sizeof(hmmBuildType),1,file);
		}
		else {
			hmmBuildType=HmmBuildType_Original;
		}

		// convenient to know # of states before format-specific stuff
		fread(&numStates,sizeof(numStates),1,file);
		Init(numStates);

		int int_doLocal;
		fread(&int_doLocal,sizeof(int_doLocal),1,file);
		SetDoLocal(int_doLocal!=0);

		if (DoLocal()) {
			fread(&(endscLinksToLeftVector.front()),sizeof(EndscLinksToLeft),numStates,file);
		}

		fread(&(otherStateInfoVector.front()),sizeof(OtherStateInfo),numStates,file);

		hmm2CmStateVector.resize(numStates);
		for (State state=IntToState(0); state<IntToState(numStates); state++) {
			int n;
			fread(&n,sizeof(n),1,file);
			for (int i=0; i<n; i++) {
				int i_state;
				fread(&i_state,sizeof(i_state),1,file);
				hmm2CmStateVector[state].push_back(IntToState(i_state));
			}
		}

		if (format>=1) {
			size_t numCmStates;
			fread(&numCmStates,sizeof(numCmStates),1,file);
			cm2HmmStateVector.resize(numCmStates);

			fread(&(cm2HmmStateVector.front()),sizeof(Cm2HmmState),numCmStates,file);
		}
	}
	else {
		// original format
		Init(numStates);
	}

	fread(cm->stid,sizeof(cm->stid[0]),numStates+1,file);
	fread(cm->sttype,sizeof(cm->sttype[0]),numStates+1,file);
	fread(cm->cfirst,sizeof(cm->cfirst[0]),numStates,file);
	fread(cm->cnum,sizeof(cm->cnum[0]),numStates,file);

	for (int state=0; state<numStates; state++) {
		fread(cm->esc[state],sizeof(cm->esc[0][0]),Alphabet_size*Alphabet_size,file);
		fread(cm->tsc[state],sizeof(cm->tsc[0][0]),MAXCONNECT,file);
	}
}
bool InfernalHmm::LoadInBinary (FILE *file)
{
	// new text format

	fscanf(file,"ProfileHMMFilter\n");

	if (fscanf(file,"FileFormatNum:\t%d\n",&loadedFileFormat)!=1) {
		throw SimpleStringException("Input profile HMM missing 'FileFormatNum'");
	}

	int numStates;
	if (fscanf(file,"states:\t%d\n",&numStates)!=1) {
		throw SimpleStringException("Input profile HMM missing 'states'");
	}
	Init(numStates);

	char cmFileNameTemp[4096];
	if (fscanf(file,"cmFileName:\t%s\n",cmFileNameTemp)!=1) {
		throw SimpleStringException("Input profile HMM missing 'cmFileName'");
	}
	fromCmFileName=cmFileNameTemp;

	fullBuildDescription.clear();
	fscanf(file,"hmmBuildCommands: ");
	while (true) {
		int ch=fgetc(file);
		if (ch=='\r' || ch=='\n' || ch==EOF) {
			break;
		}
		if (ch=='\t') {
			ch='\n'; // I'm finicky about restoring the '\n's because by making things exactly the same, it helps me to show that this text format works like the binary format by loading a binary-format file, saving in text format, loading that, saving in binary, and verifying that the file is the same as the original.  Other than that, it doesn't matter.
		}
		fullBuildDescription += ch;
	}

	char hmmTypeText[256];
	if (fscanf(file,"hmmType:\t%s\n",hmmTypeText)!=1) {
		throw SimpleStringException("Input profile HMM missing 'hmmType'");
	}
	hmmBuildType=GetHmmBuildTypeByText(hmmTypeText);

	int doLocalInt;
	if (fscanf(file,"isLocal:\t%d\n",&doLocalInt)!=1) {
		throw SimpleStringException("Input profile HMM missing 'isLocal'");
	}
	SetDoLocal(doLocalInt!=0);

	if (DoLocal()) {
		if (loadedFileFormat>=5) {
			if (fscanf(file,"localEndSelfLoopScore:\t%g\n",&localEndSelfLoopScore)!=1) {
				throw SimpleStringException("Input profile HMM missing 'localEndSelfLoopScore'");
			}
		}
		else {
			localEndSelfLoopScore=0;
		}
	}

	for (State state=GetFirstState(); state!=GetLastState(); state++) {
		int i_state=StateToInt(state);

		char sttypeText[256];
		int stid;
		if (fscanf(file,"%d\t%s\t%d\t%d\t",&(stid),sttypeText,&(cm->cfirst[i_state]),&(cm->cnum[i_state]))!=4) {
			throw SimpleStringException("Input profile HMM state #%d didn't have the expected first 4 fields (unique id,state type,first child,# children)",i_state);
		}
		cm->stid[i_state]=(char)stid;
		if (strcmp(sttypeText,"PASSTHRU")==0) {
			cm->sttype[i_state]=PASSTHRU_st;
		}
		else {
			cm->sttype[i_state]=StateCode(sttypeText);
		}

		int isLeftStateInt,isRightStateInt;
		if (fscanf(file,"%d\t%d\t%g\t%g\t",&isLeftStateInt,&isRightStateInt,&(otherStateInfoVector[state].leftwardBeginsc),&(otherStateInfoVector[state].rightwardBeginsc))!=4) {
			throw SimpleStringException("Input profile HMM state #%d didn't have the expected vestigal isLocal fields (isLeft,isRight,leftwardBeginsc,rightwardBegins)",i_state);
		}
		otherStateInfoVector[state].isLeftState=isLeftStateInt!=0;
		otherStateInfoVector[state].isRightState=isRightStateInt!=0;

		size_t numRelatedCmStates;
		if (fscanf(file,"%lu\t",&numRelatedCmStates)!=1) {
			throw SimpleStringException("Input profile HMM state #%d didn't have the specified of the # of related CM states",i_state);
		}
		for (size_t i=0; i<numRelatedCmStates; i++) {
			int cmState;
			if (fscanf(file,"%d\t",&cmState)!=1) {
				printf("%d,%lu,%lu\n",i_state,numRelatedCmStates,i);
				throw SimpleStringException("Input profile HMM state #%d supposedly had %lu related CM states, but they weren't all there (I could only read %lu)",i_state,numRelatedCmStates,i);
			}
			hmm2CmStateVector[state].push_back(CovarianceModel::IntToState(cmState));
		}

		if (IsEmittingState(state)) {
			for (int nuc=0; nuc<Alphabet_size; nuc++) {
				if (fscanf(file,"%g\t",&(cm->esc[i_state][nuc]))!=1) {
					throw SimpleStringException("Input profile HMM state #%d was supposedly emitting, but I couldn't read all its 4 emit scores",i_state);
				}
			}
		}

		// tsc
		for (int child=0; child<GetNumChildren(state); child++) {
			float tsc;
			if (fscanf(file,"%g\t",&tsc)!=1) {;
				throw SimpleStringException("Input profile HMM state #%d supposedly has %d children, but I couldn't reaed the transition score for its %dth child",i_state,GetNumChildren(state),child);
			}
			cm->tsc[i_state][child]=tsc;
		}

		if (DoLocal()) {
			// other info for local
			int numLinks;
			if (fscanf(file,"%d\t",&numLinks)!=1) {
				throw SimpleStringException("Input profile HMM state #%d didn't have the #of local end links to left, even though it's supposedly a local alignment profile HMM",i_state);
			}
			SetNumEndscLinksToLeft(state,numLinks);
			for (int link=0; link<numLinks; link++) {
				int i_linkState;
				float endsc;
				fscanf(file,"%d\t%g\t",&i_linkState,&endsc);
				SetEndscLinkToLeft_State(state,link,InfernalHmm::IntToState(i_linkState));
				SetEndscLinkToLeft_Endsc(state,link,endsc);
			}
		}

		fscanf(file,"\n");
	}

	size_t numCmStates;
	if (fscanf(file,"%lu\t",&numCmStates)!=1) {
		throw SimpleStringException("Input profile HMM didn't specify the # of CM states in the underlying CM");
	}
	cm2HmmStateVector.resize(numCmStates);
	for (size_t i=0; i<numCmStates; i++) {
		CovarianceModel::State cmState=CovarianceModel::IntToState((int)i);
		int l,r;
		fscanf(file,"%d\t%d\n",&l,&r);
		cm2HmmStateVector[cmState].hmmLeftState=InfernalHmm::IntToState(l);
		cm2HmmStateVector[cmState].hmmRightState=InfernalHmm::IntToState(r);
	}

	return true;
}
bool InfernalHmm::LoadInBinary (const char *hmmFileName)
{
	loadedHmmFileName=hmmFileName;

	// WARNING: only does what HMM needs -- no nodes, null model

	FILE *file=fopen(hmmFileName,"rt");
	if (file==NULL) {
		return false;
	}

	bool result=LoadInBinary(file);
	
	fclose(file);

	return result;
}
void InfernalHmm::SaveInDeprecatedBinaryFormat (const char *fileName)
{
	FILE *file=ThrowingFopen(fileName,"wb");
	SaveInDeprecatedBinaryFormat (file);
	fclose(file);
}
void InfernalHmm::SaveInDeprecatedBinaryFormat (FILE *file)
{
	int zero=0;
	fwrite(&zero,sizeof(zero),1,file);

	int format=3;
	fwrite(&format,sizeof(format),1,file);

	size_t s;
	s=fromCmFileName.size();
	fwrite(&s,sizeof(s),1,file);
	fwrite(fromCmFileName.c_str(),1,s,file);
	s=fullBuildDescription.size();
	fwrite(&s,sizeof(s),1,file);
	fwrite(fullBuildDescription.c_str(),1,s,file);

	fwrite(&hmmBuildType,sizeof(hmmBuildType),1,file);

	fwrite(&(cm->M),sizeof(cm->M),1,file);

	int int_doLocal=DoLocal()?1:0;
	fwrite(&int_doLocal,sizeof(int_doLocal),1,file);
	if (DoLocal()) {
		fwrite(&(endscLinksToLeftVector.front()),sizeof(EndscLinksToLeft),GetNumStates(),file);
	}

	fwrite(&(otherStateInfoVector.front()),sizeof(OtherStateInfo),GetNumStates(),file);
	for (State state=IntToState(0); state<IntToState(GetNumStates()); state++) {
		size_t n=hmm2CmStateVector[state].size();
		fwrite(&n,sizeof(n),1,file);
		for (std::list<State>::const_iterator i=hmm2CmStateVector[state].begin(); i!=hmm2CmStateVector[state].end(); i++) {
			int i_state=StateToInt(*i);
			fwrite(&i_state,sizeof(i_state),1,file);
		}
	}

	size_t numCmStates=cm2HmmStateVector.size();
	fwrite(&numCmStates,sizeof(numCmStates),1,file);
	fwrite(&(cm2HmmStateVector.front()),sizeof(Cm2HmmState),numCmStates,file);

	fwrite(cm->stid,sizeof(cm->stid[0]),GetNumStates()+1,file);
	fwrite(cm->sttype,sizeof(cm->sttype[0]),GetNumStates()+1,file);
	fwrite(cm->cfirst,sizeof(cm->cfirst[0]),GetNumStates(),file);
	fwrite(cm->cnum,sizeof(cm->cnum[0]),GetNumStates(),file);

	for (int state=0; state<GetNumStates(); state++) {
		fwrite(cm->esc[state],sizeof(cm->esc[0][0]),Alphabet_size*Alphabet_size,file);
		fwrite(cm->tsc[state],sizeof(cm->tsc[0][0]),MAXCONNECT,file);
	}
}
void InfernalHmm::SaveInBinary (FILE *file)
{
	fprintf(file,"ProfileHMMFilter\n");
	fprintf(file,"FileFormatNum:\t%d\n",GetCurrFileFormatNum());
	fprintf(file,"states:\t%d\n",GetNumStates());
	fprintf(file,"cmFileName:\t%s\n",fromCmFileName.c_str());
	fprintf(file,"hmmBuildCommands:\t%s\n",GetBuildDescription().c_str());
	fprintf(file,"hmmType:\t%s\n",GetHmmBuildTypeDescription(hmmBuildType));
	fprintf(file,"isLocal:\t%d\n",DoLocal()?1:0);
	if (DoLocal()) {
		fprintf(file,"localEndSelfLoopScore:\t%g\n",GetLocalEndSelfLoopScore());
	}

	for (State state=GetFirstState(); state!=GetLastState(); state++) {

		int i_state=StateToInt(state);

		// basic info
		const char *sttypeText;
		if (cm->sttype[i_state]==PASSTHRU_st) {
			sttypeText="PASSTHRU";
		}
		else {
			sttypeText=Statetype(cm->sttype[i_state]);
		}
		fprintf(file,"%d\t%s\t%d\t%d\t",0,sttypeText,cm->cfirst[i_state],cm->cnum[i_state]);

		// vestigal info for local alignments
		fprintf(file,"%d\t%d\t%.9g\t%.9g\t",otherStateInfoVector[state].isLeftState?1:0,otherStateInfoVector[state].isRightState?1:0,otherStateInfoVector[state].leftwardBeginsc,otherStateInfoVector[state].rightwardBeginsc);

		// hmm2cm state mapping (list of CM states that relate to this HMM state)
		size_t numRelatedCmStates=hmm2CmStateVector[state].size();
		fprintf(file,"%lu\t",numRelatedCmStates);
		for (std::list<State>::const_iterator i=hmm2CmStateVector[state].begin(); i!=hmm2CmStateVector[state].end(); i++) {
			int i_cmState=StateToInt(*i);
			fprintf(file,"%d\t",i_cmState);
		}

		// esc
		if (IsEmittingState(state)) {
			for (int nuc=0; nuc<Alphabet_size; nuc++) {
				fprintf(file,"%.9g\t",GetSingletEmissionScore(state,nuc));
			}
		}
		// tsc
		for (int child=0; child<GetNumChildren(state); child++) {
			fprintf(file,"%.9g\t",GetNthChildTsc(state,child));
		}

		if (DoLocal()) {
			// other info for local
			fprintf(file,"%d\t",GetNumEndscLinksToLeft(state));
			for (int link=0; link<GetNumEndscLinksToLeft(state); link++) {
				fprintf(file,"%d\t%.9g\t",StateToInt(GetEndscLinkToLeft_State(state,link)),GetEndscLinkToLeft_Endsc(state,link));
			}
		}

		fprintf(file,"\n");
	}

	// mapping of CM states to HMM left & right states
	fprintf(file,"%lu\n",cm2HmmStateVector.size());
	for (size_t i=0; i<cm2HmmStateVector.size(); i++) {
		CovarianceModel::State cmState=CovarianceModel::IntToState((int)i);
		fprintf(file,"%d\t%d\n",InfernalHmm::StateToInt(cm2HmmStateVector[cmState].hmmLeftState),InfernalHmm::StateToInt(cm2HmmStateVector[cmState].hmmRightState));
	}
}
void InfernalHmm::SaveInBinary (const char *cmFileName)
{
	FILE *file=fopen(cmFileName,"wt");
	if (file==NULL) {
		throw SimpleStringException("Could not open '%s' for writing.  %s:%d",cmFileName,__FILE__,__LINE__);
	}

	SaveInBinary(file);

	fclose(file);
}
void InfernalHmm::SetLocalEndSelfLoopScore (float sc)
{
	localEndSelfLoopScore=sc;
}
void InfernalHmm::CopyReverseOf(const InfernalHmm& t)
{
	// WARNING: only for HMMs
	cm=CreateCM(t.GetNumNodes(),t.GetNumStates());

	localEndSelfLoopScore=t.GetLocalEndSelfLoopScore();

	cm->name=(char *)MallocOrDie(strlen(t.GetName())+1);
	strcpy(cm->name,t.GetName());

	COPY_ARRAY_2D(esc,t.GetNumStates(),Alphabet_size*Alphabet_size);

	COPY_STATE_ARRAY(beginsc);
	endscLinksToLeftVector=t.endscLinksToLeftVector;
	otherStateInfoVector=t.otherStateInfoVector;
	hmm2CmStateVector=t.hmm2CmStateVector;
	cm2HmmStateVector=t.cm2HmmStateVector;
	cm->flags=t.cm->flags;

	// initially no children
	State myState;
	for (myState=GetFirstState(); myState!=GetLastState(); myState++) {
		SetStateType(myState,t.GetStateType(myState));
		SetNoChildren(myState);
	}

	State srcState;
	for (srcState=t.GetFirstState(); srcState!=t.GetLastState(); srcState++) {
		for (int childNum=0; childNum<t.GetNumChildren(srcState); childNum++) {
			State srcToState=t.GetNthChildState(srcState,childNum);
			float tsc=t.GetNthChildTsc(srcState,childNum);
			if (GetNumChildren(srcToState)==0) {
				SetFirstChild(srcToState,srcState);
				SetNumChildren(srcToState,1);
				SetTransitionLogScore(srcToState,0,tsc);
			}
			else {
				int int_srcToState=StateToInt(srcToState);
#ifdef _DEBUG
				int int_srcState=cm->cfirst[int_srcToState] + cm->cnum[int_srcToState];
				assert(srcState==IntToState(int_srcState));
#endif
				cm->tsc[int_srcToState][cm->cnum[int_srcToState]]=tsc;
				cm->cnum[int_srcToState]++;
			}
		}
	}
}
int InfernalHmm::GetNumStatesForCmNode (const CovarianceModel& cm,CovarianceModel::Node cmNode) const
{
	// the cheezy, slow implementation
	StateFromNodeList hmmStates;
	GetHmmStatesByCmNode(hmmStates,cm,cmNode);
	return (int)(hmmStates.size());
}
void InfernalHmm::GetHmmStatesByCmNode (StateFromNodeList& get_hmmStates,const CovarianceModel& cm,CovarianceModel::Node cmNode) const
{
	get_hmmStates.clear();

	BoolVectorByCmState hmmStateUsed;
	hmmStateUsed.assign(GetNumStates(),false);

	CovarianceModel::State cmState;
	for (cmState=cm.GetFirstStateOfNode(cmNode); cmState!=cm.GetLastStateOfNode(cmNode); cmState++) {

		InfernalHmm::State hmmLeftState=cm2HmmStateVector[cmState].hmmLeftState;
		if (hmmLeftState!=GetInvalidState()) {
			if (!hmmStateUsed[hmmLeftState]) {
				hmmStateUsed[hmmLeftState]=true;
				StateFromNode add;
				add.state=hmmLeftState;
				add.isRight=false;
				add.isNormal=!cm.IsInsertState(cmState);
				get_hmmStates.push_back(add);
			}
		}

		InfernalHmm::State hmmRightState=cm2HmmStateVector[cmState].hmmRightState;
		if (hmmRightState!=GetInvalidState()) {
			if (!hmmStateUsed[hmmRightState]) {
				hmmStateUsed[hmmRightState]=true;
				StateFromNode add;
				add.state=hmmRightState;
				add.isRight=true;
				add.isNormal=!cm.IsInsertState(cmState);
				get_hmmStates.push_back(add);
			}
		}
	}
}
void InfernalHmm::GetNormalHmmStatesOfLeftOrRightByCmNode (StateList& get_stateList,const CovarianceModel& cm,CovarianceModel::Node cmNode,bool getIsRight) const
{
	get_stateList.clear();

	StateFromNodeList stateFromNodeList;
	GetHmmStatesByCmNode (stateFromNodeList,cm,cmNode);
	for (StateFromNodeList::iterator i=stateFromNodeList.begin(); i!=stateFromNodeList.end(); i++) {

		if (i->isNormal && i->isRight==getIsRight) {
			get_stateList.push_back(i->state);
		}
	}
}
void InfernalHmm::GetHmmEdgesByCmNode (std::list<EdgeInfo>& get_edgeInfoList,const StateFromNodeList& hmmStateList,const CovarianceModel& cm,CovarianceModel::Node cmNode)
{
	BuildNonSavedInfoIfNecessary();

	get_edgeInfoList.clear();

	for (StateFromNodeList::const_iterator i=hmmStateList.begin(); i!=hmmStateList.end(); i++) {

		InfernalHmm::State fromState=i->state;

		bool isRight=i->isRight;

		if (!isRight) {
			assert(IsLeftState(fromState)); // could also be right, but I don't care
			for (int childNum=0; childNum<GetNumChildren(fromState); childNum++) {
				InfernalHmm::State toState=GetNthChildState(fromState,childNum);
				EdgeInfo edge;
				edge.fromState=fromState;
				edge.toState=toState;
				edge.childNum=childNum;
				edge.isRightSide=false;
				get_edgeInfoList.push_back(edge);
			}
		}
		else {
			assert(IsRightState(fromState));

			for (ParentAndMyChildNumVector::const_iterator parentIter=nonSavedInfoVector[fromState].parentAndMyChildNumVector.begin(); parentIter!=nonSavedInfoVector[fromState].parentAndMyChildNumVector.end(); parentIter++) {
				EdgeInfo edge;
				edge.fromState=parentIter->parentState;
				edge.toState=fromState;
				edge.childNum=parentIter->myChildNum;
				edge.isRightSide=true;
				get_edgeInfoList.push_back(edge);
			}
		}
	}
}
void InfernalHmm::CopyProbabilitiesForNode(const InfernalHmm& sourceHmm,const CovarianceModel& cm,const CovarianceModel::Node cmNode)
{
	throw SimpleStringException("Called deprecated function.  I believe this isn't used & it doesn't support local stuff.  %s:%d.",__FILE__,__LINE__);
	StateFromNodeList hmmStateList;
	GetHmmStatesByCmNode(hmmStateList,cm,cmNode);
	EdgeInfoList edgeInfoList;
	GetHmmEdgesByCmNode (edgeInfoList,hmmStateList,cm,cmNode);

	// emissions
	for (StateFromNodeList::iterator i=hmmStateList.begin(); i!=hmmStateList.end(); i++) {

		InfernalHmm::State fromState=i->state;

		if (IsEmitting(fromState)) {
			for (int nuc=0; nuc<Alphabet_size; nuc++) {
				SetSingletEmissionLogScore(fromState,nuc,sourceHmm.GetSingletEmissionScore(fromState,nuc));
			}
		}
	}

	// transitions (i.e. edges)
	for (EdgeInfoList::iterator i=edgeInfoList.begin(); i!=edgeInfoList.end(); i++) {
		State fromState=i->fromState;
		//State toState=i->toState;
		int childNum=i->childNum;

		SetTransitionLogScore(fromState,childNum,sourceHmm.GetNthChildTsc(fromState,childNum));
	}
}
void InfernalHmm::DumpProbabilitiesCsvByNode(FILE *nodeDump,const CovarianceModel& cm,bool doValues)
{
	BuildNonSavedInfoIfNecessary();

	CovarianceModel::Node cmNode;
	StateFromNodeList hmmStateList;
	EdgeInfoList edgeInfoList;
	for (cmNode=cm.GetFirstNode(); cmNode!=cm.GetLastNode(); cmNode++) {

		if (doValues) {
			fprintf(nodeDump,",,");
		}
		else {
			fprintf(nodeDump,",node #%d,",CovarianceModel::NodeToInt(cmNode));
		}

		GetHmmStatesByCmNode(hmmStateList,cm,cmNode);
		GetHmmEdgesByCmNode (edgeInfoList,hmmStateList,cm,cmNode);

		// emissions
		for (StateFromNodeList::iterator i=hmmStateList.begin(); i!=hmmStateList.end(); i++) {

			InfernalHmm::State fromState=i->state;
			if (IsEmitting(fromState)) {
				for (int nuc=0; nuc<Alphabet_size; nuc++) {
					if (doValues) {
						fprintf(nodeDump,"%f,",GetSingletEmissionScore(fromState,nuc));
					}
					else {
						fprintf(nodeDump,"%d @ %c,",StateToInt(fromState),nucs[nuc]);
					}
				}
			}
		}

		// transitions (i.e. edges)
		for (EdgeInfoList::iterator i=edgeInfoList.begin(); i!=edgeInfoList.end(); i++) {
			State fromState=i->fromState;
			State toState=i->toState;
			int childNum=i->childNum;

			if (doValues) {
				fprintf(nodeDump,"%f,",GetNthChildTsc(fromState,childNum));
			}
			else {
				if (i->isRightSide) {
					std::swap(fromState,toState); // just for printing's sake
				}
				fprintf(nodeDump,"%d --> %d,",StateToInt(fromState),StateToInt(toState));
			}
		}
	}
	fprintf(nodeDump,"\n");
}
bool InfernalHmm::AreSingletEmissionScoresSameForAllNucs(const State& state) const
{
	assert(IsEmitting(state));
	for (int nuc=0; nuc<Alphabet_size; nuc++) {
		if (GetSingletEmissionScore(state,nuc)!=GetSingletEmissionScore(state,0)) {
			return false;
		}
	}
	return true;
}
void InfernalHmm::HmmNormalize (void)
{
	for (State state=GetFirstState(); state!=GetLastState(); state++) {

		if (IsEmitting(state)) {
			double sum=0;
			for (int nuc=0; nuc<Alphabet_size; nuc++) {
				double e=GetSingletEmissionScore(state,nuc);
				sum += pow2(e);
			}
			double logSum=sreLOG2(sum);
			for (int nuc=0; nuc<Alphabet_size; nuc++) {
				SetSingletEmissionLogScore(state,nuc,(float)(GetSingletEmissionScore(state,nuc) - logSum));
			}
		}

		double sum=0;
		for (int child=0; child<GetNumChildren(state); child++) {
			double e=GetNthChildTsc(state,child);
			sum += pow2(e);
		}
		double logSum=sreLOG2(sum);
		for (int child=0; child<GetNumChildren(state); child++) {
			SetTransitionLogScore(state,child,(float)(GetNthChildTsc(state,child)-logSum));
		}
	}
}
void InfernalHmm::ConvertToInfernalSavableFormat (const CovarianceModel& sourceCm)
{
	throw SimpleStringException("Sorry InfernalHmm::ConvertToInfernalSavableFormat doesn't work; Infernal is highly dependent on node type, and it'd take a fair amount of work to make my HMMs conform to the nodes defined by Infernal, & this info is used in Infernal's alignment code.");

	// alg: go thru HMM states in order.  Whenever the corresponding CM node or the left/right side of the HMM state changes, that's a new node

	// first alloc node stuff.  GetNumStates() is an upper bound on # of nodes we'll need
	cm->ndtype = (char *)MallocOrDie(GetNumStates()  * sizeof(char));
	cm->nodemap= (int *)MallocOrDie(GetNumStates()  * sizeof(int));

	int nextNode=-1; // first HMM state will cause an increment, so we'll start at 0
	CovarianceModel::Node prevCmNode=CovarianceModel::GetInvalidNode();
	bool prevIsLeft=false; // just init to avoid compiler warnings
	for (InfernalHmm::State hmmState=GetFirstState(); hmmState!=GetLastState(); hmmState++) {

		// PASSTHRU_st -> S_st
		if (GetStateType(hmmState)==PASSTHRU_st) {
			cm->sttype[InfernalHmm::StateToInt(hmmState)]=S_st;
		}

		// what node are we in
		CovarianceModel::Node cmNode=CovarianceModel::GetInvalidNode();
		for (std::list<CovarianceModel::State>::iterator cmStateIter=hmm2CmStateVector[hmmState].begin(); cmStateIter!=hmm2CmStateVector[hmmState].end(); cmStateIter++) {
			CovarianceModel::State cmState=*cmStateIter;
			CovarianceModel::Node thisNode=sourceCm.GetNode(cmState);
			if (cmNode!=CovarianceModel::GetInvalidNode()) {
				//assert(cmNode==thisNode); // they should all map to the same node, or something's weird
				// this happens with states that are from splicing after bifurcation.  Just take the lowest-# state; it odesn't really matter
				if (thisNode<cmNode) {
					cmNode=thisNode;
				}
			}
			cmNode=thisNode;
		}
		assert(cmNode!=CovarianceModel::GetInvalidNode());
		//printf("%d\n",CovarianceModel::NodeToInt(cmNode));

		// did it change?
		if (prevCmNode!=cmNode || prevIsLeft!=IsLeftState(hmmState)) {
			// a change
			nextNode++;
			cm->nodemap[nextNode]=InfernalHmm::StateToInt(hmmState);
			if (GetStateType(hmmState)==S_st) {
				cm->ndtype[nextNode]=MATL_nd;
			}
			else {
				if (IsEndState(hmmState)) {
					cm->ndtype[nextNode]=END_nd;
				}
				else {
					//crap: on the right side, the nodes are in the reverse order.  Maybe I can just call it MATL_nd... assert(GetStateType(hmmState)==ML_st); // what other type of HMM node is there?
					cm->ndtype[nextNode]=MATL_nd;
				}
			}
		}
		cm->ndidx[InfernalHmm::StateToInt(hmmState)]=nextNode;

		// remember previous
		prevCmNode=cmNode;
		prevIsLeft=IsLeftState(hmmState);
	}

	cm->nodes=nextNode;
	cm->null[0]=cm->null[1]=cm->null[2]=cm->null[3]=1;

	for (int s=0; s<GetNumStates(); s++) {
		if (IsEmitting(IntToState(s))) {
			for (int nuc=0; nuc<Alphabet_size; nuc++) {
				cm->e[s][nuc]=(float)(pow2(cm->esc[s][nuc]));
			}
		}
		for (int child=0; child<cm->cnum[s]; child++) {
			cm->t[s][child]=(float)(pow2(cm->tsc[s][child]));
		}
	}

	// set plast/pnum
	for (InfernalHmm::State hmmState=GetFirstState(); hmmState!=GetLastState(); hmmState++) {
		int i_hmmState=InfernalHmm::StateToInt(hmmState);
		cm->pnum[i_hmmState]=0;
		cm->plast[i_hmmState]=-1;
	}
	for (InfernalHmm::State hmmState=GetFirstState(); hmmState!=GetLastState(); hmmState++) {
		for (int child=0; child<GetNumChildren(hmmState); child++) {
			int i_hmmState=InfernalHmm::StateToInt(hmmState);
			int i_childState=InfernalHmm::StateToInt(GetNthChildState(hmmState,child));
			cm->plast[i_childState]=i_hmmState;
			cm->pnum[i_childState]++;
		}
	}
}
void InfernalHmm::DumpCmAndHmmCorrespondingScores (const char *outFileName,const CovarianceModel& cm)
{
	FILE *out=ThrowingFopen(outFileName,"wt");

	fprintf(out,"CM node,node type,CM state,state type,emit/transition?,emit left nuc OR transition to state#,emit right nuc OR transition to state type,CM score,HMM left score,HMM right score\n");


	for (CovarianceModel::State cmState=cm.GetFirstState(); cmState!=cm.GetLastState(); cmState++) {

		if (cm.IsBifurcation(cmState)) {
			continue;
		}

		fprintf(out,"\n");

		InfernalHmm::State leftHmmState=GetHmmLeftStateOfCmState(cmState);
		InfernalHmm::State rightHmmState=GetHmmRightStateOfCmState(cmState);

		if (cm.IsEmittingState(cmState)) {
			int numLeftNucs=cm.EmitsLeft(cmState) ? Alphabet_size : 1;
			int numRightNucs=cm.EmitsRight(cmState) ? Alphabet_size : 1;
			for (int leftNuc=0; leftNuc<numLeftNucs; leftNuc++) {
				for (int rightNuc=0; rightNuc<numRightNucs; rightNuc++) {
					fprintf(out,"%d,%s,%d,%s,emit,",CovarianceModel::NodeToInt(cm.GetNode(cmState)),Nodetype(cm.GetNodeType(cm.GetNode(cmState))),CovarianceModel::StateToInt(cmState),Statetype(cm.GetStateType(cmState)));
					if (cm.EmitsLeft(cmState)) {
						fprintf(out,"%c",nucs[leftNuc]);
					}
					fprintf(out,",");
					if (cm.EmitsRight(cmState)) {
						fprintf(out,"%c",nucs[rightNuc]);
					}
					fprintf(out,",");
					if (cm.EmitsLeftAndRight(cmState)) {
						fprintf(out,"%f,",cm.GetPairEmissionScore(cmState,leftNuc,rightNuc));
					}
					else {
						const int nuc=cm.EmitsLeft(cmState) ? leftNuc : rightNuc;
						fprintf(out,"%f,",cm.GetSingletEmissionScore(cmState,nuc));
					}
					if (cm.EmitsLeft(cmState)) {
						fprintf(out,"%f",GetSingletEmissionScore(leftHmmState,leftNuc));
					}
					fprintf(out,",");
					if (cm.EmitsRight(cmState)) {
						fprintf(out,"%f",GetSingletEmissionScore(rightHmmState,rightNuc));
					}
					fprintf(out,",");
					fprintf(out,"\n");
				}
			}
		}

		for (int child=0; child<cm.GetNumChildren(cmState); child++) {
			CovarianceModel::State cmToState=cm.GetNthChildState(cmState,child);
			fprintf(out,"%d,%s,%d,%s,transition,%d,%s,%f,",CovarianceModel::NodeToInt(cm.GetNode(cmState)),Nodetype(cm.GetNodeType(cm.GetNode(cmState))),CovarianceModel::StateToInt(cmState),Statetype(cm.GetStateType(cmState)),CovarianceModel::StateToInt(cmToState),Statetype(cm.GetStateType(cmToState)),cm.GetNthChildTsc(cmState,child));

			InfernalHmm::State leftToHmmState=GetHmmLeftStateOfCmState(cmToState);
			if (leftToHmmState!=InfernalHmm::GetInvalidState() && leftHmmState!=InfernalHmm::GetInvalidState()) {
				int hmmChild=GetChildNum_Slow(leftHmmState,leftToHmmState);
				float tsc=GetNthChildTsc(leftHmmState,hmmChild);
				fprintf(out,"%f",tsc);
			}
			fprintf(out,",");

			InfernalHmm::State rightToHmmState=GetHmmRightStateOfCmState(cmToState);
			if (rightToHmmState!=InfernalHmm::GetInvalidState() && rightHmmState!=InfernalHmm::GetInvalidState()) {
				int hmmChild=GetChildNum_Slow(rightToHmmState,rightHmmState);
				float tsc=GetNthChildTsc(rightToHmmState,hmmChild);
				fprintf(out,"%f",tsc);
			}
			fprintf(out,",");
			fprintf(out,"\n");
		}
	}

	fclose(out);
}
CovarianceModel::NodeList InfernalHmm::GetCmNode (const CovarianceModel& cm,State hmmState) const
{
	CovarianceModel::NodeList nodeList;
	const CovarianceModel::StateList& stateList=GetCmState(hmmState);
	for (CovarianceModel::StateList::const_iterator i=stateList.begin(); i!=stateList.end(); i++) {
		CovarianceModel::Node node=cm.GetNode(*i);
		if (std::find(nodeList.begin(),nodeList.end(),node)==nodeList.end()) {
			nodeList.push_back(node);
		}
	}

	return nodeList;
}
void InfernalHmm::ComputePairInflationMatrix (PairInflationMatrix& inflationPerPair,const CovarianceModel& cm,CovarianceModel::State cmConsensusState) const
{
	assert(cm.GetStateType(cmConsensusState)==MP_st); // otherwise it's weird to call this function

	const InfernalHmm::State leftHmmConsensusState=GetHmmLeftStateOfCmState(cmConsensusState);
	const InfernalHmm::State rightHmmConsensusState=GetHmmRightStateOfCmState(cmConsensusState);

	float minInflation=+FLT_MAX; // okay really, the min inflation should be 0, but I'm cheating by calculating this under the assumption that the minimum "inflation" (HMM score minus CM score) is like inflation of 0.  So, later I subtract this out
	for (int leftNuc=0; leftNuc<Alphabet_size; leftNuc++) {
		for (int rightNuc=0; rightNuc<Alphabet_size; rightNuc++) {
			float hmmScore=GetSingletEmissionScore(leftHmmConsensusState,leftNuc) + GetSingletEmissionScore(rightHmmConsensusState,rightNuc);
			inflationPerPair[leftNuc][rightNuc] = hmmScore - cm.GetPairEmissionScore(cmConsensusState,leftNuc,rightNuc);
			minInflation=std::min(minInflation,inflationPerPair[leftNuc][rightNuc]);
		}
	}
	for (int leftNuc=0; leftNuc<Alphabet_size; leftNuc++) {
		for (int rightNuc=0; rightNuc<Alphabet_size; rightNuc++) {
			inflationPerPair[leftNuc][rightNuc] -= minInflation;
		}
	}
}
Cm2Hmm_HmmBuildType InfernalHmm::GetHmmBuildType (void) const
{
	return hmmBuildType;
}
void InfernalHmm::Save (char *f)
{
	assertr(false); // Sorry, but this function is only for CMs.  Sorry about the whole public inheritance stuff.
}
void InfernalHmm::SaveInFormat (const char *hmmFileName,HmmFileFormat hmmFileFormat)
{
	switch (hmmFileFormat) {
		case HmmFileFormat_Binary:
			SaveInDeprecatedBinaryFormat(hmmFileName); // really save in binary
			break;
		case HmmFileFormat_Text:
			SaveInBinary(hmmFileName);
			break;
		case HmmFileFormat_HMMER:
			throw SimpleStringException("Save format HMMER not implemented yet.");
			if (hmmBuildType!=HmmBuildType_Original) {
				throw SimpleStringException("This HMM is not a compact-type profile HMM, so it cannot be saved in HMMER format.");
			}
			break;
		default:
			throw SimpleStringException("Invalid save file format set, or it wasn't implemented.  Sorry, it's undoubtably my fault.");
	}
}
HmmFileFormat InfernalHmm::GetDefaultHmmFileFormat (void)
{
	return HmmFileFormat_Text;
}
void InfernalHmm::SetAdjustedEndsc(State hmmLeftState,State hmmRightState,const CovarianceModel& cm,CovarianceModel::State cmState,int linkNum)
{
	// endsc happens after emits, and the HMM emit scores may have a constant added to them relative to the CM (which is taken care of in transition scores)
	float cmEndsc=cm.GetEndsc(cmState);
	float hmmEndsc;

	if (cm.IsEmittingState(cmState)) {
		if (cm.EmitsLeft(cmState)) {
			if (cm.EmitsRight(cmState)) {
				// emits pair
				float minInflation=+FLT_MAX;
				for (int leftNuc=0; leftNuc<Alphabet_size; leftNuc++) {
					for (int rightNuc=0; rightNuc<Alphabet_size; rightNuc++) {
						float cmScore=cm.GetPairEmissionScore(cmState,leftNuc,rightNuc);
						float hmmScore=GetSingletEmissionScore(hmmLeftState,leftNuc) + GetSingletEmissionScore(hmmRightState,rightNuc);
						float inflation=hmmScore-cmScore;
						minInflation=std::min(minInflation,inflation);
					}
				}
				hmmEndsc=cmEndsc-minInflation;
			}
			else {
				// only left
				int arbitraryNuc=0; // they should all be the same (approximately)
				float adjustment=cm.GetSingletEmissionScore(cmState,arbitraryNuc) - GetSingletEmissionScore(hmmLeftState,arbitraryNuc);
				hmmEndsc=cmEndsc+adjustment;
			}
		}
		else {
			// only right
			int arbitraryNuc=0; // they should all be the same (approximately)
			float adjustment=cm.GetSingletEmissionScore(cmState,arbitraryNuc) - GetSingletEmissionScore(hmmRightState,arbitraryNuc);
			hmmEndsc=cmEndsc+adjustment;
		}
	}
	else {
		hmmEndsc=cmEndsc;
	}

	SetEndscLinkToLeft_Endsc(hmmRightState,linkNum,hmmEndsc);
}
void InfernalHmm::SetAdjustedEndscGlobally(const CovarianceModel& cm)
{
	CovarianceModel::State cmState;
	for (cmState=cm.GetFirstState(); cmState!=cm.GetLastState(); cmState++) {

		if (cm.GetNode(cmState)!=cm.GetFirstNode()) { // don't do anything for the whole first node, since that gets clobbered by the ConfigLocal function in modelconfig.c. 
			InfernalHmm::State hmmLeftState=GetHmmLeftStateOfCmState(cmState);
			InfernalHmm::State hmmRightState=GetHmmRightStateOfCmState(cmState);

			if (cm.GetEndsc(cmState)!=(float)IMPOSSIBLE) {

				// link right state to this left state, for simulating endsc
				if (hmmLeftState==InfernalHmm::GetInvalidState() || hmmRightState==InfernalHmm::GetInvalidState()) {
					assert(false);
					throw SimpleStringException("For making an HMM local, it's pretty tricky if the insert states (IL,IR) in the CM have non-IMPOSSIBLE endsc scores, because in my mapping, IL,IR states don't have both a left&right hmmState, and it's not easy to see what it should be.  So, I just assumed it wouldn't happen, making an ass out of me and you, the anonymous user of this program.  Sorry.");
				}

				// set the link related via the current cmState
				bool foundState=false;
				for (int link=0; link<GetNumEndscLinksToLeft(hmmRightState); link++) {
					if (GetEndscLinkToLeft_State(hmmRightState,link)==hmmLeftState) {
						assert(!foundState);
						SetAdjustedEndsc(hmmLeftState,hmmRightState,cm,cmState,link);
						foundState=true;
					}
				}
				assert(foundState);
			}
		}
	}
}

char *GetHmmBuildTypeDescription (Cm2Hmm_HmmBuildType hmmType)
{
	switch (hmmType) {
	case HmmBuildType_Original:
		return "compact";
	case HmmBuildType_separateMPandMLMR:
		return "expanded";
	case HmmBuildType_separateMPMLMRD:
		return "overexpanded"; // deprecated too
	default:
		assertr(false);
		return "unknown type";
	}
}
Cm2Hmm_HmmBuildType GetHmmBuildTypeByText (const char *text)
{
	if (strcmp(text,"compact")==0) {
		return HmmBuildType_Original;
	}
	if (strcmp(text,"expanded")==0) {
		return HmmBuildType_separateMPandMLMR;
	}
	throw SimpleStringException("HMM type text string (e.g. compact or expanded) was not recognized -- is some input file in the wrong format?");
}
void InfernalHmm::AddStates(int numNewStates)
{
	int newNumStates=numNewStates+GetNumStates();

	cm->sttype=(char *)ReallocOrDie(cm->sttype,sizeof(cm->sttype[0])*(newNumStates+1));
	cm->stid=(char *)ReallocOrDie(cm->stid,sizeof(cm->stid[0])*(newNumStates+1));
	cm->ndidx=(int *)ReallocOrDie(cm->ndidx,sizeof(cm->ndidx[0])*newNumStates);
	cm->cfirst=(int *)ReallocOrDie(cm->cfirst,sizeof(cm->cfirst[0])*newNumStates);
	cm->cnum=(int *)ReallocOrDie(cm->cnum,sizeof(cm->cnum[0])*newNumStates);
	cm->plast=(int *)ReallocOrDie(cm->plast,sizeof(cm->plast[0])*newNumStates);
	cm->pnum=(int *)ReallocOrDie(cm->pnum,sizeof(cm->pnum[0])*newNumStates);

	cm->begin=(float *)ReallocOrDie(cm->begin,sizeof(cm->begin[0])*newNumStates);
	cm->end=(float *)ReallocOrDie(cm->end,sizeof(cm->end[0])*newNumStates);
	cm->beginsc=(float *)ReallocOrDie(cm->beginsc,sizeof(cm->beginsc[0])*newNumStates);
	cm->endsc=(float *)ReallocOrDie(cm->endsc,sizeof(cm->endsc[0])*newNumStates);

	Realloc2d(cm->tsc,GetNumStates(),newNumStates,MAXCONNECT);
	Realloc2d(cm->esc,GetNumStates(),newNumStates,Alphabet_size*Alphabet_size);

	Realloc2d(cm->t,GetNumStates(),newNumStates,MAXCONNECT);
	Realloc2d(cm->e,GetNumStates(),newNumStates,MAXCONNECT);

	cm->M=newNumStates;
}
void InfernalHmm::MoveStatesHigher (State st_firstState,State st_destOfFirstState)
{
	int firstState=StateToInt(st_firstState);
	int destOfFirstState=StateToInt(st_destOfFirstState);

	assert(destOfFirstState>=firstState);
	int increase=destOfFirstState-firstState;
	int newNumStates=GetNumStates() + increase;

	AddStates(increase);

	// do the move, careful about the overlap
	for (int state=newNumStates-1; state>=destOfFirstState; state--) {
		cm->sttype[state]=cm->sttype[state-increase];
		cm->cfirst[state]=cm->cfirst[state-increase];
		cm->cnum[state]=cm->cnum[state-increase];
		for (int i=0; i<MAXCONNECT; i++) {
			cm->tsc[state][i]=cm->tsc[state-increase][i];
		}
		for (int i=0; i<Alphabet_size*Alphabet_size; i++) {
			cm->esc[state][i]=cm->esc[state-increase][i];
		}
	}

	// now adjust 'cfirst' members of the moved thingies
	for (int state=destOfFirstState; state<newNumStates; state++) {
		if (cm->cfirst[state]!=-1) {
			assert(cm->cfirst[state]>=firstState);
			cm->cfirst[state] += increase;
		}
	}
}
void InfernalHmm::CopyStatesVerbatimFrom(const InfernalHmm& t,State st_firstState,State st_lastState)
{
	int firstState=StateToInt(st_firstState);
	int lastState=StateToInt(st_lastState);
	int state;
	for (state=firstState; state<lastState; state++) {

		cm->sttype[state]=t.cm->sttype[state];
		cm->cfirst[state]=t.cm->cfirst[state];
		cm->cnum[state]=t.cm->cnum[state];
		for (int i=0; i<MAXCONNECT; i++) {
			cm->tsc[state][i]=t.cm->tsc[state][i];
		}
		for (int i=0; i<Alphabet_size*Alphabet_size; i++) {
			cm->esc[state][i]=t.cm->esc[state][i];
		}
	}
}
void InfernalHmm::SetLeftwardBeginProbDirectly(State state,float p)
{
	otherStateInfoVector[state].leftBeginProb=p;
}
void InfernalHmm::SetRightwardBeginProbDirectly(State state,float p)
{
	otherStateInfoVector[state].rightBeginProb=p;
}
void InfernalHmm::SetEndProbDirectly(State state,int link,float p)
{
	endscLinksToLeftVector[state][link].endProb=p;
}
float InfernalHmm::GetLeftwardBeginProbDirectly(State state)
{
	return otherStateInfoVector[state].leftBeginProb;
}
float InfernalHmm::GetRightwardBeginProbDirectly(State state)
{
	return otherStateInfoVector[state].rightBeginProb;
}
float InfernalHmm::GetEndProbDirectly(State state,int link)
{
	return endscLinksToLeftVector[state][link].endProb;
}
void InfernalHmm::BuildReverseLocalEndByStateVector (ReverseLocalEndByStateVector& v)
{
	v.resize(GetNumStates());
	for (State state=GetFirstState(); state!=GetLastState(); state++) {
		for (int link=0; link<GetNumEndscLinksToLeft(state); link++) {
			ReverseLocalEnd t;
			t.rightState=state;
			t.linkNum=link;
			v[GetEndscLinkToLeft_State(state,link)].push_back(t);
		}
	}
}
void InfernalHmm::NormalizeTransitionsToStrictProbabilitiesViaProbabilities (void)
{
	// linking local ends from Right to Left makes sense with the reversed HmmType1, but in this orientation, we should do the link-to-left with its left state.  The transitions on the left state are mutually exclusive with the local end.  The transitions on the right state go out of the states that hte local end skips, so are not mutually exclusive
	// similarly, in this case, right local begins compete with transitions of the same state, while left local begins compete only with each other.  Sum over all left local begins should be 1.  Each state's right local begin, plus transitions (plus relevant local ends) should be 1. 
	if (DoLocal()) {
		// left local begins
		int numEvents=0;
		double totalProb=0;
		for (State state=GetFirstState(); state!=GetLastState(); state++) {
			if (GetLeftwardBeginProbDirectly(state)!=0) {
				numEvents++;
				totalProb += GetLeftwardBeginProbDirectly(state);
			}
		}
		for (State state=GetFirstState(); state!=GetLastState(); state++) {
			float prob=ToProb(GetLeftwardBeginProbDirectly(state),totalProb,numEvents);
			SetLeftwardBeginProbDirectly(state,prob);
			SetLeftwardBeginsc(state,(float)(sreLOG2(prob)));
		}
	}

	ReverseLocalEndByStateVector reverseLocalEndByStateVector;
	if (DoLocal()) {
		BuildReverseLocalEndByStateVector(reverseLocalEndByStateVector);
	}

	for (State state=GetFirstState(); state!=GetLastState(); state++) {

		if (state==IntToState(56)) {
		  //int q=9;
		}

		int numEvents=0;
		double totalProb=0;
		if (DoLocal()) {
			if (GetRightwardBeginProbDirectly(state)!=0) {
				totalProb += GetRightwardBeginProbDirectly(state);
				numEvents++;
			}
			for (ReverseLocalEndList::iterator i=reverseLocalEndByStateVector[state].begin(); i!=reverseLocalEndByStateVector[state].end(); i++) {
				totalProb += GetEndProbDirectly(i->rightState,i->linkNum);
			}
			numEvents += (int)(reverseLocalEndByStateVector[state].size());
		}
		for (int child=0; child<GetNumChildren(state); child++) {
			totalProb += GetTransitionProbDirectly(state,child);
		}
		numEvents += GetNumChildren(state);

		if (DoLocal()) {

			float prob;
			if (GetRightwardBeginProbDirectly(state)!=0) {
				prob=ToProb(GetRightwardBeginProbDirectly(state),totalProb,numEvents);
				SetRightwardBeginProbDirectly(state,prob);
				SetRightwardBeginsc(state,(float)(sreLOG2(prob)));
			}
			else {
				SetRightwardBeginsc(state,(float)IMPOSSIBLE);
			}

			for (ReverseLocalEndList::iterator i=reverseLocalEndByStateVector[state].begin(); i!=reverseLocalEndByStateVector[state].end(); i++) {
				float prob=ToProb(GetEndProbDirectly(i->rightState,i->linkNum),totalProb,numEvents);
				SetEndProbDirectly(i->rightState,i->linkNum,prob);
				SetEndscLinkToLeft_Endsc(i->rightState,i->linkNum,(float)(sreLOG2(prob)));
			}
		}
		for (int child=0; child<GetNumChildren(state); child++) {
			float prob=ToProb(GetTransitionProbDirectly(state,child),totalProb,numEvents);
			SetTransitionProbDirectly(state,child,prob);
			SetTransitionLogScore(state,child,(float)(sreLOG2(prob)));
		}
	}
}
void InfernalHmm::SetLocalBeginsLikeGlobal(void)
{
	for (State state=GetFirstState(); state!=GetLastState(); state++) {
		SetLeftwardBeginsc(state,(float)IMPOSSIBLE);
		SetRightwardBeginsc(state,(float)IMPOSSIBLE);
	}
	SetLeftwardBeginsc(GetFirstState(),0);
	SetRightwardBeginsc(GetActualLastState(),0);
}
void InfernalHmm::SetLocalEndsToImpossible(void)
{
	for (State state=GetFirstState(); state!=GetLastState(); state++) {
		SetNumEndscLinksToLeft(state,0);
	}
}
void InfernalHmm::SetLocalEndsToImpossibleButPreserveLinks(void)
{
	for (State state=GetFirstState(); state!=GetLastState(); state++) {
		for (int link=0; link<GetNumEndscLinksToLeft(state); link++) {
			SetEndscLinkToLeft_Endsc(state,link,(float)IMPOSSIBLE);
			SetEndProbDirectly(state,link,0);
		}
	}
}
void InfernalHmm::AddToLocalBeginsExceptRoot (float add)
{
	for (State state=GetFirstState(); state!=GetLastState(); state++) {
		if (state!=GetFirstState() && state!=GetActualLastState()) { // avoid root state
			if (GetLeftwardBeginsc(state)!=(float)IMPOSSIBLE) {
				SetLeftwardBeginsc(state,add+GetLeftwardBeginsc(state));
			}
			if (GetRightwardBeginsc(state)!=(float)IMPOSSIBLE) {
				SetRightwardBeginsc(state,add+GetRightwardBeginsc(state));
			}
		}
	}
}
void InfernalHmm::AddToLocalEnds (float add)
{
	for (State state=GetFirstState(); state!=GetLastState(); state++) {
		for (int link=0; link<GetNumEndscLinksToLeft(state); link++) {
			if (GetEndscLinkToLeft_Endsc(state,link)!=(float)IMPOSSIBLE) {
				SetEndscLinkToLeft_Endsc(state,link,add+GetEndscLinkToLeft_Endsc(state,link));
			}
		}
	}
}

void ConvertLinearCMToHmm (InfernalHmm& hmm,const CovarianceModel& cm,bool doLocalAlignment,char *cmFileNameJustForRecordKeeping,std::string programParams)
{
	Cm2Hmm_HmmBuildType hmmType=HmmBuildType_Original;
	InfernalHmm infernalHmm;
	ExtraCm2HmmInfo extraCm2HmmInfo;
	extraCm2HmmInfo.actuallySolveScores=false;

	Cm2Hmm_WithWeighting_NoCaching (infernalHmm,hmmType,cm,cmFileNameJustForRecordKeeping,NULL,&extraCm2HmmInfo);
	infernalHmm.BuildNonSavedInfoIfNecessary();
	infernalHmm.AddBuildDescription(programParams);
	infernalHmm.SetHmmBuildType(hmmType);

	for (InfernalHmm::State v=infernalHmm.GetFirstState(); v!=infernalHmm.GetLastState(); v++) {
		infernalHmm.SetRightwardBeginsc(v,(float)IMPOSSIBLE); // start off here, since later we'll want to set it to the max
	}

	for (CovarianceModel::State cmState=cm.GetFirstState(); cmState!=cm.GetLastState(); cmState++) {
		if (cm.GetNode(cmState)==cm.GetFirstNode() && doLocalAlignment) {
			// skip this first node for local
			/*
			if (cm.GetStateType(cmState)==IR_st) {
				// it's the cursed IR state
				// murder it's unholy offspring in the HMM
				InfernalHmm::State hellspawnState=infernalHmm.GetHmmRightStateOfCmState(cmState);
				for (int nuc=0; nuc<MAXABET; nuc++) {
					infernalHmm.SetSingletEmissionLogScore(hellspawnState,nuc,(float)IMPOSSIBLE);
				}
				*/
		}
		else {
			InfernalHmm::State hmmState=infernalHmm.GetHmmLeftStateOfCmState(cmState);
			InfernalHmm::State rightHmmState=infernalHmm.GetHmmRightStateOfCmState(cmState);

			// copy children transitions
			for (int child=0; child<cm.GetNumChildren(cmState); child++) {
				float sc=cm.GetNthChildTsc(cmState,child);
				infernalHmm.SetTransitionLogScore(hmmState,child,sc);
			}

			// copy beginsc
			float beginsc=cm.GetBeginsc(cmState);
			infernalHmm.SetLeftwardBeginsc(hmmState,beginsc);
			infernalHmm.SetRightwardBeginsc(hmmState,(float)IMPOSSIBLE); // left states never need it

			// rightward
			if (rightHmmState!=InfernalHmm::GetInvalidState()) {
				float oldBeginsc=infernalHmm.GetRightwardBeginsc(rightHmmState);
				if (beginsc>oldBeginsc) {
					infernalHmm.SetRightwardBeginsc(rightHmmState,beginsc);
				}

				// copy endsc
				float endsc=cm.GetEndsc(cmState);
				for (int link=0; link<infernalHmm.GetNumEndscLinksToLeft(rightHmmState); link++) {
					InfernalHmm::State linkState=infernalHmm.GetEndscLinkToLeft_State(rightHmmState,link);
					if (linkState==hmmState) {
						infernalHmm.SetEndscLinkToLeft_Endsc(rightHmmState,link,endsc);
					}
				}
			}

			// copy emissions, if any
			if (cm.IsEmitting(cmState)) {
				assertr(infernalHmm.IsEmittingState(hmmState));
				for (int nuc=0; nuc<MAXABET; nuc++) {
					float sc=cm.GetSingletEmissionScore(cmState,nuc);
					infernalHmm.SetSingletEmissionLogScore(hmmState,nuc,sc);
				}
			}
		}
	}

	hmm.CopyFrom(infernalHmm);
}

void ConvertLinearCMToHmm (InfernalHmm& hmm,char *cmFileName,bool doLocalAlignment,std::string programParams)
{
	CovarianceModel cm;
	cm.Load(cmFileName,doLocalAlignment);

	ConvertLinearCMToHmm (hmm,cm,doLocalAlignment,cmFileName,programParams);
}
