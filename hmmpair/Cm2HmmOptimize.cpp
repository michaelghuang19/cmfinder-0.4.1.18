#include "hmmpair.h"

void MakeCmHmmForOptimization (CovarianceModel& cm,InfernalHmm& infernalHmm,ExtraCm2HmmInfo& extraCm2HmmInfo,char *cmFileName,bool doLocalAlignment,const std::string& programParams,Cm2Hmm_HmmBuildType hmmType,bool isRigorous)
{
	cm.Load(cmFileName,doLocalAlignment);

#ifdef DISABLE_LPSOLVE
	if (isRigorous) {
		throw SimpleStringException("LPSOLVE is required for rigorous filters.  Please enable it and re-compile cmzasha.  (This process is described in the reference manual.)");
	}
#endif
	ExtraCm2HmmInfo extraInfoTemp;
	extraInfoTemp.actuallySolveScores=isRigorous;
	Cm2Hmm_WithWeighting_NoCaching(infernalHmm,hmmType,cm,cmFileName,NULL,&extraInfoTemp);
	infernalHmm.BuildNonSavedInfoIfNecessary();
	infernalHmm.AddBuildDescription(programParams);
	infernalHmm.SetHmmBuildType(hmmType);

	// get inequalities
	extraCm2HmmInfo.actuallySolveScores=false;
	InfernalHmm dummyInfernalHmm;
	Cm2Hmm_WithWeighting_NoCaching (dummyInfernalHmm,hmmType,cm,cmFileName,NULL,&extraCm2HmmInfo);
}
