#include "hmmpair.h"

////////////////////////
// ObjectiveFunc

ObjectiveFunc::~ObjectiveFunc ()
{
}
double ObjectiveFunc::EvalActualObjectiveFuncLog2 (const vector<double>& problemVars)
{
	double fx;
	vector<double> gradient;
	vector2d<double> hessian;
	Eval(fx,gradient,hessian,problemVars,false);
	return fx/log(2.0);
}
bool ObjectiveFunc::SorryNoGradients (double& get_suggestedFiniteDifferenceAmount)
{
	return false;
}
SymbolicMath::Expression ObjectiveFunc::GetObjFuncExpression ()
{
	assertr(false); // ObjectiveFunc objects might not use SymbolicMath
}
double ObjectiveFunc::EvalValueOnly (const vector<double>& problemVars)
{
	double fx;
	vector<double> gradient;
	vector2d<double> hessian;
	Eval(fx,gradient,hessian,problemVars,false,false);
	return fx;
}
void ObjectiveFunc::GlobalToProblemVars (vector<double>& problemVars,const vector<double>& globalVars)
{
	assertr(false); // not implemented -- derived class needs to implement, or function shouldn't have been called
}
void ObjectiveFunc::UpdateGlobalVarsFromProblemVars (vector<double>& globalVars,const vector<double>& problemVars)
{
	assertr(false); // not implemented -- derived class needs to implement, or function shouldn't have been called
}
void ObjectiveFunc::LocalToProblemVars (vector<double>& problemVars,const vector<float>& localVars)
{
	assertr(false); // not implemented -- derived class needs to implement, or function shouldn't have been called
}
void ObjectiveFunc::ProblemToLocalVars (vector<float>& localVars,const vector<double>& problemVars)
{
	assertr(false); // not implemented -- derived class needs to implement, or function shouldn't have been called
}
const InequalityList& ObjectiveFunc::GetInequalityList (void)
{
	return emptyInequalityList;
}
const NonLinearConstraintList& ObjectiveFunc::GetNonLinearConstraintList (void)
{
	return emptyNonLinearConstraintList;
}
const InequalityList ObjectiveFunc::emptyInequalityList;
const NonLinearConstraintList ObjectiveFunc::emptyNonLinearConstraintList;


///////////////////////
// SolverWrapper

SolverWrapper::SolverWrapper (void)
{
}
SolverWrapper::~SolverWrapper ()
{
}
void SolverWrapper::SetMaxIters (int maxIters)
{
}

SolverWrapper::MessageReceiver::~MessageReceiver ()
{
}
bool SolverWrapper::MessageReceiver::EvaluatedObjectiveFunc (double functionValue,const vector<double>& problemVars)
{
	return false;
}
void SolverWrapper::MessageReceiver::PreEvaluateObjectiveFunc (const vector<double>& problemVars)
{
}
bool SolverWrapper::MessageReceiver::CarryOn (void)
{
	return true;
}



////////////////////
// GenericSymbolicObjectiveFunc

GenericSymbolicObjectiveFunc::GenericSymbolicObjectiveFunc (SymbolicMath& master_,const InequalityList& inequalityList_,const NonLinearConstraintList& nonLinearConstraintList_,int numProblemVars_)
: master(master_)
, inequalityList(inequalityList_)
, nonLinearConstraintList(nonLinearConstraintList_)
, numProblemVars(numProblemVars_)
{
}
GenericSymbolicObjectiveFunc::~GenericSymbolicObjectiveFunc ()
{
}
SymbolicMath::Expression GenericSymbolicObjectiveFunc::GetObjFuncExpression ()
{
	return master.GetExpression();
}
void GenericSymbolicObjectiveFunc::Eval (double& f,vector<double>& gradient,vector2d<double>& hessian,const vector<double>& problemVars,bool calculateHessian,bool calculateGradient)
{
	master.Eval(numProblemVars,f,gradient,hessian,problemVars,calculateHessian,calculateGradient);
}
int GenericSymbolicObjectiveFunc::GetNumProblemVars (void)
{
	return numProblemVars;
}
void GenericSymbolicObjectiveFunc::LocalToProblemVars (vector<double>& problemVars,const vector<float>& localVars)
{
	throw SimpleStringException("Not implemented %s:%d",__FILE__,__LINE__);
}
void GenericSymbolicObjectiveFunc::ProblemToLocalVars (vector<float>& localVars,const vector<double>& problemVars)
{
	throw SimpleStringException("Not implemented %s:%d",__FILE__,__LINE__);
}
const InequalityList& GenericSymbolicObjectiveFunc::GetInequalityList (void)
{
	return inequalityList;
}
const NonLinearConstraintList& GenericSymbolicObjectiveFunc::GetNonLinearConstraintList (void)
{
	return nonLinearConstraintList;
}


const NonLinearConstraintList& CachedObjectiveFunc::GetNonLinearConstraintList (void)
{
  assertr(false); // not implemented
}
