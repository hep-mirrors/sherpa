#include "Primitive_Observable.H"
#include "MyStrStream.H"
#include "MathTools.H"

using namespace APHYTOOLS;

void Differential_Jetrate::Evaluate(double value) {
  Evaluate(nout,moms,flavs,value);
}


void Differential_Jetrate::Evaluate(int,AMATOOLS::Vec4D *,
				    APHYTOOLS::Flavour *,double value) 
{
  histo->Insert(sel->ActualValue()[0],value);
}

void Differential_Jetrate::Evaluate(const APHYTOOLS::Parton_List &,double value) 
{
  histo->Insert(sel->ActualValue()[0],value);
}

