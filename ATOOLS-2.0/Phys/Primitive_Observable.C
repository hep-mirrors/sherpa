#include "Primitive_Observable.H"
#include "MyStrStream.H"
#include "MathTools.H"

using namespace ATOOLS;

void Differential_Jetrate::Evaluate(double value) {
  Evaluate(nout,moms,flavs,value);
}


void Differential_Jetrate::Evaluate(int,Vec4D *,
				    Flavour *,double value) 
{
  histo->Insert(sel->ActualValue()[0],value);
}

void Differential_Jetrate::Evaluate(const Parton_List &,double value) 
{
  histo->Insert(sel->ActualValue()[0],value);
}

