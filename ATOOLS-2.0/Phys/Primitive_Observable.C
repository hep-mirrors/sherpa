#include "Primitive_Observable.H"
#include "MyStrStream.H"
#include "MathTools.H"

using namespace ATOOLS;


void Differential_Jetrate::Evaluate(int,const Vec4D *,
				    const Flavour *,double value, int ncount) 
{
  p_histo->Insert(p_sel->ActualValue()[0],value, ncount);
}

void Differential_Jetrate::Evaluate(const Particle_List &,double value, int ncount) 
{
  p_histo->Insert(p_sel->ActualValue()[0],value, ncount);
}

