#include "Primitive_Observable.H"
#include "MyStrStream.H"
#include "MathTools.H"

using namespace ANALYSIS;
using namespace ATOOLS;



Multiplicity::Multiplicity(int type,double xmin,double xmax,int nbins,
			   const std::string & listname) :
  Primitive_Observable_Base(type,xmin,xmax,nbins,NULL)
{
  if (listname!="") {
    m_listname = listname;
    m_name = listname+"_multi.dat";
  }
  else
    m_name = "multi.dat";
}


void Multiplicity::Evaluate(const ATOOLS::Particle_List & pl,
			    double weight, int ncount)
{
  p_histo->Insert(pl.size(),weight,ncount); 
}


Primitive_Observable_Base * Multiplicity::Copy() const {
  return new Multiplicity(m_type,m_xmin,m_xmax,m_nbins,m_listname);
}





void Differential_Jetrate::Evaluate(int,const Vec4D *,
				    const Flavour *,double value, int ncount) 
{
  p_histo->Insert(p_sel->ActualValue()[0],value, ncount);
}


void Differential_Jetrate::Evaluate(const Particle_List &,
				    double value, int ncount) 
{
  p_histo->Insert(p_sel->ActualValue()[0],value, ncount);
}
