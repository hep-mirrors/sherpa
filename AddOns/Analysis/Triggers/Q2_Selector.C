#include "AddOns/Analysis/Triggers/Trigger_Base.H"
#include "AddOns/Analysis/Main/Primitive_Analysis.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"
namespace ANALYSIS {

  class Q2_Selector: public Trigger_Base {  
  protected:

    double m_xmin, m_xmax;
    int    m_mode;

  public:

    inline Q2_Selector(const double min,const double max,const int mode,
		       const std::string &inlist,const std::string &outlist):
      Trigger_Base(inlist,outlist), m_xmin(min), m_xmax(max), m_mode(mode) {}
    
    void Evaluate(const ATOOLS::Particle_List &inlist,
		  ATOOLS::Particle_List &outlist,double value,double ncount)
    {
      const ATOOLS::Blob_List *p_bl(p_ana->GetBlobList());
      ATOOLS::Vec4D l(ATOOLS::rpa->gen.PBeam(0)), lp;
      ATOOLS::Blob *me(p_bl->FindFirst(ATOOLS::btp::Signal_Process));
      for (int i(0);i<me->NOutP();++i)
	if (me->OutParticle(i)->Flav().IsLepton()) {
	  lp=me->OutParticle(i)->Momentum();
	  break;
	}
      double Q2(-(l-lp).Abs2());
      if (Q2>=m_xmin && Q2<=m_xmax) {
	outlist.resize(inlist.size());
	for (size_t i=0;i<inlist.size();++i) 
	  outlist[i] = new ATOOLS::Particle(*inlist[i]);
      }
    }

    Analysis_Object *GetCopy() const 
    {
      return new Q2_Selector(m_xmin,m_xmax,m_mode,m_inlist,m_outlist);
    }

  };// end of class Q2_Selector

}// end of namespace ANALYSIS

#include "AddOns/Analysis/Main/Primitive_Analysis.H"

using namespace ANALYSIS;

#include "ATOOLS/Org/MyStrStream.H"

DECLARE_GETTER(Q2_Selector,"Q2Sel",Analysis_Object,Argument_Matrix);

Analysis_Object *ATOOLS::Getter
<Analysis_Object,Argument_Matrix,Q2_Selector>::operator()
(const Argument_Matrix &parameters) const
{
  if (parameters.size()<1 || parameters[0].size()<4) return NULL;
  int mode(parameters[0].size()>4?ATOOLS::ToType<int>(parameters[0][4]):0);
  return new Q2_Selector(ATOOLS::ToType<double>(parameters[0][0]),
			 ATOOLS::ToType<double>(parameters[0][1]),mode,
			 parameters[0][2],parameters[0][3]);
}									

void ATOOLS::Getter
<Analysis_Object,Argument_Matrix,Q2_Selector>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"min max inlist outlist [mode]";
}
