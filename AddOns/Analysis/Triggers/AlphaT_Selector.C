#include "AddOns/Analysis/Triggers/Trigger_Base.H"
#include "ATOOLS/Org/Message.H"
namespace ANALYSIS {

  class AlphaT_Selector: public Trigger_Base {  
  protected:

    double m_xmin, m_xmax;
    int    m_mode;

    void MinDeltaHT(const std::vector<double> &et,
		    const double &ht,double &etj1,double &dmin,
		    const double &etj=0.0,const size_t &i=0)
    {
      if (i<et.size()) {
	MinDeltaHT(et,ht,etj1,dmin,etj+et[i],i+1);
	MinDeltaHT(et,ht,etj1,dmin,etj,i+1);
	return;
      }
      double etjt(2.0*etj>ht?etj:ht-etj);
      if (2.0*etjt-ht>=dmin) return;
      dmin=2.0*etjt-ht;
      etj1=etjt;
    }

  public:

    inline AlphaT_Selector(const double min,const double max,const int mode,
			   const std::string &inlist,const std::string &outlist):
      Trigger_Base(inlist,outlist), m_xmin(min), m_xmax(max), m_mode(mode) {}

    void Evaluate(const ATOOLS::Particle_List &inlist,
		  ATOOLS::Particle_List &outlist,double value,double ncount)
    {
      ATOOLS::Vec4D pmiss;
      double ht(0.0);
      std::vector<double> pt(inlist.size());
      for (size_t i(0);i<inlist.size();++i) {
	ATOOLS::Vec4D p=inlist[i]->Momentum();
	pmiss-=p;
	ht+=pt[i]=m_mode?p.EPerp():p.PPerp();
      }
      double etj1=0.0, dmin=ht;
      MinDeltaHT(pt,ht,etj1,dmin);
      double alphat=(1.0-etj1/ht)/sqrt(1.0-pmiss.PPerp2()/(ht*ht));
      msg_Debugging()<<"\\alpha_T = "<<alphat<<"\n";
      if (alphat>=m_xmin && alphat<=m_xmax) {
	outlist.resize(inlist.size());
	for (size_t i=0;i<inlist.size();++i) 
	  outlist[i] = new ATOOLS::Particle(*inlist[i]);
      }
    }

    Analysis_Object *GetCopy() const 
    {
      return new AlphaT_Selector(m_xmin,m_xmax,m_mode,m_inlist,m_outlist);
    }

  };// end of class AlphaT_Selector

}// end of namespace ANALYSIS

#include "AddOns/Analysis/Main/Primitive_Analysis.H"

using namespace ANALYSIS;

#include "ATOOLS/Org/MyStrStream.H"

DECLARE_GETTER(AlphaTSelector_Getter,"AlphaTSel",Analysis_Object,Argument_Matrix);

Analysis_Object *AlphaTSelector_Getter::operator()
(const Argument_Matrix &parameters) const
{
  if (parameters.size()<1 || parameters[0].size()<4) return NULL;
  int mode(parameters[0].size()>4?ATOOLS::ToType<int>(parameters[0][4]):1);
  return new AlphaT_Selector(ATOOLS::ToType<double>(parameters[0][0]),
			     ATOOLS::ToType<double>(parameters[0][1]),mode,
			     parameters[0][2],parameters[0][3]);
}									

void AlphaTSelector_Getter::PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"min max inlist outlist [mode]";
}
