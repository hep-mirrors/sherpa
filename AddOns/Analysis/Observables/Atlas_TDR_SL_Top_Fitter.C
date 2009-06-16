#include "AddOns/Analysis/Observables/Primitive_Observable_Base.H"

#include "ATOOLS/Math/Variable.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Math/Histogram.H"
#include "ATOOLS/Org/Shell_Tools.H"
#include "ATOOLS/Org/Exception.H"
#include <iomanip>
#include <limits>

using namespace ATOOLS;

namespace ANALYSIS {

  class Atlas_TDR_SL_Top_Fitter: public Primitive_Observable_Base {
  private:
    Histogram *p_mjj, *p_mjjb;
  public:
    Atlas_TDR_SL_Top_Fitter(const std::string &inlist,
			    Primitive_Analysis *const ana);
    ~Atlas_TDR_SL_Top_Fitter();
    void Evaluate(const ATOOLS::Particle_List &list, 
		  double weight, double ncount);
    Analysis_Object &operator+=(const Analysis_Object &obj);
    void Restore(double scale=1.0);
    void EndEvaluation(double scale=1.0);
    void Output(const std::string & pname);
    Primitive_Observable_Base *Copy() const;    
  };// end of class Atlas_TDR_SL_Top_Fitter

} // namespace ANALYSIS

using namespace ANALYSIS;

DECLARE_GETTER(Atlas_TDR_SL_Top_Fitter_Getter,"AtlasTDRslTopFit",
 	       Primitive_Observable_Base,Argument_Matrix);

void Atlas_TDR_SL_Top_Fitter_Getter::PrintInfo
(std::ostream &str,const size_t width) const
{
  str<<"[list]";
}

Primitive_Observable_Base *
Atlas_TDR_SL_Top_Fitter_Getter::operator()
  (const Argument_Matrix &parameters) const
{
  std::string inlist("FinalState");
  if (parameters.size()<1) return NULL;
  if (parameters.size()==1) {
    if (parameters[0].size()<1) return NULL;
    return new Atlas_TDR_SL_Top_Fitter(parameters[0][0],parameters());
  }
  return NULL;
}

#include "AddOns/Analysis/Main/Primitive_Analysis.H"

using namespace ATOOLS;

Atlas_TDR_SL_Top_Fitter::Atlas_TDR_SL_Top_Fitter
(const std::string &inlist,Primitive_Analysis *const ana)
{
  p_ana=ana;
  m_listname=inlist;
  m_name="AtlasTDRslTopFit_"+m_listname;
  p_mjj = new Histogram(0,0.0,150.0,50);
  p_mjjb = new Histogram(0,0.0,400.0,100);
}

Analysis_Object &Atlas_TDR_SL_Top_Fitter::operator+=
(const Analysis_Object &obj)
{
  const Atlas_TDR_SL_Top_Fitter *vob((const Atlas_TDR_SL_Top_Fitter*)&obj);
  *p_mjj+=*vob->p_mjj;
  *p_mjjb+=*vob->p_mjjb;
  return *this;
}

void Atlas_TDR_SL_Top_Fitter::EndEvaluation(double scale) 
{
  p_mjj->Finalize();
  p_mjjb->Finalize();
  if (scale!=1.0) {
    p_mjj->Scale(scale);
    p_mjjb->Scale(scale);
  }
}

void Atlas_TDR_SL_Top_Fitter::Restore(double scale) 
{
  if (scale!=1.0) {
    p_mjj->Scale(scale);
    p_mjjb->Scale(scale);
  }
  p_mjj->Restore();
  p_mjjb->Restore();
}

void Atlas_TDR_SL_Top_Fitter::Output(const std::string & pname) 
{
  msg_Debugging()<<METHOD<<"(): {\n";
  std::string bname(pname+"/"+m_name);
  ATOOLS::MakeDir(pname); 
  p_mjj->Output((bname+"_mjj.dat").c_str());
  p_mjjb->Output((bname+"_mjjb.dat").c_str());
}

Atlas_TDR_SL_Top_Fitter::~Atlas_TDR_SL_Top_Fitter()
{
  delete p_mjj;
  delete p_mjjb;
}

void Atlas_TDR_SL_Top_Fitter::Evaluate
(const ATOOLS::Particle_List &list,double weight,double ncount)
{
  msg_Debugging()<<METHOD<<"(): {\n";
  size_t imin(0), jmin(0);
  double min(std::numeric_limits<double>::max());
  for (size_t i(0);i<list.size();++i) {
    if (list[i]->Flav().Kfcode()==kf_jet)
      for (size_t j(i+1);j<list.size();++j) 
	if (list[j]->Flav().Kfcode()==kf_jet) {
	  double cur((list[i]->Momentum()+list[j]->Momentum()).Mass());
	  msg_Debugging()<<"  i = "<<i<<", j = "<<j
			 <<" -> m_jj = "<<cur<<"\n";
	  if (dabs(cur-Flavour(kf_Wplus).Mass())<min) {
	    min=cur-Flavour(kf_Wplus).Mass();
	    imin=i;
	    jmin=j;
	  }
	}
  }
  if (imin==jmin) {
    p_mjj->Insert(0.0,0.0,ncount);
    p_mjjb->Insert(0.0,0.0,ncount);
    return;
  }
  msg_Debugging()<<"  min: i = "<<imin<<", j = "<<jmin<<" -> m_jj = "
		 <<(list[imin]->Momentum()+
		    list[jmin]->Momentum()).Mass()<<"\n";
  p_mjj->Insert((list[imin]->Momentum()+
		 list[jmin]->Momentum()).Mass(),weight,ncount);
  if (min>20.0) {
    p_mjjb->Insert(0.0,0.0,ncount);
    return;
  }
  size_t kmax(0);
  double max(0.0);
  for (size_t i(0);i<list.size();++i) {
    if (list[i]->Flav().Kfcode()==kf_bjet) {
      double cur((list[imin]->Momentum()+list[jmin]->Momentum()
		  +list[i]->Momentum()).PPerp2());
      msg_Debugging()<<"  k = "<<i<<" -> pt_jjb = "<<sqrt(cur)<<"\n";
      if (cur>max) {
	max=cur;
	kmax=i;
      }
    }
  }
  msg_Debugging()<<"  max: k = "<<kmax<<" -> m_jjb = "
		 <<(list[imin]->Momentum()+list[jmin]->Momentum()
		    +list[kmax]->Momentum()).Mass()<<"\n";
  p_mjjb->Insert((list[imin]->Momentum()+list[jmin]->Momentum()
		  +list[kmax]->Momentum()).Mass(),weight,ncount);
  msg_Debugging()<<"} done\n";
}

Primitive_Observable_Base *Atlas_TDR_SL_Top_Fitter::Copy() const
{
  return new Atlas_TDR_SL_Top_Fitter(m_listname,p_ana);
}

