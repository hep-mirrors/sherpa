#include "AddOns/Analysis/Observables/Primitive_Observable_Base.H"

#include "ATOOLS/Org/MyStrStream.H"
#include <iomanip>
#include "ATOOLS/Org/Message.H"

using namespace ATOOLS;

namespace ANALYSIS {


  class Conditional_One_Particle_Observable_Base :
    public Primitive_Observable_Base {
  protected:

    ATOOLS::Flavour m_flavour;
    double          m_cutoff;

    virtual void Evaluate(const Particle_Vector& pvec,
                          double weight=1., double ncount=1.) = 0;
  public:

    Conditional_One_Particle_Observable_Base
    (const ATOOLS::Flavour flav, const double cutoff,
     const int type, const double min, const double max, const int bins,
     const std::string &inlist, const std::string &name);

    void Evaluate(const ATOOLS::Particle_List &particlelist,
                  double weight=1.,double ncount=1);

  };// end of class Conditional_One_Particle_Observable_Base

  class Conditional_One_Particle_Multi_Emin :
    public Conditional_One_Particle_Observable_Base {
  private:
    bool m_on;
    void Evaluate(const Particle_Vector& pvec,
                  double weight=1., double ncount=1.);
  public:
    Conditional_One_Particle_Multi_Emin(const ATOOLS::Flavour flav,
                                        const double cutoff, const int type,
                                        const double min, const double max,
                                        const int bins,
                                        const std::string &inlist);

    Primitive_Observable_Base * Copy() const;

  };// end of class Multi_Emin_Distribution

}// end of namespace ANALYSIS

using namespace ANALYSIS;

template <class Class>
Primitive_Observable_Base *
GetObservable(const Analysis_Key& key)
{
  ATOOLS::Scoped_Settings s{ key.m_settings };
  const auto parameters = s.SetDefault<std::string>({}).GetVector<std::string>();
  if (parameters.size() < 6)
    THROW(missing_input, "Missing parameter values.");
  const auto list = parameters.size() > 6 ? parameters[6] : "FinalState";
  const auto kf = s.Interprete<int>(parameters[0]);
  ATOOLS::Flavour flav((kf_code)abs(kf));
  if (kf<0) flav=flav.Bar();
  return new Class(flav,s.Interprete<double>(parameters[1]),
                   HistogramType(parameters[5]),
                   s.Interprete<double>(parameters[2]),
                   s.Interprete<double>(parameters[3]),
                   s.Interprete<int>(parameters[4]),list);
}

#define DEFINE_GETTER_METHOD(CLASS,NAME)        \
  Primitive_Observable_Base *         \
  ATOOLS::Getter<Primitive_Observable_Base,Analysis_Key,CLASS>::operator()(const Analysis_Key& key) const \
  { return GetObservable<CLASS>(key); }

#define DEFINE_PRINT_METHOD(NAME)         \
  void ATOOLS::Getter<Primitive_Observable_Base,Analysis_Key,NAME>::PrintInfo(std::ostream &str,const size_t width) const \
  { str<<"[kf, Emin, min, max, bins, Lin|LinErr|Log|LogErr, list] ... list is optional"; }

#define DEFINE_OBSERVABLE_GETTER(CLASS,NAME,TAG)      \
  DECLARE_GETTER(CLASS,TAG,Primitive_Observable_Base,Analysis_Key); \
  DEFINE_GETTER_METHOD(CLASS,NAME)          \
  DEFINE_PRINT_METHOD(CLASS)

#include "AddOns/Analysis/Main/Primitive_Analysis.H"

Conditional_One_Particle_Observable_Base::
Conditional_One_Particle_Observable_Base
(const ATOOLS::Flavour flav, const double cutoff,
 const int type, const double min, const double max, const int bins,
 const std::string &inlist, const std::string &name):
  Primitive_Observable_Base(type,min,max,bins), 
  m_flavour(flav),
  m_cutoff(cutoff)
{
  m_listname=inlist;
  m_name = name + "_" + m_flavour.ShellName()
                + "_" + ToString(m_cutoff) + ".dat";
}

void Conditional_One_Particle_Observable_Base::Evaluate
(const ATOOLS::Particle_List& inlist, double weight, double ncount)
{
  Particle_Vector pvec;
  for (size_t i=0;i<inlist.size();++i) {
    if (inlist[i]->Flav()==m_flavour) {
      pvec.push_back(inlist[i]);
    }
  }
  Evaluate(pvec,weight,ncount);
}

//==============================================================================

DEFINE_OBSERVABLE_GETTER(Conditional_One_Particle_Multi_Emin,
                         Conditional_One_Particle_Multi_Emin_Getter,"MultiEmin")

Conditional_One_Particle_Multi_Emin::
Conditional_One_Particle_Multi_Emin
(const ATOOLS::Flavour flav, const double cutoff,
 const int type, const double min, const double max, const int bins,
 const std::string& inlist) :
  Conditional_One_Particle_Observable_Base(flav,cutoff,type,min,max,bins,
                                           inlist,"MultiEmin"), m_on(true)
{
  if (m_flavour.HadMass()==0. && m_cutoff==0.) {
    msg_Error()<<METHOD<<"() {\n"
               <<"  MultiEmin "<<flav.Kfcode()<<" "<<cutoff<<" "<<min<<" "<<max
               <<" "<<bins<<" Lin|LinErr|Log|LogErr "<<inlist<<"\n"
               <<"  Analysis: Multiplicity of "<<flav<<" with E > "<<cutoff
               <<".\n"
               <<"  The energy cut-off for massless particles must not be zero"
               <<" for this analysis to be sensible.\n"
               <<"  Switching it off ...\n"
               <<"}\n";
    m_name+=".meaningless_analysis";
    m_on=false;
  }
}

void Conditional_One_Particle_Multi_Emin::Evaluate
(const Particle_Vector& pvec, double weight, double ncount)
{
  if (m_on) {
    size_t n(0);
    for (size_t i=0;i<pvec.size();++i) if (pvec[i]->Momentum()[0]>m_cutoff) ++n;
    p_histo->Insert((double)n,weight,ncount);
  }
}

Primitive_Observable_Base * Conditional_One_Particle_Multi_Emin::Copy() const
{
  return new Conditional_One_Particle_Multi_Emin(m_flavour,m_cutoff,m_type,
                                                 m_xmin,m_xmax,m_nbins,
                                                 m_listname);
}

