#include "ATOOLS/Org/MyStrStream.H"
#include "AddOns/Analysis/Main/Primitive_Analysis.H"
#include "AddOns/Analysis/Observables/Primitive_Observable_Base.H"
#include "ATOOLS/Org/Shell_Tools.H"
#include "ATOOLS/Math/Poincare.H"

using namespace ANALYSIS;

#include "ATOOLS/Org/MyStrStream.H"

template <class Class>
Primitive_Observable_Base* GetObservable(const Analysis_Key& key)
{
  ATOOLS::Scoped_Settings s{ key.m_settings };
  const auto min = s["Min"].SetDefault(0.0).Get<double>();
  const auto max = s["Max"].SetDefault(1.0).Get<double>();
  const auto bins = s["Bins"].SetDefault(100).Get<size_t>();
  const auto scale = s["Scale"].SetDefault("Lin").Get<std::string>();
  const auto list = s["List"].SetDefault(std::string(finalstate_list)).Get<std::string>();
  const auto reflist = s["Ref"].SetDefault("").Get<std::string>();
  std::vector<ATOOLS::Flavour> flavs;
  flavs.reserve(2);
  for (size_t i{0}; i < 2; ++i) {
    const auto flavparamname = "Flav" + ATOOLS::ToString(i + 1);
    if (!s[flavparamname].IsCustomised())
      THROW(missing_input, flavparamname + "must be set.");
    const auto rawflav = s[flavparamname].SetDefault(0).GetScalar<int>();
    flavs.push_back(ATOOLS::Flavour((kf_code)std::abs(rawflav)));
    if (rawflav < 0)
      flavs.back() = flavs.back().Bar();
  }
  return new Class(flavs[0], flavs[1],
                   HistogramType(scale), min, max, bins,
                   list, reflist);
}

#define DEFINE_GETTER_METHOD(CLASS)				\
  Primitive_Observable_Base *					\
  ATOOLS::Getter<Primitive_Observable_Base,Analysis_Key,CLASS>::operator()(const Analysis_Key& key) const \
  { return GetObservable<CLASS>(key); }

#define DEFINE_PRINT_METHOD(NAME)					\
  void ATOOLS::Getter<Primitive_Observable_Base,Analysis_Key,NAME>::PrintInfo(std::ostream &str,const size_t width) const \
  { str<<"e.g. {Flav1: kf1, Flav2: kf2, Min: 0, Max: 1, Bins: 100, Scale: Lin, List: <list>, Ref: <list>}"; }

#define DEFINE_OBSERVABLE_GETTER(CLASS,TAG)			\
  DECLARE_GETTER(CLASS,TAG,Primitive_Observable_Base,Analysis_Key);	\
  DEFINE_GETTER_METHOD(CLASS)					\
  DEFINE_PRINT_METHOD(CLASS)

using namespace ATOOLS;

  class EV_SC_Observables : public Primitive_Observable_Base {
  protected:
    std::string m_reflistname;
    ATOOLS::Flavour      m_flav1,m_flav2;
  public:
    EV_SC_Observables(const ATOOLS::Flavour & flav1, const ATOOLS::Flavour & flav2,
		     unsigned int type,double xmin,double xmax,int nbins,
		     const std::string &,
		     const std::string &);
    void Evaluate(const ATOOLS::Blob_List & blobs,double weight, double ncount);
    void EvaluateNLOcontrib(double weight, double ncount);
    void EvaluateNLOevt();
    double CalcSCth(const Vec4D& m1,const Vec4D& m2);
    virtual double Calc(const Vec4D &m) =0;
  };


EV_SC_Observables::EV_SC_Observables(const Flavour & flav1,const Flavour & flav2,
				   unsigned int type,double xmin,double xmax,int nbins,
				   const std::string & listname,
				   const std::string & reflistname) :
  Primitive_Observable_Base(type,xmin,xmax,nbins), 
  m_flav1(flav1), m_flav2(flav2)
{
  m_listname=listname;
  m_reflistname=reflistname;
  m_name="EV_SC_"+listname+"_"+reflistname+"_";
}

void EV_SC_Observables::Evaluate(const Blob_List & blobs,double weight, double ncount)
{
  double cth(0.);
  Particle_List * rlist=p_ana->GetParticleList(m_reflistname);
  for (Particle_List::const_iterator plit1=rlist->begin();plit1!=rlist->end();++plit1) {
    if ((*plit1)->Flav()==m_flav1) {
      for (Particle_List::const_iterator plit2=rlist->begin();plit2!=rlist->end();++plit2) {
	if ((*plit2)->Flav()==m_flav2 && plit1!=plit2) {
	  cth=CalcSCth((*plit1)->Momentum(),(*plit2)->Momentum());
	}
      }
    }
  }

  Particle_List * pl=p_ana->GetParticleList(m_listname);
  Vec4D fmom(0.,0.,0.,0.);
  for (Particle_List::const_iterator it=pl->begin();it!=pl->end();++it)
    fmom+=(*it)->Momentum();
  
  p_histo->Insert(Calc(fmom),cth*weight,ncount); 
}


void EV_SC_Observables::EvaluateNLOcontrib(double weight, double ncount)
{
  double cth(0.);
  Particle_List * rlist=p_ana->GetParticleList(m_reflistname);
  for (Particle_List::const_iterator plit1=rlist->begin();plit1!=rlist->end();++plit1) {
    if ((*plit1)->Flav()==m_flav1) {
      for (Particle_List::const_iterator plit2=rlist->begin();plit2!=rlist->end();++plit2) {
	if ((*plit2)->Flav()==m_flav2 && plit1!=plit2) {
	  cth=CalcSCth((*plit1)->Momentum(),(*plit2)->Momentum());
	}
      }
    }
  }

  Particle_List * pl=p_ana->GetParticleList(m_listname);
  Vec4D fmom(0.,0.,0.,0.);
  for (Particle_List::const_iterator it=pl->begin();it!=pl->end();++it)
    fmom+=(*it)->Momentum();
  
  p_histo->InsertMCB(Calc(fmom),cth*weight,ncount); 
}
 
double EV_SC_Observables::CalcSCth(const Vec4D& m1,const Vec4D& m2) {
  if (m1[0]>m2[0]) return 1.;
  if (m1[0]<m2[0]) return -1.;
  return 0.; 
}

void EV_SC_Observables::EvaluateNLOevt()
{
  p_histo->FinishMCB();
}

////////////////////////////////////////////////////////////////////////////////

 class EV_SC_ET : public EV_SC_Observables {
 protected:
 public:
   EV_SC_ET(const Flavour & flav1,const Flavour & flav2,
	   unsigned int type,double xmin,double xmax,int nbins,
	   const std::string & listname,
	   const std::string & reflistname);
   
   Primitive_Observable_Base * Copy() const;
   double Calc(const Vec4D &m1);
 };

DEFINE_OBSERVABLE_GETTER(EV_SC_ET,"EVSC_ET")

  EV_SC_ET::EV_SC_ET(const Flavour & flav1,const Flavour & flav2,
		   unsigned int type,double xmin,double xmax,int nbins,
		   const std::string & lname,
		   const std::string & rname) :
  EV_SC_Observables(flav1,flav2,type,xmin,xmax,nbins,lname,rname)
{
  m_name+="ET.dat";
}

Primitive_Observable_Base * EV_SC_ET::Copy() const 
{
  EV_SC_ET * cpo =
    new EV_SC_ET(m_flav1,m_flav2,m_type,m_xmin,m_xmax,m_nbins,m_listname,m_reflistname);
  return cpo;
}

double EV_SC_ET::Calc(const Vec4D &mom)
{
  return mom.EPerp();
}

//----------------------------------------------------------------------

  class EV_SC_PT : public EV_SC_Observables {
  protected:
  public:
    EV_SC_PT(const Flavour & flav1,const Flavour & flav2,
	    unsigned int type,double xmin,double xmax,int nbins,
	    const std::string & listname,
	    const std::string & reflistname);

    Primitive_Observable_Base * Copy() const;
    double Calc(const Vec4D &m1);
  };

DEFINE_OBSERVABLE_GETTER(EV_SC_PT,"EVSC_PT")

  EV_SC_PT::EV_SC_PT(const Flavour & flav1,const Flavour & flav2,
		   unsigned int type,double xmin,double xmax,int nbins,
		   const std::string & lname,
		   const std::string & rname) :
  EV_SC_Observables(flav1,flav2,type,xmin,xmax,nbins,lname,rname)
{
  m_name+="PT.dat";
}

Primitive_Observable_Base * EV_SC_PT::Copy() const 
{
  EV_SC_PT * cpo =
    new EV_SC_PT(m_flav1,m_flav2,m_type,m_xmin,m_xmax,m_nbins,m_listname,m_reflistname);
  return cpo;
}

double EV_SC_PT::Calc(const Vec4D &mom)
{
  return mom.PPerp();
}


//----------------------------------------------------------------------

  class EV_SC_Eta : public EV_SC_Observables {
  protected:
  public:
    EV_SC_Eta(const Flavour & flav1,const Flavour & flav2,
	     unsigned int type,double xmin,double xmax,int nbins,
	     const std::string & listname,
	     const std::string & reflistname);

    Primitive_Observable_Base * Copy() const;
    double Calc(const Vec4D &m1);
  };

DEFINE_OBSERVABLE_GETTER(EV_SC_Eta,"EVSC_Eta")

  EV_SC_Eta::EV_SC_Eta(const Flavour & flav1,const Flavour & flav2,
		     unsigned int type,double xmin,double xmax,int nbins,
		     const std::string & lname,
		     const std::string & rname) :
  EV_SC_Observables(flav1,flav2,type,xmin,xmax,nbins,lname,rname)
{
  m_name+="Eta.dat";
}

Primitive_Observable_Base * EV_SC_Eta::Copy() const 
{
  EV_SC_Eta * cpo =
    new EV_SC_Eta(m_flav1,m_flav2,m_type,m_xmin,m_xmax,m_nbins,m_listname,m_reflistname);
  return cpo;
}

double EV_SC_Eta::Calc(const Vec4D &mom)
{
  return mom.Eta();
}

//----------------------------------------------------------------------


  class EV_SC_Y : public EV_SC_Observables {
  protected:
  public:
    EV_SC_Y(const Flavour & flav1,const Flavour & flav2,
	   unsigned int type,double xmin,double xmax,int nbins,
	   const std::string & listname,
	   const std::string & reflistname);
    
    Primitive_Observable_Base * Copy() const;
    double Calc(const Vec4D &m1);
  };
 
DEFINE_OBSERVABLE_GETTER(EV_SC_Y,"EVSC_Y")
  
  EV_SC_Y::EV_SC_Y(const Flavour & flav1,const Flavour & flav2,
		 unsigned int type,double xmin,double xmax,int nbins,
		 const std::string & lname,
		 const std::string & rname) :
  EV_SC_Observables(flav1,flav2,type,xmin,xmax,nbins,lname,rname)
{
  m_name+="Y.dat";
}

Primitive_Observable_Base * EV_SC_Y::Copy() const 
{
  EV_SC_Y * cpo =
    new EV_SC_Y(m_flav1,m_flav2,m_type,m_xmin,m_xmax,m_nbins,m_listname,m_reflistname);
  return cpo;
}

double EV_SC_Y::Calc(const Vec4D &mom)
{
  return mom.Y();
}

//----------------------------------------------------------------------


