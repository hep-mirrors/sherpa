#include "AddOns/Analysis/Observables/Scaled_Observables.H"
#include "AddOns/Analysis/Main/Primitive_Analysis.H"

using namespace ANALYSIS;

#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Run_Parameter.H"

template <class Class>
Primitive_Observable_Base *GetObservable(const Analysis_Key& key)
{
  ATOOLS::Scoped_Settings s{ key.m_settings };
  const auto min = s["Min"].SetDefault(0.0).Get<double>();
  const auto max = s["Max"].SetDefault(1.0).Get<double>();
  const auto ref = s["Ref"].SetDefault(ATOOLS::rpa->gen.Ecms()).Get<double>();
  const auto bins = s["Bins"].SetDefault(100).Get<size_t>();
  const auto scale = s["Scale"].SetDefault("Lin").Get<std::string>();
  const auto list = s["List"]
    .SetDefault(std::string(finalstate_list))
    .Get<std::string>();
  return new Class(HistogramType(scale),min,max,bins,list,ref);
}

#define DEFINE_GETTER_METHOD(CLASS,NAME)				\
  Primitive_Observable_Base *					\
  ATOOLS::Getter<Primitive_Observable_Base,Analysis_Key,CLASS>::operator()(const Analysis_Key& key) const \
  { return GetObservable<CLASS>(key); }

#define DEFINE_PRINT_METHOD(NAME)					\
  void ATOOLS::Getter<Primitive_Observable_Base,Analysis_Key,NAME>::PrintInfo(std::ostream &str,const size_t width) const \
  { str<<"e.g. {Min: 1, Max: 10, Bins: 100, Scale: Lin, List: FinalState, Ref: <energy>}"; }

#define DEFINE_OBSERVABLE_GETTER(CLASS,NAME,TAG)			\
  DECLARE_GETTER(CLASS,TAG,Primitive_Observable_Base,Analysis_Key);	\
  DEFINE_GETTER_METHOD(CLASS,NAME)					\
  DEFINE_PRINT_METHOD(CLASS)

using namespace ATOOLS;
using namespace std;

Scaled_Observable_Base::Scaled_Observable_Base(int type,double xmin,double xmax,int nbins,
					       const std::string & listname, const std::string & name,
					       double ecms) :
  Primitive_Observable_Base(type,xmin,xmax,nbins), m_ecms(ecms)
{
  m_name=listname+"_"+name+".dat";

  if (listname!=std::string("")) m_listname = listname;
  m_blobtype = std::string("");
  m_blobdisc = false;
}

void Scaled_Observable_Base::Evaluate(double value,double weight, double ncount) 
{
  p_histo->Insert(value,weight,ncount); 
}

 
void Scaled_Observable_Base::Evaluate(int nout,const ATOOLS::Vec4D * moms,
					    double weight, double ncount) 
{
  for (int i=0;i<nout;i++) Evaluate(moms[i],weight,ncount);
}


void Scaled_Observable_Base::Evaluate(const Particle_List & plist,double weight,double ncount )
{
  for (Particle_List::const_iterator plit=plist.begin();plit!=plist.end();++plit) {
    Evaluate((*plit)->Momentum(),weight, ncount);
  }
}



//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DEFINE_OBSERVABLE_GETTER(Scaled_Momentum,Scaled_Momentum_Getter,"XP")

Scaled_Momentum::Scaled_Momentum(int type,double xmin,double xmax,int nbins,
				 const std::string & listname, double ecms) :
  Scaled_Observable_Base(type,xmin,xmax,nbins,listname,"ScaledMomentum",ecms) { }


void Scaled_Momentum::Evaluate(const Vec4D & mom,double weight,double ncount) 
{
  double xp = 2.*Vec3D(mom).Abs()/m_ecms;

  p_histo->Insert(xp,weight,ncount); 
} 

Primitive_Observable_Base * Scaled_Momentum::Copy() const
{
  return new Scaled_Momentum(m_type,m_xmin,m_xmax,m_nbins,m_listname,m_ecms);
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DEFINE_OBSERVABLE_GETTER(Log_Scaled_Momentum,Log_Scaled_Momentum_Getter,"LogXP")

Log_Scaled_Momentum::Log_Scaled_Momentum(int type,double xmin,double xmax,int nbins,
					 const std::string & listname, double ecms) :
  Scaled_Observable_Base(type,xmin,xmax,nbins,listname,"LogScaledMomentum", ecms) { }


void Log_Scaled_Momentum::Evaluate(const Vec4D & mom,double weight,double ncount) 
{
  double xp = 2.*Vec3D(mom).Abs()/m_ecms;
  double xi = - log(xp);

  p_histo->Insert(xi,weight,ncount); 
} 

Primitive_Observable_Base * Log_Scaled_Momentum::Copy() const
{
  return new Log_Scaled_Momentum(m_type,m_xmin,m_xmax,m_nbins,m_listname,m_ecms);
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DEFINE_OBSERVABLE_GETTER(Scaled_Energy,Scaled_Energy_Getter,"XE")

Scaled_Energy::Scaled_Energy(int type,double xmin,double xmax,int nbins,
			     const std::string & listname, double ecms) :
  Scaled_Observable_Base(type,xmin,xmax,nbins,listname,"ScaledEnergy",ecms) { }


void Scaled_Energy::Evaluate(const Vec4D & mom,double weight, double ncount) 
{
  double E = 2.*mom[0]/m_ecms;
  p_histo->Insert(E,weight,ncount); 
} 

Primitive_Observable_Base * Scaled_Energy::Copy() const
{
  return new Scaled_Energy(m_type,m_xmin,m_xmax,m_nbins,m_listname,m_ecms);
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DEFINE_OBSERVABLE_GETTER(EtaTracks,EtaTracks_Getter,"EtaTracks")

EtaTracks::EtaTracks(int type,double xmin,double xmax,int nbins,
			   const std::string & listname, double ecms) :
  Scaled_Observable_Base(type,xmin,xmax,nbins,listname,"EtaTracks",ecms) { }


void EtaTracks::Evaluate(const Vec4D & mom,double weight,double ncount) 
{
  
  double eta = 0.;
  eta=mom.Eta();
  
  if (eta<0.) {
    p_histo->Insert(eta,weight,ncount); 
  }
  else {
    p_histo->Insert(eta,weight,ncount); 
  }
} 

Primitive_Observable_Base * EtaTracks::Copy() const
{
  return new EtaTracks(m_type,m_xmin,m_xmax,m_nbins,m_listname,m_ecms);
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DEFINE_OBSERVABLE_GETTER(EtaTracksAsym,EtaTracksAsym_Getter,"EtaTracksAsym")

EtaTracksAsym::EtaTracksAsym(int type,double xmin,double xmax,int nbins,
			 const std::string & listname, double ecms) :
  Scaled_Observable_Base(type,xmin,xmax,nbins,listname,"EtaTracksAsym",ecms) { }


void EtaTracksAsym::Evaluate(const Vec4D & mom,double weight,double ncount) 
{
  
  double eta = 0.;
  eta=mom.Eta();
  
  if (eta<0.) {
    p_histo->Insert(-eta,-weight,ncount); 
  }
  else {
    p_histo->Insert(eta,weight,ncount); 
  }
} 

Primitive_Observable_Base * EtaTracksAsym::Copy() const
{
  return new EtaTracksAsym(m_type,m_xmin,m_xmax,m_nbins,m_listname,m_ecms);
}
