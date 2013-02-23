#include "AddOns/Analysis/Observables/Scaled_Observables.H"
#include "AddOns/Analysis/Main/Primitive_Analysis.H"

using namespace ANALYSIS;

#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Run_Parameter.H"

template <class Class>
Primitive_Observable_Base *GetObservable(const Argument_Matrix &parameters)
{									
  if (parameters.size()<1) return NULL;
  if (parameters.size()==1) {
    if (parameters[0].size()<4) return NULL;
    std::string list=parameters[0].size()>4?parameters[0][4]:finalstate_list;
    double ref=parameters[0].size()>5?ATOOLS::ToType<double>(parameters[0][5]):
      ATOOLS::rpa->gen.Ecms();
    return new Class(HistogramType(parameters[0][3]),
		     ATOOLS::ToType<double>(parameters[0][0]),
		     ATOOLS::ToType<double>(parameters[0][1]),
		     ATOOLS::ToType<int>(parameters[0][2]),list,ref);
  }
  else if (parameters.size()<4) return NULL;
  double min=0.0, max=1.0, ref=ATOOLS::rpa->gen.Ecms();
  size_t bins=100;
  std::string list=finalstate_list, scale="Lin";
  for (size_t i=0;i<parameters.size();++i) {
    if (parameters[i].size()<2) continue;
    if (parameters[i][0]=="MIN") min=ATOOLS::ToType<double>(parameters[i][1]);
    else if (parameters[i][0]=="MAX") max=ATOOLS::ToType<double>(parameters[i][1]);
    else if (parameters[i][0]=="REF") ref=ATOOLS::ToType<double>(parameters[i][1]);
    else if (parameters[i][0]=="BINS") bins=ATOOLS::ToType<int>(parameters[i][1]);
    else if (parameters[i][0]=="SCALE") scale=parameters[i][1];
    else if (parameters[i][0]=="LIST") list=parameters[i][1];
  }
  return new Class(HistogramType(scale),min,max,bins,list,ref);
}									

#define DEFINE_GETTER_METHOD(CLASS,NAME)				\
  Primitive_Observable_Base *					\
  ATOOLS::Getter<Primitive_Observable_Base,Argument_Matrix,CLASS>::operator()(const Argument_Matrix &parameters) const \
  { return GetObservable<CLASS>(parameters); }

#define DEFINE_PRINT_METHOD(NAME)					\
  void ATOOLS::Getter<Primitive_Observable_Base,Argument_Matrix,NAME>::PrintInfo(std::ostream &str,const size_t width) const \
  { str<<"min max bins Lin|LinErr|Log|LogErr [list [ref]]"; }

#define DEFINE_OBSERVABLE_GETTER(CLASS,NAME,TAG)			\
  DECLARE_GETTER(CLASS,TAG,Primitive_Observable_Base,Argument_Matrix);	\
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
