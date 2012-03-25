#include "ATOOLS/Org/MyStrStream.H"
#include "AddOns/Analysis/Main/Primitive_Analysis.H"
#include "AddOns/Analysis/Observables/Primitive_Observable_Base.H"
#include "ATOOLS/Org/Shell_Tools.H"
#include "ATOOLS/Math/Poincare.H"

using namespace ANALYSIS;

#include "ATOOLS/Org/MyStrStream.H"

template <class Class>
Primitive_Observable_Base *const GetObservable(const Argument_Matrix &parameters)
{									
  if (parameters.size()<1) return NULL;
  if (parameters.size()==1) {
    if (parameters[0].size()<8) return NULL;
    ATOOLS::Flavour f[2];
    for (short unsigned int i=0;i<2;++i) {
      int kf=ATOOLS::ToType<int>(parameters[0][i]);
      f[i]=ATOOLS::Flavour((kf_code)abs(kf));
      if (kf<0) f[i]=f[i].Bar();
    }
    std::string list=parameters[0][6];
    std::string reflist=parameters[0][7];
    return new Class(f[0],f[1],HistogramType(parameters[0][5]),
		     ATOOLS::ToType<double>(parameters[0][2]),
		     ATOOLS::ToType<double>(parameters[0][3]),
		     ATOOLS::ToType<int>(parameters[0][4]),list,reflist);
  }
  else if (parameters.size()<8) return NULL;
  double min=0.0, max=1.0;
  size_t bins=100;
  std::string list=finalstate_list, scale="Lin";
  std::string reflist="";
  ATOOLS::Flavour f[2];
  for (size_t i=0;i<parameters.size();++i) {
    if (parameters[i].size()<2) continue;
    for (short unsigned int j=0;j<2;++j) {
      if (parameters[i][0]==std::string("FLAV")+ATOOLS::ToString(j+1)) {
	int kf=ATOOLS::ToType<int>(parameters[i][1]);
	f[j]=ATOOLS::Flavour((kf_code)abs(kf));
	if (kf<0) f[j]=f[j].Bar();
      }
    }
    if (parameters[i].size()<2) continue;
    if (parameters[i][0]=="MIN") min=ATOOLS::ToType<double>(parameters[i][1]);
    else if (parameters[i][0]=="MAX") max=ATOOLS::ToType<double>(parameters[i][1]);
    else if (parameters[i][0]=="BINS") bins=ATOOLS::ToType<int>(parameters[i][1]);
    else if (parameters[i][0]=="SCALE") scale=parameters[i][1];
    else if (parameters[i][0]=="LIST") list=parameters[i][1];
    else if (parameters[i][0]=="REF")  reflist=parameters[i][1];
  }
  return new Class(f[0],f[1],HistogramType(scale),min,max,bins,list,reflist);
} 

#define DEFINE_GETTER_METHOD(CLASS,NAME)				\
  Primitive_Observable_Base *					\
  NAME::operator()(const Argument_Matrix &parameters) const		\
  { return GetObservable<CLASS>(parameters); }

#define DEFINE_PRINT_METHOD(NAME)					\
  void NAME::PrintInfo(std::ostream &str,const size_t width) const	\
  { str<<"kf1 kf2 min max bins Lin|LinErr|Log|LogErr list reflist"; }

#define DEFINE_OBSERVABLE_GETTER(CLASS,NAME,TAG)			\
  DECLARE_GETTER(NAME,TAG,Primitive_Observable_Base,Argument_Matrix);	\
  DEFINE_GETTER_METHOD(CLASS,NAME)					\
  DEFINE_PRINT_METHOD(NAME)

using namespace ATOOLS;

  class EV_C_Observables : public Primitive_Observable_Base {
  protected:
    std::string m_reflistname;
    ATOOLS::Flavour      m_flav1,m_flav2;
  public:
    EV_C_Observables(const ATOOLS::Flavour & flav1, const ATOOLS::Flavour & flav2,
		     unsigned int type,double xmin,double xmax,int nbins,
		     const std::string &,
		     const std::string &);
    void Evaluate(const ATOOLS::Blob_List & blobs,double weight, double ncount);
    void EvaluateNLOcontrib(double weight, double ncount);
    void EvaluateNLOevt();
    double CalcCth(const Vec4D& m1,const Vec4D& m2);
    virtual double Calc(const Vec4D &m) =0;
  };


EV_C_Observables::EV_C_Observables(const Flavour & flav1,const Flavour & flav2,
				   unsigned int type,double xmin,double xmax,int nbins,
				   const std::string & listname,
				   const std::string & reflistname) :
  Primitive_Observable_Base(type,xmin,xmax,nbins), 
  m_flav1(flav1), m_flav2(flav2)
{
  m_listname=listname;
  m_reflistname=reflistname;
  m_name="EV_C_"+listname+"_"+reflistname+"_";
}

void EV_C_Observables::Evaluate(const Blob_List & blobs,double weight, double ncount)
{
  double cth(0.);
  Particle_List * rlist=p_ana->GetParticleList(m_reflistname);
  for (Particle_List::const_iterator plit1=rlist->begin();plit1!=rlist->end();++plit1) {
    if ((*plit1)->Flav()==m_flav1) {
      for (Particle_List::const_iterator plit2=rlist->begin();plit2!=rlist->end();++plit2) {
	if ((*plit2)->Flav()==m_flav2 && plit1!=plit2) {
	  cth=CalcCth((*plit1)->Momentum(),(*plit2)->Momentum());
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


void EV_C_Observables::EvaluateNLOcontrib(double weight, double ncount)
{
  double cth(0.);
  Particle_List * rlist=p_ana->GetParticleList(m_reflistname);
  for (Particle_List::const_iterator plit1=rlist->begin();plit1!=rlist->end();++plit1) {
    if ((*plit1)->Flav()==m_flav1) {
      for (Particle_List::const_iterator plit2=rlist->begin();plit2!=rlist->end();++plit2) {
	if ((*plit2)->Flav()==m_flav2 && plit1!=plit2) {
	  cth=CalcCth((*plit1)->Momentum(),(*plit2)->Momentum());
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
 
double EV_C_Observables::CalcCth(const Vec4D& m1,const Vec4D& m2) {
  Vec4D mother=m1+m2;
  Vec4D hm1(m1);
  Poincare boost(mother);
  boost.Boost(hm1);
  return mother.CosTheta(hm1);
}

void EV_C_Observables::EvaluateNLOevt()
{
  p_histo->FinishMCB();
}

////////////////////////////////////////////////////////////////////////////////

 class EV_C_ET : public EV_C_Observables {
 protected:
 public:
   EV_C_ET(const Flavour & flav1,const Flavour & flav2,
	   unsigned int type,double xmin,double xmax,int nbins,
	   const std::string & listname,
	   const std::string & reflistname);
   
   Primitive_Observable_Base * Copy() const;
   double Calc(const Vec4D &m1);
 };

DEFINE_OBSERVABLE_GETTER(EV_C_ET,EV_C_ET_Getter,"EVC_ET")

  EV_C_ET::EV_C_ET(const Flavour & flav1,const Flavour & flav2,
		   unsigned int type,double xmin,double xmax,int nbins,
		   const std::string & lname,
		   const std::string & rname) :
  EV_C_Observables(flav1,flav2,type,xmin,xmax,nbins,lname,rname)
{
  m_name+="ET.dat";
}

Primitive_Observable_Base * EV_C_ET::Copy() const 
{
  EV_C_ET * cpo =
    new EV_C_ET(m_flav1,m_flav2,m_type,m_xmin,m_xmax,m_nbins,m_listname,m_reflistname);
  return cpo;
}

double EV_C_ET::Calc(const Vec4D &mom)
{
  return mom.EPerp();
}

//----------------------------------------------------------------------

  class EV_C_PT : public EV_C_Observables {
  protected:
  public:
    EV_C_PT(const Flavour & flav1,const Flavour & flav2,
	    unsigned int type,double xmin,double xmax,int nbins,
	    const std::string & listname,
	    const std::string & reflistname);

    Primitive_Observable_Base * Copy() const;
    double Calc(const Vec4D &m1);
  };

DEFINE_OBSERVABLE_GETTER(EV_C_PT,EV_C_PT_Getter,"EVC_PT")

  EV_C_PT::EV_C_PT(const Flavour & flav1,const Flavour & flav2,
		   unsigned int type,double xmin,double xmax,int nbins,
		   const std::string & lname,
		   const std::string & rname) :
  EV_C_Observables(flav1,flav2,type,xmin,xmax,nbins,lname,rname)
{
  m_name+="PT.dat";
}

Primitive_Observable_Base * EV_C_PT::Copy() const 
{
  EV_C_PT * cpo =
    new EV_C_PT(m_flav1,m_flav2,m_type,m_xmin,m_xmax,m_nbins,m_listname,m_reflistname);
  return cpo;
}

double EV_C_PT::Calc(const Vec4D &mom)
{
  return mom.PPerp();
}


//----------------------------------------------------------------------

  class EV_C_Eta : public EV_C_Observables {
  protected:
  public:
    EV_C_Eta(const Flavour & flav1,const Flavour & flav2,
	     unsigned int type,double xmin,double xmax,int nbins,
	     const std::string & listname,
	     const std::string & reflistname);

    Primitive_Observable_Base * Copy() const;
    double Calc(const Vec4D &m1);
  };

DEFINE_OBSERVABLE_GETTER(EV_C_Eta,EV_C_Eta_Getter,"EVC_Eta")

  EV_C_Eta::EV_C_Eta(const Flavour & flav1,const Flavour & flav2,
		     unsigned int type,double xmin,double xmax,int nbins,
		     const std::string & lname,
		     const std::string & rname) :
  EV_C_Observables(flav1,flav2,type,xmin,xmax,nbins,lname,rname)
{
  m_name+="Eta.dat";
}

Primitive_Observable_Base * EV_C_Eta::Copy() const 
{
  EV_C_Eta * cpo =
    new EV_C_Eta(m_flav1,m_flav2,m_type,m_xmin,m_xmax,m_nbins,m_listname,m_reflistname);
  return cpo;
}

double EV_C_Eta::Calc(const Vec4D &mom)
{
  return mom.Eta();
}

//----------------------------------------------------------------------


  class EV_C_Y : public EV_C_Observables {
  protected:
  public:
    EV_C_Y(const Flavour & flav1,const Flavour & flav2,
	   unsigned int type,double xmin,double xmax,int nbins,
	   const std::string & listname,
	   const std::string & reflistname);
    
    Primitive_Observable_Base * Copy() const;
    double Calc(const Vec4D &m1);
  };
 
DEFINE_OBSERVABLE_GETTER(EV_C_Y,EV_C_Y_Getter,"EVC_Y")
  
  EV_C_Y::EV_C_Y(const Flavour & flav1,const Flavour & flav2,
		 unsigned int type,double xmin,double xmax,int nbins,
		 const std::string & lname,
		 const std::string & rname) :
  EV_C_Observables(flav1,flav2,type,xmin,xmax,nbins,lname,rname)
{
  m_name+="Y.dat";
}

Primitive_Observable_Base * EV_C_Y::Copy() const 
{
  EV_C_Y * cpo =
    new EV_C_Y(m_flav1,m_flav2,m_type,m_xmin,m_xmax,m_nbins,m_listname,m_reflistname);
  return cpo;
}

double EV_C_Y::Calc(const Vec4D &mom)
{
  return mom.Y();
}

//----------------------------------------------------------------------


