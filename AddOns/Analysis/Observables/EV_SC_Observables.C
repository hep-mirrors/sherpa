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

DEFINE_OBSERVABLE_GETTER(EV_SC_ET,EV_SC_ET_Getter,"EVSC_ET")

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

DEFINE_OBSERVABLE_GETTER(EV_SC_PT,EV_SC_PT_Getter,"EVSC_PT")

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

DEFINE_OBSERVABLE_GETTER(EV_SC_Eta,EV_SC_Eta_Getter,"EVSC_Eta")

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
 
DEFINE_OBSERVABLE_GETTER(EV_SC_Y,EV_SC_Y_Getter,"EVSC_Y")
  
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


