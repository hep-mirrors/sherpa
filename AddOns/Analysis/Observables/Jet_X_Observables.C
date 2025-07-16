#include "AddOns/Analysis/Observables/Jet_Observables.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "AddOns/Analysis/Main/Primitive_Analysis.H"
#include "ATOOLS/Org/Shell_Tools.H"

using namespace ANALYSIS;

#include "ATOOLS/Org/MyStrStream.H"

template <class Class>
Primitive_Observable_Base *GetObservable(const Analysis_Key& key)
{
  ATOOLS::Scoped_Settings s{ key.m_settings };
  const auto min = s["Min"].SetDefault(0.0).Get<double>();
  const auto max = s["Max"].SetDefault(1.0).Get<double>();
  const auto bins = s["Bins"].SetDefault(100).Get<size_t>();
  const auto nmin = s["NMin"].SetDefault(  1).Get<size_t>();
  const auto nmax = s["NMax"].SetDefault( 10).Get<size_t>();
  const auto mode = s["Mode"].SetDefault(  1).Get<size_t>();
  const auto list = s["List"].SetDefault(std::string(finalstate_list)).Get<std::string>();
  const auto reflist = s["RefList"].SetDefault("").Get<std::string>();
  const auto scale = s["Scale"].SetDefault("Lin").Get<std::string>();
  return new Class(HistogramType(scale),min,max,bins,mode,nmin,nmax,list,reflist);
}

#define DEFINE_GETTER_METHOD(CLASS,NAME)				\
  Primitive_Observable_Base *					\
  ATOOLS::Getter<Primitive_Observable_Base,Analysis_Key,CLASS>::operator()(const Analysis_Key& key) const \
  { return GetObservable<CLASS>(key); }

#define DEFINE_PRINT_METHOD(NAME)					\
  void ATOOLS::Getter<Primitive_Observable_Base,Analysis_Key,NAME>::PrintInfo(std::ostream &str,const size_t width) const \
  { str<<"e.g. {Min: 1, Max: 10, Bins: 100, NMin: 1, NMax: 10, Mode: 1, Scale: Lin, List: FinalState, RefList: <list>}"; }

#define DEFINE_OBSERVABLE_GETTER(CLASS,NAME,TAG)			\
  DECLARE_GETTER(CLASS,TAG,Primitive_Observable_Base,Analysis_Key);	\
  DEFINE_GETTER_METHOD(CLASS,NAME)					\
  DEFINE_PRINT_METHOD(CLASS)

using namespace ATOOLS;

  class Jet_X_Observables : public Jet_Observable_Base {
  protected:
    std::string m_reflistname;
  public:
    Jet_X_Observables(unsigned int type,double xmin,double xmax,int nbins,
			   unsigned int mode, unsigned int minn,unsigned int maxn, 
			   const std::string & =std::string("FinalState"),
			   const std::string & =std::string(""));
    void Evaluate(const ATOOLS::Blob_List & blobs,double weight, double ncount);
    void EvaluateNLOcontrib(double weight, double ncount);
    double Calc(const Particle * p);
    virtual double Calc(const Vec4D &m1,const Vec4D &m2) =0;
  };


Jet_X_Observables::Jet_X_Observables(unsigned int type,double xmin,double xmax,int nbins,
					       unsigned int mode,unsigned int minn,unsigned int maxn, 
					       const std::string & listname,
					       const std::string & reflistname) :
  Jet_Observable_Base(type,xmin,xmax,nbins,mode,minn,maxn,listname) 
{
  if (reflistname=="") {
    m_reflistname = listname;
    m_name=listname+"_";
  }
  else {
    m_reflistname = reflistname;
    m_name=listname+"_"+reflistname+"_";
  } 
  if (m_minn!=0) {
    MyStrStream str;
    str<<m_name<<m_mode<<"_"<<m_minn<<"_";
    str>>m_name;
  }
}

double Jet_X_Observables::Calc(const Particle * p)
{
  return 0.;
}

void Jet_X_Observables::Evaluate(const Blob_List & blobs,double weight, double ncount)
{
  Particle_List * pl=p_ana->GetParticleList(m_listname);
  if ((m_mode==1 && pl->size()>=m_minn) ||
      (m_mode==2 && pl->size()==m_minn)) {
    Particle_List * pl2=p_ana->GetParticleList(m_reflistname);
    // fill
    size_t i=1;
    m_histos[0]->Insert(0.,0.,ncount);
    Vec4D rmom(0.,0.,0.,0.);
    for (Particle_List::const_iterator it=pl2->begin();it!=pl2->end();++it) {
      rmom+=(*it)->Momentum();
    }

    for (Particle_List::const_iterator it=pl->begin();it!=pl->end();++it,++i) {
      double value=Calc((*it)->Momentum(),rmom);
      m_histos[0]->Insert(value,weight,0);
      if (i<=m_maxn) m_histos[i]->Insert(value,weight,ncount);
    }
    for (; i<m_histos.size();++i) { 
      m_histos[i]->Insert(0.,0.,ncount);
    }
  }
  else {
    // fill with 0
    m_histos[0]->Insert(0.,0.,ncount);
    for (size_t i=1; i<m_histos.size();++i) {
      m_histos[i]->Insert(0.,0.,ncount);
    }
  }
}

void Jet_X_Observables::EvaluateNLOcontrib(double weight, double ncount)
{
  Particle_List * pl=p_ana->GetParticleList(m_listname);
  if ((m_mode==1 && pl->size()>=m_minn) ||
      (m_mode==2 && pl->size()==m_minn)) {
    Particle_List * pl2=p_ana->GetParticleList(m_reflistname);
    // fill
    size_t i=1;
    m_histos[0]->InsertMCB(0.,0.,ncount);
    Vec4D rmom(0.,0.,0.,0.);
    for (Particle_List::const_iterator it=pl2->begin();it!=pl2->end();++it)
      rmom+=(*it)->Momentum();

    for (Particle_List::const_iterator it=pl->begin();it!=pl->end();++it,++i) {
      double value=Calc((*it)->Momentum(),rmom);
      m_histos[0]->InsertMCB(value,weight,0);
      if (i<=m_maxn) m_histos[i]->InsertMCB(value,weight,ncount);
    }
    for (; i<m_histos.size();++i) { 
      m_histos[i]->InsertMCB(0.,0.,ncount);
    }
  }
  else {
    // fill with 0
    m_histos[0]->InsertMCB(0.,0.,ncount);
    for (size_t i=1; i<m_histos.size();++i) {
      m_histos[i]->InsertMCB(0.,0.,ncount);
    }
  }
}


////////////////////////////////////////////////////////////////////////////////

  class JetX_DeltaR_Distribution : public Jet_X_Observables {
  protected:
  public:
    JetX_DeltaR_Distribution(unsigned int type,double xmin,double xmax,int nbins,
			      unsigned int mode, unsigned int minn, unsigned int maxn, 
			      const std::string & =std::string("Jets"),
			      const std::string & =std::string("FinalState"));

    Primitive_Observable_Base * Copy() const;
    double Calc(const Vec4D &m1,const Vec4D &m2);
  };

DEFINE_OBSERVABLE_GETTER(JetX_DeltaR_Distribution,
			 JetX_DeltaR_Distribution_Getter,"JetXDR")

JetX_DeltaR_Distribution::JetX_DeltaR_Distribution(unsigned int type,double xmin,double xmax,int nbins,
						    unsigned int mode,unsigned int minn,unsigned int maxn, 
						    const std::string & lname,const std::string & rname) :
  Jet_X_Observables(type,xmin,xmax,nbins,mode,minn,maxn,lname,rname)
{
  m_name+="dR2_";
}

Primitive_Observable_Base * JetX_DeltaR_Distribution::Copy() const 
{
  JetX_DeltaR_Distribution * jdr =
    new JetX_DeltaR_Distribution(m_type,m_xmin,m_xmax,m_nbins,m_mode,m_minn,m_maxn,m_listname,m_reflistname);
  return jdr;
}

double JetX_DeltaR_Distribution::Calc(const Vec4D &mom1,const Vec4D &mom2)
{
  double pt1  = mom1.PPerp();
  double pt2  = mom2.PPerp();
  double dphi = acos((mom1[1]*mom2[1]+mom1[2]*mom2[2])/(pt1*pt2));
  double deta = mom1.Eta()-mom2.Eta();
  return sqrt(sqr(deta) + sqr(dphi)); 
}

//----------------------------------------------------------------------

  class JetX_DeltaEta_Distribution : public Jet_X_Observables {
  protected:
  public:
    JetX_DeltaEta_Distribution(unsigned int type,double xmin,double xmax,int nbins,
			      unsigned int mode, unsigned int minn, unsigned int maxn, 
			      const std::string & =std::string("Jets"),
			      const std::string & =std::string("FinalState"));

    Primitive_Observable_Base * Copy() const;
    double Calc(const Vec4D &m1,const Vec4D &m2);
  };

DEFINE_OBSERVABLE_GETTER(JetX_DeltaEta_Distribution,
			 JetX_DeltaEta_Distribution_Getter,"JetXDEta")

JetX_DeltaEta_Distribution::JetX_DeltaEta_Distribution(unsigned int type,double xmin,double xmax,int nbins,
							 unsigned int mode,unsigned int minn,unsigned int maxn, 
							 const std::string & lname,const std::string & rname) :
  Jet_X_Observables(type,xmin,xmax,nbins,mode,minn,maxn,lname,rname)
{
  m_name+="deta2_";
}

Primitive_Observable_Base * JetX_DeltaEta_Distribution::Copy() const 
{
  JetX_DeltaEta_Distribution * jde =
    new JetX_DeltaEta_Distribution(m_type,m_xmin,m_xmax,m_nbins,m_mode,m_minn,m_maxn,m_listname,m_reflistname);
  return jde;
}

double JetX_DeltaEta_Distribution::Calc(const Vec4D &mom1,const Vec4D &mom2)
{
  return dabs((mom1.Eta()-mom2.Eta()));
}
//----------------------------------------------------------------------


  class JetX_DeltaY_Distribution : public Jet_X_Observables {
  protected:
  public:
    JetX_DeltaY_Distribution(unsigned int type,double xmin,double xmax,int nbins,
			     unsigned int mode, unsigned int minn, unsigned int maxn, 
			     const std::string & =std::string("Jets"),
			     const std::string & =std::string("FinalState"));

    Primitive_Observable_Base * Copy() const;
    double Calc(const Vec4D &m1,const Vec4D &m2);
  };

DEFINE_OBSERVABLE_GETTER(JetX_DeltaY_Distribution,
			 JetX_DeltaY_Distribution_Getter,"JetXDY")

JetX_DeltaY_Distribution::JetX_DeltaY_Distribution(unsigned int type,double xmin,double xmax,int nbins,
						   unsigned int mode,unsigned int minn,unsigned int maxn, 
						   const std::string & lname,const std::string & rname) :
  Jet_X_Observables(type,xmin,xmax,nbins,mode,minn,maxn,lname,rname)
{
  m_name+="dY_";
}

Primitive_Observable_Base * JetX_DeltaY_Distribution::Copy() const 
{
  JetX_DeltaY_Distribution * jde =
    new JetX_DeltaY_Distribution(m_type,m_xmin,m_xmax,m_nbins,m_mode,m_minn,m_maxn,m_listname,m_reflistname);
  return jde;
}

double JetX_DeltaY_Distribution::Calc(const Vec4D &mom1,const Vec4D &mom2)
{
  return dabs((mom1.Y()-mom2.Y()));
}
//----------------------------------------------------------------------


  class JetX_DeltaPhi_Distribution : public Jet_X_Observables {
  protected:
  public:
    JetX_DeltaPhi_Distribution(unsigned int type,double xmin,double xmax,int nbins,
			       unsigned int mode, unsigned int minn, unsigned int maxn, 
			       const std::string & =std::string("Jets"),
			       const std::string & =std::string("FinalState"));
    
    Primitive_Observable_Base * Copy() const;
    double Calc(const Vec4D &m1,const Vec4D &m2);
  };

DEFINE_OBSERVABLE_GETTER(JetX_DeltaPhi_Distribution,
			 JetX_DeltaPhi_Distribution_Getter,"JetXDPhi")
  
JetX_DeltaPhi_Distribution::JetX_DeltaPhi_Distribution(unsigned int type,double xmin,double xmax,int nbins,
						       unsigned int mode,unsigned int minn,unsigned int maxn, 
						       const std::string & lname,const std::string & rname) :
  Jet_X_Observables(type,xmin,xmax,nbins,mode,minn,maxn,lname,rname)
{
  m_name+="dPhi_";
}

Primitive_Observable_Base * JetX_DeltaPhi_Distribution::Copy() const 
{
  JetX_DeltaPhi_Distribution * jde =
    new JetX_DeltaPhi_Distribution(m_type,m_xmin,m_xmax,m_nbins,m_mode,m_minn,m_maxn,m_listname,m_reflistname);
  return jde;
}

double JetX_DeltaPhi_Distribution::Calc(const Vec4D &mom1,const Vec4D &mom2)
{
  return mom1.DPhi(mom2);
}
//----------------------------------------------------------------------

