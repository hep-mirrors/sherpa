#include "AddOns/Analysis/Observables/PSM_Observables.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "AddOns/Analysis/Main/Primitive_Analysis.H"
#include "ATOOLS/Org/Shell_Tools.H"
#include "ATOOLS/Org/Message.H"

using namespace ANALYSIS;

#include "ATOOLS/Org/MyStrStream.H"

template <class Class>
Primitive_Observable_Base *GetObservable(const Analysis_Key& key)
{
  ATOOLS::Scoped_Settings s{ key.m_settings };

  const auto min = s["Min"].SetDefault(0.0).Get<double>();
  const auto max = s["Max"].SetDefault(1.0).Get<double>();
  const auto bins = s["Bins"].SetDefault(100).Get<size_t>();
  const auto scale = s["Scale"].SetDefault("Lin").Get<std::string>();
  const auto list = s["List"].SetDefault(std::string(finalstate_list)).Get<std::string>();

  std::vector<int> p;
  p.reserve(4);
  for (size_t i{0}; i < 4; ++i) {
    const auto pkey = "PN" + ATOOLS::ToString(i);
    p.push_back(s[pkey].SetDefault(-1).Get<int>());
  }

  return new Class(HistogramType(scale),min,max,bins,
                   p[0],p[1],p[2],p[3],
                   list);
}

#define DEFINE_GETTER_METHOD(CLASS)				\
  Primitive_Observable_Base *					\
  ATOOLS::Getter<Primitive_Observable_Base,Analysis_Key,CLASS>::operator()(const Analysis_Key& key) const \
  { return GetObservable<CLASS>(key); }

#define DEFINE_PRINT_METHOD(NAME)					\
  void ATOOLS::Getter<Primitive_Observable_Base,Analysis_Key,NAME>::PrintInfo(std::ostream &str,const size_t width) const \
  { str<<"e.g. {Min: 1, Max: 10, Bins: 100, PN0: p0, ..., PN3: p3, Scale: Lin, List: FinalState}"; }

#define DEFINE_OBSERVABLE_GETTER(CLASS,TAG)			\
  DECLARE_GETTER(CLASS,TAG,Primitive_Observable_Base,Analysis_Key);	\
  DEFINE_GETTER_METHOD(CLASS)					\
  DEFINE_PRINT_METHOD(CLASS)

using namespace ATOOLS;

DEFINE_OBSERVABLE_GETTER(PSM_Observable,"PSM")

PSM_Observable::PSM_Observable(unsigned int type,double xmin,double xmax,int nbins,
						 int p0,int p1,int p2, int p3,
						 const std::string & lname) :
  Primitive_Observable_Base(type,xmin,xmax,nbins)
{
  m_pnb.clear();
  m_pnb.push_back(p0);
  m_pnb.push_back(p1);
  m_pnb.push_back(p2);
  m_pnb.push_back(p3);
  m_listname = lname;
  m_name     = std::string("psm_");
  if (lname!=finalstate_list) m_name=lname+std::string("_")+m_name;
  if (m_pnb.size()==0) {
    MyStrStream str;
    str<<m_name<<"_";
    for (size_t i=0;i<m_pnb.size();i++) str<<m_pnb[i];
    str>>m_name;
  }

  p_histo = new Histogram(type,m_xmin,m_xmax,m_nbins);

}

void PSM_Observable::Evaluate(const Particle_List & pl,double weight, double ncount)
{
  std::vector<Vec4D> moms;
  Vec4D smom(0.,0.,0.,0.);
  for (Particle_List::const_iterator it=pl.begin();it!=pl.end();++it) {
    smom+=(*it)->Momentum();
  }
  moms.push_back(Vec4D(0.5*(smom[0]+smom[3]),0.,0.,0.5*(smom[0]+smom[3])));
  moms.push_back(Vec4D(0.5*(smom[0]-smom[3]),0.,0.,-0.5*(smom[0]-smom[3])));
  for (Particle_List::const_iterator it=pl.begin();it!=pl.end();++it) {
    moms.push_back((*it)->Momentum());
  }
  
  Vec4D ms=Vec4D(0.,0.,0.,0.);
  if (m_pnb.size()>0) {
    for (size_t i=0;i<moms.size();i++){
      int hit=0;
      for(size_t j=0;j<m_pnb.size();j++) {
	if (m_pnb[j]==(int)i) hit = 1;
      }
      if (hit) {
	if (i<2) ms -= moms[i];
	else ms += moms[i];
      }
    } 
    double st=ms.Abs2();
    if (st<0.) st=-sqrt(-st);
    else st=sqrt(st);
    p_histo->Insert(st,weight,ncount);
  }
  else {
    ms = moms[0]+moms[1];
    double y = 0.5 * log( (ms[0]+ms[3])/(ms[0]-ms[3]) );
    p_histo->Insert(y,weight,ncount);
  }
}

void PSM_Observable::Evaluate(const Blob_List & blobs,double value, double ncount)
{
  Particle_List * pl=p_ana->GetParticleList(m_listname);
  Evaluate(*pl,value, ncount);
}

void PSM_Observable::EndEvaluation(double scale) {
    p_histo->MPISync();
    p_histo->Finalize();
}

void PSM_Observable::Restore(double scale) {
    p_histo->Restore();
}

void PSM_Observable::Output(const std::string & pname) {
  p_histo->Output((pname + std::string("/") + m_name+std::string(".dat")).c_str());
}

Primitive_Observable_Base & PSM_Observable::operator+=(const Primitive_Observable_Base & ob)
{
 PSM_Observable * cob = ((PSM_Observable*)(&ob));
 if (p_histo) {
    (*p_histo)+=(cob->p_histo);
  }
  else {
    msg_Out()<<" warning "<<Name()<<" has not overloaded the operator+="<<std::endl;
  }
 
  return *this;
}

void PSM_Observable::Reset()
{
  p_histo->Reset();
}


Primitive_Observable_Base * PSM_Observable::Copy() const 
{
  return new PSM_Observable(m_type,m_xmin,m_xmax,m_nbins,
			    m_pnb[0],m_pnb[1],m_pnb[2],m_pnb[3],
			    m_listname);
}
