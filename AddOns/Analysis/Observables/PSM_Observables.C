#include "AddOns/Analysis/Observables/PSM_Observables.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "AddOns/Analysis/Main/Primitive_Analysis.H"
#include "ATOOLS/Org/Shell_Tools.H"
#include "ATOOLS/Org/Message.H"

using namespace ANALYSIS;

#include "ATOOLS/Org/MyStrStream.H"

template <class Class>
Primitive_Observable_Base *const GetObservable(const Argument_Matrix &parameters)
{									
  if (parameters.size()<1) return NULL;
  if (parameters.size()==1) {
    if (parameters[0].size()<8) return NULL;
    std::string list=parameters[0].size()>8?parameters[0][8]:finalstate_list;
    return new Class(HistogramType(parameters[0][7]),
		     ATOOLS::ToType<double>(parameters[0][0]),
		     ATOOLS::ToType<double>(parameters[0][1]),
		     ATOOLS::ToType<int>(parameters[0][2]),
		     ATOOLS::ToType<int>(parameters[0][3]),
		     ATOOLS::ToType<int>(parameters[0][4]),
		     ATOOLS::ToType<int>(parameters[0][5]),
		     ATOOLS::ToType<int>(parameters[0][6]),list);
  }
  else if (parameters.size()<8) return NULL;
  double min=0.0, max=1.0;
  size_t bins=100;
  int p0=-1,p1=-1,p2=-1,p3=-1;
  std::string list=finalstate_list, scale="Lin";
  for (size_t i=0;i<parameters.size();++i) {
    if (parameters[i].size()<2) continue;
    if (parameters[i][0]=="MIN") min=ATOOLS::ToType<double>(parameters[i][1]);
    else if (parameters[i][0]=="MAX") max=ATOOLS::ToType<double>(parameters[i][1]);
    else if (parameters[i][0]=="BINS") bins=ATOOLS::ToType<int>(parameters[i][1]);
    else if (parameters[i][0]=="PN0")  p0=ATOOLS::ToType<int>(parameters[i][1]);
    else if (parameters[i][0]=="PN1")  p1=ATOOLS::ToType<int>(parameters[i][1]);
    else if (parameters[i][0]=="PN2")  p2=ATOOLS::ToType<int>(parameters[i][1]);
    else if (parameters[i][0]=="PN3")  p3=ATOOLS::ToType<int>(parameters[i][1]);
    else if (parameters[i][0]=="LIST") list=parameters[i][1];
  }
  return new Class(HistogramType(scale),min,max,bins,p0,p1,p2,p3,list);
}									

#define DEFINE_GETTER_METHOD(CLASS,NAME)				\
  Primitive_Observable_Base *					\
  NAME::operator()(const Argument_Matrix &parameters) const		\
  { return GetObservable<CLASS>(parameters); }

#define DEFINE_PRINT_METHOD(NAME)					\
  void NAME::PrintInfo(std::ostream &str,const size_t width) const	\
  { str<<"min max bins pn0 pn1 pn2 pn3 Lin|LinErr|Log|LogErr [list]"; }

#define DEFINE_OBSERVABLE_GETTER(CLASS,NAME,TAG)			\
  DECLARE_GETTER(NAME,TAG,Primitive_Observable_Base,Argument_Matrix);	\
  DEFINE_GETTER_METHOD(CLASS,NAME)					\
  DEFINE_PRINT_METHOD(NAME)

using namespace ATOOLS;

DEFINE_OBSERVABLE_GETTER(PSM_Observable,PSM_Observable_Getter,"PSM")

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
  ATOOLS::MakeDir(pname); 
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
