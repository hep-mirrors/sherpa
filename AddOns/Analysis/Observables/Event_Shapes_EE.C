#include "AddOns/Analysis/Observables/Event_Shapes_EE.H"

using namespace ANALYSIS;

#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include <algorithm>
#include <iomanip>


DECLARE_GETTER(Event_Shapes_EE_Getter,"EEShapes",
	       Analysis_Object,Argument_Matrix);

void Event_Shapes_EE_Getter::PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"{\n"
     <<std::setw(width+7)<<" "<<"InList  list\n"
     <<std::setw(width+7)<<" "<<"OutList list\n"
     <<std::setw(width+7)<<" "<<"Qual    qualifier\n"
     <<std::setw(width+4)<<" "<<"}";
}

Analysis_Object * 
Event_Shapes_EE_Getter::operator()(const Argument_Matrix &parameters) const
{
  std::string inlist="FinalState", outlist="EEShapes";
  ATOOLS::Particle_Qualifier_Base *qualifier=NULL;
  for (size_t i=0;i<parameters.size();++i) {
    const std::vector<std::string> &cur=parameters[i];
    if (cur[0]=="InList" && cur.size()>1) inlist=cur[1];
    else if (cur[0]=="OutList" && cur.size()>1) outlist=cur[1];
    else if (cur[0]=="Qual" && cur.size()>1) {
      if (ATOOLS::rpa->gen.Beam1().IsLepton() && 
	  ATOOLS::rpa->gen.Beam2().IsLepton()) {
	qualifier = ATOOLS::Particle_Qualifier_Getter::GetObject(cur[1],cur[1]);
      }
    }
  }
  if (!qualifier) qualifier=new ATOOLS::Is_Not_Lepton(); 
  return new Event_Shapes_EE(inlist,outlist,qualifier);
}

#include "AddOns/Analysis/Main/Primitive_Analysis.H"

using namespace ATOOLS;
using namespace std;

static bool bigger(const ATOOLS::Vec3D & lhs,const ATOOLS::Vec3D & rhs) {
  return (lhs.Sqr()>rhs.Sqr()); 
}

Event_Shape_EE_Data::Event_Shape_EE_Data(double thru,double m1,double m2,double oblate,
					 ATOOLS::Vec3D ta,ATOOLS::Vec3D m1a,ATOOLS::Vec3D m2a) 
{
  //thrust(thru), major(m1), minor(m2), oblateness(oblate),
  //thrustaxis(ta), majoraxis(m1a), minoraxis(m2a) {}
  thrust = thru; major = m1; minor = m2; oblateness = oblate;
  thrustaxis = ta; majoraxis = m1a; minoraxis = m2a;
}

namespace ANALYSIS {
  std::ostream& operator<<( std::ostream& ostr, const Event_Shape_EE_Data & evt) {
    ostr<<"Event_Shape_Data : "<<evt.thrust<<"/"<<evt.thrustaxis;
    return ostr;
  }
}

namespace ATOOLS {

template <> Blob_Data<Event_Shape_EE_Data>::~Blob_Data() { }
template class Blob_Data<Event_Shape_EE_Data>;

//std::ostream & operator>>(std::ostream &) const;

}

Event_Shapes_EE::Event_Shapes_EE(const std::string & _inlistname,
				 const std::string & _outlistname,
				 SP(Particle_Qualifier_Base) _quali) :
  Final_Selector(_inlistname,_outlistname,-1,_quali),
  m_startaxes(4), m_maxidentaxes(2), m_accuracy(1.e-4),
  m_key(std::string("EvtShapeData"))
{ 
  m_isobs=false;
  m_name        = std::string("Event_Shapes_EE");
}


void Event_Shapes_EE::Evaluate(const Blob_List & blobs,double value,double ncount)
{
  Particle_List * pl=p_ana->GetParticleList(m_inlistname);
  Select(*pl,value,ncount);
}

void Event_Shapes_EE::Evaluate(const Particle_List & pl,double value,double ncount)
{
  Select(pl,value,ncount);
}

void Event_Shapes_EE::Select(const Particle_List & pl_in,double value,double ncount)
{
  Vec3D mom;
  Particle_List * pl_out = new Particle_List;
  for (Particle_List::const_iterator pit=pl_in.begin();pit!=pl_in.end();++pit) {
    if ((*p_qualifier)(*pit)) {
      pl_out->push_back(new Particle(**pit));
      mom = Vec3D((*pit)->Momentum());
      m_vectors.push_back(mom); m_vectors_save.push_back(mom);
    }
  }
  CalculateLinears();

  p_ana->AddData(m_key,new Blob_Data<Event_Shape_EE_Data>(Event_Shape_EE_Data(m_thrust,m_major,m_minor,
									      m_oblateness,
									      m_thrustaxis,m_majoraxis,
									      m_minoraxis)));
  p_ana->AddData(m_key+"_ThrustAxis",new Blob_Data<Vec4D>(Vec4D(0.0,m_thrustaxis)));
  p_ana->AddData(m_key+"_MajorAxis",new Blob_Data<Vec4D>(Vec4D(0.0,m_majoraxis)));
  p_ana->AddData(m_key+"_MinorAxis",new Blob_Data<Vec4D>(Vec4D(0.0,m_minoraxis)));
  m_vectors.clear(); m_vectors_save.clear();
  p_ana->AddParticleList(m_outlistname,pl_out);
}

void Event_Shapes_EE::CalculateLinears()
{
  m_thrust = m_major = m_minor = 0.;
  Vec3D maxthrustaxis, lastaxis, curraxis;
  double maxthrust=0., lastthrust , currthrust;
  unsigned int min_generators = m_startaxes < m_vectors.size() ? m_startaxes : m_vectors.size();

  vector<Vec3D> initialaxes;
  int addsign;

  for (int pass=0; pass<2; pass++) {
    initialaxes.clear();
    if (pass==1) RotateMoms(m_vectors,m_thrustaxis);
    sort(m_vectors.begin(),m_vectors.end(),&bigger);
    
    for(unsigned int i=1;i<=ipow(2,min_generators-1);++i) {
      Vec3D axis;
      for(unsigned int j=1;j<=min_generators;++j) {
	addsign = -1;
	if (ipow(2,j)*((i+ipow(2,j-1)-1)/ipow(2,j)) >= i) addsign = 1;
	axis = axis+addsign*m_vectors[j-1];
      }
      initialaxes.push_back(axis);
    }
    // sort the initial axes with respect to their size ( which corresponds 
    // to their thrust because of the common denominator) 
    sort(initialaxes.begin(),initialaxes.end(), &bigger);
    for(unsigned int j=0;j<initialaxes.size();j++) 
      initialaxes[j] = initialaxes[j]/initialaxes[j].Abs();

    unsigned int ident = 0;
    double sump        = SumP(m_vectors);
    maxthrust          = 0.;
    for(unsigned int j=0; (j<initialaxes.size()) && (ident<m_maxidentaxes); j++) {
      curraxis         = initialaxes[j];
      currthrust       = SumNP(m_vectors,curraxis)/sump;
      lastthrust       = 0.;
      while (currthrust > lastthrust+m_accuracy) {
	lastthrust     = currthrust;
	lastaxis       = curraxis;
	curraxis       = NewAxis(m_vectors,curraxis);
	currthrust     = SumNP(m_vectors,curraxis)/sump;
      }
      // if it gets worse then skip this axis alltogether
      if (lastthrust < maxthrust-m_accuracy) break;
      // if it is a better solution then keep this one
      if (lastthrust > maxthrust+m_accuracy) {
	ident          = 0;
	maxthrustaxis  = lastaxis;
	maxthrust      = lastthrust;
      }
      ident++;
    }
    if (pass==0) { 
      m_thrustaxis = maxthrustaxis; 
      m_thrust     = maxthrust; 
    }
    else { 
      m_majoraxis  = maxthrustaxis; 
      m_major      = CalculateThrust(m_vectors_save,m_majoraxis); 
      m_minoraxis  = cross(m_thrustaxis,m_majoraxis);
      m_minor      = CalculateThrust(m_vectors_save,m_minoraxis);
      m_oblateness = m_major-m_minor;     
    }
  }
}

void Event_Shapes_EE::RotateMoms(vector<Vec3D> & p,const Vec3D & ref) {
  for(vector<Vec3D>::iterator i=p.begin(); i!=p.end(); ++i) (*i) = (*i)-ref*(ref*(*i));
}

Vec3D Event_Shapes_EE::NewAxis(const vector<Vec3D> & p,const Vec3D & ref) {
  Vec3D nextref = Vec3D(0.,0.,0.);
  int addsign;
  for (unsigned int i=0;i<p.size();++i) {
    addsign = 1;
    if (ref*p[i]<0.) addsign = -1;
    nextref = nextref+addsign*p[i];
  }
  return nextref/nextref.Abs();  
}

double Event_Shapes_EE::CalculateThrust(const vector<Vec3D> & p, const Vec3D & n) { 
  double sum_np = 0, sum_p  = 0;
  for (unsigned int i=0; i<p.size(); i++) {
    sum_np += dabs(p[i]*n); sum_p += p[i].Abs();
  }
  return sum_np/sum_p;
}

double Event_Shapes_EE::SumP(const vector<Vec3D> & p) { 
  double sum_p = 0.;
  for (unsigned int i=0; i<p.size(); i++) sum_p += p[i].Abs();
  return sum_p;
}

double Event_Shapes_EE::SumNP(const vector<Vec3D> & p,const Vec3D & n) { 
  double sum_np = 0.;
  for (unsigned int i=0; i<p.size(); i++) sum_np += dabs(p[i]*n);
  return sum_np;
}

unsigned int Event_Shapes_EE::ipow(int base,int exponent) { 
  int result=1;
  if (exponent>0) 
    do {
      result*=base;
    } while(--exponent);
  return result;
}

Analysis_Object * Event_Shapes_EE::GetCopy() const {
  return new Event_Shapes_EE(m_inlistname,m_outlistname,p_qualifier);
}
