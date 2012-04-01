#include "AddOns/Analysis/Observables/Primitive_Observable_Base.H"

#include "ATOOLS/Org/MyStrStream.H"
#include <iomanip>

using namespace ATOOLS;

namespace ANALYSIS {

  class SFour_Particle_Observable_Base: public Primitive_Observable_Base {  
  protected:

    ATOOLS::Flavour m_flav[4];
    size_t          m_item[4];

  public:

    SFour_Particle_Observable_Base
    (const ATOOLS::Flavour flav[4],const size_t item[4],
     const int type,const double min,const double max,const int bins,
     const std::string &inlist,const std::string &name);
    
    void Evaluate(const ATOOLS::Particle_List &particlelist,
		  double weight=1.,double ncount=1);
    
    virtual bool Evaluate(const Particle *p1,const Particle *p2,
			  const Particle *p3,const Particle *p4,
			  double weight=1.,double ncount=1) const = 0;

  };// end of class SFour_Particle_Observable_Base

  class Four_Particle_Plane_Angle: public SFour_Particle_Observable_Base {  
  public:

    Four_Particle_Plane_Angle
    (const ATOOLS::Flavour flav[4],const size_t item[4],
     const int type,const double min,const double max,const int bins,
     const std::string &inlist);
    
    virtual bool Evaluate(const Particle *p1,const Particle *p2,
			  const Particle *p3,const Particle *p4,
			  double weight=1.,double ncount=1) const;

    Primitive_Observable_Base *Copy() const;
    
  };// end of class Four_Particle_Plane_Angle

}// end of namespace ANALYSIS

using namespace ANALYSIS;

template <class Class>
Primitive_Observable_Base *
GetFourParticleSelector(const Argument_Matrix &parameters) 
{									
  if (parameters.size()<1) return NULL;
  if (parameters.size()==1) {
    if (parameters[0].size()<13) return NULL;
    size_t item[4];
    ATOOLS::Flavour flav[4];
    for (size_t i(0);i<4;++i) {
      int kf=ATOOLS::ToType<int>(parameters[0][2*i]);
      flav[i]=ATOOLS::Flavour((kf_code)abs(kf));
      if (kf<0) flav[i]=flav[i].Bar();
      item[i]=ATOOLS::ToType<size_t>(parameters[0][2*i+1]);
    }
    return new Class(flav,item,
		     HistogramType(parameters[0][11]),
		     ATOOLS::ToType<double>(parameters[0][8]),
		     ATOOLS::ToType<double>(parameters[0][9]),
		     ATOOLS::ToType<int>(parameters[0][10]),
		     parameters[0][12]);
  }
  return NULL;
}									

#define DEFINE_FOUR_OBSERVABLE_GETTER_METHOD(CLASS,NAME)		\
  Primitive_Observable_Base *					\
  NAME::operator()(const Argument_Matrix &parameters) const		\
  { return GetFourParticleSelector<CLASS>(parameters); }

#define DEFINE_FOUR_OBSERVABLE_PRINT_METHOD(NAME)		\
  void NAME::PrintInfo(std::ostream &str,const size_t width) const	\
  { str<<"flav1 item1 ... flav4 item4 min max bins Lin|LinErr|Log|LogErr list"; }

#define DEFINE_FOUR_OBSERVABLE_GETTER(CLASS,NAME,TAG)		\
  DECLARE_GETTER(NAME,TAG,Primitive_Observable_Base,Argument_Matrix);	\
  DEFINE_FOUR_OBSERVABLE_GETTER_METHOD(CLASS,NAME)		\
  DEFINE_FOUR_OBSERVABLE_PRINT_METHOD(NAME)

#include "AddOns/Analysis/Main/Primitive_Analysis.H"

SFour_Particle_Observable_Base::SFour_Particle_Observable_Base    
(const ATOOLS::Flavour flav[4],const size_t item[4],
 const int type,const double min,const double max,const int bins,
 const std::string &inlist,const std::string &name):
  Primitive_Observable_Base(type,min,max,bins)
{
  m_name=name;
  for (size_t i(0);i<4;++i) {
    m_flav[i]=flav[i];
    m_item[i]=item[i];
    m_name+="_"+ToString(m_flav[i])+"-"+ToString(m_item[i]);
  }
  m_name+=".dat";
  m_xmin=min;
  m_xmax=max;
  m_listname=inlist;
}

void SFour_Particle_Observable_Base::Evaluate
(const ATOOLS::Particle_List &inlist,double weight,double ncount)
{
  int no(-1);
  size_t pos[4]={std::string::npos,std::string::npos,
		 std::string::npos,std::string::npos};
  for (size_t k(0);k<4;++k) {
    no=-1;
    for (size_t i(0);i<inlist.size();++i) {
      if (inlist[i]->Flav()==m_flav[k] || 
	  m_flav[k].Kfcode()==kf_none) {
	++no;
	if (no==(int)m_item[k]) {
	  pos[k]=i;
	  break;
	}
      }
    }
  }
  for (size_t k(0);k<4;++k) if (pos[k]==std::string::npos) return;
  Evaluate(inlist[pos[0]],inlist[pos[1]],
	   inlist[pos[2]],inlist[pos[3]],weight,ncount);
}

DEFINE_FOUR_OBSERVABLE_GETTER(Four_Particle_Plane_Angle,
			      Four_Particle_Plane_Angle_Getter,"FourPlaneAngle")
  
Four_Particle_Plane_Angle::Four_Particle_Plane_Angle
(const ATOOLS::Flavour flav[4],const size_t item[4],
 const int type,const double min,const double max,const int bins,
 const std::string &inlist):
  SFour_Particle_Observable_Base(flav,item,type,min,max,bins,
				 inlist,"FourPlaneAngle") {}

bool Four_Particle_Plane_Angle::Evaluate(const Particle *p1,const Particle *p2,
					 const Particle *p3,const Particle *p4,
					 double weight,double ncount) const
{
  Vec4D n1(1.0,cross((Vec3D)p1->Momentum(),(Vec3D)p2->Momentum()));
  Vec4D n2(1.0,cross((Vec3D)p3->Momentum(),(Vec3D)p4->Momentum()));
  p_histo->Insert(n1.Theta(n2),weight,ncount);
  return true;
}

Primitive_Observable_Base *Four_Particle_Plane_Angle::Copy() const
{
  return new Four_Particle_Plane_Angle(m_flav,m_item,m_type,m_xmin,m_xmax,
				       m_nbins,m_listname);
}

