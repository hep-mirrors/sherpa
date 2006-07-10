#include "Three_Particle_Observables.H"
#include "Primitive_Analysis.H"
#include "MyStrStream.H"

using namespace ANALYSIS;
using namespace ATOOLS;
using namespace std;

#define DEFINE_GETTER_METHOD(CLASS,NAME)                                \
  Primitive_Observable_Base *                                           \
  NAME::operator()(const String_Matrix &parameters) const               \
  { return new CLASS(parameters); }

#define DEFINE_PRINT_METHOD(NAME)                                       \
  void NAME::PrintInfo(std::ostream &str,const size_t width) const      \
  { str<<"kf1 kf2 kf3 min max bins Lin|LinErr|Log|LogErr [list]"; }

#define DEFINE_CONSTRUCTOR_METHODS(CLASS,TAG)                           \
  CLASS::CLASS(const String_Matrix & parameters) :                      \
      Three_Particle_Observable_Base(parameters, TAG) { }               \
                                                                        \
  CLASS::CLASS(const CLASS * old) :                                     \
      Three_Particle_Observable_Base(*old) { }                          \
                                                                        \
  Primitive_Observable_Base * CLASS::Copy() const { return new CLASS(this); }

#define DEFINE_OBSERVABLE_GETTER(CLASS,NAME,TAG)                        \
  DECLARE_GETTER(NAME,TAG,Primitive_Observable_Base,String_Matrix);     \
  DEFINE_GETTER_METHOD(CLASS,NAME)                                      \
  DEFINE_PRINT_METHOD(NAME)                                             \
  DEFINE_CONSTRUCTOR_METHODS(CLASS,TAG)

#ifdef USING__ROOT
#include "Scaling.H"
#include "My_Root.H"
#include "TH2D.h"
#endif 

Three_Particle_Observable_Base::Three_Particle_Observable_Base(const String_Matrix & parameters,
                                                           const std::string & obsname):
    Primitive_Observable_Base(parameters)
{
  if (parameters.size()==1) {
    if (parameters[0].size()<7) {
      msg.Error()<<"Error in "<<METHOD<<": Not enough parameters for "
          <<"observable "<<obsname<<" in Analysis.dat";
      abort();
    }
    ATOOLS::Flavour f[3];
    for (short unsigned int i=0;i<3;++i) {
      int kf=ATOOLS::ToType<int>(parameters[0][i]);
      f[i]=ATOOLS::Flavour((ATOOLS::kf::code)abs(kf));
      if (kf<0) f[i]=f[i].Bar();
    }
    m_flav1 = f[0];
    m_flav2 = f[1];
    m_flav3 = f[2];

    m_xmin  = ATOOLS::ToType<double>(parameters[0][3]);
    m_xmax  = ATOOLS::ToType<double>(parameters[0][4]);
    int nbins = ATOOLS::ToType<int>(parameters[0][5]);
    m_type  = HistogramType(parameters[0][6]);
    p_histo = new Histogram(m_type,m_xmin,m_xmax,nbins);

    MyStrStream str;
    str<<obsname<<m_flav1<<m_flav2<<m_flav3<<".dat";
    str>>m_name;

    m_listname = parameters[0].size()>7?parameters[0][7]:"Analysed";
  }
  else {
    for (size_t i=0;i<parameters.size();++i) {
      if (parameters[i].size()<2) continue;
      for (short unsigned int j=0;j<3;++j) {
        if (parameters[i][0]==std::string("FLAV")+ATOOLS::ToString(j+1)) {
          int kf=ATOOLS::ToType<int>(parameters[i][1]);
          if(j==0) m_flav1=ATOOLS::Flavour((ATOOLS::kf::code)abs(kf),kf<0);
          if(j==1) m_flav2=ATOOLS::Flavour((ATOOLS::kf::code)abs(kf),kf<0);
          if(j==2) m_flav3=ATOOLS::Flavour((ATOOLS::kf::code)abs(kf),kf<0);
        }
      }
    }

    if (m_name=="SherpaDefault") {
      MyStrStream str;
      str<<obsname<<m_flav1<<m_flav2<<m_flav3;
      str>>m_name;
    }
    if (p_histo->Title()=="SherpaDefault") {
      std::string title = "";
      MyStrStream str;
      str<<obsname<<" of "<<m_flav1<<", "<<m_flav2<<", "<<m_flav3;
      str>>title;
      p_histo->SetTitle(title);
    }
    if (m_listname=="") m_listname="Analysed";
  }
}

Three_Particle_Observable_Base::Three_Particle_Observable_Base(const Three_Particle_Observable_Base * old) :
    Primitive_Observable_Base(*old)
{
  m_flav1=old->m_flav1; m_flav2=old->m_flav2; m_flav3=old->m_flav3;
}

void Three_Particle_Observable_Base::Evaluate(const Particle_List & plist,double weight, int ncount)
{
  for (Particle_List::const_iterator plit1=plist.begin();plit1!=plist.end();++plit1) {
    if ((*plit1)->Flav()==m_flav1) {
      for (Particle_List::const_iterator plit2=plist.begin();plit2!=plist.end();++plit2) {
	if ((*plit2)->Flav()==m_flav2 && plit1!=plit2) {
	  for (Particle_List::const_iterator plit3=plist.begin();plit3!=plist.end();++plit3) {
	    if ((*plit3)->Flav()==m_flav3 && plit1!=plit3 && plit2!=plit3) {
	      Evaluate((*plit1)->Momentum(),(*plit2)->Momentum(),(*plit3)->Momentum(),weight,ncount);
	      return;
	    }
	  }
	}
      }
    }
  }
  p_histo->Insert(0.0,0.0,ncount);
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DEFINE_OBSERVABLE_GETTER(Three_Particle_Y,Three_Particle_Y_Getter,"Y3")

void Three_Particle_Y::Evaluate(const Vec4D & mom1,const Vec4D & mom2,const Vec4D & mom3,double weight, int ncount) 
{    
  double y = (mom1+mom2+mom3).Y();
  p_histo->Insert(y,weight,ncount); 
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DEFINE_OBSERVABLE_GETTER(Three_Particle_DEta,Three_Particle_DEta_Getter,"DEta3")

void Three_Particle_DEta::Evaluate(const Vec4D & mom1,const Vec4D & mom2,const Vec4D & mom3,double weight, int ncount) 
{    
   Vec4D mother = mom1+mom2;
  double deta = abs((mother.Eta()-mom3.Eta()));
  p_histo->Insert(deta,weight,ncount); 
} 

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DEFINE_OBSERVABLE_GETTER(Three_Particle_DPhi,Three_Particle_DPhi_Getter,"DPhi3")

void Three_Particle_DPhi::Evaluate(const Vec4D & mom1,const Vec4D & mom2,const Vec4D & mom3,double weight, int ncount) 
{ 
  
  Vec4D mother = mom1+mom2;

  double pt1=sqrt(mother[1]*mother[1]+mother[2]*mother[2]);
  double pt2=sqrt(mom3[1]*mom3[1]+mom3[2]*mom3[2]);
  double dphi=acos((mother[1]*mom3[1]+mother[2]*mom3[2])/(pt1*pt2));
  p_histo->Insert(dphi,weight,ncount); 
} 

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DEFINE_OBSERVABLE_GETTER(Three_Particle_DR,Three_Particle_DR_Getter,"DR3")

void Three_Particle_DR::Evaluate(const Vec4D & mom1,const Vec4D & mom2,const Vec4D & mom3,double weight, int ncount)
{ 
  
  Vec4D mother = mom1+mom2;
  
  double pt1=sqrt(mother[1]*mother[1]+mother[2]*mother[2]);
  double pt2=sqrt(mom3[1]*mom3[1]+mom3[2]*mom3[2]);
  double dphi=acos((mother[1]*mom3[1]+mother[2]*mom3[2])/(pt1*pt2));
  double c1=mother[3]/Vec3D(mother).Abs();
  double c2=mom3[3]/Vec3D(mom3).Abs();
  double deta=0.5 *log( (1 + c1)*(1 - c2)/((1-c1)*(1+c2)));
  double dr= sqrt(sqr(deta) + sqr(dphi)); 
  p_histo->Insert(dr,weight,ncount); 
} 

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DEFINE_OBSERVABLE_GETTER(Three_Particle_3Mass2,Three_Particle_3Mass2_Getter,"3Mass2")

void Three_Particle_3Mass2::Evaluate(const Vec4D & mom1,const Vec4D & mom2,const Vec4D & mom3,double weight, int ncount)
{ 
  
  double mass = (mom1+mom2+mom3).Abs2();
  p_histo->Insert(mass,weight,ncount); 
} 

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DEFINE_OBSERVABLE_GETTER(Three_Particle_3Mass,Three_Particle_3Mass_Getter,"3Mass")

void Three_Particle_3Mass::Evaluate(const Vec4D & mom1,const Vec4D & mom2,const Vec4D & mom3,double weight, int ncount) 
{ 
  
  double mass = sqrt( (mom1+mom2+mom3).Abs2() );
  p_histo->Insert(mass,weight,ncount); 
} 
