#include "Two_Particle_Observables.H"
#include "Primitive_Analysis.H"
#include "MyStrStream.H"

using namespace ANALYSIS;
using namespace ATOOLS;
using namespace std;

#define DEFINE_GETTER_METHOD(CLASS,NAME)				\
  Primitive_Observable_Base *					\
  NAME::operator()(const String_Matrix &parameters) const		\
  { return new CLASS(parameters); }

#define DEFINE_PRINT_METHOD(NAME)					\
  void NAME::PrintInfo(std::ostream &str,const size_t width) const	\
  { str<<"kf1 kf2 min max bins Lin|LinErr|Log|LogErr [list]"; }

#define DEFINE_CONSTRUCTOR_METHODS(CLASS,TAG)                           \
  CLASS::CLASS(const String_Matrix & parameters) :                      \
      Two_Particle_Observable_Base(parameters, TAG) { }               \
                                                                        \
  CLASS::CLASS(const CLASS * old) :                                     \
      Two_Particle_Observable_Base(*old) { }                            \
                                                                        \
  Primitive_Observable_Base * CLASS::Copy() const { return new CLASS(this); }

#define DEFINE_OBSERVABLE_GETTER(CLASS,NAME,TAG)			\
  DECLARE_GETTER(NAME,TAG,Primitive_Observable_Base,String_Matrix);	\
  DEFINE_GETTER_METHOD(CLASS,NAME)					\
  DEFINE_PRINT_METHOD(NAME)                                             \
  DEFINE_CONSTRUCTOR_METHODS(CLASS,TAG)

// #ifdef USING__ROOT
// #include "Scaling.H"
// #include "My_Root.H"
// #include "TH2D.h"
// #endif 

Two_Particle_Observable_Base::Two_Particle_Observable_Base(const String_Matrix & parameters,
                                                           const std::string & obsname):
    Primitive_Observable_Base(parameters)
{
  if (parameters.size()==1) {
    if (parameters[0].size()<6) {
      msg.Error()<<"Error in "<<METHOD<<": Not enough parameters for "
                 <<"observable "<<obsname<<" in Analysis.dat";
      abort();
    }
    Flavour f[2];
    for (short unsigned int i=0;i<2;++i) {
      int kf=ToType<int>(parameters[0][i]);
      f[i]=Flavour((kf::code)abs(kf));
      if (kf<0) f[i]=f[i].Bar();
    }
    m_flav1 = f[0];
    m_flav2 = f[1];

    m_xmin  = ToType<double>(parameters[0][2]);
    m_xmax  = ToType<double>(parameters[0][3]);
    int nbins = ToType<int>(parameters[0][4]);
    m_type  = HistogramType(parameters[0][5]);
    p_histo = new Histogram(m_type,m_xmin,m_xmax,nbins);

    m_listname = parameters[0].size()>6?parameters[0][6]:"Analysed";

    MyStrStream str;
    str<<obsname<<m_flav1<<m_flav2<<".dat";
    str>>m_name;
  }
  else {
    for (size_t i=0;i<parameters.size();++i) {
      if (parameters[i].size()<2) continue;
      for (short unsigned int j=0;j<2;++j) {
        if (parameters[i][0]==std::string("FLAV")+ToString(j+1)) {
          int kf=ToType<int>(parameters[i][1]);
          if(j==0) m_flav1=Flavour((kf::code)abs(kf),kf<0);
          if(j==1) m_flav2=Flavour((kf::code)abs(kf),kf<0);
        }
      }
    }

    if (m_name=="SherpaDefault") {
      MyStrStream str;
      str<<obsname<<m_flav1<<m_flav2;
      str>>m_name;
    }
    if (p_histo->Title()=="SherpaDefault") {
      std::string title = "";
      MyStrStream str;
      str<<obsname<<" of "<<m_flav1<<", "<<m_flav2;
      str>>title;
      p_histo->SetTitle(title);
    }
    if (m_listname=="") m_listname="Analysed";
  }
}

Two_Particle_Observable_Base::Two_Particle_Observable_Base(const Two_Particle_Observable_Base * old) :
    Primitive_Observable_Base(*old)
{
  m_flav1=old->m_flav1; m_flav2=old->m_flav2;
}

void Two_Particle_Observable_Base::Evaluate(const Particle_List & plist,double weight, int ncount)
{
  for (Particle_List::const_iterator plit1=plist.begin();plit1!=plist.end();++plit1) {
    if ((*plit1)->Flav()==m_flav1) {
      for (Particle_List::const_iterator plit2=plist.begin();plit2!=plist.end();++plit2) {
	if ((*plit2)->Flav()==m_flav2 && plit1!=plit2) {
	  Evaluate((*plit1)->Momentum(),(*plit2)->Momentum(),weight,ncount);
	  return;
	}
      }
    }
  }
  Evaluate(Vec4D(1.,0,0,1.),Vec4D(1.,0,0,-1.),0, ncount);
}



//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DEFINE_OBSERVABLE_GETTER(Two_Particle_Mass,Two_Particle_Mass_Getter,"Mass")

void Two_Particle_Mass::Evaluate(const Vec4D & mom1,const Vec4D & mom2,double weight, int ncount)
{
  double mass = sqrt((mom1+mom2).Abs2());
  p_histo->Insert(mass,weight,ncount);
  if (weight!=0) {
    p_ana->AddData(m_name,new Blob_Data<double>(mass));
  }
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DEFINE_OBSERVABLE_GETTER(Two_Particle_PT,Two_Particle_PT_Getter,"PT2")

void Two_Particle_PT::Evaluate(const Vec4D & mom1,const Vec4D & mom2,double weight, int ncount)
{
  double pt = sqrt(sqr(mom1[1]+mom2[1]) + sqr(mom1[2]+mom2[2]));
  p_histo->Insert(pt,weight,ncount); 
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DEFINE_OBSERVABLE_GETTER(Two_Particle_Scalar_PT,
			 Two_Particle_Scalar_PT_Getter,"SPT2")

void Two_Particle_Scalar_PT::Evaluate(const Vec4D & mom1,const Vec4D & mom2,double weight, int ncount)
{
  p_histo->Insert((mom1.PPerp()+mom2.PPerp())/2.,weight,ncount); 
} 

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

// #ifdef USING__ROOT
// DEFINE_OBSERVABLE_GETTER(Two_Particle_Angles,
// 			 Two_Particle_Angles_Getter,"Theta_2D")
// 
// Two_Particle_Angles::Two_Particle_Angles(const Flavour & _flav1,const Flavour & _flav2,
// 					 int _type,double _xmin,double _xmax,int _nbins,
// 					 const std::string & _name, const std::string & _lname) :
//   Two_Particle_Observable_Base(_flav1,_flav2,_type,_xmin,_xmax,_nbins,"Angles",_lname) 
// { 
//   (*MYROOT::myroot)(new TH2D(ToString(this).c_str(),
// 			     (m_flav1.Name()+std::string("_")+m_flav2.Name()+
// 			      std::string("_Angles")).c_str(),
// 			     64,0.,M_PI,64,0.,M_PI),
// 		    m_flav1.Name()+std::string("_")+m_flav2.Name()+
// 		    std::string("_Angles"));
// }
// 
// void Two_Particle_Angles::Evaluate(const Vec4D & mom1,const Vec4D & mom2,double weight, int ncount) 
// {
//   ((TH2D*)(*MYROOT::myroot)[m_flav1.Name()+std::string("_")+m_flav2.Name()+
// 			    std::string("_Angles")])->Fill(mom1.Theta(),mom2.Theta(),weight);
// } 
// 
// Primitive_Observable_Base * Two_Particle_Angles::Copy() const 
// {
//   return new Two_Particle_Angles(m_flav1,m_flav2,m_type,m_xmin,m_xmax,Nbins(),m_name,m_listname);
// }
// #endif

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DEFINE_OBSERVABLE_GETTER(Two_Particle_Eta,Two_Particle_Eta_Getter,"Eta2")

void Two_Particle_Eta::Evaluate(const Vec4D & mom1,const Vec4D & mom2,double weight, int ncount) 
{
  double eta = (mom1+mom2).Eta();
  p_histo->Insert(eta,weight,ncount); 
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DEFINE_OBSERVABLE_GETTER(Two_Particle_DEta,Two_Particle_DEta_Getter,"DEta")

void Two_Particle_DEta::Evaluate(const Vec4D & mom1,const Vec4D & mom2,double weight, int ncount) 
{    
    double deta = abs(mom1.Eta()-mom2.Eta());
    p_histo->Insert(deta,weight,ncount); 
} 

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DEFINE_OBSERVABLE_GETTER(Two_Particle_Y,Two_Particle_Y_Getter,"Y")

void Two_Particle_Y::Evaluate(const Vec4D & mom1,const Vec4D & mom2,double weight, int ncount) 
{
  double y = (mom1+mom2).Y();
  p_histo->Insert(y,weight,ncount); 
} 

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DEFINE_OBSERVABLE_GETTER(Two_Particle_DY,Two_Particle_DY_Getter,"DY")

void Two_Particle_DY::Evaluate(const Vec4D & mom1,const Vec4D & mom2,double weight, int ncount) 
{    
  double deta = abs((mom1.Y()-mom2.Y()));
  p_histo->Insert(deta,weight,ncount); 
} 

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DEFINE_OBSERVABLE_GETTER(Two_Particle_Angle,Two_Particle_Angle_Getter,"Angle")

void Two_Particle_Angle::Evaluate(const Vec4D & mom1,const Vec4D & mom2,double weight, int ncount) 
{ 
    double pt1=sqrt(mom1[1]*mom1[1]+mom1[2]*mom1[2]+mom1[3]*mom1[3]);
    double pt2=sqrt(mom2[1]*mom2[1]+mom2[2]*mom2[2]+mom2[3]*mom2[3]);
    double phi=acos((mom1[1]*mom2[1]+mom1[2]*mom2[2]+mom1[3]*mom2[3])/(pt1*pt2));
    p_histo->Insert(phi,weight,ncount); 
} 

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DEFINE_OBSERVABLE_GETTER(Two_Particle_DPhi,Two_Particle_DPhi_Getter,"DPhi")

void Two_Particle_DPhi::Evaluate(const Vec4D & mom1,const Vec4D & mom2,double weight, int ncount) 
{ 
    double pt1=sqrt(mom1[1]*mom1[1]+mom1[2]*mom1[2]);
    double pt2=sqrt(mom2[1]*mom2[1]+mom2[2]*mom2[2]);
    double dphi=acos((mom1[1]*mom2[1]+mom1[2]*mom2[2])/(pt1*pt2));
    p_histo->Insert(dphi,weight,ncount); 
} 

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DEFINE_OBSERVABLE_GETTER(Two_Particle_DR,Two_Particle_DR_Getter,"DR")

void Two_Particle_DR::Evaluate(const Vec4D & mom1,const Vec4D & mom2,double weight, int ncount) 
{ 
    double pt1=sqrt(mom1[1]*mom1[1]+mom1[2]*mom1[2]);
    double pt2=sqrt(mom2[1]*mom2[1]+mom2[2]*mom2[2]);
    double dphi=acos((mom1[1]*mom2[1]+mom1[2]*mom2[2])/(pt1*pt2));
    double c1=mom1[3]/Vec3D(mom1).Abs();
    double c2=mom2[3]/Vec3D(mom2).Abs();
    double deta=0.5 *log( (1 + c1)*(1 - c2)/((1-c1)*(1+c2)));
    double dr= sqrt(sqr(deta) + sqr(dphi)); 
    //cout<<"Deat in DR "<<deta<<" DR is :  "<<dr<<endl;
    //if(dr<0.4) {
    //  std::cout<<"\n>>>>>>>>>>>>>>>>>> DR = "<<dr<<"\n";
    //  std::cout<<m_flav1<<"\n"<<mom1<<"\n";
    //  std::cout<<m_flav2<<"\n"<<mom2<<"\n\n";
    //}
    p_histo->Insert(dr,weight,ncount); 
} 

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DEFINE_OBSERVABLE_GETTER(Two_Particle_CMS_Angle,Two_Particle_CMS_Angle_Getter,"CMSAngle")

void Two_Particle_CMS_Angle::Evaluate(const Vec4D & mom1,const Vec4D & mom2,double weight, int ncount) 
{
  Vec4D sum=mom1+mom2;
  Poincare boost(sum);
  Vec4D p1=boost*mom1;

  Vec3D a(sum), b(p1);
  double costh=a*b/(a.Abs()*b.Abs());

  p_histo->Insert(costh,weight,ncount); 
  if (weight!=0) {
    p_ana->AddData(m_name,new Blob_Data<double>(costh));
  }
} 

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DEFINE_OBSERVABLE_GETTER(Two_Particle_Mass2,Two_Particle_Mass2_Getter,"Mass2")

void Two_Particle_Mass2::Evaluate(const Vec4D & mom1,const Vec4D & mom2,double weight, int ncount) 
{
  double mass = (mom1+mom2).Abs2();
  p_histo->Insert(mass,weight,ncount); 
  if (weight!=0) {
    p_ana->AddData(m_name,new Blob_Data<double>(mass));
  }
}
