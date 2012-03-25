#include "AddOns/Analysis/Observables/Two_Particle_Observables.H"
#include "AddOns/Analysis/Main/Primitive_Analysis.H"
#include "ATOOLS/Org/Message.H"

using namespace ANALYSIS;

#include "ATOOLS/Org/MyStrStream.H"

template <class Class>
Primitive_Observable_Base *const GetObservable(const Argument_Matrix &parameters)
{									
  if (parameters.size()<1) return NULL;
  if (parameters.size()==1) {
    if (parameters[0].size()<6) return NULL;
    ATOOLS::Flavour f[2];
    for (short unsigned int i=0;i<2;++i) {
      int kf=ATOOLS::ToType<int>(parameters[0][i]);
      f[i]=ATOOLS::Flavour((kf_code)abs(kf));
      if (kf<0) f[i]=f[i].Bar();
    }
    std::string list=parameters[0].size()>6?parameters[0][6]:finalstate_list;
    return new Class(f[0],f[1],HistogramType(parameters[0][5]),
		     ATOOLS::ToType<double>(parameters[0][2]),
		     ATOOLS::ToType<double>(parameters[0][3]),
		     ATOOLS::ToType<int>(parameters[0][4]),list);
  }
  else if (parameters.size()<6) return NULL;
  double min=0.0, max=1.0;
  size_t bins=100;
  ATOOLS::Flavour f[2];
  std::string list=finalstate_list, scale="Lin";
  for (size_t i=0;i<parameters.size();++i) {
    if (parameters[i].size()<2) continue;
    for (short unsigned int j=0;j<2;++j) {
      if (parameters[i][0]==std::string("FLAV")+ATOOLS::ToString(j+1)) {
	int kf=ATOOLS::ToType<int>(parameters[i][1]);
	f[j]=ATOOLS::Flavour((kf_code)abs(kf));
	if (kf<0) f[j]=f[j].Bar();
      }
    }
    if (parameters[i][0]=="MIN") min=ATOOLS::ToType<double>(parameters[i][1]);
    else if (parameters[i][0]=="MAX") max=ATOOLS::ToType<double>(parameters[i][1]);
    else if (parameters[i][0]=="BINS") bins=ATOOLS::ToType<int>(parameters[i][1]);
    else if (parameters[i][0]=="SCALE") scale=parameters[i][1];
    else if (parameters[i][0]=="LIST") list=parameters[i][1];
  }
  return new Class(f[0],f[1],HistogramType(scale),min,max,bins,list);
}									

#define DEFINE_GETTER_METHOD(CLASS,NAME)				\
  Primitive_Observable_Base *					\
  NAME::operator()(const Argument_Matrix &parameters) const		\
  { return GetObservable<CLASS>(parameters); }

#define DEFINE_PRINT_METHOD(NAME)					\
  void NAME::PrintInfo(std::ostream &str,const size_t width) const	\
  { str<<"kf1 kf2 min max bins Lin|LinErr|Log|LogErr [list]"; }

#define DEFINE_OBSERVABLE_GETTER(CLASS,NAME,TAG)			\
  DECLARE_GETTER(NAME,TAG,Primitive_Observable_Base,Argument_Matrix);	\
  DEFINE_GETTER_METHOD(CLASS,NAME)					\
  DEFINE_PRINT_METHOD(NAME)

#ifdef USING__ROOT
#include "ATOOLS/Math/Scaling.H"
#include "AddOns/Analysis/Tools/My_Root.H"
#include "TH2D.h"
#endif 

using namespace ATOOLS;
using namespace std;

Two_Particle_Observable_Base::Two_Particle_Observable_Base(const Flavour & flav1,const Flavour & flav2,
							   int type,double xmin,double xmax,int nbins,
							   const std::string & listname,const std::string & name) :
  Primitive_Observable_Base(type,xmin,xmax,nbins), 
  m_flav1(flav1), m_flav2(flav2)
{
  m_listname=listname;
  MyStrStream str;
  str<<name<<m_flav1<<m_flav2<<".dat";
  str>>m_name;
  m_blobtype = std::string("");
  m_blobdisc = false;
}

/*
void Two_Particle_Observable_Base::Evaluate(double value,double weight, double ncount) 
{
  p_histo->Insert(value,weight,ncount); 
}
*/
 


void Two_Particle_Observable_Base::Evaluate(const Particle_List & plist,double weight, double ncount)
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

void Two_Particle_Observable_Base::EvaluateNLOcontrib(const Vec4D & mom1,const Vec4D & mom2,double weight, double ncount) 
{
  msg_Error()<<"ERROR virtual function Two_Particle_Observable_Base::EvaluateNLOcontrib called "<<m_name<<std::endl
     <<" not NLO-ready!!"<<m_name<<std::endl;
} 

void Two_Particle_Observable_Base::EvaluateNLOcontrib(double weight, double ncount)
{
  Particle_List * plist=p_ana->GetParticleList(m_listname);
  for (Particle_List::const_iterator plit1=plist->begin();plit1!=plist->end();++plit1) {
    if ((*plit1)->Flav()==m_flav1) {
      for (Particle_List::const_iterator plit2=plist->begin();plit2!=plist->end();++plit2) {
	if ((*plit2)->Flav()==m_flav2 && plit1!=plit2) {
	  EvaluateNLOcontrib((*plit1)->Momentum(),(*plit2)->Momentum(),weight,ncount);
	  return;
	}
      }
    }
  }
  EvaluateNLOcontrib(Vec4D(1.,0,0,1.),Vec4D(1.,0,0,-1.),0, ncount);
}

void Two_Particle_Observable_Base::EvaluateNLOevt()
{
  p_histo->FinishMCB();
}


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DEFINE_OBSERVABLE_GETTER(Two_Particle_Mass,Two_Particle_Mass_Getter,"Mass")

Two_Particle_Mass::Two_Particle_Mass(const Flavour & flav1, const Flavour & flav2,
				     int type, double xmin, double xmax, int nbins,
				     const std::string & listname) :
  Two_Particle_Observable_Base(flav1,flav2,type,xmin,xmax,nbins,listname,"Mass") { }


void Two_Particle_Mass::Evaluate(const Vec4D & mom1,const Vec4D & mom2,double weight, double ncount) 
{
  double mass = sqrt((mom1+mom2).Abs2());
  p_histo->Insert(mass,weight,ncount); 
  if (weight!=0) {
    p_ana->AddData(m_name,new Blob_Data<double>(mass));
  }
} 

void Two_Particle_Mass::EvaluateNLOcontrib(const Vec4D & mom1,const Vec4D & mom2,double weight, double ncount) 
{
  double mass = sqrt((mom1+mom2).Abs2());
  p_histo->InsertMCB(mass,weight,ncount); 
  if (weight!=0) {
    p_ana->AddData(m_name,new Blob_Data<double>(mass));
  }
} 

Primitive_Observable_Base * Two_Particle_Mass::Copy() const
{
  return new Two_Particle_Mass(m_flav1,m_flav2,m_type,m_xmin,m_xmax,m_nbins,m_listname);
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DEFINE_OBSERVABLE_GETTER(Two_Particle_PT,Two_Particle_PT_Getter,"PT2")

Two_Particle_PT::Two_Particle_PT(const Flavour & flav1,const Flavour & flav2,
				 int type,double xmin,double xmax,int nbins,
				 const std::string & listname) :
  Two_Particle_Observable_Base(flav1,flav2,type,xmin,xmax,nbins,listname,"PT") { }


void Two_Particle_PT::Evaluate(const Vec4D & mom1,const Vec4D & mom2,double weight, double ncount) 
{
  double pt = sqrt(sqr(mom1[1]+mom2[1]) + sqr(mom1[2]+mom2[2]));
  p_histo->Insert(pt,weight,ncount); 
} 

void Two_Particle_PT::EvaluateNLOcontrib(const Vec4D & mom1,const Vec4D & mom2,double weight, double ncount) 
{
  double pt = sqrt(sqr(mom1[1]+mom2[1]) + sqr(mom1[2]+mom2[2]));
  p_histo->InsertMCB(pt,weight,ncount); 
} 

Primitive_Observable_Base * Two_Particle_PT::Copy() const 
{
  return new Two_Particle_PT(m_flav1,m_flav2,m_type,m_xmin,m_xmax,m_nbins,m_listname);
}

DEFINE_OBSERVABLE_GETTER(Two_Particle_ETW,Two_Particle_ETW_Getter,"ET2W")

Two_Particle_ETW::Two_Particle_ETW(const Flavour & flav1,const Flavour & flav2,
				   int type,double xmin,double xmax,int nbins,
				   const std::string & listname) :
  Two_Particle_Observable_Base(flav1,flav2,type,xmin,xmax,nbins,listname,"ETW") { }

void Two_Particle_ETW::Evaluate(const Vec4D & mom1,const Vec4D & mom2,double weight, double ncount) 
{
  double etw = sqrt(sqr(mom1[1]+mom2[1]) + sqr(mom1[2]+mom2[2]) + sqr(Flavour(kf_Wplus).Mass()));
  p_histo->Insert(etw,weight,ncount); 
} 

void Two_Particle_ETW::EvaluateNLOcontrib(const Vec4D & mom1,const Vec4D & mom2,double weight, double ncount) 
{
  double etw = sqrt(sqr(mom1[1]+mom2[1]) + sqr(mom1[2]+mom2[2]) + sqr(Flavour(kf_Wplus).Mass()));
  p_histo->InsertMCB(etw,weight,ncount); 
} 

Primitive_Observable_Base * Two_Particle_ETW::Copy() const 
{
  return new Two_Particle_ETW(m_flav1,m_flav2,m_type,m_xmin,m_xmax,m_nbins,m_listname);
}

DEFINE_OBSERVABLE_GETTER(Two_Particle_Scalar_PT,
			 Two_Particle_Scalar_PT_Getter,"SPT2")

Two_Particle_Scalar_PT::Two_Particle_Scalar_PT(const Flavour & flav1,const Flavour & flav2,
				 int type,double xmin,double xmax,int nbins,
				 const std::string & listname) :
  Two_Particle_Observable_Base(flav1,flav2,type,xmin,xmax,nbins,listname,"SPT") { }

void Two_Particle_Scalar_PT::Evaluate(const Vec4D & mom1,const Vec4D & mom2,double weight, double ncount) 
{
  p_histo->Insert((mom1.PPerp()+mom2.PPerp())/2.,weight,ncount); 
} 

void Two_Particle_Scalar_PT::EvaluateNLOcontrib(const Vec4D & mom1,const Vec4D & mom2,double weight, double ncount) 
{
  p_histo->InsertMCB((mom1.PPerp()+mom2.PPerp())/2.,weight,ncount); 
} 

Primitive_Observable_Base * Two_Particle_Scalar_PT::Copy() const 
{
  return new Two_Particle_Scalar_PT(m_flav1,m_flav2,m_type,m_xmin,m_xmax,m_nbins,m_listname);
}

#ifdef USING__ROOT
DEFINE_OBSERVABLE_GETTER(Two_Particle_Angles,
			 Two_Particle_Angles_Getter,"Theta_2D")

Two_Particle_Angles::Two_Particle_Angles(const Flavour & _flav1,const Flavour & _flav2,
					 int _type,double _xmin,double _xmax,int _nbins,
					 const std::string & _name, const std::string & _lname) :
  Two_Particle_Observable_Base(_flav1,_flav2,_type,_xmin,_xmax,_nbins,_lname,"Angles") 
{ 
  (*MYROOT::myroot)(new TH2D(ATOOLS::ToString(this).c_str(),
			     (m_flav1.IDName()+std::string("_")+m_flav2.IDName()
                              +std::string("_Angles")).c_str(),
			     64,0.,M_PI,64,0.,M_PI),
		    m_flav1.IDName()+std::string("_")+m_flav2.IDName()+
		    std::string("_Angles"));
}

void Two_Particle_Angles::Evaluate(const Vec4D & mom1,const Vec4D & mom2,double weight, double ncount)
{
  ((TH2D*)(*MYROOT::myroot)[m_flav1.IDName()+std::string("_")+m_flav2.IDName()+
          std::string("_Angles")])->Fill(mom1.Theta(),mom2.Theta(),weight);
}

void Two_Particle_Angles::EvaluateNLOcontrib(const Vec4D & mom1,const Vec4D & mom2,double weight, double ncount)
{
  msg_Out()<<METHOD<<": not defined yet, evaluating nothing ...\n";
//   ((TH2D*)(*MYROOT::myroot)[m_flav1.IDName()+std::string("_")+m_flav2.IDName()+
//           std::string("_Angles")])->Fill(mom1.Theta(),mom2.Theta(),weight);
}

Primitive_Observable_Base * Two_Particle_Angles::Copy() const
{
  return new Two_Particle_Angles(m_flav1,m_flav2,m_type,m_xmin,m_xmax,m_nbins,m_name,m_listname);
}
#endif

DEFINE_OBSERVABLE_GETTER(Two_Particle_Eta,Two_Particle_Eta_Getter,"Eta2")

Two_Particle_Eta::Two_Particle_Eta(const Flavour & flav1,const Flavour & flav2,
				 int type,double xmin,double xmax,int nbins,
				 const std::string & listname) :
    Two_Particle_Observable_Base(flav1,flav2,type,xmin,xmax,nbins,listname,"Eta") 
{
}


void Two_Particle_Eta::Evaluate(const Vec4D & mom1,const Vec4D & mom2,double weight, double ncount) 
{
  double eta = (mom1+mom2).Eta();
  p_histo->Insert(eta,weight,ncount); 
} 

void Two_Particle_Eta::EvaluateNLOcontrib(const Vec4D & mom1,const Vec4D & mom2,double weight, double ncount) 
{
  double eta = (mom1+mom2).Eta();
  p_histo->InsertMCB(eta,weight,ncount); 
} 

Primitive_Observable_Base * Two_Particle_Eta::Copy() const 
{
  return new Two_Particle_Eta(m_flav1,m_flav2,m_type,m_xmin,m_xmax,m_nbins,m_listname);
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DEFINE_OBSERVABLE_GETTER(Two_Particle_DEta,Two_Particle_DEta_Getter,"DEta")

Two_Particle_DEta::Two_Particle_DEta(const Flavour & flav1,const Flavour & flav2,
				     int type,double xmin,double xmax,int nbins,
				     const std::string & listname) :
    Two_Particle_Observable_Base(flav1,flav2,type,xmin,xmax,nbins,listname,"deta") 

{ 
}


void Two_Particle_DEta::Evaluate(const Vec4D & mom1,const Vec4D & mom2,double weight, double ncount) 
{    
    double deta = abs((mom1.Eta()-mom2.Eta()));
    p_histo->Insert(deta,weight,ncount); 
} 

void Two_Particle_DEta::EvaluateNLOcontrib(const Vec4D & mom1,const Vec4D & mom2,double weight, double ncount) 
{    
    double deta = abs((mom1.Eta()-mom2.Eta()));
    p_histo->InsertMCB(deta,weight,ncount); 
} 

Primitive_Observable_Base * Two_Particle_DEta::Copy() const 
{
    return new Two_Particle_DEta(m_flav1,m_flav2,m_type,m_xmin,m_xmax,m_nbins,m_listname);
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DEFINE_OBSERVABLE_GETTER(Two_Particle_Y,Two_Particle_Y_Getter,"Y2")

Two_Particle_Y::Two_Particle_Y(const Flavour & flav1,const Flavour & flav2,
                               int type,double xmin,double xmax,int nbins,
                               const std::string & listname) :
  Two_Particle_Observable_Base(flav1,flav2,type,xmin,xmax,nbins,listname,"Y")
{
}

void Two_Particle_Y::Evaluate(const Vec4D & mom1,const Vec4D & mom2,
                              double weight, double ncount)
{
  double y = (mom1+mom2).Y();
  p_histo->Insert(y,weight,ncount);
}

void Two_Particle_Y::EvaluateNLOcontrib(const Vec4D & mom1,const Vec4D & mom2,double weight, double ncount) 
{
  double y = (mom1+mom2).Y();
  p_histo->InsertMCB(y,weight,ncount); 
} 

Primitive_Observable_Base * Two_Particle_Y::Copy() const
{
  return new Two_Particle_Y(m_flav1,m_flav2,m_type,m_xmin,m_xmax,m_nbins,m_listname);
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DEFINE_OBSERVABLE_GETTER(Two_Particle_DY,Two_Particle_DY_Getter,"DY")

Two_Particle_DY::Two_Particle_DY(const Flavour & flav1,const Flavour & flav2,
				     int type,double xmin,double xmax,int nbins,
				     const std::string & listname) :
    Two_Particle_Observable_Base(flav1,flav2,type,xmin,xmax,nbins,listname,"DY") 

{ 
}


void Two_Particle_DY::Evaluate(const Vec4D & mom1,const Vec4D & mom2,double weight, double ncount) 
{    
  double deta = abs((mom1.Y()-mom2.Y()));
  p_histo->Insert(deta,weight,ncount); 
} 

void Two_Particle_DY::EvaluateNLOcontrib(const Vec4D & mom1,const Vec4D & mom2,double weight, double ncount) 
{    
  double deta = abs((mom1.Y()-mom2.Y()));
  p_histo->InsertMCB(deta,weight,ncount); 
} 

Primitive_Observable_Base * Two_Particle_DY::Copy() const 
{
    return new Two_Particle_DY(m_flav1,m_flav2,m_type,m_xmin,m_xmax,m_nbins,m_listname);
}
    
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DEFINE_OBSERVABLE_GETTER(Two_Particle_Angle,Two_Particle_Angle_Getter,"Angle")

Two_Particle_Angle::Two_Particle_Angle(const Flavour & flav1,const Flavour & flav2,
				     int type,double xmin,double xmax,int nbins,
				     const std::string & listname) :
    Two_Particle_Observable_Base(flav1,flav2,type,xmin,xmax,nbins,listname,"Angle") 
{ 
}


void Two_Particle_Angle::Evaluate(const Vec4D & mom1,const Vec4D & mom2,double weight, double ncount) 
{ 
    double pt1=sqrt(mom1[1]*mom1[1]+mom1[2]*mom1[2]+mom1[3]*mom1[3]);
    double pt2=sqrt(mom2[1]*mom2[1]+mom2[2]*mom2[2]+mom2[3]*mom2[3]);
    double phi=acos((mom1[1]*mom2[1]+mom1[2]*mom2[2]+mom1[3]*mom2[3])/(pt1*pt2));
    p_histo->Insert(phi,weight,ncount); 
} 

void Two_Particle_Angle::EvaluateNLOcontrib(const Vec4D & mom1,const Vec4D & mom2,double weight, double ncount) 
{ 
    double pt1=sqrt(mom1[1]*mom1[1]+mom1[2]*mom1[2]+mom1[3]*mom1[3]);
    double pt2=sqrt(mom2[1]*mom2[1]+mom2[2]*mom2[2]+mom2[3]*mom2[3]);
    double phi=acos((mom1[1]*mom2[1]+mom1[2]*mom2[2]+mom1[3]*mom2[3])/(pt1*pt2));
    p_histo->InsertMCB(phi,weight,ncount); 
} 

Primitive_Observable_Base * Two_Particle_Angle::Copy() const 
{
    return new Two_Particle_Angle(m_flav1,m_flav2,m_type,m_xmin,m_xmax,m_nbins,m_listname);
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DEFINE_OBSERVABLE_GETTER(Two_Particle_Phi,Two_Particle_Phi_Getter,"Phi2")

Two_Particle_Phi::Two_Particle_Phi(const Flavour & flav1,const Flavour & flav2,
                                   int type,double xmin,double xmax,int nbins,
                                   const std::string & listname) :
    Two_Particle_Observable_Base(flav1,flav2,type,xmin,xmax,nbins,listname,"Phi")
{
}


void Two_Particle_Phi::Evaluate(const Vec4D & mom1,const Vec4D & mom2,double weight, double ncount)
{
    p_histo->Insert((mom1+mom2).Phi(),weight,ncount);
}

void Two_Particle_Phi::EvaluateNLOcontrib(const Vec4D & mom1,const Vec4D & mom2,double weight, double ncount)
{
    p_histo->InsertMCB((mom1+mom2).Phi(),weight,ncount);
}

Primitive_Observable_Base * Two_Particle_Phi::Copy() const
{
    return new Two_Particle_Phi(m_flav1,m_flav2,m_type,m_xmin,m_xmax,m_nbins,m_listname);
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DEFINE_OBSERVABLE_GETTER(Two_Particle_DPhi,Two_Particle_DPhi_Getter,"DPhi")

Two_Particle_DPhi::Two_Particle_DPhi(const Flavour & flav1,const Flavour & flav2,
				     int type,double xmin,double xmax,int nbins,
				     const std::string & listname) :
    Two_Particle_Observable_Base(flav1,flav2,type,xmin,xmax,nbins,listname,"dphi") 
{ 
}


void Two_Particle_DPhi::Evaluate(const Vec4D & mom1,const Vec4D & mom2,double weight, double ncount) 
{ 
    double pt1=sqrt(mom1[1]*mom1[1]+mom1[2]*mom1[2]);
    double pt2=sqrt(mom2[1]*mom2[1]+mom2[2]*mom2[2]);
    double dphi=acos((mom1[1]*mom2[1]+mom1[2]*mom2[2])/(pt1*pt2));
    p_histo->Insert(dphi,weight,ncount); 
} 

void Two_Particle_DPhi::EvaluateNLOcontrib(const Vec4D & mom1,const Vec4D & mom2,double weight, double ncount) 
{ 
    double pt1=sqrt(mom1[1]*mom1[1]+mom1[2]*mom1[2]);
    double pt2=sqrt(mom2[1]*mom2[1]+mom2[2]*mom2[2]);
    double dphi=acos((mom1[1]*mom2[1]+mom1[2]*mom2[2])/(pt1*pt2));
    p_histo->InsertMCB(dphi,weight,ncount); 
} 

Primitive_Observable_Base * Two_Particle_DPhi::Copy() const 
{
    return new Two_Particle_DPhi(m_flav1,m_flav2,m_type,m_xmin,m_xmax,m_nbins,m_listname);
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DEFINE_OBSERVABLE_GETTER(Two_Particle_DR,Two_Particle_DR_Getter,"DR")

Two_Particle_DR::Two_Particle_DR(const Flavour & flav1,const Flavour & flav2,
				 int type, double xmin, double xmax, int nbins,
				 const std::string & listname) :
    Two_Particle_Observable_Base(flav1,flav2,type,xmin,xmax,nbins,listname,"dr") 
{
 }


void Two_Particle_DR::Evaluate(const Vec4D & mom1,const Vec4D & mom2,double weight, double ncount) 
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

void Two_Particle_DR::EvaluateNLOcontrib(const Vec4D & mom1,const Vec4D & mom2,double weight, double ncount) 
{ 
    double pt1=sqrt(mom1[1]*mom1[1]+mom1[2]*mom1[2]);
    double pt2=sqrt(mom2[1]*mom2[1]+mom2[2]*mom2[2]);
    double dphi=acos((mom1[1]*mom2[1]+mom1[2]*mom2[2])/(pt1*pt2));
    double c1=mom1[3]/Vec3D(mom1).Abs();
    double c2=mom2[3]/Vec3D(mom2).Abs();
    double deta=0.5 *log( (1 + c1)*(1 - c2)/((1-c1)*(1+c2)));
    double dr= sqrt(sqr(deta) + sqr(dphi)); 
    p_histo->InsertMCB(dr,weight,ncount); 
} 

Primitive_Observable_Base * Two_Particle_DR::Copy() const 
{
    return new Two_Particle_DR(m_flav1,m_flav2,m_type,m_xmin,m_xmax,m_nbins,m_listname);
}


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DEFINE_OBSERVABLE_GETTER(Two_Particle_CMS_Angle,Two_Particle_CMS_Angle_Getter,"CMSAngle")
// angle of particle in CMS system of particles 1+2 relative to boost direction

Two_Particle_CMS_Angle::Two_Particle_CMS_Angle(const Flavour & flav1, const Flavour & flav2,
				     int type, double xmin, double xmax, int nbins,
				     const std::string & listname) :
  Two_Particle_Observable_Base(flav1,flav2,type,xmin,xmax,nbins,listname,"CMSAngle") { }


void Two_Particle_CMS_Angle::Evaluate(const Vec4D & mom1,const Vec4D & mom2,double weight, double ncount) 
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

void Two_Particle_CMS_Angle::EvaluateNLOcontrib(const Vec4D & mom1,const Vec4D & mom2,double weight, double ncount) 
{
  Vec4D sum=mom1+mom2;
  Poincare boost(sum);
  Vec4D p1=boost*mom1;

  Vec3D a(sum), b(p1);
  double costh=a*b/(a.Abs()*b.Abs());

  p_histo->InsertMCB(costh,weight,ncount); 
  if (weight!=0) {
    p_ana->AddData(m_name,new Blob_Data<double>(costh));
  }
} 

Primitive_Observable_Base * Two_Particle_CMS_Angle::Copy() const
{
  return new Two_Particle_CMS_Angle(m_flav1,m_flav2,m_type,m_xmin,m_xmax,m_nbins,m_listname);
}
 
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DEFINE_OBSERVABLE_GETTER(Two_Particle_Mass2,Two_Particle_Mass2_Getter,"Mass2")

Two_Particle_Mass2::Two_Particle_Mass2(const Flavour & flav1, const Flavour & flav2,
				     int type, double xmin, double xmax, int nbins,
				     const std::string & listname) :
  Two_Particle_Observable_Base(flav1,flav2,type,xmin,xmax,nbins,listname,"Mass2") { }


void Two_Particle_Mass2::Evaluate(const Vec4D & mom1,const Vec4D & mom2,double weight, double ncount) 
{
  double mass = (mom1+mom2).Abs2();
  p_histo->Insert(mass,weight,ncount); 
  if (weight!=0) {
    p_ana->AddData(m_name,new Blob_Data<double>(mass));
  }
} 

void Two_Particle_Mass2::EvaluateNLOcontrib(const Vec4D & mom1,const Vec4D & mom2,double weight, double ncount) 
{
  double mass = (mom1+mom2).Abs2();
  p_histo->InsertMCB(mass,weight,ncount); 
  if (weight!=0) {
    p_ana->AddData(m_name,new Blob_Data<double>(mass));
  }
} 

Primitive_Observable_Base * Two_Particle_Mass2::Copy() const
{
  return new Two_Particle_Mass2(m_flav1,m_flav2,m_type,m_xmin,m_xmax,m_nbins,m_listname);
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DEFINE_OBSERVABLE_GETTER(Two_Particle_MT2,Two_Particle_MT2_Getter,"MT2")

Two_Particle_MT2::Two_Particle_MT2(const Flavour & flav1, const Flavour & flav2,
				   int type, double xmin, double xmax, int nbins,
				   const std::string & listname) :
  Two_Particle_Observable_Base(flav1,flav2,type,xmin,xmax,nbins,listname,"MT2") { }


void Two_Particle_MT2::Evaluate(const Vec4D & mom1,const Vec4D & mom2,double weight, double ncount) 
{
  double mass = sqrt(2.*(mom1.PPerp()*mom2.PPerp()-mom1[1]*mom2[1]-mom1[2]*mom2[2]));
  p_histo->Insert(mass,weight,ncount); 
  if (weight!=0) {
    p_ana->AddData(m_name,new Blob_Data<double>(mass));
  }
} 

void Two_Particle_MT2::EvaluateNLOcontrib(const Vec4D & mom1,const Vec4D & mom2,double weight, double ncount) 
{
  double mass = sqrt(2.*(mom1.PPerp()*mom2.PPerp()-mom1[1]*mom2[1]-mom1[2]*mom2[2]));
  p_histo->InsertMCB(mass,weight,ncount); 
  if (weight!=0) {
    p_ana->AddData(m_name,new Blob_Data<double>(mass));
  }
} 

Primitive_Observable_Base * Two_Particle_MT2::Copy() const
{
  return new Two_Particle_MT2(m_flav1,m_flav2,m_type,m_xmin,m_xmax,m_nbins,m_listname);
}

