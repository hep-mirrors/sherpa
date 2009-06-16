#include "AddOns/Analysis/Observables/One_Particle_Observables.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Message.H"
#include "AddOns/Analysis/Main/Primitive_Analysis.H"

using namespace ANALYSIS;

#include "ATOOLS/Org/MyStrStream.H"

template <class Class>
Primitive_Observable_Base *const GetObservable(const Argument_Matrix &parameters)
{									
  if (parameters.size()<1) return NULL;
  if (parameters.size()==1) {
    if (parameters[0].size()<5) return NULL;
    int kf=ATOOLS::ToType<int>(parameters[0][0]);
    ATOOLS::Flavour flavour((kf_code)abs(kf));
    if (kf<0) flavour=flavour.Bar();
    std::string list=parameters[0].size()>5?parameters[0][5]:finalstate_list;
    return new Class(flavour,HistogramType(parameters[0][4]),
		     ATOOLS::ToType<double>(parameters[0][1]),
		     ATOOLS::ToType<double>(parameters[0][2]),
		     ATOOLS::ToType<int>(parameters[0][3]),list);
  }
  else if (parameters.size()<5) return NULL;
  double min=0.0, max=1.0;
  size_t bins=100;
  std::string list=finalstate_list, scale="Lin";
  ATOOLS::Flavour flavour;
  for (size_t i=0;i<parameters.size();++i) {
    if (parameters[i].size()<2) continue;
    if (parameters[i][0]=="FLAV") {
      int kf=ATOOLS::ToType<int>(parameters[i][1]);
      flavour=ATOOLS::Flavour((kf_code)abs(kf));
      if (kf<0) flavour=flavour.Bar();
    }
    else if (parameters[i][0]=="MIN") min=ATOOLS::ToType<double>(parameters[i][1]);
    else if (parameters[i][0]=="MAX") max=ATOOLS::ToType<double>(parameters[i][1]);
    else if (parameters[i][0]=="BINS") bins=ATOOLS::ToType<int>(parameters[i][1]);
    else if (parameters[i][0]=="SCALE") scale=parameters[i][1];
    else if (parameters[i][0]=="LIST") list=parameters[i][1];
  }
  return new Class(flavour,HistogramType(scale),min,max,bins,list);
}									

#define DEFINE_GETTER_METHOD(CLASS,NAME)				\
  Primitive_Observable_Base *					\
  NAME::operator()(const Argument_Matrix &parameters) const		\
  { return GetObservable<CLASS>(parameters); }

#define DEFINE_PRINT_METHOD(NAME)					\
  void NAME::PrintInfo(std::ostream &str,const size_t width) const	\
  { str<<"kf min max bins Lin|LinErr|Log|LogErr [list]"; }

#define DEFINE_OBSERVABLE_GETTER(CLASS,NAME,TAG)			\
  DECLARE_GETTER(NAME,TAG,Primitive_Observable_Base,Argument_Matrix);	\
  DEFINE_GETTER_METHOD(CLASS,NAME)					\
  DEFINE_PRINT_METHOD(NAME)

using namespace ATOOLS;
using namespace std;

One_Particle_Observable_Base::One_Particle_Observable_Base(const Flavour & flav,
							   int type,double xmin,double xmax,int nbins,
							   const std::string & listname, const std::string & name) :
  Primitive_Observable_Base(type,xmin,xmax,nbins), 
  m_flav(flav)
{
  MyStrStream str;
  str<<name<<m_flav<<".dat";
  str>>m_name;

  if (listname!=std::string("")) m_listname = listname;
  m_blobtype = std::string("");
  m_blobdisc = false;
}

void One_Particle_Observable_Base::Evaluate(double value,double weight, double ncount) 
{
  p_histo->Insert(value,weight,ncount); 
}

 
void One_Particle_Observable_Base::Evaluate(int nout,const ATOOLS::Vec4D * moms,const ATOOLS::Flavour * flavs,
					    double weight, double ncount) 
{
  for (int i=0;i<nout;i++) { if (flavs[i]==m_flav) Evaluate(moms[i],weight,ncount); }
}


void One_Particle_Observable_Base::Evaluate(const Particle_List & plist,double weight,double ncount )
{
  for (Particle_List::const_iterator plit=plist.begin();plit!=plist.end();++plit) {
    if ((*plit)->Flav()==m_flav) {
      Evaluate((*plit)->Momentum(),weight, ncount);
      return;
    }
  }
  
  Evaluate(Vec4D(1.,0,0,1.),0, ncount);
}

void One_Particle_Observable_Base::EvaluateNLOcontrib(const ATOOLS::Vec4D & mom,
							      double weight, double ncount)
{
  msg_Error()<<"ERROR virtual function One_Particle_Observable_Base::EvaluateNLOcontrib called "<<m_name<<std::endl
     <<" not NLO-ready!!"<<m_name<<std::endl;
}

void One_Particle_Observable_Base::EvaluateNLOcontrib(double weight,double ncount )
{
  Particle_List * plist=p_ana->GetParticleList(m_listname);
  for (Particle_List::const_iterator plit=plist->begin();plit!=plist->end();++plit) {
    if ((*plit)->Flav()==m_flav) {
      EvaluateNLOcontrib((*plit)->Momentum(),weight, ncount);
      return;
    }
  }
  
  EvaluateNLOcontrib(Vec4D(1.,0,0,1.),0, ncount);
}

void One_Particle_Observable_Base::EvaluateNLOevt()
{
  p_histo->FinishMCB();
}


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DEFINE_OBSERVABLE_GETTER(One_Particle_ET,One_Particle_ET_Getter,"ET")

One_Particle_ET::One_Particle_ET(const Flavour & flav,
				 int type,double xmin,double xmax,int nbins,
				 const std::string & listname) :
  One_Particle_Observable_Base(flav,type,xmin,xmax,nbins,listname,"ET") { }


void One_Particle_ET::Evaluate(const Vec4D & mom,double weight,double ncount) 
{
  double pt2 = sqr(mom[1])+sqr(mom[2]);
  double p2  = sqr(mom[3])+pt2;
  double net = mom[0]*sqrt(pt2/p2);

  p_histo->Insert(net,weight,ncount); 
} 

void One_Particle_ET::EvaluateNLOcontrib(const Vec4D & mom,double weight,double ncount) 
{
  double pt2 = sqr(mom[1])+sqr(mom[2]);
  double p2  = sqr(mom[3])+pt2;
  double net = mom[0]*sqrt(pt2/p2);
  p_histo->InsertMCB(net,weight,ncount); 
} 

Primitive_Observable_Base * One_Particle_ET::Copy() const
{
  return new One_Particle_ET(m_flav,m_type,m_xmin,m_xmax,m_nbins,m_listname);
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DEFINE_OBSERVABLE_GETTER(One_Particle_PT,One_Particle_PT_Getter,"PT")

One_Particle_PT::One_Particle_PT(const Flavour & flav,
				 int type,double xmin,double xmax,int nbins,
				 const std::string & listname) :
  One_Particle_Observable_Base(flav,type,xmin,xmax,nbins,listname,"PT") { }


void One_Particle_PT::Evaluate(const Vec4D & mom,double weight, double ncount) 
{
  double pt = sqrt(sqr(mom[1])+sqr(mom[2]));
  p_histo->Insert(pt,weight,ncount); 
} 

void One_Particle_PT::EvaluateNLOcontrib(const Vec4D & mom,double weight, double ncount) 
{
  double pt = sqrt(sqr(mom[1])+sqr(mom[2]));
  p_histo->InsertMCB(pt,weight,ncount); 
} 

Primitive_Observable_Base * One_Particle_PT::Copy() const 
{
  return new One_Particle_PT(m_flav,m_type,m_xmin,m_xmax,m_nbins,m_listname);
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DEFINE_OBSERVABLE_GETTER(One_Particle_Eta,One_Particle_Eta_Getter,"Eta")

One_Particle_Eta::One_Particle_Eta(const Flavour & flav,
				   int type,double xmin,double xmax,int nbins,
				   const std::string & listname) :
  One_Particle_Observable_Base(flav,type,xmin,xmax,nbins,listname,"Eta") { }


void One_Particle_Eta::Evaluate(const Vec4D & mom,double weight, double ncount) 
{
  double pt2=sqr(mom[1])+sqr(mom[2]);
  double pp =sqrt(pt2+sqr(mom[3]));
  double pz =dabs(mom[3]);
  double sn =mom[3]/pz;
  double value= sn*20.;
  if (pt2>1.e-10*pp*pp) {
    value = sn*0.5*log(sqr(pp+pz)/pt2);
  }
  p_histo->Insert(value,weight,ncount);
} 

void One_Particle_Eta::EvaluateNLOcontrib(const Vec4D & mom,double weight, double ncount) 
{
  double pt2=sqr(mom[1])+sqr(mom[2]);
  double pp =sqrt(pt2+sqr(mom[3]));
  double pz =dabs(mom[3]);
  double sn =mom[3]/pz;
  double value= sn*20.;
  if (pt2>1.e-10*pp*pp) {
    value = sn*0.5*log(sqr(pp+pz)/pt2);
  }
  p_histo->InsertMCB(value,weight,ncount);
} 

Primitive_Observable_Base * One_Particle_Eta::Copy() const
{
  return new One_Particle_Eta(m_flav,m_type,m_xmin,m_xmax,m_nbins,m_listname);
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DEFINE_OBSERVABLE_GETTER(One_Particle_Y,One_Particle_Y_Getter,"Y")

One_Particle_Y::One_Particle_Y(const Flavour & flav,
				   int type,double xmin,double xmax,int nbins,
				   const std::string & listname) :
  One_Particle_Observable_Base(flav,type,xmin,xmax,nbins,listname,"Y") { }


void One_Particle_Y::Evaluate(const Vec4D & mom,double weight, double ncount) 
{
  p_histo->Insert(mom.Y(),weight,ncount);
} 

void One_Particle_Y::EvaluateNLOcontrib(const Vec4D & mom,double weight, double ncount) 
{
  p_histo->InsertMCB(mom.Y(),weight,ncount);
} 

Primitive_Observable_Base * One_Particle_Y::Copy() const
{
  return new One_Particle_Y(m_flav,m_type,m_xmin,m_xmax,m_nbins,m_listname);
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DEFINE_OBSERVABLE_GETTER(One_Particle_E,One_Particle_E_Getter,"E")

One_Particle_E::One_Particle_E(const Flavour & flav,
			       int type,double xmin,double xmax,int nbins,
			       const std::string & listname) :
  One_Particle_Observable_Base(flav,type,xmin,xmax,nbins,listname,"E") { }


void One_Particle_E::Evaluate(const Vec4D & mom,double weight, double ncount) 
{
  double E = mom[0];
  p_histo->Insert(E,weight,ncount); 
} 

void One_Particle_E::EvaluateNLOcontrib(const Vec4D & mom,double weight, double ncount) 
{
  double E = mom[0];
  p_histo->InsertMCB(E,weight,ncount); 
} 

Primitive_Observable_Base * One_Particle_E::Copy() const
{
  return new One_Particle_E(m_flav,m_type,m_xmin,m_xmax,m_nbins,m_listname);
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DEFINE_OBSERVABLE_GETTER(One_Particle_M,One_Particle_M_Getter,"M")

One_Particle_M::One_Particle_M(const Flavour & flav,
				     int type,double xmin,double xmax,int nbins,
				     const std::string & listname) :
  One_Particle_Observable_Base(flav,type,xmin,xmax,nbins,listname,"M") { }


void One_Particle_M::Evaluate(const Vec4D & mom,double weight, double ncount) 
{
  double mass = mom.Mass();
  p_histo->Insert(mass,weight,ncount); 
} 

Primitive_Observable_Base * One_Particle_M::Copy() const
{
  return new One_Particle_M(m_flav,m_type,m_xmin,m_xmax,m_nbins,m_listname);
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DEFINE_OBSERVABLE_GETTER(One_Particle_BeamAngle,One_Particle_BeamAngle_Getter,"BeamAngle")

One_Particle_BeamAngle::One_Particle_BeamAngle(const Flavour & flav,
					       int type,double xmin,double xmax,int nbins,
					       const std::string & listname) :
  One_Particle_Observable_Base(flav,type,xmin,xmax,nbins,listname,"BeamAngle") { }


void One_Particle_BeamAngle::Evaluate(const Vec4D & mom,double weight, double ncount) 
{
  double ct = mom.CosTheta();
  p_histo->Insert(ct,weight,ncount); 
} 

void One_Particle_BeamAngle::EvaluateNLOcontrib(const Vec4D & mom,double weight, double ncount) 
{
  double ct = mom.CosTheta();
  p_histo->InsertMCB(ct,weight,ncount); 
} 

Primitive_Observable_Base * One_Particle_BeamAngle::Copy() const
{
  return new One_Particle_BeamAngle(m_flav,m_type,m_xmin,m_xmax,m_nbins,m_listname);
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DEFINE_OBSERVABLE_GETTER(One_Particle_EVis,One_Particle_EVis_Getter,"EVis")

One_Particle_EVis::One_Particle_EVis(const Flavour & flav,
			       int type,double xmin,double xmax,int nbins,
			       const std::string & listname) :
  One_Particle_Observable_Base(flav,type,xmin,xmax,nbins,listname,"EVis") { }


void One_Particle_EVis::Evaluate(const Vec4D & mom,double weight, double ncount) 
{ } 

void One_Particle_EVis::Evaluate(int nout,const ATOOLS::Vec4D * moms,const ATOOLS::Flavour * flavs,
					    double weight, double ncount) 
{
  ATOOLS::Vec4D momsum = Vec4D(0.,0.,0.,0.);
  for (int i=0;i<nout;i++) {
    momsum += moms[i];
  }
  p_histo->Insert(momsum.Abs(),weight,ncount); 
}


void One_Particle_EVis::Evaluate(const Particle_List & plist,double weight,double ncount )
{
  ATOOLS::Vec4D momsum = Vec4D(0.,0.,0.,0.);
  for (Particle_List::const_iterator plit=plist.begin();plit!=plist.end();++plit) {
    momsum += (*plit)->Momentum();
  }
  p_histo->Insert(momsum.Abs(),weight,ncount); 
}
Primitive_Observable_Base * One_Particle_EVis::Copy() const
{
  return new One_Particle_EVis(m_flav,m_type,m_xmin,m_xmax,m_nbins,m_listname);
}
