#include "One_Particle_Observables.H"

using namespace ANALYSIS;

#include "MyStrStream.H"

template <class Class>
Primitive_Observable_Base *const GetObservable(const String_Matrix &parameters)
{									
  if (parameters.size()<1) return NULL;
  if (parameters.size()==1) {
    if (parameters[0].size()<5) return NULL;
    int kf=ATOOLS::ToType<int>(parameters[0][0]);
    ATOOLS::Flavour flavour((ATOOLS::kf::code)abs(kf));
    if (kf<0) flavour=flavour.Bar();
    std::string list=parameters[0].size()>5?parameters[0][5]:"Analysed";
    return new Class(flavour,10*(int)(parameters[0][4]=="Log"),
		     ATOOLS::ToType<double>(parameters[0][1]),
		     ATOOLS::ToType<double>(parameters[0][2]),
		     ATOOLS::ToType<int>(parameters[0][3]),list);
  }
  else if (parameters.size()<5) return NULL;
  double min=0.0, max=1.0;
  size_t bins=100;
  std::string list="Analysed", scale="Lin";
  ATOOLS::Flavour flavour;
  for (size_t i=0;i<parameters.size();++i) {
    if (parameters[i].size()<2) continue;
    if (parameters[i][0]=="FLAV") {
      int kf=ATOOLS::ToType<int>(parameters[i][1]);
      flavour=ATOOLS::Flavour((ATOOLS::kf::code)abs(kf));
      if (kf<0) flavour=flavour.Bar();
    }
    else if (parameters[i][0]=="MIN") min=ATOOLS::ToType<double>(parameters[i][1]);
    else if (parameters[i][0]=="MAX") max=ATOOLS::ToType<double>(parameters[i][1]);
    else if (parameters[i][0]=="BINS") bins=ATOOLS::ToType<int>(parameters[i][1]);
    else if (parameters[i][0]=="SCALE") scale=parameters[i][1];
    else if (parameters[i][0]=="LIST") list=parameters[i][1];
  }
  return new Class(flavour,(scale=="Log")*10,min,max,bins,list);
}									

#define DEFINE_GETTER_METHOD(CLASS,NAME)				\
  Primitive_Observable_Base *const					\
  NAME::operator()(const String_Matrix &parameters) const		\
  { return GetObservable<CLASS>(parameters); }

#define DEFINE_PRINT_METHOD(NAME)					\
  void NAME::PrintInfo(std::ostream &str,const size_t width) const	\
  { str<<"kf min max bins Lin|Log [list]"; }

#define DEFINE_OBSERVABLE_GETTER(CLASS,NAME,TAG)			\
  DECLARE_GETTER(NAME,TAG,Primitive_Observable_Base,String_Matrix);	\
  DEFINE_GETTER_METHOD(CLASS,NAME);					\
  DEFINE_PRINT_METHOD(NAME)

using namespace ATOOLS;
using namespace std;

One_Particle_Observable_Base::One_Particle_Observable_Base(const Flavour & _flav,
							   int _type,double _xmin,double _xmax,int _nbins,
							   const std::string & _name) :
  Primitive_Observable_Base(_type,_xmin,_xmax,_nbins,NULL), 
  m_flav(_flav)
{
  m_name     = _name + m_flav.Name()+std::string(".dat");
  m_blobtype = std::string("");
  m_blobdisc = false;
}

void One_Particle_Observable_Base::Evaluate(double value,double weight, int ncount) 
{
  p_histo->Insert(value,weight,ncount); 
}

 
void One_Particle_Observable_Base::Evaluate(int nout,const ATOOLS::Vec4D * moms,const ATOOLS::Flavour * flavs,
					    double weight, int ncount) 
{
  for (int i=0;i<nout;i++) { if (flavs[i]==m_flav) Evaluate(moms[i],weight,ncount); }
}


void One_Particle_Observable_Base::Evaluate(const Particle_List & plist,double weight,int ncount )
{
  for (Particle_Const_Iterator plit=plist.begin();plit!=plist.end();++plit) {
    if ((*plit)->Flav()==m_flav) {
      Evaluate((*plit)->Momentum(),weight, ncount);
      return;
    }
  }
  
  Evaluate(Vec4D(1.,0,0,1.),0, ncount);
}



//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DEFINE_OBSERVABLE_GETTER(One_Particle_ET,One_Particle_ET_Getter,"ET");

One_Particle_ET::One_Particle_ET(const Flavour & _flav,
				 int _type,double _xmin,double _xmax,int _nbins,
				 const std::string & _name) :
  One_Particle_Observable_Base(_flav,_type,_xmin,_xmax,_nbins,"ET") { }


void One_Particle_ET::Evaluate(const Vec4D & mom,double weight,int ncount) 
{
  double pt2 = sqr(mom[1])+sqr(mom[2]);
  double p2  = sqr(mom[3])+pt2;
  double net = mom[0]*sqrt(pt2/p2);

  p_histo->Insert(net,weight,ncount); 
} 

Primitive_Observable_Base * One_Particle_ET::Copy() const
{
  return new One_Particle_ET(m_flav,m_type,m_xmin,m_xmax,m_nbins,m_name);
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DEFINE_OBSERVABLE_GETTER(One_Particle_PT,One_Particle_PT_Getter,"PT");

One_Particle_PT::One_Particle_PT(const Flavour & _flav,
				 int _type,double _xmin,double _xmax,int _nbins,
				 const std::string & _name) :
  One_Particle_Observable_Base(_flav,_type,_xmin,_xmax,_nbins,"PT") { }


void One_Particle_PT::Evaluate(const Vec4D & mom,double weight, int ncount) 
{
  double pt = sqrt(sqr(mom[1])+sqr(mom[2]));
  p_histo->Insert(pt,weight,ncount); 
} 

Primitive_Observable_Base * One_Particle_PT::Copy() const 
{
  return new One_Particle_PT(m_flav,m_type,m_xmin,m_xmax,m_nbins,m_name);
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DEFINE_OBSERVABLE_GETTER(One_Particle_Eta,One_Particle_Eta_Getter,"Eta");

One_Particle_Eta::One_Particle_Eta(const Flavour & _flav,
				   int _type,double _xmin,double _xmax,int _nbins,
				   const std::string & _name) :
  One_Particle_Observable_Base(_flav,_type,_xmin,_xmax,_nbins,"Eta") { }


void One_Particle_Eta::Evaluate(const Vec4D & mom,double weight, int ncount) 
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

Primitive_Observable_Base * One_Particle_Eta::Copy() const
{
  return new One_Particle_Eta(m_flav,m_type,m_xmin,m_xmax,m_nbins,m_name);
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DEFINE_OBSERVABLE_GETTER(One_Particle_E,One_Particle_E_Getter,"E");

One_Particle_E::One_Particle_E(const Flavour & _flav,
			       int _type,double _xmin,double _xmax,int _nbins,
			       const std::string & _name) :
  One_Particle_Observable_Base(_flav,_type,_xmin,_xmax,_nbins,"E") { }


void One_Particle_E::Evaluate(const Vec4D & mom,double weight, int ncount) 
{
  double E = mom[0];
  p_histo->Insert(E,weight,ncount); 
} 

Primitive_Observable_Base * One_Particle_E::Copy() const
{
  return new One_Particle_E(m_flav,m_type,m_xmin,m_xmax,m_nbins,m_name);
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DEFINE_OBSERVABLE_GETTER(One_Particle_EVis,One_Particle_EVis_Getter,"EVis");

One_Particle_EVis::One_Particle_EVis(const Flavour & _flav,
			       int _type,double _xmin,double _xmax,int _nbins,
			       const std::string & _name) :
  One_Particle_Observable_Base(_flav,_type,_xmin,_xmax,_nbins,"EVis") { }


void One_Particle_EVis::Evaluate(const Vec4D & mom,double weight, int ncount) 
{ } 

void One_Particle_EVis::Evaluate(int nout,const ATOOLS::Vec4D * moms,const ATOOLS::Flavour * flavs,
					    double weight, int ncount) 
{
  ATOOLS::Vec4D momsum = Vec4D(0.,0.,0.,0.);
  for (int i=0;i<nout;i++) {
    momsum += moms[i];
  }
  p_histo->Insert(momsum.Abs(),weight,ncount); 
}


void One_Particle_EVis::Evaluate(const Particle_List & plist,double weight,int ncount )
{
  ATOOLS::Vec4D momsum = Vec4D(0.,0.,0.,0.);
  for (Particle_Const_Iterator plit=plist.begin();plit!=plist.end();++plit) {
    momsum += (*plit)->Momentum();
  }
  p_histo->Insert(momsum.Abs(),weight,ncount); 
}
Primitive_Observable_Base * One_Particle_EVis::Copy() const
{
  return new One_Particle_EVis(m_flav,m_type,m_xmin,m_xmax,m_nbins,m_name);
}
