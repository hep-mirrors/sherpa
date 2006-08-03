#include "One_Particle_Observables.H"
#include "MyStrStream.H"

using namespace ANALYSIS;
using namespace ATOOLS;
using namespace std;

#define DEFINE_GETTER_METHOD(CLASS,NAME)				\
  Primitive_Observable_Base *					        \
  NAME::operator()(const String_Matrix &parameters) const		\
  { return new CLASS(parameters); }

#define DEFINE_PRINT_METHOD(NAME)					\
  void NAME::PrintInfo(std::ostream &str,const size_t width) const	\
  { str<<"kf min max bins Lin|LinErr|Log|LogErr [list]"; }

#define DEFINE_CONSTRUCTOR_METHODS(CLASS,TAG)                           \
  CLASS::CLASS(const String_Matrix & parameters) :                      \
      One_Particle_Observable_Base(parameters, TAG) { }                 \
                                                                        \
  CLASS::CLASS(const CLASS * old) :                                     \
      One_Particle_Observable_Base(*old) { }                            \
                                                                        \
  Primitive_Observable_Base * CLASS::Copy() const { return new CLASS(this); }

#define DEFINE_OBSERVABLE_GETTER(CLASS,NAME,TAG)			\
  DECLARE_GETTER(NAME,TAG,Primitive_Observable_Base,String_Matrix);	\
  DEFINE_GETTER_METHOD(CLASS,NAME)					\
  DEFINE_PRINT_METHOD(NAME)                                             \
  DEFINE_CONSTRUCTOR_METHODS(CLASS,TAG)

One_Particle_Observable_Base::One_Particle_Observable_Base(const String_Matrix & parameters,
                                                           const std::string & obsname):
    Primitive_Observable_Base(parameters)
{
  if (parameters.size()==1) {
    if (parameters[0].size()<5) {
      msg.Error()<<"Error in "<<METHOD<<": Not enough parameters for "
          <<"observable "<<obsname<<" in Analysis.dat";
      abort();
    }
    int kf = ATOOLS::ToType<int>(parameters[0][0]);
    m_flav = Flavour( (ATOOLS::kf::code)abs(kf), kf<0 );

    m_xmin  = ATOOLS::ToType<double>(parameters[0][1]);
    m_xmax  = ATOOLS::ToType<double>(parameters[0][2]);
    int nbins = ATOOLS::ToType<int>(parameters[0][3]);
    m_type  = HistogramType(parameters[0][4]);
    p_histo = new Histogram(m_type,m_xmin,m_xmax,nbins);

    m_listname = parameters[0].size()>5?parameters[0][5]:"Analysed";

    MyStrStream str;
    str<<obsname<<m_flav<<".dat";
    str>>m_name;
  }
  else {
    for (size_t i=0;i<parameters.size();++i) {
      if (parameters[i].size()<2) continue;
      if (parameters[i][0]=="FLAV") {
        int kf=ATOOLS::ToType<int>(parameters[i][1]);
        m_flav=ATOOLS::Flavour((ATOOLS::kf::code)abs(kf),kf<0);
      }
    }

    if (m_name=="SherpaDefault") {
      MyStrStream str;
      str<<obsname<<m_flav;
      str>>m_name;
    }
    if (p_histo->Title()=="SherpaDefault") {
      std::string title = "";
      MyStrStream str;
      str<<obsname<<" of "<<m_flav;
      str>>title;
      p_histo->SetTitle(title);
    }
    if (m_listname=="") m_listname="Analysed";
  }
}

One_Particle_Observable_Base::One_Particle_Observable_Base(const One_Particle_Observable_Base * old) :
    Primitive_Observable_Base(*old)
{
  m_flav = old->m_flav;
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
  for (Particle_List::const_iterator plit=plist.begin();plit!=plist.end();++plit) {
    if ((*plit)->Flav()==m_flav) {
      Evaluate((*plit)->Momentum(),weight, ncount);
      return;
    }
  }
  
  Evaluate(Vec4D(1.,0,0,1.),0, ncount);
}



//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DEFINE_OBSERVABLE_GETTER(One_Particle_ET,One_Particle_ET_Getter,"ET")

void One_Particle_ET::Evaluate(const Vec4D & mom,double weight,int ncount) 
{
  double pt2 = sqr(mom[1])+sqr(mom[2]);
  double p2  = sqr(mom[3])+pt2;
  double net = mom[0]*sqrt(pt2/p2);

  p_histo->Insert(net,weight,ncount); 
} 

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DEFINE_OBSERVABLE_GETTER(One_Particle_PT,One_Particle_PT_Getter,"PT")

void One_Particle_PT::Evaluate(const Vec4D & mom,double weight, int ncount) 
{
  double pt = sqrt(sqr(mom[1])+sqr(mom[2]));
  p_histo->Insert(pt,weight,ncount); 
} 

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DEFINE_OBSERVABLE_GETTER(One_Particle_Eta,One_Particle_Eta_Getter,"Eta")

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

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DEFINE_OBSERVABLE_GETTER(One_Particle_E,One_Particle_E_Getter,"E")

void One_Particle_E::Evaluate(const Vec4D & mom,double weight, int ncount) 
{
  double E = mom[0];
  p_histo->Insert(E,weight,ncount); 
} 

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DEFINE_OBSERVABLE_GETTER(One_Particle_BeamAngle,One_Particle_BeamAngle_Getter,"BeamAngle")

void One_Particle_BeamAngle::Evaluate(const Vec4D & mom,double weight, int ncount) 
{
  double ct = mom.CosTheta();
  p_histo->Insert(ct,weight,ncount); 
} 

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DEFINE_OBSERVABLE_GETTER(One_Particle_EVis,One_Particle_EVis_Getter,"EVis")

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
  for (Particle_List::const_iterator plit=plist.begin();plit!=plist.end();++plit) {
    momsum += (*plit)->Momentum();
  }
  p_histo->Insert(momsum.Abs(),weight,ncount); 
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DEFINE_OBSERVABLE_GETTER(One_Particle_Px,One_Particle_Px_Getter,"Px")

void One_Particle_Px::Evaluate(const Vec4D & mom,double weight, int ncount) 
{
  double px = mom[1];
  p_histo->Insert(px,weight,ncount); 
} 

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DEFINE_OBSERVABLE_GETTER(One_Particle_Py,One_Particle_Py_Getter,"Py")

void One_Particle_Py::Evaluate(const Vec4D & mom,double weight, int ncount) 
{
  double py = mom[2];
  p_histo->Insert(py,weight,ncount); 
} 

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DEFINE_OBSERVABLE_GETTER(One_Particle_Pz,One_Particle_Pz_Getter,"Pz")

void One_Particle_Pz::Evaluate(const Vec4D & mom,double weight, int ncount) 
{
  double pz = mom[3];
  p_histo->Insert(pz,weight,ncount); 
} 
