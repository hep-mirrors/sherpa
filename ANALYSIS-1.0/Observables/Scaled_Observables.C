#include "Scaled_Observables.H"
#include "Run_Parameter.H"
#include "MyStrStream.H"

using namespace ANALYSIS;
using namespace ATOOLS;
using namespace std;

#define DEFINE_GETTER_METHOD(CLASS,NAME)                                \
  Primitive_Observable_Base *                                   \
  NAME::operator()(const String_Matrix &parameters) const               \
  { return new CLASS(parameters); }

#define DEFINE_PRINT_METHOD(NAME)                                       \
  void NAME::PrintInfo(std::ostream &str,const size_t width) const      \
  { str<<"min max bins Lin|LinErr|Log|LogErr [list [ref]]"; }

#define DEFINE_CONSTRUCTOR_METHODS(CLASS,TAG)                           \
  CLASS::CLASS(const String_Matrix & parameters) :                      \
      Scaled_Observable_Base(parameters, TAG) { }               \
                                                                        \
  CLASS::CLASS(const CLASS * old) :                                     \
      Scaled_Observable_Base(*old) { }                            \
                                                                        \
  Primitive_Observable_Base * CLASS::Copy() const { return new CLASS(this); }

#define DEFINE_OBSERVABLE_GETTER(CLASS,NAME,TAG)                        \
  DECLARE_GETTER(NAME,TAG,Primitive_Observable_Base,String_Matrix);     \
  DEFINE_GETTER_METHOD(CLASS,NAME)                                      \
  DEFINE_PRINT_METHOD(NAME)                                             \
  DEFINE_CONSTRUCTOR_METHODS(CLASS,TAG)

Scaled_Observable_Base::Scaled_Observable_Base(const String_Matrix & parameters,
                                               const std::string & obsname):
    Primitive_Observable_Base(parameters)
{
  if (parameters.size()==1) {
    if (parameters[0].size()<4) {
      msg.Error()<<"Error in "<<METHOD<<": Not enough parameters for "
                 <<"observable "<<obsname<<" in Analysis.dat";
      abort();
    }

    m_xmin  = ToType<double>(parameters[0][0]);
    m_xmax  = ToType<double>(parameters[0][1]);
    int nbins = ToType<int>(parameters[0][2]);
    m_type  = HistogramType(parameters[0][3]);
    p_histo = new Histogram(m_type,m_xmin,m_xmax,nbins);

    m_listname = parameters[0].size()>4?parameters[0][4]:"Analysed";

    m_ecms=parameters[0].size()>5 ? ToType<double>(parameters[0][5]) : rpa.gen.Ecms();

    m_name="";
    if (m_listname!="Analysed") m_name=m_listname+string("_");
    m_name+=obsname+".dat";
  }
  else {
    if (m_listname=="") m_listname="Analysed";
    m_ecms = rpa.gen.Ecms();
    for (size_t i=0;i<parameters.size();++i) {
      if (parameters[i].size()<2) continue;
      if (parameters[i][0]==std::string("REF")) m_ecms = ToType<double>(parameters[i][1]);
    }

    if (m_name=="SherpaDefault") {
      m_name="";
      if (m_listname!="Analysed") m_name=m_listname+string("_");
      m_name+=obsname;
    }
    if (p_histo->Title()=="SherpaDefault") {
      p_histo->SetTitle(obsname);
    }
  }
}

Scaled_Observable_Base::Scaled_Observable_Base(const Scaled_Observable_Base * old) :
    Primitive_Observable_Base(*old)
{
  m_ecms=old->m_ecms;
}

void Scaled_Observable_Base::Evaluate(double value,double weight, int ncount) 
{
  p_histo->Insert(value,weight,ncount); 
}

 
void Scaled_Observable_Base::Evaluate(int nout,const ATOOLS::Vec4D * moms,
					    double weight, int ncount) 
{
  for (int i=0;i<nout;i++) Evaluate(moms[i],weight,ncount);
}


void Scaled_Observable_Base::Evaluate(const Particle_List & plist,double weight,int ncount )
{
  for (Particle_List::const_iterator plit=plist.begin();plit!=plist.end();++plit) {
    Evaluate((*plit)->Momentum(),weight, ncount);
  }
}



//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DEFINE_OBSERVABLE_GETTER(Scaled_Momentum,Scaled_Momentum_Getter,"XP")

void Scaled_Momentum::Evaluate(const Vec4D & mom,double weight,int ncount) 
{
  double xp = 2.*Vec3D(mom).Abs()/m_ecms;

  p_histo->Insert(xp,weight,ncount); 
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DEFINE_OBSERVABLE_GETTER(Log_Scaled_Momentum,Log_Scaled_Momentum_Getter,"LogXP")

void Log_Scaled_Momentum::Evaluate(const Vec4D & mom,double weight,int ncount)
{
  double xp = 2.*Vec3D(mom).Abs()/m_ecms;
  double xi = - log(xp);

  p_histo->Insert(xi,weight,ncount); 
} 

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DEFINE_OBSERVABLE_GETTER(Scaled_Energy,Scaled_Energy_Getter,"XE")

void Scaled_Energy::Evaluate(const Vec4D & mom,double weight, int ncount)
{
  double E = 2.*mom[0]/m_ecms;
  p_histo->Insert(E,weight,ncount); 
} 
