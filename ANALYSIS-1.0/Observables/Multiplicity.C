#include "Multiplicity.H"
#include "Message.H"
#include "MathTools.H"
#include "Particle_Qualifier.H"
#include "MyStrStream.H"

using namespace ANALYSIS;
using namespace ATOOLS;
using namespace std;

#define DEFINE_GETTER_METHOD(CLASS,NAME)                                \
  Primitive_Observable_Base *                                   \
  NAME::operator()(const String_Matrix &parameters) const               \
{ return new CLASS(parameters); }

#define DEFINE_PRINT_METHOD(NAME)                                       \
  void NAME::PrintInfo(ostream &str,const size_t width) const   \
{ str<<"min max bins Lin|LinErr|Log|LogErr [list]"; }

#define DEFINE_OBSERVABLE_GETTER(CLASS,NAME,TAG)                        \
  DECLARE_GETTER(NAME,TAG,Primitive_Observable_Base,String_Matrix);     \
  DEFINE_GETTER_METHOD(CLASS,NAME)                                      \
  DEFINE_PRINT_METHOD(NAME)

DEFINE_OBSERVABLE_GETTER(Multiplicity,Multiplicity_Getter,"Multi")

Multiplicity::Multiplicity(const String_Matrix & parameters):
    Primitive_Observable_Base(parameters)
{
  if (parameters.size()==1) {
    if (parameters[0].size()<4) {
      msg.Error()<<"Error in "<<METHOD<<": Not enough parameters for "
          <<"observable Multi in Analysis.dat";
      abort();
    }
    m_xmin  = ToType<double>(parameters[0][0]);
    m_xmax  = ToType<double>(parameters[0][1]);
    int nbins = ToType<int>(parameters[0][2]);
    m_type  = HistogramType(parameters[0][3]);
    p_histo = new Histogram(m_type,m_xmin,m_xmax,nbins);

    m_listname = parameters[0].size()>4?parameters[0][4]:"Analysed";
    m_name="";
    if (m_listname!="Analysed") m_name=m_listname+string("_");
    m_name+="Multi.dat";
  }
  else {
    if (m_listname=="") m_listname="Analysed";
    if (m_name=="SherpaDefault") {
      m_name="";
      if (m_listname!="Analysed") m_name=m_listname+string("_");
      m_name+="Multi";
    }
    if (p_histo->Title()=="SherpaDefault") {
      string title = "Multiplicity";
      p_histo->SetTitle(title);
    }
  }
}

Multiplicity::Multiplicity(const Multiplicity * old) :
  Primitive_Observable_Base(*old)
{
  m_listname = old->m_listname; // ?
  m_name     = old->m_name;     // ?
}

Primitive_Observable_Base * Multiplicity::Copy() const { return new Multiplicity(this); }

void Multiplicity::Evaluate(const Particle_List & pl,
			    double weight, int ncount)
{
  p_histo->Insert(pl.size(),weight,ncount);
}





