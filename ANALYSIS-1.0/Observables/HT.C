#include "HT.H"
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
  void NAME::PrintInfo(ostream &str,const size_t width) const	\
  { str<<"min max bins Lin|LinErr|Log|LogErr [list]"; }

#define DEFINE_OBSERVABLE_GETTER(CLASS,NAME,TAG)			\
  DECLARE_GETTER(NAME,TAG,Primitive_Observable_Base,String_Matrix);	\
  DEFINE_GETTER_METHOD(CLASS,NAME)					\
  DEFINE_PRINT_METHOD(NAME)


DEFINE_OBSERVABLE_GETTER(HT,HT_Getter,"HT")

HT::HT(const String_Matrix & parameters):
    Primitive_Observable_Base(parameters)
{
  if (parameters.size()==1) {
    if (parameters[0].size()<4) {
      msg.Error()<<"Error in "<<METHOD<<": Not enough parameters for "
          <<"observable HT in Analysis.dat";
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
    m_name+="HT.dat";
  }
  else {
    if (m_listname=="") m_listname="Analysed";
    if (m_name=="SherpaDefault") {
      m_name="";
      if (m_listname!="Analysed") m_name=m_listname+string("_");
      m_name+="HT";
    }
    if (p_histo->Title()=="SherpaDefault") {
      p_histo->SetTitle("HT");
    }
  }
}

HT::HT(const HT * old) :
    Primitive_Observable_Base(*old) {}

// HT::HT(Histogram_Base * histo,const string & listname) :
//   Primitive_Observable_Base(histo)
// {
//   if (listname!="") {
//     m_listname = listname;
//     m_name = listname+"_HT.dat";
//   }
//   else
//     m_name = "HT.dat";
// }

void HT::Evaluate(const Particle_List& pl,
		  double weight, int ncount)
{
  Particle_List* jets=p_ana->GetParticleList(m_listname);
  double HT=0.0;
  if(jets->size()==0) {
    p_histo->Insert(0.0,0.0,ncount);
    return;
  }
  for (Particle_List::const_iterator pit=jets->begin();
       pit!=jets->end();++pit) {
    HT+=(*pit)->Momentum().PPerp();
  }
  p_histo->Insert(HT,weight,ncount);
}


Primitive_Observable_Base * HT::Copy() const 
{
  return new HT(this);
}
