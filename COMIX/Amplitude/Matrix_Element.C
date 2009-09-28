#include "COMIX/Amplitude/Matrix_Element.H"

#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Math/Random.H"
#include "PHASIC++/Main/Color_Integrator.H"
#include "PHASIC++/Main/Helicity_Integrator.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/STL_Tools.H"
#include <iomanip>

using namespace COMIX;
using namespace ATOOLS;

Matrix_Element::Matrix_Element():
  m_nin(0), m_nout(0),
  p_colint(NULL), p_helint(NULL),
  m_mode(0), m_tests(0) {}

Matrix_Element::~Matrix_Element()
{
}

bool Matrix_Element::Initialize
(const size_t &nin,const size_t &nout,const std::vector<Flavour> &flavs,
 const double &isf,const double &fsf,Model *const model,
 const size_t &oew,const size_t &oqcd,
 const size_t &maxoew,const size_t &maxoqcd)
{
  m_nin=nin;
  m_nout=nout;
  m_ampl.SetOrderEW(oew);
  m_ampl.SetOrderQCD(oqcd);
  m_ampl.SetMaxOrderEW(Min(oew,maxoew));
  m_ampl.SetMaxOrderQCD(Min(oqcd,maxoqcd));
  Int_Vector incs(m_nin,1);
  incs.resize(flavs.size(),-1);
  if (!m_ampl.Construct(incs,flavs,model)) return false;
  std::map<Flavour,size_t> fc;
  for (size_t i(nin);i<flavs.size();++i) {
    std::map<Flavour,size_t>::iterator fit(fc.find(flavs[i]));
    if (fit==fc.end()) {
      fc[flavs[i]]=0;
      fit=fc.find(flavs[i]);
    }
    ++fit->second;
  }
  m_fsf=fsf;
  m_sf=m_fsf*isf;
  return true;
}

bool Matrix_Element::Map(const Matrix_Element *me,Flavour_Map &fmap)
{
  return m_ampl.Map(me->m_ampl,fmap);
}

void Matrix_Element::PrintStatistics
(std::ostream &str,const int mode) const
{
  m_ampl.PrintStatistics(str,mode);
}

double Matrix_Element::Differential(const std::vector<Vec4D> &momenta)
{
  return Differential(momenta,p_colint->I(),p_colint->J());
}

double Matrix_Element::Differential
(const std::vector<Vec4D> &momenta,
 const Int_Vector &ci,const Int_Vector &cj,const bool set)
{
  m_ampl.SetColors(ci,cj,set);
  m_ampl.SetMomenta(momenta);
  if (p_helint!=NULL && p_helint->On()) {
    m_ampl.Evaluate(p_helint->Chiralities());
    return m_ampl.Result()/m_sf;
  }
  m_ampl.EvaluateAll();
  return m_ampl.Result()/m_sf;
}

bool Matrix_Element::GaugeTest(std::vector<Vec4D> momenta)
{
  size_t nt(0);
  bool cnt(true);
  while (cnt) {
    while (!p_colint->GeneratePoint());
    m_ampl.SetColors(p_colint->I(),p_colint->J());
    if (!m_ampl.GaugeTest(momenta)) return false;
    if (m_ampl.Result()!=0.0) cnt=false;
    if (cnt) {
      if (++nt>100) 
	msg_Error()<<METHOD<<"(): Zero result. Redo gauge test."<<std::endl;
      while (!p_colint->GeneratePoint());
    }
  }
  return true;
}

void Matrix_Element::PrintGraphs(const std::string &gpath) const
{
  msg_Tracking()<<METHOD<<"(): Write diagrams to '"<<gpath<<"'.\n";
  m_ampl.WriteOutGraphs(gpath);
}

