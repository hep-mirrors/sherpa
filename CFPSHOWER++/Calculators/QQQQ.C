#include "CFPSHOWER++/Calculators/Gauge_Base.H"
#include "CFPSHOWER++/Shower/Kernel.H"
#include "MODEL/Main/Single_Vertex.H"
#include "ATOOLS/Phys/Flow.H"
#include "ATOOLS/Org/Message.H"

namespace CFPSHOWER {
  class QQQQ : public Gauge_Base {
    bool m_splitisanti;
  public:
    QQQQ(const Kernel_Info & info);
    virtual const double Scale(const Splitting & split) const;
    bool SetColours(Splitting & split);
  };
}


using namespace CFPSHOWER;
using namespace ATOOLS;

QQQQ::QQQQ(const Kernel_Info & info) :
  Gauge_Base(info)
{
  m_charge      = m_CF * m_TR;
  m_splitisanti = info.GetSplit().IsAnti();
  SetName("3 3 3 3");
}

const double QQQQ::Scale(const Splitting & split) const {
  double scale = split.t();
  switch (m_type) {
  case kernel_type::IF:
    break;
  case kernel_type::FI:
    break;
  case kernel_type::FF:
  default:
    break;
  }
  return scale;
}

bool QQQQ::SetColours(Splitting & split) {
  m_colors.clear();
  m_colors.resize(3);
  unsigned int newcol = Flow::Counter();
  unsigned int oldcol = split.GetSplitter()->GetColor()[int(m_splitisanti)]; 
  // first take care of parton j in splitter(aij) -> j + ai -> j[0] + a[1] + i[2]
  // assume ai has colour quantum numbers of gluon
  if (!m_splitisanti) {
    m_colors[0] = Color(newcol,0);
    m_colors[1] = Color(0,newcol);
    m_colors[2] = Color(oldcol,0);
  }
  else {
    m_colors[0] = Color(0,newcol);
    m_colors[1] = Color(newcol,0);
    m_colors[2] = Color(0,oldcol);
  }
  Flavour splitter     = split.GetKernel()->GetSplit();
  Flavour_Vector flavs = split.GetKernel()->GetFlavs();
  //msg_Out()<<"   *** "<<METHOD
  //	   <<"("<<splitter<<" -> "
  //	   <<flavs[0]<<" ("<<m_colors[0]<<") "
  //	   <<flavs[1]<<" ("<<m_colors[1]<<") "
  //	   <<flavs[2]<<" ("<<m_colors[2]<<")\n";
  return true;
}


DECLARE_GETTER(QQQQ,"QQQQ",Gauge_Base,Kernel_Info);

Gauge_Base * ATOOLS::Getter<Gauge_Base,Kernel_Info,QQQQ>::
operator()(const Parameter_Type & info) const
{
  if (info.GetFlavs().size()!=3) return NULL;
  if (abs(info.GetSplit().StrongCharge())==3 &&
      abs(info.GetFlavs()[0].StrongCharge())==3 &&
      abs(info.GetFlavs()[1].StrongCharge())==3 &&
      abs(info.GetFlavs()[2].StrongCharge())==3)
    return new QQQQ(info);
  return NULL;
}

void ATOOLS::Getter<Gauge_Base,Kernel_Info,QQQQ>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"gauge part for 3-3-3-3 in colour space";
}
