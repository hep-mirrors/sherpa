#include "CFPSHOWER++/SplittingFunctions/Gauge_Base.H"
#include "CFPSHOWER++/Shower/Kernel.H"
#include "MODEL/Main/Single_Vertex.H"
#include "ATOOLS/Phys/Flow.H"
#include "ATOOLS/Org/Message.H"

namespace CFPSHOWER {
  class QQG : public Gauge_Base {
    int m_mode;
  public:
    QQG(const Kernel_Info & info) : Gauge_Base(info) {
      m_charge = m_CF;
      m_colors.resize(2);
      SetName("3-3-8");
    }
    inline const double Charge(const double & scale) const { return m_charge; }    
    const double Scale(const Splitting & split) const;
    bool         SetColours(Splitting & split);
  };
}


using namespace CFPSHOWER;
using namespace ATOOLS;

const double QQG::Scale(const Splitting & split) const {
  double kt2 = (split.Mom(0)*split.Mom(1))*(split.Mom(1)*split.Mom(2))/(split.Mom(0)*split.Mom(2));
  switch (m_type) {
  case kernel_type::FF:
    if (m_muRscheme==muR_scheme::KT2) return kt2;
    break;
  default:
    break;
  }
  return split.T();
}

bool QQG::SetColours(Splitting & split) {
  unsigned int newcol  = Flow::Counter();
  Flavour_Vector flavs = split.GetKernel()->GetFlavs();
  if (!flavs[0].IsAnti()) {
    m_colors[0] = Color(newcol,0);
    m_colors[1] = Color(split.GetSplitter()->GetColor()[0],newcol);
  }
  else if (flavs[0].IsAnti()) {
    m_colors[0] = Color(0,newcol);
    m_colors[1] = Color(newcol,split.GetSplitter()->GetColor()[1]);
  }
  return true;
}

DECLARE_GETTER(QQG,"QQG",Gauge_Base,Kernel_Info);

Gauge_Base * ATOOLS::Getter<Gauge_Base,Kernel_Info,QQG>::
operator()(const Parameter_Type & info) const
{
  if (abs(info.GetSplit().StrongCharge())==3 &&
      abs(info.GetFlavs()[0].StrongCharge())==3 &&
      abs(info.GetFlavs()[1].StrongCharge())==8)
    return new QQG(info);
  return NULL;
}

void ATOOLS::Getter<Gauge_Base,Kernel_Info,QQG>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"gauge part for 3-3-8 in colour space";
}
