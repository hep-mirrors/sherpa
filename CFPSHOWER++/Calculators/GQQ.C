#include "CFPSHOWER++/Calculators/Gauge_Base.H"
#include "CFPSHOWER++/Shower/Kernel.H"
#include "MODEL/Main/Single_Vertex.H"
#include "ATOOLS/Phys/Flow.H"
#include "ATOOLS/Org/Message.H"

namespace CFPSHOWER {
  class GQQ : public Gauge_Base {
    int m_mode;
  public:
    GQQ(const Kernel_Info & info) : Gauge_Base(info) {
      m_charge = m_type==kernel_type::IF?m_CF/2.:m_TR/2.;
      m_colors.resize(2);
      SetName("8-3-3");
    }
    
    const double Scale(const Splitting & split) const;
    bool SetColours(Splitting & split);
  };
}


using namespace CFPSHOWER;
using namespace ATOOLS;

const double GQQ::Scale(const Splitting & split) const {
  switch (m_type) {
  case kernel_type::IF:
    return split.t(0)/split.z(0);
    break;
  case kernel_type::FI:
    return split.t(0)/split.y();
    break;
  case kernel_type::FF:
  default:
    break;
  }
  return split.t(0);
}

bool GQQ::SetColours(Splitting & split) {
  Flavour_Vector flavs = split.GetKernel()->GetFlavs();
  if (!flavs[0].IsAnti() && flavs[1].IsAnti()) {
    m_colors[0] = Color(split.GetSplitter()->GetColor()[0],0);
    m_colors[1] = Color(0,split.GetSplitter()->GetColor()[1]);
  }
  else if (flavs[0].IsAnti() && !flavs[1].IsAnti()) {
    m_colors[0] = Color(0,split.GetSplitter()->GetColor()[1]);
    m_colors[1] = Color(split.GetSplitter()->GetColor()[0],0);
  }
  else return false;
  return true;
}

DECLARE_GETTER(GQQ,"GQQ",Gauge_Base,Kernel_Info);

Gauge_Base * ATOOLS::Getter<Gauge_Base,Kernel_Info,GQQ>::
operator()(const Parameter_Type & info) const
{
  if (abs(info.GetSplit().StrongCharge())==8 &&
      abs(info.GetFlavs()[0].StrongCharge())==3 &&
      abs(info.GetFlavs()[1].StrongCharge())==3)
    return new GQQ(info);
  return NULL;
}

void ATOOLS::Getter<Gauge_Base,Kernel_Info,GQQ>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"gauge part for 8-3-3 in colour space";
}
