#include "CFPSHOWER++/SplittingFunctions/Gauge_Base.H"
#include "MODEL/Main/Single_Vertex.H"
#include "ATOOLS/Phys/Flow.H"
#include "ATOOLS/Org/Message.H"

namespace CFPSHOWER {
  class GGG : public Gauge_Base {
  public:
    GGG(const Kernel_Info & info) : Gauge_Base(info) {
      m_charge = 3./2.;
      SetName("8-8-8");
      m_colors.resize(2);
    }
    const double Scale(const Splitting & split) const;
    bool SetColours(Splitting & split);
  };
}

using namespace CFPSHOWER;
using namespace ATOOLS;

const double GGG::Scale(const Splitting & split) const {
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

bool GGG::SetColours(Splitting & split) {
  Parton * splitter = split.GetSplitter(), * spec = split.GetSpectator();
  // check for colour connections between splitter and spectator, make sure
  // there is at least one.
  bool same0 = splitter->GetColor()[0]==spec->GetColor()[1];
  bool same1 = splitter->GetColor()[1]==spec->GetColor()[0];
  if (!same0 && !same1) return false;
  // if double colour-connected (gluon) pair, choose one connection.
  if (same0 && same1) {
    if (ran->Get()>0.5) same0 = false; else same1 = false;
  }
  unsigned int newcol = Flow::Counter();
  // starting position: splitter anti-colour = spectator colour:
  // soft gluon inherits anti-colour and hard gluon keeps colour
  m_colors[0] = Color(splitter->GetColor()[0],newcol);
  m_colors[1] = Color(newcol,splitter->GetColor()[1]);
  // if splitter colour = spectator anti-colour:
  // soft gluon inherits colour and hard gluon keep anti-colour,
  // this is equivalent to swapping the default.
  if (same0)  swap(m_colors[0],m_colors[1]);
  // if we use the "swapped" splitting function, the roles of soft and
  // hard gluon are interchange - in other words, colours have to swap.
  if (m_tagsequence[0]==1) swap(m_colors[0],m_colors[1]);
  return true;
}

DECLARE_GETTER(GGG,"GGG",Gauge_Base,Kernel_Info);

Gauge_Base * ATOOLS::Getter<Gauge_Base,Kernel_Info,GGG>::
operator()(const Parameter_Type & info) const
{
  if (abs(info.GetSplit().StrongCharge())==8 &&
      abs(info.GetFlavs()[0].StrongCharge())==8 &&
      abs(info.GetFlavs()[1].StrongCharge())==8) return new GGG(info);
  return NULL;
}

void ATOOLS::Getter<Gauge_Base,Kernel_Info,GGG>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"gauge part for 8-8-8 in colour space";
}
