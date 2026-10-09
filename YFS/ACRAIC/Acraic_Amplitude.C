#include "YFS/ACRAIC/Acraic_Amplitude.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/MyStrStream.H"

using namespace YFS;
using namespace ATOOLS;

int Amplitude::s_nlegs = 4;

std::vector<int> Amplitude::s_bit;

void Amplitude::SetLegs(const ATOOLS::Flavour_Vector &flavs)
{
  s_bit.assign(flavs.size(), -1);
  int nb(0);
  for (size_t i(0); i < flavs.size(); ++i) {
    // METOOLS' rule: massless vector 2 states, otherwise 2s+1
    const int ns(flavs[i].IsVector() && !flavs[i].IsMassive()
                 ? 2 : flavs[i].IntSpin() + 1);
    if (ns > 2)
      THROW(not_implemented, "ACRAIC: external "+flavs[i].IDName()+" has "
            +ATOOLS::ToString(ns)+" helicity states; the amplitudes hold two"
            " per leg. Let massive vectors decay in the process.");
    if (ns > 1) s_bit[i] = nb++;
  }
  SetLegs(nb);
}

void Amplitude::SetLegs(int n) {
  if (n < 1 || n > s_maxlegs)
    THROW(not_implemented, "ACRAIC: "+ATOOLS::ToString(n)+" Born legs with two"
          " helicity states; the capacity is "+ATOOLS::ToString(s_maxlegs)
          +" (2 -> 8). Raise Amplitude::s_maxlegs.");
  if (n > s_warnlegs && n != s_nlegs)
    msg_Error()<<"ACRAIC: "<<n<<" Born legs with two helicity states, beyond the"
               <<" validated "<<s_warnlegs<<" (2 -> 4 plus a Born photon): every"
               <<" amplitude has 2^"<<n<<" = "<<NumLegSubsets<int>(n)
               <<" helicity entries, the cost grows accordingly, and the"
               <<" multiplicity is untested."<<std::endl;
  s_nlegs = n;
}

Amplitude::Amplitude(const Amplitude &o) {
  const int n(NHel());
  for (int i = 0; i < n; ++i) m_A[i] = o.m_A[i];
}

Amplitude &Amplitude::operator=(const Amplitude &o) {
  if (this == &o) return *this;
  const int n(NHel());
  for (int i = 0; i < n; ++i) m_A[i] = o.m_A[i];
  return *this;
}

Amplitude::Amplitude() {
  // only the active extent: one Amplitude is built per partition
  const int n(NHel());
  for (int i = 0; i < n; ++i) m_A[i] = Complex(0, 0);
}
