#include "METOOLS/HadronCurrents/Photon_Emission.H"
#include "METOOLS/Main/Polarization_Tools.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/MyStrStream.H"
#include <cmath>

using namespace METOOLS;
using namespace ATOOLS;

void Photon_Emission::SetAlpha(const double alpha)
{
  if (!(alpha>0.))
    THROW(fatal_error,"Photon_Emission: alpha must be positive, got "
          +ToString(alpha)+".");
  m_e = sqrt(4.*M_PI*alpha);
}

double Photon_Emission::E() const
{
  if (m_e<0.)
    THROW(fatal_error,"Photon_Emission: a current was asked to emit a photon "
          "before its coupling was set. Exchange_ME sets it when it builds "
          "the current; a current built some other way has to call SetAlpha.");
  return m_e;
}

// XYZFunc gives a leg the mass sqrt(p^2) whenever p^2 > 0, as evaluated by
// Vec4D::Abs2(). Setting k0 = |k| leaves p^2 = +-1e-16 E^2 after rounding,
// and when it lands positive the pseudo-leg acquires a "mass" of order
// 1e-8 E: k-slash u(k) is then not zero but ~1e-8, which showed up as a
// second population of gauge-check failures at 1e-9..1e-7, at hard,
// wide-angle points where nothing else is delicate. Step k0 down by ulps
// until p^2 is not positive -- tested as l*l, the dot product XYZFunc
// itself evaluates, AND as Abs2(), which rounds differently. Testing only
// Abs2() removed most of the effect but not all of it.
Vec4D Photon_Emission::LightLike(const Vec4D & k) const
{
  Vec4D l(Vec3D(k).Abs(),k[1],k[2],k[3]);
  while (l*l>0. || l.Abs2()>0.) l[0] = std::nextafter(l[0],0.);
  return l;
}

Vec4C Photon_Emission::EpsStar(const Vec4D & k,const int lambda) const
{
  if (m_gaugetest) return Vec4C(Complex(k[0],0.),Complex(k[1],0.),
                                Complex(k[2],0.),Complex(k[3],0.));
  // Outgoing photon: the Feynman rule wants the complex conjugate. The same
  // vector, from the same constructor, in every current -- see the class
  // comment for why that matters.
  return conj(Polarization_Vector(k)[lambda]);
}

Vec4C Photon_Emission::Emit(const Flavour & f,const Vec4D & p,
                            const bool incoming,const Vec4D & k,
                            const Vec4C & eps,const Vec4C & J0,
                            const Vec4C & S) const
{
  const double Q = f.Charge();
  if (Q==0.) return Vec4C(0.,0.,0.,0.);
  const double  eta = incoming ? -1. : 1.;
  const double  pk  = p*k;
  const Complex pe  = p*eps;
  return (E()*Q)*((eta*pe/pk)*J0 + S/(2.*pk));
}
