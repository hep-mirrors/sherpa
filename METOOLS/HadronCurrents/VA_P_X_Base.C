#include "METOOLS/HadronCurrents/VA_P_X_Base.H"
#include "METOOLS/Main/Polarization_Tools.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"

using namespace METOOLS;
using namespace ATOOLS;
using namespace std;

namespace {
  /// Strange mesons that can appear as X in a charm semileptonic
  /// decay. Used only to pick the DEFAULT CKM element; the "Vcq" YAML
  /// key always wins, and every channel in the shipped Decaydata.yaml
  /// sets it explicitly, so this is a convenience, not a dependency.
  bool IsStrangeMeson(const int kf) {
    switch (kf) {
    case 321: case 311: case 310: case 130:      // K+, K0, K_S, K_L
    case 323: case 313:                          // K*(892)
    case 325: case 315:                          // K*_2(1430)
    case 10321: case 10311:                      // K*_0(1430)
    case 10323: case 10313:                      // K_1(1270)
    case 20323: case 20313:                      // K_1(1400)
    case 100323: case 100313:                    // K*(1410)
      return true;
    default: return false;
    }
  }
}

VA_P_X_Base::VA_P_X_Base(const ATOOLS::Flavour_Vector & flavs,
                         const std::vector<int> & indices,
                         const std::string & name) :
  Current_Base(flavs,indices,name),
  m_M(0.), m_M2(0.), m_mX(0.),
  m_norm(1.), m_Vcq(1.), m_isospin(1.),
  m_isKSKL(false), m_ckmdressed(false),
  m_epssign(1.), m_dynamicmass(false), m_q2min(1.e-8), m_nX(0)
{
  // Index 0 is the decaying meson; everything after it is hadronic.
  // The base class only requires that there IS a parent and at least
  // one recoiling hadron - the "exactly two" check that used to live
  // here predated VA_P_PP and fired in its constructor before the
  // derived class could say anything. Each derived current asserts its
  // own arity in CheckArity(), so the error now names the right class
  // and the right expectation.
  if (p_i.size()<2)
    THROW(fatal_error,"P->X current "+name+" needs at least two indices "
          "(decaying meson and at least one recoiling hadron), got "
          +ToString(p_i.size())+".");
  m_nX = p_i.size()-1;
  m_M  = m_flavs[p_i[0]].HadMass();
  m_M2 = m_M*m_M;
  // Only meaningful for a single recoiling hadron; multi-hadron
  // currents carry their own daughter masses and ignore m_mX.
  m_mX = m_flavs[p_i[1]].HadMass();
  int kfX = int(m_flavs[p_i[1]].Kfcode());
  m_isKSKL = (m_nX==1 && (kfX==310 || kfX==130));
}

void VA_P_X_Base::CheckArity(const size_t & nX) const {
  if (m_nX!=nX)
    THROW(fatal_error,"Current "+m_name+" expects "+ToString(nX)+
          " recoiling hadron(s) plus the decaying meson, i.e. "
          +ToString(nX+1)+" indices, but was given "
          +ToString(p_i.size())+".");
}

VA_P_X_Base::~VA_P_X_Base() {
  for (map<string,FF_P_X*>::iterator it=m_ffs.begin();it!=m_ffs.end();++it)
    if (it->second) delete it->second;
  m_ffs.clear();
}

double VA_P_X_Base::DefaultVcq() const {
  const bool parentstrange = IsStrangeMeson(int(m_flavs[p_i[0]].Kfcode()))
    || int(m_flavs[p_i[0]].Kfcode())==431;
  const bool daughterstrange = IsStrangeMeson(int(m_flavs[p_i[1]].Kfcode()));
  // c -> s if the |Delta S|=1 current is the one connecting the two;
  // for a D_s parent the roles are reversed because the sbar is the
  // spectator there.
  if (parentstrange) return (daughterstrange ? Tools::Vcd : Tools::Vcs);
  return (daughterstrange ? Tools::Vcs : Tools::Vcd);
}

FF_P_X * VA_P_X_Base::MakeFF(const std::string & name,GeneralModel & model) {
  // FF_P_X is a TWO-index object by construction: it exists to supply
  // an invariant form factor F(q^2) of a P -> X transition, and its
  // getter rejects anything else. Passing p_i wholesale therefore
  // worked for every quasi-two-body current and failed for VA_P_PP,
  // whose p_i has three entries.
  //
  // The right two indices are the parent and ONE hadron: the form
  // factor never sees the hadronic system as such. For a multi-hadron
  // current the second index only fixes a nominal daughter mass, which
  // is unused there because VA_P_PP forces DYNAMIC_MASS on for every
  // form factor and the kernels pass sqrt(s) explicitly to Value().
  std::vector<int> fi;
  fi.push_back(p_i[0]);
  fi.push_back(p_i[1]);
  map<string,double> pmap;
  FF_Parameters params(ff_model::none,m_flavs,fi,pmap,name,&model);
  FormFactor_Base * ff = FF_Getter::GetObject("FF_P_X",params);
  if (ff==NULL)
    THROW(fatal_error,"Could not build form factor '"+name+"' for current "
          +m_name+". FF_Getter returned NULL for tag 'FF_P_X' with "
          +ToString(fi.size())+" indices - check that FF_P_X.C is built "
          "and registered.");
  FF_P_X * ffx = dynamic_cast<FF_P_X*>(ff);
  if (ffx==NULL) THROW(fatal_error,"Form factor '"+name+"' is not an FF_P_X.");
  m_ffs[name] = ffx;
  return ffx;
}

void VA_P_X_Base::ReadCommonParameters(GeneralModel & model) {
  m_ckmdressed  = (model("CKM_DRESSED",0.)>0.5);
  m_Vcq         = model("Vcq",DefaultVcq());
  if (m_ckmdressed) m_Vcq = 1.;   // |V_cq| already inside F(0)
  m_isospin     = model("ISOSPIN",1.);
  m_epssign     = model("EPS_SIGN",1.);
  m_dynamicmass = (model("DYNAMIC_MASS",0.)>0.5);
  m_q2min       = model("Q2MIN",1.e-8);
  m_norm        = m_Vcq*m_isospin;
  // K_S/K_L are strangeness-eigenstate MIXTURES; the weak current only
  // ever produces one flavour eigenstate, so projecting onto K_S or
  // K_L costs an extra 1/sqrt(2) in the amplitude. Exactly the same
  // m_isKSKL logic as in FF_0_PP.C.
  if (m_isKSKL) m_norm *= SQRT_05;
}

double VA_P_X_Base::
RegularisedOverQ2(const std::function<double(double)> & N,
                  const double & q2) const {
  if (std::abs(q2)>m_q2min) return N(q2)/q2;
  // N(0)=0 holds analytically by the q^2=0 kinematic constraint
  // (f+(0)=f0(0), A_0(0)=A_3(0), V_0(0)=V_3(0)), so the limit is
  // N'(0). Symmetric difference: exact to O(h^2) and, unlike a
  // closed-form derivative, valid for every shape in FF_P_X.
  const double h = 1.e-4;
  return (N(h)-N(-h))/(2.*h);
}

void VA_P_X_Base::
PolarizationVectors(const ATOOLS::Vec4D & p,
                    std::vector<ATOOLS::Vec4C> & eps) const {
  // Polarization_Vector from METOOLS/Main/Polarization_Tools.H, used
  // with the same call signature as the legacy HADRONS++ currents. The
  // basis MUST be the one the rest of METOOLS uses, otherwise the spin
  // density matrix passed to the subsequent strong decay of X sits in
  // a different frame and the spin correlations silently break. Do not
  // replace this by a hand-rolled helicity basis.
  eps.clear();
  // Legacy call signature, confirmed against HADRONS++ VA_P_V/VA_P_A:
  // Polarization_Vector(p, m^2). The mass argument is the GENERATED
  // p^2, not a nominal pole mass, so the basis is exactly orthonormal
  // for the momentum actually produced.
  Polarization_Vector pol(p,p.Abs2());
  for (size_t h=0;h<3;h++) {
    const Vec4C & e = pol[h];
    // Outgoing particle -> conjugated polarization vector.
    eps.push_back(Vec4C(std::conj(e[0]),std::conj(e[1]),
                        std::conj(e[2]),std::conj(e[3])));
  }
}

void VA_P_X_Base::
EffectiveTensorPolarizations(const ATOOLS::Vec4D & pX,
                             const ATOOLS::Vec4D & pD,
                             std::vector<ATOOLS::Vec4C> & epsT) const {
  // eps_T^mu = eps*^{mu nu} p_nu / M, Eq.(effectiveT).
  //
  // This used to assemble the rank-2 tensor by hand from the spin-1
  // basis with 1 (x) 1 -> 2 Clebsch-Gordan coefficients, because the
  // interface of METOOLS::Polarization_Tensor had not been confirmed.
  // The legacy HADRONS++ VA_P_T settles it: the tensor is obtained as
  // Polarization_Tensor pol(p,m^2); CMatrix eps = pol[h]; and the
  // contraction with a four-vector is eps.Conjugate()*p0. Using that
  // directly is not merely shorter - it guarantees the SAME helicity
  // labelling and phase convention as everything else in METOOLS, so
  // the spin density matrix handed to the strong decay of the tensor
  // is expressed in the basis the decay ME expects. A hand-built
  // basis could have differed by a permutation of the five states
  // without any error being raised.
  //
  // The explicit 1/M is the note's own normalisation of eps_T and is
  // part of the definition of the T_i form factors; the legacy current
  // absorbs it into its own h, k, b_+, b_- instead, which is why it
  // does not appear there.
  epsT.clear();
  Polarization_Tensor pol(pX,pX.Abs2());
  for (size_t h=0;h<5;h++) {
    CMatrix eps = pol[h];
    epsT.push_back((eps.Conjugate()*pD)/m_M);
  }
}

//////////////////////////////////////////////////////////////////////////////
// Shared Lorentz kernels - see VA_P_X_Base.H
//////////////////////////////////////////////////////////////////////////////

Vec4C VA_P_X_Base::ScalarKernel(const Vec4D & p,const Vec4D & Q,
                                FF_P_X * fp,FF_P_X * f0) const {
  const Vec4D  q  = p-Q, P = p+Q;
  const double q2 = q.Abs2(), s = Q.Abs2();
  const double mX = sqrt(Max(0.,s));
  // Delta = p^2 - Q^2 always from the momenta: it is the identity
  // behind q_mu V^mu = f_0 Delta.
  const double Delta = p.Abs2()-s;
  const double fpq = fp->Value(q2,mX);
  std::function<double(double)> N =
    [&](double x)->double { return Delta*(f0->Value(x,mX)-fp->Value(x,mX)); };
  const double fm = RegularisedOverQ2(N,q2);
  return Complex(fpq,0.)*Vec4C(P)+Complex(fm,0.)*Vec4C(q);
}

Vec4C VA_P_X_Base::VectorKernel(const Vec4D & p,const Vec4D & Q,
                                const Vec4C & pol,
                                FF_P_X * V,FF_P_X * A0,
                                FF_P_X * A1,FF_P_X * A2,
                                const double & csign) const {
  const Vec4D  q  = p-Q, P = p+Q;
  const double q2 = q.Abs2(), s = Q.Abs2();
  const double m  = sqrt(Max(0.,s)), M = sqrt(Max(0.,p.Abs2()));
  const double Mpm = M+m, Mmm = M-m;
  if (m<1.e-6 || Mpm<1.e-6) return Vec4C(0.,0.,0.,0.);

  // A_0(0) from Eq.(A0constraint) with the GENERATED m, so that
  // G(0)=0 stays exact however the hadronic mass runs.
  const double A1_0 = A1->Value(0.,m), A2_0 = A2->Value(0.,m);
  const double A0_0 = (Mpm*A1_0-Mmm*A2_0)/(2.*m);
  const double A0s0 = A0->Value(0.,m);
  std::function<double(double)> A0f =
    [&](double x)->double {
      return (std::abs(A0s0)>1.e-300 ? A0_0*A0->Value(x,m)/A0s0 : A0_0); };
  std::function<double(double)> G =
    [&](double x)->double {
      return 2.*m*A0f(x) - Mpm*A1->Value(x,m) + Mmm*A2->Value(x,m); };
  const double GoverQ2 = RegularisedOverQ2(G,q2);

  const Complex poldotq = CDot(pol,q);
  const Vec4C Vmu = Complex(0.,2.*V->Value(q2,m)/Mpm)*csign*cross(pol,p,Q);
  const Vec4C Amu =
      poldotq*Complex(GoverQ2,0.)*Vec4C(q)
    + Complex(Mpm*A1->Value(q2,m),0.)*pol
    - poldotq*Complex(A2->Value(q2,m)/Mpm,0.)*Vec4C(P);
  return Vmu-Amu;
}

Vec4C VA_P_X_Base::AxialKernel(const Vec4D & p,const Vec4D & Q,
                               const Vec4C & pol,
                               FF_P_X * A,FF_P_X * V0,
                               FF_P_X * V1,FF_P_X * V2,
                               const double & csign,
                               const double & phase) const {
  const Vec4D  q  = p-Q, P = p+Q;
  const double q2 = q.Abs2(), s = Q.Abs2();
  const double m  = sqrt(Max(0.,s)), M = sqrt(Max(0.,p.Abs2()));
  const double Mpm = M+m, Mmm = M-m;
  // M-m sits in a denominator here, unlike D->V.
  if (m<1.e-6 || Mmm<1.e-4) return Vec4C(0.,0.,0.,0.);

  const double V1_0 = V1->Value(0.,m), V2_0 = V2->Value(0.,m);
  const double V0_0 = (Mmm*V1_0-Mpm*V2_0)/(2.*m);      // Eq.(V0constraint)
  const double V0s0 = V0->Value(0.,m);
  std::function<double(double)> V0f =
    [&](double x)->double {
      return (std::abs(V0s0)>1.e-300 ? V0_0*V0->Value(x,m)/V0s0 : V0_0); };
  std::function<double(double)> G =
    [&](double x)->double {
      return 2.*m*V0f(x) - Mmm*V1->Value(x,m) + Mpm*V2->Value(x,m); };
  const double GoverQ2 = RegularisedOverQ2(G,q2);

  const Complex poldotq = CDot(pol,q);
  const Complex ph(cos(phase),sin(phase));
  const Vec4C Vmu = ph*(
      poldotq*Complex(GoverQ2,0.)*Vec4C(q)
    + Complex(Mmm*V1->Value(q2,m),0.)*pol
    - poldotq*Complex(V2->Value(q2,m)/Mmm,0.)*Vec4C(P) );
  const Vec4C Amu = Complex(-2.*A->Value(q2,m)/Mmm,0.)*csign*cross(pol,p,Q);
  return Vmu-Amu;
}

void VA_P_X_Base::PrintCommonParameters() const {
  msg_Tracking()<<"### "<<m_name<<": "<<m_flavs[p_i[0]]<<" -> "
                <<m_flavs[p_i[1]]<<" l nu\n"
                <<"###   M = "<<m_M<<" GeV, m_X(nominal) = "<<m_mX<<" GeV\n"
                <<"###   V_cq = "<<m_Vcq
                <<(m_ckmdressed?" (CKM-dressed F(0), V_cq set to 1)":"")
                <<", isospin = "<<m_isospin
                <<(m_isKSKL?", extra 1/sqrt(2) for K_S/K_L":"")
                <<"\n###   total prefactor = "<<m_norm
                <<", eps sign = "<<m_epssign
                <<", dynamic m_X = "<<(m_dynamicmass?"yes":"no")<<"\n";
}
