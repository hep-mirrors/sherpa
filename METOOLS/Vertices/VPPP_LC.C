#include "METOOLS/Explicit/Lorentz_Calculator.H"
#include "METOOLS/Currents/C_Scalar.H"
#include "METOOLS/Currents/C_Vector.H"
#include "METOOLS/Explicit/Vertex.H"
#include "MODEL/Main/Single_Vertex.H"
#include "ATOOLS/Phys/Spinor.H"

using namespace ATOOLS;

namespace METOOLS {

  /*!
    gamma* -> pi+ pi- pi0 as a SINGLE contact current,

      A^mu = eps^{mu,al,be,ga} p_{+,al} p_{-,be} p_{0,ga}

    with the whole dynamics - the isoscalar tower in q^2 and the rho tower in
    the two-pion sub-masses - carried by the accompanying form factor (the
    "EE3PI" mode of Meson_FormFactors.C), not by this calculator.

    WHY A CONTACT VERTEX.  The alternative, and what this model did before, is
    to factorise into gamma -> rho pi (AVVP) followed by rho -> pi pi (SSV) and
    let Comix build the chain.  That is a different object, and three separate
    problems come with it:
      - the relative weight of the three charge channels is then whatever the
        epsilon contraction and the SSV momentum difference happen to produce.
        Measured, that is (+2,-1,+1) relative to eps(e,p_+,p_-,p_0), whereas an
        I = 0 three-pion state requires an ALL-PLUS sum: |A|^2 has to be totally
        symmetric under permuting the pions, because the isospin state and the
        momentum epsilon are both totally antisymmetric;
      - the overall normalisation is split across two couplings and an internal
        propagator, so the Wess-Zumino-Witten anomaly - which PREDICTS this
        amplitude at threshold and leaves nothing to fit - cannot be imposed
        directly and ends up as a tuned number;
      - the sum-of-products structure of hep-ph/0512180 eq. (10), where an
        excited rho pairs with one specific isoscalar, is inexpressible as a
        product of two independent sums.
    Writing the current in one place makes all three explicit and checkable.

    Follows the construction of tau-decays-in-sherpa (METOOLS-hadron-currents),
    METOOLS/HadronCurrents/V_0_Isoscalar3Pi.C and FF_0_Isoscalar3Pi.C, which
    implement Hoferichter, Hoid and Kubis, arXiv:1907.01556 - but written
    natively for this model's Comix vertex interface rather than imported.

    LEG ORDER is fixed by the vertex declaration in
    LowEnergy_Model::InitThreePionVertices: { gamma, pi+, pi-, pi0 }, i.e.
    leg 0 = photon, 1 = pi+, 2 = pi-, 3 = pi0.
  */
  template <typename SType>
  class VPPP_Calculator: public Lorentz_Calculator {
  public:

    typedef std::complex<SType> SComplex;
    typedef CVec4<SType>   CVec4Type;
    typedef CScalar<SType> CScalarType;

    VPPP_Calculator(const Vertex_Key &key):
      Lorentz_Calculator(key) {}

    std::string Label() const { return "VPPP"; }

  private:

    static inline unsigned int I0() { return 0; }
    static inline unsigned int I1() { return Spinor<SType>::R1(); }
    static inline unsigned int I2() { return Spinor<SType>::R2(); }
    static inline unsigned int I3() { return Spinor<SType>::R3(); }

    template <class A, class B, class C>
    static SComplex Minor(const A &r0,const B &r1,const C &r2,
			  const unsigned int &c0,const unsigned int &c1,
			  const unsigned int &c2) {
      return ( r0[c0]*(r1[c1]*r2[c2]-r1[c2]*r2[c1])
	      -r0[c1]*(r1[c0]*r2[c2]-r1[c2]*r2[c0])
	      +r0[c2]*(r1[c0]*r2[c1]-r1[c1]*r2[c0]) );
    }

    // J^mu = eps^{mu,al,be,ga} a_al b_be c_ga.  Same index and sign
    // convention as AVVP_LC.C::EpsVec, which this model already uses for the
    // anomalous three-leg vertices.
    template <class A, class B, class C>
    static CVec4Type EpsVec(const A &a,const B &b,const C &c) {
      CVec4Type j;
      j[I0()] =  Minor(a,b,c,I1(),I2(),I3());
      j[I1()] =  Minor(a,b,c,I0(),I2(),I3());
      j[I2()] = -Minor(a,b,c,I0(),I1(),I3());
      j[I3()] =  Minor(a,b,c,I0(),I1(),I2());
      return j;
    }

    // eps_{mu,al,be,ga} v^mu a^al b^be c^ga
    template <class V, class A, class B, class C>
    static SComplex Eps4(const V &v,const A &a,const B &b,const C &c) {
      return ( v[I0()]*Minor(a,b,c,I1(),I2(),I3())
	      -v[I1()]*Minor(a,b,c,I0(),I2(),I3())
	      +v[I2()]*Minor(a,b,c,I0(),I1(),I3())
	      -v[I3()]*Minor(a,b,c,I0(),I1(),I2()) );
    }

  public:

    CObject *Evaluate(const CObject_Vector &jj)
    {
      // jj[i] and p_v->J(i) belong to declared leg id[i]; momenta outgoing
      // from the vertex. The output current's flavour is the conjugate of
      // its leg, so legs are identified by id, not by flavour.
      const std::vector<size_t> &id(p_v->V()->id);
      const size_t nl(4), out(id.back());
      if (out>=nl || jj.size()!=nl-1) return NULL;
      size_t slot[4];
      Vec4D p[4];
      for (size_t i(0);i<nl-1;++i) {
	slot[id[i]]=i;
	p[id[i]]=p_v->J(i)->P();
      }
      p[out]=-p_v->JC()->P();
      if (out==0) {
	SComplex sc(1.,0.);
	for (size_t l(1);l<nl;++l) sc*=(*jj[slot[l]]->Get<CScalarType>())[0];
	CVec4Type *j(CVec4Type::New(EpsVec(p[1],p[2],p[3])*sc));
	j->SetS(jj[0]->S()|jj[1]->S()|jj[2]->S());
	return j;
      }
      const CVec4Type &v(*jj[slot[0]]->Get<CVec4Type>());
      SComplex sc(1.,0.);
      for (size_t l(1);l<nl;++l)
	if (l!=out) sc*=(*jj[slot[l]]->Get<CScalarType>())[0];
      CScalarType *j(CScalarType::New(Eps4(v,p[1],p[2],p[3])*sc));
      j->SetS(jj[0]->S()|jj[1]->S()|jj[2]->S());
      return j;
    }

  };// end of class VPPP_Calculator

  template class VPPP_Calculator<double>;

}// end of namespace METOOLS

using namespace METOOLS;

DECLARE_GETTER(VPPP_Calculator<double>,"DVPPP",
	       Lorentz_Calculator,Vertex_Key);
Lorentz_Calculator *ATOOLS::Getter
<Lorentz_Calculator,Vertex_Key,VPPP_Calculator<double> >::
operator()(const Vertex_Key &key) const
{ return new VPPP_Calculator<double>(key); }

void ATOOLS::Getter<Lorentz_Calculator,Vertex_Key,
		    VPPP_Calculator<double> >::
PrintInfo(std::ostream &str,const size_t width) const
{ str<<"anomalous gamma* -> pi pi pi contact vertex"; }
