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
    leg 0 = photon, 1 = pi+, 2 = pi-, 3 = pi0.  Momenta are always taken from
    p_v->J(i)->P(), which is in ASCENDING leg order; jj[] arrives in cyclic
    order after the outgoing leg and is used only for the current values.
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
      // Momenta of the three pions, identified by FLAVOUR rather than by leg
      // index.  p_v->J() holds the INCOMING currents - one fewer than the
      // number of legs - so on this four-leg vertex it has three entries and
      // the outgoing leg is reached through p_v->JC() instead.  Indexing
      // J(0..3) directly reads past the end and segfaults.
      Vec4D pp(0.,0.,0.,0.), pm(0.,0.,0.,0.), p0(0.,0.,0.,0.);
      bool hp(false), hm(false), h0(false);
      for (size_t i(0);i<p_v->J().size();++i) {
	const Flavour & fl(p_v->J(i)->Flav());
	if (!fl.IsHadron()) continue;
	const Vec4D p(p_v->J(i)->P());
	if      (fl.Kfcode()==kf_pi)                          { p0=p; h0=true; }
	else if (fl.Kfcode()==kf_pi_plus && fl.IntCharge()>0) { pp=p; hp=true; }
	else if (fl.Kfcode()==kf_pi_plus && fl.IntCharge()<0) { pm=p; hm=true; }
      }
      if (!(hp&&hm&&h0)) {
	const Flavour & fc(p_v->JC()->Flav());
	const Vec4D pc(p_v->JC()->P());
	if      (!h0 && fc.Kfcode()==kf_pi) { p0=pc; h0=true; }
	else if (!hp && fc.IntCharge()>0)   { pp=pc; hp=true; }
	else if (!hm && fc.IntCharge()<0)   { pm=pc; hm=true; }
      }
      if (!(hp&&hm&&h0)) return NULL;

      // Leg 0 outgoing: solve for the photon.  This is the branch e+e- -> 3pi
      // actually uses.  The three pion currents are scalars, so they multiply
      // in as plain numbers and the vector index is carried by the epsilon.
      if (p_v->V()->id.back()==0) {
	const CScalarType &a(*jj[0]->Get<CScalarType>());
	const CScalarType &b(*jj[1]->Get<CScalarType>());
	const CScalarType &c(*jj[2]->Get<CScalarType>());
	CVec4Type *j(CVec4Type::New(EpsVec(pp,pm,p0)*(a[0]*b[0]*c[0])));
	j->SetS(a.S()|b.S()|c.S());
	return j;
      }
      // A pion is outgoing: contract the photon current into the epsilon and
      // return a scalar.  Which pion is being solved for does not change the
      // structure - the epsilon already carries all three pion momenta from
      // p_v->J(i)->P(), including the outgoing leg's own (supplied there by
      // momentum conservation).
      //
      // jj[] is CYCLIC after the outgoing leg, so for legs 1, 2 and 3 the
      // photon - leg 0 - sits at cyclic position 3-id, and the two remaining
      // pions fill the others.  Spelling the three cases out rather than
      // searching by type: CObject carries no type tag, and guessing the slot
      // is how leg-ordering bugs get in.
      const size_t out(p_v->V()->id.back());
      if (out==1 || out==2 || out==3) {
	const size_t ig(3-out);          // cyclic slot holding the photon
	const CVec4Type &v(*jj[ig]->Get<CVec4Type>());
	SComplex sc(1.,0.);
	for (size_t i(0);i<jj.size();++i)
	  if (i!=ig) sc*=(*jj[i]->Get<CScalarType>())[0];
	CScalarType *j(CScalarType::New(Eps4(v,pp,pm,p0)*sc));
	j->SetS(jj[0]->S()|jj[1]->S()|jj[2]->S());
	return j;
      }
      return NULL;
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
