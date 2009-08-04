#ifndef METOOLS_Main_XYZFuncs_H
#define METOOLS_Main_XYZFuncs_H

#include <vector>
#include "ATOOLS/Math/MyComplex.H"
#include "ATOOLS/Math/Vector.H"
#include "ATOOLS/Phys/Flavour.H"

namespace METOOLS {
  class XYZFunc {
  protected:
    int                 m_N;
    int                 m_k0n;
    bool                m_anti;
    ATOOLS::Vec4D     * p_mom;
    ATOOLS::Flavour   * p_flav;
    const int         * p_i;
    std::vector<Complex>     m_eta, m_mu;

    void CalcEtaMu();
    int MToL(int m);

    Complex S( const int s, const int i, const int j );
    Complex S( const int s, const int i, const ATOOLS::Vec4C p, Complex eta );
    Complex S( const int s, const ATOOLS::Vec4C p, Complex eta, const int j );

    Complex Z( const int t1, const int t2, const int t3, const int t4,
               const int hel_comb,
	       const Complex cR1, const Complex cL1,
	       const Complex cR2, const Complex cL2 );
    Complex Y( const int t1, const int t2, const int hel_comb, 
               const Complex cR, const Complex cL );
    Complex X( const int t1, const ATOOLS::Vec4C p2, const int t3,
	       const int hel_comb, const Complex cR, const Complex cL );
    ATOOLS::Vec4C L( const int t1, const int t2,
		      const int hel_comb, const Complex cR, const Complex cL );

  public:
    XYZFunc( int n, const ATOOLS::Vec4D *p, const ATOOLS::Flavour *fl,
             int k0_n, bool anti=false, const int *indices=NULL );
    XYZFunc( const int n, const ATOOLS::Flavour *fl, const int k0_n,
             const int *indices=NULL );
    ~XYZFunc();

    void Prepare( const ATOOLS::Vec4D *p, const bool anti=false );

    Complex Z(
	      const int t1, const int l1,
	      const int t2, const int l2,
	      const int t3, const int l3,
	      const int t4, const int l4,
	      const Complex cR1, const Complex cL1,
	      const Complex cR2, const Complex cL2 );
    Complex Y(
	      const int t1, const int l1,
	      const int t2, const int l2,
	      const Complex cR, const Complex cL );
    Complex X(
	      const int t1, const int l1,
	      const ATOOLS::Vec4C p2,
	      const int t3, const int l3,
	      const Complex cR, const Complex cL );
    Complex G(
	      const int t1, const int l1,
	      const ATOOLS::Vec4C p2,
	      const int t3, const int l3 );
    ATOOLS::Vec4C L(
	      const int t1, const int l1,
	      const int t2, const int l2,
	      const Complex cR, const Complex cL );

    ATOOLS::Vec4C Y31(const int t1, const int l1,
                       const int t2, const int l2,
                       Complex cR, Complex cL );

    ATOOLS::Vec4C Y13(const int t1, const int l1,
                       const int t2, const int l2,
                       Complex cR, Complex cL );

    ATOOLS::Vec4C X31(const int t1, const int l1,
                       const ATOOLS::Vec4C p2,
                       const int t3, const int l3,
                       Complex cR, Complex cL );

    ATOOLS::Vec4C X13(const int t1, const int l1,
                       const ATOOLS::Vec4C p2,
                       const int t3, const int l3,
                       Complex cR, Complex cL );

    ATOOLS::Vec4C L31(const int t1, const int l1,
                       const ATOOLS::Vec4C p,
                       const int t2, const int l2,
                       Complex cR, Complex cL );

    ATOOLS::Vec4C L13(const int t1, const int l1,
                       const ATOOLS::Vec4C p,
                       const int t2, const int l2,
                       Complex cR, Complex cL );
  };


  /*!
    \file XYZFuncs.H
    \brief Declares the class ATOOLS::XYZFunc

    This file can be found in the directory \c Helicities.
  */

  /*!
    \class XYZFunc
    \brief Tools to calculate X, Y, and Z functions

    This class contains everything that is necessery to calculate the value
    of an X, Y, or Z function, which can be used to calculate a decay matrix
    element.
    <b>Note!</b> Helicity combinations are coded in like a binary number such that
    - \f$++++ \equiv 0\f$
    - \f$+++- \equiv 1\f$
    - \f$++-+ \equiv 2\f$
    - ...
    - \f$---- \equiv 15\f$
    .
  */

  /*!
    \fn XYZFunc::XYZFunc( int n, const Vec4D *p, const Flavour *fl )
    \brief Constructor to set up the calculation

    This method saves the momenta and flavours so that they are available throughout the class. Furthermore, it
    calls XYZFunc::CalcEtaMu which is tasked with the calculation of \f$\eta_i\f$ and \f$\mu_i\f$ for each
    particle.
  */
  /*!
    \fn XYZFunc::CalcEtaMu()
    \brief Calculates \f$\eta_i\f$ and \f$\mu_i\f$ for each particle

    Eta's and Mu's are calcuted and stored in private attributes. For Eta's \f$k_0\f$ is set to
    \f$k_0=(1,0,\frac{1}{\sqrt{2}},\frac{1}{\sqrt{2}})\f$ if XYZFunc::m_k0n is 1 (default).
  */

  /*!
    \var XYZFunc::p_mom
    Pointer on involved 4-momenta.
  */

  /*!
    \var XYZFunc::p_flav
    Pointer on involved flavours.
  */

  /*!
    \var XYZFunc::N
    Number of involved particles (in- and out-particles)
  */
  /*!
    \var XYZFunc::m_k0n
    It is either 0, 1 or 2. The used vector \f$k_0\f$ is then either 
    \f$k_0=(1,\frac{1}{\sqrt{2}},0,\frac{1}{\sqrt{2}})\f$,
    \f$k_0=(1,0,\frac{1}{\sqrt{2}},\frac{1}{\sqrt{2}})\f$ or
    \f$k_0=(1,\frac{1}{\sqrt{2}},\frac{1}{\sqrt{2}},0)\f$, respectively.
  */
}

#endif

