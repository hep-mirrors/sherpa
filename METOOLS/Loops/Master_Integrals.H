#ifndef METOOLS_Loops_Master_Integrals_H
#define METOOLS_Loops_Master_Integrals_H

#include "ATOOLS/Math/MyComplex.H"
#include "METOOLS/Loops/Divergence_Array.H"

// set global renormalisation scale as backup to m_Z^2
#define GLOBAL_RENORMALISATION_SCALE 8315.251344

namespace METOOLS {

  //! dimension D=4-2epsilon
  inline DivArrD D() { return DivArrD(0.,0.,0.,4.,-2.,0.); }
#define D D()
  // Helper functions
  Complex CLog(const double&, const int& = 1);
  Complex CLog(const Complex&, const int& = 1);
  Complex Eta(const Complex&, const int&, const Complex&, const int&); 
  double Theta(const double&);
  Complex K(const double&, const Complex&, const Complex&);
  Complex Dilog(const Complex&, const int& = 1);
  int ReSign(const Complex&);
  int ImSign(const Complex&);
  int Prod2Sign(const Complex&,const int&, const Complex&,const int&);
  int Prod3Sign(const Complex&,const int&, const Complex&,const int&, const Complex&,const int&);
  int DivSign(const Complex&,const int&, const Complex&,const int&);
  int ProdDivSign(const Complex&,const int&, const Complex&,const int&, const Complex&,const int&);
  int DivProdSign(const Complex&,const int&, const Complex&,const int&, const Complex&,const int&);
  void rfun(Complex&, Complex&, const Complex&);
  void solabc(Complex&, Complex&, Complex&, const Complex&, const Complex&,
		       const Complex&);
  /*! Convention on input:
      - \f$ p_i^2 \f$ labels the squared outer momenta
      - \f$ m_i^2 \f$ labels the masses of the internal lines
      - \f$ s_{ij}=(p_i + p_j)^2 \f$
      - \f$ \mu^2 \f$ labels the renormalisation scale
      Convention on output:
      - array with the following entries
        [0] - \f$ \sim \frac{1}{\epsilon_{UV}} \f$
        [1] - \f$ \sim \frac{1}{\epsilon_{IR}} \f$
        [2] - \f$ \sim \frac{1}{\epsilon_{IR}^2} \f$
        [3] - \f$ \sim 1 \f$
        [4] - \f$ \sim \epsilon \f$
        [5] - \f$ \sim \epsilon^2 \f$
      All terms \f$ \sim \epsilon/\epsilon^2 \f$ are currently not calculated.
      This is suitable for 1-loop integrals, but nothing beyond.
   */
  //! A_0(m2;mu2)
  DivArrC
  Master_Tadpole(const Complex&,
                 double);
  //! B_0(s12;m02,m12;mu2)
  DivArrC
  Master_Bubble(const double&,
                const Complex&, const Complex&,
                double);
  //! B_0_prime(s12;m02,m12;mu2)
  DivArrC
  Master_Bubble_Prime(const double&,
		      const Complex&, const Complex&,
		      double);
  //! C_0(p12,p22,s12;m02,m12,m22;mu2)
  DivArrC
  Master_Triangle(const double&,  const double&,  const double&,
                  const Complex&, const Complex&, const Complex&,
                  double);
  //! D_0(p12,p22,p32,p43,s12,s23;m12,m22,m32,m42;mu2)
  DivArrC
  Master_Box(const double&,  const double&,  const double&,  const double&,
             const double&,  const double&,
             const Complex&, const Complex&, const Complex&, const Complex&,
             double);


}

#endif
