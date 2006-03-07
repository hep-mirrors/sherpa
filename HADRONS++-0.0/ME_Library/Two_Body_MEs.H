#ifndef Two_Body_MEs_H
#define Two_Body_MEs_H

#include "HD_ME_Base.H"

namespace HADRONS {
  class P_2Gamma : public HD_ME_Base {
  public:
    P_2Gamma(int,ATOOLS::Flavour *);
    void   operator()( 
        const ATOOLS::Vec4D  * _p, 
        std::vector<Complex> * _ampls_tensor, 
        std::vector<std::pair<int,int> > * _indices,
        int                    k0_n );
    double operator()(const ATOOLS::Vec4D *);
  }; // end of class

  class F_VF : public HD_ME_Base {
    int m_fermion, m_boson;
    double m_global;                        // GF/sqrt(2)
    Complex m_cR, m_cL;                     // couplings for Z function
    double m_MB, m_GB, m_MB2;               // mass, width of W
    public:
      F_VF( int nout, ATOOLS::Flavour *fl );
      void   operator()( 
          const ATOOLS::Vec4D  * _p, 
          std::vector<Complex> * _ampls_tensor, 
          std::vector<std::pair<int,int> > * _indices,
          int k0_n );
      void   SetModelParameters( GeneralModel _md );
  };
   
  class V_FF : public HD_ME_Base {
    int m_ferm1, m_ferm2;
    double m_global;                        // GF/sqrt(2)
    Complex m_cR, m_cL;                     // couplings for Z function
    double m_MB, m_GB, m_MB2;               // mass, width of W
    public:
      V_FF( int nout, ATOOLS::Flavour *fl );
      void   operator()( 
          const ATOOLS::Vec4D  * _p, 
          std::vector<Complex> * _ampls_tensor, 
          std::vector<std::pair<int,int> > * _indices,
          int k0_n );
      void   SetModelParameters( GeneralModel _md );
  };
  /*!
	\file Two_Body_MEs.H
	\brief Declares the ME classes for two-body decays

	This code can be found in the directory \c ME_Library.

    Implemented classes are
    - HADRONS::P_2Gamma
    - HADRONS::V_FF
    - HADRONS::F_VF
    .
  */

  /*!
	\class P_2Gamma 
	\brief Tool to calculate ME for \f$Pseudoscalar \to 2\gamma\f$ decay.

	This class is subclass of HADRONS::HD_ME_Base.
  */

  /*!
	\fn P_2Gamma::operator()(const ATOOLS::Vec4D *)
	\brief Returns amplitude squared.

	Up to now it does nothing else than to return \f$\frac{M}{4}\f$.
  */	

  /*!
	\class F_VF
	\brief Tool to calculate ME for Fermion -> Vector Fermion decay

	This class is subclass of HADRONS::HD_ME_Base.
  */

  /*!
	\class V_FF
	\brief Tool to calculate ME for Vector -> Fermion Fermion decay

	This class is subclass of HADRONS::HD_ME_Base.
  */

	
} // end of namespace


#endif
