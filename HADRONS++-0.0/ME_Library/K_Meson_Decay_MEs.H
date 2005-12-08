#include "HD_ME_Base.H"
#include "Flavour.H"
#include "Vector.H"
#include "MathTools.H"
#include "Histogram.H"
#include "Message.H"
#include "Tools.H"

namespace HADRONS {

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  class K_Meson_Lepton : public HD_ME_Base {
    int m_lep, m_nulep;
    double m_GF2;
    Complex m_cR, m_cL;         // couplings for X function
   public:
     K_Meson_Lepton( int nout, ATOOLS::Flavour *fl );
     double operator()(const ATOOLS::Vec4D *_p );
     double Using_Hels( const ATOOLS::Vec4D *_p );
     void SetModelParameters( GeneralModel _md );
  };

  class K_Meson_SemiLeptonic : public HD_ME_Base {
    int m_lep, m_nulep, m_pi;
    double m_GF2;
    Complex m_cR_K, m_cL_K;         // couplings for X-function with mom. of K
    Complex m_cR_pi, m_cL_pi;         // couplings for X-function with mom. of pi
    public:
      K_Meson_SemiLeptonic( int nout, ATOOLS::Flavour *fl );
      double operator()(const ATOOLS::Vec4D *_p );
      double Using_Hels( const ATOOLS::Vec4D *_p );
      void SetModelParameters( GeneralModel _md );
  };
}
