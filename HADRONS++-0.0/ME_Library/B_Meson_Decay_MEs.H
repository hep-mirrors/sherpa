#include "HD_ME_Base.H"
#include "Flavour.H"
#include "Vector.H"
#include "MathTools.H"
#include "Histogram.H"
#include "Message.H"
#include "Tools.H"

namespace HADRONS {

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  class Semileptonic_B_Meson : public HD_ME_Base {
    int m_had, m_nu, m_lep;
    double m_hadmass, m_Bmass;
    Complex Formfactor_plus(const ATOOLS::Vec4D *_p);
    Complex Formfactor_minus(const ATOOLS::Vec4D *_p);
  public:
    Semileptonic_B_Meson( int nout, ATOOLS::Flavour *fl );
    void   operator()( 
        const ATOOLS::Vec4D      * _p, 
        std::vector<Complex>     * _ampls_tensor, 
        std::vector<std::pair<int,int> > * _indices,
        int                        k0_n );
    void SetModelParameters( GeneralModel _md );
  };
  
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

};
