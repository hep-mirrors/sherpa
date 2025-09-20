#ifndef EXTRAXS_One2Three_H_to_bbg_NLO_Complete_H
#define EXTRAXS_One2Three_H_to_bbg_NLO_Complete_H

#include "METOOLS/Main/Spin_Structure.H"
#include "ATOOLS/Math/Vector.H"       
#include "ATOOLS/Phys/Flavour.H"           
#include <vector>   


namespace EXTRAXS {
  class H_to_bbg_Real;
  class H_to_bb_Virtual;
}

namespace EXTRAXS {
  
  class H_to_bbg_NLO_Complete : public METOOLS::Spin_Amplitudes {
  private:
    H_to_bbg_Real* p_real;
    H_to_bb_Virtual* p_virtual;
    
  public:
    H_to_bbg_NLO_Complete(const std::vector<ATOOLS::Flavour>& flavs,
                          const ATOOLS::Flavour& prop,
                          size_t non_prop, size_t gluon, size_t propj);
    
    ~H_to_bbg_NLO_Complete();
    
    void Calculate(const ATOOLS::Vec4D_Vector& momenta, bool anti);
    
    // Debugging/Testing
    double GetRealContribution() const;
    double GetVirtualContribution() const;
  };

}

#endif