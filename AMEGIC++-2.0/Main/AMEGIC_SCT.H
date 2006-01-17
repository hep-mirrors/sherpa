#ifndef AMEGIC_SCT_H
#define AMEGIC_SCT_H

#include "Spin_Correlation_Tensor.H"
#include "MyComplex.H"
#include "CFColor.H"

using namespace ATOOLS;
namespace AMEGIC {
  class AMEGIC_SCT : public Spin_Correlation_Tensor
  {
  private: 
    //! Converts the HADRONS-convention for spins to a pseudo-AMEGIC convention.
    int HADRONS_to_AMEGIC(size_t pol, size_t maxPol);
  public:
    AMEGIC_SCT(std::vector< std::vector<Complex> > *Ampls, 
	       CFColor *ColorMatrix, Helicity *Hel,
	       std::vector<int> comb1 = std::vector<int>(0), 
	       std::vector<int> comb2 = std::vector<int>(0));
    ~AMEGIC_SCT();
  }; // end of class AMEGIC_SCT
} // end of namespace AMEGIC
#endif
