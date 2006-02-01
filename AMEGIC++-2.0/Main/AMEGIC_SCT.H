#ifndef AMEGIC_SCT_H
#define AMEGIC_SCT_H

#include "Spin_Correlation_Tensor.H"
#include "MyComplex.H"
#include "CFColor.H"
#include "Amplitude_Base.H"

using namespace ATOOLS;
using namespace std;

namespace AMEGIC {

  //! This is an evil trick...
  class SCT_DUMMY : public Spin_Correlation_Tensor
  {
    friend class AMEGIC_SCT;
  };

  class AMEGIC_SCT : public Spin_Correlation_Tensor
  {
  private: 
    //! Stores the numbers of the associated amplitudes.
    std::vector<int> m_A1, m_A2;
    //! Converts the HADRONS-convention for spins to a pseudo-AMEGIC convention
    int HADRONS_to_AMEGIC(size_t pol, size_t maxPol);
    //! Convers the AMEGIC-convention for spins to the HADRONS-convention.
    size_t AMEGIC_to_HADRONS(int pol, size_t maxPol);
  public:
    /*    AMEGIC_SCT(std::vector< std::vector<Complex> > *Ampls, 
	       CFColor *ColorMatrix, Helicity *Hel,
	       std::vector<int> comb1 = std::vector<int>(0), 
	       std::vector<int> comb2 = std::vector<int>(0)); */
    AMEGIC_SCT(std::vector<int> A1, std::vector<int> A2,
	       Helicity* Hel, 
	       std::vector<int>* particles, size_t pPos =0);

    Spin_Correlation_Tensor* CreateSCT(std::vector<Amplitude_Base*> *graphs,CFColor *col,
				       Helicity* hel);
    ~AMEGIC_SCT();

  }; // end of class AMEGIC_SCT


} // end of namespace AMEGIC

#endif
