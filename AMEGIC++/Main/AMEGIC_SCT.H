#ifndef AMEGIC_Main_AMEGIC_SCT_H
#define AMEGIC_Main_AMEGIC_SCT_H

#include "ATOOLS/Phys/Spin_Correlation_Tensor.H"
#include "ATOOLS/Math/MyComplex.H"
#include "AMEGIC++/Amplitude/CFColor.H"
#include "AMEGIC++/Amplitude/Amplitude_Base.H"

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
    // Stores the numbers of the associated amplitudes.
    std::vector<int> m_A1, m_A2;

    // Methods to convert AMEGIC convention for spins <-> HADRONS convention.
    int HADRONS_to_AMEGIC(size_t pol, size_t maxPol);
    size_t AMEGIC_to_HADRONS(int pol, size_t maxPol);

  public:
    AMEGIC_SCT(std::vector<int> A1, std::vector<int> A2,
	       Helicity* Hel, 
	       std::vector<int>* particles, size_t pPos =0);

    Spin_Correlation_Tensor* CreateSCT(std::vector<Amplitude_Base*> *graphs,CFColor *col,
				       Helicity* hel);
    ~AMEGIC_SCT();

  }; // end of class AMEGIC_SCT

  /*!
    \class AMEGIC_SCT
    \brief A pre-structure for a Spin_Correlation_Tensor
    
    AMEGIC_SCT is the class which is responsible for transporting spin informations for a
    hard event from AMEGIC to the overall framework SHERPA. The class has two main purposes:
    
    1) Creating the complete tree-structure for spin correlations for each event results in
    a huge performance hit. The class AMEGIC_SCT has been developed to reduce this performance
    loss. The idea is to create the structure of the tree only once at start of the program
    and only to fill in the data at the end-nodes for each event.

    2) AMEGIC and HADRONS use different conventions for labeling spin orientations. AMEGIC_SCT
    takes AMEGIC-system variables as input but creates a ATOOLS::Spin_Correlation_Tensor which
    is compatible with the HADRONS-convention.
  */  
  /*!
    \var std::vector<int> AMEGIC_SCT::m_A1
    Endnode-only variable. The two fields m_A1 and m_A2 store the numbers of the amplitudes needed 
    to calculate the value for the endnode.
    Upon actual creation of a Spin_Correlation_Tensor from the AMEGIC_SCT via AMEGIC_SCT::CreateSCT,
    the values will be calculated.
  */
  /*!
    \var std::vector<int> AMEGIC_SCT::m_A2
    See AMEGIC_SCT::m_A1.
  */
  /*!
    \fn int AMEGIC_SCT::HADRONS_TO_AMEGIC(size_t pol, size_t maxPol)
    The modules HADRONS and AMEGIC use different conventions for labeling spin orientations. This
    method converts a spin orientation from the HADRONS convention to the AMEGIC convention.
    
    Note that the HADRONS convention is the one used in the class Spin_Correlation_Tensor.
  */
  /*!
    \fn size_t AMEGIC_SCT::AMEGIC_TO_HARONS(int pol, size_t maxPol)
    The modules HADRONS and AMEGIC use different conventions for labeling spin orientations. This
    method converts a spin orientation from the AMEGIC convention to the HADRONS convention.  
    
    Note that the HADRONS convention is the one used in the class Spin_Correlation_Tensor.
  */
  /*!
    \fn AMEGIC_SCT::AMEGIC_SCT(std::vector<int> A1, std::vector<int> A2,Helicity* Hel,std::vector<int>* particles, size_t pPos =0);

    Constructs AMEGIC_SCT, the pre-structure for the Spin_Correlation_Tensor. The constructor is
    given a list of particle numbers for which spin informations shall be constructed (parameter
    "particle") and the list of amplitude numbers which contribute to the process ("A1" and "A2" 
    these will simply be the list of all amplitudes.

    The contructor recursively builds up the pre-structure, then. In each step, the entries of the
    two amplitude lists are divided into different sets, each corresponding to a certain spin
    orientation of the current particle "particles[pPos]".

    Then, the follow-up nodes belonging to the different certain spin orientation are created. Each 
    follow-up node is given the corresponding lists; the counter pPos is increased by one.

    This recursion is done until pPos hits the end of the list particles. The two amplitude lists 
    arriving at an endnode are stored there so that the value of the endnode can be calculated
    from them later (on calling AMEGIC_SCT::CreateSCT).
  */
  /*!
    \fn Spin_Correlation_Tensor* AMEGIC_SCT::CreateSCT(std::vector<Amplitude_Base*> *graphs,CFColor *col, Helicity* hel)
    
    This methods extracts spin informations about the currently selected event (the phasespace point
    currently active in Basic_Sfuncs), creates a Spin_Correlation_Tensor and stores them there.

    The method recursively parses the AMEGIC_SCT-tree and copies its structure. At endnodes, the
    value is calculated from the list of amplitude numbers AMEGIC_SCT::m_A1 and AMEGIC_SCT::m_A2
    using the color-matrix color, the diagrams graphs and polarization informations obtained from
    hel.
  */
} // end of namespace AMEGIC
#endif
