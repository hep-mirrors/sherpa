#ifndef Beam_Remnant_Handler_H
#define Beam_Remnant_Handler_H

#include "Beam_Spectra_Handler.H"
#include "Primordial_KPerp.H"
#include "ISR_Handler.H"
#include "Blob_List.H"
#include "Primordial_KPerp.H"
#include "Remnant_Base.H"

namespace SHERPA {

  class Beam_Remnant_Handler {
  private:

    PDF::ISR_Handler            * p_isr;
    BEAM::Beam_Spectra_Handler  * p_beam;

    SHERPA::Remnant_Base        * p_beampart[2];
    Primordial_KPerp            * p_kperp;

    ATOOLS::Blob                * p_beamblob[2];

    std::string  m_path, m_file;
    bool         m_fill;

  public :

    // constructor
    Beam_Remnant_Handler(std::string,std::string,
			 PDF::ISR_Handler *,BEAM::Beam_Spectra_Handler *,double);

    // destructor
    ~Beam_Remnant_Handler();
    
    // member functions
    bool FillBunchBlobs(ATOOLS::Blob_List *,ATOOLS::Particle_List * = NULL);
    bool FillBeamBlobs(ATOOLS::Blob_List *,ATOOLS::Particle_List * = NULL);

  };// end of class Beam_Remnant_Handler

  /*! 
    \class Beam_Remnant_Handler
    \brief This class fills in the beam remnants after the hard scatters are defined.
  */

  /*!
    \var bool Beam_Remnant_Handler::m_fill
    This flag equals 1, if any of the two bunch-particles is a hadron and therefore has to be filled. 
  */

  /*!
    \var double Beam_Remnant_Handler::m_q2min
    The minimal \f$Q^2\f$ covered by the PDF's. At the moment it is set hard-wired to
    \f$Q^2 = 1 \mbox{\rm GeV}^2\f$, but this will be changed to match the actual PDF.
  */

  /*!
    \var std::string Beam_Remnant_Handler::m_path
    The directory where the data file can be found.
  */

  /*!
    \var std::string Beam_Remnant_Handler::m_file
    The name of the data file.
  */

  /*!
    \var PDF::ISR_Handler * Beam_Remnant_Handler::p_isr
    A pointer to the ISR_Handler to have access to the PDFs.
  */

  /*!
    \var BEAM::Beam_Spectra_Handler * Beam_Remnant_Handler::p_beam
    A pointer to the beam handler, mainly to fill also the blobs for the transition of beams to
    bunch-particles.
  */

  /*!
    \var ATOOLS::Flavour ** Beam_Remnant_Handler::p_constituents
    A pointer to the incoming hadrons' constituents. It is filled by the method 
    Constituents(ATOOLS::Flavour,ATOOLS::Flavour *&) that analyses the kf-codes of the hadrons to
    extract the two or three constituent quarks.
  */

  /*!
    \var int * Beam_Remnant_Handler::p_numberofconstituents
    The respective number of constituents for every incoming hadron, either 2 or 3 for mesons and
    baryons.
  */

  /*!
    \fn int Beam_Remnant_Handler::Constituents(ATOOLS::Flavour,ATOOLS::Flavour *&) 
    Fills in the pointers containing the constituent quakrs of the incoming hadrons.
    By analysing their kf-code. It should be stressed here that we still have to thing
    about some of the neutral mesons 9like for instance the \f$\pi^0\f$ which has
    a \f$u\bar u\f$- and a \f$d\bar d\f$-component that is not explicit in the kf-code.
  */

  /*!
    \fn bool Beam_Remnant_Handler::FillHadron(ATOOLS::Blob *,int,ATOOLS::Parton_List * = NULL)
    This is the central method to fill in the hadron remnant, after one hard scatter has been taken out.
    Some cases are to be considered:
    - The hard scatter is a constituent quark:</br>
      Then the colour flow is pretty simple. The remainder is just the other constituent (either a quark
      for incoming mesons or a diquark for baryons), and its colour is a simple triplet or anti-triplet
      depending on whether hte hard scatter constituent is an anti-triplet or a triplet.
    - The hard scatter is a gluon:</br>
      Again, the colour flow is pretty simple. The gluons triplet colour component connects to a
      corresponding anti-triplet (either an anti-quark or a diquark), its anti-triplet colour
      connects to a triplet (quark or anti-diquark). The only tricky point occurs for baryons,
      because then it has to be decided which two constituent (anti-)quarks form the (anti-)diquark.
      This is done in the following way: The single quark or anti-quark is chosen with equal
      probability among the three constituents, for the resulting diquark the eventually emerging
      different spin-states are chosen at random with a ratio of 3:1 for the spin-1 versus the spin-0
      state.
    - The hard scatter is a seaquark:</br>
      For simplicitly, focus on quarks, anti-quarks can be obtained from thre treatment described by
      simple replacements. For the hard scatter seaquark, an anti-quark has to be kicked out of the
      sea, too. It will
  */

  /*!
    \fn double Beam_Remnant_Handler::GetXLund(ATOOLS::Flavour,ATOOLS::Flavour,double)
  */

  /*!
    \fn double Beam_Remnant_Handler::GetXPDF(ATOOLS::Flavour,ATOOLS::Flavour,double,int)
  */

  /*!
    \fn double Beam_Remnant_Handler::PrimordialKPerp(unsigned int,std::vector<ATOOLS::Vec3D)
  */

  /*!
    \fn Beam_Remnant_Handler::Beam_Remnant_Handler(std::string,std::string,PDF::ISR_Handler *,BEAM::Beam_Spectra_Handler *)
  */

  /*!
    \fn Beam_Remnant_Handler::~Beam_Remnant_Handler()
  */

  /*!
    \fn bool Beam_Remnant_Handler::FillBunchBlobs(ATOOLS::Blob_List *,ATOOLS::Parton_List * = NULL)
  */

  /*!
    \fn bool Beam_Remnant_Handler::FillBeamBlobs(ATOOLS::Blob_List *,ATOOLS::Parton_List * = NULL)
  */

}// end of namespace SHERPA

#endif

