#ifndef Simple_String_H
#define Simple_String_H

#include "MI_Base.H"
#include "Environment.H"
#ifdef USING__Sherpa
#include "Remnant_Base.H"
#endif

namespace AMISIC {

  class Reggeon_Trajectory;

  class Simple_String: public MI_Base {
  private:

#ifdef USING__Sherpa
    SHERPA::Remnant_Base  *p_remnants[2];
#endif

    std::vector<Reggeon_Trajectory*> m_reggeons;

    void CleanUp();
    void PrepareTerminate();

    bool FillBlob(ATOOLS::Blob *blob);

  public:

    // constructor
    Simple_String();
    Simple_String(MODEL::Model_Base *_p_model,
		  BEAM::Beam_Spectra_Handler *_p_beam,
		  PDF::ISR_Handler *_p_isr);

    // destructor
    ~Simple_String();

    // member functions
    bool Initialize();
    void Update(const MI_Base *mibase);
    void Reset();

    bool VetoProcess(ATOOLS::Blob *blob);
    bool DiceProcess();

  }; // end of class Simple_String

} // end of namespace AMISIC

#endif
