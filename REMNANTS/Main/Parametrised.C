#include "REMNANTS/Main/Parametrised.H"
#include "PDF/Main/ISR_Handler.H"
#include "BEAM/Main/Beam_Spectra_Handler.H"
#include "ATOOLS/Org/Message.H"

#ifdef PROFILE__all
#define PROFILE__Parametrised_Beam_Remnants
#endif
#ifdef PROFILE__Parametrised_Beam_Remnants
#include "prof.hh" 
#else
#define PROFILE_HERE
#endif

using namespace REMNANTS;
using namespace ATOOLS;

Parametrised::Parametrised(const std::string path,const std::string file,
			   Remnant_Handler * const remnants,
			   BEAM::Beam_Spectra_Handler *const beam):
  p_rhandler(remnants)
{
  p_rhandler->InitializeKinematics(path,file);
}

Parametrised::~Parametrised() {}


  

