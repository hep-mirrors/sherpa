#ifndef Beam_Remnant_Handler_h
#define Beam_Remnant_Handler_h

#include <set> 
#include <iterator> 
#include <iostream> 

#include "Beam_Spectra_Handler.H"
#include "ISR_Handler.H"
#include "Particle_List.H"
#include "Particle.H"
#include "Blob_List.H"
#include "Blob.H"

namespace SHERPA {
  class Beam_Remnant_Handler {
    bool                          m_fill;
    double                        m_q2min;
    std::string                   m_dir,m_file;
    PDF::ISR_Handler            * p_isr;
    BEAM::Beam_Spectra_Handler  * p_beam;
    ATOOLS::Flavour         ** p_constituents;
    int                         * p_numberofconstituents;

    int    Constituents(ATOOLS::Flavour,ATOOLS::Flavour *&); 
    bool   FillHadron(ATOOLS::Blob *,int,ATOOLS::Particle_List * = NULL);
    double GetX_Lund(ATOOLS::Flavour,ATOOLS::Flavour,double);
    double GetX_PDF(ATOOLS::Flavour,ATOOLS::Flavour,double,int);
    double GetPT(ATOOLS::Flavour,ATOOLS::Flavour,double);
  public :
    Beam_Remnant_Handler(std::string,std::string,PDF::ISR_Handler *,BEAM::Beam_Spectra_Handler *);
    ~Beam_Remnant_Handler();

    bool FillBunchBlobs(ATOOLS::Blob_List *,ATOOLS::Particle_List * = NULL);
    bool FillBeamBlobs(ATOOLS::Blob_List *,ATOOLS::Particle_List * = NULL);
  };
}
#endif

