#ifndef Beam_Remnant_Handler_h
#define Beam_Remnant_Handler_h

#include <set> 
#include <iterator> 
#include <iostream> 

#include "Beam_Spectra_Handler.H"
#include "ISR_Handler.H"
#include "Parton_List.H"
#include "Parton.H"
#include "Blob_List.H"
#include "Blob.H"

namespace SHERPA {
  class Beam_Remnant_Handler {
    bool                          m_fill;
    double                        m_q2min;
    std::string                   m_dir,m_file;
    ISR::ISR_Handler            * p_isr;
    BEAM::Beam_Spectra_Handler  * p_beam;
    APHYTOOLS::Flavour         ** p_constituents;
    int                         * p_numberofconstituents;

    int    Constituents(APHYTOOLS::Flavour,APHYTOOLS::Flavour *&); 
    bool   FillHadron(APHYTOOLS::Blob *,int,APHYTOOLS::Parton_List * = NULL);
    double GetX_Lund(APHYTOOLS::Flavour,APHYTOOLS::Flavour,double);
    double GetX_PDF(APHYTOOLS::Flavour,APHYTOOLS::Flavour,double,int);
    double GetPT(APHYTOOLS::Flavour,APHYTOOLS::Flavour,double);
  public :
    Beam_Remnant_Handler(string,string,ISR::ISR_Handler *,BEAM::Beam_Spectra_Handler *);
    ~Beam_Remnant_Handler();

    bool FillBunchBlobs(APHYTOOLS::Blob_List *,APHYTOOLS::Parton_List * = NULL);
    bool FillBeamBlobs(APHYTOOLS::Blob_List *,APHYTOOLS::Parton_List * = NULL);
  };
}
#endif

