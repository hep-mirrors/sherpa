#ifndef Apacic_H
#define Apacic_H

#include "Hard_Processes.H"
#include "Hard_Interface.H"

#include "Parton_List.H"
#include "Blob_List.H"

namespace APACIC {
  class Apacic {
  private:
    bool             test,ini,fin,had;
    ISR::ISR_Handler         * isr;
    BEAM::Beam_Spectra_Handler       * beam;
    APHYTOOLS::Selector_Data * seldata;
    APHYTOOLS::Parton_List   * partons;
    APHYTOOLS::Blob_List     * blobs;

    Hard_Processes           * AP_proc;
    Hard_Interface           * me_ps_interface;
  public:
    Apacic(bool);
    ~Apacic();
    void Init();
    void CrossSections();
    void GenerateEvents();
    void CleanUpEvent();
    void FillBlob(APHYTOOLS::Blob *);
    void TestShower();
  };
}

#endif
