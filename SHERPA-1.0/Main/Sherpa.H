#ifndef Sherpa_h
#define Sherpa_h


#include "ME_PS_Interface.H"
#include "Soft_Interface.H"
#include "ISR_Handler.H"


#include "Amegic.H"
#include "Process_Base.H"


#include "Hard_Processes.H"

#include "Vector.H"
#include "Blob.H"
#include "Blob_List.H"
#include "Parton.H"
#include "Parton_List.H"

namespace SHERPA {
  class Sherpa {
    bool tune;

    AMEGIC::Amegic         * AM;
    ME_PS_Interface        * hard_interface;
    Soft_Interface         * soft_interface;
    
    ISR::ISR_Handler       * isr;
    
    APHYTOOLS::Parton_List * partons;
    APHYTOOLS::Blob_List   * blobs;
    
    // init a blob with a Hard event
    void FillBlob(APHYTOOLS::Blob *,AMEGIC::Process_Base *);
    // empty parton and blob buffers (to make room for a new event generation)
    void CleanUpEvent();

    bool RescaleJetrates();
    void RunPythiaTest();

    // routines for a simple parallelisation
    bool MPIInit();
    bool MPIGenerateEvents();

    bool SendBlobs();
    bool RecieveBlobs();
  public:
    Sherpa(): tune(0),AM(0), hard_interface(0), soft_interface(0),
	      isr(0),partons(0),blobs(0) {  } 
    ~Sherpa();
    // init all essentials (e.g. rpa, Flavour, msg, aqed, ..., Amegic )
    // returns 1 for successful
    bool Init();
    // calculate total crosssections (preparation for event selection)
    // returns 1 if process libs were existing and could be integrated () otherwise
    bool CrossSections();
    // generate event (as initialised above)
    // returns 1 for successful
    bool GenerateEvents();
    // GenerateOneEvent;
    void OneEvent();
    // Finish Event Simulation
    void Finalize();
  };
}
#endif

