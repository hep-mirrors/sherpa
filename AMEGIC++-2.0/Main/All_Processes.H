#ifndef Process_h
#define Process_h

#include "Process_Group.H"
#include "Single_Process.H"

#include "Topology.H"

namespace AMEGIC {
  //! The container for a list of processes.

  class All_Processes : public Process_Group {
  public:

    int  InitAllProcesses(Interaction_Model_Base *,Topology *,ATOOLS::Vec4D *&);

    void Add(Process_Base *);
    void SelectOne();
    void RescaleXSec(double);
    void SetupEnhance();

    bool CalculateTotalXSec(std::string _resdir=std::string(""));
    bool OneEvent(double = -1.);
    bool SameEvent();

    int  NIn()  { return p_selected->NIn();  }
    int  NOut() { return p_selected->NOut(); } 

    const ATOOLS::Flavour *Flavours() { return p_selected->Flavours(); } 
    const ATOOLS::Vec4D   *Momenta()  { return p_selected->Momenta();  }

  };

}

#endif



    /*!
      Adds a process to the list. Here, no care is taken for any number of legs etc..
    */
    /*!
      InitAllProcesses iterates over all constituents of the process list and calls their
      InitAmplitudes method. This of course might be either a Process_Group or a 
      Single_Process. Then, via InitLibraries, all libraries needed are linked or
      written to a file. Following that the integrators of all Single_Processes that
      have no partner are set up, i.e. all potentiall occuring integrators are either loaded in
      or constructed. Finally, AllProcesses is in the position to set up the integrators of 
      all constituents.
    */
    /*!
      To calculate total x-sections. 
      If this is successful a 'true' is returned.
    */
