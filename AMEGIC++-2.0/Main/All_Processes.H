#ifndef Process_h
#define Process_h

#include "Process_Group.H"
#include "Single_Process.H"
#include "Topology.H"

namespace AMEGIC {
  //! The container for a list of processes.
  class All_Processes : public Process_Group {
  public:

    int  InitAllProcesses(Topology *,AMATOOLS::Vec4D *&,
			  std::vector<double> &,std::vector<Single_Process *> &);

    void Add(Process_Base *);
    void SelectOne();
    void RescaleXSec(double);

    bool CalculateTotalXSec();
    bool OneEvent();
    bool SameEvent();

    int                  Nin()     { return selected->Nin(); }
    int                  Nout()    { return selected->Nout(); } 
    APHYTOOLS::Flavour * Flavs()   { return selected->Flavs(); } 
    AMATOOLS::Vec4D    * Momenta() { return selected->Momenta(); }
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
