#ifndef Decay_h
#define Decay_h

#include "Process_Group.H"
#include "Single_Process.H"
#include "Topology.H"



namespace AMEGIC {
  /*
  //! The DecayTable implements the partial widths and branching ratios.
  struct DecayTable {
    //! The width of the particular particle (fl[0])
    double width;
    //! The branching ratio
    double BR;
    //! fl[0] is the incoming flavour fl[1..fn] are the outgoings.
    ATOOLS::Flavour* fl;
    //! The number of flavours.
    int fn;
    //! Particle is decaying into the MORE list.
    DecayTable* More;
    //! Next particle in the row.
    DecayTable* Next;

  };
  
  //! Handles everything with decay widths and chains.
  class Decay_Handler {
    //! The decay table, where all informations are stored.
    DecayTable* dtfirst;
    //! Checks, if a certain flavour is included in the list of vertices.
    int  CheckInVertex(ATOOLS::Flavour);
    //! Checks, if a vertex with this flavour combination exists.
    int  CheckVertex(ATOOLS::Flavour,ATOOLS::Flavour,ATOOLS::Flavour);
    //! Adds a new DecayTable to the list.
    void Add(DecayTable*,ATOOLS::Flavour,ATOOLS::Flavour,ATOOLS::Flavour);
    //! Finds every 1->2 decay for the unstable particles.
    void FindDecayProducts2();
    //! Finds every 1->3 decay.
    //! Finds every 1->3 decay, where no 1->2 decay exists or 
    //!could not pass the mass constraint.
    void FindDecayProducts3();
    //! Used for the determination of further decays.
    void NextDecay(DecayTable*);
    //! Calculates the brnaching ratios recursively.
    void RecBR(DecayTable*,double);
    //! Calculates the branching ratios.
    void BranchingRatios();
    //! Calculates the widths recursively.
    void RecCalc(Topology*,DecayTable*,int);
  public:
    //! Prints the table of decays.
    void Print();
    //! Determines all decays of unstable particles.
    void FindUnstable();
    //! Calculates the width of the different decay modes.
    void Calculate(Topology*);
  };
  */
}

#endif
