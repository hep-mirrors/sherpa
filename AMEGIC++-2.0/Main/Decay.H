//Decay tables

#ifndef Decay_h
#define Decay_h

#include "Topology.H"

namespace AMEGIC {

  //! The DecayTable implements the partial widths and branching ratios.(ready)
  struct DecayTable {
    //! The width of the particular particle (fl[0])
    double width;
    //! The branching ratio
    double BR;
    //! fl[0] is the incoming flavour fl[1..fn] are the outgoings.
    APHYTOOLS::Flavour* fl;
    //! The number of flavours.
    int fn;
    //! Particle is decaying into the MORE list.
    DecayTable* More;
    //! Next particle in the row.
    DecayTable* Next;

  };
  
  //! Handles everything with decay widths and chains.(ready, to be extended)
  class Decay_Handler {
    //! The decay table, where all informations are stored.
    DecayTable* dt_first;
    //! Saves the Silent variable of Run_Parameter
    int save_silent;
    //! Checks, if a certain flavour is included in the list of vertices.
    int  Check_In_Vertex(Flavour);
    //! Checks, if a vertex with this flavour combination exists.
    int  Check_Vertex(Flavour,Flavour,Flavour);
    //! Adds a new DecayTable to the list.
    void Add(DecayTable*,Flavour,Flavour,Flavour);
    //! Finds every 1->2 decay for the unstable particles.
    void Find_Decay_Products_2();
    //! Finds every 1->3 decay.
    /*! Finds every 1->3 decay, where no 1->2 decay exists or 
      could not pass the mass constraint.
     */
    void Find_Decay_Products_3();
    //! Used for the determination of further decays.
    void Next_Decay(DecayTable*);
    //! Calculates the brnaching ratios recursively.
    void Rec_BR(DecayTable*,double);
    //! Calculates the branching ratios.
    void Branching_Ratios();
    //! Calculates the widths recursively.
    void Rec_Calc(Topology*,DecayTable*,int);
  public:
    //! Prints the table of decays.
    void Print();
    //! Determines all decays of unstable particles.
    void Find_Unstable();
    //! Calculates the width of the different decay modes.
    void Calculate(Topology*);
  };
}

#endif
