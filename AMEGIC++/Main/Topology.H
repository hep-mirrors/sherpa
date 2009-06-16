#ifndef AMEGIC_Main_Topology_H
#define AMEGIC_Main_Topology_H

#include "AMEGIC++/Main/Point.H"

namespace AMEGIC {
  //! A list of topologies with the same number of legs.
  struct Single_Topology {
    //! Number of topologies.
    int number;
    //! Depth of the topology.
    int depth;
    //! The list of point lists.
    Point** p;//[number][depth]
  };
  
  //! Creates and handles topologies.
  class Topology {
    //! Maximum number of topologies
    int              ntop;
    //! The list of topology lists, one entry equals a certain number of legs.
    Single_Topology* top;
    //! Build all topologies to a maximum number of legs.
    void   Build_All(int);
    //! Build all topologies for a certain number of legs.
    void   Build_Single(int,Single_Topology*);
  public:
    //! Zero Constructor.
    Topology(); 
    //! Constructs all topologies.
    Topology(int); 
    //! Destructs all topologies.
    ~Topology();
    //! Returns a list of topologies with a certain number of legs.
    inline Single_Topology* Get(int N) {return &top[N];}
    //! Copies a given point list into another.
    Point* Copy(Point*,Point*,int&);
    //! Prints a topology.
    void   Print(Point*);
  };
}
#endif











