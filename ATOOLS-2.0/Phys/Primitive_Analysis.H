#ifndef Primitive_Analysis_H
#define Primitive_Analysis_H

#include "Primitive_Observable_Base.H"
#include <vector>

#include "Particle_List.H"
#include "Blob_List.H"

namespace ATOOLS {
  class Primitive_Analysis {
  private :
    std::vector<Primitive_Observable_Base *>   observables;
    std::string                                name;
  public :
    Primitive_Analysis(std::string);
    Primitive_Analysis();
    ~Primitive_Analysis();

    void AddObservable(Primitive_Observable_Base *);
    void DoAnalysis(double);
    void DoAnalysis(const Particle_List &, double);
    void DoAnalysis(const Blob_List &, double);
    void FinishAnalysis(std::string,int=0);
    void SetBlobType(std::string);
  };
}

    // add a single Observable to the list
    // most simple analysis, assume the observables a know how to handle a single value
    // simple analysis based on a list of particles (basically flavs, moms)
    // advanced analysis, working with a bloblist, thus event generation properties can be exploit
    // finialise and write out results of analysis
#endif
