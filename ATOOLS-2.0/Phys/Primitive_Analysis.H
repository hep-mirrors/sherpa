#ifndef Primitive_Analysis_H
#define Primitive_Analysis_H

#include "Primitive_Observable_Base.H"
#include <vector>

#include "Parton_List.H"
#include "Blob_List.H"

namespace APHYTOOLS {
  class Primitive_Analysis {
  private :
    std::vector<Primitive_Observable_Base *>   observables;
    std::string                                name;
  public :
    Primitive_Analysis(std::string);
    Primitive_Analysis();
    ~Primitive_Analysis();

    // add a single Observable to the list
    void AddObservable(Primitive_Observable_Base *);
    // most simple analysis, assume the observables a know how to handle a single value
    void DoAnalysis(double);
    // simple analysis based on a list of partons (basically flavs, moms)
    void DoAnalysis(const APHYTOOLS::Parton_List &, double);
    // advanced analysis, working with a bloblist, thus event generation properties can be exploit
    void DoAnalysis(const APHYTOOLS::Blob_List &, double);
    // finialise and write out results of analysis
    void FinishAnalysis(std::string,int);
  };
}

#endif
