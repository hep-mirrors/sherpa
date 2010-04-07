#ifndef PHASIC_Process_Single_Process_H
#define PHASIC_Process_Single_Process_H

#include "PHASIC++/Process/Process_Base.H"

namespace PHASIC {

  class Single_Process: public Process_Base {
  public:

    Single_Process();
    virtual ~Single_Process();
    
    size_t Size() const;
    Process_Base *operator[](const size_t &i);

    void DeSelect();
    bool SelectOne();

    Weight_Info *OneEvent();
    Weight_Info *WeightedEvent(const int mode=0);

    double KFactor();
    double KFactor2();

    bool CalculateTotalXSec(const std::string &resultpath,
			    const bool create);
    void SetLookUp(const bool lookup);

    virtual bool Combinable(const size_t &idi,
			    const size_t &idj);
    virtual const ATOOLS::Flavour_Vector &
    CombinedFlavour(const size_t &idij);
    double BeamWeight(const double& Q2);

  };// end of class Single_Process

}// end of namespace ATOOLS

#endif
