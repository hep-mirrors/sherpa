#ifndef PHASIC__EWSudakov__Clustered__Calculator_H
#define PHASIC__EWSudakov__Clustered__Calculator_H

#include "ATOOLS/Phys/Flavour.H"
#include "AddOns/EWSudakov/EWSudakov_Calculator.H"

namespace EWSud {

  class Clustered_EWSudakov_Calculator {

  public:

    Clustered_EWSudakov_Calculator(PHASIC::Process_Base*);
    ~Clustered_EWSudakov_Calculator();
    ATOOLS::EWSudakov_Log_Corrections_Map CorrectionsMap(ATOOLS::Vec4D_Vector);
    const std::set<ATOOLS::EWSudakov_Log_Type>& ActiveLogTypes() const
    { return calculators.begin()->second->ActiveLogTypes(); }

  private:

    void AddCalculators(const ATOOLS::Cluster_Amplitude_UP&, size_t clusterings);
    void AddCalculator(const ATOOLS::Cluster_Amplitude_UP&, size_t clusterings);
    double CalcIClustered(const std::map<double, std::vector<long int>> restab,
                          const ATOOLS::Vec4D_Vector &mom,
                          const ATOOLS::Flavour_Vector &flavs);

    // the objects that will do the actual calculation, but do not need to know
    // anything about the clustering, which will be handled by
    // Clustered_EWSudakov_Calculator, before it calls its calculator to do the
    // actual work; note that there might be several possible clusterings (e.g.
    // eenunu -> ZZ/WW), such that we allow for the possibility of having
    // several calculators; the first calculator will always be the "base"
    // calculator, i.e. the one for the completely unclustered process
    std::map<ATOOLS::Flavour_Vector, std::unique_ptr<EWSudakov_Calculator>> calculators;

    PHASIC::Process_Base* proc;

    /// the default COMIX interface used to create clustered processes
    Comix_Interface* p_comixinterface;

    double m_resdist;
    bool m_disabled;

  };

}

#endif
