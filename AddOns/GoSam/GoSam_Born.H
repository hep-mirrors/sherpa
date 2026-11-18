#ifndef GoSam_Born_H
#define GoSam_Born_H

#include "PHASIC++/Process/Tree_ME2_Base.H"

#include "AddOns/GoSam/GoSam_Interface.H"

#include <map>
#include <set>
#include <utility>
#include <vector>

namespace GoSam {

  class GoSam_Born: public PHASIC::Tree_ME2_Base { 
  private:
    size_t m_pn, m_nin;
    bool m_active;
    int m_proc_id, m_oqcd, m_oqed;
    std::vector<double> m_momenta;
    std::array<double,4> m_result;
    double m_Norm;
    double m_g2qcd, m_g2qed;
    
    static ATOOLS::Flavour CombineQCD(const ATOOLS::Flavour &a,
                                      const ATOOLS::Flavour &b);

    void FillCombinations(const std::vector<std::pair<size_t,ATOOLS::Flavour> > &legs,
                          std::set<std::pair<size_t,size_t> > &combs,
                          std::map<size_t,ATOOLS::Flavour_Vector> &fls,
                          std::set<std::vector<size_t> > &visited) const;

  public:

    GoSam_Born(const PHASIC::External_ME_Args& args,
               bool active, int proc_id=-1);
    ~GoSam_Born() {}

    double Calc(const ATOOLS::Vec4D_Vector& momenta);
    int OrderQCD(const int &id=-1) const;
    int OrderEW(const int &id=-1) const;

    void FillCombinations(std::set<std::pair<size_t,size_t> > &combs,
                          std::map<size_t,ATOOLS::Flavour_Vector> &fls);
#ifdef VALIDATE_GOSAM_INTERFACE
    PHASIC::Tree_ME2_Base* p_check;
    bool m_errorOnFailedCheck;
    bool m_printAllChecks;
#endif
  };
}

#endif
