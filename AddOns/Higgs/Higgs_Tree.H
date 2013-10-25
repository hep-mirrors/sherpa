#ifndef AddOns_Higgs_Higgs_Tree_H
#define AddOns_Higgs_Higgs_Tree_H

#include "PHASIC++/Process/Tree_ME2_Base.H"

namespace MODEL {
  class Model_Base;
  class Coupling_Data;
}

namespace AMEGIC { class Basic_Sfuncs; }

namespace HIGGS {

  class Higgs_Tree: public PHASIC::Tree_ME2_Base {
  private:

    static MODEL::Model_Base *s_model;

    size_t m_oqcd, m_oew, m_n, m_mode, m_int, m_io, m_spin;
    double m_mh, m_gh, m_kg, m_kq;

    std::vector<Map_Info> m_hmap;
    std::vector<Complex>  m_resa, m_rest, m_resb, m_ress;

    AMEGIC::Basic_Sfuncs *p_bs;
    std::vector<int>      m_b;

  public:

    Higgs_Tree(const PHASIC::Process_Info& pi,
	       const ATOOLS::Flavour_Vector& flavs,
	       int mode,int io,int spin,
	       double kg,double kq);

    ~Higgs_Tree();

    double Calc(const ATOOLS::Vec4D_Vector& momenta);

    int OrderQCD(const int &id=-1);
    int OrderEW(const int &id=-1);

    double TR() const;

    std::vector<Complex> GetAmplitudes(const size_t &id);

    Complex GetPhase(const size_t &id);

    Complex GetHelicityPhase(const ATOOLS::Vec4D &pijt,const ATOOLS::Vec4D &eps1);

    void FillCombinations
    (std::set<std::pair<size_t,size_t> > &combs,
     std::map<size_t,ATOOLS::Flavour_Vector> &fls);

    std::vector<Map_Info> GetFlavourHelicityMap();

    inline static void SetModel(MODEL::Model_Base *const md) { s_model=md; }

  };// end of class Higgs_Tree

}// end of namespace HIGGS

#endif
