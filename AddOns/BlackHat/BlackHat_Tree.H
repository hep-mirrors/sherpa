#ifndef AddOns_BlackHat_BlackHat_Tree_H
#define AddOns_BlackHat_BlackHat_Tree_H

#include "PHASIC++/Process/Tree_ME2_Base.H"
#include "blackhat/BH_interface.h"
#include "blackhat/BH_Ampl.h"

namespace MODEL {
  class Model_Base;
  class Coupling_Data;
}

namespace PHASIC {
  struct External_ME_Args;
}

namespace BLACKHAT {

  class BlackHat_Tree: public PHASIC::Tree_ME2_Base {
  private:

    BH::BH_Ampl *p_ampl;

    static BH::BH_interface  *s_interface;

    size_t m_oqcd, m_oew, m_mode;
    double m_asfac, m_afac;

    void WarnForMassiveFlavours(const ATOOLS::Flavour_Vector&) const;

  public:

    static MODEL::Model_Base *s_model;

    BlackHat_Tree(const PHASIC::External_ME_Args& args,
		  BH::BH_Ampl *ampl,const int mode);

    ~BlackHat_Tree();

    double Calc(const ATOOLS::Vec4D_Vector& momenta);

    static void AddCouplings
    (const PHASIC::External_ME_Args &args,
     std::vector<std::vector<std::pair<std::string,int> > > &couplings,
     std::vector<std::pair<std::string,int> > cpls,size_t i=0);

    int OrderQCD(const int &id=-1) const;
    int OrderEW(const int &id=-1) const;

    inline static void SetInterface(BH::BH_interface *const bh) { s_interface=bh; }
    inline static void SetModel(MODEL::Model_Base *const md) { s_model=md; }

    inline static BH::BH_interface *Interface() { return s_interface; }

  };

}// end of namespace WHITEHAT

#endif
