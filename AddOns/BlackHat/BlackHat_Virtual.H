#ifndef AddOns_BlackHat_BlackHat_Virtual_H
#define AddOns_BlackHat_BlackHat_Virtual_H

#include "PHASIC++/Process/Virtual_ME2_Base.H"
#include "blackhat/BH_interface.h"
#include "blackhat/BH_Ampl.h"

namespace MODEL { class Model_Base; }

namespace BLACKHAT {

  class BlackHat_Virtual: public PHASIC::Virtual_ME2_Base {
  private:

    BH::BH_Ampl *p_ampl;

    static BH::BH_interface  *s_interface;
    static MODEL::Model_Base *s_model;

  public:

    BlackHat_Virtual(const PHASIC::Process_Info& pi,
		     const ATOOLS::Flavour_Vector& flavs,
                     BH::BH_Ampl *ampl);

    ~BlackHat_Virtual();

    void Calc(const ATOOLS::Vec4D_Vector& momenta);
    double Eps_Scheme_Factor(const ATOOLS::Vec4D_Vector& mom);
    double ScaleDependenceCoefficient(const int i);

    static void AddCouplings
    (const PHASIC::Process_Info &pi,
     std::vector<std::vector<std::pair<std::string,int> > > &couplings,
     std::vector<std::pair<std::string,int> > cpls,size_t i=0);

    inline static void SetInterface(BH::BH_interface *const bh) { s_interface=bh; }
    inline static void SetModel(MODEL::Model_Base *const md) { s_model=md; }

    inline static BH::BH_interface *Interface() { return s_interface; }

  };

}// end of namespace WHITEHAT

#endif
