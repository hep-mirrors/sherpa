#include "EXTRA_XS/Main/ME2_Base.H"
#include "Exception.H"

using namespace EXTRAXS;
using namespace ATOOLS;

namespace EXTRAXS {
  class Drell_Yan_EW : public ME2_Base {
  public:
    Drell_Yan_EW(const Single_Process& proc);
    double operator()(const Vec4D_Vector& momenta);
  };
}


Drell_Yan_EW::Drell_Yan_EW(const Single_Process& proc) :
  ME2_Base(proc)
{
}


double Drell_Yan_EW::operator()(const Vec4D_Vector& momenta) {
  return 1.0;
}

DECLARE_ME2_GETTER(0,Drell_Yan_EW, "d_db_e-_e+_2_0_1_0")
DECLARE_ME2_GETTER(1,Drell_Yan_EW, "u_ub_e-_e+_2_0_1_0")
DECLARE_ME2_GETTER(2,Drell_Yan_EW, "s_sb_e-_e+_2_0_1_0")
DECLARE_ME2_GETTER(3,Drell_Yan_EW, "c_cb_e-_e+_2_0_1_0")
DECLARE_ME2_GETTER(4,Drell_Yan_EW, "b_bb_e-_e+_2_0_1_0")

DECLARE_ME2_GETTER(5,Drell_Yan_EW, "d_db_mu-_mu+_2_0_1_0")
DECLARE_ME2_GETTER(6,Drell_Yan_EW, "u_ub_mu-_mu+_2_0_1_0")
DECLARE_ME2_GETTER(7,Drell_Yan_EW, "s_sb_mu-_mu+_2_0_1_0")
DECLARE_ME2_GETTER(8,Drell_Yan_EW, "c_cb_mu-_mu+_2_0_1_0")
DECLARE_ME2_GETTER(9,Drell_Yan_EW, "b_bb_mu-_mu+_2_0_1_0")

DECLARE_ME2_GETTER(10,Drell_Yan_EW, "d_db_tau-_tau+_2_0_1_0")
DECLARE_ME2_GETTER(11,Drell_Yan_EW, "u_ub_tau-_tau+_2_0_1_0")
DECLARE_ME2_GETTER(12,Drell_Yan_EW, "s_sb_tau-_tau+_2_0_1_0")
DECLARE_ME2_GETTER(13,Drell_Yan_EW, "c_cb_tau-_tau+_2_0_1_0")
DECLARE_ME2_GETTER(14,Drell_Yan_EW, "b_bb_tau-_tau+_2_0_1_0")
