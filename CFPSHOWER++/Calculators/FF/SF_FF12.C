#include "CFPSHOWER++/Calculators/FF/SF_FF12.H"
#include "CFPSHOWER++/Calculators/FF/Kinematics_FF2.H"
#include "PHASIC++/Channels/CSS_Kinematics.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Org/Message.H"

using namespace CFPSHOWER;
using namespace PHASIC;
using namespace ATOOLS;

SF_FF12::SF_FF12(const Kernel_Info & info) :
  SF_Base(info) {
  switch (m_logtype) {
  case log_type::soft:
    p_kinematics = new Kinematics_FF2_soft();
    break;
  case log_type::coll:
    p_kinematics = new Kinematics_FF2_coll();
    break;
  default:
    msg_Error()<<"Error in "<<METHOD<<":\n"
	       <<"   No suitable type found to define kinematics mapping.\n"
	       <<"   Will exit the run.\n";
    exit(0);
  }
  m_moms.resize(2);
}

double SF_FF12::Jacobean(const Splitting & split) const {
  return p_kinematics->Weight();
}

bool SF_FF12::Construct(Splitting * split,Configuration * config,const int & mode) {
  return (*p_kinematics)(split,config,mode);
}
