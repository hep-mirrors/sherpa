#include "CFPSHOWER++/Calculators/FF/SF_FF13.H"
#include "PHASIC++/Channels/CSS_Kinematics.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Org/Message.H"

using namespace CFPSHOWER;
using namespace PHASIC;
using namespace ATOOLS;

SF_FF13::SF_FF13(const Kernel_Info & info) : SF_Base(info) {}

double SF_FF13::Jacobean(const Splitting & split) const {
}

void SF_FF13::GeneratePoint(Splitting & split) const {
}

bool SF_FF13::InitKinematics(Splitting & split) const {
}

int SF_FF13::Construct(Splitting & split) const {
}
