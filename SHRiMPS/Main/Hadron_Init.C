#include "SHRiMPS/Main/Hadron_Init.H"
#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Org/Message.H"

using namespace SHRIMPS;
using namespace ATOOLS;
using namespace std;

void Hadron_Init::Init() {
  if(s_kftable.find(kf_pomeron)==s_kftable.end()) // if not initialized in amisic
    s_kftable[kf_pomeron]=new Particle_Info(kf_pomeron,0.0,0.0,0,0,0,1,0,"pomeron","pomeron");
  if(s_kftable.find(kf_reggeon)==s_kftable.end()) // if not initialized in amisic
    s_kftable[kf_reggeon]=new Particle_Info(kf_reggeon,0.0,0.0,0,0,0,1,0,"reggeon","reggeon");
}
