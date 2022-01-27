#include "SHRiMPS/Main/Hadron_Init.H"
#include "ATOOLS/Phys/KF_Table.H"
#include "ATOOLS/Org/Message.H"

using namespace SHRIMPS;
using namespace ATOOLS;
using namespace std;

void Hadron_Init::Init() {
  if(s_kftable.find(kf_pomeron)==s_kftable.end()) // if not initialized in amisic
    s_kftable[kf_pomeron]=new Particle_Info(kf_pomeron,0.0,0.0,0,0,1,0,"pomeron","pomeron");
  if(s_kftable.find(kf_reggeon)==s_kftable.end()) // if not initialized in amisic
    s_kftable[kf_reggeon]=new Particle_Info(kf_reggeon,0.0,0.0,0,0,1,0,"reggeon","reggeon");
  
  // Assume pp/ppbar collisions only
  if(s_kftable.find(kf_N_1440)==s_kftable.end()) // if not initialised
    s_kftable[kf_N_1440]=
      new Particle_Info(kf_N_1440,1.44,0.35,0,1,1,0,"N(1440)","N(1440)");
  if(s_kftable.find(kf_N_1440_plus)==s_kftable.end()) // if not initialised
    s_kftable[kf_N_1440_plus]=
      new Particle_Info(kf_N_1440_plus,1.44,0.35,3,1,1,0,"N(1440)+","N(1440)+");
  if(s_kftable.find(kf_N_1710)==s_kftable.end()) // if not initialised
    s_kftable[kf_N_1710]=
      new Particle_Info(kf_N_1710,1.71,0.12,0,1,1,0,"N(1710)","N(1710)");
  if(s_kftable.find(kf_N_1710_plus)==s_kftable.end()) // if not initialised
    s_kftable[kf_N_1710_plus]=
      new Particle_Info(kf_N_1710_plus,1.71,0.12,3,1,1,0,"N(1710)+","N(1710)+");
}

// we will need to add decay tables for the N1710's ....
// https://pdglive.lbl.gov/Particle.action?init=0&node=B014&home=BXXX005
