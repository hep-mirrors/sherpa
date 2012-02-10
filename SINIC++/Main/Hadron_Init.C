#include "SINIC++/Main/Hadron_Init.H"
#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Org/Message.H"

using namespace SINIC;
using namespace ATOOLS;
using namespace std;

void Hadron_Init::Init() {
  if(s_kftable.find(990)==s_kftable.end()) // if not initialized in amisic
    s_kftable[990]=new Particle_Info(990,0.0,0.0,0,0,0,1,0,"pomeron","pomeron");
  if(s_kftable.find(110)==s_kftable.end()) // if not initialized in amisic
    s_kftable[110]=new Particle_Info(110,0.0,0.0,0,0,0,1,0,"reggeon","reggeon");
}
