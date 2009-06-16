#include"SHERPA/PerturbativePhysics/Hard_Decay_Handler.H"

#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Phys/Decay_Map.H"
#include "ATOOLS/Phys/Decay_Table.H"

#include <iostream>
#include <iomanip>
#include <sstream>
#include <cstring>

using namespace SHERPA;
using namespace ATOOLS;

Hard_Decay_Handler::Hard_Decay_Handler(std::string _path,std::string _file)
{
}

Hard_Decay_Handler::~Hard_Decay_Handler() 
{
  delete p_decaymap;
}

void Hard_Decay_Handler::InitializeDecayMap()
{
  p_decaymap = new Decay_Map();
}
