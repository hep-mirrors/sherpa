#include "AMEGIC++/Main/ColorSC.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Default_Reader.H"

using namespace ATOOLS;
using namespace AMEGIC;

AMEGIC::ColorSC AMEGIC::CSC;

ColorSC::ColorSC()
{
  init = false;
  Nc = 3.;
  CF = 4./3.;
  CA = 3.;
  TR = 0.5;
}


void ColorSC::Init()
{
  if (init) return;
  init = true;
  Default_Reader reader;
  if (reader.Read(Nc, "N_COLOR", Nc)) {
    msg_Out()<<"Set N_color="<<Nc<<"."<<std::endl;
    CF = 0.5*(Nc-1./Nc);
    CA = Nc;
  }
}
