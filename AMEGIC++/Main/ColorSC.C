#include "AMEGIC++/Main/ColorSC.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Data_Reader.H"

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
  double helpd;
  Data_Reader reader(" ",";","!","=");
  if (reader.ReadFromFile(helpd,"N_COLOR")) {
    Nc = helpd;
    msg_Out()<<"Set N_color="<<Nc<<"."<<std::endl;
    CF = 0.5*(Nc-1./Nc);
    CA = Nc;
  }
}
