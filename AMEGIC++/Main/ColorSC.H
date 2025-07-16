#ifndef AMEGIC_Main_ColorSC_H
#define AMEGIC_Main_ColorSC_H

namespace AMEGIC {

  class ColorSC {
  public:
    double Nc;
    double CF;
    double CA;
    double TR;
    bool init;
    ColorSC();
    void Init();
  };
  extern ColorSC CSC;


}// end of namespace PHASIC

#endif
