#ifndef METOOLS_HadronCurrents_Line_Shapes_H
#define METOOLS_HadronCurrents_Line_Shapes_H

#include "METOOLS/HadronCurrents/Resonance_Base.H"
#include <map>

namespace METOOLS {
  class Line_Shapes {
  private:
    std::map<long int,Total_Width_Base *> m_lineshapes;
  public:
    Line_Shapes();
    ~Line_Shapes();

    void Init();
    void Add(const ATOOLS::Flavour & inflav,
	     Total_Width_Base * lineshape);
    Total_Width_Base * Get(const ATOOLS::Flavour & inflav);
  };

  extern Line_Shapes * LineShapes;
}



#endif

