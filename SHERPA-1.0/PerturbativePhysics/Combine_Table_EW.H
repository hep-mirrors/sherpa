#ifndef Combine_Table_EW_H
#define Combine_Table_EW_H

#include "Combine_Table_Base.H"


namespace SHERPA {

  struct Vertex_Info {
    ATOOLS::Flavour m_fl;
    Complex * p_cpl;
    AMEGIC::Color_Function * p_color;
    unsigned int m_mode;
    inline Vertex_Info();
  };

  inline Vertex_Info::Vertex_Info() : p_cpl(NULL), p_color(NULL), m_mode(0) {}

  class Combine_Table_EW {
  public:
    Combine_Table_EW(ATOOLS::Jet_Finder * jf, ATOOLS::Vec4D * moms, 
		  Combine_Table * up, int isrmode, int isrshoweron, int mode);
    ~Combine_Table_EW();


    void            FillTable(Leg **,int, int); 
    Combine_Table * CalcJet(int, double, double, ATOOLS::Vec4D * moms=0);
    CD_Iterator     CalcPropagator(CD_Iterator & cit);
  };
}


#endif



