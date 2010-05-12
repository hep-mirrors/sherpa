#ifndef AMEGIC_Main_Point_H
#define AMEGIC_Main_Point_H

#include "MODEL/Interaction_Models/Vertex.H"

namespace AMEGIC {

  class Point {
    
    void ResetProps(int&);
  public:
    int             number;
    int             b;
    int             t;  // t=10: intermediate particle in decay treatment
    int             zwf;
    int             m;
    int             propid;
    ATOOLS::Flavour   fl;
    ATOOLS::Flavour * extrafl;
    int                  nextra;
    Point              * left;
    Point              * right;
    Point              * middle;
    Point              * prev;
    MODEL::Single_Vertex * v;
    std::vector<Complex> cpl;
    MODEL::Color_Function     * Color;
    MODEL::Lorentz_Function   * Lorentz;
    
    Point(int extra = 0);
    Point(const Point& copy);
    
    Point& operator=(const Point& p);
    
    ~Point() {
      if (nextra>0) delete[] extrafl;
      delete Color;
      if (Lorentz) delete Lorentz;
    }
    void Print();
    void ResetExternalNumbers(int);
    void ResetProps();
    void ResetFlag();
    Point* CopyList(Point* p);
    size_t Ncpl() const {return cpl.size();}
    int CountKK();
    void GeneratePropID();
    std::string GetPropID() const;
    int FindQCDOrder(int &);
    int FindQEDOrder(int &);
  };
}
#endif












