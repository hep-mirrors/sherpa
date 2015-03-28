#ifndef AMEGIC_Main_Point_H
#define AMEGIC_Main_Point_H

#include "MODEL/Main/Model_Base.H"
#include "AMEGIC++/Amplitude/Single_Vertex.H"
#include "AMEGIC++/Amplitude/Lorentz_Function.H"

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
    AMEGIC::Single_Vertex * v;
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
    bool CountT(int & tchan,const long unsigned int & kfcode=0);
    void GeneratePropID();
    std::string GetPropID() const;
    void FindOrder(std::vector<int> &);
  };

  std::ostream &operator<<(std::ostream &str,const Point &p);
}
#endif












