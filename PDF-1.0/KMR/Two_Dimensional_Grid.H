#ifndef Two_Dimensional_Grid_H
#define Two_Dimensional_Grid_H

#include "Node.H"

namespace PDF {

  class Two_Dimensional_Grid {
  private:

    ATOOLS::Node<double> *p_root, *p_xl, *p_xr;
    ATOOLS::Node<double> *p_yll, *p_ylr, *p_yrl, *p_yrr;

    std::vector<double> m_intx;
    double              m_inty;

    size_t FindX(const double x);
    size_t FindY(const double y);

    void InterpolateX(const double x) const;
    void InterpolateY(const double y) const;

  public:

    // constructor
    Two_Dimensional_Grid(const size_t inipoints);

    //destructor
    ~Two_Dimensional_Grid();

    // member functions
    double Value(const double x1,const double x2);
    double SetValue(const double x1,const double x2,const double y);

    // inline functions
    inline double operator()(const double x1,const double x2) 
    { return Value(x1,x2); }

  };// end of class Two_Dimensional_Grid_H

}// end of namespace PDF

#endif
