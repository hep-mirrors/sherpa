#ifndef Point_h
#define Point_h

#include "Vertex.H"

namespace AMEGIC {

class Point {

  void ResetProps(int&);
public:
  int             number;
  int             b;
  int             t;  // t=10: intermediate particle in decay treatment
  int             m;
  ATOOLS::Flavour   fl;
  ATOOLS::Flavour * extrafl;
  int                  nextra;
  Point              * left;
  Point              * right;
  Point              * middle;
  Point              * prev;
  Single_Vertex      * v;
  Complex            * cpl;
  int                  ncpl;
  Color_Function     * Color;
  Lorentz_Function   * Lorentz;

  Point(int extra = 0);
  Point(const Point& copy);

  Point& operator=(const Point& p);

  ~Point() {
    if (nextra>0) delete[] extrafl;
    delete[] cpl; 
    delete Color;
    delete Lorentz;
  }
  void Print();
  void ResetExternalNumbers(int);
  void ResetProps();
  void ResetFlag();
  Point* CopyList(Point* p);
};

 // Format of a Point to transport via MPI:
  struct MPI_Point {
    
    int                   m_fl;  
    MPI_Lorentz_Function  m_lf;
    MPI_Color_Function    m_cf;
    MPI_Single_Vertex     m_v;
    double               m_cpl[8];
    int                  m_left, m_middle, m_right, m_num, m_b;
  };
  
  //for MPI usage
  
  void Point2MPI(const Point * , MPI_Point &);
  void Point2MPI(const Point * , MPI_Point *);
  int Point2MPI(const Point * , MPI_Point *,int &);
  
  Point * MPI2Point(const MPI_Point &, Point &);
  Point * MPI2Point(const MPI_Point *, Point *);
  Point * MPI2Point(const MPI_Point *, Point *, int);


  std::ostream & operator<<(std::ostream &, const MPI_Point &);
  std::ostream & operator<<(std::ostream &, const Point *);
}
#endif












