#ifndef Point_h
#define Point_h

#include "Vertex.H"

namespace AMEGIC {

class Point {

public:
  int             number;
  int             b;
  int             t;
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

  Point(int _extra = 0) : nextra(_extra)  { 
    extrafl = 0;
    v       = 0;
    Color   = new Color_Function;
    Lorentz = new Lorentz_Function;
    middle  = 0;
    ncpl = 4;
    cpl = new Complex[4];
    if (nextra>0) extrafl = new ATOOLS::Flavour[nextra]; 
  }
  Point(const Point& copy) { 
    extrafl = 0;
    Color   = new Color_Function;
    Lorentz = new Lorentz_Function;
    middle  = 0;
    ncpl = 4;
    cpl = new Complex[4];
    if (nextra>0) extrafl = new ATOOLS::Flavour[nextra]; 

    *this = copy;
  } 

  ~Point() {
    if (nextra>0) delete[] extrafl;
    delete[] cpl; 
    delete Color;
    delete Lorentz;
  }
  Point& operator=(const Point& p) {
    if (this!=&p) {
      number = p.number;
      b      = p.b;
      t      = p.t;
      m      = p.m;
      fl     = p.fl;
      
      *Color = *p.Color;
      *Lorentz = *p.Lorentz;
 
      if (nextra>0) delete[] extrafl;
      nextra = p.nextra;
      if (nextra>0) {
	extrafl = new ATOOLS::Flavour[nextra]; 
	for(short int i=0;i<nextra;i++) extrafl[i] = p.extrafl[i];
      }
      left   = p.left;
      right  = p.right;
      middle = p.middle;
      prev  = p.prev;
      v = p.v;
      //cpl's
      if (ncpl!=p.ncpl) {
	delete[] cpl;
	ncpl = p.ncpl;
	cpl = new Complex[ncpl];
      } 

      for(short int i=0;i<ncpl;i++) cpl[i] = p.cpl[i];
    }
    return *this;
  }
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
}
#endif












