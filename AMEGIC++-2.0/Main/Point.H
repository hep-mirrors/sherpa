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
  APHYTOOLS::Flavour   fl;
  APHYTOOLS::Flavour * extrafl;
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
    Color   = new Color_Function;
    Lorentz = new Lorentz_Function;
    middle  = 0;
    ncpl = 4;
    cpl = new Complex[4];
    if (nextra>0) extrafl = new APHYTOOLS::Flavour[nextra]; 
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
      
      if (p.Color->Next!=0) {
	
	Color_Function* cforig = p.Color->Next;
	Color_Function* cfcopy = Color;
	while (cforig) {
	  cfcopy->Next = new Color_Function(*cforig);
	  cfcopy = cfcopy->Next;
	  cforig = cforig->Next;
	}
      }

      *Lorentz = *p.Lorentz;
 
      if (nextra>0) delete[] extrafl;
      nextra = p.nextra;
      if (nextra>0) {
	extrafl = new APHYTOOLS::Flavour[nextra]; 
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
}
#endif












