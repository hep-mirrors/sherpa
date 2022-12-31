#include <iostream>
#include "Grid.h"
#include <stdlib.h> 

namespace Grids{

class Grids{

  static const int ReIm=2;
  static const int nLoop=2;
  static const int nEhel=2;
  static const int nEcoeff=9;
  int m_ngrids;

  std::vector<Grid> grid;

  public:

  Grids();
  ~Grids()
  {};
  
  double E[nEhel][nEcoeff+1][nLoop][ReIm];

  void Init_Grids(int ngrids);
  void Fill_Grids();
  void Fill_Grid(int);
  void Initialize_Interpolator();
  void Build_Interpolator(int);
  int Interpolate(double, double, double, double);
    
};
}

