/**
 *  * NNPDF tutorial: Loading and using a NNPDF grid with the NNPDF Driver.
 *   *
 *    * Compile with:
 *     * g++ -O3 test.cc NNPDFDriver.cc -o prog
 *      *
 *       * Author: The NNPDF Collaboration - 2014
 *        */
 
#include "iostream"
#include "cmath"
#include "NNPDFDriver.h"
using namespace std;
 
int main() {
  NNPDFDriver *nnpdf = new NNPDFDriver("share/SHERPA-MC/NNPDF30_nlo_as_0118", 0);
  double x = 0.1, Q2 = 2; // x and Q^2
  for (int i=-6;i<7;++i){
    cout << i<< ": " << nnpdf->xfx(x,Q2,i) << endl; // prints  x*PDF
  }
  return 0;
}
