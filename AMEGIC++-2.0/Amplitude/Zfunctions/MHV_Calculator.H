#ifndef MHV_Calculator_h
#define MHV_Calculator_h

#include "Basic_Sfuncs.H"

namespace AMEGIC {

  class MHV_Calculator {
    int n_part;
    int * m_arg;
    int * m_dummyarg;
    Basic_Sfuncs * p_BS;
    int CountSign(int* perm,int* signlist,int part,int sg);
    Complex Elementary_MHV_Amplitude(int* perm,int* signlist,int part,int psign=1);
    Complex Elementary_MHVbar_Amplitude(int* perm,int* signlist,int part,int psign=1);
    Complex MHV_Amplitude(int* perm,int* signlist,int vhel);
    Complex MHVbar_Amplitude(int* perm,int* signlist,int vhel);
  public:
    MHV_Calculator(int,Basic_Sfuncs *);
    ~MHV_Calculator(); 
    
    double Differential(int* perm,int* signlist);
  };  

}
#endif
