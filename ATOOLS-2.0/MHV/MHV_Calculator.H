#ifndef MHV_Calculator_h
#define MHV_Calculator_h

#include "Mom.H"

#ifndef Mom_H
#include "Basic_Sfuncs.H"
#endif



namespace AMEGIC {

#if defined(Mom_H)
class MomentumList;
typedef MomentumList Basic_Sfuncs;
#endif

  class MHV_Calculator {
    int n_part;
    int * m_arg;
    int * m_dummyarg;
  
    Basic_Sfuncs* p_BS;

    int * m_signlist;
    int * m_plist;
    int m_qlist[5];
    int print;
 

    int CountSign(int* signlist,int start,int end,int sg);
    Complex Elementary_MHV_Amplitude(int* perm,int* signlist,int part);
    Complex Elementary_MHVbar_Amplitude(int* perm,int* signlist,int part);
    Complex MHV_Amplitude(int* perm,int* signlist,int vhel);
    Complex MHVbar_Amplitude(int* perm,int* signlist,int vhel);
    
    Complex Elementary_MHVQ2_Amplitude(int* perm,int* signlist,int* qlist,int part);
    Complex Elementary_MHVQ2bar_Amplitude(int* perm,int* signlist,int* qlist,int part);
    Complex MHVQ2_Amplitude(int* perm,int* signlist,int vhel);
    Complex MHVQ2bar_Amplitude(int* perm,int* signlist,int vhel);
    Complex MHVQ4_Amplitude(int* perm,int* signlist,int vhel);
    Complex MHVQ4bar_Amplitude(int* perm,int* signlist,int vhel);
    Complex Elementary_MHVQ4_Amplitude(int* perm,int* signlist,int* qlist,int part);
    Complex Elementary_MHVQ4bar_Amplitude(int* perm,int* signlist,int* qlist,int part);
    
    void Make_Qlist(int* perm,int* plist,int* qlist,int part);
    bool Check_Qlist(int* perm,int* signlist,int* plist,int* qlist);
      
  public:

    inline MHV_Calculator() {}

    MHV_Calculator(int part,Basic_Sfuncs *BS,int* signlist,int* plist,int print=0);

    ~MHV_Calculator();

    Complex Differential(int* perm);
    int* GetQlist();
    int* GetPlist();

  };  
  
}
#endif
