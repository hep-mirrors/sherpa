#ifndef AMEGIC_Amplitude_Zfunctions_MHVCalculator_H
#define AMEGIC_Amplitude_Zfunctions_MHVCalculator_H

#define Basic_Sfuncs_In_MHV
#ifdef Basic_Sfuncs_In_MHV
#include "AMEGIC++/Amplitude/Zfunctions/Basic_Sfuncs.H"
#else
#include "AMEGIC++/Amplitude/Zfunctions/Mom.H"
#endif

#include "ATOOLS/Math/MyComplex.H"


namespace AMEGIC {
  
#ifdef Basic_Sfuncs_In_MHV
  class Basic_Sfuncs;
  typedef  Basic_Sfuncs MomentumList;
#else
  class MomentumList;
#endif
  
  class MHVCalculator {
    int n_part;
    int *m_dummyarg, *m_ndummyarg; 
    int *m_dummysl, *m_ndummysl; 
    int *m_plist, *m_signlist;
    int m_qlist[9];
    
    MomentumList *p_BS;
    
    inline int CountSign(int* signlist,int start,int end,int sg);
    void Make_Qlist(int* perm,int* plist,int* qlist,int part);
    bool Check_Qlist(int* perm,int* signlist,int* qlist);
    
    Complex Elementary_MHV_Amplitude(int* perm,int* signlist,int part);
    Complex NMHV_Amplitude(int* perm,int* signlist,int part,int vhel);
    Complex NNMHV_Amplitude(int* perm,int* signlist,int part,int vhel);	
    Complex Elementary_MHVbar_Amplitude(int* perm,int* signlist,int part);
    Complex NMHVbar_Amplitude(int* perm,int* signlist,int part,int vhel);
    Complex NNMHVbar_Amplitude(int* perm,int* signlist,int part,int vhel);
    
    Complex Elementary_MHVQ2_Amplitude(int* perm,int* signlist,int* qlist,int part);
    Complex NMHVQ2_Amplitude(int* perm,int* signlist,int* qlist,int part,int vhel);
    Complex NNMHVQ2_Amplitude(int* perm,int* signlist,int* qlist,int part,int vhel);
    Complex Elementary_MHVQ2bar_Amplitude(int* perm,int* signlist,int* qlist,int part);
    Complex NMHVQ2bar_Amplitude(int* perm,int* signlist,int* qlist,int part,int vhel);
    Complex NNMHVQ2bar_Amplitude(int* perm,int* signlist,int* qlist,int part,int vhel);
    
    Complex Elementary_MHVQ4_Amplitude(int* perm,int* signlist,int* qlist,int part);	
    Complex NMHVQ4_Amplitude(int* perm,int* signlist,int* qlist,int part,int vhel);
    Complex NNMHVQ4_Amplitude(int* perm,int* signlist,int* qlist,int part,int vhel);
    Complex Elementary_MHVQ4bar_Amplitude(int* perm,int* signlist,int* qlist,int part);
    Complex NMHVQ4bar_Amplitude(int* perm,int* signlist,int* qlist,int part,int vhel);
    Complex NNMHVQ4bar_Amplitude(int* perm,int* signlist,int* qlist,int part,int vhel);

    
  public:
    
    inline MHVCalculator() {}
    MHVCalculator(int part,int* plist);
    
    ~MHVCalculator();
    
    Complex Differential(int* perm, int* signlist);
    
    inline const int * GetQlist() { return m_qlist;}
    inline const int * GetPlist() { return m_plist;}
    void SetMomentumList(MomentumList* BS) { p_BS=BS; }
  };
  



////////////////////////////////////////////////////////////////////////////////////////////////

  inline int MHVCalculator::CountSign(int* signlist,int start,int end,int sg)
  {
    int sum = 0;
    for (int i=start;i<end;i++) sum+=signlist[i];
    if (sg==-1) return (end-start-sum)/2;
    return (end-start+sum)/2;
  }


}
#endif
