#ifndef RamboKK_h
#define RamboKK_h

#include "Single_Channel.H"

namespace PHASIC {
  class RamboKK : public Single_Channel {
    int     Nin,Nout;
    double  Z_N;
    double* xm2;
    double* p2;
    double* E;
    double* ms_out;
    int     massflag;
    int     kkp,ed;
    double  R2,maxM2,maxN,mpss;
  public:
    RamboKK(int,int,APHYTOOLS::Flavour*);
    void Set_KKmass();
    void Generate_Point(AMATOOLS::Vec4D *,APHYTOOLS::Cut_Data *);
    void Generate_Point(AMATOOLS::Vec4D *,APHYTOOLS::Cut_Data *,double*);
    void Generate_Weight(AMATOOLS::Vec4D *,APHYTOOLS::Cut_Data *);  
    void Massive_Weight(AMATOOLS::Vec4D *,double);
    void Massive_Point(AMATOOLS::Vec4D *,double);
    std::string Name() {return std::string("RamboKK");}
  };  
}
#endif
