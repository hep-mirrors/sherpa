#ifndef Sarge_h
#define Sarge_h

#include "Single_Channel.H"

namespace PHASIC {
  class Sarge : public Single_Channel {
    int      nin,nout,n_xi;
    double   s0;
    double*  x;
    double** rot;
    
    //temporarily
    void boost(int,ATOOLS::Vec4D,ATOOLS::Vec4D&,ATOOLS::Vec4D&);
    void rotat(int,ATOOLS::Vec4D&,ATOOLS::Vec4D,double**);
    
    void QcdAntenna(ATOOLS::Vec4D* &,double);
    void BasicAntenna(ATOOLS::Vec4D,ATOOLS::Vec4D,ATOOLS::Vec4D &,double*,double);
    void PermP(int,int* &);
    void Polytope(int,double*&);
  public:
    Sarge(int,int);
    ~Sarge();

    void GenerateWeight(ATOOLS::Vec4D *,ATOOLS::Cut_Data *);
    void GeneratePoint(ATOOLS::Vec4D *,ATOOLS::Cut_Data *);
    void GeneratePoint(ATOOLS::Vec4D *,ATOOLS::Cut_Data *,double *);
    void CalculateS0(ATOOLS::Cut_Data *);
    std::string Name() {return std::string("Sarge");}
    int CountResonances(ATOOLS::Flavour*&) { return 0; }
  };
}
#endif
