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
    void boost(int,AMATOOLS::Vec4D,AMATOOLS::Vec4D&,AMATOOLS::Vec4D&);
    void rotat(int,AMATOOLS::Vec4D&,AMATOOLS::Vec4D,double**);
    
    void QcdAntenna(AMATOOLS::Vec4D* &,double);
    void BasicAntenna(AMATOOLS::Vec4D,AMATOOLS::Vec4D,AMATOOLS::Vec4D &,double*,double);
    void PermP(int,int* &);
    void Polytope(int,double*&);
  public:
    Sarge(int,int);
    ~Sarge();

    void GenerateWeight(AMATOOLS::Vec4D *,APHYTOOLS::Cut_Data *);
    void GeneratePoint(AMATOOLS::Vec4D *,APHYTOOLS::Cut_Data *);
    void GeneratePoint(AMATOOLS::Vec4D *,APHYTOOLS::Cut_Data *,double *);
    void CalculateS0(APHYTOOLS::Cut_Data *);
    std::string Name() {return std::string("Sarge");}
    int CountResonances(APHYTOOLS::Flavour*&) { return 0; }
  };
}
#endif
