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
    void boost(int,AMATOOLS::vec4d,AMATOOLS::vec4d&,AMATOOLS::vec4d&);
    void rotat(int,AMATOOLS::vec4d&,AMATOOLS::vec4d,double**);
    
    void QcdAntenna(AMATOOLS::vec4d* &,double);
    void BasicAntenna(AMATOOLS::vec4d,AMATOOLS::vec4d,AMATOOLS::vec4d &,double*,double);
    void PermP(int,int* &);
    void Polytope(int,double*&);
  public:
    Sarge(int,int);
    ~Sarge();

    void GenerateWeight(AMATOOLS::vec4d *,APHYTOOLS::Cut_Data *);
    void GeneratePoint(AMATOOLS::vec4d *,APHYTOOLS::Cut_Data *);
    void GeneratePoint(AMATOOLS::vec4d *,APHYTOOLS::Cut_Data *,double *);
    void CalculateS0(APHYTOOLS::Cut_Data *);
    std::string Name() {return std::string("Sarge");}
    int CountResonances(APHYTOOLS::Flavour*&) { return 0; }
  };
}
#endif
