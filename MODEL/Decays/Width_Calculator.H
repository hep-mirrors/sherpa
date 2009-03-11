#ifndef Width_Calculator_h
#define Width_Calculator_h

#include "Flavour.H"
#include "Decay_Table.H"
#include "Model_Base.H"
#include "Vertex.H"

namespace MODEL {
  class Width_Calculator {
  private :
    double m_M,m_m[3];

    Vertex_Table * p_vertextable;
    Single_Vertex * FindVertex(const ATOOLS::Flavour &,
			       const ATOOLS::Flavour &,const ATOOLS::Flavour &);

    double TwoBodyPref(const double M,const double m1,const double m2);
    double Lambda2(const double M,const double m1,const double m2);
    double ColorFactor(Single_Vertex *);
    double Norm(ATOOLS::Flavour);

    double SFF(ATOOLS::Decay_Channel * dec);
    double SVV(ATOOLS::Decay_Channel * dec);
    double SVS(ATOOLS::Decay_Channel * dec);
    double SSS(ATOOLS::Decay_Channel * dec);
    double VFF(ATOOLS::Decay_Channel * dec);
    double VSS(ATOOLS::Decay_Channel * dec);
    double VVS(ATOOLS::Decay_Channel * dec);
    double VVV(ATOOLS::Decay_Channel * dec);
    double FFS(ATOOLS::Decay_Channel * dec);
    double FFV(ATOOLS::Decay_Channel * dec);
    double FSF(ATOOLS::Decay_Channel * dec);
    double FVF(ATOOLS::Decay_Channel * dec);
  public:
    Width_Calculator(Vertex_Table * vertextable) : p_vertextable(vertextable) {}

    double Width(ATOOLS::Decay_Channel * dec);
  };
}

#endif
