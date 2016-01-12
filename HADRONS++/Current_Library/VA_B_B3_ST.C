namespace HADRONS {
namespace VA_B_B3_FFs {

  class ST : public FormFactor_Base {
    int m_mode;
    double mq, mQ, msigma;
    double a, ap;
    double IS(double w, double factor, double exp1, double exp2);
  public:
    ST(struct GeneralModel model, double* masses, string variant);
    void CalcFFs(ATOOLS::Vec4D p0, ATOOLS::Vec4D p1);
  };

  ST::ST(GeneralModel model, double* masses, string variant) : 
    FormFactor_Base(model, masses)
  {
    m_mode = int(model("mode",0.0)+0.5);
    string prefix=variant+"_";
    mq = model(prefix+"mq",1.0);
    mQ = model(prefix+"mQ",1.0);
    msigma = model(prefix+"msigma",1.0);
    a = model(prefix+"alpha",1.0);
    ap = model(prefix+"alpha_prime",1.0);
  }

  void ST::CalcFFs(Vec4D p0, Vec4D p1)
  {
    Vec4D v0 = p0/m0;
    Vec4D v1 = p1/m1;
    double w = v0*v1;

    double aaprime = 0.5*(a+ap);
    double msigma2 = sqr(msigma);

    double I;
    switch(m_mode) {
    case 0:
      // Lambda_x (1/2+) -> Lambda_y or N (3/2-)
      I = IS(w, -sqrt(2.0)/3.0, 2.5, 3.0);
      m_V1 = I*3.0*msigma/a*(1.0+msigma/aaprime*(ap/mq+a/mQ));
      m_V2 = -I*(3.0*msigma2/mq*ap/aaprime/a-msigma/4.0/mq/mQ*(a-3.0*ap));
      m_V3 = -I*(3.0*msigma2/mQ/aaprime+aaprime/4.0/mQ);
      m_V4 = I*aaprime/2.0/mQ;
      m_A1 = I*(3.0*msigma/a-aaprime/4.0/mQ+msigma/60.0/mq/mQ*(5.0*a-23.0*ap));
      m_A2 = -I*(3.0*msigma2/mq*ap/a/aaprime-msigma/60.0/mq/mQ*(5.0*a-11.0*ap)+18.0*msigma2*msigma*ap/7.0/sqr(aaprime)/mq/mQ);
      m_A3 = I/mQ*(3.0*msigma2/aaprime+aaprime/4.0+msigma*ap/5.0/mq+18.0*msigma2*msigma*ap/7.0/sqr(aaprime)/mq);
      m_A4 = -I/mQ*(aaprime/2.0+2.0*msigma*ap/5.0/mq);
      break;
    case 1:
      // Lambda_x (1/2+) -> Lambda_y or N (3/2+)
      I = IS(w, sqrt(6.0)/5.0, 3.5, 4.0);
      m_V1 = I*msigma*aaprime/2.0/a*(1.0/mQ-5.0/3.0/mq);
      m_V2 = I*msigma/a*(6.0*msigma/a-5.0*aaprime/6.0/mq+6.0*msigma2/aaprime/mQ-msigma/6.0/mq/mQ*(5*a-ap));
      m_V3 = -I*msigma/3.0/a/mQ*(aaprime+18.0*msigma2/aaprime);
      m_V4 = I*2.0*msigma*aaprime/3.0/mQ/a;
      m_A1 = -I*msigma/a*(6.0*msigma/a-aaprime/6.0/mQ+msigma/6.0/mq/mQ*(5.0*a-ap));
      m_A2 = I*aaprime/a*(6.0*msigma2/a/aaprime-5.0*msigma/6.0/mq-2.0*msigma/3.0/mQ+aaprime/72.0/mq/mQ*(5.0*a+ap));
      m_A3 = -I*aaprime/3.0/a/mQ*(msigma-msigma2/2.0/mq/aaprime*(5.0*a-ap)+aaprime/24.0/mq*(5.0*a+ap));
      m_A4 = -I*sqr(aaprime)/36.0/mq/mQ/a*(ap+5*a);
      break;
    case 10:
      // Omega_x (1/2+) -> Omega_y or Xi (1/2+ |200000>)
      //I = IS(w, 1.0/3.0, 3.0/2.0, 1.0);
      abort();
      break;
    case 11:
      // Omega_x (1/2+) -> Omega_y or Xi (1/2+1 |200010>)
      //I = IS(w, 1.0/sqrt(6.0), 2.5, 1.0);
      abort();
      break;
    case 12:
      // Omega_x (1/2+) -> Omega_y or Xi (1/2+2 |320002>)
      //I = IS(w, -sqrt(10.0)/3.0/sqrt(3.0), 3.5, 1.0);
      abort();
      break;
    case 13:
      // Omega_x (1/2+) -> Omega_y or Xi (1/2- |210001>)
      //I = IS(w, -1.0/3.0, 2.5, 1.0);
      abort();
      break;
    case 14:
      // Omega_x (1/2+) -> Omega_y or Xi (1/2-1 |310001>)
      //I = IS(w, -sqrt(2.0)/3.0, 2.5, 1.0);
      abort();
      break;
    case 15:
      // Lambda_x (1/2+) -> Lambda_y or N (1/2+1)
      //I = IS(w, sqrt(3.0/2.0), 1.5, 1.0);
      abort();
      break;
    default:
      THROW(fatal_error,"Mode not implemented in ST::CalcFFs.");
    }
    m_calced = true;
  }

  double ST::IS(double w, double factor, double exp1, double exp2)
  {
    double aaprime2 = sqr(0.5*(a+ap));
    return factor*pow(a*ap/aaprime2, exp1)/
      pow(1.0+1.5*sqr(msigma)*(sqr(w)-1.0)/aaprime2, exp2);
  }

}
}
