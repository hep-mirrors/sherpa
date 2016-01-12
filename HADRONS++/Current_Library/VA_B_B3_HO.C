namespace HADRONS {
namespace VA_B_B3_FFs {

  class HO : public FormFactor_Base {
    int m_mode;
    double mq, mQ, msigma;
    double a, ap;
    double IH(double w, double factor, double exponent);
  public:
    HO(struct GeneralModel model, double* masses, string variant);
    void CalcFFs(ATOOLS::Vec4D p0, ATOOLS::Vec4D p1);
  };

  HO::HO(GeneralModel model, double* masses, string variant) : 
    FormFactor_Base(model, masses)
  {
    m_mode = int(model("mode",0.0)+0.5);
    string prefix=variant+"_";
    mq = model(prefix+"mq",1.0);
    mQ = model(prefix+"mQ",1.0);
    msigma = model(prefix+"msigma",1.0);
    a = model(prefix+"alpha",0.5);
    ap = model(prefix+"alpha_prime",0.5);
  }

  void HO::CalcFFs(Vec4D p0, Vec4D p1)
  {
    Vec4D v0 = p0/m0;
    Vec4D v1 = p1/m1;
    double w = v0*v1;

    double a2 = sqr(a);
    double ap2 = sqr(ap);
    double aaprime2 = 0.5*(a2+ap2);
    double msigma2 = sqr(msigma);

    double I;
    switch(m_mode) {
    case 0:
      // Lambda_x (1/2+) -> Lambda_y or N (3/2-)
      I = IH(w, -1.0/sqrt(3.0), 2.5);
      m_V1 = I*3.0*msigma/a*(1.0+msigma/aaprime2*(ap2/mq+a2/mQ));
      m_V2 = -I*(3.0*sqr(msigma)/mq*ap2/aaprime2/a -
             5.0*a*ap2*msigma/(4.0*aaprime2*mq*mQ));
      m_V3 = -I*(3.0*sqr(msigma)/mQ*a/aaprime2 + a/(2.0*mQ));
      m_V4 = I*a/mQ;
      m_A1 = I*(3.0*msigma/a-a/(2.0*mQ)*(1.0+3.0*msigma*ap2/(2.0*mq*aaprime2)));
      m_A2 = -I*(3.0*sqr(msigma)/mq*ap2/aaprime2/a +
             msigma*a*ap2/(4.0*mq*mQ*sqr(aaprime2))*(aaprime2+12.0*sqr(msigma)));
      m_A3 = I*a/mQ/aaprime2*(aaprime2/2.0+3.0*sqr(msigma) +
             ap2*msigma/mq/aaprime2*(aaprime2+6.0*sqr(msigma)));
      m_A4 = -I*(a/mQ+msigma/mq/mQ*ap2*a/aaprime2);
      break;
    case 1:
      // Lambda_x (1/2+) -> Lambda_y or N (3/2+)
      I = IH(w, 1.0/sqrt(5.0), 3.5);
      m_V1 = -I*msigma/2.0*(5.0/mq-3.0/mQ);
      m_V2 = I*msigma/a*(6.0*msigma/a-5.0*a/2.0/mq+6.0*msigma2*a/aaprime2/mQ-msigma*a/2.0/aaprime2/mq/mQ*(a2-2.0*ap2));
      m_V3 = -I*msigma/mQ*(1.0+6.0*msigma2/aaprime2);
      m_V4 = I*2.0*msigma/mQ;
      m_A1 = -I*(6.0*msigma2/a2-msigma/2.0/mQ+msigma2/6.0/aaprime2/mq/mQ*(11.0*a2-6.0*ap2));
      m_A2 = I*(6.0*msigma2/a2-5.0*msigma/2.0/mq-2.0*msigma/mQ+5.0*a2/12.0/mq/mQ-2.0*msigma2/a2/3.0/aaprime2/mq/mQ);
      m_A3 = -I*(msigma/2.0/mQ-5.0*a2/24.0/mq/mQ-msigma2/4.0/mq/mQ/aaprime2*(5.0*a2-2.0*ap2));
      m_A4 = -I*5.0*a2/6.0/mq/mQ;
      break;
    case 10:
      // Omega_x (1/2+) -> Omega_y or Xi (1/2+ |200000>)
      //I = IH(w, 1.0/3.0, 3.0/2.0);
      abort();
      break;
    case 11:
      // Omega_x (1/2+) -> Omega_y or Xi (1/2+1 |200010>)
      //I = IH(w, 1.0/sqrt(6.0), 2.5);
      abort();
      break;
    case 12:
      // Omega_x (1/2+) -> Omega_y or Xi (1/2+2 |320002>)
      //I = IH(w, -sqrt(10.0)/3.0/sqrt(3.0), 3.5);
      abort();
      break;
    case 13:
      // Omega_x (1/2+) -> Omega_y or Xi (1/2- |210001>)
      //I = IH(w, -1.0/3.0, 2.5);
      abort();
      break;
    case 14:
      // Omega_x (1/2+) -> Omega_y or Xi (1/2-1 |310001>)
      //I = IH(w, -sqrt(2.0)/3.0, 2.5);
      abort();
      break;
    case 15:
      // Lambda_x (1/2+) -> Lambda_y or N (1/2+1)
      //I = IH(w, sqrt(3.0/2.0), 1.5);
      abort();
      break;
    default:
      THROW(fatal_error,"Mode not implemented in HO::CalcFFs.");
    }
    m_calced = true;
  }

  double HO::IH(double w, double factor, double exponent)
  {
    double aaprime2 = 0.5*(sqr(a)+sqr(ap));
    return factor*pow(a*ap/aaprime2, exponent)*
      exp(-3.0*sqr(msigma)*(sqr(w)-1.0)/(2.0*aaprime2));
  }

}
}
