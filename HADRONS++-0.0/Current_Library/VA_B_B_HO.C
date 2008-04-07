namespace HADRONS {
namespace VA_B_B_FFs {

class HO : public FormFactor_Base {
  int m_mode;
  double mq, mQ, msigma;
  double a, ap;
  double IH(double w, double factor, double exponent);
public:
  HO(struct GeneralModel model, double* masses, std::string variant);
  void CalcFFs(ATOOLS::Vec4D p0, ATOOLS::Vec4D p1);
};

HO::HO(GeneralModel model, double* masses, std::string variant) : 
  FormFactor_Base(model, masses)
{
  m_mode = int(model("mode",0.0)+0.5);
  std::string prefix=variant+"_";
  mq = model(prefix+"mq",1.0);
  mQ = model(prefix+"mQ",1.0);
  msigma = model(prefix+"msigma",1.0);
  a = model(prefix+"alpha",1.0);
  ap = model(prefix+"alpha_prime",1.0);
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
      // Lambda_x (1/2+) -> Lambda_y or N (1/2+)
      I = IH(w, 1.0, 1.5);
      m_V1 = I*(1.0 + msigma/aaprime2*(ap2/mq+a2/mQ));
      m_V2 = -I*(msigma/mq*ap2/aaprime2 - a2*ap2/(4.0*aaprime2*mq*mQ));
      m_V3 = -I*(msigma/mQ*a2/aaprime2);
      m_A1 = I*(1.0-a2*ap2/(12.0*aaprime2*mq*mQ));
      m_A2 = -I*(msigma/mq*ap2/aaprime2+
          a2*ap2/(12.0*mq*mQ*aaprime2)*(1.0+12.0*msigma2/aaprime2));
      m_A3 = I*(msigma/mQ*a2/aaprime2+msigma2*a2*ap2/(mq*mQ*sqr(aaprime2)));
      break;
    case 1:
      // Lambda_x (1/2+) -> Lambda_y or N (1/2-)
      I = IH(w, 1.0, 2.5);
      m_V1 = I*a/6.0*(3.0/mq-1.0/mQ);
      m_V2 = -I*(2.0*msigma/a-a/2.0/mq+2.0*msigma2*a/mQ/aaprime2-msigma*a*(3.0*a2-2.0*ap2));
      m_V3 = I*2.0*msigma2*a/(mQ*aaprime2);
      m_A1 = I*(2.0*msigma/a-a/6.0/mQ+msigma*a/(6.0*mq*mQ*aaprime2)*(3.0*a2-2.0*ap2));
      m_A2 = I*(-2.0*msigma/a+a/2.0/mq+a/3.0/mQ);
      m_A3 = I*a/3.0/mQ*(1.0-msigma/(2.0*mq*aaprime2)*(3.0*a2-2.0*ap2));
      break;
    case 2:
      // Lambda_x (1/2+) -> Lambda_y or N (1/2+1)
      I = IH(w, sqrt(3.0/2.0), 1.5);
      m_V1 = I/2.0/aaprime2*
          (a2-ap2-msigma/3.0/aaprime2*(ap2/mq*(7.0*a2-3.0*ap2)+a2/mQ*(7.0*ap2-3.0*a2)));
      m_V2 = -I*ap2/6.0/mq/sqr(aaprime2)*(7.0*a2-3.0*ap2)*(msigma-a2/4.0/mQ);
      m_V3 = I*a2*msigma/6.0/mQ/sqr(aaprime2)*(7.0*ap2-3.0*a2);
      m_A1 = I*((a2-ap2)/2.0/aaprime2-a2*ap2/72.0/sqr(aaprime2)/mq/mQ*(7.0*a2-3.0*ap2));
      m_A2 = -I*ap2/6.0/mq/sqr(aaprime2)*
          ((7.0*a2-3.0*ap2)*(msigma+a2/6.0/mQ)+7.0*msigma2*a2/mQ/aaprime2*(a2-ap2));
      m_A3 = -I*a2*msigma/6.0/mQ/sqr(aaprime2)*
          ((7.0*ap2-3.0*a2)-7.0*msigma*ap2/mq/aaprime2*(a2-ap2));
      break;
    case 10:
      // Omega_x (1/2+) -> Omega_y or Xi (1/2+ |200000>)
      I = IH(w, 1.0/3.0, 3.0/2.0);
      m_V1 = -I*(1.0+msigma/aaprime2*(ap2/mq+a2/mQ));
      m_V2 = 2.0*I*(1.0-msigma/aaprime2*(ap2/2.0/mq-a2/mQ));
      m_V3 = 2.0*I*(1.0+msigma/aaprime2*(ap2/mq-a2/2.0/mQ));
      m_A1 = -I;
      m_A2 = I*msigma*ap2/mq/aaprime2;
      m_A3 = -I*msigma*a2/mQ/aaprime2;
      break;
    case 11:
      // Omega_x (1/2+) -> Omega_y or Xi (1/2+1 |200010>)
      I = IH(w, 1.0/sqrt(6.0), 2.5);
      m_V1 = -I/2.0/a/ap*((a2-ap2)+msigma/3.0/aaprime2*(ap2*(7.0*a2-3.0*ap2)/mq+a2*(3.0*a2-7.0*ap2)/mQ));
      m_V2 = I/a/ap*((a2-ap2)-msigma/3.0/aaprime2*(ap2*(7.0*a2-3.0*ap2)/2.0/mq-a2*(3.0*a2-7.0*ap2)/mQ));
      m_V3 = I/a/ap*((a2-ap2)+msigma/3.0/aaprime2*(ap2*(7.0*a2-3.0*ap2)/mq-a2*(3.0*a2-7.0*ap2)/2.0/mQ));
      m_A1 = -I*(a2-ap2)/2.0/a/ap;
      m_A2 = I*msigma*ap*(7.0*a2-3.0*ap2)/6.0/mq/a/aaprime2;
      m_A3 = -I*msigma*a*(3.0*a2-7.0*ap2)/6.0/mQ/ap/aaprime2;
      break;
    case 12:
      // Omega_x (1/2+) -> Omega_y or Xi (1/2+2 |320002>)
      I = IH(w, -sqrt(10.0)/3.0/sqrt(3.0), 3.5);
      m_V1 = -I*msigma*(1.0/mq-1.0/mQ);
      m_V2 = -0.5*m_V1;
      m_V3 = m_V2;
      m_A1 = 0.0;
      m_A2 = I*msigma/a*(18.0*msigma/5.0/a-a/2.0*(4.0/mq+3.0/mQ));
      m_A3 = -I*msigma/a*(18.0*msigma/5.0/a+a/2.0/mQ);
      break;
    case 13:
      // Omega_x (1/2+) -> Omega_y or Xi (1/2- |210001>)
      I = IH(w, -1.0/3.0, 2.5);
      m_V1 = -I*a/6.0*(1.0/mq+5.0/mQ);
      m_V2 = I*(2.0*msigma/a*(1.0-msigma/aaprime2*(2.0*ap2/mq-a2/mQ))-a/6.0/mq);
      m_V3 = I*4.0*msigma/a*(1.0+msigma/aaprime2*(ap2/mq-a2/2.0/mQ));
      m_A1 = I*(2.0*msigma/a-5.0*a/6.0/mQ);
      m_A2 = -I*(2.0*msigma/a-a*(1.0/2.0/mq+2.0/3.0/mQ));
      m_A3 = I*2.0*a/3.0/mQ;
      break;
    case 14:
      // Omega_x (1/2+) -> Omega_y or Xi (1/2-1 |310001>)
      I = IH(w, -sqrt(2.0)/3.0, 2.5);
      m_V1 = -I*a/3.0*(1.0/mq-1.0/mQ);
      m_V2 = I*(msigma/a*(1.0+msigma/aaprime2*(ap2/mq+a2/mQ))-a/3.0/mQ);
      m_V3 = -I*msigma/a*(1.0+msigma/aaprime2*(ap2/mq+a2/mQ));
      m_A1 = -I*(2.0*msigma/a-a/3.0/mQ);
      m_A2 = -I*(msigma/a-a/3.0*(3.0/mq+1.0/mQ));
      m_A3 = I*(3.0*msigma/a+a/3.0/mQ);
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
