namespace HADRONS {
namespace VA_B_B_FFs {

class ST : public FormFactor_Base {
  int m_mode;
  double mq, mQ, msigma;
  double a, ap;
  double IS(double w, double factor, double exp1, double exp2);
public:
  ST(struct GeneralModel model, double* masses, std::string variant);
  void CalcFFs(ATOOLS::Vec4D p0, ATOOLS::Vec4D p1);
};

ST::ST(GeneralModel model, double* masses, std::string variant) : 
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

void ST::CalcFFs(Vec4D p0, Vec4D p1)
{
  Vec4D v0 = p0/m0;
  Vec4D v1 = p1/m1;
  double w = v0*v1;

  double aaprime = 0.5*(a+ap);

  double I;
  double amap(0.0);
  switch(m_mode) {
    case 0:
      // Lambda_x (1/2+) -> Lambda_y or N (1/2+)
      I = IS(w, 1.0, 1.5, 2.0);
      m_V1 = I*(1.0+msigma/aaprime*(ap/mq+a/mQ));
      m_V2 = -I*(msigma*ap/mq/aaprime-a*ap/6.0/mq/mQ);
      m_V3 = -I*msigma*a/mQ/aaprime;
      m_A1 = I*(1.0-a*ap/18.0/mq/mQ);
      m_A2 = -I*(msigma*ap/mq/aaprime+4.0*sqr(msigma)*a*ap/5.0/mq/mQ/sqr(aaprime)+a*ap/18.0/mq/mQ);
      m_A3 = I*(msigma*a/mQ/aaprime+4.0*sqr(msigma)*a*ap/5.0/mq/mQ/sqr(aaprime));
      break;
    case 1:
      // Lambda_x (1/2+) -> Lambda_y or N (1/2-)
      I = IS(w, sqrt(2.0), 2.5, 3.0);
      m_V1 = I*aaprime/12.0*(3.0/mq-1.0/mQ);
      m_V2 = -I*(2.0*msigma/a-aaprime/4.0/mq+2.0*sqr(msigma)/aaprime/mQ-msigma/12.0/mq/mQ*(5.0*a-3.0*ap));
      m_V3 = I*2.0*sqr(msigma)/mQ/aaprime;
      m_A1 = I*(2.0*msigma/a-aaprime/12.0/mQ+msigma/36.0/mq/mQ*(11.0*a-5.0*ap));
      m_A2 = -I*(2.0*msigma/a-aaprime/4.0/mq-aaprime/6.0/mQ+msigma/18.0/mq/mQ*(a-ap));
      m_A3 = I*aaprime/6.0/mQ*(1.0+msigma/2.0/mq/aaprime*(ap-3.0*a));
      break;
    case 2:
      // Lambda_x (1/2+) -> Lambda_y or N (1/2+1)
      I = IS(w, sqrt(3.0)/2.0, 2.5, 3.0);
      m_V1 = I/2.0/ap/a*
          (sqr(a)-sqr(ap)-2.0*msigma/3.0*(a/mQ*(5.0*ap-3.0*a)-ap/mq*(5.0*a-3.0*ap)));
      m_V2 = -I*(5.0*a-3.0*ap)/3.0/mq*(msigma/a-aaprime/3.0/mQ);
      m_V3 = I*msigma/6.0/mQ/ap*(5.0*ap-3.0*a);
      m_A1 = I*((sqr(a)-sqr(ap))/2.0/a/ap-aaprime/54.0/mq/mQ*(5.0*a-3.0*ap));
      m_A2 = -I*msigma/3.0/mq/a*
          ((5.0*a-3.0*ap)+4.0*msigma*a/mQ/aaprime*(a-ap)+aaprime/18.0/mQ*(5.0*a-ap));
      m_A3 = -I*msigma/3.0/mQ/ap*
          ((5.0*ap-3.0*a)-4.0*msigma*ap/mQ/aaprime*(a-ap));
      break;
    case 10:
      // Omega_x (1/2+) -> Omega_y or Xi (1/2+ |200000>)
      I = IS(w, 1.0/3.0, 1.5, 2.0);
      m_V1 = -I*(1.0+msigma/aaprime*(ap/mq+a/mQ));
      m_V2 = 2.0*I*(1.0-msigma/aaprime*(ap/2.0/mq-a/mQ));
      m_V3 = 2.0*I*(1.0+msigma/aaprime*(ap/mq-a/2.0/mQ));
      m_A1 = -I;
      m_A2 = I*msigma*ap/mq/aaprime;
      m_A3 = -I*msigma*a/mQ/aaprime;
      break;
    case 11: {
      // Omega_x (1/2+) -> Omega_y or Xi (1/2+1 |200010>)
      I = IS(w, 1.0/2.0/sqrt(3.0), 2.5, 3.0);
      amap = sqr(a)-sqr(ap);
      m_V1 = -I/2.0/a/ap*(amap+2.0*msigma/3.0*(ap*(5.0*a-3.0*ap)/mq+a*(3.0*a-5.0*ap)/mQ));
      m_V2 = I/a/ap*(amap-2.0*msigma/3.0*(ap*(5.0*a-3.0*ap)/2.0/mq-a*(3.0*a-5.0*ap)/mQ));
      m_V3 = I/a/ap*(amap+2.0*msigma/3.0*(ap*(5.0*a-3.0*ap)/mq-a*(3.0*a-5.0*ap)/2.0/mQ));
      m_A1 = -I*amap/2.0/a/ap;
      m_A2 = I*msigma/3.0/mq*(5.0*a-3.0*ap)/a;
      m_A3 = -I*msigma/3.0/mQ*(3.0*a-5.0*ap)/ap;
      break;
	}
    case 12:
      // Omega_x (1/2+) -> Omega_y or Xi (1/2+2 |320002>)
      I = IS(w, -4.0/9.0*sqrt(5.0), 3.5, 4.0);
      m_V1 = -I*msigma*aaprime/2.0/a*(1.0/mq-1.0/mQ);
      m_V2 = -0.5*m_V1;
      m_V3 = m_V2;
      m_A1 = 0.0;
      m_A2 = I*msigma/a*(27.0*msigma/5.0/a-aaprime/4.0*(4.0/mq+3.0/mQ));
      m_A3 = -I*msigma/a*(27.0*msigma/5.0/a+aaprime/4.0/mQ);
      break;
    case 13:
      // Omega_x (1/2+) -> Omega_y or Xi (1/2- |210001>)
      I = IS(w, -sqrt(2.0)/3.0, 2.5, 3.0);
      m_V1 = -I*aaprime/12.0*(1.0/mq+5.0/mQ);
      m_V2 = I*(2.0*msigma/a*(1.0-msigma/aaprime*(2.0*ap/mq-a/mQ))-aaprime/12.0/mq);
      m_V3 = I*4.0*msigma/a*(1.0+msigma/aaprime*(ap/mq-a/2.0/mQ));
      m_A1 = I*(2.0*msigma/a-5.0*aaprime/12.0/mQ);
      m_A2 = -I*(2.0*msigma/a-aaprime*(1.0/4.0/mq+1.0/3.0/mQ));
      m_A3 = I*aaprime/3.0/mQ;
      break;
    case 14:
      // Omega_x (1/2+) -> Omega_y or Xi (1/2-1 |310001>)
      I = IS(w, -2.0/3.0, 2.5, 3.0);
      m_V1 = -I*aaprime/6.0*(1.0/mq-1.0/mQ);
      m_V2 = I*(msigma/a*(1.0+msigma/aaprime*(ap/mq+a/mQ))-aaprime/6.0/mq);
      m_V3 = -I*msigma/a*(1.0+msigma/aaprime*(ap/mq+a/mQ));
      m_A1 = -I*(2.0*msigma/a-aaprime/6.0/mQ);
      m_A2 = -I*(msigma/a-aaprime/6.0*(3.0/mq+1.0/mQ));
      m_A3 = I*(3.0*msigma/a+aaprime/6.0/mQ);
      break;
    default:
      THROW(fatal_error,"Mode not implemented in HO::CalcFFs.");
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
