#include "XYZFuncs.H"
#include "Exception.H"

#include "Current_Base.H"

using namespace HADRONS;
using namespace ATOOLS;
using namespace HELICITIES;
using namespace std;

namespace HADRONS {

  class VA_B_B_hepph9409272 : public Current_Base
  {
    double m_V_CKM;
    double m_v, m_a;
    double m0, m1;
    double m_V1, m_V2, m_V3, m_A1, m_A2, m_A3;
    double m_V1_0, m_V2_0, m_V3_0, m_A1_0, m_A2_0, m_A3_0;
    double m_V1_Lambda12, m_V2_Lambda12, m_V3_Lambda12, m_A1_Lambda12,
      m_A2_Lambda12, m_A3_Lambda12;
    double m_V1_Lambda24, m_V2_Lambda24, m_V3_Lambda24, m_A1_Lambda24,
      m_A2_Lambda24, m_A3_Lambda24;
    
  public:
    VA_B_B_hepph9409272(ATOOLS::Flavour* flavs, int n, int* indices,
                        std::string name) :
      Current_Base(flavs, n, indices, name) {};
    void SetModelParameters( struct GeneralModel _md );
    void Calc(const ATOOLS::Vec4D * moms);
    double Fit(double q2, double f0, double Lambda12, double Lambda24);
  };
}

void VA_B_B_hepph9409272::SetModelParameters( struct GeneralModel model )
{
  double V_CKM(1.0),
    V1_0(0.0), V1_Lambda1(1.0), V1_Lambda2(1.0),
    A1_0(0.0), A1_Lambda1(1.0), A1_Lambda2(1.0),
    V2_0_M(0.0), A2_0_M(0.0);

  if ((m_flavs[0].Kfcode()==kf_Sigma_plus||m_flavs[0].Kfcode()==kf_Sigma_minus)
      && m_flavs[1].Kfcode()==kf_Lambda) {
    V_CKM = Tools::Vud;
    V1_0 = 0.0;
    V1_Lambda1 = -0.32;
    V1_Lambda2 = -1.72;
    A1_0 = 0.60;
    A1_Lambda1 = 0.77;
    A1_Lambda2 = 1.05;
    V2_0_M = 1.04;
    A2_0_M = 0.0;
  }
  else if (m_flavs[0].Kfcode()==kf_Lambda && m_flavs[1].Kfcode()==kf_p_plus) {
    V_CKM = Tools::Vus;
    V1_0 = -1.19;
    V1_Lambda1 = 0.71;
    V1_Lambda2 = 0.98;
    A1_0 = -0.99;
    A1_Lambda1 = 0.81;
    A1_Lambda2 = 1.12;
    V2_0_M = -0.85;
    A2_0_M = -0.025;
  }
  else if (m_flavs[0].Kfcode()==kf_Sigma_minus && m_flavs[1].Kfcode()==kf_n) {
    V_CKM = Tools::Vus;
    V1_0 = -0.97;
    V1_Lambda1 = 0.64;
    V1_Lambda2 = 0.90;
    A1_0 = 0.27;
    A1_Lambda1 = 0.83;
    A1_Lambda2 = 1.16;
    V2_0_M = 0.62;
    A2_0_M = 0.0061;
  }
  else if (m_flavs[0].Kfcode()==kf_Xi_minus && m_flavs[1].Kfcode()==kf_Lambda) {
    V_CKM = Tools::Vus;
    V1_0 = 1.19;
    V1_Lambda1 = 0.68;
    V1_Lambda2 = 0.89;
    A1_0 = 0.33;
    A1_Lambda1 = 0.81;
    A1_Lambda2 = 1.10;
    V2_0_M = 0.07;
    A2_0_M = 0.0076;
  }
  else if (m_flavs[0].Kfcode()==kf_Xi_minus && m_flavs[1].Kfcode()==kf_Sigma) {
    V_CKM = Tools::Vus;
    V1_0 = 0.69;
    V1_Lambda1 = 0.75;
    V1_Lambda2 = 1.05;
    A1_0 = 0.94;
    A1_Lambda1 = 0.81;
    A1_Lambda2 = 1.12;
    V2_0_M = 0.98;
    A2_0_M = 0.022;
  }
  else if (m_flavs[0].Kfcode()==kf_Xi && m_flavs[1].Kfcode()==kf_Sigma_plus) {
    V_CKM = Tools::Vus;
    V1_0 = 0.98;
    V1_Lambda1 = 0.75;
    V1_Lambda2 = 1.05;
    A1_0 = 1.33;
    A1_Lambda1 = 0.81;
    A1_Lambda2 = 1.12;
    V2_0_M = 1.38;
    A2_0_M = 0.0306;
  }

  m_V_CKM    = model("V_CKM",V_CKM);
  m_v        = model("v", 1.0);
  m_a        = model("a", -1.0);

  m_V1_0     = model("V1_0",V1_0);
  m_V1_Lambda12     = sqr(model("V1_Lambda1",V1_Lambda1));
  m_V1_Lambda24     = sqr(sqr(model("V1_Lambda2",V1_Lambda2)));

  m_A1_0     = model("A1_0",A1_0);
  m_A1_Lambda12     = sqr(model("A1_Lambda1",A1_Lambda1));
  m_A1_Lambda24     = sqr(sqr(model("A1_Lambda2",A1_Lambda2)));

  m_V2_0     = model("V2_0",p_masses[0]*V2_0_M);
  m_V2_Lambda12     = sqr(model("V2_Lambda1",1.0));
  m_V2_Lambda24     = sqr(sqr(model("V2_Lambda2",1.0)));

  m_A2_0     = model("A2_0",p_masses[0]*A2_0_M);
  m_A2_Lambda12     = sqr(model("A2_Lambda1",1.0));
  m_A2_Lambda24     = sqr(sqr(model("A2_Lambda2",1.0)));

  m_V3_0     = model("V3_0",0.0);
  m_V3_Lambda12     = sqr(model("V3_Lambda1",1.0));
  m_V3_Lambda24     = sqr(sqr(model("V3_Lambda2",1.0)));

  m_A3_0     = model("A3_0",0.0);
  m_A3_Lambda12     = sqr(model("A3_Lambda1",1.0));
  m_A3_Lambda24     = sqr(sqr(model("A3_Lambda2",1.0)));
}

void VA_B_B_hepph9409272::Calc(const ATOOLS::Vec4D * moms)
{
  double q2=(moms[p_i[0]]-moms[p_i[1]]).Abs2();

  double V1=Fit(q2,m_V1_0,m_V1_Lambda12,m_V1_Lambda24);
  double V2=Fit(q2,m_V2_0,m_V2_Lambda12,m_V2_Lambda24);
  double V3=Fit(q2,m_V3_0,m_V3_Lambda12,m_V3_Lambda24);
  double A1=Fit(q2,m_A1_0,m_A1_Lambda12,m_A1_Lambda24);
  double A2=Fit(q2,m_A2_0,m_A2_Lambda12,m_A2_Lambda24);
  double A3=Fit(q2,m_A3_0,m_A3_Lambda12,m_A3_Lambda24);

  Complex v1=m_v*(V1+(p_masses[0]+p_masses[1])*V2/p_masses[0]);
  Complex a1=m_a*(A1+(p_masses[0]+p_masses[1])*A2/p_masses[0]);

  Complex v2=m_v*(-V2+V3)/p_masses[0];
  Complex a2=m_a*(-A2+A3)/p_masses[0];

  Complex v3=m_v*(-V2-V3)/p_masses[0];
  Complex a3=m_a*(-A2-A3)/p_masses[0];

  XYZFunc F(2, moms, m_flavs, Tools::k0, m_anti, p_i);
  for(int h0=0; h0<2; h0++) {
    for(int h1=0; h1<2; h1++) {
      // 1 is "the barred spinor" in the current, 0 is the not-barred one
      // (other way around than VA_F_F, maybe unify?)
      Vec4C amp(Vec4D(0.0, 0.0, 0.0, 0.0));
      amp += F.L(1,h1, 0,h0, v1+a1, v1-a1);
      amp += moms[p_i[0]]*F.Y(1,h1, 0,h0, v2+a2, v2-a2);
      amp += moms[p_i[1]]*F.Y(1,h1, 0,h0, v3+a3, v3-a3);
      vector<pair<int,int> > spins;
      spins.push_back(make_pair(0,h0));
      spins.push_back(make_pair(1,h1));
      Insert( m_V_CKM*amp ,spins );
    }
  }
}

double VA_B_B_hepph9409272::Fit(double q2, double f0,
                                double Lambda12, double Lambda24)
{
  return f0/(1-q2/Lambda12+sqr(q2)/Lambda24);
}

DEFINE_CURRENT_GETTER(VA_B_B_hepph9409272,VA_B_B_hepph9409272_Getter,"VA_B_B_hepph9409272")

void VA_B_B_hepph9409272_Getter::PrintInfo(std::ostream &st,const size_t width) const {
  st<<std::endl;
}
