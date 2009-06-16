#include "ATOOLS/Org/CXXFLAGS_PACKAGES.H"
#ifdef USING__GOLEM95

#include "EXTRA_XS/NLO/Loop_ME_Base.H"
#include "HELICITIES/Loops/Box_4.H"
#include "HELICITIES/Main/Polarization_Tools.H"

using namespace EXTRAXS;
using namespace PHASIC;
using namespace ATOOLS;
using namespace HELICITIES;

namespace EXTRAXS {
  class PPPP_QED_Loop : public Loop_ME_Base {
    std::vector<std::vector<double> > m_smat;
    std::vector<ATOOLS::Vec4D> m_delta_a1;
    std::vector<ATOOLS::Vec4D> m_delta_a2;
    std::vector<ATOOLS::Vec4D> m_delta_a3;
    std::vector<ATOOLS::Vec4D> m_delta_a4;
    
    HELICITIES::Box_4 m_box4;

    void FillKinematics(const Vec4D_Vector& momenta);

  public:
    PPPP_QED_Loop(const Process_Info& pi, const Flavour_Vector& flavs) :
      Loop_ME_Base(pi, flavs),
      m_smat(4, std::vector<double>(4, 0.0)),
      m_delta_a1(4), m_delta_a2(4), m_delta_a3(4), m_delta_a4(4),
      m_box4(m_smat, m_delta_a1, m_delta_a2, m_delta_a3, m_delta_a4)
    {
    }

    ~PPPP_QED_Loop() {
    }

    void Calc(const ATOOLS::Vec4D_Vector& momenta);
  };
}


void PPPP_QED_Loop::Calc(const Vec4D_Vector& momenta) {
  DEBUG_FUNC(momenta.size());
  FillKinematics(momenta);
  m_box4.Calc();

  vector<Polarization_Vector> eps;
  for (size_t i=0; i<4; ++i) {
    eps.push_back(Polarization_Vector(momenta[i], m_flavs[i].Mass()));
  }
  vector<int> spins(4);
  DivArrC amp;
  for (spins[0]=0; spins[0]<2; ++spins[0]) {
    for (spins[1]=0; spins[1]<2; ++spins[1]) {
      for (spins[2]=0; spins[2]<2; ++spins[2]) {
        for (spins[3]=0; spins[3]<2; ++spins[3]) {
          Complex e01(eps[0][spins[0]]*eps[1][spins[1]]);
          Complex e02(eps[0][spins[0]]*eps[2][spins[2]]);
          Complex e03(eps[0][spins[0]]*eps[3][spins[3]]);
          Complex e12(eps[1][spins[1]]*eps[2][spins[2]]);
          Complex e13(eps[1][spins[1]]*eps[3][spins[3]]);
          Complex e23(eps[2][spins[2]]*eps[3][spins[3]]);
          Lorentz_Ten4<Complex> ten=
            BuildTensor(eps[0][spins[0]], eps[1][spins[1]],
                        eps[2][spins[2]], eps[3][spins[3]]);
          
          for (size_t i=0; i<amp.GetResult().size(); ++i) {
            amp[i]=Contraction(ten, 1, 2, 3, 4,
                               m_box4[i], 1, 2, 3, 4);
          }
          //amp*= DivArrC(trace(...));
          m_res.Insert(amp, spins);
        }
      }
    }
  }
  exit(0);
}

void PPPP_QED_Loop::FillKinematics(const Vec4D_Vector& momenta) {
  Vec4D_Vector r(4);
  r[0]=Vec4D(0.0, 0.0, 0.0, 0.0);
  r[1]=r[0]+momenta[1];
  r[2]=r[1]+momenta[2];
  r[3]=r[2]+momenta[3];
  
  for (size_t i=0; i<4; ++i) {
    m_delta_a1[i]=r[i]-r[0];
    m_delta_a2[i]=r[i]-r[1];
    m_delta_a3[i]=r[i]-r[2];
    m_delta_a4[i]=r[i]-r[3];
  }
  
  for (size_t i=0; i<4; ++i) {
    DEBUG_VAR(m_smat[i].size());
    for (size_t j=0; j<4; ++j) {
      m_smat[i][j]=(r[i]-r[j]).Abs2();
      if (IsZero(m_smat[i][j],1e-10)) m_smat[i][j]=0.0;
    }
  }
}


DECLARE_LOOPME_GETTER(PPPP_QED_Loop_Getter,"PPPP_QED_Loop")
Loop_ME_Base *PPPP_QED_Loop_Getter::operator()(const Process_Info &pi) const
{
  if (pi.m_loopgenerator!="Internal") return NULL;
  if (pi.m_fi.m_nloqcdtype!=nlo_type::lo) return NULL;
  if (pi.m_fi.m_nloewtype==nlo_type::loop) {
    Flavour_Vector fl=pi.ExtractFlavours();
    if (fl.size()!=4) return NULL;
    for (size_t i=0; i<4; ++i) if (!fl[i].IsPhoton()) return NULL;
    return new PPPP_QED_Loop(pi, fl);
  }
  return NULL;
}

#endif
