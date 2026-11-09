#include "METOOLS/HadronCurrents/VA_N_NPi.H"
#include "METOOLS/HadronCurrents/FormFactors/FF_N_N.H"
#include "METOOLS/HadronCurrents/FormFactors/FF_N_NPi.H"
#include "METOOLS/Main/XYZFuncs.H"
#include "MODEL/Main/Model_Base.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Message.H"

using namespace METOOLS;
using namespace ATOOLS;
using namespace std;

VA_N_NPi::VA_N_NPi(const ATOOLS::Flavour_Vector& flavs,
                   const std::vector<int>& indices,
                   const std::string& name) :
  Current_Base(flavs,indices,name),
  p_F1V(NULL), p_gpiNN(NULL), p_FPiPole(NULL),
  m_boson(kf_photon), m_pipolearg(0),
  m_mN(0.), m_norm(1.), m_s2(0.), m_gA(0.), m_fpi(0.),
  m_mNP(0.), m_mCNP(0.), m_mex2(0.), m_mpole2(0.),
  m_cNP(0.), m_cCNP(0.), m_cPF(0.), m_cX(0.), m_cCT(0.),
  m_onPF(false), m_wardcheck(true)
{
  for (size_t k(0);k<4;++k) p_NP[k] = p_CNP[k] = NULL;
  if (p_i.size()!=3)
    THROW(fatal_error,"Current "+m_name+" needs three legs [N_out, N_in, "
          "pion], got "+ToString(p_i.size())+".");
  for (size_t k(0);k<2;++k) {
    const kf_code kf(m_flavs[p_i[k]].Kfcode());
    if (kf!=kf_p_plus && kf!=kf_n)
      THROW(fatal_error,"Current "+m_name+" needs nucleon legs 0 and 1, got "
            +ToString(m_flavs[p_i[k]])+".");
  }
  const kf_code kfpi(m_flavs[p_i[2]].Kfcode());
  if (kfpi!=kf_pi && kfpi!=kf_pi_plus)
    THROW(fatal_error,"Current "+m_name+" needs a pion as leg 2, got "
          +ToString(m_flavs[p_i[2]])+".");
}

VA_N_NPi::~VA_N_NPi()
{
  for (size_t k(0);k<4;++k) {
    if (p_NP[k])  delete p_NP[k];
    if (p_CNP[k]) delete p_CNP[k];
  }
  if (p_F1V) delete p_F1V;
  if (p_gpiNN) delete p_gpiNN;
  if (p_FPiPole) delete p_FPiPole;
}

void VA_N_NPi::CheckTopology() const
{
  if (NIncomingLegs()!=1 || !IsIncoming(p_i[1]))
    THROW(fatal_error,"Current "+m_name+" is the N -> N' pi current: Indices "
          "must be [outgoing nucleon, incoming nucleon, pion], got ["
          +ToString(p_i[0])+", "+ToString(p_i[1])+", "+ToString(p_i[2])+"].");
}

double VA_N_NPi::CPi(const int Qin,const int Qout) const
{
  // pi-N-N isospin factor, N(Qin) -> N(Qout) + pi(Qin-Qout)
  if (Qin==Qout) return Qin==1 ? 1. : -1.;
  return sqrt(2.);
}

void VA_N_NPi::MakeVertexFF(FF_N_N * ff[4],const int iout,const int iin,
                            const int vff,const int aff,GeneralModel & md)
{
  static const std::string names[4] = { "F1", "F2", "FA", "FP" };
  const std::vector<int> idx = { iout, iin };
  for (size_t k(0);k<4;++k) {
    std::map<std::string,double> pmap;
    FF_Parameters params(k<2 ? vff : aff,m_ffflavs,idx,pmap,names[k],&md);
    ff[k] = new FF_N_N(params);
  }
}

void VA_N_NPi::SetModelParameters(struct GeneralModel model)
{
  // Exchange_ME calls this a second time with only the default v and a.
  if (p_gpiNN) return;

  const int QN4(m_flavs[p_i[0]].Charge()); 
  const int QN2(m_flavs[p_i[1]].Charge()); 
  const int Qpi(m_flavs[p_i[2]].Charge()); 
  const int QX(QN4+Qpi-QN2); // charge of the exchanged boson calculated from hadron current
  if (std::abs(QX)>1)
    THROW(fatal_error,"Current "+m_name+": no single boson connects these legs.");
  m_boson = static_cast<kf_code>(std::abs(int(model("BOSON",QX!=0 ? 24. : 22.))));
  if (m_boson!=kf_photon && m_boson!=kf_Z && m_boson!=kf_Wplus)
    THROW(fatal_error,"Current "+m_name+": BOSON must be 22, 23 or 24.");
  if ((QX!=0)!=(m_boson==kf_Wplus))
    THROW(fatal_error,"Current "+m_name+": BOSON "+ToString(m_boson)+
          " does not match the charges of the legs.");

  m_mN   = 0.5*(p_masses[0]+p_masses[1]);
  m_s2   = std::abs(MODEL::s_model->ComplexConstant(string("csin2_thetaW")));
  m_norm = (m_boson==kf_Wplus) ?
    std::abs(MODEL::s_model->ComplexConstant(string("CKM_0_0"))) : 1.;
  m_gA   = model("GA",1.2754);
  m_fpi  = model("fpi",0.093);
  m_pipolearg = int(model("PIPOLE_FF_ARG",0.)+0.5);

  // Internal nucleons are appended, so FF_N_N::FixMode sees the vertex.
  m_ffflavs = m_flavs;
  m_ffflavs.push_back(Flavour(kf_p_plus));
  m_ffflavs.push_back(Flavour(kf_n));
  const int ip(m_flavs.size()), in(ip+1);
  auto nuc = [ip,in](const int Q) { return Q==1 ? ip : in; };

  GeneralModel mX(model), mW(model);
  mX["BOSON"] = m_boson;
  mW["BOSON"] = 24.;
  const int vff(int(model("VECTOR_FF",10.)+0.5));
  const int aff(m_boson==kf_photon ? ff_zero : int(model("AXIAL_FF",20.)+0.5));

  // NP: X absorbed by N(Q2) -> N(Q2+QX), then pi emitted.
  const int QR(QN2+QX); // charge of resonance 
  if (QR==0 || QR==1) {
    MakeVertexFF(p_NP,nuc(QR),p_i[1],vff,aff,mX);
    m_mNP = m_ffflavs[nuc(QR)].Mass(true);
    m_cNP = CPi(QR,QN4);
  }
  // CNP: pi emitted from N(Q2) -> N(Q2-Qpi), then X absorbed.
  const int QRC(QN2 - Qpi);
  if (QRC==0 || QRC==1) {
    MakeVertexFF(p_CNP,p_i[0],nuc(QRC),vff,aff,mX);
    m_mCNP = m_ffflavs[nuc(QRC)].Mass(true);
    m_cCNP = CPi(QN2,QRC);
  }
  // PF: exchanged pion of charge Q2-Q4; no X pi0 pi0 vertex.
  const int Qex(QN2-QN4); 
  m_cPF  = CPi(QN2,QN4);
  m_mex2 = sqr(Flavour(Qex==0 ? kf_pi : kf_pi_plus).HadMass());
  if      (m_boson==kf_photon) m_cX = Qpi;
  else if (m_boson==kf_Z)      m_cX = Qpi*(1.-2.*m_s2);
  else                         m_cX = sqrt(2.)*(Qex!=0 ? 1. : -1.);
  m_onPF = (m_cX!=0.);
  // CT: Kroll-Ruderman (+ axial and pion-pole terms for Z, W).
  if (m_boson==kf_Wplus)
    m_cCT = Qpi==0 ? 1. : ((QN2==0)==(QX>0) ? 1. : -1.)/sqrt(2.);
  else m_cCT = Qpi/sqrt(2.);
  m_mpole2 = sqr(Flavour(m_boson==kf_Wplus ? kf_pi_plus : kf_pi).HadMass());

  if (!p_NP[0] && !p_CNP[0] && !m_onPF && m_cCT==0.)
    THROW(fatal_error,"Current "+m_name+": no diagram is charge allowed.");

  // F1V for PF and CT, pi-N-N and F_PiPole couplings: built on the appended
  // p, n pair with BOSON 24, so that FF_N_N gives F^p - F^n and bare F_A.
  const std::vector<int> pn = { ip, in };
  std::map<std::string,double> pmap;
  FF_Parameters pv(vff,m_ffflavs,pn,pmap,"F1",&mW); // params input
  p_F1V = new FF_N_N(pv);
  FF_Parameters pg(int(model("PION_NN_FF",1.)+0.5),m_ffflavs,pn,pmap,"FpiNN",&mW);
  p_gpiNN = new FF_N_NPi(pg);
  FF_Parameters pr(int(model("PIPOLE_FF",11.)+0.5),m_ffflavs,pn,pmap,"FPiPole",&mW);
  p_FPiPole = new FF_N_NPi(pr);

  msg_Tracking()<<"### "<<m_name<<": "<<m_flavs[p_i[1]]<<" -> "
                <<m_flavs[p_i[0]]<<" "<<m_flavs[p_i[2]]<<" via "
                <<Flavour(m_boson)<<", diagrams NP "<<(p_NP[0]!=NULL)
                <<", CNP "<<(p_CNP[0]!=NULL)<<", PF "<<m_onPF
                <<", CT "<<(m_cCT!=0.)<<", g_piNN = "<<p_gpiNN->FF_Assigner(0.)
                <<", norm = "<<m_norm<<"\n";
}

void VA_N_NPi::Calc(const ATOOLS::Vec4D_Vector& moms, bool anti)
{
  if (anti)
    THROW(not_implemented,"Current "+m_name+" has no anti-nucleon version.");
  const Vec4D  p4(moms[p_i[0]]), p2(moms[p_i[1]]), p5(moms[p_i[2]]);
  const Vec4D  q(p4+p5-p2), p(p4+p5), pc(p2-p5), K(p2-p4+p5);
  const double m4(p_masses[0]), m2(p_masses[1]), mpi2(p5.Abs2());
  const double Q2(-q.Abs2()), t((p2-p4).Abs2()), m2N(2.*m_mN);
  const double p4p5(p4*p5), p2p5(p2*p5);

  // NP: A + B p5slash, G D = c_pi g_piNN/(2m_N)/(p^2-m^2)
  double A(0.), B(0.), GD(0.), F1(0.), F2(0.), FA(0.), FP(0.);
  if (p_NP[0]) {
    A  = p.Abs2()-sqr(m4);  B = m4+m_mNP;
    GD = m_cNP*p_gpiNN->FF_Assigner(-mpi2)/m2N/(p.Abs2()-sqr(m_mNP));
    F1 = p_NP[0]->FF_Assigner(Q2); F2 = p_NP[1]->FF_Assigner(Q2);
    FA = p_NP[2]->FF_Assigner(Q2); FP = p_NP[3]->FF_Assigner(Q2);
  }
  // CNP: A' - B' p5slash
  double Ac(0.), Bc(0.), GDc(0.), F1c(0.), F2c(0.), FAc(0.), FPc(0.);
  if (p_CNP[0]) {
    Ac  = pc.Abs2()-sqr(m2);  Bc = m_mCNP+m2;
    GDc = m_cCNP*p_gpiNN->FF_Assigner(-mpi2)/m2N/(pc.Abs2()-sqr(m_mCNP));
    F1c = p_CNP[0]->FF_Assigner(Q2); F2c = p_CNP[1]->FF_Assigner(Q2);
    FAc = p_CNP[2]->FF_Assigner(Q2); FPc = p_CNP[3]->FF_Assigner(Q2);
  }
  // PF: everything but Y(g5) and (k1+k2)^mu
  const double F1V((m_onPF || m_cCT!=0.) ? p_F1V->FF_Assigner(Q2) : 0.);
  const double GPF(m_onPF ? m_cPF*p_gpiNN->FF_Assigner(-t)/m2N*(m4+m2)
                            /(t-m_mex2)*m_cX*F1V : 0.);
  // CT: Kroll-Ruderman, and for Z, W the axial + pion-pole terms (F_PiPole)
  const double KR(m_cCT*m_gA/m_fpi*F1V*(m_boson==kf_Z ? 1.-2.*m_s2 : 1.));
  const double FR((m_cCT!=0. && m_boson!=kf_photon) ?
                  m_cCT/m_fpi*p_FPiPole->FF_Assigner(m_pipolearg ? -t : Q2) : 0.);

  XYZFunc F(moms,m_flavs,anti,p_i);
  Vec4C amp[2][2];
  for (int h0=0;h0<2;h0++) {
    for (int h1=0;h1<2;h1++) {
      const Complex Y(F.Y(0,h0,1,h1,1.,1.)),       Y5(F.Y(0,h0,1,h1,1.,-1.));
      const Complex X(F.X(0,h0,p5,1,h1,1.,1.)),    X5(F.X(0,h0,p5,1,h1,1.,-1.));
      const Vec4C   L(F.L(0,h0,1,h1,1.,1.)),       L5(F.L(0,h0,1,h1,1.,-1.));
      const Vec4C   W(F.W(0,h0,p5,1,h1,1.,1.)),    W5(F.W(0,h0,p5,1,h1,1.,-1.));
      Vec4C J(Complex(0.,0.),Complex(0.,0.),Complex(0.,0.),Complex(0.,0.));
      if (p_NP[0]) {
        // P(chi) = ubar4 g5 (A + B p5slash) g^mu chi u2
        auto P = [&](const Complex & y,const Vec4C & l,const Vec4C & w)
          { return -A*l+(2.*B*y)*p5-B*w; };
        // sigma term: ubar4 g5 (A + B p5slash)(q^mu + m2 g^mu - g^mu p4slash
        //                                      - g^mu p5slash) g5-chirality u2
        const Vec4C sig = (A*Y5-B*X5)*q + m2*P(Y5,L5,W5)
          - ( A*((2.*Y5)*p4-m4*L5)
              + B*((-2.*X5)*p4+2.*p4p5*L5-m4*((2.*Y5)*p5-W5)) )
          - ( A*W5-(2.*B*X5)*p5+B*mpi2*L5 );
        J += Complex(0.,-GD)*( F1*P(Y5,L5,W5)-FA*P(Y,L,W)+(F2/m2N)*sig
                       -(FP/m2N*(A*Y-B*X))*q );
      }
      if (p_CNP[0]) {
        // sigma term: ubar4 (-q^mu + qslash g^mu)(A' - B' p5slash) g5 u2,
        // ubar4 qslash = ubar4 (m4 + p5slash - p2slash)
        const Vec4C sig = (-(Ac*Y5-Bc*X5))*q + m4*(Ac*L5-Bc*W5)
          + Ac*((2.*Y5)*p5-W5) - Bc*((2.*X5)*p5-mpi2*L5)
          - ( Ac*((2.*Y5)*p2+m2*L5)
              - Bc*((2.*X5)*p2-2.*p2p5*L5-m2*W5) );
        J += Complex(0.,-GDc)*( F1c*(Ac*L5-Bc*W5)-FAc*(Ac*L+Bc*W)+(F2c/m2N)*sig
                        -(FPc/m2N*(Ac*Y+Bc*X))*q );
      }
      if (m_onPF) J += (Complex(0.,-GPF)*Y5)*K;
      if (m_cCT!=0.) {
        J += Complex(0.,-KR)*L5;
        if (m_boson!=kf_photon)
          J += Complex(0.,FR)*( L+(((m4-m2)*Y+2.*X)/(2.*(Q2+m_mpole2)))*q );
      }
      amp[h0][h1] = m_norm*J;
      vector<pair<int,int> > spins;
      spins.push_back(make_pair(0,h0));
      spins.push_back(make_pair(1,h1));
      spins.push_back(make_pair(2,0));
      Insert(amp[h0][h1],spins);
    }
  }
  if (m_wardcheck && m_boson==kf_photon) CheckWard(q,amp);
}

void VA_N_NPi::CheckWard(const ATOOLS::Vec4D & q,const ATOOLS::Vec4C amp[2][2])
{
  // Gauge invariance, unit test 3 of implementation_1.tex: q.J = 0 for the
  // photon. Exact for m_p = m_n; otherwise O((m_n-m_p)/m_N) in the charged
  // channels. Done once.
  m_wardcheck = false;
  double qJ(0.), norm(0.);
  for (int h0=0;h0<2;h0++)
    for (int h1=0;h1<2;h1++) {
      qJ += std::abs(Vec4C(q)*amp[h0][h1]);
      for (size_t k(0);k<4;++k) norm += std::abs(q[k])*std::abs(amp[h0][h1][k]);
    }
  const double ratio(norm>0. ? qJ/norm : 0.);
  msg_Info()<<"### "<<m_name<<": Ward check |q.J|/(|q||J|) = "<<ratio<<" for "
            <<m_flavs[p_i[1]]<<" -> "<<m_flavs[p_i[0]]<<" "<<m_flavs[p_i[2]]<<".\n";
  if (ratio>1.e-2)
    THROW(fatal_error,"Current "+m_name+" fails the photon Ward identity, "
          "|q.J|/(|q||J|) = "+ToString(ratio)+".");
}

DECLARE_GETTER(METOOLS::VA_N_NPi,"VA_N_NPi",
               METOOLS::Current_Base,METOOLS::ME_Parameters);
METOOLS::Current_Base *ATOOLS::Getter<METOOLS::Current_Base,
                                      METOOLS::ME_Parameters,
                                      METOOLS::VA_N_NPi>::
operator()(const ME_Parameters &parameters) const
{
  // A nucleon beam with a PDF hands the hard process partons, never the
  // nucleon this current describes.
  for (size_t b(0);b<2;++b) {
    const Flavour bunch(rpa->gen.Bunch(b));
    if (rpa->gen.PDF(b)==nullptr) continue;
    for (size_t k(0);k<parameters.indices.size();++k)
      if (parameters.flavs[parameters.indices[k]]==bunch)
        THROW(fatal_error,"Current VA_N_NPi needs an intact "+ToString(bunch)
              +" beam, but beam "+ToString(b+1)+" has a PDF. Set PDF_SET: "
              "Intact for that beam.");
  }
  return new METOOLS::VA_N_NPi(parameters.flavs,parameters.indices,"VA_N_NPi");
}

void ATOOLS::Getter<METOOLS::Current_Base,
                    METOOLS::ME_Parameters,METOOLS::VA_N_NPi>::
PrintInfo(std::ostream &st,const size_t width) const {
  st<<"Example: the nucleon-pion side of $e^-p\\rightarrow e^-p\\pi^0$ \n\n"
    <<"Order: 0 = outgoing nucleon, 1 = incoming nucleon, 2 = pion \n\n"
    <<"Nucleon Born terms: nucleon pole, crossed nucleon pole, pion in "
    <<"flight, contact term. \n\n"
    <<"Parameters: {\\tt BOSON}, {\\tt VECTOR\\_FF}, {\\tt AXIAL\\_FF}, "
    <<"{\\tt PION\\_NN\\_FF}, {\\tt PIPOLE\\_FF}, {\\tt PIPOLE\\_FF\\_ARG}. \n"
    <<"The nucleon beam needs {\\tt PDF\\_SET: Intact}. \n"
    <<std::endl;
}
