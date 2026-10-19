#include "METOOLS/HadronCurrents/VA_N_N.H"
#include "METOOLS/Main/XYZFuncs.H"
#include "MODEL/Main/Model_Base.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Message.H"

using namespace METOOLS;
using namespace ATOOLS;
using namespace std;

VA_N_N::VA_N_N(const ATOOLS::Flavour_Vector& flavs,
               const std::vector<int>& indices,
               const std::string& name) :
  Current_Base(flavs,indices,name),
  p_F1p(NULL), p_F2p(NULL), p_F1n(NULL), p_F2n(NULL), p_FA(NULL),
  m_boson(kf_photon), m_mN(0.), m_mpi(0.),
  m_s2(0.), m_T3(0.), m_norm(1.), m_isp(true)
{
  for (size_t k(0);k<p_i.size();++k) {
    const kf_code kf(m_flavs[p_i[k]].Kfcode());
    if (kf!=kf_p_plus && kf!=kf_n)
      THROW(fatal_error,"Current "+m_name+" needs nucleon legs, got "
            +ToString(m_flavs[p_i[k]])+".");
  }
}

VA_N_N::~VA_N_N()
{
  if (p_F1p) delete p_F1p;
  if (p_F2p) delete p_F2p;
  if (p_F1n) delete p_F1n;
  if (p_F2n) delete p_F2n;
  if (p_FA)  delete p_FA;
}

void VA_N_N::CheckTopology() const
{
  if (p_i.size()!=2)
    THROW(fatal_error,"Current "+m_name+" needs exactly two legs, got "
          +ToString(p_i.size())+".");
  if (NIncomingLegs()!=1 || IsIncoming(p_i[0]))
    THROW(fatal_error,"Current "+m_name+" is the quasi-elastic N -> N' "
          "current: Indices must be [outgoing, incoming] nucleon, got ["
          +ToString(p_i[0])+", "+ToString(p_i[1])+"].");
}

FormFactor_Base * VA_N_N::MakeFF(const std::string & name,const int ffmodel,
                                 GeneralModel & model)
{
  std::map<std::string,double> pmap;
  FF_Parameters params(ffmodel,m_flavs,p_i,pmap,name,&model);
  FormFactor_Base * ff = FF_Getter::GetObject("FF_N_N",params);
  if (ff==NULL)
    THROW(fatal_error,"Could not build form factor '"+name+"' for current "
          +m_name+". Check that FF_N_N.C is built and registered.");
  return ff;
}

void VA_N_N::SetModelParameters(struct GeneralModel model)
{
  // Exchange_ME calls this a second time with only the default v and a.
  // Those are not used here, and must not rebuild the form factors with
  // default settings.
  if (p_FA) return;

  const kf_code kf0(m_flavs[p_i[0]].Kfcode()), kf1(m_flavs[p_i[1]].Kfcode());
  const bool cc(kf0!=kf1);
  m_boson = static_cast<kf_code>(std::abs(int(model("BOSON", cc ? 24. : 22.))));
  if (m_boson != kf_photon && m_boson != kf_Z && m_boson != kf_Wplus)
    THROW(fatal_error,"Current "+m_name+": BOSON must be 22, 23 or 24.");
  if (cc!=(m_boson==kf_Wplus))
    THROW(fatal_error,"Current "+m_name+": BOSON "+ToString(m_boson)+
          " does not match the nucleon legs "+ToString(m_flavs[p_i[0]])+
          ", "+ToString(m_flavs[p_i[1]])+".");

  m_isp  = (kf0==kf_p_plus);
  m_T3   = m_isp ? 0.5 : -0.5;
  m_mN   = 0.5*(p_masses[0]+p_masses[1]);
  m_mpi  = Flavour(m_boson==kf_Z ? kf_pi : kf_pi_plus).HadMass();
  m_s2   = std::abs(MODEL::s_model->ComplexConstant(string("csin2_thetaW")));
  m_norm = (m_boson==kf_Wplus) ?
    std::abs(MODEL::s_model->ComplexConstant(string("CKM_0_0"))) : 1.;

  const int vff(int(model("VECTOR_FF",10.)+0.5));
  const int aff(m_boson==kf_photon ? ff_zero : int(model("AXIAL_FF",20.)+0.5));
  p_F1p = MakeFF("F1p",vff,model);
  p_F2p = MakeFF("F2p",vff,model);
  p_F1n = MakeFF("F1n",vff,model);
  p_F2n = MakeFF("F2n",vff,model);
  p_FA  = MakeFF("FA", aff,model);

  msg_Tracking()<<"### "<<m_name<<": "<<m_flavs[p_i[1]]<<" -> "
                <<m_flavs[p_i[0]]<<" via "<<Flavour(m_boson)
                <<", VECTOR_FF = "<<vff<<", AXIAL_FF = "<<aff
                <<", norm = "<<m_norm<<"\n";
}

void VA_N_N::Calc(const ATOOLS::Vec4D_Vector& moms, bool anti)
{
  const Vec4D  q(moms[p_i[0]]-moms[p_i[1]]), P(moms[p_i[0]]+moms[p_i[1]]);
  const double Q2(-q.Abs2()), m2N(2.*m_mN);

  const double F1p((*p_F1p)(moms).real()), F2p((*p_F2p)(moms).real());
  const double F1n((*p_F1n)(moms).real()), F2n((*p_F2n)(moms).real());
  const double F1N(m_isp ? F1p : F1n), F2N(m_isp ? F2p : F2n);
  double F1(F1p-F1n), F2(F2p-F2n), FA((*p_FA)(moms).real());
  if (m_boson==kf_photon) { F1 = F1N; F2 = F2N; }
  else if (m_boson==kf_Z) {
    F1 = m_T3*F1-2.*m_s2*F1N;
    F2 = m_T3*F2-2.*m_s2*F2N;
    FA = m_T3*FA;
  }
  const double FP(4.*sqr(m_mN)/(Q2+sqr(m_mpi))*FA);

  // Gordon: gamma^mu coefficient F1 + (m'+m)/(2m_N) F2, plus -P^mu F2/(2m_N).
  const double V(F1+(p_masses[0]+p_masses[1])/m2N*F2);
  const Complex cR(0.,V+FA), cL(0.,V-FA);      // i gamma^mu (V - F_A g5)
  const Complex sR(0.,1.),   sL(0.,1.);        // i
  const Complex pR(0.,FP),   pL(0.,-FP);       // i (-F_P g5)

  XYZFunc F(moms,m_flavs,anti,p_i);
  for (int h0=0;h0<2;h0++) {
    for (int h1=0;h1<2;h1++) {
      const Vec4C amp = Complex(m_norm,0.)*
        ( F.L(0,h0,1,h1,cR,cL)
        - (F.Y(0,h0,1,h1,sR,sL)*Complex(F2/m2N,0.))*Vec4C(P)
        + (F.Y(0,h0,1,h1,pR,pL)*Complex(1./m2N,0.))*Vec4C(q) );
      vector<pair<int,int> > spins;
      spins.push_back(make_pair(0,h0));
      spins.push_back(make_pair(1,h1));
      Insert(amp,spins);
    }
  }
}

DECLARE_GETTER(METOOLS::VA_N_N,"VA_N_N",
               METOOLS::Current_Base,METOOLS::ME_Parameters);
METOOLS::Current_Base *ATOOLS::Getter<METOOLS::Current_Base,
                                      METOOLS::ME_Parameters,
                                      METOOLS::VA_N_N>::
operator()(const ME_Parameters &parameters) const
{
  // A nucleon beam with a PDF hands the hard process partons, never the
  // nucleon this current describes.
  for (size_t b(0);b<2;++b) {
    const Flavour bunch(rpa->gen.Bunch(b));
    if (rpa->gen.PDF(b)==nullptr) continue;
    for (size_t k(0);k<parameters.indices.size();++k)
      if (parameters.flavs[parameters.indices[k]]==bunch)
        THROW(fatal_error,"Current VA_N_N needs an intact "+ToString(bunch)
              +" beam, but beam "+ToString(b+1)+" has a PDF. Set PDF_SET: "
              "Intact for that beam.");
  }
  return new METOOLS::VA_N_N(parameters.flavs,parameters.indices,"VA_N_N");
}

void ATOOLS::Getter<METOOLS::Current_Base,
                    METOOLS::ME_Parameters,METOOLS::VA_N_N>::
PrintInfo(std::ostream &st,const size_t width) const {
  st<<"Example: the proton line of $e^-p\\rightarrow e^-p$ \n\n"
    <<"Order: 0 = outgoing nucleon, 1 = incoming nucleon \n\n"
    <<"\\[\\bar{u}(p') [\\gamma^\\mu F_1 + i\\sigma^{\\mu\\nu}q_\\nu F_2/2m_N"
    <<" - (\\gamma^\\mu F_A + q^\\mu F_P/2m_N)\\gamma_5 ] u(p) \\] \n\n"
    <<"Parameters: {\\tt BOSON}, {\\tt VECTOR\\_FF}, {\\tt AXIAL\\_FF}. \n"
    <<"The nucleon beam needs {\\tt PDF\\_SET: Intact}. \n"
    <<std::endl;
}
