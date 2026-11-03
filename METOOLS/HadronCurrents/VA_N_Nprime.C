#include "METOOLS/HadronCurrents/VA_N_Nprime.H"
#include "METOOLS/Main/XYZFuncs.H"
#include "MODEL/Main/Model_Base.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Phys/Flavour_Tags.H" 
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Message.H"

using namespace METOOLS;
using namespace ATOOLS;
using namespace std;

VA_N_Nprime::VA_N_Nprime(const ATOOLS::Flavour_Vector& flavs,
               const std::vector<int>& indices,
               const std::string& name) :
  Current_Base(flavs,indices,name),
  p_F1prime(NULL), p_F2prime(NULL),
  m_boson(kf_photon), m_N(0.), m_R(0.), m_sum(0.), m_diff(0.), m_norm(1.), m_vanishes(false)
{
  for (size_t k(0);k<p_i.size();++k) {
    const kf_code kf(m_flavs[p_i[k]].Kfcode());
    if (kf!=kf_p_plus && kf!=kf_n 
        && kf!=kf_N_1440 && kf!=kf_N_1440_plus
        && kf!=kf_N_1535 && kf!=kf_N_1535_plus)
      THROW(fatal_error,"Current "+m_name+" needs nucleon legs, got "
            +ToString(m_flavs[p_i[k]])+".");
  }
}


VA_N_Nprime::~VA_N_Nprime()
{
  if (p_F1prime) delete p_F1prime;
  if (p_F2prime) delete p_F2prime;
}


void VA_N_Nprime::CheckTopology() const
{
  if (p_i.size()!=2)
    THROW(fatal_error,"Current "+m_name+" needs exactly two legs, got "
          +ToString(p_i.size())+".");
  if (NIncomingLegs()!=1 || IsIncoming(p_i[0]))
    THROW(fatal_error,"Current "+m_name
          +": indices must be [outgoing resonance, incoming nucleon].");

  const kf_code kf_out(m_flavs[p_i[0]].Kfcode());
  const kf_code kf_in(m_flavs[p_i[1]].Kfcode());
  const bool outgoing_resonance =
    kf_out==kf_N_1440 || kf_out==kf_N_1440_plus ||
    kf_out==kf_N_1535 || kf_out==kf_N_1535_plus;
  const bool incoming_nucleon =
    kf_in==kf_p_plus || kf_in==kf_n;
  if (!outgoing_resonance || !incoming_nucleon)
    THROW(fatal_error,"Current "+m_name
          +": requires an outgoing N(1440)/N(1535) resonance "
          "and an incoming proton/neutron.");
  const bool out_charged =
    kf_out==kf_N_1440_plus || kf_out==kf_N_1535_plus;
  const bool in_charged = kf_in==kf_p_plus;
  if (out_charged != in_charged)
    THROW(fatal_error,"Current "+m_name
          +": photon transitions must conserve electric charge.");
}


FormFactor_Base * VA_N_Nprime::MakeFFprime(const std::string & name,const int ffmodel,
                                 GeneralModel & model)
{
  std::map<std::string,double> pmap;
  FF_Parameters params(ffmodel,m_flavs,p_i,pmap,name,&model);
  FormFactor_Base * ff = FF_Getter::GetObject("FF_N_Nprime",params);
  if (ff==NULL)
    THROW(fatal_error,"Could not build form factor '"+name+"' for current "
          +m_name+". Check that FF_N_Nprime.C is built and registered.");
  return ff;
}


void VA_N_Nprime::SetModelParameters(struct GeneralModel model)
{
  // Exchange_ME may call this a second time with default parameters.
  // Do not rebuild the transition form factors.
  if (p_F1prime) return;

  // This current describes electromagnetic transitions only.
  m_boson = kf_photon;

  // Average mass of the incoming nucleon and outgoing resonance.
  m_R   = p_masses[0];  // Outgoing resonance mass
  m_N   = p_masses[1];  // Incoming nucleon mass
  m_sum = m_N + m_R;    // Sum of masses
  m_diff = m_R - m_N;    // Difference of masses

  // No CKM factor for photon exchange.
  m_norm = 1.;

  // Select the transition vector form-factor model.
  const int vff(int(model("VECTOR_FF",10.) + 0.5));

  p_F1prime = MakeFFprime("F1prime",vff,model);
  p_F2prime = MakeFFprime("F2prime",vff,model);

  msg_Tracking()
    << "### " << m_name << ": "
    << m_flavs[p_i[1]] << " -> " << m_flavs[p_i[0]]
    << " via " << Flavour(m_boson)
    << ", VECTOR_FF = " << vff
    << ", norm = " << m_norm << "\n";
}


void VA_N_Nprime::Calc(const ATOOLS::Vec4D_Vector& moms, bool anti)
{
  const Vec4D q(moms[p_i[0]]-moms[p_i[1]]);
  const Vec4D P(moms[p_i[0]]+moms[p_i[1]]);
  const double Q2(-q.Abs2());
  if (Q2 == 0.0)
    THROW(fatal_error,
          "Current "+m_name+": Q2 = 0 encountered in VA_N_Nprime::Calc().");

  const double F1prime((*p_F1prime)(moms).real());
  const double F2prime((*p_F2prime)(moms).real());
  const kf_code kf_out(m_flavs[p_i[0]].Kfcode());
  const bool positive_parity(
    kf_out==kf_N_1440 || kf_out==kf_N_1440_plus
  );

  // Parity-dependent coefficient from the transition-current convention.
  const double cP(positive_parity ? 1. : (m_diff)/m_sum);
  // Gordon-reduced coefficient of gamma^mu.
  const double V(F1prime+cP*F2prime);
  const Complex cR(0.,V), cL(0.,positive_parity ? V : -V);
  const Complex sR(0.,1.), sL(0.,positive_parity ? 1. : -1.);
  const Complex pR(0.,1.), pL(0.,positive_parity ? 1. : -1.);

  XYZFunc F(moms,m_flavs,anti,p_i);

  for (int h0=0;h0<2;h0++) {
    for (int h1=0;h1<2;h1++) {
      const Vec4C amp = Complex(m_norm,0.)*
        ( F.L(0,h0,1,h1,cR,cL)
        + (F.X(0,h0,q,1,h1,sR,sL)*Complex(F1prime/(-q.Abs2()),0.))*Vec4C(q)
        - (F.Y(0,h0,1,h1,pR,pL)*Complex(F2prime/m_sum,0.))*Vec4C(P) );

      vector<pair<int,int> > spins;
      spins.push_back(make_pair(0,h0));
      spins.push_back(make_pair(1,h1));
      Insert(amp,spins);
    }
  }
}


DECLARE_GETTER(METOOLS::VA_N_Nprime,"VA_N_Nprime",
               METOOLS::Current_Base,METOOLS::ME_Parameters);
METOOLS::Current_Base *ATOOLS::Getter<METOOLS::Current_Base,
                                      METOOLS::ME_Parameters,
                                      METOOLS::VA_N_Nprime>::
operator()(const ME_Parameters &parameters) const
{
  // A nucleon beam with a PDF hands the hard process partons, never the
  // nucleon this current describes.
  for (size_t b(0);b<2;++b) {
    const Flavour bunch(rpa->gen.Bunch(b));
    if (rpa->gen.PDF(b)==nullptr) continue;
    for (size_t k(0);k<parameters.indices.size();++k)
      if (parameters.flavs[parameters.indices[k]]==bunch)
        THROW(fatal_error,"Current VA_N_Nprime needs an intact "+ToString(bunch)
              +" beam, but beam "+ToString(b+1)+" has a PDF. Set PDF_SET: "
              "Intact for that beam.");
  }
  return new METOOLS::VA_N_Nprime(parameters.flavs,parameters.indices,"VA_N_Nprime");
}

void ATOOLS::Getter<METOOLS::Current_Base,
                    METOOLS::ME_Parameters,METOOLS::VA_N_Nprime>::
PrintInfo(std::ostream &st,const size_t width) const {
  st << "Example: the resonance transition current "
   << "of a nucleon to N(1440) or N(1535).\n\n"
   << "Order: 0 = outgoing resonance, 1 = incoming nucleon.\n\n"
   << "Parameters: {\\tt VECTOR\\_FF}.\n\n"
   << "The nucleon beam needs {\\tt PDF\\_SET: Intact}.\n"
   << std::endl;
}
