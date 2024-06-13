#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "PHASIC++/Main/Phase_Space_Handler.H"
#include "PHASIC++/Main/Phase_Space_Point.H"
#include "PHASIC++/Main/Process_Integrator.H"

using namespace PHASIC;
using namespace ATOOLS;

Phase_Space_Point::Phase_Space_Point(Phase_Space_Handler * psh)
    : p_pshandler(psh), p_beamhandler(NULL), p_isrhandler(NULL),
      p_moms(p_pshandler->Momenta()), p_cuts(NULL), p_beamchannels(NULL), p_isrchannels(NULL),
      p_fsrchannels(NULL), m_Ecms(ATOOLS::rpa->gen.Ecms()), m_smin(0.),
      m_weight(0.), m_ISsymmetryfactor(1.) {}

Phase_Space_Point::~Phase_Space_Point() {
  if (p_fsrchannels)
    delete p_fsrchannels;
  if (p_isrchannels)
    delete p_isrchannels;
  if (p_beamchannels)
    delete p_beamchannels;
  if (p_cuts)
    delete p_cuts;
}

///////////////////////////////////////////////////////////////////////////////
//
// Initialization
//
//////////////////////////////////////////////////////////////////////////////

void Phase_Space_Point::Init() {
  p_beamhandler = p_pshandler->GetBeamSpectra();
  p_isrhandler = p_pshandler->GetISRHandler();
  p_yfshandler = p_pshandler->GetYFSHandler();
  m_nin = p_pshandler->Process()->Process()->NIn();
  m_nout = p_pshandler->Process()->Process()->NOut();
  m_nvec = m_nin + m_nout;
  Flavour_Vector flavours = p_pshandler->Flavs();
  m_masses.resize(flavours.size());
  for (size_t i = 0; i < flavours.size(); i++) {
    m_masses[i] = flavours[i].Mass();
    if (i < m_nin)
      m_masses2[i] = sqr(m_masses[i]);
  }
  m_osmass =
      (m_nout == 1 && flavours[2].Kfcode() != kf_instanton ? m_masses[m_nin] : 0.0);
  m_beamspkey.Assign(std::string("BEAM::s'"), 5, 0, p_pshandler->GetInfo());
  m_beamykey.Assign(std::string("BEAM::y"), 3, 0, p_pshandler->GetInfo());
  p_beamhandler->AssignKeys(p_pshandler->GetInfo());
  m_isrspkey.Assign(std::string("ISR::s'"), 5, 0, p_pshandler->GetInfo());
  m_isrykey.Assign(std::string("ISR::y"), 3, 0, p_pshandler->GetInfo());
  p_isrhandler->AssignKeys(p_pshandler->GetInfo());
  InitFixedIncomings();
  msg_Tracking() << "================================================\n"
                 << METHOD << " for "
                 << (p_beamchannels ? p_beamchannels->NChannels() : 0)
                 << " beam channels, "
                 << (p_isrchannels ? p_isrchannels->NChannels() : 0)
                 << " ISR channels, "
                 << (p_fsrchannels ? p_fsrchannels->NChannels() : 0)
                 << " FSR channels:\n";
  if (p_beamchannels) {
    p_beamchannels->Reset();
    p_beamchannels->Print();
  }
  if (p_isrchannels) {
    p_isrchannels->Reset();
    p_isrchannels->Print();
  }
  if (p_fsrchannels) {
    p_fsrchannels->Reset();
    p_fsrchannels->Print();
  }
}

void Phase_Space_Point::InitFixedIncomings() {
  if (m_nin == 1) {
    m_Eprime = m_masses[0];
    m_sprime = m_fixedsprime = sqr(m_Eprime);
    m_y = m_fixedy = 0.;
    m_ISmoms[0] = Vec4D(m_Eprime, 0., 0., 0.);
  } else if (m_nin == 2) {
    if (p_beamhandler->On() || p_isrhandler->On() != 0)
      return;
    for (int i = 0; i < 2; ++i) {
      m_ISmoms[i] = p_beamhandler->GetBeam(i)->OutMomentum();
    }
    m_sprime = m_fixedsprime = (m_ISmoms[0] + m_ISmoms[1]).Abs2();
    m_Eprime = sqrt(m_sprime);
    m_y = m_fixedy = (m_ISmoms[0] + m_ISmoms[1]).Y();
  }
}

void Phase_Space_Point::InitCuts(Process_Integrator *const process) {
  if (p_cuts != NULL)
    delete p_cuts;
  p_cuts = new Cut_Data();
  process->Process()->InitCuts(p_cuts);
  process->Process()->FillOnshellConditions();
  process->Process()->BuildCuts(p_cuts);
  if (process->NIn() > 1) {
    m_smin = ATOOLS::Max(sqr(process->ISRThreshold()), p_cuts->Smin());
    process->ISR()->SetFixedSprimeMin(m_smin);
    process->Beam()->SetSprimeMin(m_smin);
  }
}

///////////////////////////////////////////////////////////////////////////////
//
// Generation of single phase space point
//
//////////////////////////////////////////////////////////////////////////////

bool Phase_Space_Point::operator()(Process_Integrator *const process,
                                   const psmode::code &mode) {
  p_pshandler->GetInfo()->ResetAll();
  Reset(mode);
  if (m_nin == 2) {
    if (!DefineBeamKinematics())
      return false;
    if (!DefineISRKinematics(process)) {
      if (p_beamchannels)
        p_beamchannels->NoGenerate();
      if (p_isrchannels)
        p_isrchannels->NoGenerate();
      p_fsrchannels->NoGenerate();
      return false;
    }
  }
  if(p_yfshandler->Mode()==YFS::yfsmode::off) DefineFSRKinematics();
  CorrectMomenta();
  return true;
}

bool Phase_Space_Point::DefineBeamKinematics() {
  // active beam handling necessary - generate a point:
  // -- s' and rapidity y for collider mode,
  // -- s' only for relic density mode (integration over frame implicit)
  // -- s', y = E1/(E1+E2), and cos(theta) for DM annihilation mode
  if (p_beamhandler->On() && p_beamchannels != NULL) {
    p_beamhandler->SetLimits();
    p_beamchannels->GeneratePoint();
    if (!p_beamhandler->MakeBeams(p_moms))
      return false;
  }
  m_sprime = p_beamhandler->Sprime();
  m_y = (p_beamhandler->GetBeam(0)->InMomentum() +
         p_beamhandler->GetBeam(1)->InMomentum())
            .Y();
  m_y += p_beamhandler->Y();
  return true;
}

bool Phase_Space_Point::DefineISRKinematics(Process_Integrator *const process) {
  // active ISR handling necessary - generate a point:
  // -- s' and rapidity y for collider mode,
  if (p_isrhandler->On() && p_isrchannels != NULL &&
      m_mode != psmode::no_gen_isr) {
    if (m_smin > m_sprime * p_isrhandler->Upper1() * p_isrhandler->Upper2()) {
      msg_Tracking() << METHOD << ": Event rejected due to insufficient s' \n"
                     << "    Upper limits for x are " << p_isrhandler->Upper1()
                     << " and " << p_isrhandler->Upper2() << "\n";
      return false;
    }
    p_isrhandler->Reset();
    if (!(m_mode & psmode::no_lim_isr)) {
      p_isrhandler->SetSprimeMax(m_sprime * p_isrhandler->Upper1() *
                                 p_isrhandler->Upper2());
      p_isrhandler->SetSprimeMin(m_smin);
    }
    p_isrhandler->SetPole(m_sprime);
    if (!(m_mode & psmode::no_gen_isr)) {
      p_isrhandler->SetMasses(process->Process()->Selected()->Flavours());
      p_isrhandler->SetLimits(m_beamykey[2]);
      if (!p_isrhandler->CheckMasses())
        return false;
      if (m_nin == 2 && m_nout == 1 &&
          process->Process()->Selected()->Flavours()[2].Kfcode() == kf_instanton) {
        if (p_pshandler->Active()->Process()->SPrimeMin() > 0.)
          m_isrspkey[0] = p_pshandler->Active()->Process()->SPrimeMin();
        if (p_pshandler->Active()->Process()->SPrimeMax() > 0.)
          m_isrspkey[1] = p_pshandler->Active()->Process()->SPrimeMax();
      }
      m_isrspkey[4] = sqr(m_osmass);
      p_isrchannels->GeneratePoint();
    }
    m_sprime = m_osmass ? m_isrspkey[4] : m_isrspkey[3];
    m_y += m_isrykey[2];
    m_ISsymmetryfactor = p_isrhandler->GenerateSwap(
                             p_pshandler->Active()->Process()->Flavours()[0],
                             p_pshandler->Active()->Process()->Flavours()[1])
                             ? 2.0
                             : 1.0;
  }
  if(p_yfshandler->HasISR()){
    p_isrchannels->GeneratePoint();
    p_yfshandler->SetLimits(m_smin);
    DefineFSRKinematics();
    p_yfshandler->SetBornMomenta(p_moms);
    // p_yfshandler->SetFlavours(p_pshandler->Active()->Process()->Flavours());
    m_sprime = m_osmass ? m_isrspkey[4] : m_isrspkey[3];
    p_yfshandler->SetMomenta(p_moms);
    p_yfshandler->SetSprime(m_sprime);
    if(!p_yfshandler->MakeYFS(p_moms)) return 0;
    DefineFSRKinematics();
    p_yfshandler->SetMomenta(p_moms);
    return(p_yfshandler->CalculateFSR());
  }
  else if(p_yfshandler->Mode()==YFS::yfsmode::fsr){
    DefineFSRKinematics();
    p_yfshandler->SetBornMomenta(p_moms);
    return(p_yfshandler->CalculateFSR(p_moms));

  }
  return p_isrhandler->MakeISR(m_sprime, m_isrykey[2], p_moms,
                               process->Process()->Selected()->Flavours());
}

bool Phase_Space_Point::DefineFSRKinematics() {
  p_fsrchannels->GeneratePoint(p_moms.data(), p_pshandler->Cuts());
  return true;
}


void Phase_Space_Point::CorrectMomenta() {
  if (m_nin!=2 || (m_nin==2 && m_nout==1 &&
       p_pshandler->Active()->Process()->Flavours()[2].Kfcode()==kf_instanton))
    return;
  if(p_yfshandler->Mode()==YFS::yfsmode::isr) return;
  Vec4D  momsum(0.,0.,0.,0.);
  size_t imax(0);
  double Emax(0.0);
  for (size_t i(0); i < m_nin; ++i)
    momsum += -p_moms[i];
  for (size_t i(m_nin); i < m_nvec; ++i) {
    momsum += p_moms[i];
    if (p_moms[i][0] > Emax) {
      Emax = p_moms[i][0];
      imax = i;
    }
    p_moms[i][0] = sqrt(p_moms[i].PSpat2() + sqr(m_masses[i]));
  }
  p_moms[imax] -= momsum;
  p_moms[imax][0] = sqrt(p_moms[imax].PSpat2() + sqr(m_masses[imax]));
  double E0tot(0);
  for (size_t i(0); i < m_nvec; ++i)
    E0tot += (i < m_nin ? -1. : 1.) * p_moms[i][0];
  double p2[2] = {p_moms[0].PSpat2(), p_moms[1].PSpat2()};
  double E0[2] = {-p_moms[0][0], -p_moms[1][0]};
  double E1[2] = {p2[0] / E0[0],
                  -(Vec3D(p_moms[0]) * Vec3D(p_moms[1])) / E0[1]};
  double E1tot = E1[0] + E1[1];
  double E2[2] = {p2[0] * m_masses2[0] / (2 * pow(E0[0], 3)),
                  (p2[0] - sqr(E1[1])) / (2 * E0[1])};
  double E2tot = E2[0] + E2[1];
  double eps1 = -E0tot / E1tot;
  double eps = eps1 * (1 - eps1 * E2tot / E1tot);
  p_moms[1] = -p_moms[1] + p_moms[0] * eps;
  p_moms[0] = -p_moms[0] - p_moms[0] * eps;
  // Indices in m_masses need to be shifted if swap happened in ISR_Handler
  int swap = p_isrhandler->Swap();
  for (int i(0); i < 2; ++i) {
    if (m_masses[i^swap] == 0.0 && p_moms[i][1] == 0.0 && p_moms[i][2] == 0.0)
      p_moms[i][0] = -std::abs(p_moms[i][3]);
    else
      p_moms[i][0] = E0[i] + E1[i] * eps + E2[i] * sqr(eps);
  }
  for (size_t i(0); i < m_nin; ++i)
    p_moms[i] = -p_moms[i];
}

double Phase_Space_Point::CalculateWeight() {
  m_weight = 0.;
  if (Check4Momentum()) {
    m_weight = 1.0;
    if (p_isrchannels && !(m_mode & psmode::no_gen_isr)) {
      p_isrchannels->GenerateWeight();
      m_weight *= p_isrchannels->Weight();
    }
    if (p_beamchannels) {
      p_beamchannels->GenerateWeight();
      m_weight *= p_beamchannels->Weight();
    }
    p_fsrchannels->GenerateWeight(p_moms.data(), p_cuts);
    m_weight *= p_fsrchannels->Weight();
  }
  return m_weight;
}

bool Phase_Space_Point::Check4Momentum() {
  Vec4D pin, pout;
  for (int i = 0; i < m_nin; i++)
    pin += p_moms[i];
  for (int i = m_nin; i < m_nin + m_nout; i++)
    pout += p_moms[i];
  double sin = pin.Abs2(), sout = pout.Abs2();
  static double accu(sqrt(Accu()));
  if (!IsEqual(pin, pout, accu) || !IsEqual(sin, sout, accu)) {
    int prec(ATOOLS::msg->Error().precision());
    msg_Error().precision(12);
    msg_Error() << "ERROR in " << METHOD << ": [accu = " << accu << "] {\n";
    for (int i = 0; i < m_nin + m_nout; ++i)
      msg_Error() << "   p_" << i << " = " << p_moms[i] << " ("
                  << p_moms[i].Abs2() << ")\n";
    msg_Error() << "   p_in  = " << pin << " (" << sin << ")\n"
                << "   p_out = " << pout << " (" << sout << ")\n"
                << "   diff  = " << pout - pin << " (" << sout - sin << ")\n}\n"
                << "   Will return 0 as phase space weight.\n";
    msg_Error().precision(prec);
    return false;
  }
  return true;
}

void Phase_Space_Point::Print(std::ostream &str) {
  str << METHOD << " is generating phase space points with:\n";
  if (m_nin > 1) {
    if (p_beamchannels)
      str << "  Beam   : " << p_beamchannels->Name() << " (" << p_beamchannels
          << ") "
          << "  (" << p_beamchannels->Number() << "," << p_beamchannels->N()
          << ")\n";
    if (p_isrchannels)
      str << "  ISR    : " << p_isrchannels->Name() << " (" << p_isrchannels
          << ") "
          << "  (" << p_isrchannels->Number() << "," << p_isrchannels->N()
          << ")\n";
  }
  str << "  FSR    : " << p_fsrchannels->Name() << " (" << p_fsrchannels << ") "
      << "  (" << p_fsrchannels->Number() << "," << p_fsrchannels->N() << ")\n";
  str << "Printing all channels:\n";
  if (p_beamchannels)
    p_beamchannels->Print();
  if (p_isrchannels)
    p_isrchannels->Print();
  p_fsrchannels->Print();
}
