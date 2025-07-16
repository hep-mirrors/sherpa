#include "AddOns/EWSud/Clustered_Calculator.H"

#include "ATOOLS/Org/Run_Parameter.H"
#include "PHASIC++/Process/ME_Generator_Base.H"
#include "PHASIC++/Process/ME_Generators.H"
#include "SHERPA/PerturbativePhysics/Matrix_Element_Handler.H"
#include "SHERPA/SoftPhysics/Resonance_Finder.H"

using namespace ATOOLS;
using namespace PHASIC;
using namespace SHERPA;
using namespace EWSud;

Clustered_Calculator::Clustered_Calculator(Process_Base* _proc)
  : proc{_proc}
{
  DEBUG_FUNC(proc->Name());

  Scoped_Settings ewsudsettings{
    Settings::GetMainSettings()["EWSUD"] };
  m_resdist =
    ewsudsettings["CLUSTERING_THRESHOLD"].SetDefault(10.0).Get<double>();
  m_disabled =
    ewsudsettings["CLUSTERING_DISABLED"].SetDefault(false).Get<bool>();
  if(Settings::GetMainSettings()["EWSUDAKOV_CLUSTERING_DISABLED"].IsSetExplicitly()){
    THROW(fatal_error, "Avoid Using old syntax, prefer the new EWSUD: CLUSTERING_DISABLED");
  }
  auto ampl = Amplitudes::CreateAmplitude(proc);
  const Flavour_Vector& flavs = ampl->Flavs();

  // Add calculator for the unclustered base process and capture its Comix
  // interface, which we will use to build clustered processes
  auto base_calculator =
    std::unique_ptr<Calculator>(new Calculator{proc});
  calculators.emplace(
      std::make_pair(flavs, std::move(base_calculator)));
  p_comixinterface = &calculators.begin()->second->GetComixInterface();

  // Add calculators for clustered processes
  if (!m_disabled) {
    AddCalculators(ampl, 0);
  }

  msg_Debugging() << "Added " << calculators.size() << " calculators\n";
}

void Clustered_Calculator::AddCalculators(const Cluster_Amplitude_UP& ampl, size_t clusterings)
{
  DEBUG_FUNC(*ampl);
  ClusterLeg_Vector& legs = ampl->Legs();
  size_t n_legs = legs.size();
  for (size_t i {ampl->NIn()}; i < n_legs; ++i) {
    if (ampl->Flav(i).IsLepton()) {
      for (size_t j {i + 1}; j < n_legs; ++j) {
        auto clustered_kfc = kf_none;
        if (ampl->Flav(i) == ampl->Flav(j).Bar()) {
          clustered_kfc = kf_Z;
        } else if (ampl->Flav(i) == ampl->Flav(j).IsoWeakPartner()) {
          clustered_kfc = kf_Wplus;
        }
        if (clustered_kfc != kf_none) {
          Flavour flav {clustered_kfc};
          if (clustered_kfc == kf_Wplus
              && (ampl->Flav(i).Charge() + ampl->Flav(j).Charge() < 0)) {
            flav = flav.Bar();
          }
          Cluster_Amplitude_UP new_ampl = CopyClusterAmpl(ampl);
          new_ampl->CombineLegs(new_ampl->Leg(i), new_ampl->Leg(j), flav);
          AddCalculators(new_ampl, clusterings + 1);
          break;
        }
      }
    }
  }
  // Add current calculator (except the base one, which has already been added
  // to calculaturs in the ctor).
  if (clusterings > 0)
    AddCalculator(ampl, clusterings);
}

Clustered_Calculator::~Clustered_Calculator()
{
}

void Clustered_Calculator::AddCalculator(const Cluster_Amplitude_UP& ampl, size_t clusterings)
{
  const Flavour_Vector& flavs = ampl->Flavs();
  auto it = calculators.find(flavs);
  if (it != calculators.end())
    return;

  // Create clustered process using the base calculator's COMIX interface
  PHASIC::Process_Base* clustered_proc = p_comixinterface->GetProcess(*ampl);
  if (!clustered_proc) {
    Process_Info pi = p_comixinterface->CreateProcessInfo(ampl.get());
    pi.m_mincpl[1] -= clusterings;
    pi.m_maxcpl[1] -= clusterings;
    clustered_proc = p_comixinterface->InitializeProcess(pi);
  }

  // add calculator
  auto calculator = std::unique_ptr<Calculator>(
      new Calculator{clustered_proc});
  calculators.emplace(
      std::make_pair(flavs, std::move(calculator)));
}

double Clustered_Calculator::CalcIClustered(
    const std::map<double, std::vector<long int>> restab,
    const Vec4D_Vector &mom, const Flavour_Vector &flavs)
{
  double ClusteredIOperator{0.0};
  const auto& EWConsts = calculators[flavs]->GetEWGroupConstants();
  for(const auto& clij: restab){
    const auto i = clij.second[0];
    const auto j = clij.second[1];
    const double Qi  {flavs[i].Charge()};
    const double Qj  {flavs[j].Charge()};
    const double sij {(mom[i] + mom[j]).Abs2()};
    const double mi2 {sqr(flavs[i].Mass())};
    const double mj2 {sqr(flavs[j].Mass())};
    ClusteredIOperator += 2. * sqr(log(EWConsts.m_mw2 / sij));
    if(mi2 != 0.0)
      ClusteredIOperator -= sqr(log(mi2/EWConsts.m_mw2));
    if(mj2 != 0.0)
      ClusteredIOperator -= sqr(log(mj2/EWConsts.m_mw2));
    ClusteredIOperator *= Qi * Qj;
  }
  return EWConsts.delta_prefactor * ClusteredIOperator;
}

EWSudakov_Log_Corrections_Map
Clustered_Calculator::CorrectionsMap(Vec4D_Vector mom)
{
  DEBUG_FUNC(mom);

  if (m_disabled) {
    return calculators.begin()->second->CorrectionsMap(mom);
  }

  Flavour_Vector flavs = proc->Flavours();
  double ClusteredIOperator{0.0};

  if (msg_LevelIsDebugging()) {
    msg_Out() << "Will use input process for EWSud: ";
    for (const auto& flav : flavs) {
      msg_Out() << flav.ShellName() << " ";
    }
    msg_Out() << '\n';
  }

  size_t nflavs {flavs.size()};
  std::map<double, std::vector<long int>> restab;
  for (long int i {0}; i < nflavs; ++i) {
    if (flavs[i].IsLepton()) {
      for (long int j {i + 1}; j < nflavs; ++j) {
        long int kf {kf_none};
        if (flavs[j] == flavs[i].Bar()) {
          kf = kf_Z;
        } else if (flavs[j] == flavs[i].IsoWeakPartner()) {
          kf = kf_Wplus;
        }
        if (kf != kf_none) {
          Flavour resonance{kf};
          double invariant_mass {(mom[i]+mom[j]).Mass()};
          double mdist{
            std::abs(invariant_mass - resonance.Mass()) / resonance.Width()};
          msg_Debugging() << "found resonance candidate " << i << ", " << j <<
            " -> " << kf << " (" << mdist << ") ";
          if (mdist < m_resdist) {
            msg_Debugging()<<"-> accept\n";
            long int ida[3]={i,j,kf};
            restab[mdist]=std::vector<long int>(ida,ida+3);
          } else {
            msg_Debugging()<<"-> reject\n";
          }
        }
      }
    }
  }

  // replace resonances starting with the least off-shell one
  std::vector<double> clusterings;
  std::unordered_set<size_t> clustered_indizes;
  for (const auto& mdist_ida_pair : restab) {
    if (clustered_indizes.find(mdist_ida_pair.second[0]) == clustered_indizes.end()
        && clustered_indizes.find(mdist_ida_pair.second[1]) == clustered_indizes.end()) {
      clusterings.push_back(mdist_ida_pair.first);
      clustered_indizes.insert(mdist_ida_pair.second[0]);
      clustered_indizes.insert(mdist_ida_pair.second[1]);
    }
  }
  std::set<size_t> removelist;
  ClusteredIOperator  = CalcIClustered(restab,mom,flavs);
  for (const auto& mdist : clusterings) {
    if (restab[mdist][2] == kf_Wplus) {
      flavs[restab[mdist][0]] =
        Flavour{kf_Wplus, (flavs[restab[mdist][1]].Charge() + flavs[restab[mdist][2]].Charge() < 0)};
    } else {
      flavs[restab[mdist][0]] =
        Flavour{kf_Z};
    }
    mom[restab[mdist][0]] += mom[restab[mdist][1]];
    removelist.insert(restab[mdist][1]);
  }
  for (auto i = removelist.rbegin(); i != removelist.rend(); ++ i) {
    flavs.erase(flavs.begin() + *i);
    mom.erase(mom.begin() + *i);
  }

  if (msg_LevelIsDebugging()) {
    msg_Out() << "Will use clustered process for EWSud: ";
    for (const auto& flav : flavs) {
      msg_Out() << flav.ShellName() << " ";
    }
    msg_Out() << '\n';
  }

  assert(calculators.find(flavs) != calculators.end());
  EWSudakov_Log_Corrections_Map CorrectionsMaps {calculators[flavs]->CorrectionsMap(mom)};
  CorrectionsMaps[EWSudakov_Log_Type::lI] = ClusteredIOperator;
  return CorrectionsMaps;
}
