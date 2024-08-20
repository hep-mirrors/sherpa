#include "PHASIC++/Selectors/Selector.H"
#include "ATOOLS/Phys/Fastjet_Helpers.H"

#include "PHASIC++/Main/Process_Integrator.H"
#include "ATOOLS/Org/Run_Parameter.H"

#include <numeric>

/*--------------------------------------------------------------------

 NJettiness arXiv:1004.2489

 --------------------------------------------------------------------*/

using namespace ATOOLS;

namespace PHASIC {
  class NJettiness_Selector: public Selector_Base {
  private:
    enum class Frame{
      Undefined = -1,
      Laboratory = 0,
      CentreOfMass,
    } m_frame;

    enum class Norm{
      Undefined = -1,
      JetEnergy,
      JetEnergy2,
      CentreOfMass
    } m_norm;

    size_t m_N;
    double m_taucut, m_delta_r, m_ptmin, m_etmin, m_etamax, m_ymax;
    fjcore::RecombinationScheme m_reco;
    fjcore::JetAlgorithm m_algo;

    inline void SetFrame(const std::string& frame);
    inline void SetNorm(const std::string& norm);

    static fjcore::JetAlgorithm GetAlgorithm(const std::string& algo);
  public:
    NJettiness_Selector(Process_Base* const proc,
                          ATOOLS::Scoped_Settings s);
    ~NJettiness_Selector() = default;

  protected:
  public:
    static void PrintCommonInfoLines(std::ostream& str, size_t width);

    bool Trigger(ATOOLS::Selector_List &sl) override;
    void BuildCuts(Cut_Data *) override;
  };

  NJettiness_Selector::NJettiness_Selector(Process_Base* const proc,
                                           ATOOLS::Scoped_Settings s)
    : Selector_Base("",proc)
    {
      DEBUG_INFO(METHOD);
      rpa->gen.AddCitation(1, "N-Jettiness is published under \\cite{Stewart:2010tn,Stewart:2015waa}");
      /// Operates on a fixed number of hard directions, N
      m_N = s["N"].SetDefault(0).Get<size_t>();

      /// Passes if the value of N-jettiness is above a given cut
      m_taucut = s["TauN"].SetDefault(1.).Get<double>();
      const auto frame = s["Frame"]        .SetDefault("Lab").Get<std::string>();
      SetFrame(frame);
      const auto norm  = s["Normalization"].SetDefault("jet2") .Get<std::string>();
      SetNorm(norm);

      /// additional parameters might be frame, which underlying algo to
      /// determine the hard-directions... for the moment use anti-kt in
      /// in the lab frame
      const auto algo = s["Algorithm"]          .SetDefault("antikt")  .Get<std::string>();
      const auto reco = s["RecombinationScheme"].SetDefault("E")       .Get<std::string>();
      m_reco = ATOOLS::GetRecombinationScheme(reco);

      m_delta_r = s["DR"]    .SetDefault(1)                                .Get<double>();
      // min/max settings
      m_ptmin   = s["PTMin"] .SetDefault("None").UseZeroReplacements()     .Get<double>();
      m_etmin   = s["ETMin"] .SetDefault("None").UseZeroReplacements()     .Get<double>();
      m_etamax  = s["EtaMax"].SetDefault("None").UseMaxDoubleReplacements().Get<double>();
      m_ymax    = s["YMax"]  .SetDefault("None").UseMaxDoubleReplacements().Get<double>();


    }

  bool NJettiness_Selector::Trigger(ATOOLS::Selector_List &sl)
  {
    bool trigger = false;
    DEBUG_FUNC(p_proc->Flavours() << " " << sl.ExtractMomenta());

    auto mom = sl.ExtractMomenta();

    switch (m_frame) {
    case Frame::Laboratory:
    default:
      break;
    case Frame::CentreOfMass:
      {
        Poincare cms(mom[0]+mom[1]);
        for(auto& m: mom){
          cms.Boost(m);
        }
      }
      break;
    }

    /// fastjet clustering of fs particles
    std::vector<fjcore::PseudoJet> input, jets;
    /// prepare input partons
    for(size_t i{m_nin}; i < sl.size();++i){
      if(Flavour(kf_jet).Includes(sl[i].Flavour())){
        input.push_back(MakePseudoJet(sl[i].Flavour(),mom[i]));
      }
    }
    if(m_N > input.size()){
      THROW(fatal_error, METHOD + " Cannot have more hard directions than partons");
    }

    if(m_N > 0){
      /// create cluster sequence and get the hard fs directions (these need to be added to beam)
      /// in a hadron-collider setup
      fjcore::ClusterSequence cs(input,
                                 fjcore::JetDefinition(m_algo,m_delta_r, m_reco));
      jets = cs.exclusive_jets(static_cast<int>(m_N));

      /// Check user selection for hard directions
      for(auto&j : jets){
        if(j.pt() < m_ptmin)                       trigger = false;
        if(sqrt(j.m2() + j.pt()*j.pt()) < m_etmin) trigger = false;
        if(j.eta() > m_etamax)                     trigger = false;
        if(j.rapidity() > m_ymax)                  trigger = false;
      }
    }

    /// add beams if necessary
    for(size_t i{0}; i< m_nin; ++i){
      if(Flavour(kf_jet).Includes(sl[i].Flavour())) {
        jets.push_back(MakePseudoJet(sl[i].Flavour(), mom[i]));
      }
    }

    /// compute TauN, i.e. the sum of the minimum for each FS partons w.r.t to the hard directions
    std::vector<double> TauN;
    for(auto& i: input){
      std::vector<double> TauN_i;
      for(auto& j : jets){
        double normalization = 1;
        switch (m_norm) {
        case Norm::Undefined:
          THROW(fatal_error, METHOD + " undefined normalisation scheme for N-Jettiness");
          break;
        case Norm::JetEnergy:
          normalization = j.E();
          break;
        case Norm::JetEnergy2:
          normalization = 2.*j.E();
        break;
        case Norm::CentreOfMass:
          normalization = 2. * sqrt(mom[0].E() * mom[1].E());
        }
        TauN_i.push_back(2.*fjcore::dot_product(j, i)/normalization); /// there are of course different choices here one can make
      }
      TauN.push_back(*std::min_element(TauN_i.begin(),TauN_i.end()));
    }
    const double TauN_Val = std::accumulate(TauN.begin(), TauN.end(), 0.0);

    trigger = (TauN_Val > m_taucut);

    if(!trigger){
      msg_Debugging() << "Point discarded by " << METHOD << std::endl;
    } else {
      msg_Debugging() << "Point passed: " << METHOD << std::endl;
    }

    return (1 - m_sel_log->Hit(1 - trigger));
  }

  fjcore::JetAlgorithm NJettiness_Selector::GetAlgorithm(const std::string& algo)
  {
    const bool leplep(rpa->gen.Bunch(0).IsLepton() and rpa->gen.Bunch(1).IsLepton());
    if(algo == "kt"){
      if(leplep) return fjcore::ee_kt_algorithm;
      else return fjcore::kt_algorithm;
    } else if(algo == "antikt") {
      return fjcore::antikt_algorithm;
    } else if(algo == "genkt"){
      if(leplep) return fjcore::ee_genkt_algorithm;
      else return fjcore::genkt_algorithm;
    } else if(algo == "cambridge") {
      return fjcore::cambridge_algorithm;
    } else return fjcore::undefined_jet_algorithm;
  }

  inline void NJettiness_Selector::SetFrame(const std::string& frame)
  {
    m_frame = Frame::Undefined;
    if(frame       == "lab" or frame == "laboratory")    m_frame = Frame::Laboratory;
    else if (frame == "com" or frame == "centreofmass")  m_frame = Frame::CentreOfMass;
    return;
  }

  inline void NJettiness_Selector::SetNorm(const std::string& norm)
  {
    m_norm = Norm::Undefined;
    if(norm       == "jet"     or norm == "jetenergy")     m_norm = Norm::JetEnergy;
    else if (norm == "jet2"    or norm == "jetenergy2")    m_norm = Norm::JetEnergy2;
    else if (norm == "com"     or norm == "centreofmass")  m_norm = Norm::CentreOfMass;
    return;
  }

  void NJettiness_Selector::BuildCuts(Cut_Data *) {return;}
}

using namespace PHASIC;

DECLARE_GETTER(NJettiness_Selector,"NJettiness",Selector_Base,Selector_Key);

Selector_Base *ATOOLS::Getter<Selector_Base,Selector_Key,NJettiness_Selector>::
operator()(const Selector_Key &key) const
{
  auto s = key.m_settings["NJettiness"];

  return new NJettiness_Selector(key.p_proc, s);
}

void ATOOLS::Getter<Selector_Base,Selector_Key,NJettiness_Selector>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"NJettiness:\n";
  str<< width << "  N: number of hard directions\n"
     << width << "  TauN: NJettiness-cut\n"
     << width << "  Frame: frame in which to evaluate NJettiness: lab, com (default: lab)"
     << width << "  Normalization: normalization for the jettiness contribution: jet, jet2, com (default: jet2)"
     << width << "  Algorithm: antikt, kt, genkt, cambridge\n"
     << width << "  RecombinationScheme: E, pt, pt2, Et, Et2, BIpt, BIpt2 (default: E) \n"
     << width << "  PTMin: minimum jet pT for each hard direction\n"
     << width << "  ETMin: minimum jet eta for each hard direction\n"
     << width << "  DR: jet distance parameter (default: 1.)\n"
     << width << "  EtaMax: maximum jet eta for each hard direction (default: None)\n"
     << width << "  YMax: maximum jet rapidity for each hard direction (default: None)\n";
  str<<"cite: arXiv:1004.2489";
}
