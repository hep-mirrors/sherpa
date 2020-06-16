#include "PHASIC++/Selectors/Fastjet_Selector_Base.H"

namespace PHASIC {
  class Fastjet_Selector: public Fastjet_Selector_Base, public ATOOLS::Tag_Replacer {
    int m_bmode;
    ATOOLS::Algebra_Interpreter m_calc;
    ATOOLS::Vec4D_Vector m_p;
    std::vector<double> m_mu2;

    std::string ReplaceTags(std::string &expr) const;
    ATOOLS::Term *ReplaceTags(ATOOLS::Term *term) const;
    void AssignId(ATOOLS::Term *term);

  public:
    Fastjet_Selector(Process_Base* const proc,
                     ATOOLS::Scoped_Settings s,
                     int bmode,
                     const std::string& expression);
    bool Trigger(ATOOLS::Selector_List&);
    void BuildCuts(Cut_Data*);
  };
}

using namespace PHASIC;
using namespace ATOOLS;


/*---------------------------------------------------------------------

  General form - flavours etc are unknown, will operate on a Particle_List.

  --------------------------------------------------------------------- */

Fastjet_Selector::Fastjet_Selector(Process_Base* const proc,
                                   Scoped_Settings s,
                                   int bmode,
                                   const std::string& expression):
  Fastjet_Selector_Base("FastjetSelector", proc, s),
  m_bmode(bmode)
{
  m_p.resize(m_nin+m_nout);
  m_mu2.resize(m_nout);

  m_calc.AddTag("H_T2","1.0");
  m_calc.AddTag("P_SUM","(1.0,0.0,0.0,0.0)");
  for (size_t i=0;i<m_p.size();++i)
    m_calc.AddTag("p["+ToString(i)+"]",ToString(m_p[i]));
  for (size_t i=0;i<m_mu2.size();++i)
    m_calc.AddTag("MU_"+ToString(i)+"2",ToString(m_mu2[i]));

  m_calc.SetTagReplacer(this);
  m_calc.Interprete(expression);

  msg_Debugging()<<METHOD<<"(): '"<<expression<<"' {\n";
  msg_Indent();
  if (msg_LevelIsDebugging()) m_calc.PrintEquation();
  msg_Debugging()<<"}\n";
}

void Fastjet_Selector::BuildCuts(Cut_Data *cuts)
{
  if (m_isnlo) return;
  if (m_nj==m_nout) { // number of jets equals number of outgoing particles
    double sumM2=0.;
    for (int i=m_nin;i<m_n;i++) {
      sumM2+=sqr(p_fl[i].SelMass());
    }

    for (int i=m_nin;i<m_n;i++) {
      // PT cut
      cuts->energymin[i] = Max(sqrt(sqr(m_ptmin)+sqr(p_fl[i].SelMass())),
                               cuts->energymin[i]);
      double Emax2 = sqr((m_smax+2.*sqr(p_fl[i].SelMass())-sumM2)
                       /(2.*sqrt(m_smax)));
      double cosmax = Min(cuts->cosmax[0][i],
                          sqrt(1.-sqr(m_ptmin)/(Emax2-sqr(p_fl[i].SelMass()))));
      cuts->cosmax[0][i] = cuts->cosmax[1][i] = cosmax;
      cuts->cosmax[i][0] = cuts->cosmax[i][1] = cosmax;
      cuts->etmin[i] = Max(sqrt(sqr(m_ptmin)+sqr(p_fl[i].SelMass())
                           *(1.-sqr(cuts->cosmax[0][i]))),cuts->etmin[i]);

      // ET cut
      cuts->energymin[i] = Max(sqrt(sqr(m_etmin)+sqr(p_fl[i].SelMass())),
                               cuts->energymin[i]);
      cosmax = Min(cuts->cosmax[0][i],
                          sqrt(1.-sqr(m_etmin)/(Emax2-sqr(p_fl[i].SelMass()))));
      cuts->cosmax[0][i] = cuts->cosmax[1][i] = cosmax;
      cuts->cosmax[i][0] = cuts->cosmax[i][1] = cosmax;
      cuts->etmin[i] = Max(sqrt(sqr(m_etmin)+sqr(p_fl[i].SelMass())
                           *(1.-sqr(cuts->cosmax[0][i]))),cuts->etmin[i]);

      // Eta cut
      cuts->cosmin[1][i] = cuts->cosmin[i][1]
        = Max(cuts->cosmin[1][i],tanh(-m_eta));
      cuts->cosmin[0][i] = cuts->cosmin[i][0]
        = Max(cuts->cosmin[0][i],tanh(-m_eta));
      cuts->cosmax[0][i] = cuts->cosmax[i][0]
        = Min(cuts->cosmax[0][i],tanh(m_eta));
      cuts->cosmax[1][i] = cuts->cosmax[i][1]
        = Min(cuts->cosmax[1][i],tanh(m_eta));

      // Y cut
      cuts->cosmax[0][i] = cuts->cosmax[i][0] =
        Min(cuts->cosmax[0][i],tanh(m_y)/sqrt(1.-sqr(p_fl[i].SelMass())
                                                 /sqr(cuts->energymin[i])));
      cuts->cosmax[1][i] = cuts->cosmax[i][1] =
        Min(cuts->cosmax[1][i],tanh(m_y)/sqrt(1.-sqr(p_fl[i].SelMass())
                                                  /sqr(cuts->energymin[i])));
    }
  }
}

std::string Fastjet_Selector::ReplaceTags(std::string &expr) const
{
  return m_calc.ReplaceTags(expr);
}

Term *Fastjet_Selector::ReplaceTags(Term *term) const
{
  if (term->Id()>=1000) {
    term->Set(m_mu2[term->Id()-1000]);
    return term;
  }
  if (term->Id()>=100) {
    term->Set(m_p[term->Id()-100]);
    return term;
  }
  else if (term->Id()==5) {
    double ht(0.0);
    for (size_t i(0);i<m_p.size();++i) ht+=m_p[i].PPerp();
    term->Set(sqr(ht));
    return term;
  }
  else if (term->Id()==6) {
    Vec4D sum(0.0,0.0,0.0,0.0);
    for (size_t i(0);i<m_p.size();++i) sum+=m_p[i];
    term->Set(sum);
    return term;
  }
  return term;
}

void Fastjet_Selector::AssignId(Term *term)
{
  if (term->Tag()=="H_T2") term->SetId(5);
  else if (term->Tag()=="P_SUM") term->SetId(6);
  else if (term->Tag().find("MU_")==0) {
    int idx(ToType<int>(term->Tag().substr(3,term->Tag().length()-4)));
    if (idx>=m_mu2.size()) THROW(fatal_error,"Invalid syntax");
    term->SetId(1000+idx);
  }
  else {
    int idx(ToType<int>(term->Tag().substr(2,term->Tag().length()-3)));
    if (idx>=m_nin+m_nout) THROW(fatal_error,"Invalid syntax");
    term->SetId(100+idx);
  }
}

bool Fastjet_Selector::Trigger(Selector_List &sl)
{
  if (m_nj<0) return true;

  m_p.clear();
  for (size_t i(0);i<m_nin;++i) m_p.push_back(sl[i].Momentum());
  std::vector<fjcore::PseudoJet> input,jets;
  for (size_t i(m_nin);i<sl.size();++i) {
    if (ToBeClustered(sl[i].Flavour(), m_bmode)) {
      input.push_back(MakePseudoJet(sl[i].Flavour(),sl[i].Momentum()));
    } else {
      m_p.push_back(sl[i].Momentum());
    }
  }
  int nj=m_p.size();

  fjcore::ClusterSequence cs(input,*p_jdef);
  jets=fjcore::sorted_by_pt(cs.inclusive_jets());

  if (m_eekt) {
    for (size_t i(0);i<input.size();++i) {
      if (cs.exclusive_dmerge_max(i)>sqr(m_ptmin)) {
        m_p.emplace_back(jets[i].E(),jets[i].px(),jets[i].py(),jets[i].pz());
      }
    }
  } else {
    for (size_t i(0);i<jets.size();++i) {
      if (m_bmode==0 || BTag(jets[i], m_bmode)) {
        Vec4D pj(jets[i].E(),jets[i].px(),jets[i].py(),jets[i].pz());
        if (pj.PPerp() > m_ptmin
            && pj.EPerp() > m_etmin
            && dabs(pj.Eta()) < m_eta
            && dabs(pj.Y()) < m_y)
          m_p.push_back(pj);
      }
    }
  }
  for (size_t i(0);i<input.size();++i)
    m_mu2[i]=cs.exclusive_dmerge_max(i);

  bool trigger((int)(m_p.size()-nj)>=m_nj);
  if (trigger) trigger=(int)m_calc.Calculate()->Get<double>();

  return (1-m_sel_log->Hit(1-trigger));
}


DECLARE_GETTER(Fastjet_Selector,"FastjetSelector",Selector_Base,Selector_Key);

Selector_Base *ATOOLS::Getter<Selector_Base,Selector_Key,Fastjet_Selector>::
operator()(const Selector_Key &key) const
{
  auto s = key.m_settings["FastjetSelector"];

  const auto expression = s["Expression"].SetDefault("")  .Get<std::string>();
  const auto bmode      = s["BMode"]     .SetDefault(0)   .Get<int>();

  return new Fastjet_Selector(key.p_proc, s, bmode, expression);
}

void ATOOLS::Getter<Selector_Base,Selector_Key,Fastjet_Selector>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"FastjetSelector:\n"
     <<width<<"  Expression: boolean expression\n";
  Fastjet_Selector_Base::PrintCommonInfoLines(str, width);
  str<<width<<"  BMode: 0|1|2";
}
