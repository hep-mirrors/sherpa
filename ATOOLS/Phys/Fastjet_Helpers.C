#include "ATOOLS/Phys/Fastjet_Helpers.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/MyStrStream.H"

namespace ATOOLS {

  bool ToBeClustered(const ATOOLS::Flavour& flav, int bmode)
  {
    return
      (bmode==0 && Flavour(kf_jet).Includes(flav)) ||
      (bmode>0 && (Flavour(kf_jet).Includes(flav) || flav.Kfcode()==kf_b));
  }

  bool ToBeClustered(const ATOOLS::Flavour& flav,
                     const ATOOLS::Vec4D& mom,
                     const ATOOLS::Jet_Inputs& jetinput)
  {
    for (size_t i(0);i<jetinput.size();++i)
      if (jetinput[i].m_fl.Includes(flav) &&
          (jetinput[i].m_ptmin==0.      || mom.PPerp()>jetinput[i].m_ptmin) &&
          (jetinput[i].m_etamax>100000. || mom.Eta()<jetinput[i].m_etamax) &&
          (jetinput[i].m_ymax>100000.   || mom.Y()<jetinput[i].m_ymax))
        return true;
    return false;
  }

  fjcore::PseudoJet MakePseudoJet(const ATOOLS::Flavour& flav,
                                   const Vec4D& mom)
  {
    fjcore::PseudoJet ret(mom[1],mom[2],mom[3],mom[0]);
    ret.set_user_index((long int)flav);
    return ret;
  }

  bool BTag(const fjcore::PseudoJet& jet, int bmode)
  {
    if (bmode==0) return false;

    int nb=0;
    std::vector<fjcore::PseudoJet> cons = jet.constituents();
    for (size_t i=0; i<cons.size(); ++i) {
      if (cons[i].user_index()==5) ++nb;
      if (cons[i].user_index()==-5) {
        if (bmode==1) ++nb;
        else if (bmode==2) --nb;
      }
    }
    return (nb!=0);
  }

  Flavour FlavourTag(const fjcore::PseudoJet& jet,
                     const Jet_Identifications& jetids,
                     const kf_code& notagkf)
  {
    DEBUG_FUNC(jet.constituents().size()<<" constituents");
    std::vector<fjcore::PseudoJet> cons = jet.constituents();
    if (msg_LevelIsDebugging()) {
      for (size_t i(0);i<cons.size();++i)
        msg_Out()<<cons[i].user_index()<<": "
                 <<"("<<cons[i].E()<<","<<cons[i].px()
                 <<","<<cons[i].py()<<","<<cons[i].pz()<<")"<<std::endl;
    }
    for (size_t i(0);i<jetids.size();++i) {
      Vec4D momid(0.,0.,0.,0.);
      for (size_t j(0);j<cons.size();++j) {
        msg_Debugging()<<cons[j].user_index()<<" <-> "
                       <<(long int)jetids[i]->Flavour()<<std::endl;
        if (cons[j].user_index()==(long int)jetids[i]->Flavour()) {
          momid+=Vec4D(cons[j].E(),cons[j].px(),cons[j].py(),cons[j].pz());
        }
        msg_Debugging()<<"momid: "<<momid<<std::endl;
        msg_Debugging()<<" E="<<momid[0]<<" => "
                       <<momid[0]/jet.E()<<" <-> "
                       <<jetids[i]->EMin()<<std::endl;
        msg_Debugging()<<"ET="<<momid.MPerp()<<" => "
                       <<momid.MPerp()/jet.Et()<<" <-> "
                       <<jetids[i]->ETMin()<<std::endl;
        msg_Debugging()<<"pT="<<momid.PPerp()<<" => "
                       <<momid.PPerp()/jet.perp()<<" <-> "
                       <<jetids[i]->PTMin()<<std::endl;
      }
      if (jetids[i]->JetIDMode()&JetIdMode::larger) {
        if (jetids[i]->JetIDMode()&JetIdMode::relative) {
          if ((jetids[i]->EMin()>0. &&
               momid[0]/jet.E()>jetids[i]->EMin()) ||
              (jetids[i]->ETMin()>0. &&
               momid.MPerp()/jet.Et()>jetids[i]->ETMin()) ||
              (jetids[i]->PTMin()>0. &&
               momid.PPerp()/jet.perp()>jetids[i]->PTMin())) {
            msg_Debugging()<<"identified as "<<jetids[i]->Flavour()<<std::endl;
            return jetids[i]->Flavour();
          }
        }
        else {
          if ((jetids[i]->EMin()>0. &&
               momid[0]>jetids[i]->EMin()) ||
              (jetids[i]->ETMin()>0. &&
               momid.MPerp()>jetids[i]->ETMin()) ||
              (jetids[i]->PTMin()>0. &&
               momid.PPerp()>jetids[i]->PTMin())) {
            msg_Debugging()<<"identified as "<<jetids[i]->Flavour()<<std::endl;
            return jetids[i]->Flavour();
          }
        }
      }
      else THROW(not_implemented,"Relation not implemented yet.");
    }
    msg_Debugging()<<"not identified"<<std::endl;
    return Flavour(notagkf);
  }


  std::ostream &operator<<(std::ostream &str,const Jet_Input &ji)
  {
    std::string out("");
    if (ji.m_ptmin>0.)     out+="pT>"+ToString(ji.m_ptmin)+",";
    if (ji.m_etamax<1000.) out+="|eta|<"+ToString(ji.m_etamax)+",";
    if (ji.m_ymax<1000.)   out+="|y|<"+ToString(ji.m_ymax);
    if (out[out.size()-1]==',') out.substr(0,out.size()-1);
    return str<<ji.m_fl<<((out!="")?"["+out+"]":"");
  }

  std::ostream &operator<<(std::ostream &str,const JetIdMode::code &idm)
  {
    std::string out="";
    if (idm&JetIdMode::larger)   out+=">";
    else                         out+="<";
    if (idm&JetIdMode::relative) out+="[rel]";
    else                         out+="[abs]";
    return str<<out;
  }

  std::ostream &operator<<(std::ostream &str,const Jet_Identification &jid)
  {
    return str<<"Jet_Identification: "<<jid.Flavour()
              <<", (pT,ET,E)"<<jid.JetIDMode()
              <<"("<<jid.PTMin()<<","<<jid.ETMin()<<","<<jid.EMin()<<")";
  }
}

using namespace ATOOLS;

Jet_Identification::Jet_Identification(const ATOOLS::Flavour& flid,
                                       const double& ptmin, const double& etmin,
                                       const double& emin,
                                       const JetIdMode::code& jetidmode) :
  m_flid(flid), m_ptmin(ptmin), m_etmin(etmin), m_emin(emin),
  m_jetidmode(jetidmode)
{
  DEBUG_FUNC("");
}

