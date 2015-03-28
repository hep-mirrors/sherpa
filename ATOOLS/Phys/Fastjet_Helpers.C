#include "ATOOLS/Phys/Fastjet_Helpers.H"
#include "ATOOLS/Org/Message.H"
#ifdef USING__FASTJET

namespace ATOOLS {

  bool BTag(const fastjet::PseudoJet& jet, int bmode)
  {
    if (bmode==0) return false;

#ifdef USING__FASTJET__3
    int nb=0;
    std::vector<fastjet::PseudoJet> cons = jet.constituents();
    for (size_t i=0; i<cons.size(); ++i) {
      if (cons[i].user_index()==5) ++nb;
      if (cons[i].user_index()==-5) {
        if (bmode==1) ++nb;
        else if (bmode==2) --nb;
      }
    }
    return (nb!=0);
#else
    return false;
#endif
  }


  bool ToBeClustered(const ATOOLS::Flavour& flav, int bmode)
  {
    return
      (bmode==0 && Flavour(kf_jet).Includes(flav)) || 
      (bmode>0 && (Flavour(kf_jet).Includes(flav) || flav.Kfcode()==kf_b));
  }


  fastjet::PseudoJet MakePseudoJet(const Flavour& flav, const Vec4D& mom)
  {
    fastjet::PseudoJet ret(mom[1],mom[2],mom[3],mom[0]);
    ret.set_user_index((long int) flav);
    return ret;
  }

}

#endif
