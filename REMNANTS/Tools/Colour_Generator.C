#include "REMNANTS/Tools/Colour_Generator.H"
#include "REMNANTS/Main/Remnant_Handler.H"
#include "REMNANTS/Main/Remnant_Base.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Message.H"
#include <list>
#include <set>

using namespace REMNANTS;
using namespace ATOOLS;
using namespace std;

    
Colour_Generator::Colour_Generator() {
  for(size_t beam=0;beam<2;beam++) p_remnants[beam] = NULL;
}

Colour_Generator::~Colour_Generator() {}

void Colour_Generator::Initialize(Remnant_Handler * rhandler) {
  for(size_t beam=0;beam<2;beam++) p_remnants[beam] = rhandler->GetRemnant(beam);
}

void Colour_Generator::ConnectColours(Blob *const showerblob) {
  // Extract incoming particles/shower initiators - they are the ones that
  // matter for the colour handling.
  for (size_t i=0;i<showerblob->NInP();i++) {
    Particle * part = showerblob->InParticle(i);
    size_t beam = part->Beam();
    if (beam==-1 || part->ProductionBlob()!=NULL) continue;
    p_inparts[beam] = part;
  }
  // First replace t-channel connected colours in IS.  Without t-channel colour
  // exchange, go for the unconnected colour flows of both beams.
  // The logic of the replacements is detailed in the methods.
  if (!TChannelColourFlows() && !SChannelColourFlows()) {
    msg_Error()<<"Warning in "<<METHOD<<": No colours in incoming partons.\n"
	       <<(*showerblob)<<"\n";
    Output();
  }
  if (!showerblob->CheckColour(true)) {
    msg_Out()<<METHOD<<" did not conserved colour in:\n"<<(*showerblob)<<"\n";
    Output();
    exit(1);
  }
  showerblob->UnsetStatus(blob_status::needs_beams);
  ResetFlags();
}

bool Colour_Generator::TChannelColourFlows() {
  // Do nothing if no t-channel exchange
  if ((p_inparts[0]->GetFlow(1)==0 ||
       (p_inparts[0]->GetFlow(1)!=p_inparts[1]->GetFlow(2)))  &&
      (p_inparts[1]->GetFlow(1)==0 ||
       (p_inparts[0]->GetFlow(2)!=p_inparts[1]->GetFlow(1)))) {
    return false;
  }
  // Do nothing if t-channel exchange but no replacement colours
  for (size_t pos=0;pos<2;pos++) {
    if (p_inparts[pos]->GetFlow(1)!=0 &&
	p_inparts[pos]->GetFlow(1)==p_inparts[1-pos]->GetFlow(2) &&
	m_cols[pos][0].empty() && m_cols[1-pos][1].empty()) {
      return false;
    }
  }
  // Logic is as follows:
  // For sea-quarks and their spectators or gluons, we have two colours -
  // triplet and anti-triplet - in each beam, one pair of triplet-anti-triplet is
  // connected through the t-channel.  We define a donor of colour, replacing either
  // the triplet or anti-triplet colour in the link, taking it from the stack in
  // one of the two beams, selected by DefineColourDonor.  If the recipient particle
  // is a valence quark, its beam will become the colour donor.  In both cases, the
  // replacement colour is put back into the stack of the recipient beam, and the
  // replaced colour is erased from both stacks.
  // Find t-channel colour flow and define beam to take replacement colour 
  int tbeam = (p_inparts[0]->GetFlow(1)!=0 &&
	       p_inparts[0]->GetFlow(1)==p_inparts[1]->GetFlow(2))?0:1;
  int donor = DefineColourDonor(tbeam);
  if (donor==-1) return true;
  // Make sure the valence side is the colour donor.
  Particle * spectator = p_remnants[1-donor]->GetSpectator();
  if (!spectator && !p_inparts[1-donor]->Flav().IsGluon()) {
    donor     = 1-donor;
    spectator = p_remnants[1-donor]->GetSpectator();
  }
  if (donor==tbeam) {
    // replace triplet colour in donor (and as anti-triplet in recipient)
    // and push new colour into recipient beam
    ReplaceBoth(donor,0);
    if (spectator || p_inparts[1-donor]->Flav().IsGluon())
      Replace(1-donor,0,spectator?spectator:p_inparts[1-donor]);
  }
  else {
    // replace anti-triplet colour (and as triplet in recipient)
    // and push new colour into recipient beam
    ReplaceBoth(donor,1);
    if (spectator || p_inparts[1-donor]->Flav().IsGluon())
      Replace(1-donor,1,spectator?spectator:p_inparts[1-donor]);
  }
  return true;
}

int Colour_Generator::DefineColourDonor(const size_t & tbeam) {
  // Define the stack to take the new t-channel colour from, according to availability.
  bool tripcol(false), anticol(false);
  for (list<int>::iterator cit=m_cols[tbeam][0].begin();
       cit!=m_cols[tbeam][0].end();cit++) {
    if (m_vetoed[tbeam][0].find(*cit)==m_vetoed[tbeam][0].end()) {
      tripcol = true;
      break;
    }
  }
  for (list<int>::iterator cit=m_cols[1-tbeam][1].begin();
       cit!=m_cols[1-tbeam][1].end();cit++) {
    if (m_vetoed[1-tbeam][1].find(*cit)==m_vetoed[1-tbeam][1].end()) {
      anticol = true;
      break;
    }
  }
  // Make sure - by direct veto - that no singlet gluons are produced by picking a colour
  // from stack of the beam that is already in the stack of the recipient.
  // It is possible that this means that the t-channel colour cannot be replaced, and that
  // colour replacements are delegated to the SChannelColourFlow method.
  if (tripcol && p_inparts[1-tbeam]->Flav().IsGluon()) {
    if (m_cols[1-tbeam][0].front()==m_cols[tbeam][0].front()) {
      tripcol = false;
    }
  }
  if (anticol && p_inparts[tbeam]->Flav().IsGluon()) {
    if (m_cols[tbeam][1].front()==m_cols[1-tbeam][1].front()) {
      anticol = false;
    }
  }
  
  if (tripcol && anticol) return (ran->Get()>0.5?tbeam:1-tbeam);
  if (tripcol)            return tbeam;
  if (anticol)            return 1-tbeam;
  return -1;
}

bool Colour_Generator::SChannelColourFlows() {
  // Check if there is any danger to draw identical colours from both stacks.  This
  // could lead to colour sinlget gluons in the final state and must be avoided.
  for (size_t beam=0;beam<2;beam++) {
    if (m_cols[beam][0].empty() || m_cols[1-beam][1].empty()) continue;
    if (m_cols[beam][0].front()==m_cols[1-beam][1].front()) {
      //msg_Out()<<"Potential problem in "<<METHOD<<": identical colour "
      //       <<"("<<m_cols[beam][0].front()<<") on stack for "
      //       <<p_inparts[0]->Flav()<<" + "<<p_inparts[1]->Flav()<<".\n";
      //Output();
      return ConstrainedColourFlows(beam);
    }
  }
  // Overall logic is the following: unless we extract a valence quark, we assume
  // that extracted particles come as colour-octets: either as gluon or as sea-quark
  // plus a spectator.  We then replace one colour index - say triplet - of the system
  // with one from the stack of previous colours and put the new anti-triplet colour
  // of the system on the triplet stack.  This is achieved in in AssignColours.
  // For valence quarks we only replacve their colour with one from the stack.
  // flag to make sure that we do not add a colour from a particle to replace
  // the anti-colour.
  bool replacecol(false);
  for (size_t beam=0;beam<2;beam++) {
    Particle * spectator = p_remnants[beam]->GetSpectator();
    if (spectator) {
      if (p_inparts[beam]->Flav().IsQuark() && p_inparts[beam]->Flav().IsAnti()) {
	AssignColours(beam,spectator,p_inparts[beam]);
      }
      else {
	AssignColours(beam,p_inparts[beam],spectator);
      }
      replacecol = true;
    }
    else if (p_inparts[beam]->Flav().IsGluon()) {
      AssignColours(beam,p_inparts[beam],p_inparts[beam]);
      replacecol = true;
    }
    else if (p_inparts[beam]->Flav().IsQuark()) {
      Replace(beam,int(p_inparts[beam]->Flav().IsAnti()),p_inparts[beam]);
      replacecol = true;
    }
  }
  return replacecol;
}

bool Colour_Generator::ConstrainedGGFlows(const size_t & tbeam) {
  // Both beam particles are gluons - the logic is to replace, if possible, one of their
  // colours with a dangerous one (the common one on stack for both beams) and to
  // replace the colour of the other gluon with a harmless one.
  int ncolt = AvailableColours(tbeam), ncola = AvailableColours(1-tbeam), replace = 0;
  if (ncolt==3 && ncola==3) replace = ran->Get()>0.5?1:2;
  else if (ncolt==3)        replace = 1;
  else if (ncola==3)        replace = 2;
  int newcolt, oldcolt, newcola, oldcola;
  switch (replace) {
  case 2: 
    // harmless triplet colour from stack for gluon in anti-triplet beam and
    // potentially dangerous triplet colour fro stack for gluon from triplet beam
    oldcolt = p_inparts[tbeam]->GetFlow(1);
    oldcola = p_inparts[1-tbeam]->GetFlow(1);
    Replace(1-tbeam,0,p_inparts[1-tbeam]);
    Replace(tbeam,0,p_inparts[tbeam]);
    newcolt = p_inparts[tbeam]->GetFlow(1);
    newcola = p_inparts[1-tbeam]->GetFlow(1);
    break;
  case 1:
    // harmless anti-triplet colour from stack for gluon in triplet beam and
    // potentially dangerous anti-triplet colour for gluon from anti-triplet beam
    oldcolt = p_inparts[tbeam]->GetFlow(2);
    oldcola = p_inparts[1-tbeam]->GetFlow(2);
    Replace(1-tbeam,1,p_inparts[1-tbeam]);
    Replace(tbeam,1,p_inparts[tbeam]);
    newcolt = p_inparts[tbeam]->GetFlow(2);
    newcola = p_inparts[1-tbeam]->GetFlow(2);
    break;
  case 0:
  default:
    // potentially dangerous triplet and anti-triplet colours for both gluons,
    // and hope for the best.
    msg_Error()<<METHOD<<" tries potentially dangerous colour replacement - hope for the best.\n";
    oldcolt = p_inparts[tbeam]->GetFlow(1);
    oldcola = p_inparts[1-tbeam]->GetFlow(2);
    Replace(1-tbeam,1,p_inparts[1-tbeam]);
    Replace(tbeam,0,p_inparts[tbeam]);
    newcolt = p_inparts[tbeam]->GetFlow(1);
    newcola = p_inparts[1-tbeam]->GetFlow(2);
    break;    
  }
  //msg_Out()<<METHOD<<"(replace = "<<replace<<"): "
  //	   <<oldcolt<<" --> "<<newcolt<<" and "<<oldcola<<" --> "<<newcola<<".\n"
  //	   <<(*p_inparts[0]->DecayBlob());
  //Output();
  return true;
}

bool Colour_Generator::ConstrainedGQFlows(const size_t & tbeam) {
  // Pretty similar to the two gluon case.  Gluon on triplet beam, with potentially
  // dangerous colour, plus a quark on the anti-triplet beam
  int ncolt = AvailableColours(tbeam);
  int ncola = AvailableColours(1-tbeam);
  int newcolt, oldcolt, newcola, oldcola;
  bool anti = p_inparts[1-tbeam]->Flav().IsAnti();
  Particle * aspec = p_remnants[1-tbeam]->GetSpectator();
  //msg_Out()<<METHOD<<"(anti = "<<anti<<", aspec = ["<<aspec<<"], ncol_t = "<<ncolt<<")\n";
  if (!anti) {
    // Quark on anti-triplet beam is a quark, therefore it cannot replace its original colour
    // with a dangerous one, and the gluon can do the replacement with a tricky colour.
    oldcolt = p_inparts[tbeam]->GetFlow(1);
    Replace(tbeam,0,p_inparts[tbeam]);
    newcolt = p_inparts[tbeam]->GetFlow(1);
    if (ncola==1 || ncola==3) {
      // There is a colour on stack for the quark - take it.
      oldcola = p_inparts[1-tbeam]->GetFlow(1);
      Replace(1-tbeam,0,p_inparts[1-tbeam]);
      newcolt = p_inparts[1-tbeam]->GetFlow(1);
    }
    else if (aspec && ncola==2) {
      // There is nocolour on stack for the quark but a anti-colour for the potential spectator.
      oldcola = aspec->GetFlow(2);
      Replace(1-tbeam,1,aspec);
      newcolt = aspec->GetFlow(2);
    }
  }
  else if (anti) {
    // Quark on anti-triplet beam is an anti-quark, therefore we may be forced to replace it -
    // try to avoid this potentially dangoerous situation, by checking for anti-triplet colour
    // replacement either for the gluon or by moving the dangerous anti-triplet colour on a
    // potential spectator.
    if (ncolt==2 || ncolt==3) {
      // gluon on triplet replaces the anti-triplet colour
      oldcolt = p_inparts[tbeam]->GetFlow(2);
      Replace(tbeam,1,p_inparts[tbeam]);
      newcolt = p_inparts[tbeam]->GetFlow(2);
      oldcola = p_inparts[1-tbeam]->GetFlow(2);
      Replace(1-tbeam,1,p_inparts[1-tbeam]);
      newcola = p_inparts[1-tbeam]->GetFlow(2);
    }
    else if (ncolt==1 && aspec) {
      // gluon on triplet and spectator on anti-triplet replace triplet colour
      oldcolt = p_inparts[tbeam]->GetFlow(1);
      Replace(tbeam,0,p_inparts[tbeam]);
      newcolt = p_inparts[tbeam]->GetFlow(1);
      oldcola = aspec->GetFlow(1);
      Replace(1-tbeam,0,aspec);
      newcola = aspec->GetFlow(1);
    }
    else {
      // dangerous: both gluon on triplet and anti-quark on anti-triplet replace
      // potentially upsetting colour
      oldcolt = p_inparts[tbeam]->GetFlow(1);
      Replace(tbeam,0,p_inparts[tbeam]);
      newcolt = p_inparts[tbeam]->GetFlow(1);
      oldcola = p_inparts[1-tbeam]->GetFlow(2);
      Replace(1-tbeam,1,p_inparts[1-tbeam]);
      newcola = p_inparts[1-tbeam]->GetFlow(2);
    }
  }
  //msg_Out()<<METHOD<<"(ncol_a = "<<ncola<<" for anti = "<<anti<<", spec = "<<aspec<<"): "
  //	   <<oldcolt<<" --> "<<newcolt<<" and "<<oldcola<<" --> "<<newcola<<".\n"
  //	   <<(*p_inparts[0]->DecayBlob());
  //if (aspec) msg_Out()<<(*aspec)<<"\n";
  //Output();
  return true;
}

bool Colour_Generator::ConstrainedQGFlows(const size_t & tbeam) {
  // Pretty similar to the two gluon case.  Gluon on triplet beam, with potentially
  // dangerous colour, plus a quark on the anti-triplet beam
  int ncolt = AvailableColours(tbeam);
  int ncola = AvailableColours(1-tbeam);
  int newcolt, oldcolt, newcola, oldcola;
  bool anti = p_inparts[tbeam]->Flav().IsAnti();
  Particle * tspec = p_remnants[tbeam]->GetSpectator();
  //msg_Out()<<METHOD<<"(anti = "<<anti<<", tspec = ["<<tspec<<", ncol_a = "<<ncola<<"])\n";
  if (anti) {
    // Quark on triplet beam is an anti-quark, therefore it cannot replace its original colour
    // with a dangerous one, and the gluon can do the replacement with a tricky colour.
    oldcola = p_inparts[1-tbeam]->GetFlow(2);
    Replace(1-tbeam,1,p_inparts[1-tbeam]);
    newcola = p_inparts[1-tbeam]->GetFlow(2);
    if (ncolt==2 || ncolt==3) {
      // There is an anti-colour on stack for the anti-quark - take it.
      oldcolt = p_inparts[tbeam]->GetFlow(2);
      Replace(tbeam,1,p_inparts[tbeam]);
      newcolt = p_inparts[tbeam]->GetFlow(2);
    }
    else if (tspec && ncolt==1) {
      // There is no anti-colour on stack for the anti-quark but a colour for the potential spectator.
      oldcolt = tspec->GetFlow(1);
      Replace(tbeam,0,tspec);
      newcolt = tspec->GetFlow(1);
    }
  }
  else if (!anti) {
    // Quark on triplet beam is a quark, therefore we may be forced to replace it - try to avoid this
    // potentially dangoerous situation, by checking for triplet colour replacement either for the
    // gluon or by moving the dangerous triplet colour on a potential spectator.
    if (ncola==1 || ncola==3) {
      // quark on triplet and gluon on anti-triplet replaces the triplet colour
      oldcolt = p_inparts[tbeam]->GetFlow(1);
      Replace(tbeam,0,p_inparts[tbeam]);
      newcolt = p_inparts[tbeam]->GetFlow(1);
      oldcola = p_inparts[1-tbeam]->GetFlow(1);
      Replace(1-tbeam,0,p_inparts[1-tbeam]);
      newcola = p_inparts[1-tbeam]->GetFlow(1);
    }
    else if (tspec) {
      // spectator on triplet and gluon on anti-triplet replace the anti-triplet colour
      oldcolt = tspec->GetFlow(2);
      Replace(tbeam,1,tspec);
      newcolt = tspec->GetFlow(2);
      oldcola = p_inparts[1-tbeam]->GetFlow(2);
      Replace(1-tbeam,1,p_inparts[1-tbeam]);
      newcola = p_inparts[1-tbeam]->GetFlow(2);
    }
    else {
      // dangerous: both partons replace the potentially upsetting colour
      oldcolt = p_inparts[tbeam]->GetFlow(1);
      Replace(tbeam,0,p_inparts[tbeam]);
      newcolt = p_inparts[tbeam]->GetFlow(1);
      oldcola = p_inparts[1-tbeam]->GetFlow(2);
      Replace(1-tbeam,1,p_inparts[1-tbeam]);
      newcola = p_inparts[1-tbeam]->GetFlow(2);
    }
  }
  //msg_Out()<<METHOD<<"(ncol_a = "<<ncola<<" for anti = "<<anti<<", spec = "<<tspec<<"): "
  //	   <<oldcolt<<" --> "<<newcolt<<" and "<<oldcola<<" --> "<<newcola<<".\n"
  //	   <<(*p_inparts[0]->DecayBlob());
  //if (tspec) msg_Out()<<(*tspec)<<"\n";
  //Output();
  return true;
}

bool Colour_Generator::ConstrainedQQFlows(const size_t & tbeam) {
  int  ncolt = AvailableColours(tbeam), ncola = AvailableColours(1-tbeam);
  bool tanti = p_inparts[tbeam]->Flav().IsAnti();
  bool aanti = p_inparts[1-tbeam]->Flav().IsAnti();
  Particle * tspec = p_remnants[tbeam]->GetSpectator();
  Particle * aspec = p_remnants[1-tbeam]->GetSpectator();
  int newcolt = 0, oldcolt = 0, newcola = 0, oldcola = 0;
  bool replacecol = false;
  //msg_Out()<<METHOD<<" for\n"
  //	   <<"   triplet (tanti = "<<tanti<<", tspec = ["<<tspec<<"]) and "
  //	   <<"anti-triplet (aanti = "<<aanti<<", aspec = ["<<aspec<<"]).\n";
  if (tanti && aanti) {
    // harmless - both quarks are anti-quarks - the anti-quark on the anti-triplet beam
    // will get the potentially dangerous anti-triplet colour, the anti-quark on the triplet beam
    // will either get a colour from stack - if available - or not.
    oldcola = p_inparts[1-tbeam]->GetFlow(2);
    Replace(1-tbeam,1,p_inparts[1-tbeam]);
    newcola = p_inparts[1-tbeam]->GetFlow(2);
    if (ncolt==3) {
      oldcolt = p_inparts[tbeam]->GetFlow(2);
      Replace(tbeam,1,p_inparts[tbeam]);
      newcolt = p_inparts[tbeam]->GetFlow(2);
    }
    else if (tspec) {
      oldcolt = tspec->GetFlow(1);
      Replace(tbeam,0,tspec);
      newcolt = tspec->GetFlow(1);
    }
    replacecol = true;
  }
  else if (tanti && !aanti) {
    // harmless - both quarks can draw colours from stack without any problem
    if (ncolt==3) {
      oldcolt = p_inparts[tbeam]->GetFlow(2);
      Replace(tbeam,1,p_inparts[tbeam]);
      newcolt = p_inparts[tbeam]->GetFlow(2);
    }
    else if (tspec) {
      oldcolt = tspec->GetFlow(1);
      Replace(tbeam,0,tspec);
      newcolt = tspec->GetFlow(1);
    }
    if (ncola==3) {
      oldcola = p_inparts[1-tbeam]->GetFlow(1);
      Replace(1-tbeam,0,p_inparts[1-tbeam]);
      newcola = p_inparts[1-tbeam]->GetFlow(1);
    }
    else if (aspec) {
      oldcola = aspec->GetFlow(2);
      Replace(1-tbeam,1,aspec);
      newcola = aspec->GetFlow(2);
    }
    replacecol = true;
  }
  else if (!tanti && aanti) {
    // this is the dangerous one and it may lead to a situation where
    // we have to take out chances.
    int replace = 0;
    if (tspec && aspec)       replace = ran->Get()>0.5?1:2;
    else if (tspec && !aspec) replace = 1;
    else if (!tspec && aspec) replace = 2;
    //msg_Out()<<"--> replace = "<<replace<<"\n";
    if (replace==1) {
      m_vetoed[tbeam][1].erase(p_inparts[tbeam]->GetFlow(1));
      oldcolt = tspec->GetFlow(2);
      Replace(tbeam,1,tspec);
      newcolt = tspec->GetFlow(2);
      oldcola = p_inparts[1-tbeam]->GetFlow(2);
      Replace(1-tbeam,1,p_inparts[1-tbeam]);
      newcola = p_inparts[1-tbeam]->GetFlow(2);
    }
    else if (replace==2) {
      oldcolt = p_inparts[tbeam]->GetFlow(1);
      Replace(tbeam,0,p_inparts[tbeam]);
      newcolt = p_inparts[tbeam]->GetFlow(1);
      m_vetoed[1-tbeam][0].erase(p_inparts[1-tbeam]->GetFlow(2));
      oldcola = aspec->GetFlow(1);
      Replace(1-tbeam,0,aspec);
      newcola = aspec->GetFlow(1);
    }
    else {
      oldcolt = p_inparts[tbeam]->GetFlow(1);
      Replace(tbeam,0,p_inparts[tbeam]);
      newcolt = p_inparts[tbeam]->GetFlow(1);
      oldcola = p_inparts[1-tbeam]->GetFlow(2);
      Replace(1-tbeam,1,p_inparts[1-tbeam]);
      newcola = p_inparts[1-tbeam]->GetFlow(2);
    }
    replacecol = true;
  }
  else if (!tanti && !aanti) {
    // harmless - both quarks are anti-quarks - the anti-quark on the anti-triplet beam
    // will get the potentially dangerous anti-triplet colour, the anti-quark on the triplet beam
    // will either get a colour from stack - if available - or not.
    oldcolt = p_inparts[tbeam]->GetFlow(1);
    Replace(tbeam,0,p_inparts[tbeam]);
    newcolt = p_inparts[tbeam]->GetFlow(1);
    if (ncola==3) {
      oldcola = p_inparts[1-tbeam]->GetFlow(1);
      Replace(1-tbeam,0,p_inparts[1-tbeam]);
      newcola = p_inparts[1-tbeam]->GetFlow(1);
    }
    else if (aspec) {
      oldcola = aspec->GetFlow(2);
      Replace(1-tbeam,1,aspec);
      newcola = aspec->GetFlow(2);
    }
    replacecol = true;
  }
  //msg_Out()<<"--> (ncol_a = "<<ncola<<" for anti = "<<aanti<<", spec = "<<aspec<<") "
  //	   <<"(ncol_t = "<<ncolt<<" for trip = "<<tanti<<", spec = "<<tspec<<"): "
  //	   <<oldcolt<<" --> "<<newcolt<<" and "<<oldcola<<" --> "<<newcola<<".\n"
  //	   <<(*p_inparts[0]->DecayBlob())<<"\n"<<"Spectators:\n";
  //if (tspec) msg_Out()<<(*tspec)<<"\n";
  //if (aspec) msg_Out()<<(*aspec)<<"\n";
  //Output();
  if (replacecol) return true;
  exit(1);
  return true;
}

bool Colour_Generator::ConstrainedColourFlows(const size_t & tbeam) {
  ran->Get();
  // tbeam is beam particle with tricky triplet colour
  int abeam = 1-tbeam;
  Flavour tflav = p_inparts[tbeam]->Flav(), aflav = p_inparts[abeam]->Flav();
  if (tflav.IsGluon() && aflav.IsGluon()) return ConstrainedGGFlows(tbeam);
  if (tflav.IsGluon() && aflav.IsQuark()) return ConstrainedGQFlows(tbeam);
  if (tflav.IsQuark() && aflav.IsGluon()) return ConstrainedQGFlows(tbeam);
  if (tflav.IsQuark() && aflav.IsQuark()) return ConstrainedQQFlows(tbeam);
  msg_Error()<<METHOD<<" no obvious solution available for "
	     <<"["<<tflav<<"]"<<"["<<aflav<<"] and we should not end up here.\n"
	     <<"   Will exit the run.\n";
  exit(1);
}

void Colour_Generator::AssignColours(const size_t & beam,Particle * trip,Particle * anti) {
  // Replace one colour (either triplet or anti-triplet) in either of the two particles,
  // depending on availability of replacement colours on stack.
  // ncol checks if there are colours available on stack:
  // 1 = triplet, 2 = anti-triplet, 3 = both
  size_t ncol = AvailableColours(beam);
  if (ncol==1 || (ncol==3 && ran->Get()>0.5))
    Replace(beam,0,trip);
  else if (ncol==2 || ncol==3)
    Replace(beam,1,anti);
}

size_t Colour_Generator::AvailableColours(const size_t & beam) {
  // Check availability of replacement colours on stack, taking into account vetoed
  // replacements, as encoded in the m_vetoed sets.
  size_t ncolours = 0;
  for (size_t index=0;index<2;index++) {
    for (list<int>::iterator cit=m_cols[beam][index].begin();
	 cit!=m_cols[beam][index].end();cit++) {
      if (m_vetoed[beam][index].find(*cit)==m_vetoed[beam][index].end()) {
	ncolours += (index+1);
	break;
      }
    }
  }
  return ncolours;
}

void Colour_Generator::ReplaceBoth(const int & beam,const size_t & index) {
  // This is the replacement method for t-channel colour exchanges.
  // If a new colour (newcol) can be taken from stack through the NextColour method, replace the old
  // colour (oldcol) in both beam particles with it.  For the recipient beam, this is handled in the
  // ReplaceInIS and ReplaceInFS methods.  Remove the old colour from both conjugate stacks,
  // as it does not need to be balanced any more, and push the new colour on the stack of the recipient
  // beam, where it must be balanced now.  Follow the replacement downstream through the blob list.
  int newcol = NextColour(beam,index);
  if (newcol!=-1) {
    Particle * part = p_inparts[beam];
    int oldcol = part->GetFlow(index+1);
    part->SetFlow(index+1,newcol);
    m_cols[beam][1-index].remove(oldcol);
    m_cols[1-beam][index].remove(oldcol);
    m_cols[1-beam][index].push_back(newcol);
    Blob * showerblob = part->DecayBlob();
    if (showerblob) {
      ReplaceInFS(oldcol,newcol,index,showerblob);
      ReplaceInIS(oldcol,newcol,index,showerblob);
      ReplaceInFS(oldcol,newcol,1-index,showerblob);
      ReplaceInIS(oldcol,newcol,1-index,showerblob);
    }
  }
}

void Colour_Generator::Replace(const int & beam,const size_t & index,ATOOLS::Particle * part) {
  // This is the replacement method for s-channel colour exchanges.
  // If a new colour (newcol) can be taken from stack through the NextColour method, replace the old
  // colour (oldcol) of the respective (beam) particle with it.  Remove oldcol from the conjugate
  // stack, as it does not need to be balanced any more.  Follow the replacement downstream
  // through the blob list.
  int newcol = NextColour(beam,index);
  if (newcol!=-1) {
    int oldcol = part->GetFlow(index+1);
    part->SetFlow(index+1,newcol);
    m_cols[beam][1-index].remove(oldcol);
    Blob * showerblob = part->DecayBlob();
    if (showerblob) {
      ReplaceInFS(oldcol,newcol,index,showerblob);
      ReplaceInIS(oldcol,newcol,index,showerblob);
    }
  }
}

void Colour_Generator::
ReplaceInFS(const int & oldcol, const int & newcol, const size_t & index,Blob * blob) {
  // Replace old colour with new colour at index, in all outgoing particles of the blob.
  // Follow ghte replacement recursively downstream.
  for (size_t j=0;j<blob->NOutP();j++) {
    Particle * part = blob->OutParticle(j);
    if (part->GetFlow(index+1)==oldcol) {
      part->SetFlow(index+1,newcol);
      Blob * dec = part->DecayBlob();
      if (dec!=NULL &&
	  !(dec->Type()==btp::Signal_Process || dec->Type()==btp::Hard_Collision)) {
	ReplaceInFS(oldcol,newcol,index,blob);
      }
    }
    if (blob->Type()==btp::Shower && part->Info()=='I' && part->GetFlow(2-index)==oldcol) {
      part->SetFlow(2-index,newcol);
    }
  }
}

void Colour_Generator::
ReplaceInIS(const int & oldcol, const int & newcol, const size_t & index,Blob * blob) {
  // Replace old colour with new colour at index, in all incoming particles of blob.
  for (size_t j=0;j<blob->NInP();j++) {
    Particle * part = blob->InParticle(j);
    if (part->GetFlow(index+1)==oldcol) part->SetFlow(index+1,newcol);
  }
}

int Colour_Generator::NextColour(const size_t & beam,const size_t & index) {
  // Pull the next colour from the stack specified by beam and (colour) index.
  // Return -1 if there is no allowed colour - disallowed colours are stored in the
  // m_vetoed sets.
  int col = -1;
  for (list<int>::iterator cit=m_cols[beam][index].begin();
       cit!=m_cols[beam][index].end();cit++) {
    if (m_vetoed[beam][index].find(*cit)==m_vetoed[beam][index].end()) {
      col = (*cit);
      m_cols[beam][index].erase(cit);
      break;
    }
  }
  return col;
}

void Colour_Generator::
AddColour(const size_t & beam,const size_t & pos,Particle * const part) {
  // Adds colour to the beam stack - triplets become anti-triplets and vice versa,
  // ready for consumption in the NextColour method.
  if (part->GetFlow(pos+1)!=0) {
    m_cols[beam][1-pos].push_back(part->GetFlow(pos+1));
    m_vetoed[beam][1-pos].insert(part->GetFlow(pos+1));
    m_vetoed[beam][pos].insert(part->GetFlow(pos+1));
  }
}

void Colour_Generator::Output() {
  for (size_t beam=0;beam<2;beam++) {
    for (size_t pos=0;pos<2;pos++) {
      msg_Out()<<"   ["<<beam<<pos<<"]: ";
      for (list<int>::iterator cit=m_cols[beam][pos].begin();
	   cit!=m_cols[beam][pos].end();cit++) msg_Out()<<" "<<(*cit);
      msg_Out()<<"--- vetoed: ";
      for (set<int>::iterator cit=m_vetoed[beam][pos].begin();
	   cit!=m_vetoed[beam][pos].end();cit++) msg_Out()<<" "<<(*cit);
      msg_Out()<<"\n";
    }
  }
}
