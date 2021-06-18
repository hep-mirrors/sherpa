#include "SHRiMPS/Ladders/Ladder.H"

using namespace SHRIMPS;
using namespace ATOOLS;


Ladder::Ladder(const Vec4D & position) :
  m_position(position)
{
  for (size_t i=0;i<2;i++) {
    m_inpart[i] = Ladder_Particle(Flavour(kf_gluon),Vec4D(0.,0.,0.,0.),m_position);
    m_inpart[i].SetIS(true);
  }
}

Ladder::~Ladder() {
  m_emissions.clear();
  m_tprops.clear();
}

Ladder_Particle * Ladder::
AddRapidity(const double y,const Flavour & flav,const Vec4D & mom) {
  return &(m_emissions[y] = Ladder_Particle(flav,mom,m_position));
}

void Ladder::DeleteRapidity(LadderMap::iterator yiter) {
  if (yiter!=m_emissions.end()) m_emissions.erase(yiter);
}

void Ladder::AddPropagator(T_Prop prop)
{
  m_tprops.push_back(prop);
}

void Ladder::UpdatePropagatorKinematics() {
  LadderMap::iterator lmit = m_emissions.begin();
  TPropList::iterator tpit = m_tprops.begin();
  Vec4D q = m_inpart[0].Momentum();
  while (tpit!=m_tprops.end()) {
    q -= lmit->second.Momentum();
    tpit->SetQ2(dabs(q.Abs2()));
    tpit->SetQ(q);
    tpit++;
    lmit++;
  }
  LadderMap::reverse_iterator lmrit = m_emissions.rbegin();
  TPropList::reverse_iterator tprit = m_tprops.rbegin();
  q = m_inpart[1].Momentum();
  while (tprit!=m_tprops.rend()) {
    if (tprit->Q()[0]<0.) {
      q -= lmrit->second.Momentum();
      tprit->SetQ2(dabs(q.Abs2()));
      tprit->SetQ(q);
      tprit++;
      lmrit++;
    }
    else break;
  }
}

bool Ladder::ExtractHardest(TPropList::iterator & winner,
			    const double & qt2min) {
  winner        = m_tprops.end();
  double qt2max = 0., qt2;
  TPropList::iterator pit;
  for (pit=m_tprops.begin();pit!=m_tprops.end();pit++) {
    qt2      = pit->QT2();	 
    if (qt2>qt2max) {
      winner = pit;
      qt2max = qt2;
    }
  }
  return (winner!=m_tprops.end());
}

void Ladder::HardestIncomingMomenta(const TPropList::iterator & winner,
				    Vec4D & q0,Vec4D & q1) {
  q0 = m_inpart[0].Momentum(); q1 = m_inpart[1].Momentum();
  TPropList::iterator pit = winner;
  pit++;
  if (pit!=m_tprops.end()) q1 = pit->Q();
  if (winner==m_tprops.begin()) return;
  pit = winner;
  pit--;
  q0 = pit->Q();
}

void Ladder::Reset(const bool & all) { 
  m_emissions.clear(); 
  m_tprops.clear();
}

void Ladder::ResetFS() {
  m_emissions.clear(); 
  m_tprops.clear();
}

void Ladder::OutputRapidities() {
  msg_Out()<<"=== - ";
  for (LadderMap::const_iterator yiter=m_emissions.begin();
       yiter!=m_emissions.end();yiter++) msg_Out()<<yiter->first<<" - ";
  msg_Out()<<"===\n";  
}

std::ostream & SHRIMPS::
operator<<(std::ostream & s,Ladder & ladder) {
  s<<"   ---------------------------------------------------------\n"
   <<"Ladder ("<<ladder.GetProps()->size()<<" props) "
   <<"at position "<<ladder.Position()<<" (b= "
   <<(sqrt(sqr(ladder.Position()[1])+sqr(ladder.Position()[2])))<<"):\n"
   <<"  in = "<<(*ladder.InPart(0))<<"\n"
   <<"       "<<(*ladder.InPart(1))<<"\n";
  int i(0);
  TPropList::const_iterator citer=ladder.GetProps()->begin();
  for (LadderMap::const_iterator yiter=ladder.GetEmissions()->begin();
       yiter!=ladder.GetEmissions()->end();yiter++) {
    s<<"  y_{"<<i<<"} = "<<yiter->first<<", k_{"<<i<<"} = "
     <<yiter->second<<"\n";
    if (citer!=ladder.GetProps()->end()) {
      s<<"   "<<(*citer)<<"\n";
      citer++;
    }
    i++;
  }
  s<<"   ---------------------------------------------------------\n";
  return s;
}

bool Ladder::CheckFourMomentum() {
  Vec4D check(m_inpart[0].Momentum()+m_inpart[1].Momentum());
  double shat(check.Abs2());
  for (LadderMap::iterator liter=m_emissions.begin();
       liter!=m_emissions.end();liter++) {
    check -= liter->second.Momentum();
  }
  if (dabs(check.Abs2())/shat>1.e-6) {
    msg_Error()<<"-------------------------------------------\n"
	       <<METHOD<<" failed: check = "<<check<<", "<<check.Abs2()<<"\n"
	       <<(*this)<<"\n";
    return false;
  }
  return true;
}
