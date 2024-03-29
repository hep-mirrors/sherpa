#include "SHRiMPS/Ladders/Ladder.H"

using namespace SHRIMPS;
using namespace ATOOLS;

Ladder::Ladder(const Vec4D & position,const bool & rescatter) :
m_position(position), m_isrescatter(rescatter), m_type(ladder_type::unknown)
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

void Ladder::DeletePropagator(TPropList::iterator piter) {
  if (piter!=m_tprops.end()) m_tprops.erase(piter); 
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
  winner = m_tprops.end();
  double qt2max = 0., qt2;
  for (TPropList::iterator pit=m_tprops.begin();pit!=m_tprops.end();pit++) {
    qt2 = pit->QT2(); 
    if (qt2>qt2max) {
      winner = pit;
      qt2max = qt2;
    }
  }
  return (winner!=m_tprops.end());
}

void Ladder::FixType(const double & ymin, const double & ymax) {
  //msg_Out()<<METHOD<<"("<<ymin<<", "<<ymax<<") for " <<m_tprops.size()
  //<<" propagators and "<<m_emissions.size()<<" emissions.\n";
  if (m_tprops.size()==1 && m_tprops.begin()->Col()==colour_type::singlet) {
    m_type = ladder_type::elastic;
  }
  else {
    //double maxSDmass2 = 2.5; // 62 GeV - ~20% mass cut
    //double maxSDmass2 = 282; // 7 TeV - 20% mass cut
    double maxSDmass2 = -1.; // no mass cut
    int i = 0;
    //for (LadderMap::iterator emissions = m_emissions.begin();
    //	 emissions != m_emissions.end(); emissions++) {
      //msg_Out()<<"i = "<<i<<" emission's y: "<<emissions->first<<" emission's second: "
    //<<emissions->second<<std::endl;
    //  i++; 
    //}
    size_t forwards=0, backwards=0;
    TPropList::iterator pit1 = m_tprops.begin();
    LadderMap::iterator lmit = m_emissions.begin();
    i = 0;
    while (lmit->first>ymax) {
      //msg_Out()<<"i = "<<i<<" positive beam zone emission's momentum: "
      //<<lmit->second.Momentum()<<" emission's y: "<<lmit->first<<std::endl;
      i++;
      forwards++;
      lmit++; pit1++;
    }
    Vec4D PmiddleAndBackwards;
    while(lmit != m_emissions.end()) {
      PmiddleAndBackwards += lmit->second.Momentum();
      //msg_Out()<<"i = "<<i<<" (BW) middle or negative beam zone emission's momentum: "
      //<<lmit->second.Momentum()<<" emission's y: "<<lmit->first<<std::endl;
      i++;
      lmit++;
    }
    i = m_emissions.size() - 1;
    TPropList::iterator pit2 = m_tprops.end(); pit2--;
    lmit = m_emissions.end(); lmit--;
    while (lmit->first<ymin) {
      //msg_Out()<<"i = "<<i<<" negative beam zone emission's momentum: "
      //<<lmit->second.Momentum()<<" emission's y: "<<lmit->first<<std::endl;
      i--;
      backwards++;
      lmit--; pit2--;
    }
    Vec4D PmiddleAndForwards;
    while(lmit != m_emissions.begin()) {
      PmiddleAndForwards += lmit->second.Momentum();
      //msg_Out()<<"i = "<<i<<" (FW) middle or positive beam zone emission's momentum: "
      //<<lmit->second.Momentum()<<" emission's y: "<<lmit->first<<std::endl;
      i--;
      lmit--;
    }
    PmiddleAndForwards += lmit->second.Momentum();
    //msg_Out()<<"i = "<<i<<" (FW) middle or positive beam zone emission's momentum: "
    //<<lmit->second.Momentum()<<" emission's y: "<<lmit->first<<std::endl;
    double mass2FW = PmiddleAndForwards.Abs2();
    double mass2BW = PmiddleAndBackwards.Abs2();
    double largestMass2 = 0.;
    if (mass2FW > mass2BW) largestMass2 = mass2FW;
    else largestMass2 = mass2BW;
    //std::ofstream massSDbw, massSDfw, massInelBW, massInelFW;
    //massSDbw.open("./massSDbw62G.txt", std::ios::app);
    //massSDfw.open("./massSDfw62G.txt", std::ios::app);
    //massInelBW.open("./massInelBW62G.txt", std::ios::app);
    //massInelFW.open("./massInelFW762G.txt", std::ios::app);
    if ((&*pit1)!=(&*pit2) ||
	pit1->Col()==colour_type::octet || pit2->Col()==colour_type::octet ||
	pit1->Col()==colour_type::triplet || pit2->Col()==colour_type::triplet) {
      if (largestMass2 > maxSDmass2) {
        m_type = ladder_type::inelastic;//THIS SHOULD BE inelastic!!!!!!! IT WAS JUST A TEST!!!!!
        //massInelBW << mass2BW << std::endl;
        //massInelFW << mass2FW << std::endl;
      }
      else if (forwards>1 && backwards>1) m_type = ladder_type::double_diffractive;
      else {
        m_type = ladder_type::single_diffractive;//THIS SHOULD BE single_diffractive!!!!!!! IT WAS JUST A TEST!!!!!
        //massSDbw << mass2BW << std::endl;
        //massSDfw << mass2FW << std::endl;
      }
    }
    else if (forwards>1 && backwards>1) {
      m_type = ladder_type::double_diffractive;
    }
    else {
      m_type = ladder_type::single_diffractive;//THIS SHOULD BE single_diffractive!!!!!!! IT WAS JUST A TEST!!!!!;
      //massSDbw << mass2BW << std::endl;
      //massSDfw << mass2FW << std::endl;
    }
    //massSDfw.close();
    //massInelFW.close();
    //massSDbw.close();
    //massInelBW.close();
  }
  //msg_Out()<<METHOD<<"(type = "<<int(m_type)<<"): " <<m_inpart[0].Momentum()<<" + "<<m_inpart[1].Momentum()<<"\n";
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
  for (LadderMap::const_iterator yiter=m_emissions.begin(); yiter!=m_emissions.end();yiter++) msg_Out()<<yiter->first<<" - ";
  msg_Out()<<"===\n"; 
}

std::ostream & SHRIMPS::
operator<<(std::ostream & s,Ladder & ladder) {
s<<" ---------------------------------------------------------\n"
<<"Ladder ("<<ladder.GetProps()->size()<<" props) "
<<"at position "<<ladder.Position()<<" (b= "
<<(sqrt(sqr(ladder.Position()[1])+sqr(ladder.Position()[2])))<<"):\n"
<<" in = "<<(*ladder.InPart(0))<<"\n"
<<" "<<(*ladder.InPart(1))<<"\n";
int i(0);
TPropList::const_iterator citer=ladder.GetProps()->begin();
for (LadderMap::const_iterator yiter=ladder.GetEmissions()->begin();
yiter!=ladder.GetEmissions()->end();yiter++) {
s<<" y_{"<<i<<"} = "<<yiter->first<<", k_{"<<i<<"} = "<<yiter->second;
if (citer!=ladder.GetProps()->end()) {
s<<" "<<(*citer);
citer++;
}
i++;
}
s<<" ---------------------------------------------------------\n";
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

