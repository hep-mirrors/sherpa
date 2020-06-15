#include "CFPSHOWER++/Tools/Parton.H"
#include "CFPSHOWER++/Tools/Splitting.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"

using namespace CFPSHOWER;
using namespace ATOOLS;
using namespace std;

namespace CFPSHOWER {
  size_t Parton::s_cid=0;
  size_t Parton::s_cnt=0;
}


Parton::Parton(const Flavour & flav,const Vec4D & mom,
	       const Color & color,const size_t & beam,const size_t & id) :
  m_flav(flav), m_mom(mom), m_color(color), m_beam(beam), m_xB(0.), m_on(true)
{
  s_cid++; s_cnt++;
  if (id>0) m_id = id;
  else m_id = s_cid;
  m_offsprings.resize(2);
  for (size_t i=0;i<3;i++) p_softpartners[i] = NULL;
}

Parton::~Parton() {
  s_cnt--;
}

void Parton::SetXB() 
{
  if (m_beam==1) m_xB = -m_mom.PPlus()/rpa->gen.PBeam(0).PPlus();
  if (m_beam==2) m_xB = -m_mom.PMinus()/rpa->gen.PBeam(1).PMinus();
}

void Parton::AddWeight(const Splitting & split,const bool & accept) {
  // add acceptance or rejection weight for specific spectator to its weight vector.
  // The weight vector consists of BranchingWeights, linking the branching scale t to
  // the weight (probability) associated with the splitting (either reject or accept).
  // The weight is cumulative, from the last splitting, i.e. since the last time
  // ClearWeights() has been called.
  // TODO: will have to upgrade this later to capture the effect of on-the-flight
  // variations.
  double weight = (accept?split.GetWeight()->Accept():split.GetWeight()->Reject());
  if (weight!=1.0) {
    BranchingWeight wt;
    Parton * spec         = split.GetSpectator();
    PWV_Map::iterator wit = m_specweights.insert(make_pair(spec,Weight_Vector())).first;
    BranchingWeight brwt  = (wit->second.empty() ?
			     BranchingWeight(0.,1.) :
			     wit->second.back());
    brwt.m_t              = split.t(0);
    brwt.m_weight        *= weight;
    wit->second.push_back(brwt);
  }
}

double Parton::GetWeight(const double & t) const {
  if (m_specweights.empty()) return 1.;
  double weight = 1.;
  // Iterate over all weights accumulated so far (with different spectators)
  // bisection in branching scales to match the scales stored with the argument t,
  // then multiply the weight with the acceptance/rejection weights.
  for (PWV_Map::const_iterator wit = m_specweights.begin();
       wit!=m_specweights.end(); wit++) {
    const Weight_Vector & wts = wit->second;
    size_t upper = 0, lower = wts.size()-1, mid = (upper+lower)/2, take;
    double scale = wts[mid].m_t;
    while (lower-upper>1) {
      ((t>scale) ? lower : upper) = mid;
      mid = (upper+lower)/2;
      scale = wts[mid].m_t;
    }
    if      (t <= wts[lower].m_t) weight *= wts[lower].m_weight;
    else if (t <= wts[upper].m_t) weight *= wts[upper].m_weight;
  }
  return weight;
}

ostream & CFPSHOWER::operator<<(ostream & s,Parton & part) {
  if (part.On()) s<<"  "; else s<<"# ";
  s<<"Parton("<<part.Id()<<", beam = "<<part.Beam()<<"): "
   <<"["<<part.Flav()<<", mom = "<<part.Mom()<<", ";
  if (part.Beam()>0) s<<"xB = "<<part.XB()<<", ";
  s<<"col = "<<part.GetColor()<<"] <--> ";
  for (Parton_List::const_iterator pit=part.GetSpectators()->begin();
       pit!=part.GetSpectators()->end(); pit++) s<<(*pit)->Id()<<" ";
  s<<"\n";
  return s;
}
