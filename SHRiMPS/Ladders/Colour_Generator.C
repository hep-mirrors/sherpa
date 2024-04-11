#include "SHRiMPS/Ladders/Colour_Generator.H"
#include "ATOOLS/Phys/Flow.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"

using namespace SHRIMPS;
using namespace ATOOLS;
using namespace std;

Colour_Generator::Colour_Generator() : m_mode(lcolours::random) {}

Colour_Generator::~Colour_Generator() {}

bool Colour_Generator::operator()(Ladder * ladder){
  p_ladder    = ladder;
  p_emissions = p_ladder->GetEmissions(); 
  p_props     = p_ladder->GetProps();
  PickStartColours();
  IterateColours(p_ladder->GetEmissions()->begin(),p_props->begin());
  PickEndColours();
  return true;
}

void Colour_Generator::PickStartColours() {
  Ladder_Particle * inpart = p_ladder->InPart(0);
  if (inpart->Flavour().IsQuark() && !inpart->Flavour().IsAnti()) {
    inpart->SetFlow(1,-1);
    m_keepindex = 0;
  }
  else if (inpart->Flavour().IsQuark() &&  inpart->Flavour().IsAnti()) {
    inpart->SetFlow(2,-1);
    m_keepindex = 1;
  }
  else if (inpart->Flavour().IsGluon()) {
    for (size_t i=0;i<2;i++) inpart->SetFlow(i+1,-1);
    m_keepindex = m_mode==(lcolours::fix_all ? 1 : (ran->Get()>0.5 ? 0 : 1));
  }
  for (size_t i=0;i<2;i++) m_propcolours[i] = inpart->GetFlow(i+1);
}

void Colour_Generator::
IterateColours(LadderMap::iterator lmit,TPropList::iterator tpit) {
  Ladder_Particle * outpart = &lmit->second;
  T_Prop          * prop    = &*tpit;
  //msg_Out()<<METHOD<<": "<<&lmit->first<<" "
  //	   <<"("<<prop->Col()<<": ["<<m_propcolours[0]<<" "<<m_propcolours[1]<<"])\n";
  if (outpart==&p_ladder->GetEmissions()->rbegin()->second) return;
  if (m_mode==lcolours::no_singlets) prop->SetCol(colour_type::octet);
  if (prop->Col()==colour_type::singlet) {
    for (size_t i=0;i<2;i++) {
      outpart->SetFlow(i+1,m_propcolours[i]);
      m_propcolours[i] = 0;
    }
  }
  else if (prop->Col()==colour_type::octet) {
    if (outpart->Flavour().IsGluon()) {
      if (m_propcolours[0]!=0 && m_propcolours[1]!=0) {
	size_t keep = ((m_mode==lcolours::fix_emissions ||
			m_mode==lcolours::fix_all) ? 
		       m_keepindex : (ran->Get()>0.5 ? 0 : 1 ));
	outpart->SetFlow(keep+1,m_propcolours[keep]);
	outpart->SetFlow(2-keep,-1);
	m_propcolours[keep] = outpart->GetFlow(2-keep);
      }
      else if (m_propcolours[0]==0 && m_propcolours[1]==0) {
	//msg_Out()<<"                           --> have to fill zeros.\n";
	for (size_t i=0;i<2;i++) {
	  outpart->SetFlow(i+1,-1);
	  m_propcolours[1-i] = outpart->GetFlow(i+1);
	}
      }
    }
    else if (outpart->Flavour().IsQuark()) {
      size_t keep = outpart->Flavour().IsAnti() ? 0 : 1;
      if (m_propcolours[keep]!=0 && m_propcolours[1-keep]==0) {
	outpart->SetFlow(2-keep);
	m_propcolours[keep] = outpart->GetFlow(2-keep);
      }
    }
  }
  //msg_Out()<<"                           --> "<<&lmit->first<<" "
  //	   <<"("<<prop->Col()<<": ["<<m_propcolours[0]<<" "<<m_propcolours[1]<<"])\n";
  lmit++; tpit++;
  IterateColours(lmit,tpit);
}

void Colour_Generator::PickEndColours() {
  Ladder_Particle * inpart  = p_ladder->InPart(1);
  Ladder_Particle * outpart = &p_ladder->GetEmissions()->rbegin()->second;
  if (m_propcolours[0]==0 && m_propcolours[1]==0) {
    if ((inpart->Flavour().IsQuark() && !inpart->Flavour().IsAnti()) ||
	inpart->Flavour().IsGluon()) inpart->SetFlow(1,-1);
    if ((inpart->Flavour().IsQuark() &&  inpart->Flavour().IsAnti()) ||
	inpart->Flavour().IsGluon()) inpart->SetFlow(2,-1);
    for (size_t i=1;i<3;i++) outpart->SetFlow(i,inpart->GetFlow(i));
  }
  else if (m_propcolours[0]!=0 && m_propcolours[1]!=0) {
    size_t incol = ( (m_mode==lcolours::fix_emissions ||
		      m_mode==lcolours::fix_all) ?
		     m_keepindex : (ran->Get()>0.5 ? 0 : 1 ));
    inpart->SetFlow(2-incol,m_propcolours[incol]);
    inpart->SetFlow(1+incol,-1);
    outpart->SetFlow(1+incol,inpart->GetFlow(incol+1));
    outpart->SetFlow(2-incol,m_propcolours[1-incol]);
  }
}
