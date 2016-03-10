#include "SHRiMPS/Beam_Remnants/Colour_Generator.H"
#include "ATOOLS/Org/Message.H"

using namespace SHRIMPS;
using namespace ATOOLS;
using namespace std;

Colour_Generator::Colour_Generator() {}

Colour_Generator::~Colour_Generator() {}

bool Colour_Generator::operator()(Ladder * ladder){
  p_ladder    = ladder;
  p_emissions = p_ladder->GetEmissions(); 
  p_props     = p_ladder->GetProps();
  if (!p_ladder->IsRescatter()) PickStartColours();
  IterateColours(p_emissions->begin(),p_props->begin(),
		 p_ladder->GetIn(0)->GetFlow(1),
		 p_ladder->GetIn(0)->GetFlow(2));
  UpdateColours();
  return true;
}

void Colour_Generator::PickStartColours() {
  Ladder_Particle * beampart(p_ladder->GetIn(0));
  int beam(beampart->m_mom[3]>0.?1:0);
  int newcols[2];
  bool hit(false);
  for (size_t pos=0;pos<2;pos++) {
    newcols[pos] = 0;
    // make sure not to fill the wrong colours for quarks
    if (beampart->m_flav.IsQuark() &&
	((!beampart->m_flav.IsAnti() && pos==0) ||
	 (beampart->m_flav.IsAnti() && pos==1))) continue;
    // take first colour from stack -- but do not take a second one
    // the stuff below is tricky, since our Flow class has its colours
    // on components/counters 1 and 2 ...
    if (!hit && !m_colours[beam][pos].empty()) {
      beampart->SetFlow(2-pos,(*m_colours[beam][pos].begin()));
      m_colours[beam][pos].erase(m_colours[beam][pos].begin());
      hit = true;
    }
    // no colour on stack - generate new one
    else {
      beampart->SetFlow(2-pos);
      newcols[pos] = beampart->GetFlow(2-pos);
    }
  }
  // put new colours on stack
  for (size_t pos=0;pos<2;pos++) {
    if (newcols[pos]!=0) m_colours[0][1-pos].insert(newcols[pos]);
  }
}

void Colour_Generator::IterateColours(LadderMap::iterator out,
				      TPropList::iterator prop,
				      int col1,int col2) {
  Ladder_Particle * lpart = &out->second;
  if (prop->m_col==colour_type::singlet) {
    lpart->SetFlow(1,col1);
    lpart->SetFlow(2,col2);
    col1 = col2 = -1;
  }
  else {
    if (col1==-1) {
      lpart->SetFlow(1);
      col2 = lpart->GetFlow(1);
    }
    else lpart->SetFlow(1,col1);
    lpart->SetFlow(2);
    col1 = lpart->GetFlow(2);
  }
  out++; prop++;
  if (prop!=p_props->end()) IterateColours(out,prop,col1,col2);
  else FinishColours(col1,col2);
}

void Colour_Generator::FinishColours(int col1,int col2) {
  Ladder_Particle * lpart(&p_ladder->GetEmissions()->rbegin()->second);
  Ladder_Particle * beampart(p_ladder->GetIn(1));
  lpart->SetFlow(1,col1);
  lpart->SetFlow(2);
  bool issinglet(p_ladder->GetProps()->rbegin()->m_col==colour_type::singlet);
  beampart->SetFlow(1,issinglet?lpart->GetFlow(1):col2);
  beampart->SetFlow(2,lpart->GetFlow(2));
}

void Colour_Generator::UpdateColours() {
  if (!(m_colours[1][0].size()==0 && m_colours[1][1].size()==0)) {
    int trial((m_colours[1][0].size()<m_colours[1][1].size())?1:2);
    bool hit(false);
    Ladder_Particle * beam0(p_ladder->GetIn(0)), * beam1(p_ladder->GetIn(1));
    if (beam0->GetFlow(3-trial)!=beam1->GetFlow(trial))
      hit = ReplaceColours(trial);
    if (!hit && beam0->GetFlow(trial)!= beam1->GetFlow(3-trial))
      hit = ReplaceColours(3-trial);
  }
  else {
    for (size_t pos=1;pos<3;pos++)
      m_colours[1][2-pos].insert(p_ladder->GetIn(1)->GetFlow(pos));
  }
}

bool Colour_Generator::ReplaceColours(const size_t & pos) {
  for (set<int>::iterator cit=m_colours[1][pos-1].begin();
       cit!=m_colours[1][pos-1].end();cit++) {
    int colnew((*cit));
    if (colnew!=p_ladder->GetIn(1)->GetFlow(3-pos)) {
      m_colours[1][pos-1].erase(cit);
      int colold(p_ladder->GetIn(1)->GetFlow(pos));
      p_ladder->GetIn(1)->SetFlow(pos,colnew);
      for (LadderMap::iterator liter=p_ladder->GetEmissions()->begin();
	   liter!=p_ladder->GetEmissions()->end();liter++) {
	if (liter->second.GetFlow(pos)==colold) {
	  liter->second.SetFlow(pos,colnew);
	  return true;
	}
      }
    }
  }
  m_colours[1][2-pos].insert(p_ladder->GetIn(1)->GetFlow(pos));
  return false;
}

void Colour_Generator::Reset() { 
  for (size_t beam=0;beam<2;beam++) { 
    for (size_t flow=0;flow<2;flow++) { 
      m_colours[beam][flow].clear();
    } 
  } 
}

void Colour_Generator::OutputStack() {
  for (int beam=0;beam<2;beam++) {
    for (int pos=0;pos<2;pos++) {
      msg_Out()<<"Colours in stack["<<beam<<"]["<<pos<<"] : {";
      for (set<int>::iterator col=m_colours[beam][pos].begin();
	   col!=m_colours[beam][pos].end();col++) {
	msg_Out()<<" "<<(*col);
      }
      msg_Out()<<" }\n";
    }
  }
}
