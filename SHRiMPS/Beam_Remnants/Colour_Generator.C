#include "SHRiMPS/Beam_Remnants/Colour_Generator.H"
#include "ATOOLS/Phys/Flow.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"

using namespace SHRIMPS;
using namespace ATOOLS;
using namespace std;

Colour_Generator::Colour_Generator() : m_dir(2*int(rpa->gen.PBeam(0)[3]>0)-1)
{}

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

/*
size_t Colour_Generator::CountSingletProps() {
  size_t nsinglet = 0;
  for (TPropList::iterator pit=p_props->begin();pit!=p_props->end();pit++) {
    if ((*pit)->Col()==colour_type::singlet) nsinglet++;
  }
  return nsinglet;
}
*/

void Colour_Generator::PickStartColours() {
  //msg_Out()<<METHOD<<" =========================================\n";
  //OutputStack();
  Ladder_Particle * inpart = p_ladder->InPart(0);
  if (!m_colours[0][0].empty()) {
    inpart->SetFlow(1,*m_colours[0][0].begin());
    m_colours[0][0].erase(inpart->GetFlow(1));
  }
  else {
    inpart->SetFlow(1,-1);
    m_colours[0][1].insert(inpart->GetFlow(1));
  }
  if (!m_colours[1][0].empty()) {
    inpart->SetFlow(2,*m_colours[1][0].begin());
  }
  else {
    inpart->SetFlow(2,-1);
  }
  m_colours[0][0].insert(inpart->GetFlow(2));
  m_propcolours[0] = inpart->GetFlow(1);
  m_propcolours[1] = inpart->GetFlow(2);
  //msg_Out()<<"---------------------------------------------\n"
  //	   <<METHOD<<" ("<<m_propcolours[0]<<", "<<m_propcolours[1]<<") for:\n"
  //<<(*inpart)<<"\n";
}

void Colour_Generator::IterateColours(LadderMap::iterator out,TPropList::iterator prop) {
  Ladder_Particle * outpart = &out->second;
  if (m_propcolours[0]==0 && m_propcolours[1]==0) {
    //msg_Out()<<"---------------------------------------------\n"
    //	     <<"---"<<METHOD<<" from singlet.  Pick new colours.\n";
    for (size_t i=0;i<2;i++) {
      outpart->SetFlow(i+1,-1);
      m_propcolours[1-i] = outpart->GetFlow(i+1);
    }
  }
  else {
    outpart->SetFlow(1,m_propcolours[0]);
    if (prop->Col()==colour_type::octet) {
      outpart->SetFlow(2,-1);
      m_propcolours[0] = outpart->GetFlow(2);
    }
    else if (prop->Col()==colour_type::singlet) {
      outpart->SetFlow(2,m_propcolours[1]);
      m_propcolours[0] = m_propcolours[1] = 0;
      if (!m_colours[1][0].empty() &&
	  outpart->GetFlow(2)==(*m_colours[1][0].begin())) {
	//msg_Out()<<"Gotcha!  Next propagator is singlet - need to update colours\n"
	//	 <<"   Get "<<outpart->GetFlow(2)<<" out of stack and initial state\n"
	//	 <<"   "<<(*p_ladder->InPart(0))<<"\n";
	//OutputStack();
	m_colours[0][0].erase(outpart->GetFlow(2));
	outpart->SetFlow(2,-1);
	p_ladder->InPart(0)->SetFlow(2,outpart->GetFlow(2));
	//msg_Out()<<"   and replace with "<<outpart->GetFlow(2)<<", also in initial state:\n"
	//	 <<"   "<<(*p_ladder->InPart(0))<<"\n"
	//	 <<"   "<<(*outpart)<<"\n";
	m_colours[0][0].insert(outpart->GetFlow(2));
	//OutputStack();
      }
    }
  }
  //msg_Out()<<"---------------------------------------------\n"
  //	   <<METHOD<<": "<<(prop!=p_props->end()?prop->Col():colour_type::none)
  //	   <<" ("<<m_propcolours[0]<<", "<<m_propcolours[1]<<") for:\n"<<(*outpart)<<"\n";
  out++;
  prop++;
  if (prop!=p_props->end()) IterateColours(out,prop);
}

void Colour_Generator::PickEndColours() {
  Ladder_Particle * outpart = &p_ladder->GetEmissions()->rbegin()->second;
  Ladder_Particle * inpart  = p_ladder->InPart(1);
  if (m_propcolours[0]!=0 && m_propcolours[0]!=0) {
    outpart->SetFlow(1,m_propcolours[0]);
    inpart->SetFlow(1,m_propcolours[1]);
    if (m_colours[1][0].find(inpart->GetFlow(1))!=m_colours[1][0].end()) {
      //msg_Out()<<"   erase old "<<inpart->GetFlow(1)<<" from stack.\n";
      m_colours[1][0].erase(inpart->GetFlow(1));
    }
    else if (!m_colours[1][0].empty() &&
	     p_ladder->InPart(0)->GetFlow(2)!=m_propcolours[1]) {
      //msg_Out()<<METHOD<<" with a singlet in between!\n"<<(*p_ladder)<<"\n";
      //OutputStack();
      inpart->SetFlow(1,*m_colours[1][0].begin());
      //msg_Out()<<"   erase "<<inpart->GetFlow(1)<<" from stack.\n";
      m_colours[1][0].erase(inpart->GetFlow(1));
      //msg_Out()<<"   need to replace final state colours "
      //	       <<m_propcolours[1]<<" --> "<<inpart->GetFlow(1)<<"\n";
      ReplaceFSColour(1,m_propcolours[1],inpart->GetFlow(1));
    }
    else {
      //msg_Out()<<"   insert new "<<inpart->GetFlow(1)<<" into stack.\n";
      m_colours[1][1].insert(inpart->GetFlow(1));
    }
    inpart->SetFlow(2,-1);
    outpart->SetFlow(2,inpart->GetFlow(2));
    m_colours[1][0].insert(inpart->GetFlow(2));
  }
  else {
    //msg_Out()<<METHOD<<" with a singlet as last propagator\n";
    if (!m_colours[1][0].empty()) {
      inpart->SetFlow(1,*m_colours[1][0].begin());
      m_colours[1][0].erase(inpart->GetFlow(1));
    }
    else {
      inpart->SetFlow(1,-1);
      m_colours[1][1].insert(inpart->GetFlow(1));
    }
    inpart->SetFlow(2,-1);
    m_colours[1][0].insert(inpart->GetFlow(2));
    for (size_t i=1;i<3;i++) outpart->SetFlow(i,inpart->GetFlow(i));
  }
  //msg_Out()<<"---------------------------------------------\n"
  //	   <<METHOD<<" for (col = "
  //	   <<((m_colours[1][0].size()>0)?(*m_colours[1][0].begin()):0)<<"):\n"
  //	   <<(*inpart)<<"\n"<<(*outpart)<<"\n";
  //OutputStack();
}

bool Colour_Generator::
ReplaceFSColour(const size_t & pos,const int & orig,const int & repl) {
  for (LadderMap::iterator lmit=p_emissions->begin();lmit!=p_emissions->end();lmit++) {
    if (lmit->second.GetFlow(pos)==orig) {
      lmit->second.SetFlow(pos,repl);
      //msg_Out()<<METHOD<<" replaced colour["<<pos<<"] = "<<orig<<" with "<<repl<<"\n"
      //       <<"   "<<lmit->second<<"\n"
      //       <<(*p_ladder)<<"\n";
      return true;
    }
  }
  return false;
}

void Colour_Generator::Reset() { 
  for (size_t beam=0;beam<2;beam++) { 
    for (size_t flow=0;flow<2;flow++) { 
      m_colours[beam][flow].clear();
    }
    m_propcolours[beam] = 0;
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
