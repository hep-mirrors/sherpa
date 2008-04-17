#include "All_Decays.H"
#include "Model_Base.H"
#include "Width_Calculator.H"
#include "Vertex.H"
#include <complex>

using namespace MODEL;
using namespace ATOOLS;


All_Decays::All_Decays(Model_Base * model) :
  p_vertextable(model->GetVertexTable()),p_decays(new DecayMap)
{  
  //std::cout<<METHOD<<" vertextable = "<<p_vertextable<<" with "<<p_vertextable->size()
  //	   <<" from model "<<model<<std::endl;
}

All_Decays::~All_Decays()
{
  delete p_decays;
}

bool All_Decays::AddToDecays(const ATOOLS::Flavour & flav) 
{
  if (CheckInVertex(flav) && p_decays->find(flav)==p_decays->end()) {
    (*p_decays)[flav] = new Decay_Table(flav);
    return true;
  }
  return false;
}

bool All_Decays::AddToDecays(ATOOLS::Decay_Channel * dc) 
{
  Flavour flav = dc->GetDecaying();
  // This is maybe a bit too much - maybe if there's no vertex, 
  // we should have isotropic decays.
  if (!CheckInVertex(flav)) return false;
  DMIterator dmit = p_decays->find(flav);
  if (dmit==p_decays->end()) {
    Decay_Table * dt = new Decay_Table(flav);
    dt->SetWidthGenerator(std::string("Sherpa"));
    dt->AddDecayChannel(dc);
    (*p_decays)[flav] = new Decay_Table(flav);
  }
  else dmit->second->AddDecayChannel(dc);
  if (dc->ProcessName()==std::string("")) dc->SetProcessName("");

  return true;
}

bool All_Decays::InitializeDecayTables() {
  if (p_decays->empty()) return false;
  BinaryDecays();
  ThreeBodyDecays();
  ArrangeDecays();
  PrintDecayTables();
  return true;
}

bool All_Decays::CalculateWidths() {
  Width_Calculator calc(p_vertextable);
  for (DMIterator dmit=p_decays->begin();dmit!=p_decays->end();dmit++) {
    for (int i=0;i<dmit->second->NumberOfDecayChannels();i++) {
      dmit->second->GetDecayChannel(i)->SetWidth(calc.Width(dmit->second->GetDecayChannel(i)));
    }
    dmit->second->UpdateWidth();
 }
  PrintDecayTables();
  return true;
}

bool All_Decays::CheckInVertex(const ATOOLS::Flavour & flav) {
  Vertex_List & vertexlist = (*p_vertextable)[flav];
  if (vertexlist.size()>0) return 1;
  return 0;
}

void All_Decays::BinaryDecays()
{
  Decay_Channel * dc(NULL);
  Decay_Table * dt(NULL);
  Single_Vertex * sv;
  Flavour inflav;
  bool    construct;
  double  maxest;
  for (DMIterator dmit=p_decays->begin();dmit!=p_decays->end();dmit++) {
    inflav = dmit->first;
    dt     = dmit->second;
    maxest = 0.;
    construct = (dt->NumberOfDecayChannels()==0);
    if (construct) {
      Vertex_List & vertexlist = (*p_vertextable)[inflav];
      for (size_t i=0;i<vertexlist.size();i++) {
	sv = vertexlist[i];
	if (sv->on && sv->nleg==3 &&
	    sv->in[0]!=Flavour(kf_shgluon) && 
	    sv->in[1]!=Flavour(kf_shgluon) && 
	    sv->in[2]!=Flavour(kf_shgluon)) {
	  if (inflav.Mass()>sv->in[1].Mass()+sv->in[2].Mass()) {
	    dc = new Decay_Channel(inflav);
	    dc->AddDecayProduct(sv->in[1]);
	    dc->AddDecayProduct(sv->in[2]);
	    dc->SetProcessName("");
	    dt->AddDecayChannel(dc);
	    if (std::norm(sv->cpl[0].Value())+std::norm(sv->cpl[1].Value())>maxest) 
	      maxest = std::norm(sv->cpl[0].Value())+std::norm(sv->cpl[1].Value());
	  }
	}
      }
    }
    else {
      for (short int i=0;i<dt->NumberOfDecayChannels();i++) {
	dc = dt->GetDecayChannel(i);
	if (dc->NumberOfDecayProducts()==2) {
	  // must check whether vertex exists, otherwise, isotropic decay to be 
	  // initialised - needs some infrastructure with metype/pstype.
	  // to be done.
	}
      }
    }
    m_maxestimates[inflav] = maxest;
  }
}

void All_Decays::ThreeBodyDecays() {
  Decay_Channel * dc(NULL);
  Decay_Table * dt(NULL);
  Single_Vertex * sv1,* sv2;
  Flavour inflav;
  double  crit,est;
  for (DMIterator dmit=p_decays->begin();dmit!=p_decays->end();dmit++) {
    inflav    = dmit->first;
    crit      = m_maxestimates.find(inflav)->second;
    dt        = dmit->second;
    Vertex_List & vertexlist = (*p_vertextable)[inflav];
    for (size_t i=0;i<vertexlist.size();i++) {
      sv1 = vertexlist[i];
      if (sv1->on && sv1->nleg==3 &&
	  sv1->in[0]!=Flavour(kf_shgluon) && 
	  sv1->in[1]!=Flavour(kf_shgluon) && sv1->in[1]!=inflav && 
	  sv1->in[2]!=Flavour(kf_shgluon) && sv1->in[2]!=inflav) {
	est = std::norm(sv1->cpl[0].Value())+std::norm(sv1->cpl[1].Value());
	if (inflav.Mass()<sv1->in[1].Mass()+sv1->in[2].Mass() &&
	    (sv1->in[1].Mass()+sv1->in[2].Mass()-inflav.Mass())/
	    sqrt(sqr(sv1->in[1].Width())+sqr(sv1->in[2].Width()))<est/crit) {
	  for (int j=1;j<3;j++) {
	    Vertex_List & vertexlist2 = (*p_vertextable)[sv1->in[j]];
	    for (size_t k=0;k<vertexlist2.size();k++) {
	      sv2 = vertexlist2[k];
	      if (sv2->on && sv2->nleg==3) {
		if (inflav.Mass()>sv1->in[3-j].Mass()+sv2->in[1].Mass()+sv2->in[2].Mass()) {
		  //std::cout<<sv1->in[0]<<" --> "<<sv1->in[1]<<" "<<sv1->in[2]<<" ("<<j<<") ==> "
		  //	   <<sv2->in[0]<<" --> "<<sv2->in[1]<<" "<<sv2->in[2]<<" + "
		  //	   <<sv1->in[3-j]<<" from "<<vertexlist2.size()<<std::endl;
		  dc = new Decay_Channel(inflav);
		  dc->AddDecayProduct(sv1->in[3-j]);
		  dc->AddDecayProduct(sv2->in[1]);
		  dc->AddDecayProduct(sv2->in[2]);
		  dc->SetProcessName("");
		  dt->AddDecayChannel(dc);
		}
	      }
	    }
	  }
	}
      }
    }
  }  
}

void All_Decays::ArrangeDecays() {
  Flavour in1;
  Decay_Table * ref, * comp;
  std::string dcname;
  for (DMIterator refit=p_decays->begin();refit!=p_decays->end();refit++) {
    in1 = refit->first;
    ref = refit->second;

    DMIterator compit(refit);
    compit++;
    while (compit!=p_decays->end()) {
      if (in1!=compit->first) compit++;
      else {
	msg_Error()<<"Potential error in "<<METHOD<<" : "<<std::endl
		   <<"   found two decay tables with identical initial particle."<<std::endl
		   <<"   Attempt rescue ... "<<std::endl;
	comp = compit->second;
	if (comp==NULL) {
	  msg_Error()<<"   2nd decay table was empty.  Continue & hope for the best."<<std::endl;
	}
	else {
	  msg_Error()<<"   Move decay channels of 2nd table into first decay table."<<std::endl
		     <<"   This should solve the problem."<<std::endl;
	  for (int i=0;i<comp->NumberOfDecayChannels();i++) 
	    ref->AddDecayChannel(comp->GetDecayChannel(i));
	  delete comp;
	}
	DMIterator delit(compit);
	compit++;
	p_decays->erase(delit);
      }
    }

    for (int i=0;i<ref->NumberOfDecayChannels()-1;i++) {
      dcname = ref->GetDecayChannel(i)->ProcessName();
      bool kick(false);
      do {
	for (int j=i+1;j<ref->NumberOfDecayChannels();j++) {
	  if (dcname==ref->GetDecayChannel(j)->ProcessName()) ref->EraseDecayChannel(j);
	}
      } while (kick);
    }
  }
}

void All_Decays::PrintDecayTables() {
  for (DMIterator dmit=p_decays->begin();dmit!=p_decays->end();dmit++) dmit->second->Output();
}
