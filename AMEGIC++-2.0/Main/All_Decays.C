#include "All_Decays.H"

using namespace AMEGIC;
using namespace ATOOLS;
using namespace std;

All_Decays::All_Decays(Interaction_Model_Base * _model,Topology * _top) :
  p_model(_model), p_top(_top), p_decay(NULL)
{
  Vertex * vertex = p_model->GetVertex();
  for (int i=0;i<vertex->MaxNumber();++i) {
    if ((*vertex)[i]->on && ((*vertex)[i]->nleg==3)) {
      m_vertextable[(*vertex)[i]->in[0]].push_back((*vertex)[i]);
    }
  }
}

bool All_Decays::AddToDecays(const Flavour & flav) 
{
  DMIterator dit = m_decays.find(flav);
  if (dit->first==flav) return 0;
  if (CheckInVertex(flav)) {
    m_particles.insert(flav);
    return 1;
  }
  msg.Error()<<"Error in All_Decays::AddToDecays("<<flav<<") :"<<endl
	     <<"   could not add flavour to list of all_decays, no vertex found. Abort run."<<endl;
  abort();
  return 0;
}

bool All_Decays::AddToDecays(ATOOLS::Decay_Channel * _dec) 
{
  Flavour flav = _dec->GetDecaying();
  DMIterator dit = m_decays.find(flav);
  if (dit->first==flav) {
    msg.Error()<<"Error in All_Decays::AddToDecays("<<flav<<") :"<<endl
	       <<"   could not add flavour to list of all_decays with specfied decay."<<endl
	       <<"   Already booked for unspecified decays. Will continue."<<endl;
    return 0;
  }
  if (CheckInVertex(flav)) {
    double mass = 0.;
    Full_Decay_Channel * dc = new Full_Decay_Channel(_dec);
    if (dc->CreateDecay()) {
      Full_Decay_Table * dt = new Full_Decay_Table(flav,false);
      dt->AddDecayChannel(dc);
      msg.Tracking()<<"Added Decay_Channel :"<<endl;
      if (msg.LevelIsTracking()) dc->Output();
      m_decays.insert(std::make_pair(flav,dt));
      return 1;
    }
    msg.Error()<<"Error in All_Decays::AddToDecays("<<flav<<") :"<<endl
	       <<"   Specified decay is impossible : ";dc->Output();
    msg.Error()<<"   Will continue and hope for the best."<<endl;
    delete dc;
    return 0;
  }
  msg.Error()<<"Error in All_Decays::AddToDecays("<<flav<<") :"<<endl
	     <<"   could not add flavour to list of all_decays, no vertex found. Abort run."<<endl;
  abort();
  return 0;
}

void All_Decays::PrintDecayings()
{
  if (msg.LevelIsTracking()) {
    msg.Out()<<"Final list : "<<m_particles.size()<<endl;
    for (FlSetIter flit=m_particles.begin();flit!=m_particles.end();++flit) msg.Out()<<(*flit)<<endl;
  }
}

bool All_Decays::CheckInVertex(Flavour flav)
{
  Vertex_List & vertexlist = m_vertextable[flav];
  if (vertexlist.size()>0) return 1;
  return 0;
}

bool All_Decays::CalculateWidths(std::string)
{
  for (DMIterator dmit=m_decays.begin();dmit!=m_decays.end();dmit++) {
    dmit->second->CalculateWidths();
  }
  return 1;
}

double All_Decays::Width(ATOOLS::Flavour _fl) 
{
  DMIterator dit = m_decays.find(_fl);
  if (dit->first==_fl) return dit->second->Width();
  msg.Error()<<"Error in All_Decays::Width("<<_fl<<")."<<endl
	     <<"   Flavour not found in internal tables. Return width = 0 GeV."<<endl;
  return 0.;
}

Full_Decay_Table * All_Decays::GetFullDecayTable(ATOOLS::Flavour _fl) 
{
  DMIterator dit = m_decays.find(_fl);
  if (dit->first==_fl) return dit->second;
  msg.Error()<<"Error in All_Decays::GetFullDecayTable("<<_fl<<")."<<endl
	     <<"   Flavour not found in internal tables. Return NULL."<<endl;
  return NULL;
}

bool All_Decays::UnweightedEvent(ATOOLS::Decay_Channel * _dec,double _mass)
{
  DMIterator dit = m_decays.find(_dec->GetDecaying());
  for (int i=0;i<dit->second->NumberOfChannels();i++) {
    if (dit->second->GetChannel(i)==_dec) {
      p_decay = dit->second->GetFullChannel(i);
      return p_decay->OneEvent(_mass);
    } 
  }
}

bool All_Decays::InitializeDecayTables() {
  BinaryDecays();
  //ThreeBodyDecays();
  ArrangeDecays();
  return InitializeDecays();
}

void All_Decays::BinaryDecays()
{
  Vertex_List          vertexlist;
  Flavour              flavs[3];
  Full_Decay_Table   * dt;
  Full_Decay_Channel * dc;
  bool                 skippit;
  DMIterator           dmit;
  for (FlSetIter flit=m_particles.begin();flit!=m_particles.end();++flit) {
    msg.Tracking()<<"Construct decays for "<<(*flit)<<endl;
    skippit = 0;
    dmit    = m_decays.find((*flit));
    if (dmit->first==(*flit)) {
      if (dmit->second->IsEvaluated()) skippit = 1;
    }
    if (!skippit) {
      vertexlist = m_vertextable[(*flit)];
      if (vertexlist.size()==0) {
	msg.Error()<<"Error in Decay_Handler::BinaryDecays()."<<endl
		   <<"   Zero-length vertex list. Abort"<<endl;
	abort();
      }
      dt = NULL;
      if (!dt) dt = new Full_Decay_Table((*flit));

      for (int i=0;i<vertexlist.size();i++) {
	for (int j=0;j<3;j++) flavs[j] = vertexlist[i]->in[j];
	if (flavs[0].Mass()>flavs[1].Mass()+flavs[2].Mass()) {
	  dc = new Full_Decay_Channel((*flit));
	  for (int j=1;j<3;j++) dc->AddDecayProduct(flavs[j]);
	  dc->CreateDecay();
	  dt->AddDecayChannel(dc);
	  msg.Tracking()<<"Added Decay_Channel :"<<endl;
	  if (msg.LevelIsTracking()) dc->Output();
	}
      }
      m_decays.insert(std::make_pair((*flit),dt));
    }
  }
}

void All_Decays::ThreeBodyDecays()
{
  Vertex_List           primarylist,secondarylist;
  Flavour              primaries[3],flavs[3];
  Full_Decay_Table   * dt;
  Full_Decay_Channel * dc;
  bool                 skippit;
  DMIterator           dmit;
  for (FlSetIter flit=m_particles.begin();flit!=m_particles.end();++flit) {
    skippit = 0;
    dmit    = m_decays.find((*flit));
    dt      = NULL;
    if (dmit->first==(*flit)) {
      if (dmit->second->IsEvaluated()) skippit = 1;
                                  else dt = dmit->second;
    }
    if (!skippit) {
      primarylist = m_vertextable[(*flit)];
      if (primarylist.size()==0) {
	msg.Error()<<"Error in Decay_Handler::ThreeBodyDecays()."<<endl
		   <<"   Zero-length vertex list. Abort"<<endl;
	abort();
      }
      if (!dt) {
	msg.Tracking()<<"Funny finding in Decay_Handler::ThreeBodyDecays()."<<endl
		      <<"   Have to instantiate new Full_Decay_Table for "<<(*flit)<<endl;
	dt = new Full_Decay_Table((*flit));
      }
      for (int i=0;i<primarylist.size();i++) {
	for (int j=0;j<3;j++) primaries[j] = primarylist[i]->in[j];
	if (primaries[0].Mass()<primaries[1].Mass()+primaries[2].Mass() &&
	    primaries[0]!=primaries[1] && primaries[0]!=primaries[2]) {
	  for (int k=1;k<3;k++) {
	    secondarylist = m_vertextable[primaries[k]];
	    flavs[2]      = primaries[3-k];
	    if (secondarylist.size()==0) continue;
	    for (int l=0;l<secondarylist.size();l++) {
	      flavs[0] = secondarylist[l]->in[1];
	      flavs[1] = secondarylist[l]->in[2];
	      if (!(primaries[0].Mass()<flavs[0].Mass()+flavs[1].Mass()+flavs[2].Mass())) {
		dc = new Full_Decay_Channel((*flit));
		for (int m=0;m<3;m++) dc->AddDecayProduct(flavs[m]);
		dc->CreateDecay();
		dt->AddDecayChannel(dc);
	      }
	    }
	  }
	}
      }
    }
  }
}

void All_Decays::ArrangeDecays()
{
  for (DMIterator dmit=m_decays.begin();dmit!=m_decays.end();dmit++) dmit->second->ArrangeDecays();
}

bool All_Decays::InitializeDecays() { 
  bool okay = 1;
  for (DMIterator dmit=m_decays.begin();dmit!=m_decays.end();dmit++) {
    if (!dmit->second->IsEvaluated()) okay = okay && dmit->second->InitAllDecays(p_model,p_top);
  }
  msg.Tracking()<<"All_Decays::InitializeDecays() ";
  if (okay) msg.Tracking()<<" successful."<<endl;
  else {
    msg.Tracking()<<" failed."<<endl;
    msg.Error()<<"Some libraries were missing ! Type make install and rerun."<<endl;
    abort();
  }
  return okay; 
}

