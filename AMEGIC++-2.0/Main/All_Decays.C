#include "All_Decays.H"

using namespace AMEGIC;
using namespace APHYTOOLS;
using namespace AORGTOOLS;
using namespace std;

All_Decays::All_Decays(Interaction_Model_Base * _model,Topology * _top) :
  p_model(_model), p_top(_top)
{
  Vertex * vertex = p_model->GetVertex();

  for (int i=0;i<vertex->MaxNumber();++i) {
    if ((*vertex)[i]->on) {
      m_vertextable[(*vertex)[i]->in[0]].push_back((*vertex)[i]);
    }
  }
  FindUnstableParticles();
  InitializeDecayTables();
}

void All_Decays::FindUnstableParticles() 
{
  int count = 0;
  Flavour flav;
  for (int i=1;i<MAX_PARTICLES;i++) {
    flav = Flavour(particles[i].kfc);
    if (!flav.IsHadron() && !flav.IsStable() && flav.IsOn() && flav.Size()==1) {
      if (CheckInVertex(flav)==1) {
	if (m_particles.size()==0) {
	  m_particles.push_back(flav);
	}
	else if (m_particles.size()==1) {
	  if (flav.PSMass()>m_particles[0].PSMass()) m_particles.push_back(flav);
	                                        else m_particles.push_front(flav);
	}
	else if (flav.PSMass()>m_particles.back().PSMass()) m_particles.push_back(flav);
	else {
	  double compare=0.;
	  for (DecayingParticleIterator dit=m_particles.begin();dit!=m_particles.end();++dit) {
	    if (flav.PSMass()>compare && flav.PSMass()<(*dit).PSMass()) {
	      m_particles.insert(dit,flav);
	      break;
	    }
	  }
	}
      }
    }
  }
  //msg.Debugging()<<"Final list :"<<endl;
  //for (DecayingParticleIterator dit=m_particles.begin();dit!=m_particles.end();++dit) msg.Debugging()<<(*dit)<<endl;
}


bool All_Decays::CheckInVertex(Flavour flav)
{
  VertexList & vertexlist = m_vertextable[flav];
  if (vertexlist.size()>0) return 1;
  return 0;
}


void All_Decays::InitializeDecayTables() {
  BinaryDecays();
}

void All_Decays::BinaryDecays()
{
  Flavour flav[3];
  VertexList vertexlist;
  DecayList  decaylist;
  for (DecayingParticleIterator dit=m_particles.begin();dit!=m_particles.end();++dit) {
    flav[0]    = (*dit);
    vertexlist = m_vertextable[flav[0]];
    if (vertexlist.size()==0) {
      msg.Error()<<"Error in Decay_Handler::BinaryDecays()."<<endl
		 <<"   Zero-length vertex list. Abort"<<endl;
      abort();
    }
    for (int i=0;i<vertexlist.size();i++) {
      flav[1] = vertexlist[i]->in[1];
      flav[2] = vertexlist[i]->in[2];
      if (flav[0].PSMass()>flav[1].PSMass()+flav[2].PSMass()) {
	decaylist.push_back(new Single_Process(1,2,flav));
      }
    }
    m_decaytable.insert(make_pair(flav[0],decaylist));
    decaylist.clear();
  }
}


bool All_Decays::InitializeDecays()
{
  DecayList                decaylist;
  DecayingParticleIterator dit = m_particles.begin();

  bool okay = 1;
  while ((*dit)) {
    decaylist = m_decaytable[(*dit)];
    for (int i=0;i<decaylist.size();i++) {
      okay = okay && decaylist[i]->InitAmplitude(p_model,p_top);
    }
    dit++;
  }
  if (!okay) {
    msg.Error()<<"Initialization of the decays failed !"<<endl
	       <<"   Try to continue with setup of integrators."<<endl; 
  }

  dit = m_particles.begin();
  while ((*dit)) {
    decaylist = m_decaytable[(*dit)];
    for (int i=0;i<decaylist.size();i++) {
      okay = okay && decaylist[i]->SetUpIntegrator();
    }
    dit++;
  }
}

bool All_Decays::CalculateBranchingWidths(string _resdir) {}
