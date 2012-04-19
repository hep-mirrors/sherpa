#include "SHRiMPS/Beam_Remnants/Colour_Reconnections.H"
#include "SHRiMPS/Tools/MinBias_Parameters.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Run_Parameter.H"

using namespace SHRIMPS;
using namespace ATOOLS;
using namespace std;

Colour_Reconnections::Colour_Reconnections() :
  m_on(MBpars.ReconnMode()!=reconn_mode::off),
  m_reconn(MBpars("ReconnProb")), 
  m_Q02(MBpars("QRC2")), m_b02(4.*m_Q02*sqr(rpa->hBar()*rpa->c())),
  m_eta(2.),
  m_ycut(MBpars("originalY")-MBpars("deltaY")),
  m_analyse(true)
{
  if (m_analyse) {
    m_histomap[string("Reconn_MassBefore")] = new Histogram(0,0.0,2000.0,2000);
    m_histomap[string("Reconn_MassAfter")]  = new Histogram(0,0.0,2000.0,2000);
  }
}

Colour_Reconnections::~Colour_Reconnections() 
{
  if (m_analyse) {
    Histogram * histo;
    string name;
    for (map<string,Histogram *>::iterator hit=m_histomap.begin();
       hit!=m_histomap.end();hit++) {
      histo = hit->second;
      name  = string("Ladder_Analysis/")+hit->first+string(".dat");
      histo->Finalize();
      histo->Output(name);
      delete histo;
    }
    m_histomap.clear();
  }
}
  
bool Colour_Reconnections::
FinishConfiguration(Blob_List * blobs,const double & smin) {
  m_shuffled = false;
  m_smin     = (smin<0. || MBpars.ReconnMode()==reconn_mode::fix)?m_Q02:smin;
  m_newcols.clear();
  m_colours.clear();
  m_trips.clear();
  m_antis.clear();
  m_links.clear();
  m_pairs.clear();
  HarvestParticles(blobs);
  FillWeightTable();
  ShuffleColours();
  blobs->push_back(AddReconnectionBlob());
  return true;
}

void Colour_Reconnections::HarvestParticles(Blob_List * blobs) {
  Blob * blob;
  Particle * part;
  for (Blob_List::iterator bit=blobs->begin();bit!=blobs->end();bit++) {
    blob = (*bit);
    if (blob->Has(blob_status::needs_hadronization)) {
      for (int i=0;i<blob->NOutP();i++) {
	part = blob->OutParticle(i);
	if (dabs(part->Momentum().Y())>m_ycut) part->SetInfo('B');
	if (part->Status()==part_status::active &&
	    part->DecayBlob()==NULL) {
	  unsigned int col1 = part->GetFlow(1);
	  unsigned int col2 = part->GetFlow(2);
	  colpair cols = colpair(col1,col2);
	  m_newcols[part] = cols;
	  if (col1!=0) {
	    m_trips.insert(part);
	    m_colours.insert(col1);
	  }
	  if (col2!=0) {
	    m_antis.insert(part);
	  }
	}
      }
      blob->UnsetStatus(blob_status::needs_beams);
      blob->UnsetStatus(blob_status::needs_showers);
      blob->UnsetStatus(blob_status::needs_harddecays);
      blob->UnsetStatus(blob_status::needs_hadronization);
    }
  }
}

void Colour_Reconnections::FillWeightTable() {
  set<Particle *>::iterator tripit, antiit;
  Particle * trip, * anti;
  for (tripit=m_trips.begin();tripit!=m_trips.end();tripit++) {
    trip = (*tripit);
    map<double,Particle *> dists;
    for (antiit=m_antis.begin();antiit!=m_antis.end();antiit++) {
      anti = (*antiit);
      if (anti==trip) continue;
      double dist(Distance(trip,anti));
      switch (ColourConnected(trip,anti)) {
      case 2:
	dist *= m_reconn;
	break;
      case 1:
	dist /= m_reconn;
	break;
      case 0:
      default:
	break;
      }
      dists[dist] = anti;
    }
    m_links[trip] = dists;
  }
}

void Colour_Reconnections::ShuffleColours() {
  OutputWeightTable();
  map<Particle *,map<double, Particle *> >::iterator mapit;
  map<double,Particle *>           dists;
  map<double,Particle *>::iterator distit;
  Particle * test1, * test2, * trip, * anti;
  while (!m_trips.empty()) {
    double maxdist = 0.;
    msg_Out()<<"Start looping to look for next colour connection: "
	     <<m_trips.size()<<" particles still to do.\n";
    for (mapit=m_links.begin();mapit!=m_links.end();mapit++) {
      test1 = mapit->first;
      if (m_trips.find(test1)==m_trips.end()) continue;
      dists = mapit->second;
      distit = dists.begin();
      while (distit!=dists.end()) {
	test2 = distit->second;
	if (m_antis.find(test2)!=m_antis.end()) {
	  if (distit->first>maxdist) {
	    trip    = test1;
	    anti    = test2;
	    maxdist = distit->first;
	    msg_Out()<<"     candidate  (dist = "<<maxdist<<") "
		     <<"["<<test1->Number()<<" "<<test2->Number()<<"]\n";
	  }
	  break;
	}
	distit++;
      }
    }
    if (trip==NULL || anti==NULL) {
      msg_Error()<<"Error in "<<METHOD<<":\n"
		 <<"   did not find a viable pair!\n";
      exit(1);
    }
    msg_Out()<<"   * want to establish connection between "
	     <<"["<<trip->Number()<<"]"
	     <<"("<<trip->GetFlow(1)<<", "<<trip->GetFlow(2)<<") and "
	     <<"["<<anti->Number()<<"]"
	     <<"("<<anti->GetFlow(1)<<", "<<anti->GetFlow(2)<<"), "
	     <<"dist = "<<maxdist<<".\n";
    m_trips.erase(trip);
    m_antis.erase(anti);    
    m_pairs.push_back(partpair(trip,anti));
    msg_Out()<<"   * now "<<m_trips.size()<<" / "<<m_antis.size()<<" "
	     <<"particles left for triplet/antitriplet.\n";
    if (m_trips.size()==1 &&
	(*m_trips.begin())==(*m_antis.begin())) {
      msg_Out()<<"Would have to save last gluon.\n"
	       <<(**m_trips.begin())<<"\n";
      SaveLastGluon((*m_trips.begin()));
      exit(1);
    }
  }
}

void Colour_Reconnections::SaveLastGluon(Particle * part) {
  partpair winner;
  double combdist(1.e99), testdist;
  Particle *test1, * test2, * trip, * anti;
  partdists pdists(m_links[part]), tdists;
  for (std::list<partpair>::iterator ppit=m_pairs.begin();
       ppit!=m_pairs.end();ppit++) {
    test1 = ppit->first;
    test2 = ppit->second;
    testdist = 0.;
    tdists = m_links[test1];
    for (partdists::iterator tit=tdists.begin();tit!=tdists.end();tit++) {
      if (tit->second==part) {
	testdist += tit->first;
	break;
      }
    }
    for (partdists::iterator pit=pdists.begin();pit!=pdists.end();pit++) {
      if (pit->second==test2) {
	testdist += pit->first;
	break;
      }
    }
    if (testdist<combdist) {
      trip = test1;
      anti = test2;
      combdist = testdist;
    }
  }
  msg_Out()<<"Would like to insert "<<part->Number()<<" between "
	   <<"["<<trip->Number()<<" and "<<anti->Number()<<"], "
	   <<"comb = "<<combdist<<".\n";
}

const double Colour_Reconnections::
Distance(Particle * part1,Particle * part2,const bool & spat) {
  double sij((part1->Momentum()+part2->Momentum()).Abs2());
  double dist(pow((m_smin+sij)/m_smin,m_eta));
  if (spat && part1->ProductionBlob()!=part2->ProductionBlob()) {
    double deltar = (part1->ProductionBlob()->Position().Perp()-
		     part2->ProductionBlob()->Position().Perp()).Abs2();
    dist  *= exp(deltar/m_b02);
  }
  return dist;
}

Blob * Colour_Reconnections::AddReconnectionBlob() {
  Blob * blob = new Blob();
  blob->SetType(btp::Soft_Collision);
  blob->SetTypeSpec("Four_Momentum_Compensation");
  blob->SetId();
  blob->SetStatus(blob_status::needs_hadronization);
  Particle * partin, * partout;
  map<Particle *,colpair>::iterator pcit;
  for (map<Particle *,map<double,Particle *> >::iterator pit=m_links.begin();
       pit!=m_links.end();pit++) {
    partin = pit->first;
    partin->SetStatus(part_status::decayed);
    blob->AddToInParticles(partin);
    pcit = m_newcols.find(partin);
    if (pcit==m_newcols.end()) {
      msg_Error()<<"Error in "<<METHOD<<":\n"
		 <<"   Did not find particle ["<<partin->Number()<<"] "
		 <<"in new colours list.\n"
		 <<"   Will exit the run.\n";
      exit(1);
    }
    partout = new Particle(0,partin->Flav(),partin->Momentum(),partin->Info());
    partout->SetFlow(1,pcit->second.first);
    partout->SetFlow(2,pcit->second.second);
    partout->SetNumber();
    blob->AddToOutParticles(partout);
  }
  return blob;
}

void Colour_Reconnections::OutputWeightTable() {
  for (map<Particle *,map<double, Particle *> >::iterator mapit=m_links.begin();
       mapit!=m_links.end();mapit++) {
    msg_Out()<<"Links for particle ["<<mapit->first->Number()<<"]"
	     <<"("<<mapit->first->GetFlow(1)<<", "
	     <<mapit->first->GetFlow(2)<<"):\n";
    map<double,Particle *> dists = mapit->second;
    for (map<double, Particle *>::iterator distit=dists.begin();
	 distit!=dists.end();distit++) {
      msg_Out()<<"   "<<distit->first<<"     "
	       <<"["<<distit->second->Number()<<"]"
	       <<"("<<distit->second->GetFlow(1)<<", "
	       <<distit->second->GetFlow(2)<<")\n";
    }
  }
}

