#include "SINIC++/Beam_Remnants/Colour_Reconnections.H"
#include "SINIC++/Tools/MinBias_Parameters.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Run_Parameter.H"

using namespace SINIC;
using namespace ATOOLS;
using namespace std;

Colour_Reconnections::Colour_Reconnections() :
  m_on(MBpars.ReconnMode()!=reconn_mode::off),
  m_colfac(1./(64.-1.)), m_reconn(MBpars("ReconnProb")), 
  m_Q02(MBpars("QRC2")), m_b02(4.*m_Q02*sqr(rpa->hBar()*rpa->c())),
  m_inveta(-2.),
  m_ycut(MBpars("originalY")-2.*MBpars("deltaY")),
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
  if (!m_on) return true;
  m_shuffled = false;
  m_smin     = (smin<0. || MBpars.ReconnMode()==reconn_mode::fix)?m_Q02:smin;
  m_pclist.clear();
  HarvestParticles(blobs);
  ShuffleColours();
  if (m_shuffled) blobs->push_back(AddReconnectionBlob());
  return true;
}

void Colour_Reconnections::HarvestParticles(Blob_List * blobs) {
  Blob * blob;
  Particle * part;
  PartList parts;
  for (Blob_List::iterator bit=blobs->begin();bit!=blobs->end();bit++) {
    blob = (*bit);
    if (blob->Has(blob_status::needs_hadronization)) {
      for (int i=0;i<blob->NOutP();i++) {
	part = blob->OutParticle(i);
	if (dabs(part->Momentum().Y())>m_ycut) part->SetInfo('B');
	if (part->Status()==part_status::active &&
	    part->DecayBlob()==NULL) {
	  parts.push_back(part);
	}
      }
      blob->UnsetStatus(blob_status::needs_beams);
      blob->UnsetStatus(blob_status::needs_showers);
      blob->UnsetStatus(blob_status::needs_harddecays);
      blob->UnsetStatus(blob_status::needs_hadronization);
    }
  }
  m_sorter.Sort(&parts,&m_pclist);
}


void Colour_Reconnections::ShuffleColours() {
  PCList::iterator pit1(m_pclist.begin()), pit2, pit3, pit4;
  while (pit1!=m_pclist.end()) {
    pit2 = pit1; pit2++;
    if (pit1->second.first!=pit2->second.second ||
	pit1->second.first==0) { 
      pit1++;
      continue;
    }
    if (pit2==m_pclist.end()) break;
    pit3 = pit1; pit3++;
    while (pit3!=m_pclist.end()) {
      pit4 = pit3; pit4++;
      if (pit4==m_pclist.end()) break;
      if ((pit3->first!=pit2->first) &&
	  !(pit1->second.second==pit4->second.first && 
	    pit1->second.second!=0) &&
	  (pit3->second.first==pit4->second.second && 
	   pit3->second.first!=0)) {
	double w12(Weight(pit1->first,pit2->first));
	double w34(Weight(pit3->first,pit4->first));
	double w14(Weight(pit1->first,pit4->first));
	double w32(Weight(pit3->first,pit2->first));
	double w1234(w12*w34),w1432(m_colfac*m_reconn*w14*w32);
	double summed(w1234+w1432);
	double m12(sqrt((pit1->first->Momentum()+
			 pit2->first->Momentum()).Abs2()));
	double m34(sqrt((pit3->first->Momentum()+
			 pit4->first->Momentum()).Abs2()));
	if (m_analyse) {
	  m_histomap[string("Reconn_MassBefore")]->Insert(m12);
	  m_histomap[string("Reconn_MassBefore")]->Insert(m34);
	}
	if (w1432>summed*ran->Get()) {
	  SkewList(pit1,pit2,pit3,pit4);
	  if (m_analyse) {
	    double m14(sqrt((pit1->first->Momentum()+
			     pit4->first->Momentum()).Abs2()));
	    double m32(sqrt((pit3->first->Momentum()+
			     pit2->first->Momentum()).Abs2()));
	    m_histomap[string("Reconn_MassAfter")]->Insert(m14);
	    m_histomap[string("Reconn_MassAfter")]->Insert(m32);
	  }
	  break;
	}
	else {
	  if (m_analyse) {
	    m_histomap[string("Reconn_MassAfter")]->Insert(m12);
	    m_histomap[string("Reconn_MassAfter")]->Insert(m34);
	  }
	}
      }
      pit3++;
    }
    pit1++;
  }
}

void Colour_Reconnections::
SkewList(PCList::iterator & pit1,PCList::iterator & pit2,
	 PCList::iterator & pit3,PCList::iterator & pit4) {
  m_shuffled = true;
  PCList::iterator start(pit2), stop(pit4), test(stop), test1;
  bool ring(true);
  while (test!=m_pclist.end()) {
    if (test->second.first==0) {
      ring = false;
      break;
    }
    test++;
  }
  //msg_Out()<<METHOD<<"(ring = "<<ring<<") for "
  //	   <<m_pclist.size()<<" particles.\n";
  unsigned int help(pit2->second.second);
  pit2->second.second = pit4->second.second;
  pit4->second.second = help;
  PCList helplist;

  if (!ring) {
    while (start!=stop) {
      helplist.push_back((*start));
      start = m_pclist.erase(start);
    } 
    size_t count(helplist.size());
    //msg_Out()<<METHOD<<" must reorder "<<count<<" elements.\n";
    while (count>0 && helplist.begin()->second.second!=0) {
      helplist.push_back((*helplist.begin()));
      helplist.pop_front();
      count--;
    }
    while (!helplist.empty()) {
      m_pclist.push_back((*helplist.begin()));
      helplist.pop_front();
    }
  }
  else {
    test  = stop;
    //msg_Out()<<"try to extract ring and rotate it, start with: "
    //	     <<test->first->Number()<<".\n";
    while (test!=m_pclist.end()) {
      //msg_Out()<<"  add (+): "<<test->first->Number()<<".\n";
      helplist.push_back((*test));
      test = m_pclist.erase(test);
      if (test==m_pclist.end() ||
	  test->second.second!=helplist.back().second.first) break;
    }
    test--;
    while (helplist.back().second.first==test->second.second &&
	   test!=start) {
      helplist.push_back((*test));
      //msg_Out()<<"  add (-): "<<test->first->Number()<<", "
      //       <<"check col: "<<helplist.back().second.first<<".\n";
      test = m_pclist.erase(test);
      test--;
    }
    //msg_Out()<<"check treatment of ring: "<<helplist.size()<<" members "
    //	     <<"for pos:"<<start->first->Number()<<".\n";
    //for (PCList::iterator pit=helplist.begin();
    //	 pit!=helplist.end();pit++) {
    //  msg_Out()<<(*(pit->first))
    //	       <<"["<<pit->second.first<<", "<<pit->second.second<<"]\n";
    //}
    //msg_Out()<<"---------------------------------------------------\n";
    m_pclist.splice(start,helplist);
    //msg_Out()<<"splicing successful before "<<start->first->Number()<<".\n";
  }
  
  //msg_Out()<<"Output "<<m_pclist.size()<<" particles.\n";
  //for (PCList::iterator pit=m_pclist.begin();
  //   pit!=m_pclist.end();pit++) {
  //msg_Out()<<(*(pit->first))
  //	     <<"["<<pit->second.first<<", "<<pit->second.second<<"]\n";
  //}
}
  
const double Colour_Reconnections::
Weight(ATOOLS::Particle * part1,ATOOLS::Particle * part2,const bool & spatial) {
  double weight(pow(1.+(part1->Momentum()+part2->Momentum()).Abs2()/m_smin,
		    m_inveta));
  double dist2(0.);
  if (spatial && part1->ProductionBlob()!=part2->ProductionBlob()) {
    dist2 = (part1->ProductionBlob()->Position().Perp()-
	     part2->ProductionBlob()->Position().Perp()).Abs2();
    weight *= exp(-dist2/m_b02);
  }
  return weight;
}

Blob * Colour_Reconnections::AddReconnectionBlob() {
  Blob * blob = new Blob();
  blob->SetType(btp::Soft_Collision);
  blob->SetTypeSpec("Four_Momentum_Compensation");
  blob->SetId();
  blob->SetStatus(blob_status::needs_hadronization);
  Particle * partin, * partout;
  for (PCList::iterator pit=m_pclist.begin();pit!=m_pclist.end();pit++) {
    partin = pit->first;
    blob->AddToInParticles(partin);
    partin->SetStatus(part_status::decayed);
    partout = new Particle(0,partin->Flav(),partin->Momentum(),partin->Info());
    partout->SetFlow(1,pit->second.first);
    partout->SetFlow(2,pit->second.second);
    partout->SetNumber();
    blob->AddToOutParticles(partout);
  }
  return blob;
}
