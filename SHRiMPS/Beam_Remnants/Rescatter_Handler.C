#include "SHRiMPS/Beam_Remnants/Rescatter_Handler.H"
#include "SHRiMPS/Tools/MinBias_Parameters.H"
#include "MODEL/Main/Model_Base.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Phys/Momenta_Stretcher.H"

using namespace SHRIMPS;
using namespace ATOOLS;
using namespace std;

Rescatter_Handler::Rescatter_Handler(Beam_Remnant_Handler * beams) :
  m_rescatter(MBpars.RescMode()!=resc_mode::off), m_mustmatch(false), 
  m_rescprob(MBpars("RescProb")),
  m_singprob(MBpars.RescOverSing()==resc_over_sing::on?MBpars("RescProb1"):0.), 
  p_beams(beams),
  p_alphaS(static_cast<MODEL::Strong_Coupling *>
	   (MODEL::s_model->GetScalarFunction(string("strong_cpl")))),
  m_Ylimit(MBpars("originalY")-MBpars("deltaY")),
  m_analyse(false)
{ 
  if (m_analyse) {
    m_histomap[string("Rescatter_B")]  = new Histogram(0,0.0,20.0,200);
    m_histomap[string("Rescatter_b1")] = new Histogram(0,0.0,20.0,200);
    m_histomap[string("Rescatter_b2")] = new Histogram(0,0.0,20.0,200);
    m_histomap[string("Rescatter_wt")] = new Histogram(0,0.0,20.0,200);
    m_histomap[string("Rescatter_Y")]  = new Histogram(0,-5.0,5.0,20);
    m_histomap[string("Rescatter_E")]  = new Histogram(0,0.0,2000.0,200);
    m_histomap[string("Rescatter_N")]  = new Histogram(0,0.0,10.0,10);
  }
}
  
Rescatter_Handler::~Rescatter_Handler() {
  if (!m_histomap.empty()) {
    Histogram * histo;
    string name;
    for (map<string,Histogram *>::iterator 
	   hit=m_histomap.begin();hit!=m_histomap.end();hit++) {
      histo = hit->second;
      name  = string("Ladder_Analysis/")+hit->first+string(".dat");
      histo->Finalize();
      histo->Output(name);
      delete histo;
    }
    m_histomap.clear();
  }
}

void Rescatter_Handler::
ResetCollision(Omega_ik * eikonal,const double & smin,const double & B) {
  if (!m_rescatter) return;
  msg_Tracking()
    <<"###########################################"<<endl
    <<"###########################################"<<endl
    <<METHOD<<": "
    <<"particles: "<<m_particles.size()<<", "
    <<"part pairs: "<<m_probpairs.size()<<", "
    <<"blobs: "<<m_treatedblobs.size()<<endl
    <<"###########################################"<<endl
    <<"###########################################"<<endl;
  m_treatedblobs.clear();
  m_partmap.clear();
  m_intervals.clear();
  ResetRescatter(true);
  
  p_eikonalweights = eikonal->GetWeights();
  m_smin = smin;
  m_B    = B;
}

void Rescatter_Handler::ResetRescatter(const bool & enforce) { 
  if (!m_rescatter) return;
  m_Nfact = 1; m_Nresc = 0;
  m_particles.clear(); 
  m_probpairs.clear(); 
}

void Rescatter_Handler::UpdateCollision(Blob_List * blobs) {
  if (!m_rescatter) return;
  Blob * blob;
  bool addparts(false);
  for (size_t i=0;i<blobs->size();++i) {
    blob = (*blobs)[i];
    if (!blob->Has(blob_status::needs_beams) ||
	m_treatedblobs.find(blob)!=m_treatedblobs.end()) continue;
    else {
      m_treatedblobs.insert(blob);
      addparts = DealWithBlob(blob);
    }
  }
  if (addparts) {
    m_Nresc++;
    m_Nfact *= m_Nresc;
  }
}

bool Rescatter_Handler::DealWithBlob(ATOOLS::Blob * blob) {
  bool stretch(false);
  for (int i=0;i<blob->NOutP();i++) {
    if (blob->OutParticle(i)->Momentum().Abs2()<-1e-8) stretch=true;
  }
  if(stretch){
    Momenta_Stretcher momstretcher;
    if (!momstretcher.StretchBlob(blob)){
      msg_Error()<<"Error in "<<METHOD<<": "
		 <<"cannot adjust momenta to put all particles on-shell.\n";
    }   
  } 
  Vec4D pos(blob->Position()/(rpa->hBar()*rpa->c()));
  m_b1 = pos.PPerp();
  m_b2 = sqrt(m_B*m_B+m_b1*m_b1-2.*m_B*m_b1*pos.CosPhi());
  PartList plist;
  Particle * part;
  for (int i=0;i<blob->NOutP();i++) {
    part = blob->OutParticle(i);
    if (!part->DecayBlob() && part->Info()!='B') plist.push_back(part);
  }
  PSetYSet singlets;
  m_sorter.Sort(&plist,&singlets);
  if (singlets.size()>1) {
    double y1, y2;
    PSetYSet::iterator plit(singlets.begin());
    PSetYSet::iterator plend(singlets.end()); plend--;
    y1 = (*(*plit)->rbegin())->Momentum().Y();
    do {
      plit++;
      y2 = (*(*plit)->begin())->Momentum().Y();
      if (y1<y2) m_intervals.push_back(make_pair(y1,y2));
      y1 = (*(*plit)->rbegin())->Momentum().Y();
    } while (plit!=plend);
  }

  for (PSetYSet::iterator plit=singlets.begin();plit!=singlets.end();plit++) {
    while (!((*plit)->empty())) {
      AddParticleToRescatters(*(*plit)->begin());
      (*plit)->erase((*plit)->begin());
    }
    delete (*plit);
  }
  return true;
}

void Rescatter_Handler::AddParticleToRescatters(Particle * part) {
  if (!m_rescatter) return;
  double y1(part->Momentum().Y()),kt12(part->Momentum().PPerp2()),y2,kt22;
  double sup,s12,prob,ybar,ymax;
  int nbeam(dabs(y1)>m_Ylimit);
  double expo;
  bool allowed, singlet;
  for (set<Particle *, partcomp>::iterator piter=m_particles.begin();
       piter!=m_particles.end();piter++) {
    allowed = true;
    y2      = (*piter)->Momentum().Y();
    singlet = false;
    for (list<pair<double,double> >::iterator sit=m_intervals.begin();
	 sit!=m_intervals.end();sit++) {
      if ((y1<=sit->first && y2>=sit->second) ||
	  (y2<=sit->first && y1>=sit->second)) {
        singlet = true;
	continue;
      }
    }
    prob  = singlet?m_singprob:1.;
    if (prob<0.00000001 || !CanRescatter((*piter),part)) continue;
    s12   = Max(0.,((*piter)->Momentum()+part->Momentum()).Abs2());
    kt22  = (*piter)->Momentum().PPerp2();
    sup   = m_rescprob * SuppressionTerm(kt12,kt22);
    ybar  = (part->Momentum()+(*piter)->Momentum()).Y();
    ymax  = Max(dabs(y1),dabs((*piter)->Momentum().Y()));
    prob *= p_eikonalweights->RescatterProbability(m_b1,m_b2,y1,y2,sup,
						   nbeam+int(dabs(y2)>m_Ylimit)); 
    expo  = (singlet?0.:1.)+p_eikonalweights->EffectiveIntercept(m_b1,m_b2,ybar);
    prob *= pow(s12/Max(s12,m_smin),expo);
    prob /= sqrt(double(m_Nfact));
    // new stuff in the line below.
    prob *= 4.*exp(-ymax/2.-dabs(ybar));
    if (m_analyse) m_histomap["Rescatter_wt"]->Insert(m_B,prob);
    PartPair partpair;
    partpair.first  = y1<y2?part:(*piter);
    partpair.second = y1<y2?(*piter):part;
    m_probpairs[prob] = partpair;
  }

  m_particles.insert(part);
}


bool Rescatter_Handler::
IsColourConnected(const Particle * part1,const Particle * part2) const
{
  if (!(part1->GetFlow(1)==part2->GetFlow(2) && part1->GetFlow(1)!=0) &&
      !(part1->GetFlow(2)==part2->GetFlow(1) && part1->GetFlow(2)!=0))
    return false;
  return true;
}

bool Rescatter_Handler::
CanRescatter(const Particle * part1,const Particle * part2) const
{
  if (part1==part2) return false;
  if (m_mustmatch && !IsColourConnected(part1,part2)) return false;
  if (part1->Flav().IsQuark() && part2->Flav().IsQuark()) return false;
  return true;
}



bool Rescatter_Handler::
SelectRescatter(Particle *& part1,Particle *& part2) {
  if (!m_rescatter || m_probpairs.empty()) return false;
  while (!m_probpairs.empty()) {
    if (m_probpairs.begin()->first>ran->Get()) {
      part1   = m_probpairs.begin()->second.first;
      part2   = m_probpairs.begin()->second.second;
      DeleteProbPairs(part1,part2);
      if (m_analyse) {
	double Y((part1->Momentum()+part2->Momentum()).Y());
	double E(sqrt((part1->Momentum()+part2->Momentum()).Abs2()));
	m_histomap[string("Rescatter_B")]->Insert(m_B);
	m_histomap[string("Rescatter_b1")]->Insert(m_b1);
	m_histomap[string("Rescatter_b2")]->Insert(m_b2);
	m_histomap[string("Rescatter_Y")]->Insert(Y);
	m_histomap[string("Rescatter_E")]->Insert(E);
      }
      return true;
    }
    m_probpairs.erase(m_probpairs.begin());
  }
  if (m_analyse) m_histomap["Rescatter_N"]->Insert(m_Nresc);
  return false;
}


bool Rescatter_Handler::ConnectBlobs(Blob_List * blobs,Blob * add) {
  if (m_partmap.empty()) return true;
  bool       found;
  Blob     * blob;
  Particle * repl, * orig, * test, * test2;
  Particle_Vector partvec;
  while (!m_partmap.empty()) {
    orig  = m_partmap.begin()->first; 
    repl  = m_partmap.begin()->second;
    for (int i=0;i<2;i++) {
      blob  = repl->ProductionBlob();
      found = false;
      for (short int i=0;i<blob->NInP();i++) {
	test = blob->InParticle(i);
	if (test->ProductionBlob()==NULL &&
	    test->Info()=='I' &&
	    test->GetFlow(1)==orig->GetFlow(1) &&
	    test->GetFlow(2)==orig->GetFlow(2)) {
	  blob->AddToInParticles(orig);
	  orig->SetInfo('R');
	  orig->SetStatus(part_status::decayed);
	  delete blob->RemoveInParticle(test);
	  found = true;
	  if (!blob->CheckColour()) {
	    msg_Tracking()<<"Problem in "<<METHOD<<":\n"
			  <<"   Scattering blob ("<<blob->Id()<<") "
			  <<"seems fishy: "
			  <<"Bad colour configuration.\n"<<(*blob)<<"\n";
	    return false;
	  }
	  break;
	}
      }
      if (found) break;
      Particle * help=orig;orig=repl;repl=help;
    }
    if (!found && add) {
      blob  = repl->ProductionBlob();
      for (short int i=0;i<add->NInP();i++) {
	test = add->InParticle(i);
	if (test->Number()==repl->Number()) {
	  for (short int i=0;i<blob->NInP();i++) {
	    test2 = blob->InParticle(i);
	    if (test2->ProductionBlob()==NULL &&
		test2->Info()=='I' &&
		(repl->GetFlow(1)==0 || 
		 test2->GetFlow(1)==repl->GetFlow(1)+100) &&
		(repl->GetFlow(2)==0 ||
		 test2->GetFlow(2)==repl->GetFlow(2)+100)) {
	      add->AddToOutParticles(test2);
	      add->AddToInParticles(orig);
	      test2->SetInfo('R');
	      orig->SetInfo('R');
	      orig->SetStatus(part_status::decayed);
	      blobs->push_back(add);
	      delete add->RemoveInParticle(test);
	      if (!add->CheckColour()) {
		msg_Tracking()<<"Problem in "<<METHOD<<":\n"
			      <<"   Extra blob ("<<add->Id()<<") seems fishy: "
			      <<"Bad colour configuration.\n"<<(*add)<<"\n";
		return false;
	      }
	      found = true;
	      break;
	    }
	  }
	  if (!found) { 
	    msg_Tracking()<<"Error in "<<METHOD<<".\n"<<(*orig)<<"\n"<<(*repl)
			  <<"\n"<<(*add)<<"\n"<<(*blob)<<"\n"; 
	    if (add) {
	      add->DeleteOutParticles();
	      add->DeleteInParticles();
	      delete add;
	    }
	    return false;
	  }
	}
	if (found) break;
      }
    }
    if (!found) {
      msg_Tracking()<<"WARNING in "<<METHOD<<":\n"
		    <<"   Failed to connect "
		    <<orig->Number()<<" <--> "<<repl->Number()<<" in "
		    <<blobs->size()<<" blobs.\n";
      if (add) {
	add->DeleteOutParticles();
	add->DeleteInParticles();
	delete add;
      }
      return false;
    }
    m_partmap.erase(m_partmap.begin());
  }
  return true;
}

void Rescatter_Handler::
DeleteProbPairs(Particle *& part1,Particle *& part2) {
  ProbPairMap::iterator piter=m_probpairs.begin(), help;
  while (piter!=m_probpairs.end()) {
    if (piter->second.first ==part1 ||
	piter->second.second==part1 ||
	piter->second.first ==part2 ||
	piter->second.second==part2) {
      help = piter; help++;
      m_probpairs.erase(piter);
      piter = help;
    }
    else piter++;
  }
  m_particles.erase(part1);
  m_particles.erase(part2);
}


void Rescatter_Handler::
FillInitialStateIntoBlob(ATOOLS::Blob * blob,Ladder * ladder) {
  blob->AddToInParticles(ladder->GetIn1()->GetParticle());
  blob->AddToInParticles(ladder->GetIn2()->GetParticle());
}

double Rescatter_Handler::
SuppressionTerm(const double & q02,const double & q12) {
  return sqrt(p_alphaS->Weight(q02,true)*p_alphaS->Weight(q12,true));
}
