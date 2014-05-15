#include "SHERPA/Tools/HepMC2_Interface.H"
#ifdef USING__HEPMC2

#include "ATOOLS/Phys/Blob_List.H"
#include "ATOOLS/Phys/Particle.H"
#include "ATOOLS/Phys/NLO_Subevt.H"
#include "ATOOLS/Math/Vector.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Exception.H"
#include "MODEL/Main/Model_Base.H"
#include "PHASIC++/Main/Phase_Space_Handler.H"

#include "HepMC/GenEvent.h"
#include "HepMC/GenVertex.h"
#include "HepMC/GenParticle.h"
#include "HepMC/SimpleVector.h"
#include "HepMC/PdfInfo.h"
#ifdef USING__HEPMC2__UNITS
#include "HepMC/Units.h"
#endif

using namespace SHERPA;
using namespace ATOOLS;

HepMC2_Interface::HepMC2_Interface():
  p_event(NULL)
{
}

HepMC2_Interface::~HepMC2_Interface()
{
  delete p_event;
  DeleteGenSubEventList();
}

bool HepMC2_Interface::Sherpa2HepMC(ATOOLS::Blob_List *const blobs,
                                    HepMC::GenEvent& event, double weight)
{
#ifdef USING__HEPMC2__UNITS
  event.use_units(HepMC::Units::GEV,
                  HepMC::Units::MM);
#endif
  event.set_event_number(ATOOLS::rpa->gen.NumberOfGeneratedEvents());
  size_t decid(11);
  std::map<size_t,size_t> decids;
  Blob *sp(blobs->FindFirst(btp::Signal_Process));
  if (sp) {
    Blob_Data_Base *info((*sp)["Decay_Info"]);
    if (info) {
      DecayInfo_Vector decs(info->Get<DecayInfo_Vector>());
      for (size_t i(0);i<decs.size();++i) decids[decs[i]->m_id]=++decid;
    }
  }
  m_blob2genvertex.clear();
  m_particle2genparticle.clear();
  HepMC::GenVertex * vertex;
  std::vector<HepMC::GenParticle*> beamparticles;
  for (ATOOLS::Blob_List::iterator blit=blobs->begin();
       blit!=blobs->end();++blit) {
    if (Sherpa2HepMC(*(blit),vertex,decids)) {
      event.add_vertex(vertex);
      if ((*blit)->Type()==ATOOLS::btp::Signal_Process) {
        if ((**blit)["NLO_subeventlist"]) {
          THROW(fatal_error,"Events containing correlated subtraction events"
                +std::string(" cannot be translated into the full HepMC event")
                +std::string(" format.\n")
                +std::string("   Try 'EVENT_OUTPUT=HepMC_Short' instead."));
        }
        event.set_signal_process_vertex(vertex);
        Blob_Data_Base *pdfinfodb((**blit)["PDFInfo"]);
        if (pdfinfodb && (*blit)->NInP()==2) {
          PHASIC::PDF_Info pi=pdfinfodb->Get<PHASIC::PDF_Info>();
          double q(sqrt(sqrt(pi.m_muf12*pi.m_muf22)));
          HepMC::PdfInfo pdfinfo(pi.m_fl1,pi.m_fl2,pi.m_x1,pi.m_x2,
                                 q,pi.m_xf1,pi.m_xf2);
          event.set_pdf_info(pdfinfo);
        }
      }
      else if ((*blit)->Type()==ATOOLS::btp::Beam || 
	       (*blit)->Type()==ATOOLS::btp::Bunch) {
        for (HepMC::GenVertex::particles_in_const_iterator 
	       pit=vertex->particles_in_const_begin();
             pit!=vertex->particles_in_const_end(); ++pit) {
          if ((*pit)->production_vertex()==NULL) {
            beamparticles.push_back(*pit);
          }
        }
      }
    }
  }
  if (beamparticles.size()==2) {
    event.set_beam_particles(beamparticles[0],beamparticles[1]);
  }
  std::vector<double> weights;
  weights.push_back(weight);
  if (sp) {
    Blob_Data_Base *info((*sp)["MEWeight"]);
    if (!info) THROW(fatal_error,"Missing weight info.");
    double meweight(info->Get<double>());
    weights.push_back(meweight);
    Blob_Data_Base *ofni((*sp)["Weight_Norm"]);
    if (!ofni) THROW(fatal_error,"Missing weight normalisation.");
    double weightnorm(ofni->Get<double>());
    weights.push_back(weightnorm);
    ofni=(*sp)["Trials"];
    if (!ofni) THROW(fatal_error,"Missing nof trials.");
    double trials(ofni->Get<double>());
    weights.push_back(trials);
  }
  event.weights()=weights;
  return true;
}

bool HepMC2_Interface::Sherpa2ShortHepMC(ATOOLS::Blob_List *const blobs,
                                         HepMC::GenEvent& event, double weight)
{
#ifdef USING__HEPMC2__UNITS
  event.use_units(HepMC::Units::GEV,
                  HepMC::Units::MM);
#endif
  event.set_event_number(ATOOLS::rpa->gen.NumberOfGeneratedEvents());
  HepMC::GenVertex * vertex=new HepMC::GenVertex();
  std::vector<HepMC::GenParticle*> beamparticles;
  std::vector<std::pair<HepMC::FourVector,int> > beamparts,
                                                 remnantparts1, remnantparts2;
  Blob *sp(blobs->FindFirst(btp::Signal_Process));
  NLO_subevtlist* subevtlist(NULL);
  ME_wgtinfo* wgtinfo(0);
  
  if (sp) {
    Blob_Data_Base* seinfo=(*sp)["ME_wgtinfo"];
    if (seinfo) wgtinfo=seinfo->Get<ME_wgtinfo*>();
   
    Blob_Data_Base * bdb((*sp)["NLO_subeventlist"]);
    if (bdb) subevtlist=bdb->Get<NLO_subevtlist*>();
  }
  for (ATOOLS::Blob_List::iterator blit=blobs->begin();
       blit!=blobs->end();++blit) {
    Blob* blob=*blit;
    for (int i=0;i<blob->NInP();i++) {
      if (blob->InParticle(i)->ProductionBlob()==NULL) {
        Particle* parton=blob->InParticle(i);
        ATOOLS::Vec4D mom  = parton->Momentum();
        HepMC::FourVector momentum(mom[1],mom[2],mom[3],mom[0]);
        HepMC::GenParticle* inpart = 
	  new HepMC::GenParticle(momentum,parton->Flav().HepEvt(),2);
        vertex->add_particle_in(inpart);
        // distinct because SHRIMPS has no bunches for some reason
        if (blob->Type()==btp::Beam || blob->Type()==btp::Bunch) {
          beamparticles.push_back(inpart);
          beamparts.push_back(std::make_pair(momentum,parton->Flav().HepEvt()));
        }
      }
    }
    for (int i=0;i<blob->NOutP();i++) {
      if (blob->OutParticle(i)->DecayBlob()==NULL) {
        Particle* parton=blob->OutParticle(i);
        ATOOLS::Vec4D mom  = parton->Momentum();
        HepMC::FourVector momentum(mom[1],mom[2],mom[3],mom[0]);
        HepMC::GenParticle* outpart = 
	  new HepMC::GenParticle(momentum,parton->Flav().HepEvt(),1);
        vertex->add_particle_out(outpart);
        if (blob->Type()==btp::Beam) {
          if      (mom[3]>0) 
	    remnantparts1.push_back(std::make_pair(momentum,
						   parton->Flav().HepEvt()));
          else if (mom[3]<0) 
	    remnantparts2.push_back(std::make_pair(momentum,
						   parton->Flav().HepEvt()));
          else THROW(fatal_error,"Ill defined beam remnants.");
        }
      }
    }
    
    if ((*blit)->Type()==ATOOLS::btp::Signal_Process) {
      Blob_Data_Base *pdfinfodb((**blit)["PDFInfo"]);
      if (pdfinfodb && (*blit)->NInP()==2) {
        PHASIC::PDF_Info pi=pdfinfodb->Get<PHASIC::PDF_Info>();
        double q(sqrt(sqrt(pi.m_muf12*pi.m_muf22)));
        HepMC::PdfInfo pdfinfo(pi.m_fl1,pi.m_fl2,pi.m_x1,pi.m_x2,
                               q,pi.m_xf1,pi.m_xf2);
        event.set_pdf_info(pdfinfo);
      }
    }
  }
  event.add_vertex(vertex);
  if (beamparticles.size()==2) {
    event.set_beam_particles(beamparticles[0],beamparticles[1]);
  }
  std::vector<double> weights; weights.push_back(weight);
  if (sp && !subevtlist) {
    Blob_Data_Base *info;
    info=((*sp)["MEWeight"]);
    if (!info) THROW(fatal_error,"Missing weight info.");
    double meweight(info->Get<double>());
    weights.push_back(meweight);
    info=((*sp)["Weight_Norm"]);
    if (!info) THROW(fatal_error,"Missing weight normalisation.");
    double weightnorm(info->Get<double>());
    weights.push_back(weightnorm);
    info=(*sp)["Trials"];
    if (!info) THROW(fatal_error,"Missing nof trials.");
    double trials(info->Get<double>());
    weights.push_back(trials);
    //alphaS value && power
    double rscale2 = (*sp)["Renormalization_Scale"]->Get<double>();
    double alphaS = MODEL::s_model->ScalarFunction("alpha_S",rscale2);
    event.set_alphaQCD(alphaS);
    double aSpower = ((*sp)["OQCD"])->Get<int>();
    weights.push_back(aSpower);
    
    if (wgtinfo) {
      //dump weight_0
      weights.push_back(wgtinfo->m_w0);
      //dump number of usr weights
      weights.push_back(wgtinfo->m_nx);
      //store xprimes
      weights.push_back(wgtinfo->m_y1);
      weights.push_back(wgtinfo->m_y2);
      //fill in usr weights
      for (int i=0;i<wgtinfo->m_nx;i++) {
	weights.push_back(wgtinfo->p_wx[i]);
      }
    }
  }
  
  event.weights()=weights;

  if (subevtlist) {
    // build GenEvent for all subevts (where only the signal is available)
    // use stored beam & remnant particles from above
    // sub->m_flip==1 -> momenta need to be inverted for analyses
    for (size_t i(0);i<subevtlist->size();++i) {
      NLO_subevt * sub((*subevtlist)[i]);
      if (sub->m_result==0.) continue;
      HepMC::GenVertex * subvertex(new HepMC::GenVertex());
      HepMC::GenEvent * subevent(new HepMC::GenEvent());
      // assume that only 2->(n-2) processes
      HepMC::GenParticle *beam[2];
      if (beamparts.size()==2) {
        for (size_t j(0);j<beamparts.size();++j) {
          beam[j] = new HepMC::GenParticle(beamparts[j].first,
                                           beamparts[j].second,2);
          subvertex->add_particle_in(beam[j]);
        }
      }
      else {
        const ATOOLS::Vec4D *mom(sub->p_mom);
        const ATOOLS::Flavour *fl(sub->p_fl);
        for (size_t j(0);j<2;++j) {
          HepMC::FourVector momentum;
          momentum.set( mom[j][1], mom[j][2], mom[j][3],mom[j][0]);
          ATOOLS::Flavour flc(fl[j]);
          HepMC::GenParticle* inpart
              = new HepMC::GenParticle(momentum,flc.HepEvt(),2);
          subvertex->add_particle_in(inpart);
        }
      }
      const ATOOLS::Vec4D *mom(sub->p_mom);
      const ATOOLS::Flavour *fl(sub->p_fl);
      for (size_t j(2);j<(*subevtlist)[i]->m_n;++j) {
        HepMC::FourVector momentum;
	momentum.set( mom[j][1], mom[j][2], mom[j][3],mom[j][0]);
        ATOOLS::Flavour flc(fl[j]);
        HepMC::GenParticle* outpart
            = new HepMC::GenParticle(momentum,flc.HepEvt(),1);
        subvertex->add_particle_out(outpart);
      }
      // if beamremnants are present:
      // scale beam remnants of real event for energy momentum conservation :
      //   flavours might not add up properly for sub events,
      //   but who cares. they go down the beam pipe.
      // also assume they are all massless :
      //   this will give momentum conservation violations
      //   which are collider dependent only
      //
      //   for real event (as constructed in BRH) :
      //
      //     P = p + \sum r_i  with  P^2 = m^2  and  r_i^2 = p^2 = 0
      //
      //     further  P  = ( P^0 , 0 , 0 , P^z )
      //              p  = (  p  , 0 , 0 ,  p  )
      //             r_i = ( r_i , 0 , 0 , r_i )
      //
      //     => P^0 = p + \sum r_i
      //        P^z = \sqrt{(P^0)^2 - m^2} <  p + \sum r_i
      //
      // in a mass-symmetric collider, the excess is the same for both beams
      // ensuring momentum conservation
      //
      //   for sub event (constructed here):
      //
      //     P = p~ + \sum r_i~ = p~ + \sum u*r_i
      //
      //     where u = ( P^0 - p^0~ ) / \sum r_i^0
      //
      //     again, the r_i~ = u*r_i are constructed such that
      //
      //     => P^0 = p~ + \sum u*r_i
      //        P^z = \sqrt{(P^0)^2 - m^2} <  p~ + \sum u*r_i =  p + \sum r_i
      //
      //     leading to the same momentum conservation violations per beam
      //
      if (remnantparts1.size()!=0 && remnantparts2.size()!=0) {
        double res1(0.),res2(0.);
        for (size_t j(0);j<remnantparts1.size();++j) {
          res1+=remnantparts1[j].first.e();
        }
        for (size_t j(0);j<remnantparts2.size();++j) {
          res2+=remnantparts2[j].first.e();
        }
        ATOOLS::Vec4D hardparton[2];
        for (size_t j(0);j<2;++j) {
	  hardparton[j]=ATOOLS::Vec4D(mom[j][0],Vec3D( mom[j]));
        }
        // incoming partons might need to be flipped due to particle sorting
        bool flip(hardparton[0][3]<0);
        double u1((beamparts[0].first.e()-hardparton[flip?1:0][0])/res1);
        double u2((beamparts[1].first.e()-hardparton[flip?0:1][0])/res2);
        // filling
        for (size_t j(0);j<remnantparts1.size();++j) {
          HepMC::FourVector momentum(u1*remnantparts1[j].first.px(),
                                     u1*remnantparts1[j].first.py(),
                                     u1*remnantparts1[j].first.pz(),
                                     u1*remnantparts1[j].first.e());
          HepMC::GenParticle* outpart
              = new HepMC::GenParticle(momentum,remnantparts1[j].second,1);
          subvertex->add_particle_out(outpart);
        }
        for (size_t j(0);j<remnantparts2.size();++j) {
          HepMC::FourVector momentum(u2*remnantparts2[j].first.px(),
                                     u2*remnantparts2[j].first.py(),
                                     u2*remnantparts2[j].first.pz(),
                                     u2*remnantparts2[j].first.e());
          HepMC::GenParticle* outpart
              = new HepMC::GenParticle(momentum,remnantparts2[j].second,1);
          subvertex->add_particle_out(outpart);
        }
        if (beamparticles.size()==2) {
          subevent->set_beam_particles(beam[0],beam[1]);
        }
      }
      subevent->add_vertex(subvertex);
      // not enough info in subevents to set PDFInfo properly,
      // so set only flavours, x1, x2, and q from the Signal_Process
      HepMC::PdfInfo subpdfinfo(*event.pdf_info());
      double q = sqrt(sub->m_mu2[stp::fac]);
      if (q!=0. && q != subpdfinfo.scalePDF()) 
        subpdfinfo.set_scalePDF(q);
      if (sub->m_xf1) subpdfinfo.set_pdf1(sub->m_xf1);
      if (sub->m_xf2) subpdfinfo.set_pdf2(sub->m_xf2);
      subevent->set_pdf_info(subpdfinfo);
      // add weight
      std::vector<double> subweight; subweight.push_back(sub->m_result);
      Blob_Data_Base *info;
      info=((*sp)["MEWeight"]);
      if (!info) THROW(fatal_error,"Missing weight info.");
      double meweight(info->Get<double>());
      subweight.push_back(meweight);
      info=((*sp)["Weight_Norm"]);
      if (!info) THROW(fatal_error,"Missing weight normalisation.");
      double weightnorm(info->Get<double>());
      subweight.push_back(weightnorm);
      info=(*sp)["Trials"];
      if (!info) THROW(fatal_error,"Missing nof trials.");
      double trials(info->Get<double>());
      subweight.push_back(trials);
      //alphaS value && power
      double alphaS = MODEL::s_model->ScalarFunction("alpha_S",sub->m_mu2[stp::ren]);
      subevent->set_alphaQCD(alphaS);
      double aSpower = ((*sp)["OQCD"])->Get<int>();
      subweight.push_back(aSpower);
      subweight.push_back(sub->m_mewgt);
      
      if (wgtinfo) {
	//dump number of usr weights
	subweight.push_back(wgtinfo->m_nx);
	//store xprimes
	subweight.push_back(wgtinfo->m_y1);
	subweight.push_back(wgtinfo->m_y2);
	//fill in usr weights
	for (int i=0;i<wgtinfo->m_nx;i++) {
	  subweight.push_back(wgtinfo->p_wx[i]);
	}
      }
      //
      subevent->weights()=subweight;
      // set the event number (could be used to identify correlated events)
      subevent->set_event_number(ATOOLS::rpa->gen.NumberOfGeneratedEvents());
      m_subeventlist.push_back(subevent);
    }
  }
  return true;
}

bool HepMC2_Interface::Sherpa2HepMC(ATOOLS::Blob_List *const blobs,
				    double weight)
{
  if (blobs->empty()) {
    msg_Error()<<"Error in "<<METHOD<<"."<<std::endl
               <<"   Empty list - nothing to translate into HepMC."<<std::endl
               <<"   Continue run ... ."<<std::endl;
    return true;
  }
  if (p_event!=NULL) delete p_event;
  for (size_t i=0; i<m_subeventlist.size();++i)
    delete m_subeventlist[i];
  m_subeventlist.clear();
  p_event = new HepMC::GenEvent();
  return Sherpa2HepMC(blobs, *p_event, weight);
}

bool HepMC2_Interface::Sherpa2ShortHepMC(ATOOLS::Blob_List *const blobs,
                                         double weight)
{
  if (blobs->empty()) {
    msg_Error()<<"Error in "<<METHOD<<"."<<std::endl
               <<"   Empty list - nothing to translate into HepMC."<<std::endl
               <<"   Continue run ... ."<<std::endl;
    return true;
  }
  if (p_event!=NULL) delete p_event;
  DeleteGenSubEventList();
  p_event = new HepMC::GenEvent();
  return Sherpa2ShortHepMC(blobs, *p_event, weight);
}

bool HepMC2_Interface::Sherpa2HepMC(ATOOLS::Blob * blob, 
				    HepMC::GenVertex *& vertex,
				    const std::map<size_t,size_t> &decids)
{
  if (m_ignoreblobs.count(blob->Type())) return false;
  int count = m_blob2genvertex.count(blob);
  if (count>0) {
    vertex = m_blob2genvertex[blob];
    return true;
  }
  else {
    ATOOLS::Vec4D pos = blob->Position();
    HepMC::FourVector position(pos[1],pos[2],pos[3],pos[0]);
    vertex = new HepMC::GenVertex(position,blob->Id());
    vertex->weights().push_back(1.);
    if (blob->Type()==btp::Signal_Process)      vertex->set_id(1);
    else if (blob->Type()==btp::Hard_Collision) vertex->set_id(2);
    else if (blob->Type()==btp::Hard_Decay)     vertex->set_id(3);
    else if (blob->Type()==btp::Shower || 
	     blob->Type()==btp::QED_Radiation)  vertex->set_id(4);
    else if (blob->Type()==btp::Fragmentation)  vertex->set_id(5);
    else if (blob->Type()==btp::Hadron_Decay)   vertex->set_id(6);
      //{  
      //if ((*blob)["Partonic"]!=NULL) vertex->set_id(-6);
      //else vertex->set_id(6);
      //}
    else vertex->set_id(0);
  }

  bool okay = 1;
  HepMC::GenParticle * _particle;
  for (int i=0;i<blob->NInP();i++) {
    if (Sherpa2HepMC(blob->InParticle(i),_particle,decids)) {
      vertex->add_particle_in(_particle);
    }
    else okay = 0;
  }
  for (int i=0;i<blob->NOutP();i++) {
    if (Sherpa2HepMC(blob->OutParticle(i),_particle,decids)) {
      vertex->add_particle_out(_particle);
    }
    else okay = 0;
  }
  m_blob2genvertex.insert(std::make_pair(blob,vertex));
  if (!okay) {
    msg_Error()<<"Error in HepMC2_Interface::Sherpa2HepMC(Blob,Vertex)."<<std::endl
		       <<"   Continue event generation with new event."<<std::endl;
  }
  if (msg_LevelIsDebugging()) {
    ATOOLS::Vec4D check = blob->CheckMomentumConservation();
    double test         = ATOOLS::Vec3D(check).Abs();
    if (ATOOLS::dabs(1.-vertex->check_momentum_conservation()/test)>1.e-5 &&
	ATOOLS::dabs(test)>1.e-5)
      {
	msg_Error()<<"ERROR in "<<METHOD<<std::endl
			   <<"   Momentum not conserved. Continue."<<std::endl
			   <<"ERROR in Blob -> Vertex : "<<vertex->check_momentum_conservation()
			   <<" <- "<<test<<" "<<check
			   <<std::endl<<(*blob)<<std::endl;
	vertex->print(std::cout);
	std::cout<<"--------------------------------------------------------"<<std::endl;
      }
  }
  return okay;
}

bool HepMC2_Interface::Sherpa2HepMC(ATOOLS::Particle * parton,HepMC::GenParticle *& particle,
				    const std::map<size_t,size_t> &decids)
{
  int count = m_particle2genparticle.count(parton);
  if (count>0) {
    particle = m_particle2genparticle[parton];
    return true;
  }

  ATOOLS::Vec4D mom  = parton->Momentum();
  HepMC::FourVector momentum(mom[1],mom[2],mom[3],mom[0]);
  int status=11;
  if (parton->DecayBlob()==NULL ||
      m_ignoreblobs.count(parton->DecayBlob()->Type())!=0) {
    status=1;
  }
  else {
    if (parton->DecayBlob()->Type()==ATOOLS::btp::Hadron_Decay ||
        parton->DecayBlob()->Type()==ATOOLS::btp::Hadron_Mixing) {
      status=2;
    }
    else if (parton->DecayBlob()->Type()==ATOOLS::btp::Signal_Process ||
             (parton->ProductionBlob() &&
              parton->ProductionBlob()->Type()==ATOOLS::btp::Signal_Process)) {
      status=3;
    }
    else if (parton->DecayBlob()->Type()==ATOOLS::btp::Bunch) {
      status=4;
    }
    else {
      status=11;
      for (std::map<size_t,size_t>::const_iterator did(decids.begin());
	   did!=decids.end();++did) if (did->first&parton->MEId()) status=did->second;
    }
  }
  particle = new HepMC::GenParticle(momentum,parton->Flav().HepEvt(),status);
  for (int i=1;i<3;i++) {
    if (parton->GetFlow(i)>0) particle->set_flow(i,parton->GetFlow(i));
  }
  m_particle2genparticle.insert(std::make_pair(parton,particle));
  return true;
}

void HepMC2_Interface::DeleteGenSubEventList()
{
  for (size_t i=0; i<m_subeventlist.size();++i)
    delete m_subeventlist[i];
  m_subeventlist.clear();
}

#endif
