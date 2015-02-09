#include "SHERPA/Tools/HepMC2_Interface.H"
#ifdef USING__HEPMC2

#include "ATOOLS/Phys/Blob_List.H"
#include "ATOOLS/Phys/Particle.H"
#include "ATOOLS/Phys/NLO_Subevt.H"
#include "ATOOLS/Phys/Weight_Info.H"
#include "ATOOLS/Math/Vector.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "MODEL/Main/Model_Base.H"
#include "PHASIC++/Main/Phase_Space_Handler.H"
#include "SHERPA/Tools/Scale_Variations.H"

#include "HepMC/GenEvent.h"
#include "HepMC/GenVertex.h"
#include "HepMC/GenParticle.h"
#include "HepMC/SimpleVector.h"
#include "HepMC/PdfInfo.h"
#include "HepMC/WeightContainer.h"
#ifdef USING__HEPMC2__UNITS
#include "HepMC/Units.h"
#endif

using namespace SHERPA;
using namespace ATOOLS;

EventInfo::EventInfo(ATOOLS::Blob * sp, const double &wgt,
                     bool namedweights) :
  m_usenamedweights(namedweights), p_sp(sp),
  m_oqcd(0), m_oew(0), m_wgt(wgt),
  m_mewgt(0.), m_wgtnorm(0.), m_ntrials(0.),
  m_mur2(0.), m_muf12(0.), m_muf22(0.),
  m_alphas(0.), m_alpha(0.), m_type(PHASIC::nlo_type::lo),
  p_wgtinfo(NULL), p_pdfinfo(NULL), p_subevtlist(NULL),
  p_nsvmap(NULL)
{
  if (p_sp) {
    DEBUG_FUNC(*p_sp);
    Blob_Data_Base *db;
    ReadIn(db,"MEWeight",true);
    m_mewgt=db->Get<double>();
    ReadIn(db,"Weight_Norm",true);
    m_wgtnorm=db->Get<double>();
    ReadIn(db,"Trials",true);
    m_ntrials=db->Get<double>();
    ReadIn(db,"Renormalization_Scale",true);
    m_mur2=db->Get<double>();
    SetAlphaS();
    SetAlpha();
    ReadIn(db,"OQCD",true);
    m_oqcd=db->Get<int>();
    ReadIn(db,"OEW",true);
    m_oew=db->Get<int>();
    ReadIn(db,"MEWeightInfo",true);
    p_wgtinfo=db->Get<ME_Weight_Info*>();
    ReadIn(db,"PDFInfo",true);
    p_pdfinfo=&db->Get<ATOOLS::PDF_Info>();
    m_muf12=p_pdfinfo->m_muf12;
    m_muf22=p_pdfinfo->m_muf22;
    ReadIn(db,"NLO_subeventlist",false);
    if (db) p_subevtlist=db->Get<NLO_subevtlist*>();
    if (p_subevtlist) m_type=p_subevtlist->Type();
    ReadIn(db,"ScaleVariations",false);
    if (db) p_nsvmap=db->Get<NamedScaleVariationMap*>();
  }
}

EventInfo::EventInfo(const EventInfo &evtinfo) :
  m_usenamedweights(evtinfo.m_usenamedweights),
  p_sp(evtinfo.p_sp),
  m_oqcd(evtinfo.m_oqcd), m_oew(evtinfo.m_oew),
  m_wgt(0.), m_mewgt(0.), m_wgtnorm(0.),
  m_ntrials(evtinfo.m_ntrials),
  m_mur2(0.), m_muf12(0.), m_muf22(0.),
  m_alphas(0.), m_alpha(0.), m_type(evtinfo.m_type),
  p_wgtinfo(NULL), p_pdfinfo(evtinfo.p_pdfinfo),
  p_subevtlist(evtinfo.p_subevtlist), p_nsvmap(evtinfo.p_nsvmap)
{
}

bool EventInfo::ReadIn(ATOOLS::Blob_Data_Base* &db,std::string name,bool abort)
{
  db=(*p_sp)[name];
  if (abort && !db) THROW(fatal_error,name+" information missing.");
}

bool EventInfo::WriteTo(HepMC::GenEvent &evt, const int& idx)
{
  DEBUG_FUNC("use named weights: "<<m_usenamedweights);
  HepMC::WeightContainer wc;
  if (m_usenamedweights) {
#ifdef HEPMC_HAS_NAMED_WEIGHTS
    // fill standard entries to ensure backwards compatability
    wc["Weight"]=m_wgt;
    wc["MEWeight"]=m_mewgt;
    wc["WeightNormalisation"]=m_wgtnorm;
    wc["NTrials"]=m_ntrials;
    // additional entries for LO/LOPS reweighting
    // x1,x2,muf2 can be found in PdfInfo; alphaS,alphaQED in their infos
    wc["MuR2"]=m_mur2;
    wc["OQCD"]=m_oqcd;
    wc["OEW"]=m_oew;
    // fill scale variations map into weight container
    msg_Debugging()<<"#named wgts: "<<p_nsvmap->size()<<std::endl;
    if (p_nsvmap) {
      for (NamedScaleVariationMap::const_iterator it(p_nsvmap->begin());
           it!=p_nsvmap->end();++it) {
        if (idx==-1) wc[it->first]=it->second->Value();
        else         wc[it->first]=it->second->Value(idx);
      }
    }
    size_t nt(0);
    if (p_wgtinfo && p_wgtinfo->m_type!=mewgttype::none) {
      wc["Reweight_B"]=p_wgtinfo->m_B;
      if (p_wgtinfo->m_B) nt|=1;
      wc["Reweight_VI"]=p_wgtinfo->m_VI;
      if (p_wgtinfo->m_VI) nt|=2;
      wc["Reweight_KP"]=p_wgtinfo->m_KP;
      if (p_wgtinfo->m_KP) nt|=4;
      wc["Reweight_KP_x1p"]=p_wgtinfo->m_y1;
      wc["Reweight_KP_x2p"]=p_wgtinfo->m_y2;
      wc["Reweight_MuR2"]=p_wgtinfo->m_mur2;
      wc["Reweight_MuF2"]=p_wgtinfo->m_muf2;
      if (p_wgtinfo->m_type&mewgttype::muR) {
        for (size_t i=0;i<p_wgtinfo->m_wren.size();++i) {
          wc["Reweight_VI_wren_"+ToString(i)]=p_wgtinfo->m_wren[i];
        }
      }
      if (p_wgtinfo->m_type&mewgttype::muF) {
        for (size_t i=0;i<p_wgtinfo->m_wfac.size();++i) {
          wc["Reweight_KP_wfac_"+ToString(i)]=p_wgtinfo->m_wfac[i];
        }
      }
      if (p_wgtinfo->m_dadsinfos.size()) {
        wc.push_back((p_wgtinfo->m_dadsinfos.size()));
        for (size_t i(0);i<p_wgtinfo->m_dadsinfos.size();++i) {
          wc["Reweight_DADS_"+ToString(i)+"_Weight"]
              =p_wgtinfo->m_dadsinfos[i].m_wgt;
          wc["Reweight_DADS_"+ToString(i)+"_x1"]
              =p_wgtinfo->m_dadsinfos[i].m_pdf.m_x1;
          wc["Reweight_DADS_"+ToString(i)+"_x2"]
              =p_wgtinfo->m_dadsinfos[i].m_pdf.m_x2;
          wc["Reweight_DADS_"+ToString(i)+"_fl1"]
              =p_wgtinfo->m_dadsinfos[i].m_pdf.m_fl1;
          wc["Reweight_DADS_"+ToString(i)+"_fl2"]
              =p_wgtinfo->m_dadsinfos[i].m_pdf.m_fl2;
          wc["Reweight_DADS_"+ToString(i)+"_MuR2"]
              =p_wgtinfo->m_dadsinfos[i].m_mur2;
          wc["Reweight_DADS_"+ToString(i)+"_MuF12"]
              =p_wgtinfo->m_dadsinfos[i].m_pdf.m_muf12;
          wc["Reweight_DADS_"+ToString(i)+"_MuF22"]
              =p_wgtinfo->m_dadsinfos[i].m_pdf.m_muf22;
        }
      }
    }
    if (p_subevtlist) {
      nt|=32;
      wc["Reweight_RS"]=m_mewgt;
    }
    wc["NLO_Type"]=nt;
#else
    THROW(fatal_error,"Asked for named weights, but HepMC version too old.");
#endif
  }
  else {
    wc.push_back(m_wgt);
    wc.push_back(m_mewgt);
    wc.push_back(m_wgtnorm);
    wc.push_back(m_ntrials);
    wc.push_back(m_mur2);
    wc.push_back(m_muf12);
    wc.push_back(m_muf22);
    wc.push_back(m_oqcd);
    wc.push_back(m_oew);
    wc.push_back(p_subevtlist?32:0);
    if (p_wgtinfo) {
      //dump weight_0
      wc.push_back(p_wgtinfo->m_B);
      wc.push_back(p_wgtinfo->m_VI);
      wc.push_back(p_wgtinfo->m_KP);
      wc.push_back(p_wgtinfo->m_RS);
      //dump number of usr weights
      size_t nentries(0);
      if (p_wgtinfo->m_type&mewgttype::muR) nentries+=2;
      if (p_wgtinfo->m_type&mewgttype::muF) nentries+=16;
      wc.push_back(nentries);
      //store xprimes
      wc.push_back(p_wgtinfo->m_y1);
      wc.push_back(p_wgtinfo->m_y2);
      //fill in usr weights
      if (p_wgtinfo->m_type&mewgttype::muR) {
        for (size_t i=0;i<p_wgtinfo->m_wren.size();++i) {
          wc.push_back(p_wgtinfo->m_wren[i]);
        }
      }
      if (p_wgtinfo->m_type&mewgttype::muF) {
        for (size_t i=0;i<p_wgtinfo->m_wfac.size();++i) {
          wc.push_back(p_wgtinfo->m_wfac[i]);
        }
      }
      if (p_wgtinfo->m_dadsinfos.size()) {
        wc.push_back((p_wgtinfo->m_dadsinfos.size()));
        for (size_t i(0);i<p_wgtinfo->m_dadsinfos.size();++i) {
          wc.push_back(p_wgtinfo->m_dadsinfos[i].m_wgt);
          wc.push_back(p_wgtinfo->m_dadsinfos[i].m_pdf.m_x1);
          wc.push_back(p_wgtinfo->m_dadsinfos[i].m_pdf.m_x2);
          wc.push_back(p_wgtinfo->m_dadsinfos[i].m_pdf.m_fl1);
          wc.push_back(p_wgtinfo->m_dadsinfos[i].m_pdf.m_fl2);
          wc.push_back(p_wgtinfo->m_dadsinfos[i].m_pdf.m_muf12);
          wc.push_back(p_wgtinfo->m_dadsinfos[i].m_pdf.m_muf22);
        }
      }
    }
  }
  evt.weights()=wc;
  if (p_pdfinfo) {
    double q(sqrt(sqrt(p_pdfinfo->m_muf12*p_pdfinfo->m_muf22)));
    HepMC::PdfInfo pdfinfo(p_pdfinfo->m_fl1,p_pdfinfo->m_fl2,
                           p_pdfinfo->m_x1,p_pdfinfo->m_x2,
                           q,p_pdfinfo->m_xf1,p_pdfinfo->m_xf2);
    evt.set_pdf_info(pdfinfo);
  }
  evt.set_alphaQCD(m_alphas);
  evt.set_alphaQED(m_alpha);
  return true;
}

void EventInfo::SetAlphaS()
{
  m_alphas=MODEL::s_model->ScalarFunction("alpha_S",m_mur2);
}

void EventInfo::SetAlpha()
{
  m_alpha=MODEL::s_model->ScalarFunction("alpha_QED");
}

HepMC2_Interface::HepMC2_Interface() :
  m_usenamedweights(false), p_event(NULL)
{
  Data_Reader reader(" ",";","!","=");
  reader.AddComment("#");
  reader.AddWordSeparator("\t");
#ifdef HEPMC_HAS_NAMED_WEIGHTS
  m_usenamedweights=reader.GetValue<int>("HEPMC_USE_NAMED_WEIGHTS",true);
#endif
}

HepMC2_Interface::~HepMC2_Interface()
{
  delete p_event;
  DeleteGenSubEventList();
}

bool HepMC2_Interface::Sherpa2HepMC(ATOOLS::Blob_List *const blobs,
                                    HepMC::GenEvent& event, double weight)
{
  DEBUG_FUNC("");
#ifdef USING__HEPMC2__UNITS
  event.use_units(HepMC::Units::GEV,
                  HepMC::Units::MM);
#endif
  event.set_event_number(ATOOLS::rpa->gen.NumberOfGeneratedEvents());
  size_t decid(11);
  std::map<size_t,size_t> decids;
  Blob *sp(blobs->FindFirst(btp::Signal_Process));
  EventInfo evtinfo(sp,weight,m_usenamedweights);
  evtinfo.WriteTo(event);
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
  return true;
}

bool HepMC2_Interface::Sherpa2ShortHepMC(ATOOLS::Blob_List *const blobs,
                                         HepMC::GenEvent& event, double weight)
{
#ifdef USING__HEPMC2__UNITS
  event.use_units(HepMC::Units::GEV,
                  HepMC::Units::MM);
#endif
  Blob *sp(blobs->FindFirst(btp::Signal_Process));
  EventInfo evtinfo(sp,weight,m_usenamedweights);
  // when subevtlist, fill hepmc-subevtlist
  if (evtinfo.SubEvtList()) return SubEvtList2ShortHepMC(evtinfo);
  event.set_event_number(ATOOLS::rpa->gen.NumberOfGeneratedEvents());
  evtinfo.WriteTo(event);
  HepMC::GenVertex * vertex=new HepMC::GenVertex();
  std::vector<HepMC::GenParticle*> beamparticles;
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
      }
    }
  }
  event.add_vertex(vertex);
  if (beamparticles.size()==2) {
    event.set_beam_particles(beamparticles[0],beamparticles[1]);
  }

  return true;
}

bool HepMC2_Interface::SubEvtList2ShortHepMC(EventInfo &evtinfo)
{
  DEBUG_FUNC("subevts: "<<evtinfo.SubEvtList()->size());
  // build GenEvent for all subevts (where only the signal is available)
  // purely partonic, no beam information, may add, if needed
  for (size_t i(0);i<evtinfo.SubEvtList()->size();++i) {
    EventInfo subevtinfo(evtinfo);
    const NLO_subevt * sub((*evtinfo.SubEvtList())[i]);
    if (sub->m_result==0.) continue;
    HepMC::GenVertex * subvertex(new HepMC::GenVertex());
    HepMC::GenEvent * subevent(new HepMC::GenEvent());
    // set the event number (could be used to identify correlated events)
    subevent->set_event_number(ATOOLS::rpa->gen.NumberOfGeneratedEvents());
    // assume that only 2->(n-2) processes
    for (size_t j(0);j<2;++j) {
      HepMC::FourVector momentum(sub->p_mom[j][1],sub->p_mom[j][2],
                                 sub->p_mom[j][3],sub->p_mom[j][0]);
      HepMC::GenParticle* inpart =
        new HepMC::GenParticle(momentum,sub->p_fl[j].HepEvt(),2);
      subvertex->add_particle_in(inpart);
    }
    for (size_t j(2);j<sub->m_n;++j) {
      HepMC::FourVector momentum(sub->p_mom[j][1],sub->p_mom[j][2],
                                 sub->p_mom[j][3],sub->p_mom[j][0]);
      HepMC::GenParticle* outpart =
        new HepMC::GenParticle(momentum,sub->p_fl[j].HepEvt(),1);
      subvertex->add_particle_out(outpart);
    }
    subevent->add_vertex(subvertex);
    // not enough info in subevents to set PDFInfo properly,
    // so set flavours and x1, x2 from the Signal_Process
    // reset muR, muF, alphaS, alpha
    subevtinfo.SetWeight(sub->m_result);
    subevtinfo.SetMEWeight(sub->m_mewgt);
    subevtinfo.SetMuR2(sub->m_mu2[stp::ren]);
    subevtinfo.SetMuF12(sub->m_mu2[stp::fac]);
    subevtinfo.SetMuF22(sub->m_mu2[stp::fac]);
    subevtinfo.SetAlphaS();
    subevtinfo.SetAlpha();
    subevtinfo.WriteTo(*subevent,i);
    if (msg_LevelIsDebugging()) subevent->print(msg_Out());
    m_subeventlist.push_back(subevent);
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
