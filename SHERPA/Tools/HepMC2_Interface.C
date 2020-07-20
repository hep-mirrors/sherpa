#include "SHERPA/Tools/HepMC2_Interface.H"
#ifdef USING__HEPMC2

#include "ATOOLS/Phys/Blob_List.H"
#include "ATOOLS/Phys/Particle.H"
#include "ATOOLS/Phys/NLO_Subevt.H"
#include "ATOOLS/Phys/Weight_Info.H"
#include "ATOOLS/Math/Vector.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Scoped_Settings.H"
#include "MODEL/Main/Model_Base.H"
#include "PHASIC++/Main/Phase_Space_Handler.H"

#include "HepMC/GenEvent.h"
#include "HepMC/GenVertex.h"
#include "HepMC/GenParticle.h"
#include "HepMC/SimpleVector.h"
#include "HepMC/PdfInfo.h"
#include "HepMC/WeightContainer.h"
#ifdef USING__HEPMC2__DEFS
#include "HepMC/HepMCDefs.h"
#ifdef HEPMC_HAS_CROSS_SECTION
#include "HepMC/GenCrossSection.h"
#endif
#endif
#ifdef USING__HEPMC2__UNITS
#include "HepMC/Units.h"
#endif

using namespace SHERPA;
using namespace ATOOLS;

EventInfo::EventInfo(ATOOLS::Blob * sp, const double &wgt,
                     bool namedweights,
                     bool extendedweights,
                     bool includemeonlyweights) :
  p_sp(sp),
  m_wgt(wgt),
  m_usenamedweights(namedweights),
  m_extendedweights(extendedweights),
  m_variationsources(1, ATOOLS::Variations_Source::all),
  m_mewgt(0.), m_wgtnorm(wgt), m_ntrials(1.),
  m_pswgt(0.), m_pwgt(0.),  m_userhook(false), m_userweight(0.),
  m_mur2(0.), m_muf12(0.), m_muf22(0.),
  m_alphas(0.), m_alpha(0.), m_type(ATOOLS::nlo_type::lo),
  p_wgtinfo(NULL), p_pdfinfo(NULL), p_subevtlist(NULL)
{
  if (p_sp) {
    DEBUG_FUNC(*p_sp);
    Blob_Data_Base *db;
    ReadIn(db,"MEWeight",false);
    if (db) m_mewgt=db->Get<double>();
    m_pswgt=m_wgt/m_mewgt;
    ReadIn(db,"Weight_Norm",true);
    m_wgtnorm=db->Get<double>();
    ReadIn(db,"Trials",true);
    m_ntrials=db->Get<double>();
    ReadIn(db,"PDFInfo",false);
    if (db) {
      p_pdfinfo=&db->Get<ATOOLS::PDF_Info>();
      m_muf12=p_pdfinfo->m_muf12;
      m_muf22=p_pdfinfo->m_muf22;
    }
    ReadIn(db,"UserHook",false);
    if (db) {
      m_userhook=true;
      m_userweight=db->Get<double>();
    }
    ReadIn(db,"Renormalization_Scale",false);
    if (db) m_mur2=db->Get<double>();
    SetAlphaS();
    SetAlpha();
    if (m_extendedweights) {
      ReadIn(db,"Orders",true);
      m_orders=db->Get<std::vector<double> >();
      ReadIn(db,"MEWeightInfo",true);
      p_wgtinfo=db->Get<ME_Weight_Info*>();
    }
    ReadIn(db,"NLO_subeventlist",false);
    if (db) p_subevtlist=db->Get<NLO_subevtlist*>();
    if (p_subevtlist) m_type=p_subevtlist->Type();

    if (includemeonlyweights)
      m_variationsources.push_back(ATOOLS::Variations_Source::main);

    ReadIn(db, "WeightsMap", false);
    if (db) {
      m_wgtmap = db->Get<Weights_Map>();
      if (m_wgtmap.HasVariations() && !m_usenamedweights)
        THROW(fatal_error,
              "Scale, AlphaS and/or PDF variations cannot be written to " +
                  std::string("HepMC without using named weights. ") +
                  std::string("Try HEPMC_USE_NAMED_WEIGHTS: true"));
    }
  }
}

EventInfo::EventInfo(const EventInfo &evtinfo) :
  m_usenamedweights(evtinfo.m_usenamedweights),
  m_extendedweights(evtinfo.m_extendedweights),
  m_variationsources(evtinfo.m_variationsources),
  p_sp(evtinfo.p_sp),
  m_orders(evtinfo.m_orders),
  m_wgt(0.), m_mewgt(0.), m_wgtnorm(0.),
  m_ntrials(evtinfo.m_ntrials), m_pswgt(evtinfo.m_pswgt), m_pwgt(0.),
  m_mur2(0.), m_muf12(0.), m_muf22(0.),
  m_alphas(0.), m_alpha(0.), m_type(evtinfo.m_type),
  p_wgtinfo(NULL), p_pdfinfo(evtinfo.p_pdfinfo),
  p_subevtlist(evtinfo.p_subevtlist),
  m_wgtmap(evtinfo.m_wgtmap)
{
}

void EventInfo::ReadIn(ATOOLS::Blob_Data_Base* &db,std::string name,bool abort)
{
  db=(*p_sp)[name];
  if (abort && !db) THROW(fatal_error,name+" information missing.");
}

bool EventInfo::WriteTo(HepMC::GenEvent &evt, const int& idx)
{
  DEBUG_FUNC("use named weights: "<<m_usenamedweights
             <<", extended weights: "<<m_extendedweights);
  HepMC::WeightContainer wc;
  if (m_usenamedweights) {
#ifdef HEPMC_HAS_NAMED_WEIGHTS
    // fill standard entries to ensure backwards compatability
    wc["Weight"]=m_wgt;
    wc["MEWeight"]=m_mewgt;
    wc["WeightNormalisation"]=m_wgtnorm;
    wc["NTrials"]=m_ntrials;
    if (m_userhook) wc["UserHook"]=m_userweight;
    if (m_extendedweights) {
      wc["PSWeight"]=m_pswgt;
      // additional entries for LO/LOPS reweighting
      // x1,x2,muf2 can be found in PdfInfo; alphaS,alphaQED in their infos
      wc["MuR2"]=m_mur2;
      wc["OQCD"]=m_orders[0];
      wc["OEW"]=m_orders[1];
      if (p_wgtinfo) {
        wc["Reweight_B"]=p_wgtinfo->m_B;
        wc["Reweight_MuR2"]=p_wgtinfo->m_mur2;
        wc["Reweight_MuF2"]=p_wgtinfo->m_muf2;
        if (p_wgtinfo->m_type&mewgttype::VI) {
          wc["Reweight_VI"]=p_wgtinfo->m_VI;
          for (size_t i=0;i<p_wgtinfo->m_wren.size();++i) {
            wc["Reweight_VI_wren_"+ToString(i)]=p_wgtinfo->m_wren[i];
          }
        }
        if (p_wgtinfo->m_type&mewgttype::KP) {
          wc["Reweight_KP"]=p_wgtinfo->m_KP;
          wc["Reweight_KP_x1p"]=p_wgtinfo->m_y1;
          wc["Reweight_KP_x2p"]=p_wgtinfo->m_y2;
          for (size_t i=0;i<p_wgtinfo->m_wfac.size();++i) {
            wc["Reweight_KP_wfac_"+ToString(i)]=p_wgtinfo->m_wfac[i];
          }
        }
        if (p_wgtinfo->m_type&mewgttype::DADS &&
            p_wgtinfo->m_dadsinfos.size()) {
          wc["Reweight_DADS_N"]=p_wgtinfo->m_dadsinfos.size();
          for (size_t i(0);i<p_wgtinfo->m_dadsinfos.size();++i) {
            wc["Reweight_DADS_"+ToString(i)+"_Weight"]
                =p_wgtinfo->m_dadsinfos[i].m_wgt;
            if (p_wgtinfo->m_dadsinfos[i].m_wgt) {
              wc["Reweight_DADS_"+ToString(i)+"_x1"]
                  =p_wgtinfo->m_dadsinfos[i].m_x1;
              wc["Reweight_DADS_"+ToString(i)+"_x2"]
                  =p_wgtinfo->m_dadsinfos[i].m_x2;
              wc["Reweight_DADS_"+ToString(i)+"_fl1"]
                  =p_wgtinfo->m_dadsinfos[i].m_fl1;
              wc["Reweight_DADS_"+ToString(i)+"_fl2"]
                  =p_wgtinfo->m_dadsinfos[i].m_fl2;
            }
          }
        }
        if (p_wgtinfo->m_type&mewgttype::METS &&
            p_wgtinfo->m_clusseqinfo.m_txfl.size()) {
          wc["Reweight_ClusterStep_N"]=p_wgtinfo->m_clusseqinfo.m_txfl.size();
          for (size_t i(0);i<p_wgtinfo->m_clusseqinfo.m_txfl.size();++i) {
            wc["Reweight_ClusterStep_"+ToString(i)+"_t"]
                =p_wgtinfo->m_clusseqinfo.m_txfl[i].m_t;
            wc["Reweight_ClusterStep_"+ToString(i)+"_x1"]
                =p_wgtinfo->m_clusseqinfo.m_txfl[i].m_xa;
            wc["Reweight_ClusterStep_"+ToString(i)+"_x2"]
                =p_wgtinfo->m_clusseqinfo.m_txfl[i].m_xb;
            wc["Reweight_ClusterStep_"+ToString(i)+"_fl1"]
                =p_wgtinfo->m_clusseqinfo.m_txfl[i].m_fla;
            wc["Reweight_ClusterStep_"+ToString(i)+"_fl2"]
                =p_wgtinfo->m_clusseqinfo.m_txfl[i].m_flb;
          }
        }
        if (p_wgtinfo->m_type&mewgttype::H) {
          wc["Reweight_RDA_N"]=p_wgtinfo->m_rdainfos.size();
          for (size_t i(0);i<p_wgtinfo->m_rdainfos.size();++i) {
            wc["Reweight_RDA_"+ToString(i)+"_Weight"]
                =p_wgtinfo->m_rdainfos[i].m_wgt;
            if (p_wgtinfo->m_rdainfos[i].m_wgt) {
              const double mur2(p_wgtinfo->m_rdainfos[i].m_mur2);
              wc["Reweight_RDA_"+ToString(i)+"_MuR2"] = mur2;
              wc["Reweight_RDA_"+ToString(i)+"_MuF12"]
                  =p_wgtinfo->m_rdainfos[i].m_muf12;
              wc["Reweight_RDA_"+ToString(i)+"_MuF22"]
                  =p_wgtinfo->m_rdainfos[i].m_muf22;
              wc["Reweight_RDA_"+ToString(i)+"_Dipole"]
                  =10000*p_wgtinfo->m_rdainfos[i].m_i
                    +100*p_wgtinfo->m_rdainfos[i].m_j
                        +p_wgtinfo->m_rdainfos[i].m_k;
              wc["Reweight_RDA_"+ToString(i)+"_AlphaS"]
                  =MODEL::s_model->ScalarFunction("alpha_S", mur2);
            }
          }
        }
        wc["Reweight_Type"]=p_wgtinfo->m_type;
        // fill EW Sudakov K factor contributions into weight container
        // TODO: port to HepMC3 interface
        for (const auto& kv : p_wgtinfo->m_ewsudakovkfacdelta) {
          MyStrStream s;
          s << "EWSud_" << kv.first;
          wc[s.str()] = kv.second;
        }
      }
      if (p_subevtlist) {
        wc["Reweight_RS"]=m_pwgt;
        wc["Reweight_Type"]=64|(p_wgtinfo?p_wgtinfo->m_type:0);
      }
    }
    else {
      // if using minimal weights still dump event type if RS need correls
      if (p_subevtlist) wc["Reweight_Type"]=64;
    }
    // fill weight variations into weight container
    if (m_wgtmap.HasVariations()) {

      // output subevent weights instead if this is a subevent
      Weights_Map& wgtmap =
          (idx == -1) ? m_wgtmap : (*p_subevtlist)[idx]->m_results;

      for (const auto& source : m_variationsources) {

        for (const auto type : s_variations->ManagedVariationTypes()) {

          // calculate contributions
          Weights weights = Weights {type};
          double relfac {1.0};
          if (source == ATOOLS::Variations_Source::all) {
            weights *= wgtmap.Combine(type);
            relfac = wgtmap.NominalIgnoringVariationType(type);
          } else {
            // calculate nominal, relfac and weights ignoring shower weights
            std::unordered_set<std::string> shower_keys {
                "PS", "PS_QCUT", "MC@NLO_PS", "MC@NLO_QCUT"};
            for (const auto& v : wgtmap) {
              if (shower_keys.find(v.first) != shower_keys.end())
                continue;
              if (v.second.Type() == type) {
                weights *= v.second;
              } else {
                relfac *= v.second.Nominal();
              }
            }
            relfac *= wgtmap.BaseWeight();
          }

          // do remaining combination and output resulting weights
          size_t num_vars = weights.Size() - 1;
          for (size_t i(0); i < num_vars; ++i) {
            const std::string varname {weights.Name(i + 1)};
            const std::string typevarname {
                (source == ATOOLS::Variations_Source::main)
                    ? "ME_ONLY_" + varname
                    : varname};
            wc[typevarname] = weights.Variation(i) * relfac;
            msg_Debugging() << typevarname << " (" << typevarname
                            << "): " << weights.Variation(i) * relfac << '\n';
          }
        }
      }
    }
#else
    THROW(fatal_error,"Asked for named weights, but HepMC version too old.");
#endif
  }
  else {
    // only offer basic event record for unnamed weights
    wc.push_back(m_wgt);
    wc.push_back(m_mewgt);
    wc.push_back(m_wgtnorm);
    wc.push_back(m_ntrials);
    if (m_extendedweights) {
      wc.push_back(m_pswgt);
      wc.push_back(m_mur2);
      wc.push_back(m_muf12);
      wc.push_back(m_muf22);
      wc.push_back(m_orders[0]);
      wc.push_back(m_orders[1]);
    }
    wc.push_back(p_subevtlist?64:0);
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
  evt.set_event_scale(m_mur2); 
  return true;
}

void EventInfo::SetAlphaS()
{
  m_alphas=MODEL::s_model->ScalarFunction("alpha_S",m_mur2);
}

void EventInfo::SetAlpha()
{
  m_alpha=MODEL::s_model->ScalarConstant("alpha_QED");
}

HepMC2_Interface::HepMC2_Interface() :
  m_usenamedweights(false),
  m_extendedweights(false),
  m_includemeonlyweights(false),
  m_hepmctree(false),
  p_event(NULL)
{
  Settings& s = Settings::GetMainSettings();
#ifdef HEPMC_HAS_NAMED_WEIGHTS
  m_usenamedweights =
    s["HEPMC_USE_NAMED_WEIGHTS"].SetDefault(false).Get<bool>();
#endif
  m_extendedweights =
    s["HEPMC_EXTENDED_WEIGHTS"].SetDefault(false).Get<bool>();
  m_includemeonlyweights =
    s["HEPMC_INCLUDE_ME_ONLY_VARIATIONS"].SetDefault(false).Get<bool>();
  // Switch for disconnection of 1,2,3 vertices from PS vertices
  m_hepmctree = s["HEPMC_TREE_LIKE"].SetDefault(false).Get<bool>();
}

HepMC2_Interface::~HepMC2_Interface()
{
  delete p_event;
  DeleteGenSubEventList();
}


bool HepMC2_Interface::Sherpa2ShortHepMC(ATOOLS::Blob_List *const blobs,
                                         HepMC::GenEvent& event)
{
  const auto weight(blobs->Weight());
#ifdef USING__HEPMC2__UNITS
  event.use_units(HepMC::Units::GEV,
                  HepMC::Units::MM);
#endif
  Blob *sp(blobs->FindFirst(btp::Signal_Process));
  if (!sp) sp=blobs->FindFirst(btp::Hard_Collision);
  Blob *mp(blobs->FindFirst(btp::Hard_Collision));  
  if (!mp) event.set_mpi(-1);
  EventInfo evtinfo(sp,weight,
                    m_usenamedweights,m_extendedweights,m_includemeonlyweights);
  // when subevtlist, fill hepmc-subevtlist
  if (evtinfo.SubEvtList()) return SubEvtList2ShortHepMC(evtinfo);
  event.set_event_number(ATOOLS::rpa->gen.NumberOfGeneratedEvents());
  evtinfo.WriteTo(event);
  HepMC::GenVertex * vertex=new HepMC::GenVertex();
  std::vector<HepMC::GenParticle*> beamparticles, inparticles;
  for (ATOOLS::Blob_List::iterator blit=blobs->begin();
       blit!=blobs->end();++blit) {
    Blob* blob=*blit;
    if (m_ignoreblobs.count(blob->Type())) continue;
    for (int i=0;i<blob->NInP();i++) {
      if (blob->InParticle(i)->ProductionBlob()==NULL) {
        Particle* parton=blob->InParticle(i);
        auto flav = parton->Flav();
        if (flav.Kfcode() == kf_lepton) {
          if (sp->InParticle(0)->Flav().IsLepton()) {
            flav = sp->InParticle(0)->Flav();
          } else {
            flav = sp->InParticle(1)->Flav();
          }
        }
        HepMC::GenParticle* inpart;
        Sherpa2ShortHepMC(parton->Momentum(), flav, true, inpart);
        vertex->add_particle_in(inpart);
	inparticles.push_back(inpart);
        // distinct because SHRIMPS has no bunches for some reason
        if (blob->Type()==btp::Beam || blob->Type()==btp::Bunch) {
          beamparticles.push_back(inpart);
        }
      }
    }
    for (int i=0;i<blob->NOutP();i++) {
      if (blob->OutParticle(i)->DecayBlob()==NULL ||
	  m_ignoreblobs.count(blob->OutParticle(i)->DecayBlob()->Type())!=0) {
        Particle* parton=blob->OutParticle(i);
        HepMC::GenParticle* outpart;
        Sherpa2ShortHepMC(parton->Momentum(), parton->Flav(), false, outpart);
        vertex->add_particle_out(outpart);
      }
    }
  }
  event.add_vertex(vertex);
  if (beamparticles.empty() && inparticles.size()==2) {
    for (size_t j(0);j<2;++j) {
      HepMC::GenVertex *beamvertex = new HepMC::GenVertex();
      event.add_vertex(beamvertex);
      auto flav = (j == 0) ? rpa->gen.Beam1() : rpa->gen.Beam2();
      if (flav.Kfcode() == kf_lepton) {
        flav = sp->InParticle(j)->Flav();
      }
      HepMC::GenParticle* beampart;
      Sherpa2ShortHepMC(rpa->gen.PBeam(j), flav, true, beampart);
      beamparticles.push_back(beampart);
      beamvertex->add_particle_in(beampart);
      beamvertex->add_particle_out(inparticles[j]);
    }
  }
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
    if (sub->m_result==0. &&
	!(sub->IsReal() && m_subeventlist.empty())) continue;
    HepMC::GenVertex * subvertex(new HepMC::GenVertex());
    HepMC::GenEvent * subevent(new HepMC::GenEvent());
    // set the event number (could be used to identify correlated events)
    subevent->set_event_number(ATOOLS::rpa->gen.NumberOfGeneratedEvents());
    // assume that only 2->(n-2) processes, flip for Comix, flavs are correct
    HepMC::GenParticle *beamparticles[2]={NULL,NULL};
    for (size_t j(0);j<2;++j) {
      HepMC::GenVertex *beamvertex = new HepMC::GenVertex();
      subevent->add_vertex(beamvertex);
      auto flav = (j == 0) ? rpa->gen.Beam1() : rpa->gen.Beam2();
      if (flav.Kfcode() == kf_lepton) {
        flav = sub->p_fl[j];
      }
      HepMC::GenParticle* beampart;
      Sherpa2ShortHepMC(rpa->gen.PBeam(j), flav, true, beampart);
      beamparticles[j] = beampart;
      beamvertex->add_particle_in(beampart);
      double flip(sub->p_mom[j][0]<0.);
      Vec4D momentum {flip ? -1.0 * sub->p_mom[j] : sub->p_mom[j]};
      HepMC::GenParticle* inpart;
      Sherpa2ShortHepMC(momentum, sub->p_fl[j], true, inpart);
      subvertex->add_particle_in(inpart);
      beamvertex->add_particle_out(inpart);
    }
    subevent->set_beam_particles(beamparticles[0],beamparticles[1]);
    for (size_t j(2);j<sub->m_n;++j) {
      HepMC::GenParticle* outpart;
      Sherpa2ShortHepMC(sub->p_mom[j], sub->p_fl[j], false, outpart);
      subvertex->add_particle_out(outpart);
    }
    subevent->add_vertex(subvertex);
    // not enough info in subevents to set PDFInfo properly,
    // so set flavours and x1, x2 from the Signal_Process
    // reset muR, muF, alphaS, alpha
    subevtinfo.SetWeight(sub->m_result);
    subevtinfo.SetPartonicWeight(sub->m_mewgt);
    subevtinfo.SetMuR2(sub->m_mu2[stp::ren]);
    subevtinfo.SetMuF12(sub->m_mu2[stp::fac]);
    subevtinfo.SetMuF22(sub->m_mu2[stp::fac]);
    subevtinfo.SetAlphaS();
    subevtinfo.SetAlpha();
    subevtinfo.WriteTo(*subevent,i);
    m_subeventlist.push_back(subevent);
  }
  return true;
}


bool HepMC2_Interface::Sherpa2ShortHepMC(ATOOLS::Blob_List *const blobs)
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
  return Sherpa2ShortHepMC(blobs, *p_event);
}

bool HepMC2_Interface::Sherpa2ShortHepMC(const Vec4D& mom,
                                         const Flavour& flav,
                                         bool incoming,
                                         HepMC::GenParticle*& particle)
{
  HepMC::FourVector momentum(mom[1], mom[2], mom[3], mom[0]);
  int status = 1;
  if (incoming)
    status = (flav.StrongCharge() == 0) ? 4 : 11;
  particle = new HepMC::GenParticle(momentum, (long int)flav, status);
  return true;
}

// HS: Short-hand that takes a blob list, creates a new GenEvent and
// calls the actual Sherpa2HepMC
bool HepMC2_Interface::Sherpa2HepMC(ATOOLS::Blob_List *const blobs)
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
  return Sherpa2HepMC(blobs, *p_event);
}

// The actual code --- calls the Blob to GenVertex code
bool HepMC2_Interface::Sherpa2HepMC(ATOOLS::Blob_List *const blobs,
                                    HepMC::GenEvent& event)
{
  const auto weight(blobs->Weight());
  DEBUG_FUNC("");
#ifdef USING__HEPMC2__UNITS
  event.use_units(HepMC::Units::GEV,
                  HepMC::Units::MM);
#endif
  // Signal Process blob --- there is only one
  Blob *sp(blobs->FindFirst(btp::Signal_Process));
  if (!sp) sp=blobs->FindFirst(btp::Hard_Collision);
  // Meta info
  event.set_event_number(ATOOLS::rpa->gen.NumberOfGeneratedEvents());
  EventInfo evtinfo(sp,weight,
                    m_usenamedweights,m_extendedweights,m_includemeonlyweights);
  evtinfo.WriteTo(event);
  
  m_blob2genvertex.clear();
  m_particle2genparticle.clear();
  HepMC::GenVertex * vertex, * psvertex(NULL);
  std::vector<HepMC::GenParticle*> beamparticles;
  for (ATOOLS::Blob_List::iterator blit=blobs->begin();
       blit!=blobs->end();++blit) {
    // Call the Blob to vertex code, changes vertex pointer above
    if (Sherpa2HepMC(*(blit),vertex)) {
      event.add_vertex(vertex);
      for (size_t i(0);i<(*blit)->NInP();++i)
	if ((*blit)->InParticle(i)->ProductionBlob()==NULL) {
	  psvertex=vertex;
	  break;
	}
      if ((*blit)->Type()==ATOOLS::btp::Signal_Process) {
        if ((**blit)["NLO_subeventlist"]) {
          THROW(fatal_error,"Events containing correlated subtraction events"
                +std::string(" cannot be translated into the full HepMC event")
                +std::string(" format.\n")
                +std::string("   Try 'EVENT_OUTPUT: HepMC_Short' instead."));
        }
        event.set_signal_process_vertex(vertex);
      }
      // Find beam particles
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
  } // End Blob_List loop
  if (beamparticles.empty()) {
    HepMC::GenParticle* inpart[2];
    for (size_t j {0}; j < 2; ++j) {
      auto flav = (j == 0) ? rpa->gen.Beam1() : rpa->gen.Beam2();
      Sherpa2ShortHepMC(rpa->gen.PBeam(j), flav, true, inpart[j]);
      psvertex->add_particle_in(inpart[j]);
    }
    event.set_beam_particles(inpart[0], inpart[1]);
  }
  if (beamparticles.size()==2) {
    event.set_beam_particles(beamparticles[0],beamparticles[1]);
  }


  // Disconnect ME, MPI and hard decay vertices from PS vertices to get a
  // tree-like record --- manipulates the final GenEvent
  // Can't use ->set_production_vertex/set_end_vertex as they are private
  // Need to use GenVertex::remove_particle(Pointer to particle)
  // But: iterator loses validity when calling GenVertex::remove_particle in the particle loop
  // Hence: fill vector with pointers and call GenVertex::remove_particle 
  if (m_hepmctree) {
    DEBUG_INFO("HEPMC_TREE_LIKE true --- straighten to "
               <<"tree enabled (disconnect 1,2,3 vertices)");
    // Iterate over all vertices to find PS vertices
    int vtx_id = -1;
    for (HepMC::GenEvent::vertex_const_iterator vit=event.vertices_begin();
       vit!=event.vertices_end(); ++vit) {

      // Is this a PS Vertex?
      if ((*vit)->id()==4) {
        std::vector<HepMC::GenParticle*> remove;
        //// Loop over outgoing particles
        for (HepMC::GenVertex::particles_out_const_iterator pout
               =(*vit)->particles_out_const_begin();
             pout!=(*vit)->particles_out_const_end(); ++pout) {
            if ( (*pout)->end_vertex() ) {
              vtx_id = (*pout)->end_vertex()->id(); //
              // Disconnect outgoing particle from end/decay vertex of type (1,2,3)
              if (vtx_id==1 || vtx_id==2 || vtx_id==3 )
                  remove.push_back((*pout));
            }
        }
        // Loop over incoming particles
        for (HepMC::GenVertex::particles_in_const_iterator pin
               =(*vit)->particles_in_const_begin();
             pin!=(*vit)->particles_in_const_end(); ++pin) {
          vtx_id = (*pin)->production_vertex()->id();
          // Disconnect incoming particle from production vertex of type (1,2,3)
          if (vtx_id==1 || vtx_id==2 || vtx_id==3 )
                  remove.push_back((*pin));
        }
        // Iterate over Genparticle pointers to remove from current vertex (*vit)
        for (unsigned int nrem=0;nrem<remove.size();++nrem) {
            (*vit)->remove_particle(remove[nrem]);
        }
      } // Close if statement (vertex id==4)
    } // Close loop over vertices
  }
  return true;
}

// HS: this converts a Blob to GenVertex
bool HepMC2_Interface::Sherpa2HepMC(ATOOLS::Blob * blob, 
				    HepMC::GenVertex *& vertex)
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
    if (blob->Type()==btp::Signal_Process)      vertex->set_id(1); // signal
    else if (blob->Type()==btp::Hard_Collision) vertex->set_id(2); // mpi
    else if (blob->Type()==btp::Hard_Decay)     vertex->set_id(3); // hard-decay
    else if (blob->Type()==btp::Shower || 
             blob->Type()==btp::QED_Radiation)  vertex->set_id(4); // PS/QED
    else if (blob->Type()==btp::Fragmentation)  vertex->set_id(5); // frag
    else if (blob->Type()==btp::Hadron_Decay)   vertex->set_id(6); // had-decay
      //{  
      //if ((*blob)["Partonic"]!=NULL) vertex->set_id(-6);
      //else vertex->set_id(6);
      //}
    else vertex->set_id(0);
  }

  bool okay = 1;
  HepMC::GenParticle * _particle;
  for (int i=0;i<blob->NInP();i++) {
    if (Sherpa2HepMC(blob->InParticle(i),_particle)) {
      vertex->add_particle_in(_particle);
    }
    else okay = 0;
  }
  for (int i=0;i<blob->NOutP();i++) {
    if (Sherpa2HepMC(blob->OutParticle(i),_particle)) {
      vertex->add_particle_out(_particle);
    }
    else okay = 0;
  }
  m_blob2genvertex.insert(std::make_pair(blob,vertex));
  if (!okay) {
    msg_Error()<<"Error in HepMC2_Interface::Sherpa2HepMC(Blob,Vertex).\n"
               <<"    Continue event generation with new event."<<std::endl;
  }
  if (msg_LevelIsDebugging()) {
    ATOOLS::Vec4D check = blob->CheckMomentumConservation();
    double test         = ATOOLS::Vec3D(check).Abs();
    if (ATOOLS::dabs(1.-vertex->check_momentum_conservation()/test)>1.e-5 &&
        ATOOLS::dabs(test)>1.e-5) {
      msg_Error()<<"ERROR in "<<METHOD<<std::endl
                 <<"    Momentum not conserved. Continue."<<std::endl
                 <<"ERROR in Blob -> Vertex : "
                 <<vertex->check_momentum_conservation()
                 <<" <- "<<test<<" "<<check
                 <<std::endl<<(*blob)<<std::endl;
      vertex->print(msg_Error());
      msg_Error()<<"-----------------------------------------------"<<std::endl;
    }
  }
  return okay;
}

// HS: Sherpa Particle to HepMC::Genparticle --- fills m_particle2genparticle
// and changes the pointer reference particle ('new')
bool HepMC2_Interface::Sherpa2HepMC(ATOOLS::Particle * parton,
                                    HepMC::GenParticle *& particle)
{
  // HS: do nothing if parton has already been converted  
  int count = m_particle2genparticle.count(parton);
  if (count>0) {
    particle = m_particle2genparticle[parton];
    return true;
  }

  // Translate momentum vector
  ATOOLS::Vec4D mom  = parton->Momentum();
  HepMC::FourVector momentum(mom[1],mom[2],mom[3],mom[0]);

  int status=11;
  // Assign status 1 to stable blobs (those without decays)
  // or for that Rivet specific bit, set the particle stable (1)
  // if its DecayBlob has been cut out
  if (parton->DecayBlob()==NULL ||
      m_ignoreblobs.count(parton->DecayBlob()->Type())!=0) {
    status=1;
    DEBUG_VAR(*parton);
  }
  // Non-stable particles --- what about Hard_Decay?
  else {
    if (parton->DecayBlob()->Type()==ATOOLS::btp::Hadron_Decay ||
        parton->DecayBlob()->Type()==ATOOLS::btp::Hadron_Mixing) {
      status=2;
    }
    // Set all particles going in/out of ME to status 3
    else if (parton->DecayBlob()->Type()==ATOOLS::btp::Signal_Process ||
             (parton->ProductionBlob() &&
              parton->ProductionBlob()->Type()==ATOOLS::btp::Signal_Process)) {
      status=3;
    }
    // E - gamma collider specific
    else if (parton->DecayBlob()->Type()==ATOOLS::btp::Bunch) {
      status=4;
    }
  }
  if (parton->Status()==part_status::documentation) status=20;
  particle = new HepMC::GenParticle(momentum,(long int)parton->Flav(),status);
  for (int i=1;i<3;i++) {
    if (parton->GetFlow(i)>0) particle->set_flow(i,parton->GetFlow(i));
  }
  m_particle2genparticle.insert(std::make_pair(parton,particle));
  return true;
}

void HepMC2_Interface::AddCrossSection(HepMC::GenEvent& event,
                                       const double &xs, const double &err)
{
#ifdef HEPMC_HAS_CROSS_SECTION
  HepMC::GenCrossSection gxs;
  gxs.set_cross_section(xs,err);
  event.set_cross_section(gxs);
#else
  msg_Info()<<METHOD<<"(): Cannot add XS info to GenEvent."<<std::endl;
#endif
}

bool HepMC2_Interface::StartsLikeVariationName(const std::string& s)
{
  return (s.find("MUR") == 0 || s.find("ME_ONLY") == 0 || s.find("QCUT") == 0);
}

void HepMC2_Interface::DeleteGenSubEventList()
{
  for (size_t i=0; i<m_subeventlist.size();++i)
    delete m_subeventlist[i];
  m_subeventlist.clear();
}

#endif
