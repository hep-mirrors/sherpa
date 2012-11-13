#include"SHERPA/PerturbativePhysics/Hard_Decay_Handler.H"

#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/Return_Value.H"
#include "ATOOLS/Org/Shell_Tools.H"
#include "ATOOLS/Phys/Blob_List.H"
#include "ATOOLS/Phys/Cluster_Amplitude.H"
#include "ATOOLS/Math/Random.H"
#include "PHASIC++/Decays/Decay_Map.H"
#include "PHASIC++/Decays/Decay_Table.H"
#include "PHASIC++/Decays/Decay_Channel.H"

#include "MODEL/Main/Model_Base.H"
#include "MODEL/Interaction_Models/Single_Vertex.H"
#include "MODEL/Interaction_Models/Color_Function.H"
#include "PHASIC++/Decays/Color_Function_Decay.H"
#include "PHASIC++/Channels/Multi_Channel.H"
#include "PHASIC++/Channels/Rambo.H"
#include "METOOLS/Main/Spin_Structure.H"
#include "ATOOLS/Org/Run_Parameter.H"

#include "EXTRA_XS/One2Two/Comix1to2.H"
#include "EXTRA_XS/One2Three/Comix1to3.H"

#include <iostream>
#include <iomanip>
#include <sstream>
#include <cstring>
#include <assert.h>

using namespace SHERPA;
using namespace ATOOLS;
using namespace MODEL;
using namespace EXTRAXS;
using namespace PHASIC;
using namespace METOOLS;
using namespace std;

Hard_Decay_Handler::Hard_Decay_Handler(std::string path, std::string file)
{
  Data_Reader dr(" ",";","!","=");
  dr.AddWordSeparator("\t");
  dr.SetInputPath(path);
  dr.SetInputFile(file);
  m_mass_smearing=dr.GetValue<int>("HARD_MASS_SMEARING",1);
  m_spincorr=rpa->gen.HardSC();
  /*
    TODO: Writing out a 1->3 channel which might have a large width in one
    resolved configuration, and a small one in another?
    */
  m_store_results=dr.GetValue<int>("STORE_DECAY_RESULTS",0);
  m_decay_tau=dr.GetValue<int>("DECAY_TAU_HARD",0);
  m_set_widths=dr.GetValue<int>("HDH_SET_WIDTHS",0);
  m_resultdir=dr.GetValue<std::string>("RESULT_DIRECTORY","Results");
  if (m_store_results) {
    MakeDir("Results/Decays/", true);
  }
  m_offshell=dr.GetValue<std::string>("RESOLVE_DECAYS", "ByWidth");

  DEBUG_FUNC("");
  p_decaymap = new Decay_Map(this);
  KFCode_ParticleInfo_Map::const_iterator it;
  for (it=s_kftable.begin();it!=s_kftable.end();it++) {
    Flavour flav(it->first);
    if (Decays(flav)) {
      Decay_Table* dt=new Decay_Table(flav, this);
      vector<Decay_Table*> decaytables;
      decaytables.push_back(dt);
      p_decaymap->insert(make_pair(flav,decaytables));
      ReadDecayTable(flav);
    }
  }
  
  // initialize them sorted by masses:
  Decay_Map::iterator dmit;
  Vertex_Table offshell;
  msg_Info()<<"Initialising hard decay tables."<<endl;
  size_t i(0);
  for (dmit=p_decaymap->begin(); dmit!=p_decaymap->end(); ++dmit) {
    offshell.insert(make_pair(dmit->first,Vertex_List()));
    msg_Info()<<"  Initialising two-body decays. Step "
              <<++i<<"/"<<p_decaymap->size()<<" ("<<dmit->first<<")            "
              <<endl;
    InitializeDirectDecays(dmit->second.at(0));
  }
  i=0;
  if (p_decaymap->size()) msg_Info()<<endl;
  for (dmit=p_decaymap->begin(); dmit!=p_decaymap->end(); ++dmit) {
    msg_Info()<<"  Initialising three-body decays. Step "
              <<++i<<"/"<<p_decaymap->size()<<" ("<<dmit->first<<")            "
              <<endl;
    if (m_offshell=="None") {
      PRINT_INFO("Warning: Ignoring offshell decays as requested.");
    }
    else if (m_offshell=="Threshold") {
      RefineDecaysThreshold(dmit->second.at(0));
    }
    else if (m_offshell=="ByWidth") {
      RefineDecaysByWidth(dmit->second.at(0));
    }
    else {
      THROW(fatal_error, "Parameter RESOLVE_DECAYS set to wrong value.")
    }
    dmit->second.at(0)->UpdateWidth();
    if (m_set_widths)
      dmit->second.at(0)->Flav().SetWidth(dmit->second.at(0)->TotalWidth());
  }
  if (p_decaymap->size()) msg_Info()<<endl<<*p_decaymap<<endl;
  WriteDecayTables();
}

Hard_Decay_Handler::~Hard_Decay_Handler()
{
}

void Hard_Decay_Handler::InitializeDirectDecays(Decay_Table* dt)
{
  DEBUG_FUNC(dt->Flav());
  Flavour inflav=dt->Flav();
  Vertex_Table::const_iterator vlit=s_model->VertexTable()->find(inflav);
  const Vertex_List& vertexlist(vlit->second);

  for (size_t i=0;i<vertexlist.size();i++) {
    Single_Vertex* sv=vertexlist[i];
    if (!ProperVertex(sv)) continue;
    DEBUG_VAR(*sv);
    Decay_Channel* dc=new Decay_Channel(inflav, this);
    for (int j=1; j<sv->nleg; ++j) dc->AddDecayProduct(sv->in[j]);

    assert(sv->Color.size()==1);
    Comix1to2* diagram=new Comix1to2(dc->Flavs());
    dc->AddDiagram(diagram,new Color_Function_Decay(sv->Color[0]));

    dc->SetChannels(new Multi_Channel(""));
    dc->Channels()->SetNin(1);
    dc->Channels()->SetNout(dc->NOut());
    Rambo* rambo = new Rambo(1,dc->NOut(),&dc->Flavs().front(),this);
    dc->Channels()->Add(rambo);
    dc->Channels()->Reset();

    if (CalculateWidth(dc)) dt->AddDecayChannel(dc);
    else delete dc;
  }
  dt->UpdateWidth();
  if (m_set_widths) dt->Flav().SetWidth(dt->TotalWidth());
}

void Hard_Decay_Handler::RefineDecaysThreshold(Decay_Table* dt) {
  DEBUG_FUNC(dt->Flav()<<" "<<dt->size());
  size_t dtsize=dt->size();
  for (size_t i=0;i<dtsize;++i) {
    Decay_Channel* dc=dt->at(i);
    DEBUG_VAR(*dc);
    double outmass=0.0;
    for (size_t j=1; j<dc->Flavs().size(); ++j)
      outmass+=dc->Flavs()[j].Mass();
    if (dt->Flav().Mass()>outmass) continue;

    vector<Decay_Channel*> new_dcs=ResolveDecay(dc);
    for (size_t j=0; j<new_dcs.size(); ++j) {
      dt->AddDecayChannel(new_dcs[j]);
    }
    if (new_dcs.size()>0) {
      dc->SetActive(-1);
      for (size_t j=0; j<new_dcs.size(); ++j) {
        DEBUG_INFO("Adding "<<*new_dcs[j]);
      }
    }
    else {
      for (size_t j=0; j<new_dcs.size(); ++j) {
        new_dcs[j]->SetActive(-1);
      }
    }
  }
}

void Hard_Decay_Handler::RefineDecaysByWidth(Decay_Table* dt) {
  DEBUG_FUNC(dt->Flav()<<" "<<dt->size());
  size_t dtsize=dt->size();
  for (size_t i=0;i<dtsize;++i) {
    Decay_Channel* dc=dt->at(i);
    DEBUG_VAR(*dc);
    vector<Decay_Channel*> new_dcs=ResolveDecay(dc);
    double sum_resolved_widths=0.0;
    for (size_t j=0; j<new_dcs.size(); ++j) {
      Decay_Channel* dup=dt->GetDecayChannel(new_dcs[j]->Flavs());
      if (dup) {
        DEBUG_INFO("Duplicate for "<<*dup);
        for (DiagColVec::const_iterator it=new_dcs[j]->GetDiagrams().begin();
             it!=new_dcs[j]->GetDiagrams().end(); ++it) {
          dup->AddDiagram(it->first, it->second);
        }
        for (vector<Single_Channel*>::const_iterator it=new_dcs[j]->Channels()->Channels().begin();
             it!=new_dcs[j]->Channels()->Channels().end(); ++it) {
          dup->AddChannel(*it);
        }
        new_dcs[j]->ResetDiagrams();
        new_dcs[j]->ResetChannels();
        delete new_dcs[j];
        new_dcs[j]=NULL;
        sum_resolved_widths-=dup->Width();
        CalculateWidth(dup);
        sum_resolved_widths+=dup->Width();
      }
      else {
        sum_resolved_widths+=new_dcs[j]->Width();
        dt->AddDecayChannel(new_dcs[j]);
      }
    }
    DEBUG_INFO("resolved="<<sum_resolved_widths<<" vs. "<<dc->Width());

    if (sum_resolved_widths>dc->Width()) {
      dc->SetActive(-1);
    }
    else {
      for (size_t j=0; j<new_dcs.size(); ++j) {
        if (new_dcs[j]) new_dcs[j]->SetActive(-1);
      }
    }
  }
}

vector<Decay_Channel*> Hard_Decay_Handler::ResolveDecay(Decay_Channel* dc1)
{
  DEBUG_FUNC(*dc1);
  vector<Decay_Channel*> new_dcs;
  const std::vector<ATOOLS::Flavour> flavs1(dc1->Flavs());
  for (size_t j=1;j<flavs1.size();++j) {
    bool ignore=false;
    for (size_t k=1; k<j; ++k) {
      // TODO Do we really have to avoid double counting e.g. in h -> Z Z?
      // Further iterations: W+ -> b t -> b W b -> b b .. ?
      if (flavs1[j]==flavs1[k]) ignore=true;
    }
    if (ignore) continue;
    Vertex_Table::const_iterator it=s_model->VertexTable()->find(flavs1[j]);
    const Vertex_List& vertexlist(it->second);
    for (size_t k=0;k<vertexlist.size();k++) {
      Single_Vertex* sv = vertexlist[k];
      if (!ProperVertex(sv)) continue;
      // TODO so far special case 1->3 only
      Decay_Channel* dc=new Decay_Channel(flavs1[0], this);
      size_t nonprop(0), propi(0), propj(0);
      dc->AddDecayProduct(flavs1[3-j]);
      dc->AddDecayProduct(sv->in[1]);
      dc->AddDecayProduct(sv->in[2]);
      DEBUG_FUNC("trying "<<*dc);
      // TODO what about W+ -> b t -> b W b -> b b ... two diagrams, factor 2, ...?
      // TODO what about W' -> b t -> b W b ... two diagrams, factor 2, ...?
      for (size_t l=1; l<4; ++l) {
        if (dc->Flavs()[l]==flavs1[3-j]) nonprop=l;
      }
      for (size_t l=1; l<4; ++l) {
        if (l!=nonprop && propi>0) propj=l;
        else if (l!=nonprop && propi==0) propi=l;
      }

      assert(dc1->GetDiagrams().size()==1);
      assert(sv->Color.size()==1);
      DEBUG_VAR(dc->Flavs());
      DEBUG_VAR(flavs1[j]);
      Comix1to3* diagram=new Comix1to3(dc->Flavs(),flavs1[j],
                                       nonprop, propi, propj);
      Color_Function_Decay* col1=new Color_Function_Decay(*dc1->GetDiagrams()[0].second);
      DEBUG_VAR(*col1);
      Color_Function_Decay col2(sv->Color[0]);
      DEBUG_VAR(col2);
      vector<int> bumps = col1->Multiply(col2);
      DEBUG_VAR(*col1);
      DEBUG_INFO("Contracting "<<0+bumps[0]<<" with "<<j);
      col1->Contract(0+bumps[0],j);
      DEBUG_VAR(*col1);
      dc->AddDiagram(diagram, col1);

      dc->SetChannels(new Multi_Channel(""));
      dc->Channels()->SetNin(1);
      dc->Channels()->SetNout(dc->NOut());
      Rambo* rambo = new Rambo(1,dc->NOut(),&dc->Flavs().front(),this);
      dc->Channels()->Add(rambo);
      dc->Channels()->Reset();

      if (CalculateWidth(dc)) new_dcs.push_back(dc);
      else delete dc;
    }
  }

  return new_dcs;
  /* TODO:
     What about cases where a 1->3 was resolved from one 1->2, but would have
     been possible from another 1->2 where it was decided not to be resolved?
     */
}

bool Hard_Decay_Handler::CalculateWidth(Decay_Channel* dc)
{
  // Integrate or use results read from decay table file
  double outmass=0.0;
  for (size_t l=1; l<dc->Flavs().size(); ++l)
    outmass+=dc->Flavs()[l].Mass();
  if (dc->Flavs()[0].Mass()>outmass) {
    DEBUG_INFO("Starting calculation of widths now.");
    DEBUG_INFO(*dc);
    if (m_store_results && m_read.find(dc->Flavs()[0])!=m_read.end()) {
      if (m_read[dc->Flavs()[0]].find(dc->IDCode())!=
          m_read[dc->Flavs()[0]].end()) {
        const vector<double>& results(m_read[dc->Flavs()[0]][dc->IDCode()]);
        dc->SetIWidth(results[0]);
        dc->SetIDeltaWidth(results[1]);
        dc->SetMax(results[2]);
        dc->SetActive(int(results[3]+0.5));
      }
      else {
        msg_Info()<<"    Integrating "<<dc->Name()<<endl;
        dc->CalculateWidth();
      }
    }
    else {
      msg_Info()<<"    Integrating "<<dc->Name()<<endl;
      dc->CalculateWidth();
    }
  }
  else {
    dc->SetIWidth(0.0);
    dc->SetIDeltaWidth(0.0);
    dc->SetMax(0.0);
  }
  dc->SetWidth(dc->IWidth());
  dc->SetDeltaWidth(dc->IDeltaWidth());
  return true;
}

bool Hard_Decay_Handler::ProperVertex(MODEL::Single_Vertex* sv)
{
  if (!sv->on || sv->dec) return false;

  for (int i(0); i<sv->nleg; ++i)
    if (sv->in[i].IsDummy()) return false;

  if (sv->nleg!=3) return false; // TODO

  // ignore radiation graphs. should we?
  for (int i=1; i<sv->nleg; ++i) {
    if (sv->in[i].Kfcode()==sv->in[0].Kfcode()) {
      return false;
    }
  }

  // what about extra particles like Z4 if Z stable?

  return true;
}


void Hard_Decay_Handler::CreateDecayBlob(ATOOLS::Particle* inpart)
{
  DEBUG_FUNC(inpart->Flav());
  if(inpart->DecayBlob()) abort();
  if(!Decays(inpart->Flav())) return;
  Blob* blob = p_bloblist->AddBlob(btp::Hard_Decay);
  blob->SetStatus(blob_status::needs_showers);
  blob->AddToInParticles(inpart);
  blob->SetTypeSpec("Sherpa");
  Decay_Table* table=p_decaymap->FindDecay(blob->InParticle(0)->Flav());
  if (table==NULL) {
    msg_Error()<<METHOD<<" decay table not found, retrying event."<<endl
               <<*blob<<endl;
    throw Return_Value::Retry_Event;
  }
  blob->AddData("dc",new Blob_Data<Decay_Channel*>(table->Select()));

  DEBUG_INFO("p_onshell="<<inpart->Momentum());
  blob->AddData("p_onshell",new Blob_Data<Vec4D>(inpart->Momentum()));
  DEBUG_INFO("succeeded.");
}

void Hard_Decay_Handler::DefineInitialConditions(Cluster_Amplitude* ampl,
                                                 Blob* initial_blob)
{
  DEBUG_FUNC(*ampl);
  for (int i=0; i<initial_blob->NOutP(); ++i) {
    ampl->Leg(initial_blob->NInP()+i)->SetMom
      (initial_blob->OutParticle(i)->Momentum());
  }
  size_t imax=ampl->Legs().size()-1;
  for (int i=0; i<initial_blob->NOutP(); ++i) {
    if (initial_blob->OutParticle(i)->DecayBlob()) {
      AddDecayClustering(ampl, initial_blob->OutParticle(i)->DecayBlob(),
                         imax, 1<<(initial_blob->NInP()+i));
    }
  }
}

void Hard_Decay_Handler::AddDecayClustering(ATOOLS::Cluster_Amplitude*& ampl,
                                            Blob* blob,
                                            size_t& imax,
                                            size_t idmother)
{
  DEBUG_FUNC("blob->Id()="<<blob->Id()<<" idmother="<<idmother);
  DEBUG_VAR(*blob);
  Particle_Vector daughters=blob->GetOutParticles();
  if (daughters.size()==2) {
    Cluster_Amplitude* copy=ampl->InitPrev();
    copy->CopyFrom(ampl);
    Cluster_Leg *lij(ampl->IdLeg(idmother));
    lij->SetStat(3);
    size_t idk(0);
    for (size_t i=0; i<copy->Legs().size(); ++i) {
      copy->Leg(i)->SetK(0);
      if (copy->Leg(i)->Col().m_i==lij->Col().m_j ||
          copy->Leg(i)->Col().m_j==lij->Col().m_i) 
        idk=copy->Leg(i)->Id();
    }
    if (lij->Col().m_i==0 && lij->Col().m_j==0) {
      if (ampl->Next()) {
        Cluster_Amplitude *next(ampl->Next());
        for (size_t i=0;i<next->Legs().size();++i) {
          if (next->Leg(i)->Id()&lij->Id()) {
            size_t id(next->Leg(i)->Id()^lij->Id());
            for (size_t j=0;j<ampl->Legs().size();++j)
              if (ampl->Leg(j)->Id()&id) {
              idk=ampl->Leg(j)->Id();
              break;
            }
            break;
          }
        }
      }
      else {
        // Ad hoc EW partner
        size_t ampl_nout=ampl->Legs().size()-ampl->NIn();
        size_t select=ampl->Legs().size();
        do {
          select=ampl->NIn()+floor(ran->Get()*ampl_nout);
        } while (ampl->Leg(select)->Id()&idmother || select>ampl->Legs().size()-1);
        idk=ampl->Leg(select)->Id();
      }
    }
    if (idk==0) THROW(fatal_error,"Colour partner not found");
    lij->SetK(idk);
    Cluster_Leg *d1(copy->IdLeg(idmother));
    d1->SetMom(daughters[0]->Momentum());
    d1->SetFlav(daughters[0]->Flav());
    copy->CreateLeg(daughters[1]->Momentum(), daughters[1]->RefFlav());
    size_t idnew=1<<(++imax);
    copy->Legs().back()->SetId(idnew);
    copy->Legs().back()->SetStat(1);
    Cluster_Amplitude::SetColours(ampl->IdLeg(idmother),
                                  copy->IdLeg(idmother),
                                  copy->Legs().back());
    
    DEBUG_VAR(*copy);
    Cluster_Amplitude* tmp=copy;
    while (tmp->Next()) {
      tmp=tmp->Next();
      for (size_t i=0; i<tmp->Legs().size(); ++i) {
        if (tmp->Leg(i)->Id()&idmother) {
          tmp->Leg(i)->SetId(tmp->Leg(i)->Id()|idnew);
	}
        if (tmp->Leg(i)->K()&idmother) {
          tmp->Leg(i)->SetK(tmp->Leg(i)->K()|idnew);
        }
      }
      DEBUG_VAR(*tmp);
    }
    if (daughters[0]->DecayBlob())
      AddDecayClustering(copy, daughters[0]->DecayBlob(), imax, idmother);
    if (daughters[1]->DecayBlob())
      AddDecayClustering(copy, daughters[1]->DecayBlob(), imax, idnew);
    ampl=copy;
  }
  else if (daughters.size()==3) {
    DEBUG_VAR("size=3");
    Cluster_Amplitude* step1=ampl->InitPrev();
    step1->CopyFrom(ampl);
    Cluster_Leg *lij(ampl->IdLeg(idmother));
    lij->SetStat(3);
    size_t idk(0);
    for (size_t i=0; i<step1->Legs().size(); ++i) {
      step1->Leg(i)->SetK(0);
      if (step1->Leg(i)->Col().m_i==lij->Col().m_j ||
          step1->Leg(i)->Col().m_j==lij->Col().m_i) 
        idk=step1->Leg(i)->Id();
    }
    if (lij->Col().m_i==0 && lij->Col().m_j==0) {
      if (ampl->Next()) {
        Cluster_Amplitude *next(ampl->Next());
        for (size_t i=0;i<next->Legs().size();++i) {
          if (next->Leg(i)->Id()&lij->Id()) {
            size_t id(next->Leg(i)->Id()^lij->Id());
            for (size_t j=0;j<ampl->Legs().size();++j)
              if (ampl->Leg(j)->Id()&id) {
              idk=ampl->Leg(j)->Id();
              break;
            }
            break;
          }
        }
      }
      else {
        // Ad hoc EW partner
        size_t ampl_nout=ampl->Legs().size()-ampl->NIn();
        size_t select=ampl->Legs().size();
        do {
          select=ampl->NIn()+floor(ran->Get()*ampl_nout);
        } while (ampl->Leg(select)->Id()&idmother || select>ampl->Legs().size()-1);
        idk=ampl->Leg(select)->Id();
      }
    }
    if (idk==0) THROW(fatal_error,"Colour partner not found");
    lij->SetK(idk);
    Cluster_Leg *d1(step1->IdLeg(idmother));
    d1->SetMom(daughters[0]->Momentum());
    d1->SetFlav(daughters[0]->Flav());
    // todo: 1->2 qcd shower with ew fs recoil partner
    // d1->SetK(idmother);// not that simple: w->qq' has color connection in fs
    Decay_Channel* dc(NULL);
    Blob_Data_Base* data = (*blob)["dc"];
    if(data) {
      dc=data->Get<Decay_Channel*>();
      DEBUG_VAR(*dc);
    }
    else THROW(fatal_error, "Internal error.");
    Comix1to3* amp=dynamic_cast<Comix1to3*>(dc->GetDiagrams()[0].first);
    if (!amp) THROW(fatal_error, "Internal error.");
    Flavour prop_flav=amp->Prop();
    Vec4D prop_mom=daughters[1]->Momentum()+daughters[2]->Momentum();
    step1->CreateLeg(prop_mom, prop_flav);
    size_t idnew1=1<<(++imax);
    step1->Legs().back()->SetId(idnew1);
    step1->Legs().back()->SetStat(0);
    Cluster_Amplitude::SetColours(ampl->IdLeg(idmother),
                                  step1->IdLeg(idmother),
                                  step1->Legs().back());
    DEBUG_VAR(*step1);
    Cluster_Amplitude* tmp=step1;
    while (tmp->Next()) {
      tmp=tmp->Next();
      for (size_t i=0; i<tmp->Legs().size(); ++i) {
        if (tmp->Leg(i)->Id()&idmother) {
          tmp->Leg(i)->SetId(tmp->Leg(i)->Id()|idnew1);
	}
        if (tmp->Leg(i)->K()&idmother) {
          tmp->Leg(i)->SetK(tmp->Leg(i)->K()|idnew1);
        }
      }
      DEBUG_VAR(*tmp);
    }

    
    Cluster_Amplitude* step2=step1->InitPrev();
    step2->CopyFrom(step1);
    step1->IdLeg(idnew1)->SetStat(0);
    for (size_t i=0; i<step2->Legs().size(); ++i) step2->Leg(i)->SetK(0);
    Cluster_Leg *d2(step2->IdLeg(idnew1));
    d2->SetMom(daughters[1]->Momentum());
    d2->SetFlav(daughters[1]->Flav());
    step2->CreateLeg(daughters[2]->Momentum(), daughters[2]->Flav());
    size_t idnew2=1<<(++imax);
    step2->Legs().back()->SetId(idnew2);
    step2->Legs().back()->SetStat(0);
    Cluster_Amplitude::SetColours(step1->IdLeg(idnew1),
                                  step2->IdLeg(idnew1),
                                  step2->Legs().back());
    DEBUG_VAR(*step2);
    tmp=step2;
    while (tmp->Next()) {
      tmp=tmp->Next();
      for (size_t i=0; i<tmp->Legs().size(); ++i) {
        if (tmp->Leg(i)->Id()&idnew1) {
          tmp->Leg(i)->SetId(tmp->Leg(i)->Id()|idnew2);
	}
        if (tmp->Leg(i)->K()&idnew1) {
          tmp->Leg(i)->SetK(tmp->Leg(i)->K()|idnew2);
        }
      }
      DEBUG_VAR(*tmp);
    }

    
    if (daughters[0]->DecayBlob())
      AddDecayClustering(step2, daughters[0]->DecayBlob(), imax, idmother);
    if (daughters[1]->DecayBlob())
      AddDecayClustering(step2, daughters[1]->DecayBlob(), imax, idnew1);
    if (daughters[2]->DecayBlob())
      AddDecayClustering(step2, daughters[2]->DecayBlob(), imax, idnew2);
    ampl=step2;
  }
  else {
    PRINT_VAR(*blob);
    THROW(fatal_error, "1 -> n not implemented yet.");
  }
}

void Hard_Decay_Handler::ReadDecayTable(Flavour decayer)
{
  DEBUG_FUNC(decayer);
  if (!m_store_results) return;
  Data_Reader reader = Data_Reader("|",";","!");
  reader.SetAddCommandLine(false);
  reader.AddComment("#");
  reader.AddComment("//");
  reader.AddWordSeparator("\t");
  reader.SetInputPath("Results/Decays/");
  reader.SetInputFile(decayer.ShellName());
  
  vector<vector<string> > file;
  if(reader.MatrixFromFile(file)) {
    for (size_t iline=0; iline<file.size(); ++iline) {
      if (file[iline].size()==5) {
        string decaychannel=file[iline][0];
        vector<double> results(4);
        for (size_t i=0; i<4; ++i) results[i]=ToType<double>(file[iline][i+1]);
        m_read[decayer].insert(make_pair(decaychannel, results));
      }
      else {
        PRINT_INFO("Wrong format in decay table in Results directory.");
      }
    }
  }
}

void Hard_Decay_Handler::WriteDecayTables()
{
  if (!m_store_results) return;
  
  Decay_Map::iterator dmit;
  for (dmit=p_decaymap->begin(); dmit!=p_decaymap->end(); ++dmit) {
    ofstream ostr(("Results/Decays/"+dmit->first.ShellName()).c_str());
    ostr<<"# Decay table for "<<dmit->first<<endl<<endl;
    Decay_Table::iterator dtit;
    for (dtit=dmit->second[0]->begin(); dtit!=dmit->second[0]->end(); ++dtit) {
      ostr<<(*dtit)->IDCode()<<"\t"<<(*dtit)->IWidth()<<"\t"
          <<(*dtit)->IDeltaWidth()<<"\t"<<(*dtit)->Max()<<"\t"
          <<(*dtit)->Active()<<endl;
    }
    ostr.close();
  }
}

bool Hard_Decay_Handler::Decays(const ATOOLS::Flavour& flav)
{
  if (!flav.IsOn()) return false;
  if (flav.IsHadron()) return false;
  if (flav.Kfcode()==kf_tau && !m_decay_tau) return false;
  if (flav.IsStable()) return false;
  return true;
}
