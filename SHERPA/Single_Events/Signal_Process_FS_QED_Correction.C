#include "SHERPA/Single_Events/Signal_Process_FS_QED_Correction.H"

#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Math/Vector.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Phys/Blob.H"
#include "ATOOLS/Phys/Blob_List.H"
#include "ATOOLS/Phys/Cluster_Amplitude.H"
#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Phys/Momenta_Stretcher.H"
#include "ATOOLS/Phys/Particle.H"
#include "MODEL/Main/Model_Base.H"
#include "MODEL/Interaction_Models/Single_Vertex.H"


using namespace SHERPA;
using namespace ATOOLS;
using namespace PHASIC;
using namespace MODEL;
using namespace std;

////////////////////////////////////////////////////////////////////////////////
////                                                                        ////
////     in the documentation LEPTON is synonymous                          ////
////     to EVERYTHING NOT STRONGLY CHARGED                                 ////
////                                                                        ////
////////////////////////////////////////////////////////////////////////////////


Signal_Process_FS_QED_Correction::Signal_Process_FS_QED_Correction
(Matrix_Element_Handler *_mehandler, Soft_Photon_Handler *_sphotons) :
  m_on(true), m_qed(true), m_findresonances(true), m_resdist(1.),
  p_mehandler(_mehandler), p_sphotons(_sphotons)
{
  DEBUG_FUNC("");
  m_name      = string("Lepton_FS_QED_Corrections:");
  m_type      = eph::Perturbative;
  // general switch
  Data_Reader reader(" ",";","!","=");
  reader.AddComment("#");
  reader.AddWordSeparator("\t");
  reader.SetInputFile(rpa->gen.Variable("ME_DATA_FILE"));
  bool impliciteon = (reader.GetValue<std::string>("ME_QED","On")=="On");
  bool expliciteon = (reader.GetValue<std::string>("ME_QED","")=="On");
  // look whether there is any hadronisation following
  // if not, do not even put them on-shell -> switch everthing off
  if (!impliciteon) {
    m_qed = false;
    Data_Reader reader1(" ",";","!","=");
    reader1.AddComment("#");
    reader1.AddWordSeparator("\t");
    reader1.SetInputFile(rpa->gen.Variable("FRAGMENTATION_DATA_FILE"));
    m_on = (reader1.GetValue<std::string>("FRAGMENTATION","")!="Off");
  }
  // switch off if there are hard decays, have their own QED corrections,
  // cannot tell here what has been corrected and what not -- OR --
  // if NLO_Mode Fixed_Order, switch off completely, unless explicitely stated
  if (!expliciteon &&
      (p_mehandler->HasNLO()==1 ||
       !(reader.GetValue<std::string>("HARD_DECAYS","Off")=="Off"))) {
    m_on = false; m_qed = false;
  }
  msg_Debugging()<<"on="<<m_on<<" ,  qed="<<m_qed<<std::endl;

  // read in resonance finding parameters
  m_findresonances = (reader.GetValue<std::string>("ME_QED_CLUSTERING","On")
                                                                        =="On");
  m_resdist        = reader.GetValue<double>("ME_QED_CLUSTERING_THRESHOLD",1.);

  if (m_on && m_qed) {
    m_name += p_sphotons->SoftQEDGenerator();
  }
  else
    m_name += "None";

  // do not do the rest if not needed
  if (!m_qed) return;

  // identify non-QCD subprocesses
  Process_Vector pvec(p_mehandler->AllProcesses());
  for (size_t i=0;i<pvec.size();++i) {
    for (size_t j=0;j<pvec[i]->Size();++j) {
      SubInfoVector siv;
      FindSubProcessInfosContainingLeptons((*pvec[i])[j]->Info(),siv);
      msg_Debugging()<<"Process: "<<(*pvec[i])[j]->Name()<<" -> "
                     <<siv.size()<<" non-QCD production subprocesses.\n";
      for (size_t k=0;k<siv.size();++k) msg_Debugging()<<*siv[k]<<endl;
      m_proc_lep_map.insert(make_pair((*pvec[i])[j]->Name(),siv));
    }
  }

  // extract non-QCD resonances of the model to find resonant unresolved
  // resonant lepton production
  for (size_t i=0;i<pvec.size();++i) {
    for (size_t j=0;j<pvec[i]->Size();++j) {
      Vertex_List vlist;
      FindProcessPossibleResonances((*pvec[i])[j]->Flavours(),vlist);
      msg_Debugging()<<"Process: "<<(*pvec[i])[j]->Name()<<" -> "
                     <<vlist.size()<<" non-QCD resonances.\n";
      for (size_t k=0;k<vlist.size();++k) msg_Debugging()<<*vlist[k]<<endl;
      m_proc_restab_map[(*pvec[i])[j]->Name()]=vlist;
    }
  }

  // initialise photons helpers
  if(s_kftable.find(kf_PhotonsHelperNeutral)==s_kftable.end())
    s_kftable[kf_PhotonsHelperNeutral]
        =new Particle_Info(kf_PhotonsHelperNeutral,0.0,0,0,0,
                           "PH0","H_P^{0}");
  if(s_kftable.find(kf_PhotonsHelperPlus)==s_kftable.end())
    s_kftable[kf_PhotonsHelperPlus]
        =new Particle_Info(kf_PhotonsHelperPlus,0.0,3,0,0,
                           "PH+","H_P^{+}");
  if(s_kftable.find(kf_PhotonsHelperPlusPlus)==s_kftable.end())
    s_kftable[kf_PhotonsHelperPlusPlus]
        =new Particle_Info(kf_PhotonsHelperPlusPlus,0.0,6,0,0,
                           "PH++","H_P^{++}");
  if(s_kftable.find(kf_PhotonsHelperPlusPlusPlus)==s_kftable.end())
    s_kftable[kf_PhotonsHelperPlusPlusPlus]
        =new Particle_Info(kf_PhotonsHelperPlusPlusPlus,0.0,9,0,0,
                           "PH+++","H_P^{+++}");
}

Signal_Process_FS_QED_Correction::~Signal_Process_FS_QED_Correction() {}


Return_Value::code Signal_Process_FS_QED_Correction::Treat
(Blob_List * bloblist, double & weight)
{
  if (!m_on) return Return_Value::Nothing;
  if (bloblist->empty()) {
    msg_Error()<<"Signal_Process_FS_QED_Correction::Treat"
	       <<"("<<bloblist<<","<<weight<<"): "<<endl
               <<"   Blob list contains "<<bloblist->size()<<" entries."<<endl
               <<"   Continue and hope for the best."<<endl;
    return Return_Value::Error;
  }
  // look for QCD corrected hard process in need for QED
  Blob * sigblob(bloblist->FindLast(btp::Shower));
  if (!sigblob) return Return_Value::Nothing;
  // if already treated -> nothing to do
  if (sigblob->TypeSpec()=="YFS-type_QED_Corrections_to_ME")
    return Return_Value::Nothing;
  if (sigblob->TypeSpec()=="setting_leptons_on-shell")
    return Return_Value::Nothing;
  // extract FS leptons
  // two vectors -> the ones from the blob and the ones to be massive
  DEBUG_FUNC(m_qed);
  Particle_Vector fslep(sigblob->GetOutParticles());
  Particle_Vector mfslep;
  for (Particle_Vector::iterator it=fslep.begin();it!=fslep.end();) {
    if ((*it)->Flav().Strong() || (*it)->Flav().IsDiQuark() || 
	(*it)->DecayBlob()!=NULL) {
      fslep.erase(it);
    }
    else {
      mfslep.push_back(new Particle(-1,(*it)->Flav(),(*it)->Momentum(),'F'));
      (*mfslep.rbegin())->SetNumber(0);
      (*mfslep.rbegin())->SetOriginalPart(*it);
      (*mfslep.rbegin())->SetFinalMass((*it)->FinalMass());
      ++it;
    }
  }
  // if no leptons, nothing to do
  // if only one lepton, cannot do anything
  if (fslep.size()<2) {
    sigblob->UnsetStatus(blob_status::needs_extraQED);
    for (Particle_Vector::iterator it=mfslep.begin();it!=mfslep.end();++it)
      delete *it;
    return Return_Value::Nothing;
  }
  // if switched off or no need for QED stop here and build a blob
  if (!m_qed || !sigblob->Has(blob_status::needs_extraQED)) {
    Blob * onshellblob = bloblist->AddBlob(btp::QED_Radiation);
    onshellblob->SetTypeSpec("setting_leptons_on-shell");
    if (sigblob->Has(blob_status::needs_extraQED))
      sigblob->UnsetStatus(blob_status::needs_extraQED);
    for (Particle_Vector::iterator it=fslep.begin();it!=fslep.end();++it) {
      (*it)->SetInfo('H');
      (*it)->SetStatus(part_status::decayed);
      onshellblob->AddToInParticles(*it);
    }
    for (Particle_Vector::iterator it=mfslep.begin();it!=mfslep.end();++it) {
      onshellblob->AddToOutParticles(*it);
    }
    onshellblob->SetStatus(blob_status::needs_hadronization);
    return Return_Value::Success;
  }
  // put them on-shell (spoils consistency of pertubative calculation,
  // but necessary for YFS)
  if (!PutOnMassShell(mfslep)) {
    msg_Error()<<"Signal_Process_FS_QED_Correction::Treat("
	       <<bloblist<<","<<weight<<"): \n"
               <<"  Leptons could not be put on their mass shell.\n"
               <<"  Trying new event.\n"
               <<"  The event contained a ";
    for (Particle_Vector::iterator it=mfslep.begin();it!=mfslep.end();++it)
       msg_Error()<<(*it)->Flav().ShellName()<<"-";
    if (mfslep.size()==2) msg_Error()<<"pair";
    else                  msg_Error()<<"set";
    msg_Error()<<" of too little invariant mass to be put\n"
	       <<"  on their mass shell. If you are sensitive to this specific"
	       <<" signature consider\n  to set the respective particles"
	       <<" massive in the perturbative calculation using\n"
	       <<"  'MASSIVE[<id>]=1' to avoid this problem.\n";
    for (Particle_Vector::iterator it=mfslep.begin();it!=mfslep.end();++it)
      delete *it;
    return Return_Value::New_Event;
  }
  // build effective verteces for resonant production
  // use subprocess infos if possible
  Blob_Vector blobs = BuildResonantBlobs(mfslep);
  // add radiation
  for (Blob_Vector::iterator it=blobs.begin();it!=blobs.end();++it) {
    // do nothing if no resonance determined
    if ((*it)->InParticle(0)->Flav().Kfcode()!=kf_none) {
      (*it)->SetStatus(blob_status::needs_extraQED);
      if (!p_sphotons->AddRadiation(*it)) {
        msg_Error()<<"Signal_Process_FS_QED_Correction::Treat("<<bloblist
                   <<","<<weight<<"): "<<endl
                   <<"  Higher order QED corrections failed."<<endl
                   <<"  Retrying event."<<endl;
        for (Particle_Vector::iterator it=mfslep.begin();it!=mfslep.end();++it)
          delete *it;
        for (Blob_Vector::iterator it=blobs.begin();it!=blobs.end();++it)
          delete *it;
        return Return_Value::Retry_Event;
      }
    }
  }
  for (Blob_Vector::iterator it=blobs.begin();it!=blobs.end();++it) {
    msg_Debugging()<<**it<<endl;
    (*it)->DeleteInParticles();
  }
  sigblob->UnsetStatus(blob_status::needs_extraQED);
  // build new QED radiation blob
  Blob * QEDblob = bloblist->AddBlob(btp::QED_Radiation);
  QEDblob->SetTypeSpec("YFS-type_QED_Corrections_to_ME");
  for (Particle_Vector::iterator it=fslep.begin();it!=fslep.end();++it) {
    // set info back to hard process, otherwise
    // check for momentum conservation does not work
    (*it)->SetInfo('H');
    (*it)->SetStatus(part_status::decayed);
    QEDblob->AddToInParticles(*it);
  }
  // first fill in all LO particles
  for (Blob_Vector::iterator it=blobs.begin();it!=blobs.end();++it) {
    while ((*it)->NOutP() && (*it)->OutParticle(0)->Info()!='S') {
      Particle * part((*it)->RemoveOutParticle(0,true));
      QEDblob->AddToOutParticles(part);
    }
  }
  // then append all photons
  for (Blob_Vector::iterator it=blobs.begin();it!=blobs.end();++it) {
    while ((*it)->NOutP()) {
      Particle * part((*it)->RemoveOutParticle(0,true));
      QEDblob->AddToOutParticles(part);
    }
  }
  QEDblob->SetStatus(blob_status::needs_hadronization);
  for (size_t i=0;i<blobs.size();++i) {
    delete blobs[i];
    blobs[i]=NULL;
  }
  return Return_Value::Success;
}

bool Signal_Process_FS_QED_Correction::PutOnMassShell
(const Particle_Vector& partvec)
{
  // if massless in ME put on mass shell for YFS
  bool allonshell(true); kf_code kfc;
  std::vector<double>masses(partvec.size(),0.);
  for (size_t i=0;i<partvec.size();++i) {
    kfc=partvec[i]->Flav().Kfcode();
    if(kfc==kf_graviton || kfc==kf_gscalar)
      masses[i]=sqrt(fabs(partvec[i]->Momentum().Abs2()));
    else masses[i]=partvec[i]->Flav().Mass(1);
    //If one of the two squared masses is zero, IsEqual always returns 0.
    if (!IsEqual(partvec[i]->Momentum().Abs2(),sqr(masses[i]),1E-4))
      allonshell=false;
  }
  if (allonshell) return true;
  Momenta_Stretcher momstretch;
  return momstretch.StretchMomenta(partvec,masses);
}

Flavour Signal_Process_FS_QED_Correction::DetermineGenericResonance
(const Particle_Vector& partvec)
{
  int chargesum(0);
  for (size_t i=0;i<partvec.size();++i)
    chargesum+=partvec[i]->Flav().IntCharge();
  if      (chargesum==0)  return Flavour(kf_PhotonsHelperNeutral);
  else if (chargesum==3)  return Flavour(kf_PhotonsHelperPlus);
  else if (chargesum==-3) return Flavour(kf_PhotonsHelperPlus).Bar();
  else if (chargesum==6)  return Flavour(kf_PhotonsHelperPlusPlus);
  else if (chargesum==-6) return Flavour(kf_PhotonsHelperPlusPlus).Bar();
  else if (chargesum==9)  return Flavour(kf_PhotonsHelperPlusPlusPlus);
  else if (chargesum==-9) return Flavour(kf_PhotonsHelperPlusPlusPlus).Bar();
  // i got no clue what this might be
  return Flavour(kf_none);
}

Vec4D Signal_Process_FS_QED_Correction::MomentumSum
(const Particle_Vector& partvec)
{
  Vec4D sum(0.,0.,0.,0.);
  for (size_t i=0;i<partvec.size();++i) {
    sum += partvec[i]->Momentum();
  }
  return sum;
}

void Signal_Process_FS_QED_Correction::FindSubProcessInfosContainingLeptons
(const Process_Info& pi, SubInfoVector& siv)
{
  // loop over FS of process -> find subprocs, that are subsequent decays
  for (size_t i=0;i<pi.m_fi.m_ps.size();++i) {
    if (pi.m_fi.m_ps[i].m_ps.size()>1)
      FindSubProcessInfosContainingLeptons(pi.m_fi.m_ps[i],siv);
  }
}

void Signal_Process_FS_QED_Correction::FindSubProcessInfosContainingLeptons
(const Subprocess_Info& spi, SubInfoVector& siv)
{
  // assume connected leptons are produced in same subprocess info
  size_t count(0), leps(0);
  for (size_t i=0;i<spi.m_ps.size();++i) {
    if (spi.m_ps[i].m_ps.size()==0) {
      count++;
      if (!spi.m_ps[i].m_fl.Strong()) leps++;
    }
    else {
      FindSubProcessInfosContainingLeptons(spi.m_ps[i],siv);
    }
  }
  if (count==spi.m_ps.size() && leps!=0) {
    // -> final subprocess info
    // if contains leptons, add to siv
    siv.push_back(&spi);
  }
}

void Signal_Process_FS_QED_Correction::FindProcessPossibleResonances
(const Flavour_Vector& fv, MODEL::Vertex_List& vlist)
{
  const Vertex_Table * vtab(s_model->GetVertexTable());
  Flavour_Vector fslep;
  for (size_t i(2);i<fv.size();++i)
    if (!fv[i].Strong()) fslep.push_back(fv[i]);
  for (Vertex_Table::const_iterator it(vtab->begin());it!=vtab->end();++it) {
    if (it->first.IsOn()      && !it->first.Strong() &&
        it->first.IsMassive() && !it->first.IsDummy()) {
      for (size_t i(0);i<it->second.size();++i) {
        bool on(true);
        double m(it->first.Mass());
        Single_Vertex * v(it->second[i]);
        for (size_t j(1);j<v->nleg;++j) {
          if (!v->on || v->dec)        { on=false; break; }
          if (v->in[j]==v->in[0])      { on=false; break; }
          if (v->in[j].IsDummy())      { on=false; break; }
          if ((m-=v->in[j].Mass())<0.) { on=false; break; }
          bool flavfound(false);
          for (size_t k(0);k<fslep.size();++k)
            if (v->in[j]==fslep[k])    { flavfound=true; break; }
          if (!flavfound)              { on=false; break; }
        }
        if (on) vlist.push_back(v);
      }
    }
  }
}

bool Signal_Process_FS_QED_Correction::FindResonances
(Particle_Vector& pv,std::vector<Particle_Vector>& rpvs,Flavour_Vector& rfl,
 const Vertex_List& vlist)
{
  if (vlist.empty()) return false;
  DEBUG_FUNC("find resonances in "<<pv.size()<<" particles");
  // find a combination in pv for which a vertex exists such that the
  // IS flavour is on-shell within m_resdist times its width
  // book-keep first to later disentangle competing resonances
  // need to book-keep i,j,k,abs(mij-mk)/wk
  std::map<double,std::vector<size_t> > restab;
  for (size_t i(0);i<pv.size();++i) {
    for (size_t j(i+1);j<pv.size();++j) {
      for (size_t k(0);k<vlist.size();++k) {
        double mdist(abs((pv[i]->Momentum()+pv[j]->Momentum()).Mass()
                         -vlist[k]->in[0].Mass())/vlist[k]->in[0].Width());
        if (vlist[k]->nleg==3 &&
            ((pv[i]->Flav()==vlist[k]->in[1] &&
              pv[j]->Flav()==vlist[k]->in[2]) ||
             (pv[i]->Flav()==vlist[k]->in[2] &&
              pv[j]->Flav()==vlist[k]->in[1])) &&
            mdist<m_resdist) {
          size_t ida[3]={i,j,k};
          restab[mdist]=std::vector<size_t>(ida,ida+3);
        }
      }
    }
  }
  if (restab.empty()) {
    msg_Debugging()<<"no resonances found"<<std::endl;
    return false;
  }
  if (msg_LevelIsDebugging()) {
    msg_Debugging()<<"resonances found:\n";
    for (std::map<double,std::vector<size_t> >::const_iterator
         it=restab.begin();it!=restab.end();++it)
      msg_Debugging()<<it->second[0]<<it->second[1]<<it->second[2]<<": "
                     <<vlist[it->second[2]]->in[0]<<" -> "
                     <<vlist[it->second[2]]->in[1]<<" "
                     <<vlist[it->second[2]]->in[2]
                     <<", |m-M|/W="<<it->first<<std::endl;
  }
  Particle_Vector usedparts;
  for (std::map<double,std::vector<size_t> >::const_iterator it=restab.begin();
       it!=restab.end();++it) {
    bool valid(true);
    for (size_t i(0);i<usedparts.size();++i)
      if (pv[it->second[0]]==usedparts[i] ||
          pv[it->second[1]]==usedparts[i]) { valid=false; break; }
    if (!valid) continue;
    usedparts.push_back(pv[it->second[0]]);
    usedparts.push_back(pv[it->second[1]]);
    msg_Debugging()<<"constructing decay: "<<vlist[it->second[2]]->in[0]<<" -> "
                                           <<vlist[it->second[2]]->in[1]<<" "
                                           <<vlist[it->second[2]]->in[2]<<"\n";
    rfl.push_back(vlist[it->second[2]]->in[0]);
    Particle_Vector parts;
    parts.push_back(pv[it->second[0]]);
    parts.push_back(pv[it->second[1]]);
    rpvs.push_back(parts);
  }
  for (Particle_Vector::iterator it=usedparts.begin();it!=usedparts.end();++it)
    for (Particle_Vector::iterator pit=pv.begin();pit!=pv.end();++pit)
      if (*it==*pit) { pv.erase(pit); break; }
  return true;
}

Blob_Vector Signal_Process_FS_QED_Correction::BuildResonantBlobs
(Particle_Vector& pv)
{
  DEBUG_FUNC("");
  // get production subprocesses for the active process
  std::string name(p_mehandler->Process()->Name());
  SubInfoVector siv(m_proc_lep_map[name]);
  // create blobs accordingly (only if lepton list is unambiguous)
  Blob_Vector blobs;
  msg_Debugging()<<siv.size()<<" subprocess infos for process "
                 <<name<<std::endl;
  msg_Debugging()<<"Particle content unambiguous? "
                 <<(ContainsNoAmbiguities(pv)?"yes":"no")<<std::endl;
  if (siv.size() && ContainsNoAmbiguities(pv)) {
    for (size_t i=0;i<siv.size();++i) {
      blobs.push_back(new Blob(Vec4D(0.,0.,0.,0.)));
      FillBlob(blobs[i],*siv[i],pv);
      msg_Debugging()<<"built decay blob for subprocess:"<<endl;
      msg_Debugging()<<*blobs[i]<<endl;
    }
  }
  // find/reconstruct possible resonances in the final state
  std::vector<Particle_Vector> rpvs;
  std::vector<Flavour>         rfl;
  const Vertex_List& vlist(m_proc_restab_map[name]);
  if (m_findresonances && pv.size()>1 && FindResonances(pv,rpvs,rfl,vlist)) {
    for (size_t i=0;i<rpvs.size();++i) {
      blobs.push_back(new Blob(Vec4D(0.,0.,0.,0.)));
      FillBlob(*blobs.rbegin(),rfl[i],rpvs[i]);
      msg_Debugging()<<"built blob for identified resonance:"<<endl;
      msg_Debugging()<<**blobs.rbegin()<<endl;
    }
  }
  // otherwise create global resonant blob
  // if there are leptons not contained in defined resonant blobs
  if (pv.size()) {
    blobs.push_back(new Blob(Vec4D(0.,0.,0.,0.)));
    FillBlob(*blobs.rbegin(),DetermineGenericResonance(pv),pv);
    msg_Debugging()<<"built generic blob:"<<endl;
    msg_Debugging()<<**blobs.rbegin()<<endl;
  }
  return blobs;
}

bool Signal_Process_FS_QED_Correction::ContainsNoAmbiguities
(const Particle_Vector& pv)
{
  // look whether any particle occurs more than once
  std::set<std::string> checklist;
  for (size_t i=0;i<pv.size();++i) {
    if (checklist.find(pv[i]->Flav().IDName()) == checklist.end())
      checklist.insert(pv[i]->Flav().IDName());
    else return false;
  }
  return true;
}

void Signal_Process_FS_QED_Correction::FillBlob
(Blob * blob, const Subprocess_Info& spi, Particle_Vector& pv)
{
  // find the leptons owned by this subprocess
  Particle_Vector localpv;
  bool onlyleptons(true);
  for (size_t i=0;i<spi.m_ps.size();++i) {
    if (spi.m_ps[i].m_fl.Strong()) onlyleptons=false;
    for (Particle_Vector::iterator it=pv.begin();it!=pv.end();) {
      if ((*it)->Flav()==spi.m_ps[i].m_fl) {
        localpv.push_back(*it);
        pv.erase(it);
      }
      else ++it;
    }
  }
  if (onlyleptons) FillBlob(blob,spi.m_fl,localpv);
  else FillBlob(blob,DetermineGenericResonance(localpv),localpv);
}

void Signal_Process_FS_QED_Correction::FillBlob
(Blob * blob, const Flavour& resflav, Particle_Vector& pv)
{
  Vec4D sum(MomentumSum(pv));
  for (Particle_Vector::iterator it=pv.begin();it!=pv.end();) {
    blob->AddToOutParticles(*it);
    pv.erase(it);
  }
  blob->AddToInParticles(new Particle(-1,resflav,sum,'R'));
  blob->InParticle(0)->SetFinalMass(blob->InParticle(0)->Momentum().Mass());
}

void Signal_Process_FS_QED_Correction::CleanUp(const size_t & mode) {}

void Signal_Process_FS_QED_Correction::Finish(const std::string &) {}

