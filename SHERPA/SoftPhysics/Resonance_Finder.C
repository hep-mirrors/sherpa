#include "SHERPA/SoftPhysics/Resonance_Finder.H"

#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Phys/Blob.H"
#include "ATOOLS/Phys/Blob_List.H"
#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Phys/Particle.H"

#include "MODEL/Main/Model_Base.H"
#include "MODEL/Main/Single_Vertex.H"

#include "PHASIC++/Process/Process_Base.H"
#include "PHASIC++/Process/Process_Info.H"
#include "PHASIC++/Process/Subprocess_Info.H"

#include "SHERPA/PerturbativePhysics/Matrix_Element_Handler.H"

using namespace SHERPA;
using namespace ATOOLS;
using namespace PHASIC;
using namespace MODEL;
using namespace std;

Resonance_Finder::Resonance_Finder(ATOOLS::Data_Reader* reader,
                                   Matrix_Element_Handler * meh) :
  m_on(true), m_resdist(1.), p_mehandler(meh)
{
  m_on = (reader->GetValue<std::string>("ME_QED_CLUSTERING","On")=="On");
  m_resdist = reader->GetValue<double>("ME_QED_CLUSTERING_THRESHOLD",1.);

  if (!m_on) return;
  ScanModelForEWResonances();
  InitialiseHelperParticles();
  IdentifyEWSubprocesses();
}

Resonance_Finder::Resonance_Finder() :
  m_on(false), m_resdist(1.)
{
}

Resonance_Finder::~Resonance_Finder()
{

}

void Resonance_Finder::ScanModelForEWResonances()
{
  Process_Vector pvec(p_mehandler->AllProcesses());
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
}

void Resonance_Finder::FindProcessPossibleResonances
(const Flavour_Vector& fv, MODEL::Vertex_List& vlist)
{
  const Vertex_Table *vtab(s_model->VertexTable());
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
        for (size_t j(1);j<v->in.size();++j) {
          if (v->dec)        { on=false; break; }
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

void Resonance_Finder::InitialiseHelperParticles()
{
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

void Resonance_Finder::IdentifyEWSubprocesses()
{
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
}

void Resonance_Finder::FindSubProcessInfosContainingLeptons
(const Process_Info& pi, SubInfoVector& siv)
{
  // loop over FS of process -> find subprocs, that are subsequent decays
  for (size_t i=0;i<pi.m_fi.m_ps.size();++i) {
    if (pi.m_fi.m_ps[i].m_ps.size()>1)
      FindSubProcessInfosContainingLeptons(pi.m_fi.m_ps[i],siv);
  }
}

void Resonance_Finder::FindSubProcessInfosContainingLeptons
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

void Resonance_Finder::BuildResonantBlobs
(Particle_Vector& pv, Blob_Vector& blobs)
{
  DEBUG_FUNC("");
  // get production subprocesses for the active process
  std::string name(p_mehandler->Process()->Name());
  SubInfoVector siv(m_proc_lep_map[name]);
  // create blobs accordingly (only if lepton list is unambiguous)
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
  if (m_on && pv.size()>1 && FindResonances(pv,rpvs,rfl,vlist)) {
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
}

bool Resonance_Finder::ContainsNoAmbiguities
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

void Resonance_Finder::FillBlob
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

void Resonance_Finder::FillBlob
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

bool Resonance_Finder::FindResonances
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
        if (vlist[k]->in.size()==3 &&
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
                     <<", |m-M|/W="<<it->first
                     <<" (max: "<<m_resdist<<")"<<std::endl;
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

Flavour Resonance_Finder::DetermineGenericResonance
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

Vec4D Resonance_Finder::MomentumSum
(const Particle_Vector& partvec)
{
  Vec4D sum(0.,0.,0.,0.);
  for (size_t i=0;i<partvec.size();++i) {
    sum += partvec[i]->Momentum();
  }
  return sum;
}

