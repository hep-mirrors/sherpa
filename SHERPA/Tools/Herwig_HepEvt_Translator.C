#include "SHERPA/Tools/Herwig_HepEvt_Translator.H"
#include "SHERPA/Tools/HepEvt_Interface.H"
#include "ATOOLS/Phys/Blob_List.H"

using namespace SHERPA;
using namespace ATOOLS;
using namespace std;

Herwig_HepEvt_Translator::Herwig_HepEvt_Translator(HepEvt_Interface * interface) :
  p_interface(interface)
{}

void Herwig_HepEvt_Translator::HepEvt2Particle(const int pos)
{
  /*
    if (abs(p_idhep[pos])==9902210) return;

    std::cout<<pos<<": stat,id "<<p_isthep[pos]<<","<<p_idhep[pos]
    <<"; mos : "<<p_jmohep[2*pos]<<","<<p_jmohep[2*pos+1]
    <<"; das : "<<p_jdahep[2*pos]<<","<<p_jdahep[2*pos+1]<<std::endl
    <<"    mom: "<<p_phep[5*pos+3]<<" "<<p_phep[5*pos+0]<<" "<<p_phep[5*pos+1]
    <<" "<<p_phep[5*pos+2]<<" "<<p_phep[5*pos+4]<<std::endl;

    Flavour flav;
    if ((m_generator==gtp::Herwig) &&
    (p_idhep[pos]==94 || p_idhep[pos]==0)) flav=Flavour(kf_none);
    else flav.FromHepEvt(p_idhep[pos]);
    Vec4D momentum     = Vec4D(p_phep[3+pos*5],p_phep[0+pos*5],p_phep[1+pos*5],p_phep[2+pos*5]);
    Particle * newpart = new Particle(pos+1,flav,momentum);
    newpart->SetFinalMass(p_phep[4+pos*5]);
    newpart->SetStatus(part_status::code(p_isthep[pos]));
    m_convertH2S[pos]=std::pair<Particle*,bool>(newpart,true);
    std::cout<<"   "<<(*newpart)<<std::endl;
  */
}

bool Herwig_HepEvt_Translator::ConstructBlobs(ATOOLS::Blob_List * const blobs)
{
  /*
    bool signalwarner=true,breakit=false,first_fsr=false;
    int helper, helper1;
    ATOOLS::Particle * part, * mother;
    ATOOLS::Blob * blob, * help, *fsr,
    * signal = new ATOOLS::Blob(), 
    * beam1  = new ATOOLS::Blob(), 
    * beam2  = new ATOOLS::Blob(), 
    * isr1   = new ATOOLS::Blob(), 
    * isr2   = new ATOOLS::Blob(), 
    * fsr1   = new ATOOLS::Blob(), 
    * fsr2   = new ATOOLS::Blob(),
    * cf     = new ATOOLS::Blob();
    signal->SetType(btp::Signal_Process);
    beam1->SetType(btp::Beam);
    beam2->SetType(btp::Beam);
    isr1->SetType(btp::IS_Shower);
    isr2->SetType(btp::IS_Shower);
    fsr1->SetType(btp::FS_Shower);
    fsr2->SetType(btp::FS_Shower);
    cf->SetType(btp::Cluster_Formation);

    blobs->push_back(beam1);  beam1->SetId();
    blobs->push_back(beam2);  beam2->SetId();
    blobs->push_back(isr1);   isr1->SetId();
    blobs->push_back(isr2);   isr2->SetId();
    blobs->push_back(signal); signal->SetId();
    blobs->push_back(fsr1);   fsr1->SetId();
    blobs->push_back(fsr2);   fsr2->SetId();
    blobs->push_back(cf);     cf->SetId();

    Translation_Map::iterator piter, miter;
    for (int i=0;i<m_nhep;i++) {
    piter = m_convertH2S.find(i);
    if (piter==m_convertH2S.end() || !piter->second.second) continue;
    part = piter->second.first;
    if (part->Status()!=part_status::active && 
    part->Status()!=part_status::decayed) 
    part->SetStatus(part_status::documentation);
    switch (p_isthep[i]) {
    case 101:
    beam1->AddToInParticles(part);
    part->SetStatus(part_status::decayed);
    piter->second.second = false;
    break;
    case 102:
    beam2->AddToInParticles(part);
    part->SetStatus(part_status::decayed);
    piter->second.second = false;
    break;
    case 121:
    isr1->AddToOutParticles(part);
    signal->AddToInParticles(part);
    part->SetStatus(part_status::decayed);
    part->SetInfo('G');
    piter->second.second = false;
    break;
    case 122:
    isr2->AddToOutParticles(part);
    signal->AddToInParticles(part);
    part->SetStatus(part_status::decayed);
    part->SetInfo('G');
    piter->second.second = false;
    break;
    case 123:
    first_fsr = true;
    case 124:
    if (p_isthep[p_jmohep[2*i]-1]==125 || 
    p_isthep[p_jmohep[2*i]-1]==195 || p_isthep[p_jmohep[2*i]-1]==155) break;
    if (part->Flav()==Flavour(kf_none) && 
    p_idhep[p_jdahep[2*i]-1]==94 && p_isthep[p_jdahep[2*i]-1]==144) {
    helper = p_jdahep[2*(p_jdahep[2*(p_jdahep[2*i]-1)]-1)]-1;
    if (p_idhep[helper]==0 && p_isthep[helper]==155) {
    signalwarner=false;
    //std::cout<<"Gotcha : "<<p_jdahep[2*i]-1<<" : "<<p_idhep[p_jdahep[2*i]-1]
    //	   <<" -> "<<helper<<", "<<p_jdahep[2*helper]-1
    //	   <<" -> "<<p_jdahep[2*(p_jdahep[2*helper]-1)]-1<<std::endl;
    for (int k=0;k<2;k++) {
    helper1 = p_jdahep[2*helper+k]-1;
    piter = m_convertH2S.find(helper1);
    if (piter==m_convertH2S.end() || !piter->second.second) continue;
    part = piter->second.first;
    piter->second.second = false;
    part->SetStatus(part_status::decayed);
    part->SetInfo('H');
    signal->AddToOutParticles(part);
    if (k==0) { 
    //std::cout<<"Add "<<helper1+1<<" to "<<first_fsr<<std::endl;
    if (first_fsr) fsr1->AddToInParticles(part);
    else fsr2->AddToInParticles(part);
    }
    else {
    fsr = new ATOOLS::Blob();
    blobs->push_back(fsr);
    fsr->SetType(btp::FS_Shower);
    fsr->SetId();
    fsr->AddToInParticles(part);
    }
    if (p_jdahep[2*helper1]>p_jdahep[2*helper1+1]) {
    helper1 = p_jdahep[2*helper1]-1;
    piter   = m_convertH2S.find(helper1);
    if (piter==m_convertH2S.end() || !piter->second.second) continue;
    part = piter->second.first;
    piter->second.second = false;
    part->SetStatus(part_status::decayed);
    part->SetInfo('f');
    if (k==0) { 
    if (first_fsr) fsr1->AddToOutParticles(part);
    else fsr2->AddToOutParticles(part);
    }
    else fsr->AddToOutParticles(part);	    
    help = new ATOOLS::Blob();
    blobs->push_back(help);
    help->SetType(btp::Hard_Decay);
    help->SetId();
    help->AddToInParticles(part);
    for (int j=p_jdahep[2*helper1]-1;j<=p_jdahep[2*helper1+1]-1;j++) {
    //std::cout<<j+1<<" in "<<p_jdahep[2*helper1]<<"..."<<p_jdahep[2*helper1+1]<<std::endl;
    piter = m_convertH2S.find(j);
    if (piter==m_convertH2S.end() || !piter->second.second) continue;
    piter->second.second = false;
    part = piter->second.first;
    help->AddToOutParticles(part);
    part->SetInfo('h');
    if (p_isthep[j]==1) part->SetStatus(part_status::active);
    else {
    blob = new ATOOLS::Blob();
    blobs->push_back(blob);
    blob->SetType(btp::FS_Shower);
    blob->SetId();
    blob->AddToInParticles(part);
    part->SetStatus(part_status::decayed);
    //std::cout<<"---------------"<<j+1<<" -> "<<p_jdahep[2*j]<<" ->"
    //	   <<p_jdahep[2*(p_jdahep[2*j]-1)]<<","<<p_jdahep[2*(p_jdahep[2*j]-1)+1]
    //	   <<"---------------------------"<<std::endl;
    for (int l=p_jdahep[2*(p_jdahep[2*j]-1)];l<=p_jdahep[2*(p_jdahep[2*j]-1)+1];l++) {
    //std::cout<<"Find this : "<<p_jdahep[2*j]<<" -> "<<l<<std::endl;
    breakit=true;
    piter = m_convertH2S.find(l-1);
    if (piter==m_convertH2S.end() || !piter->second.second) continue;
    part = piter->second.first;
    piter->second.second = false;
    blob->AddToOutParticles(part);
    cf->AddToInParticles(part);
    part->SetStatus(part_status::decayed);
    part->SetInfo('f');
    }
    }
    }
    }
    //if (k==0) std::cout<<(*signal)<<(*fsr2)<<(*help)<<std::endl;
    //if (k==1) std::cout<<(*signal)<<(*fsr)<<(*help)<<std::endl;
    }
    //if (breakit) {
    // std::cout<<(*blobs)<<std::endl;
    // abort();
    //}
    first_fsr = false;
    break;
    }
    else{
    msg_Error()<<"Error in HepEvt_Interface::ConstructBlobsFromHerwig : "<<std::endl
    <<"   Unexpected feature in HepEvt, will continue & hope for the best."<<std::endl;
    first_fsr = false;
    break;
    }
    }
    //std::cout<<"Check124 "<<i+1<<" -> "<<p_jdahep[2*i]<<" ("<<p_isthep[p_jdahep[2*i]-1]<<")"
    //	       <<" <- "<<p_jmohep[2*i]<<" ("<<p_isthep[p_jmohep[2*i]-1]<<")"<<std::endl;
    signal->AddToOutParticles(part);
    if (first_fsr) fsr1->AddToInParticles(part);
    else fsr2->AddToInParticles(part);
    part->SetStatus(part_status::decayed);
    piter->second.second = false;
    part->SetInfo('H');
    if (p_isthep[p_jdahep[2*i]-1]==195 || p_isthep[p_jdahep[2*i]-1]==3) {
    //helper = i;
    helper = p_jdahep[2*i]-1;
    //std::cout<<"Try 1 : "<<i+1<<"-> "<<helper+1<<std::endl;
    if (p_isthep[p_jdahep[2*i]-1]==3) helper = (p_jdahep[2*helper]-1);
    //std::cout<<"Try 2 : "<<i+1<<"-> "<<helper+1<<std::endl;
    piter = m_convertH2S.find(helper);
    if (piter==m_convertH2S.end() || !piter->second.second) continue;
    piter->second.second = false;
    part = piter->second.first;
    part->SetStatus(part_status::decayed);
    part->SetInfo('f');
    if (first_fsr) fsr1->AddToOutParticles(part);
    else fsr2->AddToOutParticles(part);
    help = new ATOOLS::Blob();
    blobs->push_back(help);
    help->SetType(btp::Hard_Decay);
    help->SetId();
    help->AddToInParticles(part);
    for (int j=p_jdahep[2*helper]-1;j<=p_jdahep[2*helper+1]-1;j++) {
    //std::cout<<helper+1<<" -> "<<j+1<<" into decay blob."<<std::endl;
    piter = m_convertH2S.find(j);
    if (piter==m_convertH2S.end() || !piter->second.second) continue;
    part = piter->second.first;
    help->AddToOutParticles(part);
    piter->second.second = false;
    part->SetInfo('h');
    if (p_isthep[j]==1) part->SetStatus(part_status::active);
    else {
    blob = new ATOOLS::Blob();
    blobs->push_back(blob);
    blob->SetType(btp::FS_Shower);
    blob->SetId();
    blob->AddToInParticles(part);
    part->SetStatus(part_status::decayed);
    //std::cout<<"---------------"<<j+1<<" -> "<<p_jdahep[2*j]<<" ->"
    //	     <<p_jdahep[2*(p_jdahep[2*j]-1)]<<","<<p_jdahep[2*(p_jdahep[2*j]-1)+1]
    //	     <<"---------------------------"<<std::endl;
    for (int l=p_jdahep[2*(p_jdahep[2*j]-1)];l<=p_jdahep[2*(p_jdahep[2*j]-1)+1];l++) {
    //std::cout<<"Find this : "<<p_jdahep[2*j]<<" -> "<<l<<std::endl;
    breakit=true;
    piter = m_convertH2S.find(l-1);
    if (piter==m_convertH2S.end() || !piter->second.second) continue;
    part = piter->second.first;
    piter->second.second = false;
    blob->AddToOutParticles(part);
    cf->AddToInParticles(part);
    part->SetStatus(part_status::decayed);
    part->SetInfo('f');
    }
    }
    }
    //std::cout<<(*blobs)<<std::endl<<std::endl<<std::endl<<std::endl;
    //abort();
    }
    first_fsr = false;
    break;
    case 125:
    case 155:
    break;
    case 141:
    beam1->AddToOutParticles(part);
    isr1->AddToInParticles(part);
    part->SetStatus(part_status::decayed);
    part->SetInfo('I');
    piter->second.second = false;
    break;
    case 142:
    beam2->AddToOutParticles(part);
    isr2->AddToInParticles(part);
    part->SetStatus(part_status::decayed);
    part->SetInfo('I');
    piter->second.second = false;
    break;
    case 143:
    // Outgoing jet 
    break;
    case 144:
    // Outgoing jet 
    break;
    case 160:
    break;
    case 158:
    case 159:
    case 161:
    case 162:
    break;
    case 181:
    case 182:
    case 183:
    case 184:
    case 185:
    case 186:
    part->SetStatus(part_status::decayed);
    part->SetInfo('C');
    cf->AddToOutParticles(part);
    if (p_jdahep[2*i]!=0) {
    blob = new ATOOLS::Blob();
    blobs->push_back(blob);
    blob->SetType(btp::Cluster_Decay);
    blob->SetId();
    blob->AddToInParticles(part);
    }
    piter->second.second = false;
    break;
    case 195:
    case 196:
    case 197:
    part->SetStatus(part_status::active);
    part->SetInfo('P');
    miter = m_convertH2S.find(p_jmohep[2*i]-1);
    if (miter==m_convertH2S.end()) continue;
    mother = miter->second.first;
    blob   = mother->DecayBlob();
    if (!blob) {
    blob = new ATOOLS::Blob();
    blobs->push_back(blob);
    blob->SetType(btp::Cluster_Decay);
    blob->SetId();
    blob->AddToInParticles(mother);
    }
    blob->AddToOutParticles(part);
    piter->second.second = false;
    blob->SetPosition(Vec4D(p_vhep[4*i+3],p_vhep[4*i+0],p_vhep[4*i+1],p_vhep[4*i+2]));
    if (p_jdahep[2*i]!=0) {
    part->SetStatus(part_status::decayed);
    part->SetInfo('p');
    blob = new ATOOLS::Blob();
    blobs->push_back(blob);
    blob->SetType(btp::Hadron_Decay);
    blob->SetId();
    blob->AddToInParticles(part);
    }
    break;
    case 198:
    part->SetStatus(part_status::active);
    part->SetInfo('D');
    miter = m_convertH2S.find(p_jmohep[2*i]-1);
    if (miter==m_convertH2S.end()) continue;
    mother = miter->second.first;
    blob   = mother->DecayBlob();
    if (!blob) {
    blob = new ATOOLS::Blob();
    blobs->push_back(blob);
    blob->SetType(btp::Hadron_Decay);
    blob->SetId();
    blob->AddToInParticles(mother);
    }
    blob->AddToOutParticles(part);
    piter->second.second = false;
    blob->SetPosition(Vec4D(p_vhep[4*i+3],p_vhep[4*i+0],p_vhep[4*i+1],p_vhep[4*i+2]));
    if (p_jdahep[2*i]!=0) {
    part->SetStatus(part_status::decayed);
    part->SetInfo('d');
    blob = new ATOOLS::Blob();
    blobs->push_back(blob);
    blob->SetType(btp::Hadron_Decay);
    blob->SetId();
    blob->AddToInParticles(part);
    }
    break;
    case 199:
    part->SetStatus(part_status::active);
    part->SetInfo('D');
    miter = m_convertH2S.find(p_jmohep[2*i]-1);
    if (miter==m_convertH2S.end()) continue;
    mother = miter->second.first;
    blob   = mother->DecayBlob();
    if (!blob) {
    blob = new ATOOLS::Blob();
    blobs->push_back(blob);
    blob->SetType(btp::Hadron_Decay);
    blob->SetId();
    blob->AddToInParticles(mother);
    }
    blob->AddToOutParticles(part);
    piter->second.second = false;
    blob->SetPosition(Vec4D(p_vhep[4*i+3],p_vhep[4*i+0],p_vhep[4*i+1],p_vhep[4*i+2]));
    if (p_jdahep[2*i]!=0) {
    part->SetStatus(part_status::decayed);
    part->SetInfo('d');
    blob = new ATOOLS::Blob();
    blobs->push_back(blob);
    blob->SetType(btp::Hadron_To_Parton);
    blob->SetId();
    blob->AddToInParticles(part);
    help = new ATOOLS::Blob();
    blobs->push_back(help);
    help->SetType(btp::Cluster_Formation);
    help->SetId();
    FollowDaughters(i,blob,help);
    }
    else {
    msg_Error()<<"Error in HepEvt_Interface::ConstructBlobsFromHerwig."<<std::endl
    <<"   Mother of heavy hadron has no decay blob yet."<<std::endl
    <<"   Will abort."<<Particle::Counter()<<" / "<<Blob::Counter()<<std::endl;
    }
    break;      
    case 200:
    part->SetStatus(part_status::active);
    part->SetInfo('D');
    miter = m_convertH2S.find(p_jmohep[2*i]-1);
    if (miter==m_convertH2S.end()) continue;
    mother = miter->second.first;
    blob   = mother->DecayBlob();
    if (!blob) {
    msg_Error()<<"Error in HepEvt_Interface::ConstructBlobsFromHerwig."<<std::endl
    <<"   Mother of heavy hadron flavour has no decay blob yet."<<std::endl
    <<"   Will create a new blob and hope for the best."
    <<Particle::Counter()<<" / "<<Blob::Counter()<<std::endl;
    blob = new ATOOLS::Blob();
    blobs->push_back(blob);
    blob->SetType(btp::Hadron_Decay);
    blob->SetId();
    blob->AddToInParticles(mother);
    }
    blob->AddToOutParticles(part);
    piter->second.second = false;
    blob->SetPosition(Vec4D(p_vhep[4*i+3],p_vhep[4*i+0],p_vhep[4*i+1],p_vhep[4*i+2]));
    if (p_jdahep[2*i]!=0) {
    part->SetStatus(part_status::decayed);
    part->SetInfo('d');
    blob = new ATOOLS::Blob();
    blobs->push_back(blob);
    blob->SetType(btp::Hadron_Mixing);
    blob->SetId();
    blob->AddToInParticles(part);
    }
    else {
    msg_Error()<<"Error in HepEvt_Interface::ConstructBlobsFromHerwig."<<std::endl
    <<"   Mother of heavy hadron has no decay blob yet."<<std::endl
    <<"   Will abort."<<Particle::Counter()<<" / "<<Blob::Counter()<<std::endl;
    }
    break;      
    case 1:
    if (!piter->second.second) continue;
    miter = m_convertH2S.find(p_jmohep[2*i]-1);
    if (miter==m_convertH2S.end()) continue;
    mother = miter->second.first;
    blob   = mother->DecayBlob();
    if (blob) {
    mother->SetStatus(part_status::decayed);
    if (mother->Info()=='P') mother->SetInfo('p');
    else if (mother->Info()=='D') mother->SetInfo('d');
    blob->AddToOutParticles(part);
    piter->second.second = false;
    blob->SetPosition(Vec4D(p_vhep[4*i+3],p_vhep[4*i+0],p_vhep[4*i+1],p_vhep[4*i+2]));
    }
    part->SetStatus(part_status::active);
    part->SetInfo('D');
    break;
    case 2:
    if (p_isthep[p_jmohep[2*i]-1]==141) isr1->AddToOutParticles(part);
    if (p_isthep[p_jmohep[2*i]-1]==142) isr2->AddToOutParticles(part);
    if (p_isthep[p_jmohep[2*i]-1]==143) fsr1->AddToOutParticles(part);
    if (p_isthep[p_jmohep[2*i]-1]==144) fsr2->AddToOutParticles(part);
    cf->AddToInParticles(part); 
    part->SetStatus(part_status::decayed); 
    piter->second.second = false; 
    default : break;
    }
    }
    if (signalwarner) {
    if (signal->NOutP()!=2 || signal->NInP()<2) {
    msg_Error()<<"Error in HepEvt_Interface::ConstructBlobsFromHerwig"<<std::endl
    <<"   Signal is funny: "<<signal->NInP()<<" -> "<<signal->NOutP()<<std::endl
    <<"   ====================================="<<std::endl
    <<(*signal)<<std::endl
    <<"   ====================================="<<std::endl
    <<"   Clear blobs, return false and hope for the best."
    <<Particle::Counter()<<" / "<<Blob::Counter()<<std::endl
    <<(*blobs)<<std::endl;

    if (!blobs->empty()) {
    for (Blob_List::iterator blit=blobs->begin();blit!=blobs->end();++blit) delete (*blit);
    blobs->clear();
    }
    DeleteObsolete(1);
    return false;
    }
    }
    DeleteObsolete(-1);
    return true;
    }

    void Herwig_HepEvt_Translator::FollowDaughters(const int & i,ATOOLS::Blob * hadron,ATOOLS::Blob * clusters) 
    {
    Translation_Map::iterator piter;
    Particle * part;
    int end = Max(p_jdahep[2*i]-1,p_jdahep[2*i+1]-1);
    for (int j=p_jdahep[2*i]-1;j<=end;j++) {
    if (j>m_nhep || j<0) {
    msg_Error()<<"Error in HepEvt_Interface::FollowDaughters("<<i<<")"<<std::endl
    <<"   counter "<<j<<" out of bounds ("<<m_nhep<<")."<<std::endl
    <<"   Will abort the run."<<std::endl;
    abort();
    }
    if (p_isthep[j]==183) {
    // Final cluster
    piter = m_convertH2S.find(j);
    if (piter==m_convertH2S.end()) {
    msg_Error()<<"Error in HepEvt_Interface::FollowDaughters."<<std::endl
    <<"   Mix up of secondary clusters."<<std::endl
    <<"   Will abort."<<std::endl;
    abort();
    }
    if (piter->second.second) {
    piter->second.second = false;
    part = piter->second.first;
    part->SetStatus(part_status::decayed); 
    part->SetInfo('C');
    clusters->AddToOutParticles(part);
    } 
    piter = m_convertH2S.find(i);
    if (piter==m_convertH2S.end()) {
    msg_Error()<<"Error in HepEvt_Interface::FollowDaughters."<<std::endl
    <<"   Mix up of secondary clusters."<<std::endl
    <<"   Will abort."<<std::endl;
    abort();
    }
    if (piter->second.second) {
    piter->second.second = false;
    part = piter->second.first;
    part->SetStatus(part_status::decayed); 
    part->SetInfo('f');
    hadron->AddToOutParticles(part);
    clusters->AddToInParticles(part);
    } 
    }
    else if (p_isthep[j]==1) {
    // Final lepton
    piter = m_convertH2S.find(j);
    if (piter==m_convertH2S.end()) {
    msg_Error()<<"Error in HepEvt_Interface::FollowDaughters."<<std::endl
    <<"   Mix up of secondary clusters."<<std::endl
    <<"   Will abort."<<std::endl;
    abort();
    }
    if (piter->second.second) {
    piter->second.second = false;
    part = piter->second.first;
    part->SetStatus(part_status::active); 
    part->SetInfo('F');
    hadron->AddToOutParticles(part);
    } 
    } 
    else if (p_isthep[j]==160) {
    // Final spectator
    piter = m_convertH2S.find(j);
    if (piter==m_convertH2S.end()) {
    msg_Error()<<"Error in HepEvt_Interface::FollowDaughters."<<std::endl
    <<"   Mix up of secondary clusters."<<std::endl
    <<"   Will abort."<<std::endl;
    abort();
    }
    if (piter->second.second) {
    piter->second.second = false;
    part = piter->second.first;
    part->SetStatus(part_status::decayed); 
    part->SetInfo('S');
    hadron->AddToOutParticles(part);
    clusters->AddToInParticles(part);
    } 
    } 
    else {
    piter = m_convertH2S.find(j);
    if (piter==m_convertH2S.end()) {
    msg_Error()<<"Error in HepEvt_Interface::FollowDaughters."<<std::endl
    <<"   Mix up of secondary clusters."<<std::endl
    <<"   Will abort."<<std::endl;
    abort();
    }
    p_isthep[j]=155;
    FollowDaughters(j,hadron,clusters);
    }
    }
    }

    bool Herwig_HepEvt_Translator::IdentifyBlobs(ATOOLS::Blob_List * const blobs) 
    {
    int                       counter;
    bool                      test;
    Blob_List::iterator       biter, biter2;
    Particle                * search, * dummy;
    Blob                    * prod, * meps, * dec;
    std::vector<Blob *>       obsoletes;
    std::vector<Particle *>   incomings;
    incomings.clear();

    //std::cout<<(*blobs)<<std::endl<<"---------------------------------------------------------"<<std::endl;

    // Beams and bunches
    for (biter=blobs->begin();biter!=blobs->end();biter++) {
    if ((*biter)->Type()!=btp::Unspecified) continue;
    if ((*biter)->NInP()==0 && (*biter)->NOutP()==1) {
    ATOOLS::Particle *incoming=(*biter)->OutParticle(0);
    (*biter)->RemoveOutParticle(incoming);
    delete *biter;
    blobs->erase(biter);
    ATOOLS::Blob * beam=incoming->DecayBlob();
    if (incoming->Flav().IsHadron()) {
    if (beam->OutParticle(0)->DecayBlob()->OutParticle(0)->Flav().IsHadron()) {
    beam->SetType(btp::Bunch);
    beam->OutParticle(0)->DecayBlob()->SetType(btp::Beam);
    }
    else {
    beam->SetType(btp::Beam);
    }
    }
    else {
    kf_code in=incoming->Flav().Kfcode();
    for (int i=0;i<beam->NOutP();++i) {
    kf_code out=beam->OutParticle(i)->Flav().Kfcode();
    if (in==kf_e && out==kf_photon) {
    beam->SetType(btp::Bunch);
    beam->OutParticle(0)->DecayBlob()->SetType(btp::Beam);
    }

    //if (in==kf_e && out==kf_e) {
    //  beam->SetType(btp::Bunch);
    //  beam->OutParticle(0)->DecayBlob()->SetType(btp::Beam);
    //}
    }
    }
    }
    }
    //std::cout<<"Beams and bunches"<<std::endl;

    // (IS) Shower Blobs
    counter = 0;
    for (biter=blobs->begin();biter!=blobs->end();biter++) {
    if ((*biter)->NInP()==1 && (*biter)->Type()==btp::Unspecified) {
    search = (*biter)->InParticle(0);
    if (search->ProductionBlob()->Type()==btp::Beam) {
    (*biter)->SetType(btp::IS_Shower); 
    counter++; 
    }
    }
    if (counter==2) { incomings.clear(); break; }
    }
    //std::cout<<"IS"<<std::endl;

  
    // ME Blob
    for (biter=blobs->begin();biter!=blobs->end();biter++) {
    if ((*biter)->NInP()==2 && (*biter)->Type()==btp::Unspecified) {
    //std::cout<<"Check for signal : "<<std::endl<<((**biter))<<std::endl;
    test = false;
    if (incomings.size()==2) {
    if (((*biter)->InParticle(0)==incomings[0] && 
    (*biter)->InParticle(1)==incomings[1] ) || 
    ((*biter)->InParticle(0)==incomings[1] && 
    (*biter)->InParticle(1)==incomings[0] ) ) {
    test = true;
    incomings.clear();
    }
    //std::cout<<"Unspecified 2 incomings: Test = true "<<test<<std::endl;
    }
    else {
    counter = 0;
    meps    = NULL;
    for (int i=0;i<2;i++) {
    search = (*biter)->InParticle(i);
    prod   = search->ProductionBlob();
    if (prod->NInP()==1 && prod->NOutP()==1 &&
    prod->Type()==btp::Unspecified) {
    if (prod->InParticle(0)->ProductionBlob()->Type()==btp::IS_Shower) {
    if (meps==NULL) {
    meps = prod;
    meps->SetType(btp::ME_PS_Interface_IS);
    counter++;
    }
    else {
    meps->AddToInParticles(prod->RemoveInParticle(0));
    meps->AddToOutParticles(prod->RemoveOutParticle(0));
    counter++;
    for (biter2=blobs->begin();biter2!=blobs->end();biter2++) {
    if ((*biter2)==prod) { 
    delete (*biter2);
    blobs->erase(biter2); 
    break; 
    }
    }
    }
    }
    }
    }
    if (counter==2) test = true;
    //std::cout<<"Unspecified mot 2 incomings: Test = true "<<test<<std::endl;
    }
    if (test) {
    (*biter)->SetType(btp::Signal_Process);
    meps     = NULL;
    for (int i=0;i<(*biter)->NOutP();i++) {
    search = (*biter)->OutParticle(i);
    dec    = search->DecayBlob();
    if (dec->NInP()==1 && dec->NOutP()==1 &&
    dec->Type()==btp::Unspecified) {
    dummy = dec->OutParticle(0); 
    if (dummy->DecayBlob()->Type()==btp::Unspecified) {
    if (meps==NULL) {
    meps = dec;
    meps->SetType(btp::ME_PS_Interface_FS);
    dummy->DecayBlob()->SetType(btp::FS_Shower);
    }
    else {
    meps->AddToInParticles(dec->RemoveInParticle(0));
    dec->RemoveOutParticle(0);
    meps->AddToOutParticles(dummy);
    dummy->DecayBlob()->SetType(btp::FS_Shower);
    for (biter2=blobs->begin();biter2!=blobs->end();biter2++) {
    if ((*biter2)==dec) { 
    delete (*biter2);
    blobs->erase(biter2); 
    break; 
    }
    }
    }
    }
    }
    else if (dec->NInP()==1 && dec->NOutP()>1 &&
    dec->Type()==btp::Unspecified) dec->SetType(btp::FS_Shower);
    }
    break;
    }
    }
    else if ((*biter)->Type()==btp::Signal_Process) {
    //std::cout<<"Is signal : "<<std::endl<<((**biter))<<std::endl;
    for (int i=0;i<(*biter)->NOutP();i++) {
    search = (*biter)->OutParticle(i);
    dec    = search->DecayBlob();
    if (dec->NInP()==1 && dec->Type()==btp::Unspecified) {
    if (dec->NOutP()!=1) 
    dec->SetType(btp::FS_Shower);
    else {
    meps = dec->OutParticle(0)->DecayBlob();
    if (meps==NULL) 
    dec->SetType(btp::FS_Shower);	    
    else if (meps->OutParticle(0)->Flav()==Flavour(kf_string))
    dec->SetType(btp::FS_Shower);	    
    }
    }
    }
    }
    }
    //std::cout<<"ME"<<std::endl;

    // Fragmentation blob
    Flavour cluster; cluster.FromHepEvt(91);
    Flavour string;  string.FromHepEvt(92);
    for (biter=blobs->begin();biter!=blobs->end();biter++) {
    if ((*biter)->NOutP()==1 && (*biter)->Type()==btp::Unspecified &&
    (*biter)->OutParticle(0)->Flav()==string) {
    search = (*biter)->OutParticle(0);
    for (biter2=blobs->begin();biter2!=blobs->end();biter2++) {
    if ((*biter2)->NInP()==1 &&
    (*biter2)->Type()==btp::Unspecified &&
    (*biter2)->InParticle(0)==search &&
    (*biter2)!=(*biter)) {
    (*biter)->RemoveOutParticle(0);
    (*biter2)->RemoveInParticle(0);
    delete search; search=NULL;
    for (int i=(*biter2)->NOutP()-1;i>=0;i--) {
    dummy = (*biter2)->RemoveOutParticle(i);
    (*biter)->AddToOutParticles(dummy);
    }
    (*biter)->SetType(btp::Fragmentation);
    delete (*biter2);
    blobs->erase(biter2);
    break;
    }
    }
    }
    }
    //std::cout<<"Frag"<<std::endl;

    // Hadron blobs
    for (biter=blobs->begin();biter!=blobs->end();biter++) {
    if ((*biter)->NInP()==1 && 
    (*biter)->Type()==btp::Unspecified &&
    (*biter)->InParticle(0)->Flav().IsHadron()) {
    (*biter)->SetType(btp::Hadron_Decay);
    }
    }
    //std::cout<<"Hadron"<<std::endl;

    Blob *nirwana = new Blob();
    nirwana->SetStatus(blob_status::inactive);
    // Nirwana particles
    for (biter=blobs->begin();biter!=blobs->end();biter++) {
    for (size_t i=0;i<(size_t)(*biter)->NOutP();++i) {
    Particle *part=(*biter)->OutParticle(i);
    if (part->Status()!=part_status::active && part->DecayBlob()==NULL) {
    nirwana->AddToInParticles(part);
    }
    }
    }
    nirwana->SetTypeSpec("Nirwana");
    blobs->push_back(nirwana);
    //std::cout<<"Nirwana"<<std::endl;


    int blobid = 0;
    for (biter=blobs->begin();biter!=blobs->end();biter++) {
    (*biter)->SetId(-blobid);
    blobid++;
    }

    return true;
    }

    void Herwig_HepEvt_Translator::DeleteObsolete(const int mode)
    {
    //std::cout<<"Delete them !! "<<mode<<std::endl;
    if (!m_convertH2S.empty()) {
    for (Translation_Map::iterator piter=m_convertH2S.begin();
    piter!=m_convertH2S.end();piter++) {
    switch(mode) {
    case 0:
    delete (piter->second.first); piter->second.first=NULL; 
    break;  
    case 1:
    if (!piter->second.first->ProductionBlob() &&
    !piter->second.first->DecayBlob() &&
    piter->second.first) {
    delete (piter->second.first); piter->second.first=NULL;
    } 
    break;  
    default:
    if (piter->second.second) {
    delete (piter->second.first); piter->second.first=NULL; 
    }
    break;
    }
    }
    m_convertH2S.clear();
    }
  */
}

