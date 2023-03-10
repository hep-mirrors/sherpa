#include "HADRONS++/Main/Hadron_Decay_Channel.H"
#include "HADRONS++/ME_Library/HD_ME_Base.H"
#include "HADRONS++/ME_Library/Generic.H"
#include "HADRONS++/ME_Library/Current_ME.H"
#include "HADRONS++/Current_Library/Current_Base.H"
#include "HADRONS++/PS_Library/HD_PS_Base.H"
#include "PHASIC++/Decays/Decay_Table.H"
#include "PHASIC++/Channels/Multi_Channel.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Math/Vector.H"
#include "ATOOLS/Org/Shell_Tools.H"
#include "ATOOLS/Phys/Blob.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/My_MPI.H"

using namespace HADRONS;
using namespace ATOOLS;
using namespace METOOLS;
using namespace std;

Hadron_Decay_Channel::Hadron_Decay_Channel(Flavour fl, const Mass_Selector* ms) :
  Decay_Channel(fl, ms),
  m_always_integrate(false),
  m_cp_asymmetry_C(0.0), m_cp_asymmetry_S(0.0)
{
}

Hadron_Decay_Channel::~Hadron_Decay_Channel()
{
}

void Hadron_Decay_Channel::Initialise(Scoped_Settings s, GeneralModel startmd)
{
  m_physicalflavours=m_flavours;
  for (size_t i=0; i<m_flavours.size(); ++i) {
    map<kf_code,kf_code>::const_iterator it = 
      Tools::aliases.find(m_flavours[i].Kfcode());
    if (it!=Tools::aliases.end())
      m_physicalflavours[i] = Flavour(it->second, m_flavours[i].IsAnti());
  }
  
  double totalmass=0.0;
  for (size_t i=1; i<m_flavours.size(); ++i) {
    totalmass+=m_flavours[i].HadMass();
  }
  if(totalmass>m_flavours[0].HadMass()) {
    msg_Error()<<"Error in "<<METHOD<<" for "<<Name()<<"\n"
	       <<"    Total outgoing mass heavier than incoming particle.\n"
	       <<"    Will return and hope for the best.\n";
    return;
  }
  SetChannels(new PHASIC::Multi_Channel(""));
  Channels()->SetNin(1);
  Channels()->SetNout(NOut());
  m_startmd=startmd;

  m_cp_asymmetry_S = s["CPAsymmetryS"].SetDefault(0.0).Get<double>();
  m_cp_asymmetry_C = s["CPAsymmetryC"].SetDefault(0.0).Get<double>();

  // convert C and S to lambda, assuming DeltaGamma=0 for the determination of C and S.
  // this allows DeltaGamma dependent terms in the asymmetry
  double Abs2 = -1.0*(m_cp_asymmetry_C-1.0)/(m_cp_asymmetry_C+1.0);
  double Im = m_cp_asymmetry_S/(m_cp_asymmetry_C+1.0);
  double Re = sqrt(Abs2-sqr(Im));
  m_cp_asymmetry_lambda = Complex(Re, Im);

  if (s["BR"].GetItemsCount()==0) ProcessPartonic(s);

  DEBUG_INFO("processing ME block");
  for (auto tmp: s["ME"].GetItems()) {
    auto sme=tmp[tmp.GetKeys()[0]];
    HD_ME_Base* me = SelectME(tmp.GetKeys()[0], sme);
    AddDiagram(me);
  }
  if(m_diagrams.size()==0) {
    DEBUG_INFO("No ME specified, adding Generic.");
    int n=NOut()+1;
    vector<int> decayindices(n);
    for(int i=0;i<n;i++) decayindices[i]=i;
    AddDiagram(new Generic(m_physicalflavours,decayindices,"Generic"));
  }

  DEBUG_INFO("processing PhaseSpace block");
  for (auto tmp: s["PhaseSpace"].GetItems()) {
    string name = tmp.GetKeys()[0];
    auto sps=tmp[tmp.GetKeys()[0]];
    
    if(!AddPSChannel(name, sps))
      THROW(fatal_error, name+" is not a valid phase space channel.");
  }
  if(Channels()->NChannels() == 0) AddPSChannel(string("Isotropic"), s["PhaseSpace"]["Isotropic"]);

  DEBUG_INFO("processing IntResults block");
  auto results = s["IntResults"].SetDefault({-1.0}).GetItems();
  if (results.size()!=3) {
    msg_Info()<<"Calculating width (PR1) for "<<Name()<<endl;
    CalculateWidth();
    msg_Info()<<"   yields: "<<m_iwidth<<endl;
    std::cerr<<"HADRON_DECAYS: { Channels: { "<<m_flavours[0].Kfcode()<<": { "<<FSIDCode()<<": { IntResults: ["
             <<m_iwidth<<", "<<m_ideltawidth<<", "<<m_max<<"]}}}}"<<std::endl;
  }
  else {
    m_iwidth = results[0].Get<double>();
    m_ideltawidth = results[1].Get<double>();
    m_max=results[2].Get<double>();
    
    if (IsNan(m_iwidth)) {
      PRINT_INFO("Found nan in "<<Name()<<". Ignoring and continuing.");
      return;
    }
    
    if (m_always_integrate) {
      double oldwidth(m_iwidth), oldmax(m_max);
      msg_Info()<<"Calculating width (PR2) for "<<Name()<<endl;
      CalculateWidth();
      msg_Info()<<"   yields: "<<m_iwidth<<endl;
      // check whether result is different from before and write out if it is
      if(oldwidth!=m_iwidth || oldmax!=m_max)
        std::cerr<<"HADRON_DECAYS: { Channels: { "<<m_flavours[0].Kfcode()<<": { "<<FSIDCode()<<": { IntResults: ["
                 <<m_iwidth<<", "<<m_ideltawidth<<", "<<m_max<<"]}}}}"<<std::endl;
    }
  }
}

void Hadron_Decay_Channel::ProcessPartonic(Scoped_Settings& s)
{
  DEBUG_FUNC(IDCode());
  std::string me("");
  if (m_flavours.size()==4) {
    if (m_flavours[1].IsPhoton() && m_flavours[2].IsGluon() && m_flavours[3].IsGluon()) {
      me = "QQ_PGG[0,1,2,3]";
    }
    else if ((m_flavours[1].IsFermion() || m_flavours[1].IsDiQuark())&&
             m_flavours[2].IsGluon() && m_flavours[3].IsQuark()) {
      me="QQ_QVQ_Spectator[0,1,2,3]";
      AddPSChannel(string("IsotropicSpectator_1"), s["PhaseSpace"]["IsotropicSpectator_1"]);
    }
  }
  else if (m_flavours.size()==5) {
    if ((m_flavours[1].IsQuark() || m_flavours[1].IsDiQuark())&&
        m_flavours[2].IsFermion() && m_flavours[3].IsFermion() && m_flavours[4].IsFermion()) {
      me="QQ_QQQQ_Spectator[0,1,2,3,4]";
      AddPSChannel(string("IsotropicSpectator_1"), s["PhaseSpace"]["IsotropicSpectator_1"]);
    }
  }
  if (me=="") THROW(fatal_error, "ME not found for "+IDCode());

  AddDiagram(SelectME(me, s["ME"][me]));
}

bool Hadron_Decay_Channel::SetColorFlow(ATOOLS::Blob* blob)
{
  int n_q(0), n_g(0);
  for(int i=0;i<blob->NOutP();i++) {
    if(blob->OutParticle(i)->Flav().IsQuark())      n_q++;
    else if(blob->OutParticle(i)->Flav().IsGluon()) n_g++;
  }
  if(n_q>0 || n_g>0) {
    blob->SetStatus(blob_status::needs_showers);
    Particle_Vector outparts=blob->GetOutParticles();
    if(m_diagrams.size()>0) {
      // try if the matrix element knows how to set the color flow
      HD_ME_Base* firstme=(HD_ME_Base*) m_diagrams[0];
      bool anti=blob->InParticle(0)->Flav().IsAnti();
      if(firstme->SetColorFlow(outparts,n_q,n_g,anti)) return true;
    }
    // otherwise try some common situations
    int n=outparts.size();
    if(n_q==2 && n_g==0 && n==2) {
      if(outparts[0]->Flav().IsAnti()) {
        outparts[0]->SetFlow(2,-1);
        outparts[1]->SetFlow(1,outparts[0]->GetFlow(2));
      }
      else {
        outparts[0]->SetFlow(1,-1);
        outparts[1]->SetFlow(2,outparts[0]->GetFlow(1));
      }
      return true;
    }
    else if(n_q==0 && n_g==2) {
      int inflow(-1), outflow(-1);
      Particle_Vector::iterator pit;
      for(pit=outparts.begin(); pit!=outparts.end(); pit++) {
        if((*pit)->Flav().IsGluon()) {
          (*pit)->SetFlow(2,inflow);
          (*pit)->SetFlow(1,outflow);
          inflow=(*pit)->GetFlow(1);
          outflow=(*pit)->GetFlow(2);
        }
      }
      return true;
    }
    else if(n_q==0 && n_g==n) {
      outparts[0]->SetFlow(2,-1);
      outparts[0]->SetFlow(1,-1);
      for(int i=1;i<n-1;++i) {
        unsigned int c=Flow::Counter();
        outparts[i]->SetFlow(2,c-1);
        outparts[i]->SetFlow(1,c);
      }
      outparts[n-1]->SetFlow(2,outparts[n-2]->GetFlow(1));
      outparts[n-1]->SetFlow(1,outparts[0]->GetFlow(2));
      return true;
    }
    else {
      msg_Error()<<METHOD<<" wasn't able to set the color flow for"<<endl
                 <<*blob<<endl;
      return false;
    }
  }
  else return true;
}

HD_ME_Base * Hadron_Decay_Channel::SelectME(string me_string, Scoped_Settings s)
{
  // TODO convert to YAML eventually
  Data_Reader reader(",",";","#","]");
  reader.AddWordSeparator("[");
  vector<string> resultstrings;
  reader.SetString(me_string);
  reader.VectorFromString(resultstrings);
  if(resultstrings.size()==1 && (resultstrings[0]=="Generic" || resultstrings[0]=="Current_ME")) {
    for(int i=0;i<NOut()+1;i++) 
      resultstrings.push_back( ToString<size_t>(i) );
  }
  if(int(resultstrings.size())!=NOut()+2)
    THROW(fatal_error, "Wrong number of indices in "+FSIDCode());

  int n=NOut()+1;
  vector<int> indices(n);
  for(int i=0; i<n; i++) indices[i] = ToType<int>(resultstrings[i+1]);
  ME_Parameters fi(m_physicalflavours, indices);

  HD_ME_Base* me = HD_ME_Getter_Function::GetObject(resultstrings[0],fi);
  if(me==NULL)
    THROW(fatal_error, "ME not recognized in "+FSIDCode());

  me->SetModelParameters(s);
  return me;
}

void Hadron_Decay_Channel::LatexOutput(std::ostream& f, double totalwidth)
{
  f<<"$"<<GetDecaying().TexName()<<"$ $\\to$ ";
  for (size_t i=1; i<m_flavours.size(); ++i)
    f<<"$"<<m_flavours[i].TexName()<<"$ ";
  f<<" & ";
  char helpstr[100];
  snprintf( helpstr, 100, "%.4f", Width()/totalwidth*100. );
  f<<helpstr;
  if( DeltaWidth() > 0. ) {
    snprintf( helpstr, 100, "%.4f", DeltaWidth()/totalwidth*100. );
    f<<" $\\pm$ "<<helpstr;
  }
  f<<" \\% ";
  if(Origin()!="") {
    f<<"[\\verb;"<<Origin()<<";]";
  }
  f<<"\\\\"<<endl;
  if((m_diagrams.size()>0 &&
      ((HD_ME_Base*) m_diagrams[0])->Name()!="Generic")) {
    snprintf( helpstr, 100, "%.4f", IWidth()/totalwidth*100. );
    f<<" & "<<helpstr;
    if( IDeltaWidth() > 0. ) {
      snprintf( helpstr, 100, "%.4f", IDeltaWidth()/totalwidth*100. );
      f<<" $\\pm$ "<<helpstr;
    }
    f<<" \\% ";
  }
  for(size_t i=0;i<m_diagrams.size();i++) {
    HD_ME_Base* me=(HD_ME_Base*) m_diagrams[i];
    if(me->Name()=="Current_ME") {
      Current_ME* cme=(Current_ME*) me;
      f<<"\\verb;"<<cme->GetCurrent1()->Name()
       <<";$\\otimes$\\verb;"<<cme->GetCurrent2()->Name()<<"; & \\\\"<<endl;
    }
    else if (me->Name()=="Generic") {
      // do nothing
    }
    else {
      f<<"\\verb;"<<me->Name()<<"; & \\\\"<<endl;
    }
  }
}

bool Hadron_Decay_Channel::AddPSChannel(string name,
                                        Scoped_Settings s)
{
  PHASIC::Single_Channel * sc=
    HD_Channel_Selector::GetChannel(1, NOut(),&m_flavours.front(),name,s,p_ms);
  if (sc!=NULL) {
    sc->SetAlpha(s["Weight"].SetDefault(1.0).Get<double>());
    Channels()->Add(sc);
    return true;
  }
  else return false;
}
