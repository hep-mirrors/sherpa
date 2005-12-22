#include "Amegic_Apacic_Interface.H"

#include "XS_Selector.H" 
#include "Data_Read.H"
#include "Random.H"
#include "Message.H"
#include "Running_AlphaS.H"
#include "Run_Parameter.H"
#include "Amegic.H"
#include "Cluster_Partons_CKKW.H"

using namespace SHERPA;
using namespace AMEGIC;
using namespace APACIC;
using namespace EXTRAXS;
using namespace MODEL;
using namespace ATOOLS;

Amegic_Apacic_Interface::Amegic_Apacic_Interface(Matrix_Element_Handler * me,
						 Shower_Handler * shower) :
  Perturbative_Interface(me,shower), m_lastshowerveto(0), 
  p_cluster(NULL), p_filler(NULL), p_xs(NULL),
  p_two2two(new XS_Group(2,2,"Core processes")),
  p_one2N(new XS_Group(1,3,"Simple decays")),
  p_blob_psme_IS(NULL), p_blob_psme_FS(NULL)
{
  p_fl      = new Flavour[4];
  p_moms    = new Vec4D[4];
  p_cluster = new Cluster_Partons_CKKW
    (p_mehandler,shower->JetFinder(),m_maxjetnumber,
     p_shower->GetISRHandler()->On(),
     p_shower->ISROn(),p_shower->FSROn());
  p_filler  = new Tree_Filler
    (p_cluster,m_maxjetnumber,p_shower->ISROn(),p_shower->FSROn());
}  

Amegic_Apacic_Interface::~Amegic_Apacic_Interface() 
{
  if (p_two2two) { delete p_two2two; p_two2two = NULL; }
  if (p_one2N)   { delete p_one2N;   p_one2N   = NULL; }
  if (p_cluster) { delete p_cluster; p_cluster = NULL; }
  if (p_filler)  { delete p_filler;  p_filler = NULL;  }
}


bool Amegic_Apacic_Interface::ClusterConfiguration(Blob * blob)
{
  if (blob->NInP()==1) {
    m_isdecay = 1;
    if (!(p_cluster->ClusterConfiguration(blob))) return 0; 
    p_fl[0]   = blob->InParticle(0)->Flav();
    p_moms[0] = blob->InParticle(0)->Momentum();
    for (int i=0;i<blob->NOutP();i++) {
      p_fl[i+1]   = blob->OutParticle(i)->Flav();
      p_moms[i+1] = blob->OutParticle(i)->Momentum();
    }
  }
  else {
    m_isdecay = 0;
    if (!(p_cluster->ClusterConfiguration
	  (blob,p_mehandler->GetISR_Handler()->X1(),
	   p_mehandler->GetISR_Handler()->X2()))) {
      return 0;
    }
    for (int i=0;i<4;i++) {
      p_fl[i]   = p_cluster->Flav(i); 
      p_moms[i] = p_cluster->Momentum(i);
    }
  }
  // prepare Blob , will be inserted later
  if (p_blob_psme_IS) { delete p_blob_psme_IS; p_blob_psme_IS = 0; }
  if (p_blob_psme_FS) { delete p_blob_psme_FS; p_blob_psme_FS = 0; }

  if (!m_isdecay && p_shower->ISROn()) {
    p_blob_psme_IS = new Blob();
    p_blob_psme_IS->SetType(btp::ME_PS_Interface_IS);
    p_blob_psme_IS->SetTypeSpec("Sherpa");
    p_blob_psme_IS->SetStatus(1);
    p_blob_psme_IS->SetId();
    for (int i=0;i<blob->NInP();++i) {
      p_blob_psme_IS->AddToOutParticles(blob->InParticle(i));
    }
  }
  Blob * dec;
  if (p_shower->FSROn()) {
    p_blob_psme_FS = new Blob();
    p_blob_psme_FS->SetType(btp::ME_PS_Interface_FS);
    p_blob_psme_FS->SetTypeSpec("Sherpa");
    p_blob_psme_FS->SetStatus(1);
    p_blob_psme_FS->SetId();
    for (int i=0;i<blob->NOutP();++i) {
      dec = NULL;
      if (blob->OutParticle(i)->DecayBlob()) {
	if (blob->OutParticle(i)->DecayBlob()->Type()==btp::Hard_Decay) 
	  dec = blob->OutParticle(i)->DecayBlob();
      }
      p_blob_psme_FS->AddToInParticles(blob->OutParticle(i));
      if (dec) blob->OutParticle(i)->SetDecayBlob(dec);
    }
  }
  return 1;
}

int Amegic_Apacic_Interface::DefineInitialConditions(ATOOLS::Blob * blob)
{
  if (!p_shower->ISROn() && !p_shower->FSROn())  {
    m_weight = 1.;
    return 1;
  }
  m_nin  = blob->NInP();
  m_nout = blob->NOutP();
  ClusterConfiguration(blob);
  if (!m_isdecay) {
    p_xs = p_cluster->GetXS(p_two2two,p_fl);
    p_cluster->SetColours(p_xs,p_moms,p_fl);
    p_cluster->CalculateWeight();
    p_blob_psme_FS->
      AddData("OrderStrong",new ATOOLS::Blob_Data<double>
	      (p_cluster->OrderStrong()));
    p_blob_psme_FS->
      AddData("OrderEWeak",new ATOOLS::Blob_Data<double>
	      (p_cluster->OrderEWeak()));
    p_blob_psme_FS->
      AddData("HardScale",new ATOOLS::Blob_Data<double>
	      (p_cluster->HardScale()));
    p_blob_psme_FS->AddData("QCDScale",new ATOOLS::Blob_Data<double>
			    (p_cluster->QCDScale()));
    p_blob_psme_FS->AddData
      ("Core_Process",new ATOOLS::Blob_Data<XS_Base*>(p_xs));
    m_weight = p_cluster->Weight();
    if (p_mehandler->Weight()==1. && p_mehandler->UseSudakovWeight()) {
      if (m_weight>ran.Get()) {
	p_shower->CleanUp();
	p_filler->FillTrees(blob,p_shower->GetIniTrees(),
			    p_shower->GetFinTree());
	blob->AddData("Sud_Weight",new Blob_Data<double>(m_weight));
	m_weight=1.;
	return 1;
      }
      if (m_lastshowerveto==-1) {
	m_weight=1.;
	return 3;
      }
      m_weight=1.;
      return 0;
    }
    else {
      if (!(p_mehandler->UseSudakovWeight())) m_weight=1.;
      else {
	// update me weight
	Blob_Data_Base * message = (*blob)["ME_Weight"];
	if (message) {
	  double weight = message->Get<double>();
	  weight*=m_weight;
	  blob->AddData("ME_Weight",new Blob_Data<double>(weight));
	  blob->AddData("Sud_Weight",new Blob_Data<double>(m_weight));
	}
	else {
	  msg.Out()<<"Amegic_Apacic_Interface::DefineInitialConditions(..): "
		   <<"Missing signal process information."<<std::endl;
	}
      }
      p_filler->FillTrees(blob,p_shower->GetIniTrees(),p_shower->GetFinTree());
      return 1;
    }
  }
  else {
    p_xs = 0;
    if (XS_Selector::FindInGroup(p_one2N,p_xs,m_nin,m_nout,p_fl)==
	std::string::npos) {
      EXTRAXS::XS_Selector selector(p_xs);
      p_xs = selector.GetXS(m_nin,m_nout,p_fl);
      if (p_xs) p_one2N->Add(p_xs);
    }
    if (p_xs) {
      if (!(p_xs->SetColours(p_moms))) {
	msg.Error()<<"Amegic_Apacic_Interface::DefineInitialConditions(..): "
		   <<"Extra_XS unable to define colour flow. Return 0."
		   <<std::endl;
	return 0;
      }
    }
    else {
      int col1 = blob->InParticle(0)->GetFlow(1);
      int col2 = blob->InParticle(0)->GetFlow(2);
      if (p_cluster->SetDecayColours(p_moms,p_fl,col1,col2)!=0) {
	msg.Error()<<"Amegic_Apacic_Interface::DefineInitialConditions(..): "
		   <<"No Extra_XS found to define colour flow. Return 0."
		   <<std::endl;
	return 0;
      }
    }
    p_filler->FillDecayTree(p_shower->GetFinTree());
  }
  return 1;
}



bool Amegic_Apacic_Interface::FillBlobs(ATOOLS::Blob_List * bl)
{
  if (p_blob_psme_IS) {
    p_blob_psme_IS->SetId();
    bl->push_back(p_blob_psme_IS);  
    p_blob_psme_IS=0;
  }
  if (p_blob_psme_FS) {
    p_blob_psme_FS->SetId();
    for (int i=0;i<p_blob_psme_FS->NInP();++i) {
      if (p_blob_psme_FS->InParticle(i)->DecayBlob()!=p_blob_psme_FS) 
	p_blob_psme_FS->InParticle(i)->SetDecayBlob(p_blob_psme_FS);
    }
    bl->push_back(p_blob_psme_FS);  
    p_blob_psme_FS=0;
  }
  return 1;
}

int Amegic_Apacic_Interface::PerformShowers()
{
  int jetveto(-1), losejv(1);
  if (p_mehandler->UseSudakovWeight()) {
    double qmin2i,qmin2f; 
    double scale(p_mehandler->FactorisationScale());
    double ycut(p_cluster->YCut());
    double yp(((Process_Base*)p_mehandler->GetAmegic()->GetProcess())->Ycut());
    if (yp!=-1.0) ycut=yp; 
    p_cluster->JetvetoPt2(qmin2i,qmin2f);
    p_shower->SetFactorisationScale(scale);
    jetveto=1;
    if (m_maxjetnumber==m_nout && p_mehandler->OrderStrong()==0) jetveto = 0;
    size_t qcdjets(0);
    Process_Base *proc((Process_Base*)p_mehandler->GetAmegic()->GetProcess());
    for (size_t i(0);i<proc->NOut();++i)
      if (proc->Flavours()[proc->NIn()+i].Strong()) ++qcdjets;
    if (p_mehandler->MinQCDJets()==qcdjets) losejv=0;
  }
  m_lastshowerveto = 
    p_shower->PerformShowers(jetveto,losejv,
			     p_mehandler->GetISR_Handler()->X1(),
			     p_mehandler->GetISR_Handler()->X2(),
			     p_cluster->YCut());
  return m_lastshowerveto;
}


int Amegic_Apacic_Interface::PerformDecayShowers()
{
  return p_shower->PerformDecayShowers(0);
}
