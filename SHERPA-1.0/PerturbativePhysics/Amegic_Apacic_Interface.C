#include "Amegic_Apacic_Interface.H"

#include "XS_Selector.H" 
#include "Data_Read.H"
#include "Random.H"
#include "Message.H"
#include "Running_AlphaS.H"
#include "Run_Parameter.H"

using namespace SHERPA;
using namespace EXTRAXS;
using namespace AMEGIC;
using namespace APACIC;
using namespace ATOOLS;
using namespace MODEL;

using namespace EXTRAXS;

Amegic_Apacic_Interface::Amegic_Apacic_Interface(Matrix_Element_Handler * me,
						 Shower_Handler * shower) :
  Perturbative_Interface(me,shower), m_lastshowerveto(0), p_blob_psme_IS(0), p_blob_psme_FS(0)
{
  p_jf      = 0;
  p_cluster = 0;
  p_xs      = 0;

  p_two2two = new XS_Group(2,2,"Core processes");
  p_one2N   = new XS_Group(1,3,"Simple decays");
  p_fl      = new Flavour[4];
  p_moms    = new Vec4D[4];

  m_ycut    = rpa.gen.Ycut();

  if (rpa.gen.Beam1().IsLepton() && rpa.gen.Beam2().IsLepton()) m_type = 1;
  else if ((!rpa.gen.Beam1().IsLepton() && !rpa.gen.Beam2().IsLepton())) m_type = 4;
  else m_type = 4;

  p_jf       = new ATOOLS::Jet_Finder(m_ycut,m_type);
  p_cluster  = new Cluster_Partons(p_mehandler,p_jf,m_maxjetnumber,
				   p_shower->GetISRHandler()->On(),p_shower->ISROn(),p_shower->FSROn());


  double dr   = rpa.gen.DeltaR();
  double ycut = rpa.gen.Ycut();
  m_qmin_i = sqrt(ycut)*rpa.gen.Ecms();
  m_qmin_f = sqrt(ycut)*dr*rpa.gen.Ecms();
  m_jetscale = rpa.gen.RenormalizationScaleFactor() * sqr(Min(m_qmin_i,m_qmin_f));
}  

Amegic_Apacic_Interface::~Amegic_Apacic_Interface() 
{
  if (p_two2two) { delete p_two2two; p_two2two = NULL; }
  if (p_one2N)   { delete p_one2N;   p_one2N   = NULL; }
  if (p_jf)      { delete p_jf;      p_jf      = NULL; }
  if (p_cluster) { delete p_cluster; p_cluster = NULL; }
  // note :
  //  p_shower and p_mehandler are deleted in Jet_Evolution
  //  p_fl an p_moms are deleted in Perturbative_Interface
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
    if (!(p_cluster->ClusterConfiguration(blob,p_mehandler->GetISR_Handler()->X1(),
					  p_mehandler->GetISR_Handler()->X2()))) {
      return 0; // Failure!
    }
    for (int i=0;i<4;i++) {
      p_fl[i]   = p_cluster->Flav(i); 
      p_moms[i] = p_cluster->Momentum(i);
    }
  }

  // prepare Blob , will be inserted later
  /*
    std::cout<<"Delete two blobs : "<<std::endl;
    if (p_blob_psme_IS) std::cout<<p_blob_psme_IS<<std::endl;
    if (p_blob_psme_FS) std::cout<<p_blob_psme_FS<<std::endl;
  */
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

  return 1;  // OK!
}

int Amegic_Apacic_Interface::DefineInitialConditions(ATOOLS::Blob * blob)
{
  if (!p_shower->ISROn() && !p_shower->FSROn())  {
    m_weight = 1.;
    return 1;
  }

  //  PROFILE_LOCAL("Amegic_Apacic_Interface::DefineInitialConditions");
  m_nin  = blob->NInP();
  m_nout = blob->NOutP();

  ClusterConfiguration(blob);

  if (!m_isdecay) {
    p_xs = p_cluster->GetXS(p_two2two,p_fl);
    
    p_cluster->SetColours(p_xs,p_moms,p_fl);

    
    if (m_type==1) { // e+ e-
      double sprime = (p_mehandler->Momenta()[0]+p_mehandler->Momenta()[1]).Abs2();
      m_jetscale    = rpa.gen.RenormalizationScaleFactor() * m_ycut * sprime;
    }
    
    double asscale;
    if (p_xs) {
      //      p_xs->Selected()->CalculateScale(p_moms);
      //      std::cout<<" scale="<<p_xs->Selected()->Scale(PHASIC::stp::as);
      asscale=m_scale = p_xs->Selected()->Scale(PHASIC::stp::as);
    }
    else {
      m_scale   = p_cluster->Scale();
      asscale = p_cluster->AsScale();
    }
    
    p_cluster->CalculateWeight(m_scale,asscale,m_jetscale,m_qmin_i,m_qmin_f);

    p_blob_psme_FS->AddData("d#hard_scale",new ATOOLS::Blob_Data<double>(sqrt(m_scale)));
    p_blob_psme_FS->AddData("d#as_scale",new ATOOLS::Blob_Data<double>(sqrt(asscale)));

    // try alternative scale for UE
    m_scale=asscale;
    p_cluster->FixJetvetoPt2(m_scale);
    p_blob_psme_FS->AddData("d#soft_scale",new ATOOLS::Blob_Data<double>(sqrt(asscale)));

    p_blob_psme_FS->AddData("Core_Process",new ATOOLS::Blob_Data<XS_Base*>(p_xs));
    p_blob_psme_FS->AddData("OrderStrong",new ATOOLS::Blob_Data<double>(double(p_cluster->OrderStrong())));
    p_blob_psme_FS->AddData("OrderEWeak",new ATOOLS::Blob_Data<double>(double(p_cluster->OrderEWeak())));



    m_weight = p_cluster->Weight();
    if (p_mehandler->Weight()==1. && p_mehandler->UseSudakovWeight()) {
      if (m_weight>ran.Get()) {
	p_shower->CleanUp();
	p_cluster->FillTrees(p_shower->GetIniTrees(),p_shower->GetFinTree());
	m_weight=1.;
	return 1;
      }
      if (m_lastshowerveto==3) {
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
	}
	else {
	  msg.Out()<<"WARNING in Amegic_Apacic_Interface::DefineInitialConditions: "<<std::endl
		   <<"   ME_Weight not found in Amegic_Apacic_Interface::DefineInitialCondition() "<<std::endl;
	}

      }
      p_shower->CleanUp();
      p_cluster->FillTrees(p_shower->GetIniTrees(),p_shower->GetFinTree());
      return 1;
    }
  }
  else {
    p_xs = 0;
    if (XS_Selector::FindInGroup(p_one2N,p_xs,m_nin,m_nout,p_fl)==std::string::npos) {
      EXTRAXS::XS_Selector selector(p_xs);
      p_xs = selector.GetXS(m_nin,m_nout,p_fl);
      if (p_xs) p_one2N->Add(p_xs);
    }
    if (p_xs) {
      if (!(p_xs->SetColours(p_moms))) {
	msg.Error()<<"ERROR in Amegic_Apacic_Interface::DefineInitialConditions."<<std::endl
		   <<"   Extra_XS was unable to define colour flow. Return 0."<<std::endl;
	return 0;
      }
    }
    else {
      int col1 = blob->InParticle(0)->GetFlow(1);
      int col2 = blob->InParticle(0)->GetFlow(2);
      if (p_cluster->SetDecayColours(p_moms,p_fl,col1,col2)!=0) {
	msg.Error()<<"ERROR in Amegic_Apacic_Interface::DefineInitialConditions."<<std::endl
		   <<"   No Extra_XS found to define colour flow. Return 0."<<std::endl;
	return 0;
      }
    }
    p_cluster->FillDecayTree(p_shower->GetFinTree(),p_xs);
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
  // PROFILE_LOCAL("Amegic_Apacic_Interface::PerformShowers");
  int jetveto=-1;
  int losejv=1;
  if (p_mehandler->UseSudakovWeight()) {
    double qmin2i,qmin2f; 
    double scale = p_mehandler->FactorisationScale();
    double ycut=p_cluster->Ycut();
    if (ycut==rpa.gen.Ycut()) {
      p_cluster->JetvetoPt2(qmin2i,qmin2f);
    }
    else {
      qmin2i=qmin2f=scale;
    }
    p_shower->SetJetvetoPt2(qmin2i,qmin2f);
    p_shower->SetFactorisationScale(scale);
    jetveto=1;
    if (m_maxjetnumber==m_nout && p_mehandler->OrderStrong()==0) {
      //      std::cout<<" no jetveto mu_f="<<scale<<std::endl;
      jetveto=0;
    }
    if (m_nout==2) losejv=0;
  }
  return m_lastshowerveto = p_shower->PerformShowers(jetveto,losejv,
						     p_mehandler->GetISR_Handler()->X1(),
						     p_mehandler->GetISR_Handler()->X2(),
						     p_cluster->Ycut());
}


int Amegic_Apacic_Interface::PerformDecayShowers()
{
  return p_shower->PerformDecayShowers(0);
}
