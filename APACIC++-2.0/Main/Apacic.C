#include "Apacic.H"
#include "Initial_State_Shower.H"
#include "Final_State_Shower.H"
#include "Tree.H"

#include "ISR_Handler.H"
#include "ISR_Info.H"
#include "Blob.H"

#include "Run_Parameter.H"
#include "Random.H"

#ifdef NO_ANALYSIS__all
#define NO_ANALYSIS__Apacic
#endif

using namespace APACIC;
using namespace PDF;
using namespace ATOOLS;


Apacic::Apacic(ISR_Handler * isr,MODEL::Model_Base * model,int maxjetnumber,
	       bool isron,bool fsron,Data_Read * dataread):
  m_isron(isron), m_fsron(fsron), m_showers(isron||fsron),
  p_inishower(NULL), p_finshower(NULL), p_initrees(NULL), p_fintree(NULL),
  m_info_cms(8), m_info_lab(8)
{
  if (m_fsron) {
    p_fintree   = new Tree();
    p_finshower = new Final_State_Shower(model,dataread);
  }
  if (m_isron) {
    p_initrees  = new Tree*[2];
    for (int i=0;i<2;i++) p_initrees[i] = new Tree();
    p_inishower = new Initial_State_Shower(isr,p_finshower,model,dataread);
  }
}
  
Apacic::~Apacic() 
{
  if (p_fintree)     { delete p_fintree; p_fintree = 0; }
  if (p_initrees)    {
    for (int i=0;i<2;i++) { delete p_initrees[i]; p_initrees[i] = 0; }
    delete p_initrees;p_initrees = 0;
  }  
  if (p_inishower)   { delete p_inishower; p_inishower = 0; }
  if (p_finshower)   { delete p_finshower; p_finshower = 0; }
}

void Apacic::PrepareTrees() {
  if (m_fsron) p_fintree->Reset(); 
  if (m_isron) for (int i=0;i<2;i++) p_initrees[i]->Reset();
}

void Apacic::SetJetvetoPt2(const double q2i, const double q2f)
{ 
//   std::cout<<" Apacic::SetJetvetoPt2("<<q2i<<","<<q2f<<")\n";
  
  if (m_fsron) 
    p_finshower->SetJetvetoPt2(q2f); 
  if (m_isron)
    p_inishower->SetJetvetoPt2(q2i,q2f); 
}

void Apacic::SetFactorisationScale(const double scale)
{
//   std::cout<<" Apacic::SetFactorisationScale("<<scale<<")\n";
  if (m_isron) {
    p_inishower->SetFactorisationScale(scale);
  }
}

int Apacic::PerformShowers(bool ini,bool fin,int jetveto,double x1,double x2, double ycut) {
  if (!m_showers) return 1;
  if (msg.LevelIsDebugging()) {
    if (m_fsron) 
      p_finshower->OutputTree(p_fintree);
    if (m_isron) {
      p_inishower->OutputTree(p_initrees[0]);
      p_inishower->OutputTree(p_initrees[1]);
    }
  }
  if (m_fsron) {
    Vec4D sum=p_fintree->GetRoot()->part->Momentum();
    Poincare cms(sum);
    p_fintree->BoRo(cms);

    int fsrstatus = p_finshower->PerformShower(p_fintree,jetveto);
    int number=0;

    if (fsrstatus!=1) {
      if (fsrstatus!=3) msg.Out()<<" "<<fsrstatus<<" FS shower failed !\n";
      else msg_Debugging()<<" FS shower asks for new event\n";
      return fsrstatus;
    }
 
    // check ME if still njet ME!
    // if isr is on, this check will be performed after the initial state shower
    if (!m_isron) {
      //      double ycut   = rpa.gen.Ycut();
      double s  = sqr(rpa.gen.Ecms());
      double dr2 = sqr(rpa.gen.DeltaR());
      SetJetvetoPt2(s*ycut,dr2*s*ycut);
      if (!p_finshower->ExtraJetCheck()) {
	return 3;
      }
    }
    p_finshower->SetAllColours(p_fintree->GetRoot());

    if (!m_isron) {
      //      Vec4D vl =Vec4D(x1+x2,0.,0.,x2-x1);
      Vec4D vl =Vec4D(sum[0],-1.*Vec3D(sum));
      Poincare lab(vl);
      p_fintree->BoRo(lab);
      Vec4D sum_fs=p_finshower->GetMomentum(p_fintree->GetRoot(),number);
      //      std::cout<<" sum_fs:"<<sum_fs<<"\n";
    }
  }

  if (m_isron) {
    p_inishower->InitShowerPT(p_initrees[0]->GetRoot()->maxpt2);
    if (!(p_inishower->PerformShower(p_initrees,jetveto))) return 0;


    // boost in Beam-System (using x1, x2 from initiator) done in is-shower

    // determine Boost+Rotate in cms(using moms from root1 und root2)
    Vec4D mom1=p_initrees[0]->GetRoot()->part->Momentum();
    Vec4D mom2=p_initrees[1]->GetRoot()->part->Momentum();

    //    std::cout<<" sum_is:"<<mom1+mom2<<std::endl;

    Vec4D vl =Vec4D(mom1[0]+mom2[0], -1.*Vec3D(mom1+mom2));
    Poincare lab(vl);
    lab.BoostBack(mom1);
    lab.BoostBack(mom2);
    Poincare rot(Vec4D::ZVEC,mom1);
    rot.RotateBack(mom1);
    rot.RotateBack(mom2);

#ifndef NO_ANALYSIS__Apacic
    mom1=p_initrees[0]->GetRoot()->part->Momentum();
    mom2=p_initrees[1]->GetRoot()->part->Momentum();
    
    m_info_lab[iic::E_1]=mom1[0];
    m_info_lab[iic::t_1]=mom1.Abs2();
    m_info_lab[iic::Em_1]=mom1[0]/mom1.Mass();
    m_info_lab[iic::E_2]=mom2[0];
    m_info_lab[iic::t_2]=mom2.Abs2();
    m_info_lab[iic::Em_2]=mom2[0]/mom2.Mass();
    m_info_lab[iic::z_1]=p_initrees[0]->GetRoot()->z;
    m_info_lab[iic::z_2]=p_initrees[1]->GetRoot()->z;

    Poincare test(mom1+mom2);
    test.Boost(mom1);
    test.Boost(mom2);

    m_info_cms[iic::E_1]=mom1[0];
    m_info_cms[iic::t_1]=mom1.Abs2();
    m_info_cms[iic::Em_1]=mom1[0]/mom1.Mass();
    m_info_cms[iic::E_2]=mom2[0];
    m_info_cms[iic::t_2]=mom2.Abs2();
    m_info_cms[iic::Em_2]=mom2[0]/mom2.Mass();
    m_info_lab[iic::mu_1]=p_initrees[0]->GetRoot()->t;
    m_info_lab[iic::mu_2]=p_initrees[1]->GetRoot()->t;
#endif

    // rotate and boost fs in beam-system "if (fin)"
    if (fin) {
      p_fintree->BoRo(rot);
      p_fintree->BoRo(lab);

      // check ME if still njet ME!
      //      double ycut   = rpa.gen.Ycut();
      double s  = sqr(rpa.gen.Ecms());
      double dr2 = sqr(rpa.gen.DeltaR());
      SetJetvetoPt2(s*ycut,dr2*s*ycut);
//       std::cout<<"extra_jet "<<s*ycut<<"/n";
      int number=0;
      Vec4D sum_fs=p_finshower->GetMomentum(p_fintree->GetRoot(),number);
      //      std::cout<<" sum_fs:"<<sum_fs<<"\n";

      //      std::cout<<" sum_left : "<<p_finshower->GetMomentum(p_initrees[0]->GetInitiator(),number)<<"\n";
      //      std::cout<<" sum_right : "<<p_finshower->GetMomentum(p_initrees[1]->GetInitiator(),number)<<"\n";
      if (!p_finshower->ExtraJetCheck()) {
	//	std::cout<<" extrajetcheck failed"<<std::endl;
	return 3;
      }
    }
  }

  if (msg.LevelIsDebugging()) {
    if (m_fsron) 
      p_finshower->OutputTree(p_fintree);
    if (m_isron) {
      p_inishower->OutputTree(p_initrees[0]);
      p_inishower->OutputTree(p_initrees[1]);
    }
  }
  return 1;
}

bool Apacic::ExtractPartons(bool ini,bool fin,Blob_List * bl,Particle_List * pl) {
  if (fin) {
    p_fintree->CheckStructure(true);
    p_finshower->ExtractPartons(p_fintree->GetRoot(),0,bl,pl);
  }
  if (ini) {
    for (int i=0;i<2;i++) {
      p_initrees[i]->CheckStructure(true);
      p_inishower->ExtractPartons(p_initrees[i]->GetInitiator(),i,0,bl,pl);
    }
  }
//   for (Blob_List::iterator bit=bl->begin();bit!=bl->end();++bit) {
//     if ((*bit)->Type()==btp::Signal_Process) {
//       (*bit)->AddData("Shower_Info_cms",new Blob_Data<std::vector<double> >(m_info_cms));
//       (*bit)->AddData("Shower_Info_lab",new Blob_Data<std::vector<double> >(m_info_lab));
//       break;
//     }
//   }
  return 1;
}

ATOOLS::Blob_Data_Base *const Apacic::Info(const int frame) const
{
  if (frame==0) return new ATOOLS::Blob_Data<std::vector<double> >(m_info_cms);
  return new ATOOLS::Blob_Data<std::vector<double> >(m_info_lab);
}
