#include "Simple_Chain.H"

#include "Phase_Space_Handler.H"
#include "Running_AlphaS.H"
#include "Single_XS.H"
#include "XS_Selector.H"
#include "Channel_Elements.H"
#include "Particle.H"
#include "Random.H"
#ifdef USING__Sherpa
#include "Matrix_Element_Handler.H"
#endif
#include <sys/stat.h>

#ifdef PROFILE__Simple_Chain
#include "prof.hh"
#else
#define PROFILE_HERE
#endif

#ifdef DEBUG__Simple_Chain
const std::string differentialfile=std::string("differential.dat");
const std::string integralfile=std::string("integral.dat");
const std::string normalizedfile=std::string("normalized.dat");
#endif

using namespace AMISIC;

Simple_Chain::Simple_Chain():
  MI_Base("Simple Chain",MI_Base::HardEvent,4,4,1),
  m_differential(std::vector<GridFunctionType*>(0)),
  m_norm((GridResultType)1.0),
  p_total(NULL),
  m_xsextension("_xs.dat"),
  m_maxextension("_max.dat"),
  m_mcextension("_mc"),
  p_processes(NULL),
  p_fsrinterface(NULL), 
  p_environment(NULL),
  p_model(NULL),
  p_beam(NULL),
  p_isr(NULL),
  m_nflavour(3),
  m_scalescheme(2),
  m_kfactorscheme(1),
  m_maxtrials(1),
  m_external(false)
{
  SetInputFile("MI.dat");
  SetInputFile("XS.dat",1);
  SetInputFile("Run.dat",2);
  SetInputFile("Model.dat",3);
  SetOutputFile("SC.log");
  SetStart(0.0,0);
  SetStop(1.6,0);
  SetStart(ATOOLS::rpa.gen.Ecms(),1);
  SetStop(0.0,1);
}

Simple_Chain::Simple_Chain(MODEL::Model_Base *_p_model,
			   BEAM::Beam_Spectra_Handler *_p_beam,PDF::ISR_Handler *_p_isr):
  MI_Base("Simple Chain",MI_Base::HardEvent,4,4,1),
  m_differential(std::vector<GridFunctionType*>(0)),
  m_norm((GridResultType)1.0),
  p_total(NULL),
  m_xsfile(std::string("XS.dat")),
  m_xsextension("_xs.dat"),
  m_maxextension("_max.dat"),
  m_mcextension("_mc"),
  p_processes(NULL),
  p_fsrinterface(NULL), 
  p_environment(NULL),
  p_model(_p_model),
  p_beam(_p_beam),
  p_isr(_p_isr),
  m_nflavour(3),
  m_scalescheme(2),
  m_kfactorscheme(1),
  m_maxtrials(100),
  m_external(true)
{
  SetInputFile("MI.dat");
  SetInputFile("XS.dat",1);
  SetInputFile("Run.dat",2);
  SetInputFile("Model.dat",3);
  SetOutputFile("SC.log");
  SetStart(0.0,0);
  SetStop(1.6,0);
  SetStart(ATOOLS::rpa.gen.Ecms(),1);
  SetStop(0.0,1);
}

Simple_Chain::~Simple_Chain()
{
  CleanUp();
}

void Simple_Chain::CleanUp() 
{
  if (p_fsrinterface!=NULL) delete p_fsrinterface;
#ifndef USING__Sherpa
  if (p_processes!=NULL) delete p_processes;
#endif
  if (!m_external) {
    if (p_environment!=NULL) delete p_environment;
    p_environment=NULL;
    p_model=NULL;
    p_beam=NULL;
    p_isr=NULL;
  }
  if (p_total!=NULL) delete p_total;
  while (m_blobs.size()>0) {
    while (m_blobs.begin()->size()>0) {
      delete *m_blobs.begin()->begin();
      m_blobs.begin()->erase(m_blobs.begin()->begin());
    }
    m_blobs.erase(m_blobs.begin());
  }
  m_filename.clear();
  while (m_differential.size()>0) {
    delete *m_differential.begin();
    m_differential.erase(m_differential.begin());
  }
}

bool Simple_Chain::HaveBlob(ATOOLS::Blob *blob) 
{
  for (std::vector<Blob_List>::iterator blit=m_blobs.begin();blit!=m_blobs.end();++blit) {
    for (Blob_Iterator bit=blit->begin();bit!=blit->end();++bit) {
      if (((*(*bit)->InParticle(0)==*blob->InParticle(0))&&(*(*bit)->InParticle(1)==*blob->InParticle(1))&&
	   (*(*bit)->OutParticle(0)==*blob->OutParticle(0))&&(*(*bit)->OutParticle(1)==*blob->OutParticle(1)))||
	  ((*(*bit)->InParticle(1)==*blob->InParticle(0))&&(*(*bit)->InParticle(0)==*blob->InParticle(1))&&
	   (*(*bit)->OutParticle(0)==*blob->OutParticle(0))&&(*(*bit)->OutParticle(1)==*blob->OutParticle(1)))||
	  ((*(*bit)->InParticle(1)==*blob->InParticle(0))&&(*(*bit)->InParticle(0)==*blob->InParticle(1))&&
	   (*(*bit)->OutParticle(1)==*blob->OutParticle(0))&&(*(*bit)->OutParticle(0)==*blob->OutParticle(1)))||
	  ((*(*bit)->InParticle(0)==*blob->InParticle(0))&&(*(*bit)->InParticle(1)==*blob->InParticle(1))&&
	   (*(*bit)->OutParticle(1)==*blob->OutParticle(0))&&(*(*bit)->OutParticle(0)==*blob->OutParticle(1)))) {
	return true;
      }
    }
  }
  return false;
}

ATOOLS::Blob *Simple_Chain::GetBlob(ATOOLS::Flavour *flavour)
{
  ATOOLS::Blob *newblob = new ATOOLS::Blob();
  for (unsigned int j=0;j<2;++j) newblob->AddToInParticles(new ATOOLS::Particle(j,flavour[j]));
  for (unsigned int k=2;k<4;++k) newblob->AddToOutParticles(new ATOOLS::Particle(k,flavour[k]));
  if (HaveBlob(newblob)) {
    ATOOLS::msg.Tracking()<<"Simple_Chain::GetBlob("<<flavour<<"): "
			  <<"Tried to add a blob twice. Ignore last attempt."<<std::endl; 
    delete newblob;
    return NULL;
  }
  return newblob;
}

void Simple_Chain::FillMode(EXTRAXS::QCD_Processes::Mode mode)
{
  PROFILE_HERE;
  ATOOLS::Flavour temp[4];
  ATOOLS::Blob *newblob;
  unsigned int i, j;
  switch (mode) {
  case EXTRAXS::QCD_Processes::All:
  case EXTRAXS::QCD_Processes::gggg:
    m_blobs.push_back(ATOOLS::Blob_List(0));
    for (i=0;i<4;++i) temp[i]=ATOOLS::Flavour((ATOOLS::kf::code)21);
    if ((newblob=GetBlob(temp))!=NULL) {
      m_blobs[m_blobs.size()-1].push_back(newblob);
    }
    m_filename.push_back("gg_to_gg__grid");
//     m_create.push_back(false);
    if (mode==EXTRAXS::QCD_Processes::gggg) break;
  case EXTRAXS::QCD_Processes::qqbgg:
    m_blobs.push_back(ATOOLS::Blob_List(0));
    temp[3]=temp[2]=ATOOLS::Flavour((ATOOLS::kf::code)21);
    for (i=1;i<=m_nflavour;++i) {
      temp[0]=ATOOLS::Flavour((ATOOLS::kf::code)i);
      temp[1]=ATOOLS::Flavour((ATOOLS::kf::code)i).Bar();
      if ((newblob=GetBlob(temp))!=NULL) {
	m_blobs[m_blobs.size()-1].push_back(newblob);
      }
    }
    m_filename.push_back("qqb_to_gg__grid");
//     m_create.push_back(false);
    if (mode==EXTRAXS::QCD_Processes::qqbgg) break;
  case EXTRAXS::QCD_Processes::ggqqb:
    m_blobs.push_back(ATOOLS::Blob_List(0));
    temp[1]=temp[0]=ATOOLS::Flavour((ATOOLS::kf::code)21);
    for (i=1;i<=m_nflavour;++i) {
      temp[2]=ATOOLS::Flavour((ATOOLS::kf::code)i);
      temp[3]=ATOOLS::Flavour((ATOOLS::kf::code)i).Bar();
      if ((newblob=GetBlob(temp))!=NULL) {
	m_blobs[m_blobs.size()-1].push_back(newblob);
      }
    }
    m_filename.push_back("gg_to_qqb__grid");
//     m_create.push_back(false);
    if (mode==EXTRAXS::QCD_Processes::ggqqb) break;
  case EXTRAXS::QCD_Processes::qgqg:
    m_blobs.push_back(ATOOLS::Blob_List(0));
    temp[2]=temp[0]=ATOOLS::Flavour((ATOOLS::kf::code)21);
    for (i=1;i<=m_nflavour;++i) {
      temp[3]=temp[1]=ATOOLS::Flavour((ATOOLS::kf::code)i);
      if ((newblob=GetBlob(temp))!=NULL) {
	m_blobs[m_blobs.size()-1].push_back(newblob);
      }
      temp[3]=temp[1]=ATOOLS::Flavour((ATOOLS::kf::code)i).Bar();
      if ((newblob=GetBlob(temp))!=NULL) {
	m_blobs[m_blobs.size()-1].push_back(newblob);
      }
    }
    m_filename.push_back("qg_to_qg__grid");
//     m_create.push_back(false);
    if (mode==EXTRAXS::QCD_Processes::qgqg) break;
  case EXTRAXS::QCD_Processes::q1q2q1q2:
    m_blobs.push_back(ATOOLS::Blob_List(0));
    for (i=1;i<=m_nflavour;++i) {
      for (j=i+1;j<=m_nflavour;++j) {
	temp[2]=temp[0]=ATOOLS::Flavour((ATOOLS::kf::code)i);
	temp[3]=temp[1]=ATOOLS::Flavour((ATOOLS::kf::code)j);
	if ((newblob=GetBlob(temp))!=NULL) {
	  m_blobs[m_blobs.size()-1].push_back(newblob);
	}
	temp[2]=temp[0]=ATOOLS::Flavour((ATOOLS::kf::code)i).Bar();
	temp[3]=temp[1]=ATOOLS::Flavour((ATOOLS::kf::code)j).Bar();
	if ((newblob=GetBlob(temp))!=NULL) {
	  m_blobs[m_blobs.size()-1].push_back(newblob);
	}
      }
    }
  case EXTRAXS::QCD_Processes::q1q2bq1q2b:
    if (mode==EXTRAXS::QCD_Processes::q1q2bq1q2b) m_blobs.push_back(ATOOLS::Blob_List(0));
    for (i=1;i<=m_nflavour;++i) {
      for (j=i+1;j<=m_nflavour;++j) {
	temp[2]=temp[0]=ATOOLS::Flavour((ATOOLS::kf::code)i);
	temp[3]=temp[1]=ATOOLS::Flavour((ATOOLS::kf::code)j).Bar();
	if ((newblob=GetBlob(temp))!=NULL) {
	  m_blobs[m_blobs.size()-1].push_back(newblob);
	}
	temp[2]=temp[0]=ATOOLS::Flavour((ATOOLS::kf::code)i).Bar();
	temp[3]=temp[1]=ATOOLS::Flavour((ATOOLS::kf::code)j);
	if ((newblob=GetBlob(temp))!=NULL) {
	  m_blobs[m_blobs.size()-1].push_back(newblob);
	}
      }
    }
    if (mode==EXTRAXS::QCD_Processes::q1q2bq1q2b) {
      m_filename.push_back("q1q2b_to_q1q2b__grid");
//       m_create.push_back(false);
      break;
    }
    m_filename.push_back("q1q2_to_q1q2__grid");
//     m_create.push_back(false);
    if (mode==EXTRAXS::QCD_Processes::q1q2q1q2) break;
  case EXTRAXS::QCD_Processes::q1q1q1q1:
    m_blobs.push_back(ATOOLS::Blob_List(0));
    for (i=1;i<=m_nflavour;++i) {
      temp[3]=temp[2]=temp[1]=temp[0]=ATOOLS::Flavour((ATOOLS::kf::code)i);
      if ((newblob=GetBlob(temp))!=NULL) {
	m_blobs[m_blobs.size()-1].push_back(newblob);
      }
      temp[3]=temp[2]=temp[1]=temp[0]=ATOOLS::Flavour((ATOOLS::kf::code)i).Bar();
      if ((newblob=GetBlob(temp))!=NULL) {
	m_blobs[m_blobs.size()-1].push_back(newblob);
      }
    }
    m_filename.push_back("q1q1_to_q1q1__grid");
//     m_create.push_back(false);
    if (mode==EXTRAXS::QCD_Processes::q1q1q1q1) break;
  case EXTRAXS::QCD_Processes::q1q1bq1q1b:
    m_blobs.push_back(ATOOLS::Blob_List(0));
    for (i=1;i<=m_nflavour;++i) {
      temp[2]=temp[0]=ATOOLS::Flavour((ATOOLS::kf::code)i);
      temp[3]=temp[1]=ATOOLS::Flavour((ATOOLS::kf::code)i).Bar();
      if ((newblob=GetBlob(temp))!=NULL) {
	m_blobs[m_blobs.size()-1].push_back(newblob);
      }
    }
    m_filename.push_back("q1q1b_to_q1q1b__grid");
//     m_create.push_back(false);
    if (mode==EXTRAXS::QCD_Processes::q1q1bq1q1b) break;
  case EXTRAXS::QCD_Processes::q1q1bq2q2b:
    m_blobs.push_back(ATOOLS::Blob_List(0));
    for (i=1;i<=m_nflavour;++i) {
      for (j=i+1;j<=m_nflavour;++j) {
	temp[0]=ATOOLS::Flavour((ATOOLS::kf::code)i);
	temp[1]=ATOOLS::Flavour((ATOOLS::kf::code)i).Bar();
	temp[2]=ATOOLS::Flavour((ATOOLS::kf::code)j);
	temp[3]=ATOOLS::Flavour((ATOOLS::kf::code)j).Bar();
	if ((newblob=GetBlob(temp))!=NULL) {
	  m_blobs[m_blobs.size()-1].push_back(newblob);
	}
      }
    }
    m_filename.push_back("q1q1b_to_q2q2b__grid");
//     m_create.push_back(false);
    if (mode==EXTRAXS::QCD_Processes::q1q1bq2q2b) break;
    break;
  case EXTRAXS::QCD_Processes::Unknown:
  case EXTRAXS::QCD_Processes::None:
    ATOOLS::msg.Error()<<"Simple_Chain::FillMode("<<mode<<"): "
		       <<"Wrong parameter. Abort."<<std::endl;
    break;
  }
}

bool Simple_Chain::ReadInData()
{
  PROFILE_HERE;
  ATOOLS::Data_Reader *reader = new ATOOLS::Data_Reader("=",";","!");
  reader->SetInputPath(InputPath());
  reader->SetInputFile(InputFile());
  reader->SetMatrixType(reader->MTransposed);
  if (!reader->ReadFromFile(m_scalescheme,"SCALE_SCHEME")) m_scalescheme=11;
  if (!reader->ReadFromFile(m_kfactorscheme,"K_FACTOR_SCHEME")) m_kfactorscheme=1;
  if (!reader->ReadFromFile(m_nflavour,"N_FLAVOUR")) m_nflavour=3;
  if (!reader->ReadFromFile(m_error,"ERROR")) m_error=1.e-2;
  std::string outputpath, help[2];
  MyStrStream converter;
  converter<<ATOOLS::rpa.gen.Bunch(0);
  converter>>help[0];
  converter.clear();
  converter<<ATOOLS::rpa.gen.Bunch(1);
  converter>>help[1];
  outputpath=std::string("MI-Grid__")+help[0]+std::string("_")+help[1]+
    std::string("__")+ATOOLS::ToString(ATOOLS::rpa.gen.Ecms())+std::string("_GeV__")+
    p_isr->PDF(0)->Type()+std::string("_")+p_isr->PDF(1)->Type()+std::string("__as_")+
    ATOOLS::ToString(static_cast<MODEL::Running_AlphaS*>
		     (p_model->GetScalarFunction("alpha_S"))->Order())+std::string("/");
  SetOutputPath(OutputPath()+outputpath);
  std::vector<std::string> comments;
  comments.push_back("->");
  comments.push_back("FOR");
  comments.push_back("INTO");
  comments.push_back("IN");
  reader->AddIgnore(comments);
  std::vector<std::vector<std::string> > temp;
  reader->MatrixFromFile(temp,"CREATE GRID");
  for (unsigned int i=0;i<temp.size();++i) {
    if (temp[i].size()>3) {
      if ((temp[i][0]==std::string("all"))&&
	  (temp[i][1]==std::string("QCD"))&&
	  (temp[i][2]==std::string("2"))&&
	  (temp[i][3]==std::string("2"))) {
	FillMode(EXTRAXS::QCD_Processes::All);
      }
      else {
	ATOOLS::Flavour flavour[4];
	int current;
	bool success=true;
	for (unsigned int j=0;j<4;++j) {
	  flavour[j]=ATOOLS::Flavour(ATOOLS::kf_table.FromString(temp[i][j]));
	  if (flavour[j].Kfcode()==ATOOLS::kf::none) {
	    reader->ReadFromString(current,"",temp[i][j]);
	    flavour[j]=ATOOLS::Flavour((kf::code)abs(current));
	    if (current<0) flavour[j]=flavour[j].Bar();
	    if (flavour[j].Kfcode()==ATOOLS::kf::none) success=false;
	  }
	}
	if (!success) continue;
	ATOOLS::Blob *newblob;
	if ((newblob=GetBlob(flavour))!=NULL) {
	  m_blobs.push_back(ATOOLS::Blob_List(0));
	  m_blobs[m_blobs.size()-1].push_back(newblob);
	  std::string filename=std::string();
	  for (int j=0;j<2;++j) {
	    if (flavour[j].IsAnti()) filename+=std::string("-");
	    filename+=ATOOLS::ToString((int)flavour[j].Kfcode())+std::string("_");
	  }
	  filename+=std::string("to_");
	  for (int j=2;j<4;++j) {
	    if (flavour[j].IsAnti()) filename+=std::string("-");
	    filename+=ATOOLS::ToString((int)flavour[j].Kfcode())+std::string("_");
	  }
	  filename+=std::string("_grid");
	  if (temp[i].size()>4) m_filename.push_back(temp[i][4]);
	  else m_filename.push_back(filename);
// 	  m_create.push_back(false);
	}
      }
    }
    else if (temp[i].size()>1) {
      if ((temp[i][0]==std::string("all"))&&
	  (temp[i][1]==std::string("gggg"))) {
	FillMode(EXTRAXS::QCD_Processes::gggg);
      }
      if ((temp[i][0]==std::string("all"))&&
	  (temp[i][1]==std::string("qqbgg"))) {
	FillMode(EXTRAXS::QCD_Processes::qqbgg);
      }
      if ((temp[i][0]==std::string("all"))&&
	  (temp[i][1]==std::string("ggqqb"))) {
	FillMode(EXTRAXS::QCD_Processes::ggqqb);
      }
      if ((temp[i][0]==std::string("all"))&&
	  (temp[i][1]==std::string("qgqg"))) {
	FillMode(EXTRAXS::QCD_Processes::qgqg);
      }
      if ((temp[i][0]==std::string("all"))&&
	  (temp[i][1]==std::string("q1q1q1q1"))) {
	FillMode(EXTRAXS::QCD_Processes::q1q1q1q1);
      }
      if ((temp[i][0]==std::string("all"))&&
	  (temp[i][1]==std::string("q1q2q1q2"))) {
	FillMode(EXTRAXS::QCD_Processes::q1q2q1q2);
      }
      if ((temp[i][0]==std::string("all"))&&
	  (temp[i][1]==std::string("q1q1bq1q1b"))) {
	FillMode(EXTRAXS::QCD_Processes::q1q1bq1q1b);
      }
      if ((temp[i][0]==std::string("all"))&&
	  (temp[i][1]==std::string("q1q2bq1q2b"))) {
	FillMode(EXTRAXS::QCD_Processes::q1q2bq1q2b);
      }
      if ((temp[i][0]==std::string("all"))&&
	  (temp[i][1]==std::string("q1q1bq2q2b"))) {
	FillMode(EXTRAXS::QCD_Processes::q1q1bq2q2b);
      }
    }
  }
  delete reader;
  return (bool)m_blobs.size();
}

bool Simple_Chain::CreateGrid(ATOOLS::Blob_List& bloblist,std::string& filename)
{
  PROFILE_HERE;
  if (!m_external) {
    p_environment = new AMEGIC::Environment(InputPath(),InputFile(2));
    p_environment->InitializeTheEnvironment();
    p_model=p_environment->Model();
    p_beam=p_environment->BeamSpectraHandler();
    p_isr=p_environment->ISRHandler();
  }
  p_processes = new EXTRAXS::Simple_XS(InputPath(),InputFile(1),p_model);
  if (p_processes->Size()>0) {
    ATOOLS::msg.Tracking()<<"Simple_Chain::CreateGrid(..): "
			  <<"Found an initialized process group."<<std::endl
			  <<"   Empty group and start with QCD_Processes."<<std::endl;
    p_processes->Clear();
  }
  p_processes->InitializeProcesses(p_beam,p_isr,false);  
  p_processes->SetScaleScheme(m_scalescheme);
  p_processes->SetKFactorScheme(m_kfactorscheme);
  ATOOLS::Flavour flavour[4];
  flavour[0]=flavour[1]=flavour[2]=flavour[3]=ATOOLS::kf::jet;
  EXTRAXS::QCD_Processes *group;
  group = new EXTRAXS::QCD_Processes(p_isr,p_beam,flavour,p_processes->SelectorData(),
				     p_processes->ScaleScheme(),p_processes->KFactorScheme(),
				     p_processes->ScaleFactor(),false);
  for (Blob_Iterator bit=bloblist.begin();bit!=bloblist.end();++bit) {
    for (size_t i=0;i<group->NIn();++i) flavour[i]=(*bit)->InParticle(i)->Flav();
    for (size_t j=0;j<group->NOut();++j) flavour[group->NIn()+j]=(*bit)->OutParticle(j)->Flav();
    EXTRAXS::Single_XS *newxs = group->XSSelector()->GetXS(group->NIn(),group->NOut(),
							   flavour,p_isr->KMROn());
    if (newxs==NULL) {
      ATOOLS::msg.Error()<<"Simple_Chain::CreateGrid(..): "
			 <<"Did not find any process! Abort calculation."<<std::endl;
      delete p_processes;
      p_processes=NULL;
      if (!m_external) {
	delete p_environment;
	p_environment=NULL;
	p_model=NULL;
	p_beam=NULL;
	p_isr=NULL;
      }
      return false;
    }
    newxs->SetScaleScheme(m_scalescheme);
    newxs->SetKFactorScheme(m_kfactorscheme);
    group->Add(newxs);
  }
  group->SetScaleScheme(m_scalescheme);
  group->SetKFactorScheme(m_kfactorscheme);
  group->PSHandler(false)->SetError(m_error);
  p_processes->PushBack(group);
  m_comments.clear();
  std::string processname=filename.substr(0,filename.length()-4);
  size_t pos=0;
  while ((pos=processname.find("_"))!=std::string::npos) processname[pos]=' ';
  m_comments.push_back(std::string("processes : ")+processname);
  GridHandlerVector gridhandler=GridHandlerVector(1+group->Size());
  for (unsigned int i=0;i<gridhandler.size();++i) gridhandler[i] = new GridHandlerType();
  bool read=true;
  gridhandler[0]->Grid()->SetMonotony(GridFunctionType::None);
  read=read&&gridhandler[0]->ReadIn(ATOOLS::Type::TFStream,OutputPath()+filename+m_xsextension);
  for (unsigned int i=1;i<gridhandler.size();++i) {
    gridhandler[i]->Grid()->SetMonotony(GridFunctionType::None);
    read=read&&gridhandler[i]->ReadIn(ATOOLS::Type::TFStream,
				      OutputPath()+(*group)[i-1]->Name()+m_maxextension);
  }
  if (!read) for (unsigned int i=0;i<gridhandler.size();++i) {
    gridhandler[i]->Grid()->Clear();
    p_gridcreator = new GridCreatorType(gridhandler,group);
    p_gridcreator->ReadInArguments(InputFile(),InputPath());
    if (mkdir(OutputPath().c_str(),448)==0) {
      ATOOLS::msg.Out()<<"Simple_Chain::CreateGrid(..): "
		       <<"Created output directory "<<OutputPath()<<"."<<std::endl;
    }
    if (mkdir((OutputPath()+filename+m_mcextension+std::string("/")).c_str(),448)==0) {
      ATOOLS::msg.Out()<<"Simple_Chain::CreateGrid(..): Created output directory "
		       <<OutputPath()+filename+m_mcextension+std::string("/")<<"."<<std::endl;
    }
    p_gridcreator->SetXSExtension(m_xsextension);
    p_gridcreator->SetMaxExtension(m_maxextension);
    p_gridcreator->SetMCExtension(m_mcextension);
    p_gridcreator->SetOutputPath(OutputPath());
    p_gridcreator->SetOutputFile(filename);
    ATOOLS::Exception_Handler::AddTerminatorObject(this);
    p_gridcreator->CreateGrid();
    p_gridcreator->WriteOutGrid(m_comments);
    ATOOLS::Exception_Handler::RemoveTerminatorObject(this);
    delete p_gridcreator;
  }
  m_comments.clear();
  m_differential.push_back(new GridFunctionType(*gridhandler[0]->Grid()));
  delete gridhandler[0];
  for (unsigned int i=1;i<gridhandler.size();++i) {
    m_maxima[(*group)[i-1]->Name()] = new GridFunctionType(*gridhandler[i]->Grid());
    delete gridhandler[i];
  }
  delete p_processes;
  p_processes=NULL;
  if (!m_external) {
    delete p_environment;
    p_environment=NULL;
    p_model=NULL;
    p_beam=NULL;
    p_isr=NULL;
  }
  return true;
}

bool Simple_Chain::InitializeBlobList()
{  
  PROFILE_HERE;
  if (!m_external) {
    p_environment = new AMEGIC::Environment(InputPath(),InputFile(2));
    p_environment->InitializeTheEnvironment();
    p_model=p_environment->Model();
    p_beam=p_environment->BeamSpectraHandler();
    p_isr=p_environment->ISRHandler();
  }
  SetStart(sqrt(p_isr->SprimeMax()),1);
  SetStop(sqrt(p_isr->SprimeMin()),1);
  p_processes = new EXTRAXS::Simple_XS(InputPath(),InputFile(1),p_model);
  if (p_processes->Size()>0) {
    ATOOLS::msg.Tracking()<<"Simple_Chain::InitializeBlobList(): "
			  <<"Found an initialized process group."<<std::endl
			  <<"   Empty group and start with QCD_Processes."<<std::endl;
    p_processes->Clear();
  }
  p_processes->InitializeProcesses(p_beam,p_isr,false);  
  p_processes->SetScaleScheme(m_scalescheme);
  p_processes->SetKFactorScheme(m_kfactorscheme);
  std::vector<EXTRAXS::QCD_Processes*> group=std::vector<EXTRAXS::QCD_Processes*>(m_blobs.size());
  ATOOLS::Flavour flavour[4];
  ATOOLS::Flow::ResetCounter();
  for (unsigned int i=0;i<m_blobs.size();++i) {
    flavour[0]=flavour[1]=flavour[2]=flavour[3]=ATOOLS::kf::jet;
    group[i]=new EXTRAXS::QCD_Processes(p_isr,p_beam,flavour,p_processes->SelectorData(),
					p_processes->ScaleScheme(),p_processes->KFactorScheme(),
					p_processes->ScaleFactor(),false);
    for (Blob_Iterator bit=m_blobs[i].begin();bit!=m_blobs[i].end();++bit) {
      for (size_t j=0;j<group[i]->NIn();++j) flavour[j]=(*bit)->InParticle(j)->Flav();
      for (size_t j=0;j<group[i]->NOut();++j) flavour[group[i]->NIn()+j]=(*bit)->OutParticle(j)->Flav();
      EXTRAXS::Single_XS *newxs = group[i]->XSSelector()->GetXS(group[i]->NIn(),group[i]->NOut(),
								flavour,p_isr->KMROn());
      if (newxs==NULL) {
	ATOOLS::msg.Error()<<"Simple_Chain::InitializeBlobList(): "
			   <<"Did not find any process! Abort."<<std::endl;
	delete group[i];
 	delete p_processes;
 	p_processes=NULL;
	if (!m_external) {
	  delete p_environment;
	  p_environment=NULL;
	  p_model=NULL;
	  p_beam=NULL;
	  p_isr=NULL;
	}
	return false;
      }
      newxs->SetScaleScheme(m_scalescheme);
      newxs->SetKFactorScheme(m_kfactorscheme);
      group[i]->Add(newxs);
    }
    group[i]->SetScaleScheme(m_scalescheme);
    group[i]->SetKFactorScheme(m_kfactorscheme);
    p_processes->PushBack(group[i]);
  }
  p_fsrinterface = new FSRChannel(2,2,flavour,p_total->XAxis()->Variable());
  p_fsrinterface->SetAlpha(1.0);
  p_fsrinterface->SetAlphaSave(1.0);
  for (unsigned int i=0;i<m_blobs.size();++i) {
    group[i]->SetFSRInterface(p_fsrinterface);
    group[i]->SetFSRMode(2);
    group[i]->CreateFSRChannels();
    group[i]->InitIntegrators();
    group[i]->PSHandler(false)->ReadIn(OutputPath()+m_filename[i]+m_mcextension,16);
  }
#ifdef USING__Sherpa
  p_mehandler = new SHERPA::Matrix_Element_Handler();
  p_mehandler->SetXS(p_processes);
#endif
  return true;
}

void Simple_Chain::CalculateSigmaND()
{
  double eps=0.0808, eta=-0.4525, X=21.70, Y=56.08, b=2.3;
  if (p_isr->Flav(0).IsAnti()^p_isr->Flav(1).IsAnti()) Y=98.39;
  double s=ATOOLS::sqr(ATOOLS::rpa.gen.Ecms());
  double mp=ATOOLS::Flavour(ATOOLS::kf::p_plus).PSMass();
  double mpi=ATOOLS::Flavour(ATOOLS::kf::pi).PSMass();
  double ap=0.25, s0=8.0, y0=log(s/(mp*mp));
  double M1res=2.0, M2res=2.0, cres=2.0;
  double M1min=mp+2.0*mpi, M2min=mp+2.0*mpi;
  double ymin=4.0*log(1.0+2.0*mpi/mp);
  double MmaxAX2=0.213*s;
  double Del0=3.2-9.0/log(s)+17.4/ATOOLS::sqr(log(s));
  double MmaxXX2=s*(0.07-0.44/log(s)+1.36/ATOOLS::sqr(log(s)));
  double BAX=-0.47+150.0/s;
  double BXX=-1.05+40.0/sqrt(s)+8000.0/(s*s);
  double JAX=0.5/ap*log((b+ap*log(s/(M2min*M2min)))/(b+ap*log(s/MmaxAX2)));
  JAX+=0.5*cres/(b+ap*log(s/(M2res*M2min))+BAX)*log(1.0+M2res*M2res/(M2min*M2min));
  double s1, s2, s3, JXX=0.5/ap*((y0-ymin)*(log((y0-ymin)/Del0)-1.0)+Del0);
  JXX+=s1=0.5*cres/ap*log(log(s*s0/(M1min*M1min*M2res*M2min))/
		       log(s*s0/(MmaxXX2*M2res*M2min)))*log(1.0+M2res*M2res/(M2min*M2min));
  JXX+=s2=0.5*cres/ap*log(log(s*s0/(M2min*M2min*M1res*M1min))/
		       log(s*s0/(MmaxXX2*M1res*M1min)))*log(1.0+M1res*M1res/(M1min*M1min));
  JXX+=s3=cres*cres/(2.0*ap*log(s*s0/(M1res*M2res*M1min*M2min))+BXX)*
    log(1.0+M1res*M1res/(M1min*M1min))*log(1.0+M2res*M2res/(M2min*M2min));
  double xstot=X*pow(s,eps)+Y*pow(s,eta);
  double xsel=0.0511*xstot*xstot/(4*(b+pow(s,eps))-4.2);
  double xssd=0.0336*X*sqrt(X)*JAX;
  double xsdd=0.0084*X*JXX;
  SetNorm(xstot-xsel-xssd-xsdd);
}

bool Simple_Chain::CalculateTotal()
{
  PROFILE_HERE;
  if (m_differential.size()==0) return false;
  GridFunctionType *differential;
  differential = new GridHandlerType::GridFunctionType();
  differential->SetMonotony(differential->None);
  differential->XAxis()->SetVariable(m_differential[0]->XAxis()->Variable());
  differential->YAxis()->SetVariable(m_differential[0]->YAxis()->Variable());
  differential->XAxis()->SetScaling(m_differential[0]->XAxis()->Scaling()->Name());
  differential->YAxis()->SetScaling(m_differential[0]->YAxis()->Scaling()->Name());
  std::vector<GridFunctionType*>::iterator diffit;
  std::set<GridArgumentType> *xvalue = new std::set<GridArgumentType>();
  for (diffit=m_differential.begin();diffit!=m_differential.end();++diffit) {
    for (unsigned int i=0;i<(*diffit)->XDataSize();++i) {
      GridArgumentType x=(*diffit)->XYData(i).first;
      if (xvalue->find(x)==xvalue->end()) xvalue->insert(x);
    }
  }
  for (std::set<GridArgumentType>::iterator xit=xvalue->begin();xit!=xvalue->end();++xit) {
    GridResultType y=(GridResultType)0.0;
    for (diffit=m_differential.begin();diffit!=m_differential.end();++diffit) {
      y+=(*diffit)->Y(*xit,(*diffit)->Interpolation);
    }
    differential->AddPoint(*xit,y);
  }
  delete xvalue;
#ifdef DEBUG__Simple_Chain
  std::vector<std::string> comments;
  comments.push_back("  Differential XS   "); 
  GridHandlerType *gridhandler = new GridHandlerType(differential);
  GridCreatorBaseType *gridcreator = new GridCreatorBaseType();
  gridcreator->SetOutputPath(OutputPath());
  gridcreator->SetOutputFile(differentialfile);
  gridcreator->WriteSingleGrid(gridhandler,comments);
  delete gridhandler;
#endif
  p_total = differential->IntegralY(0.0,0.0,ATOOLS::nullstring,ATOOLS::nullstring,false);
  m_sigmahard=p_total->YMax();
  p_total->ScaleY(1.0/m_norm);
#ifdef DEBUG__Simple_Chain
  comments.clear();
  comments.push_back("   Integrated XS    "); 
  GridFunctionType *total = differential->IntegralY(0.0,0.0,"Id","Id",false);
  gridhandler = new GridHandlerType(total);
  gridcreator->SetOutputFile(integralfile);
  gridcreator->WriteSingleGrid(gridhandler,comments);
  delete gridhandler;
  delete total;
  comments.clear();
  comments.push_back("   Normalized XS    "); 
  gridhandler = new GridHandlerType(p_total);
  gridcreator->SetOutputFile(normalizedfile);
  gridcreator->WriteSingleGrid(gridhandler,comments);
  delete gridcreator;
  delete gridhandler;
#endif
  delete differential;
  return true;
}

bool Simple_Chain::Initialize()
{
  PROFILE_HERE;
  if (!CheckInputPath()) return false;
  if (!CheckInputFile()) return false;
  CleanUp();
  CalculateSigmaND();
  ATOOLS::Data_Reader *reader = new ATOOLS::Data_Reader("=",";","!");
  reader->SetInputPath(InputPath());
  if (!m_external) {
    ATOOLS::ParticleInit(InputPath());
    reader->SetInputFile(InputFile());
    std::string initfile=std::string("Run.dat");
    reader->ReadFromFile(initfile,"ENVIRONMENT");
    SetInputFile(initfile,2);
    ATOOLS::rpa.Init(InputPath(),InputFile(2));
  }
  if (!ReadInData()) return false;
  reader->SetInputFile(InputFile());
  std::string xsfile=std::string("XS.dat");
  reader->ReadFromFile(xsfile,"XS_FILE");
  SetInputFile(xsfile,1);
  double stop;
  if (!reader->ReadFromFile(stop,"EVENT_X_MIN")) stop=Stop(0);
  delete reader;
  for (unsigned int i=0;i<m_blobs.size();++i) {
    if (!CreateGrid(m_blobs[i],m_filename[i])) {
      CleanUp();
      throw(ATOOLS::Exception(ATOOLS::ex::critical_error,
			      std::string("Grid creation failed for ")+OutputPath()+m_filename[i],
			      "Simple_Chain","Initialize"));
    }
  }
  if (!CalculateTotal()) {
    throw(ATOOLS::Exception(ATOOLS::ex::critical_error,"Determination of \\sigma_{tot} failed.",
			    "Simple_Chain","Initialize"));
  }
  SetStart(p_total->XMax(),0);
  SetStop(ATOOLS::Max(p_total->XMin(),stop),0);
  if (!InitializeBlobList()) {
    throw(ATOOLS::Exception(ATOOLS::ex::critical_error,"Cannot initialize selected processes.",
			    "Simple_Chain","Initialize"));
  }  
  return true;
}

bool Simple_Chain::FillBlob(ATOOLS::Blob *blob)
{
  PROFILE_HERE;
  m_filledblob=false;
  if (p_processes==NULL) {
    ATOOLS::msg.Error()<<"Simple_Chain::FillBlob(..): "
		       <<"Processes are not initialized! Abort."<<std::endl;
    return false;
  }
  if (blob==NULL) {
    ATOOLS::msg.Error()<<"Simple_Chain::FillBlob(..): "
		       <<"Blob is not initialized! Abort."<<std::endl;
    return false;
  }
  blob->DeleteOwnedParticles();
  if (m_selected<(unsigned int)p_processes->Size()) {
#ifdef DEBUG__Simple_Chain
    // std::cout<<"Simple_Chain::FillBlob(..): Generating one event."<<std::endl;
#endif
    size_t trials=0;
    bool test=false;
    while (!test && trials<m_maxtrials) {
      ++trials;
      dynamic_cast<EXTRAXS::XS_Group*>((*p_processes)[m_selected])->SetAtoms(true);
      (*p_processes)[m_selected]->SelectOne();
      if (!(*p_processes)[m_selected]->OneEvent(-1.,1)) {
	ATOOLS::msg.Error()<<"Simple_Chain::FillBlob(): "
			   <<"Could not select any process! Retry."<<std::endl;
	return false;
      }
#ifdef DEBUG__Simple_Chain
      // std::cout<<"   Completed one event."<<std::endl;
#endif
      p_xs=static_cast<EXTRAXS::XS_Base*>((*p_processes)[m_selected]->Selected());
      ATOOLS::Vec4D ptot;
      double x[2];
      test=true;
      for (size_t j=0;j<p_xs->NIn();++j) {
	ptot+=p_xs->Momenta()[j];
	x[j]=2.0*p_xs->Momenta()[j][0]/ATOOLS::rpa.gen.Ecms();
	if (m_last[j+2]<=x[j]) test=false;
      }
      if (test==false) continue;
      m_last[1]-=sqrt(ptot.Abs2());
      m_last[2]-=x[0];
      m_last[3]-=x[1];
      if (m_last[1]<=0.0) return true;
      p_xs->SetColours(p_xs->Momenta());
      ATOOLS::Particle *particle;
      for (size_t j=0;j<p_xs->NIn();++j) {
	particle = new ATOOLS::Particle(0,p_xs->Flavours()[j]);
	particle->SetMomentum(p_xs->Momenta()[j]);
 	particle->SetFlow(1,p_xs->Colours()[j][0]);
 	particle->SetFlow(2,p_xs->Colours()[j][1]);
	particle->SetStatus(1);
	blob->AddToInParticles(particle);
      }
      for (size_t j=0;j<p_xs->NOut();++j) {
	particle = new ATOOLS::Particle(0,p_xs->Flavours()[p_xs->NIn()+j]);
	particle->SetMomentum(p_xs->Momenta()[p_xs->NIn()+j]);
 	particle->SetFlow(1,p_xs->Colours()[p_xs->NIn()+j][0]);
 	particle->SetFlow(2,p_xs->Colours()[p_xs->NIn()+j][1]);
	particle->SetStatus(1);
	blob->AddToOutParticles(particle);
      }
      m_filledblob=true;
      return true;
    }
    ATOOLS::msg.Tracking()<<"Simple_Chain::FillBlob(..): "
			  <<"Cannot create momentum configuration."<<std::endl;
  }
  return false;
}

bool Simple_Chain::DiceProcess()
{
  PROFILE_HERE;
  if (m_differential.size()==0) return false;
  if (m_dicedparameter) m_dicedparameter=false;
  else {
    m_dicedprocess=false;
    return true;
  }
  p_fsrinterface->SetValue(m_last[0]);
  int criterion=ATOOLS::Variable::TypeToSelectorID(m_differential[0]->XAxis()->Variable().Type());
  ATOOLS::Mom_Data initialdata=p_processes->SelectorData()->RemoveData(criterion);
  p_processes->SelectorData()->AddData(criterion,initialdata.flavs,initialdata.help,
 				       m_last[0],initialdata.max);
  p_processes->ResetSelector(p_processes->SelectorData());
  Data_To_Function<GridResultType,unsigned int> sorter;
  sorter.SetMonotony(sorter.None);
  GridResultType cur, norm=(GridResultType)0.0;
  for (unsigned int i=0;i<m_differential.size();++i) {
    cur=(*m_differential[i])(m_last[0]);
    sorter.AddPoint(cur,i);
    norm+=cur;
  }
  double rannr=ATOOLS::ran.Get();
  cur=(GridResultType)0.0;
  for (int i=sorter.XDataSize()-1;i>=0;--i) {
    if ((cur+=sorter.XData(i)/norm)>rannr) {
      m_selected=sorter.XYData(i).second;
      p_processes->SetSelected((*p_processes)[m_selected]);
      if (m_last[1]<(*p_processes)[m_selected]->ISR()->SprimeMin()) {
	ATOOLS::msg.Error()<<"Simple_Chain::DiceProcess(): s' out of bounds: "
			   <<m_last[1]<<" vs. "<<(*p_processes)[m_selected]->ISR()->SprimeMin()
			   <<std::endl<<"   Cannot create any process."<<std::endl;
	return false;
      }
      PDF::ISR_Handler *isr=(*p_processes)[m_selected]->ISR();
      double sprimemin=isr->SprimeMin(), sprimemax=isr->SprimeMax();
      isr->SetSprimeMax(m_last[1]*m_last[1]);
      isr->SetSprimeMin(4.0*m_last[0]*m_last[0]);
      SetMaximum((*p_processes)[m_selected]);
      FillBlob(p_blob);
      isr->SetSprimeMax(sprimemax);
      isr->SetSprimeMin(sprimemin);
      m_dicedprocess=true;
      return m_filledblob;
    }
  }
  throw(ATOOLS::Exception(ATOOLS::ex::critical_error,"Internal Error. Could not select any process.",
			  "Simple_Chain","DiceProcess"));
  return false;
}

bool Simple_Chain::DiceOrderingParameter()
{ 
  PROFILE_HERE;
  if (m_last[0]<=m_stop[0]) {
    ATOOLS::msg.Tracking()<<"Simple_Chain::DiceOrderingParameter(): "
			  <<"Ordering parameter exceeded allowed range."<<std::endl
			  <<"   Cannot proceed in blob creation."<<std::endl;
    m_stophard=true;
    return false;
  }
  m_last[0]=(*p_total)[(*p_total)(m_last[0])-log(ATOOLS::ran.Get())]; 
  if (m_last[0]<=m_stop[0]) { 
    m_dicedparameter=false;
    m_stophard=true;
    return true;
  }
  m_dicedparameter=true;
  m_stophard=false;
  return true;
}

void Simple_Chain::SetMaximum(EXTRAXS::XS_Base *process) 
{
  process->SetTotalXS(0.);
  if ((*process)[0]==NULL) process->SetMax((*m_maxima[process->Name()])(m_last[0]),true);
  else {
    EXTRAXS::XS_Group *node=dynamic_cast<EXTRAXS::XS_Group*>(process);
    for (size_t i=0;i<node->Size();++i) SetMaximum((*node)[i]);
    node->SetMax(0.,false);
    node->SetTotalXS(1.);
  }
}

void Simple_Chain::Reset()
{
  for (unsigned int i=0;i<4;++i) m_last[i]=m_start[i];
}

void Simple_Chain::Update(const MI_Base *mibase)
{
  return;
}

void Simple_Chain::PrepareTerminate() 
{
  p_gridcreator->WriteOutGrid(m_comments);
}

