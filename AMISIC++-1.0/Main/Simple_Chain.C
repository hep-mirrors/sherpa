#include "Simple_Chain.H"

#include "Particle.H"
#include "Random.H"
#include "Channel_Elements.H"

#ifdef DEBUG__Simple_Chain
const std::string integralfile=std::string("integral.dat");
#endif

using namespace AMISIC;

Simple_Chain::Simple_Chain():
  MI_Base("Simple Chain",MI_Base::HardEvent,2),
  m_differential(std::vector<GridFunctionType*>(0)),
  p_total(NULL),
  m_environmentfile(std::string("Run.dat")),
  m_xsfile(std::string("XS.dat")),
  m_xsextension("_xs.dat"),
  m_maxextension("_max.dat"),
  p_processes(NULL),
  p_fsrinterface(NULL), 
  p_environment(NULL),
  p_model(NULL),
  p_beam(NULL),
  p_isr(NULL),
  m_external(false) 
{
  SetInputFile("MI.dat");
  SetOutputFile("SC.log");
  SetStart(0.0,0);
  SetStop(1.6,0);
  SetStart(ATOOLS::rpa.gen.Ecms(),1);
  SetStop(0.0,1);
}

Simple_Chain::Simple_Chain(MODEL::Model_Base *_p_model,
			   BEAM::Beam_Spectra_Handler *_p_beam,PDF::ISR_Handler *_p_isr):
  MI_Base("Simple Chain",MI_Base::HardEvent,2),
  m_differential(std::vector<GridFunctionType*>(0)),
  p_total(NULL),
  m_environmentfile(std::string("Run.dat")),
  m_xsfile(std::string("XS.dat")),
  m_xsextension("_xs.dat"),
  m_maxextension("_max.dat"),
  p_processes(NULL),
  p_fsrinterface(NULL), 
  p_environment(NULL),
  p_model(_p_model),
  p_beam(_p_beam),
  p_isr(_p_isr),
  m_external(true)
{
  SetInputFile("MI.dat");
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
  if (p_processes!=NULL) delete p_processes;
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
  m_create.clear();
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

void Simple_Chain::FillMode(EXTRAXS::QCD_Processes_C::Mode mode)
{
  const unsigned int nflavour=4;
  ATOOLS::Flavour temp[4];
  ATOOLS::Blob *newblob;
  unsigned int i, j;
  switch (mode) {
  case EXTRAXS::QCD_Processes_C::All:
  case EXTRAXS::QCD_Processes_C::gggg:
    m_blobs.push_back(ATOOLS::Blob_List(0));
    for (i=0;i<4;++i) temp[i]=ATOOLS::Flavour((ATOOLS::kf::code)21);
    if ((newblob=GetBlob(temp))!=NULL) {
      m_blobs[m_blobs.size()-1].push_back(newblob);
    }
    m_filename.push_back("gg_to_gg__grid");
    m_processname.push_back("g g -> g g");
    m_create.push_back(false);
    if (mode==EXTRAXS::QCD_Processes_C::gggg) break;
  case EXTRAXS::QCD_Processes_C::qqbgg:
    m_blobs.push_back(ATOOLS::Blob_List(0));
    temp[3]=temp[2]=ATOOLS::Flavour((ATOOLS::kf::code)21);
    for (i=1;i<=nflavour;++i) {
      temp[0]=ATOOLS::Flavour((ATOOLS::kf::code)i);
      temp[1]=ATOOLS::Flavour((ATOOLS::kf::code)i).Bar();
      if ((newblob=GetBlob(temp))!=NULL) {
	m_blobs[m_blobs.size()-1].push_back(newblob);
      }
    }
    m_filename.push_back("qqb_to_gg__grid");
    m_processname.push_back("q qb -> g g");
    m_create.push_back(false);
    if (mode==EXTRAXS::QCD_Processes_C::qqbgg) break;
  case EXTRAXS::QCD_Processes_C::ggqqb:
    m_blobs.push_back(ATOOLS::Blob_List(0));
    temp[1]=temp[0]=ATOOLS::Flavour((ATOOLS::kf::code)21);
    for (i=1;i<=nflavour;++i) {
      temp[2]=ATOOLS::Flavour((ATOOLS::kf::code)i);
      temp[3]=ATOOLS::Flavour((ATOOLS::kf::code)i).Bar();
      if ((newblob=GetBlob(temp))!=NULL) {
	m_blobs[m_blobs.size()-1].push_back(newblob);
      }
    }
    m_filename.push_back("gg_to_qqb__grid");
    m_processname.push_back("g g -> q qb");
    m_create.push_back(false);
    if (mode==EXTRAXS::QCD_Processes_C::ggqqb) break;
  case EXTRAXS::QCD_Processes_C::qgqg:
    m_blobs.push_back(ATOOLS::Blob_List(0));
    temp[2]=temp[0]=ATOOLS::Flavour((ATOOLS::kf::code)21);
    for (i=1;i<=nflavour;++i) {
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
    m_processname.push_back("q g -> q g");
    m_create.push_back(false);
    if (mode==EXTRAXS::QCD_Processes_C::qgqg) break;
  case EXTRAXS::QCD_Processes_C::q1q2q1q2:
    m_blobs.push_back(ATOOLS::Blob_List(0));
    for (i=1;i<=nflavour;++i) {
      for (j=i+1;j<=nflavour;++j) {
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
  case EXTRAXS::QCD_Processes_C::q1q2bq1q2b:
    if (mode==EXTRAXS::QCD_Processes_C::q1q2bq1q2b) m_blobs.push_back(ATOOLS::Blob_List(0));
    for (i=1;i<=nflavour;++i) {
      for (j=i+1;j<=nflavour;++j) {
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
    if (mode==EXTRAXS::QCD_Processes_C::q1q2bq1q2b) {
      m_filename.push_back("q1q2b_to_q1q2b__grid");
      m_processname.push_back("q1 q2b -> q1 q2b");
      m_create.push_back(false);
      break;
    }
    m_filename.push_back("q1q2_to_q1q2__grid");
    m_processname.push_back("q1 q2 -> q1 q2");
    m_create.push_back(false);
    if (mode==EXTRAXS::QCD_Processes_C::q1q2q1q2) break;
  case EXTRAXS::QCD_Processes_C::q1q1q1q1:
    m_blobs.push_back(ATOOLS::Blob_List(0));
    for (i=1;i<=nflavour;++i) {
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
    m_processname.push_back("q1 q1 -> q1 q1");
    m_create.push_back(false);
    if (mode==EXTRAXS::QCD_Processes_C::q1q1q1q1) break;
  case EXTRAXS::QCD_Processes_C::q1q1bq1q1b:
    m_blobs.push_back(ATOOLS::Blob_List(0));
    for (i=1;i<=nflavour;++i) {
      temp[2]=temp[0]=ATOOLS::Flavour((ATOOLS::kf::code)i);
      temp[3]=temp[1]=ATOOLS::Flavour((ATOOLS::kf::code)i).Bar();
      if ((newblob=GetBlob(temp))!=NULL) {
	m_blobs[m_blobs.size()-1].push_back(newblob);
      }
    }
    m_filename.push_back("q1q1b_to_q1q1b__grid");
    m_processname.push_back("q1 q1b -> q1 q1b");
    m_create.push_back(false);
    if (mode==EXTRAXS::QCD_Processes_C::q1q1bq1q1b) break;
  case EXTRAXS::QCD_Processes_C::q1q1bq2q2b:
    m_blobs.push_back(ATOOLS::Blob_List(0));
    for (i=1;i<=nflavour;++i) {
      for (j=i+1;j<=nflavour;++j) {
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
    m_processname.push_back("q1 q1b -> q2 q2b");
    m_create.push_back(false);
    if (mode==EXTRAXS::QCD_Processes_C::q1q1bq2q2b) break;
    break;
  case EXTRAXS::QCD_Processes_C::Unknown:
  case EXTRAXS::QCD_Processes_C::None:
    ATOOLS::msg.Error()<<"Simple_Chain::FillMode("<<mode<<"): "
		       <<"Wrong parameter. Abort."<<std::endl;
    break;
  }
}

bool Simple_Chain::ReadInData()
{
  std::vector<std::string> comments;
  comments.push_back("->");
  comments.push_back("FOR");
  comments.push_back("INTO");
  comments.push_back("IN");
  ATOOLS::Data_Reader *reader = new ATOOLS::Data_Reader("=",";","!");
  reader->SetFileName(m_inputpath+m_inputfile);
  reader->AddIgnore(comments);
  std::string outputpath;
  reader->ReadFromFile(outputpath,"GRID DIRECTORY");
  m_outputpath+=outputpath;
  std::vector<std::vector<std::string> > temp;
  reader->ArrayFromFile(temp,"CREATE GRID",ATOOLS::noinputtag,reader->MTransposed);
  for (unsigned int i=0;i<temp.size();++i) {
    if (temp[i].size()>3) {
      if ((temp[i][0]==std::string("all"))&&
	  (temp[i][1]==std::string("QCD"))&&
	  (temp[i][2]==std::string("2"))&&
	  (temp[i][3]==std::string("2"))) {
	FillMode(EXTRAXS::QCD_Processes_C::All);
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
	  std::string filename=std::string(), processname=std::string();
	  for (int j=0;j<2;++j) {
	    if (flavour[j].IsAnti()) {
	      filename+=std::string("-");
	      processname+=std::string("-");
	    }
	    filename+=ATOOLS::ToString((int)flavour[j].Kfcode())+std::string("_");
	    processname+=ATOOLS::ToString((int)flavour[j].Kfcode())+std::string(" ");
	  }
	  filename+=std::string("to_");
	  processname+=std::string("-> ");
	  for (int j=2;j<4;++j) {
	    if (flavour[j].IsAnti()) {
	      filename+=std::string("-");
	      processname+=std::string("-");
	    }
	    filename+=ATOOLS::ToString((int)flavour[j].Kfcode())+std::string("_");
	    processname+=ATOOLS::ToString((int)flavour[j].Kfcode())+std::string(" ");
	  }
	  filename+=std::string("_grid");
	  if (temp[i].size()>4) m_filename.push_back(temp[i][4]);
	  else m_filename.push_back(filename);
	  m_processname.push_back(processname);
	  m_create.push_back(false);
	}
      }
    }
    else if (temp[i].size()>1) {
      if ((temp[i][0]==std::string("all"))&&
	  (temp[i][1]==std::string("gggg"))) {
	FillMode(EXTRAXS::QCD_Processes_C::gggg);
      }
      if ((temp[i][0]==std::string("all"))&&
	  (temp[i][1]==std::string("qqbgg"))) {
	FillMode(EXTRAXS::QCD_Processes_C::qqbgg);
      }
      if ((temp[i][0]==std::string("all"))&&
	  (temp[i][1]==std::string("ggqqb"))) {
	FillMode(EXTRAXS::QCD_Processes_C::ggqqb);
      }
      if ((temp[i][0]==std::string("all"))&&
	  (temp[i][1]==std::string("qgqg"))) {
	FillMode(EXTRAXS::QCD_Processes_C::qgqg);
      }
      if ((temp[i][0]==std::string("all"))&&
	  (temp[i][1]==std::string("q1q1q1q1"))) {
	FillMode(EXTRAXS::QCD_Processes_C::q1q1q1q1);
      }
      if ((temp[i][0]==std::string("all"))&&
	  (temp[i][1]==std::string("q1q2q1q2"))) {
	FillMode(EXTRAXS::QCD_Processes_C::q1q2q1q2);
      }
      if ((temp[i][0]==std::string("all"))&&
	  (temp[i][1]==std::string("q1q1bq1q1b"))) {
	FillMode(EXTRAXS::QCD_Processes_C::q1q1bq1q1b);
      }
      if ((temp[i][0]==std::string("all"))&&
	  (temp[i][1]==std::string("q1q2bq1q2b"))) {
	FillMode(EXTRAXS::QCD_Processes_C::q1q2bq1q2b);
      }
      if ((temp[i][0]==std::string("all"))&&
	  (temp[i][1]==std::string("q1q1bq2q2b"))) {
	FillMode(EXTRAXS::QCD_Processes_C::q1q1bq2q2b);
      }
    }
  }
  delete reader;
  return (bool)m_blobs.size();
}

bool Simple_Chain::CreateGrid(ATOOLS::Blob_List& bloblist,std::string& filename,std::string& processname)
{
  if (!m_external) {
    p_environment = new AMEGIC::Environment(m_inputpath,m_environmentfile);
    p_environment->InitializeTheEnvironment();
    p_model=p_environment->Model();
    p_beam=p_environment->BeamSpectraHandler();
    p_isr=p_environment->ISRHandler();
  }
  p_processes = new EXTRAXS::SimpleXSecs(m_inputpath,m_xsfile,p_model);
  if (p_processes->Size()>0) {
    ATOOLS::msg.Tracking()<<"Simple_Chain::CreateGrid(..): "
			  <<"Found an initialized process group."<<std::endl
			  <<"   Empty group and start with QCD_Processes_C."<<std::endl;
    p_processes->Clear();
  }
  p_processes->InitializeProcesses(p_beam,p_isr);  
  ATOOLS::Flavour flavour[4];
  flavour[0]=flavour[1]=flavour[2]=flavour[3]=ATOOLS::kf::jet;
  EXTRAXS::QCD_Processes_C *group;
  group = new EXTRAXS::QCD_Processes_C(p_isr,p_beam,flavour,p_processes->SelectorData(),
				       p_processes->ScaleScheme(),p_processes->KFactorScheme(),
				       p_processes->ScaleFactor(),false);
  for (Blob_Iterator bit=bloblist.begin();bit!=bloblist.end();++bit) {
    for (int i=0;i<group->Nin();++i) flavour[i]=(*bit)->InParticle(i)->Flav();
    for (int j=0;j<group->Nout();++j) flavour[group->Nin()+j]=(*bit)->OutParticle(j)->Flav();
    EXTRAXS::Single_XS *newxs = group->XSSelector()->GetXS(group->Nin(),group->Nout(),flavour);
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
    group->Add(newxs);
  }
  group->SetName(processname);
  p_processes->PushBack(group);
  std::vector<std::string> comments;
  comments.push_back(std::string("processes : ")+processname);
  GridHandlerVector gridhandler=GridHandlerVector(2);
  for (unsigned int i=0;i<gridhandler.size();++i) gridhandler[i] = new GridHandlerType();
  GridCreatorType *gridcreator = new GridCreatorType(gridhandler,p_processes);
  gridcreator->ReadInArguments(m_inputfile,m_inputpath);
  if (mkdir(m_outputpath.c_str(),448)==0) {
    ATOOLS::msg.Out()<<"Simple_Chain::CreateGrid(..): "
		     <<"Created output directory "<<m_outputpath<<"."<<std::endl;
  }
  gridcreator->SetXSExtension(m_xsextension);
  gridcreator->SetMaxExtension(m_maxextension);
  gridcreator->SetOutputPath(m_outputpath);
  gridcreator->SetOutputFile(filename);
  gridcreator->CreateGrid();
  gridcreator->WriteOutGrid(comments);
  delete gridcreator;
  for (unsigned int i=0;i<gridhandler.size();++i) delete gridhandler[i];
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
  if (!m_external) {
    p_environment = new AMEGIC::Environment(m_inputpath,m_environmentfile);
    p_environment->InitializeTheEnvironment();
    p_model=p_environment->Model();
    p_beam=p_environment->BeamSpectraHandler();
    p_isr=p_environment->ISRHandler();
  }
  SetStart(sqrt(p_isr->SprimeMax()),1);
  SetStop(sqrt(p_isr->SprimeMin()),1);
  p_processes = new EXTRAXS::SimpleXSecs(m_inputpath,m_xsfile,p_model);
  if (p_processes->Size()>0) {
    ATOOLS::msg.Tracking()<<"Simple_Chain::InitializeBlobList(): "
			  <<"Found an initialized process group."<<std::endl
			  <<"   Empty group and start with QCD_Processes_C."<<std::endl;
    p_processes->Clear();
  }
  p_processes->InitializeProcesses(p_beam,p_isr);  
  std::vector<EXTRAXS::QCD_Processes_C*> group=std::vector<EXTRAXS::QCD_Processes_C*>(m_blobs.size());
  ATOOLS::Flavour flavour[4];
  ATOOLS::Flow::ResetCounter();
  for (unsigned int i=0;i<m_blobs.size();++i) {
    flavour[0]=flavour[1]=flavour[2]=flavour[3]=ATOOLS::kf::jet;
    group[i]=new EXTRAXS::QCD_Processes_C(p_isr,p_beam,flavour,p_processes->SelectorData(),
					  p_processes->ScaleScheme(),p_processes->KFactorScheme(),
					  p_processes->ScaleFactor(),false);
    for (Blob_Iterator bit=m_blobs[i].begin();bit!=m_blobs[i].end();++bit) {
      for (int j=0;j<group[i]->Nin();++j) flavour[j]=(*bit)->InParticle(j)->Flav();
      for (int j=0;j<group[i]->Nout();++j) flavour[group[i]->Nin()+j]=(*bit)->OutParticle(j)->Flav();
      EXTRAXS::Single_XS *newxs = group[i]->XSSelector()->GetXS(group[i]->Nin(),group[i]->Nout(),flavour);
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
      group[i]->Add(newxs);
    }
    group[i]->SetName(m_processname[i]);
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
  }
  return true;
}

bool Simple_Chain::CalculateTotal()
{
  if (m_differential.size()==0) return false;
  GridFunctionType *differential;
  GridResultType total;
  differential = new GridHandlerType::GridFunctionType();
  differential->SetMonotony(differential->None);
  differential->XAxis()->SetVariable(m_differential[0]->XAxis()->Variable());
  differential->YAxis()->SetVariable(m_differential[0]->YAxis()->Variable());
  differential->XAxis()->SetScaling(m_differential[0]->XAxis()->Scaling()->Name());
  differential->YAxis()->SetScaling(m_differential[0]->YAxis()->Scaling()->Name());
  std::vector<GridFunctionType*>::iterator diffit=m_differential.begin();
  for (unsigned int i=0;i<(*diffit)->XDataSize();++i) {
    differential->AddPoint((*diffit)->XYData(i).first,(*diffit)->XYData(i).second);
  }
  for (++diffit;diffit!=m_differential.end();++diffit) {
    for (unsigned int i=0;i<(*diffit)->XDataSize();++i) {
      GridArgumentType x=(*diffit)->XYData(i).first;
      differential->ReplaceXPoint(x,differential->Y(x,differential->Interpolation)+(*diffit)->XYData(i).second);
    }
  }
  p_total = differential->IntegralY();
  total=p_total->YData(p_total->YDataSize()-1);
  p_total->SetMonotony(p_total->None);
  for (unsigned int i=0;i<p_total->XDataSize();++i) {
    p_total->ReplaceXPoint(p_total->XYData(i).first,total-p_total->XYData(i).second);
  }
  p_total->SetMonotony(p_total->MUnknown);
  p_total->MoveY(-p_total->YMin());
  p_total->ScaleY((GridArgumentType)1.0/p_total->YMax());
#ifdef DEBUG__Simple_Chain
  std::vector<std::string> comments;
  comments.push_back("   Integrated XS    "); 
  GridHandlerType *gridhandler = new GridHandlerType(p_total);
  GridCreatorBaseType *gridcreator = new GridCreatorBaseType();
  gridcreator->SetOutputPath(m_outputpath);
  gridcreator->SetOutputFile(integralfile);
  gridcreator->WriteSingleGrid(gridhandler,comments);
  delete gridcreator;
  delete gridhandler;
#endif
  delete differential;
  return true;
}

bool Simple_Chain::Initialize()
{
  if (!CheckInputPath()) return false;
  if (!CheckInputFile()) return false;
  CleanUp();
  ATOOLS::Data_Reader *reader = new ATOOLS::Data_Reader("=",";","!");
  if (!m_external) {
    std::string initfile;
    ATOOLS::ParticleInit(m_inputpath);
    reader->SetFileName(m_inputpath+m_inputfile);
    if (!reader->ReadFromFile(initfile,"ENVIRONMENT")) initfile=std::string("Run.dat");
    ATOOLS::rpa.Init(m_inputpath,initfile);
  }
  if (!ReadInData()) return false;
  reader->SetFileName(m_inputpath+m_environmentfile);
  if (!reader->ReadFromFile(m_xsfile,"XS_FILE")) m_xsfile=std::string("XS.dat");
  delete reader;
  ATOOLS::Data_Writer *writer = new ATOOLS::Data_Writer("=",";","!");
  writer->SetFileName(m_inputpath+m_xsfile);
  writer->SetBlank(32);
  writer->WriteComment("========================");
  writer->WriteComment("     Dummy XS File      ");
  writer->WriteComment("========================");
  writer->WriteToFile(std::string(" "));
  writer->WriteToFile(std::string(" Init QCD 2->2"));
  delete writer;
  for (unsigned int i=0;i<m_blobs.size();++i) {
    GridHandlerType *xsgridhandler = new GridHandlerType();
    GridHandlerType *maxgridhandler = new GridHandlerType();
    if (xsgridhandler->ReadIn(ATOOLS::Type::TFStream,m_outputpath+m_filename[i]+m_xsextension)&&
	maxgridhandler->ReadIn(ATOOLS::Type::TFStream,m_outputpath+m_filename[i]+m_maxextension)) {
      m_differential.push_back(new GridFunctionType(*xsgridhandler->Grid()));
      m_maximum.push_back(new GridFunctionType(*maxgridhandler->Grid()));
      delete xsgridhandler;
      delete maxgridhandler;
    }
    else {
      ATOOLS::msg.Error()<<"Simple_Chain::Initialize(): "
			 <<"File "<<m_filename[i]+m_xsextension<<" does not exist "<<std::endl
			 <<"   or does not contain any grid information."<<std::endl
			 <<"   Scheduling corresponding blob for grid creation."<<std::endl;
      m_create[i]=true;
    }
    if (m_create[i]) {
      if (!CreateGrid(m_blobs[i],m_filename[i],m_processname[i])) {
	ATOOLS::msg.Error()<<"Simple_Chain::Initialize(): "
			   <<"Grid creation for "<<m_filename[i]<<" failed! "<<std::endl
			   <<"   Abort initialization."<<std::endl;
	CleanUp();
	delete xsgridhandler;
	delete maxgridhandler;
	return false;
      }
      if (xsgridhandler->ReadIn(ATOOLS::Type::TFStream,m_outputpath+m_filename[i]+m_xsextension)&&
	  maxgridhandler->ReadIn(ATOOLS::Type::TFStream,m_outputpath+m_filename[i]+m_maxextension)) {
	m_differential.push_back(new GridFunctionType(*xsgridhandler->Grid()));
	m_maximum.push_back(new GridFunctionType(*maxgridhandler->Grid()));
	delete xsgridhandler;
	delete maxgridhandler;
      }
      else {
	ATOOLS::msg.Error()<<"Simple_Chain::Initialize(): "
			   <<"Grid creation failed for "<<m_outputpath+m_filename[i]<<std::endl
			   <<"   Run cannot continue."<<std::endl;
	abort();
      }
    }
  }
  if (!CalculateTotal()) {
    ATOOLS::msg.Error()<<"Simple_Chain::Initialize(): "
		       <<"Determination of \\sigma_{tot} failed. "<<std::endl
		       <<"   Run cannot continue."<<std::endl;
    abort();
  }
  SetStart(p_total->XMax(),0);
  SetStop(p_total->XMin(),0);
  if (!InitializeBlobList()) {
    ATOOLS::msg.Error()<<"Simple_Chain::Initialize(): "
		       <<"Cannot initialize selected processes. "<<std::endl
		       <<"   Run cannot continue."<<std::endl;
    abort();
  }  
  return true;
}

bool Simple_Chain::FillBlob(ATOOLS::Blob *blob)
{
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
  if (m_selected<(unsigned int)p_processes->Size()) {
#ifdef DEBUG__Simple_Chain
    std::cout<<"Simple_Chain::FillBlob(..): Generating one event."<<std::endl;
#endif
    if ((*p_processes)[m_selected]->OneEvent()) {
#ifdef DEBUG__Simple_Chain
      std::cout<<"   Completed one event."<<std::endl;
#endif
      p_xs=(*p_processes)[m_selected]->Selected();
      ATOOLS::Vec4D ptot;
      for (int j=0;j<p_xs->Nin();++j) ptot+=p_xs->Momenta()[j];
      m_last[1]-=sqrt(ptot.Abs2());
      p_xs->SetColours(p_xs->Momenta());
      ATOOLS::Particle *particle;
      for (int j=0;j<p_xs->Nin();++j) {
	particle = new ATOOLS::Particle(0,p_xs->Flavs()[j]);
	particle->SetMomentum(p_xs->Momenta()[j]);
 	particle->SetFlow(1,p_xs->Colours()[j][0]);
 	particle->SetFlow(2,p_xs->Colours()[j][1]);
	particle->SetStatus(1);
	blob->AddToInParticles(particle);
      }
      for (int j=0;j<p_xs->Nout();++j) {
	particle = new ATOOLS::Particle(0,p_xs->Flavs()[p_xs->Nin()+j]);
	particle->SetMomentum(p_xs->Momenta()[p_xs->Nin()+j]);
 	particle->SetFlow(1,p_xs->Colours()[p_xs->Nin()+j][0]);
 	particle->SetFlow(2,p_xs->Colours()[p_xs->Nin()+j][1]);
	particle->SetStatus(1);
	blob->AddToOutParticles(particle);
      }
      return true;
    }
    ATOOLS::msg.Error()<<"Simple_Chain::FillBlob(): "
		       <<"Could not select any process! Abort."<<std::endl;
    return false;
  }
  ATOOLS::msg.Error()<<"Simple_Chain::FillBlob(): "
		     <<"Internal error!"
		     <<"   Please submit a bug report."<<std::endl;
  return false;
}

bool Simple_Chain::DiceProcess()
{
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
      if (m_last[1]<(*p_processes)[m_selected]->ISR()->SprimeMin()) {
	ATOOLS::msg.Tracking()<<"Simple_Chain::DiceProcess(): s' out of bounds."<<std::endl
			      <<"   Cannot create any process. Abort."<<std::endl;
	return false;
      }
      double sprimemax=(*p_processes)[m_selected]->ISR()->SprimeMax();
      (*p_processes)[m_selected]->SetMax((*m_maximum[m_selected])(m_last[0]*0.2),1);
      (*p_processes)[m_selected]->ISR()->SetSprimeMax(m_last[1]*m_last[1]);
      m_dicedprocess=true;
      bool result=FillBlob(p_blob);
      (*p_processes)[m_selected]->ISR()->SetSprimeMax(sprimemax);
      return result;
    }
  }
  ATOOLS::msg.Tracking()<<"Simple_Chain::DiceProcess(): "
			<<"Could not select any process. "<<std::endl
			<<"   Returning most likely instead."<<std::endl;
  if (sorter.XDataSize()!=0) {
    m_selected=sorter.XYData(0).second;
    (*p_processes)[m_selected]->SetMax((*m_maximum[m_selected])(m_last[0]*0.1),1);
    m_dicedprocess=true;
    return FillBlob(p_blob);
  }
  m_dicedprocess=false;
  return false;
}

bool Simple_Chain::DiceOrderingParameter()
{ 
  if (m_last[0]<=m_stop[0]) {
    ATOOLS::msg.Error()<<"Simple_Chain::DiceOrderingParameter(): "
		       <<"Ordering parameter exceeded allowed range."<<std::endl
		       <<"   Cannot proceed in blob creation."<<std::endl;
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

void Simple_Chain::Reset()
{
  m_last[0]=m_start[0];
  m_last[1]=m_start[1];
}

void Simple_Chain::Update(MI_Base *mibase)
{
  UpdateAll(this);
  return;
}
