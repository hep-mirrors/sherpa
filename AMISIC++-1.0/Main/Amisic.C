#include "Amisic.H"

#include "Particle.H"
#include "Random.H"
#include "Channel_Elements.H"

#ifdef DEBUG__Amisic
const std::string integralfile=std::string("integral.dat");
#endif

using namespace AMISIC;

template <class Argument_Type,class Result_Type>
Amisic::Grid_Creator<Argument_Type,Result_Type>::Grid_Creator(GridHandlerType *_p_gridhandler,
							      EXTRAXS::SimpleXSecs *_p_processes):
  GridCreatorBaseType(_p_gridhandler),
  p_processes(_p_processes)
{
  if (p_processes==NULL) {
      ATOOLS::msg.Error()<<"Grid_Creator::Grid_Creator("<<_p_gridhandler<<","<<_p_processes<<"): "
			 <<"Process handler is not initialized! Abort."<<std::endl;
      abort();
  }
}

template <class Argument_Type,class Result_Type>
bool Amisic::Grid_Creator<Argument_Type,Result_Type>::InitializeCalculation()
{
  p_xaxis=p_gridhandler->Grid()->XAxis();
  m_criterion=ATOOLS::Variable::TypeToSelectorID(p_xaxis->Variable().Type());
  m_initialdata=p_processes->SelectorData()->RemoveData(m_criterion);
  return true;
}

template <class Argument_Type,class Result_Type>
Result_Type Amisic::Grid_Creator<Argument_Type,Result_Type>::CalculateSingleValue(GridArgumentType nextleft,
										  GridArgumentType nextright)
{
  GridArgumentType lower, upper, middle;
  GridResultType result;
  lower=ATOOLS::Max((*p_xaxis)[(*p_xaxis)(nextleft)*(GridArgumentType)(3.0/4.0)
			       +(*p_xaxis)(nextright)*(GridArgumentType)(1.0/4.0)],(*p_xaxis)[GridXMin()]);
  upper=ATOOLS::Min((*p_xaxis)[(*p_xaxis)(nextleft)*(GridArgumentType)(1.0/4.0)
			       +(*p_xaxis)(nextright)*(GridArgumentType)(3.0/4.0)],(*p_xaxis)[GridXMax()]);
  middle=(*p_xaxis)[((*p_xaxis)(nextleft)+(*p_xaxis)(nextright))/(GridArgumentType)2.0];
  p_processes->SelectorData()->SetData(m_criterion,m_initialdata.flavs,m_initialdata.help,lower,upper);
  p_processes->ResetSelector(p_processes->SelectorData());
  p_processes->CalculateTotalXSec();
  result=(GridResultType)p_processes->Total()/(upper-lower);
  ATOOLS::msg.Out()<<"Amisic::Grid_Creator::CalculateSingleValue(): Got value for "<<middle<<" GeV"<<std::endl
		   <<"   Calculation for "<<lower<<" GeV < "<<p_xaxis->Variable().Name()
		   <<" < "<<upper<<" GeV yielded "<<result*rpa.Picobarn()<<" pb/GeV"<<std::endl;
  return result;
}

Amisic::Amisic():
  m_differential(std::vector<GridFunctionType*>(0)),
  p_total(NULL),
  m_start((GridArgumentType)0.0),
  m_stop((GridArgumentType)1.6),
  m_inputdirectory(std::string("./")),
  m_inputfile(std::string("MI.dat")),
  m_environmentfile(std::string("Run.dat")),
  m_xsfile(std::string("XS.dat")),
  m_outputdirectory(std::string("./Grid")),
  p_environment(NULL),
  p_processes(NULL),
  p_fsrinterface(NULL) {}

Amisic::~Amisic()
{
  CleanUp();
}

void Amisic::CleanUp() 
{
  if (p_fsrinterface!=NULL) delete p_fsrinterface;
  if (p_processes!=NULL) delete p_processes;
  if (p_environment!=NULL) delete p_environment;
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

bool Amisic::HaveBlob(ATOOLS::Blob *blob) 
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

ATOOLS::Blob *Amisic::GetBlob(ATOOLS::Flavour flavour[4])
{
  ATOOLS::Blob *newblob = new ATOOLS::Blob();
  for (unsigned int j=0;j<2;++j) newblob->AddToInParticles(new ATOOLS::Particle(j,flavour[j]));
  for (unsigned int k=2;k<4;++k) newblob->AddToOutParticles(new ATOOLS::Particle(k,flavour[k]));
  if (HaveBlob(newblob)) {
    ATOOLS::msg.Tracking()<<"Amisic::CreateBlob(..): "
			  <<"Tried to add a blob twice. Ignore last attempt."<<std::endl; 
    delete newblob;
    return NULL;
  }
  return newblob;
}

void Amisic::FillMode(EXTRAXS::QCD_Processes::Mode mode)
{
  const unsigned int nflavour=4;
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
    m_filename.push_back("gg_to_gg__grid.dat");
    m_processname.push_back("g g -> g g");
    m_create.push_back(false);
    if (mode==EXTRAXS::QCD_Processes::gggg) break;
  case EXTRAXS::QCD_Processes::qqbgg:
    m_blobs.push_back(ATOOLS::Blob_List(0));
    temp[3]=temp[2]=ATOOLS::Flavour((ATOOLS::kf::code)21);
    for (i=1;i<=nflavour;++i) {
      temp[0]=ATOOLS::Flavour((ATOOLS::kf::code)i);
      temp[1]=ATOOLS::Flavour((ATOOLS::kf::code)i).Bar();
      if ((newblob=GetBlob(temp))!=NULL) {
	m_blobs[m_blobs.size()-1].push_back(newblob);
      }
    }
    m_filename.push_back("qqb_to_gg__grid.dat");
    m_processname.push_back("q qb -> g g");
    m_create.push_back(false);
    if (mode==EXTRAXS::QCD_Processes::qqbgg) break;
  case EXTRAXS::QCD_Processes::ggqqb:
    m_blobs.push_back(ATOOLS::Blob_List(0));
    temp[1]=temp[0]=ATOOLS::Flavour((ATOOLS::kf::code)21);
    for (i=1;i<=nflavour;++i) {
      temp[2]=ATOOLS::Flavour((ATOOLS::kf::code)i);
      temp[3]=ATOOLS::Flavour((ATOOLS::kf::code)i).Bar();
      if ((newblob=GetBlob(temp))!=NULL) {
	m_blobs[m_blobs.size()-1].push_back(newblob);
      }
    }
    m_filename.push_back("gg_to_qqb__grid.dat");
    m_processname.push_back("g g -> q qb");
    m_create.push_back(false);
    if (mode==EXTRAXS::QCD_Processes::ggqqb) break;
  case EXTRAXS::QCD_Processes::qgqg:
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
    m_filename.push_back("qg_to_qg__grid.dat");
    m_processname.push_back("q g -> q g");
    m_create.push_back(false);
    if (mode==EXTRAXS::QCD_Processes::qgqg) break;
  case EXTRAXS::QCD_Processes::q1q2q1q2:
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
  case EXTRAXS::QCD_Processes::q1q1q1q1:
    if (mode==EXTRAXS::QCD_Processes::q1q1q1q1) m_blobs.push_back(ATOOLS::Blob_List(0));
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
    if (mode==EXTRAXS::QCD_Processes::q1q1q1q1) {
      m_filename.push_back("q1q1_to_q1q1__grid.dat");
      m_processname.push_back("q1 q1 -> q1 q1");
      m_create.push_back(false);
      break;
    }
    m_filename.push_back("q1q2_to_q1q2__grid.dat");
    m_processname.push_back("q1 q2 -> q1 q2");
    m_create.push_back(false);
    if (mode==EXTRAXS::QCD_Processes::q1q2q1q2) break;
  case EXTRAXS::QCD_Processes::q1q2bq1q2b:
    m_blobs.push_back(ATOOLS::Blob_List(0));
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
    m_filename.push_back("q1q2b_to_q1q2b__grid.dat");
    m_processname.push_back("q1 q2b -> q1 q2b");
    m_create.push_back(false);
    if (mode==EXTRAXS::QCD_Processes::q1q2bq1q2b) break;
  case EXTRAXS::QCD_Processes::q1q1bq1q1b:
    m_blobs.push_back(ATOOLS::Blob_List(0));
    for (i=1;i<=nflavour;++i) {
      temp[2]=temp[0]=ATOOLS::Flavour((ATOOLS::kf::code)i);
      temp[3]=temp[1]=ATOOLS::Flavour((ATOOLS::kf::code)i).Bar();
      if ((newblob=GetBlob(temp))!=NULL) {
	m_blobs[m_blobs.size()-1].push_back(newblob);
      }
    }
    m_filename.push_back("q1q1b_to_q1q1b__grid.dat");
    m_processname.push_back("q1 q1b -> q1 q1b");
    m_create.push_back(false);
    if (mode==EXTRAXS::QCD_Processes::q1q1bq1q1b) break;
  case EXTRAXS::QCD_Processes::q1q1bq2q2b:
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
    m_filename.push_back("q1q1b_to_q2q2b__grid.dat");
    m_processname.push_back("q1 q1b -> q2 q2b");
    m_create.push_back(false);
    if (mode==EXTRAXS::QCD_Processes::q1q1bq2q2b) break;
    break;
  case EXTRAXS::QCD_Processes::Unknown:
  case EXTRAXS::QCD_Processes::None:
    ATOOLS::msg.Error()<<"Amisic::FillMode("<<mode<<"): "
		       <<"Wrong parameter. Abort."<<std::endl;
    break;
  }
}

bool Amisic::ReadInData()
{
  std::vector<std::string> comments;
  comments.push_back("->");
  comments.push_back("FOR");
  comments.push_back("INTO");
  comments.push_back("IN");
  ATOOLS::Data_Reader *reader = new ATOOLS::Data_Reader("=",";","!");
  reader->SetFileName(m_inputdirectory+m_inputfile);
  reader->AddIgnore(comments);
  if (!reader->ReadFromFile(m_outputdirectory,"GRID DIRECTORY")) m_outputdirectory=std::string("./");
  std::vector<std::vector<std::string> > temp;
  reader->ArrayFromFile(temp,"CREATE GRID",ATOOLS::noinputtag,reader->MTransposed);
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
	  filename+=std::string("_grid.dat");
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

bool Amisic::CreateGrid(ATOOLS::Blob_List& bloblist,std::string& filename,std::string& processname)
{
  p_environment = new AMEGIC::Environment(m_inputdirectory,m_environmentfile);
  p_environment->InitializeTheEnvironment();
  p_processes = new EXTRAXS::SimpleXSecs(m_inputdirectory,m_xsfile,p_environment->Model());
  if (p_processes->Size()>0) {
    ATOOLS::msg.Tracking()<<"Amisic::CreateGrid(..): "
			  <<"Found an initialized process group."<<std::endl
			  <<"   Empty group and start with QCD_Processes."<<std::endl;
    p_processes->Clear();
  }
  p_processes->InitializeProcesses(p_environment->BeamSpectraHandler(),p_environment->ISRHandler());  
  ATOOLS::Flavour flavour[4];
  flavour[0]=flavour[1]=flavour[2]=flavour[3]=ATOOLS::kf::jet;
  EXTRAXS::QCD_Processes *group;
  group = new EXTRAXS::QCD_Processes(p_environment->ISRHandler(),p_environment->BeamSpectraHandler(),
				     flavour,p_processes->SelectorData(),p_processes->ScaleScheme(),
				     p_processes->KFactorScheme(),p_processes->ScaleFactor(),false);
  for (Blob_Iterator bit=bloblist.begin();bit!=bloblist.end();++bit) {
    for (int i=0;i<group->Nin();++i) flavour[i]=(*bit)->InParticle(i)->Flav();
    for (int j=0;j<group->Nout();++j) flavour[group->Nin()+j]=(*bit)->OutParticle(j)->Flav();
    EXTRAXS::Single_XS *newxs = group->XSSelector()->GetXS(group->Nin(),group->Nout(),flavour);
    if (newxs==NULL) {
      ATOOLS::msg.Error()<<"Amisic::CreateGrid(..): "
			 <<"Did not find any process! Abort calculation."<<std::endl;
      delete p_processes;
      p_processes=NULL;
      delete p_environment;
      p_environment=NULL;
      return false;
    }
    group->Add(newxs);
  }
  group->SetName(processname);
  p_processes->PushBack(group);
  std::vector<std::string> comments;
  comments.push_back(std::string("processes : ")+processname);
  GridHandlerType *gridhandler = new GridHandlerType();
  GridCreatorType *gridcreator = new GridCreatorType(gridhandler,p_processes);
  gridcreator->ReadInArguments(m_inputdirectory+m_inputfile);
  if (mkdir(m_outputdirectory.c_str(),448)==0) {
    ATOOLS::msg.Out()<<"Amisic::CreateGrid(..): "
		     <<"Created output directory "<<m_outputdirectory<<"."<<std::endl;
  }
  gridcreator->CreateGrid(m_outputdirectory+filename);
  gridcreator->WriteOutGrid(m_outputdirectory+filename,comments);
  delete gridcreator;
  delete gridhandler;
  delete p_processes;
  p_processes=NULL;
  delete p_environment;
  p_environment=NULL;
  return true;
}

bool Amisic::InitializeBlobList()
{  
  p_environment = new AMEGIC::Environment(m_inputdirectory,m_environmentfile);
  p_environment->InitializeTheEnvironment();
  p_processes = new EXTRAXS::SimpleXSecs(m_inputdirectory,m_xsfile,p_environment->Model());
  if (p_processes->Size()>0) {
    ATOOLS::msg.Tracking()<<"Amisic::InitializeBlobList(): "
			  <<"Found an initialized process group."<<std::endl
			  <<"   Empty group and start with QCD_Processes."<<std::endl;
    p_processes->Clear();
  }
  p_processes->InitializeProcesses(p_environment->BeamSpectraHandler(),p_environment->ISRHandler());  
  std::vector<EXTRAXS::QCD_Processes*> group=std::vector<EXTRAXS::QCD_Processes*>(m_blobs.size());
  ATOOLS::Flavour flavour[4];
  ATOOLS::Flow::ResetCounter();
  for (unsigned int i=0;i<m_blobs.size();++i) {
    flavour[0]=flavour[1]=flavour[2]=flavour[3]=ATOOLS::kf::jet;
    group[i]=new EXTRAXS::QCD_Processes(p_environment->ISRHandler(),p_environment->BeamSpectraHandler(),
					flavour,p_processes->SelectorData(),p_processes->ScaleScheme(),
					p_processes->KFactorScheme(),p_processes->ScaleFactor(),false);
    for (Blob_Iterator bit=m_blobs[i].begin();bit!=m_blobs[i].end();++bit) {
      for (int j=0;j<group[i]->Nin();++j) flavour[j]=(*bit)->InParticle(j)->Flav();
      for (int j=0;j<group[i]->Nout();++j) flavour[group[i]->Nin()+j]=(*bit)->OutParticle(j)->Flav();
      EXTRAXS::Single_XS *newxs = group[i]->XSSelector()->GetXS(group[i]->Nin(),group[i]->Nout(),flavour);
      if (newxs==NULL) {
	ATOOLS::msg.Error()<<"Amisic::InitializeBlobList(): "
			   <<"Did not find any process! Abort."<<std::endl;
	delete group[i];
 	delete p_processes;
 	p_processes=NULL;
	delete p_environment;
	p_environment=NULL;
	return false;
      }
      group[i]->Add(newxs);
    }
    group[i]->SetName(m_processname[i]);
    p_processes->PushBack(group[i]);
  }
  p_processes->CalculateTotalXSec();
  p_fsrinterface = new FSRChannel(2,2,flavour,p_total->XAxis()->Variable());
  p_fsrinterface->SetAlpha(1.0);
  p_fsrinterface->SetAlphaSave(1.0);
  for (unsigned int i=0;i<m_blobs.size();++i) {
    group[i]->SetFSRInterface(p_fsrinterface);
    group[i]->SetFSRMode(2);
    group[i]->CreateFSRChannels();
  }
  return true;
}

bool Amisic::CalculateTotal()
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
  p_total->Monotony();
  p_total->MoveY(-p_total->YMin());
  p_total->ScaleY((GridArgumentType)1.0/p_total->YMax());
#ifdef DEBUG__Amisic
  EXTRAXS::SimpleXSecs *processes;
  GridHandlerType *gridhandler = new GridHandlerType(p_total);
  GridCreatorType *gridcreator = new GridCreatorType(gridhandler,processes);
  gridcreator->WriteOutGrid(integralfile);
  delete gridcreator;
  delete gridhandler;
#endif
  delete differential;
  return true;
}

bool Amisic::Initialize(std::string tempidir,std::string tempifile,bool creategrid)
{
  CleanUp();
  if (tempidir!=ATOOLS::nullstring) SetInputDirectory(tempidir);
  if (tempifile!=ATOOLS::nullstring) SetInputFile(tempifile);
  CleanUp();
  ATOOLS::ParticleInit(m_inputdirectory);
  ATOOLS::rpa.Init(m_inputdirectory);
  if (!ReadInData()) return false;
  ATOOLS::Data_Reader *reader = new ATOOLS::Data_Reader("=",";","!");
  reader->SetFileName(m_inputdirectory+m_environmentfile);
  if (!reader->ReadFromFile(m_xsfile,"XS_FILE")) m_xsfile=std::string("XS.dat");
  delete reader;
  ATOOLS::Data_Writer *writer = new ATOOLS::Data_Writer("=",";","!");
  writer->SetFileName(m_inputdirectory+m_xsfile);
  writer->SetBlank(32);
  writer->WriteComment("========================");
  writer->WriteComment("     Dummy XS File      ");
  writer->WriteComment("========================");
  writer->WriteToFile(std::string(" "));
  writer->WriteToFile(std::string(" Init QCD 2->2"));
  delete writer;
  for (unsigned int i=0;i<m_blobs.size();++i) {
    GridHandlerType *newgridhandler = new GridHandlerType();
    if (!newgridhandler->ReadIn(ATOOLS::Type::TFStream,m_outputdirectory+m_filename[i])) {
      ATOOLS::msg.Error()<<"Amisic::Initialize("<<tempidir<<","<<tempifile<<"): "
			 <<"File "<<m_outputdirectory+m_filename[i]<<" does not exist "<<std::endl
			 <<"   or does not contain any grid information."<<std::endl
			 <<"   Scheduling corresponding blob for grid creation."<<std::endl;
      m_create[i]=true;
    }
    else {
      m_differential.push_back(new GridFunctionType(*newgridhandler->Grid()));
      delete newgridhandler;
    }
    if (m_create[i]) {
      if (creategrid) {
	if (!CreateGrid(m_blobs[i],m_filename[i],m_processname[i])) {
	  ATOOLS::msg.Error()<<"Amisic::Initialize("<<tempidir<<","<<tempifile<<"): Grid creation for "
			     <<m_outputdirectory+m_filename[i]<<" failed! "<<std::endl
			     <<"   Abort initialization."<<std::endl;
	  CleanUp();
	  delete newgridhandler;
	  return false;
	}
	if (!newgridhandler->ReadIn(ATOOLS::Type::TFStream,m_outputdirectory+m_filename[i])) {
	  ATOOLS::msg.Error()<<"Amisic::Initialize("<<tempidir<<","<<tempifile<<"): Grid creation for "
			     <<m_outputdirectory+m_filename[i]<<" failed! "<<std::endl
			     <<"   Abort initialization."<<std::endl;
	  CleanUp();
	  delete newgridhandler;
	  return false;
	}
	else {
	  m_differential.push_back(new GridFunctionType(*newgridhandler->Grid()));
	  delete newgridhandler;
	}
      }
      else {
	ATOOLS::msg.Error()<<"Amisic::Initialize("<<tempidir<<","<<tempifile<<"): "
			   <<"Grid files are missing!"<<std::endl
			   <<"   Run cannot continue."<<std::endl;
	delete newgridhandler;
	abort();
      }
    }
  }
  if (!CalculateTotal()) {
    ATOOLS::msg.Error()<<"Amisic::Initialize("<<tempidir<<","<<tempifile<<"): "
		       <<"Determination of \\sigma_{tot} failed. "<<std::endl
		       <<"   Run cannot continue."<<std::endl;
    abort();
  }
  if (!InitializeBlobList()) {
    ATOOLS::msg.Error()<<"Amisic::Initialize("<<tempidir<<","<<tempifile<<"): "
		       <<"Cannot initialize selected processes. "<<std::endl
		       <<"   Run cannot continue."<<std::endl;
    abort();
  }  
  return true;
}

ATOOLS::Blob *Amisic::CreateBlob(unsigned int i,double ran)
{
  ATOOLS::Blob *blob=NULL;
  if (p_processes==NULL) {
    ATOOLS::msg.Error()<<"Amisic::CreateBlob("<<i<<","<<ran<<"): "
		       <<"Processes not initialized! Abort."<<std::endl;
    return blob;
  }
  if (i<(unsigned int)p_processes->Size()) {
    if ((*p_processes)[i]->OneEvent()) {
      EXTRAXS::XS_Base *selected=(*p_processes)[i]->Selected();
      blob = new ATOOLS::Blob();
      for (int j=0;j<selected->Nin();++j) {
	blob->AddToInParticles(new ATOOLS::Particle(selected->Flavs()[j]));
	blob->InParticle(j)->SetMomentum(selected->Momenta()[j]);
      }
      for (int j=0;j<selected->Nout();++j) {
	blob->AddToOutParticles(new ATOOLS::Particle(selected->Flavs()[selected->Nin()+j]));
	blob->OutParticle(j)->SetMomentum(selected->Momenta()[selected->Nin()+j]);
      }
      return blob;
    }
  }
  if (blob==NULL) ATOOLS::msg.Error()<<"Amisic::SelectBlob("<<i<<","<<ran<<"): "
				     <<"Could not select any process! Abort."<<std::endl;
  return blob;
}

ATOOLS::Blob *Amisic::DiceProcess(GridArgumentType parameter,double ran[2])
{
  if (m_differential.size()==0) return NULL;
  p_fsrinterface->SetValue(parameter);
  int criterion=ATOOLS::Variable::TypeToSelectorID(m_differential[0]->XAxis()->Variable().Type());
  ATOOLS::Mom_Data m_initialdata=p_processes->SelectorData()->RemoveData(criterion);
  p_processes->SelectorData()->SetData(criterion,m_initialdata.flavs,m_initialdata.help,
				       parameter,m_initialdata.max);
  p_processes->ResetSelector(p_processes->SelectorData());
  Data_To_Function<GridResultType,unsigned int> sorter;
  GridResultType cur, norm=(GridResultType)0.0;
  for (unsigned int i=0;i<m_differential.size();++i) {
    cur=(*m_differential[i])(parameter);
    sorter.AddPoint(cur,i);
    norm+=cur;
  }
  cur=(GridResultType)0.0;
  for (int i=sorter.XDataSize()-1;i>=0;--i) {
    if ((cur+=sorter.XData(i)/norm)>ran[0]) {
      return CreateBlob(sorter.XYData(i).second,ran[1]);
    }
  }
  ATOOLS::msg.Tracking()<<"Amisic::DiceProcess("<<parameter<<","<<ran[0]<<","<<ran[1]<<"): "
			<<"Could not select any process. "<<std::endl
			<<"   Returning most likely instead."<<std::endl;
  if (sorter.XDataSize()!=0) return CreateBlob(sorter.XYData(0).second,ran[1]);
  return NULL;
}

Amisic::GridArgumentType Amisic::DiceParameter(GridArgumentType &last,double ran)
{ 
  return last=(*p_total)[(*p_total)(last)-log(ran)]; 
}

ATOOLS::Blob_List *Amisic::CreateProcesses()
{
  ATOOLS::Blob_List *blobs = new ATOOLS::Blob_List();
  double rans[3];
  GridArgumentType last=m_start;
  if (last<=m_stop) last=p_total->XMax();
  ATOOLS::Blob *blob;
  do {
    for (unsigned int i=0;i<3;++i) do { rans[i]=ATOOLS::ran.Get(); } while (rans[i]==0.0);  
    if ((blob=DiceProcess(DiceParameter(last,rans[2]),rans))!=NULL) blobs->push_back(blob);
  } while (last>m_stop);
  return blobs;
}
