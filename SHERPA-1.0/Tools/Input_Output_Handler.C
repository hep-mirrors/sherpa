#include "Input_Output_Handler.H"
#include "Message.H"
#include "Run_Parameter.H"
#include "Data_Read.H"
#include "Exception.H"
#include "MyStrStream.H"
#include <stdio.h>

using namespace SHERPA;
using namespace ATOOLS;

extern "C" {
  void outhepevt_();
}

const iotype::code SHERPA::operator|(const iotype::code code1,const iotype::code code2)
{
  return (iotype::code)((int)code1|(int)code2);
}
  
const iotype::code SHERPA::operator&(const iotype::code code1,const iotype::code code2)
{
  return (iotype::code)((int)code1&(int)code2);
}  


Input_Output_Handler::Input_Output_Handler(const std::string mode,
                                           const std::vector<std::string> & outfiles,
                                           const std::vector<std::string> & infiles,
                                           const std::string _path, const int _filesize,
                                           const int precision):
  m_on(true), m_io(1), m_precision(precision),
  m_outtype(iotype::Unknown), m_screenout(iotype::Unknown), m_intype(iotype::Unknown), 
#ifdef USING__CLHEP
  p_hepmc(NULL), 
#endif
  p_hepevt(NULL),
  p_instream(NULL),
  m_path(_path), 
  m_filesize(_filesize), 
  m_evtnumber(0), 
  m_evtcount(0)
{
  if (InitialiseOutput(mode,outfiles)) m_io += 2;
  if (InitialiseInput(infiles))   m_io += 4;
}

Input_Output_Handler::~Input_Output_Handler() 
{
#ifdef USING__CLHEP
  if (p_hepmc)        { delete p_hepmc;  p_hepmc  = NULL; }
#endif
  if (p_hepevt!=NULL) { delete p_hepevt; p_hepevt=NULL; }
  if (!m_outmap.empty()) {
    for (std::map<iotype::code,NameStream *>::iterator oit=m_outmap.begin();
         oit!=m_outmap.end();oit++) {
      if (oit->second!=NULL) { 
        if (oit->second->outstream.good()) oit->second->outstream.close();
        delete oit->second;
        oit->second=NULL;
      }
    }
    m_outmap.clear();
  }
  if (p_instream) {
    p_instream->close();
    delete p_instream; p_instream=NULL;
  }
}

bool Input_Output_Handler::InitialiseInput(const std::vector<std::string> & infiles) {
  return false;

//   for (size_t i=0;i<infiles.size();++i) {
//     m_intype=m_intype|(iotype::code)(pow(2.,int(i))*(infiles[i]!=std::string("")));
//   }
//   std::string gentype;
//   switch (m_intype) {
//   case iotype::Sherpa:
//     m_filename = m_path+std::string("/")+infiles[0]+std::string(".evts"); 
//     p_instream = new std::ifstream(m_filename.c_str()); 
//     if (!p_instream->good()) {
//       msg.Error()<<"ERROR in Input_Output_Handler."<<std::endl
//                  <<"   Event file "<<m_filename<<" not found."<<std::endl
//                  <<"   Will abort the run."<<std::endl;
//       abort();
//     }
//     (*p_instream)>>gentype>>m_filesize;
//     if (gentype!=std::string("Sherpa")) {
//       msg.Error()<<"ERROR in HepEvt_Interface."<<std::endl
//                  <<"   Generator type:"<<gentype<<" not SHERPA."<<std::endl
//                  <<"   Cannot guarantee i/o operations. Will abort the run."<<std::endl;
//       abort();
//     }
//     break;
// #ifdef USING__CLHEP
//   case iotype::HepMC: 
//     if (p_hepmc==NULL) p_hepmc  = new HepMC_Interface();
//     break;
// #endif
//   case iotype::HepEvt:
//     // Obacht !!!!!
//     if (p_hepevt==NULL) p_hepevt = new HepEvt_Interface(false,1,m_path,infiles[2]);
//     else abort();  // Schlamassel !!!!
//     break;
//   case iotype::D0HepEvt:
//     msg.Error()<<"D0_HEPEVT_INPUT not implemented yet. Aborting."<<std::endl;
//     abort();
//     break;
//   default :
//     break;
//   }
}


bool Input_Output_Handler::InitialiseOutput(const std::string mode, 
                                            const std::vector<std::string> & outfiles) {
  iotype::code test;
  NameStream * ns;
  for (size_t i=0;i<outfiles.size();++i) {
    if (outfiles[i]!="") {
      test = (iotype::code)(pow(2.,int(i)));
      switch(test) {
      case 1:
        m_outtype   = m_outtype|test;
        ns          = new NameStream(m_path+"/"+outfiles[0],".evts",m_precision);
        m_outmap[iotype::Sherpa] = ns;
        break;
      case 2:
#ifdef USING__CLHEP
        m_outtype   = m_outtype|test;
        ns          = new NameStream(m_path+"/"+outfiles[1],".hepmc",m_precision);
        m_outmap[iotype::HepMC] = ns;
        if (p_hepmc==NULL) p_hepmc = new HepMC_Interface();
        break;
#else
        msg.Error()<<"Error in "<<METHOD<<": HepMC format can only be created when Sherpa "
                   <<"was linked with CLHEP, please read our Howto to fix this."<<std::endl;
        abort();
#endif
      case 4:
#ifdef USING__CLHEP
        m_outtype   = m_outtype|test;
        ns          = new NameStream(m_path+"/"+outfiles[2],".old.hepmc",m_precision);
        m_outmap[iotype::OldHepMC] = ns;
        if (p_hepmc==NULL) p_hepmc = new HepMC_Interface();
        break;
#else
        msg.Error()<<"Error in "<<METHOD<<": HepMC format can only be created when Sherpa "
                   <<"was linked with CLHEP, please read our Howto to fix this."<<std::endl;
        abort();
#endif
      case 8:
        m_outtype   = m_outtype|test;
        ns          = new NameStream(m_path+"/"+outfiles[3],".hepevt",m_precision);
        m_outmap[iotype::HepEvt]     = ns;
        if (p_hepevt==NULL) p_hepevt = new HepEvt_Interface();
        break;
      case 16:
        m_outtype   = m_outtype|test;
        ns          = new NameStream(m_path+"/"+outfiles[4],".d0.hepevt",m_precision);
        m_outmap[iotype::D0HepEvt]   = ns;
        if (p_hepevt==NULL) p_hepevt = new HepEvt_Interface();
        break;
      default :
        msg.LogFile()<<"ERROR in Input_Output_Handler::Input_Output_Handler("<<test<<")"<<std::endl
                     <<"   No valid output format specified. Continue run."<<std::endl;
        break;
      }
    }
  }

  if (mode=="Sherpa") m_screenout=iotype::Sherpa;
  else if (mode=="HepMC") {
#ifdef USING__CLHEP
    m_screenout=iotype::HepMC;
    if (p_hepmc==NULL) p_hepmc   = new HepMC_Interface(); 
#else
    msg.Error()<<"Error in "<<METHOD<<": HepMC format can only be created when Sherpa "
                <<"was linked with CLHEP, please read our Howto to fix this."<<std::endl;
    abort();
#endif
  }
  else if (mode=="HepEvt") {
    m_screenout=iotype::HepEvt;
    if (p_hepevt==NULL) p_hepevt = new HepEvt_Interface(  );
  }
  else m_io-=1;

  if (m_outtype==iotype::Unknown) return false;
  return true;
}

void Input_Output_Handler::AddOutputMode(const iotype::code c1) 
{
  if (m_on==false) return;
  msg.Error()<<"Error in Input_Output_Handler::AddOutputMode("<<c1<<")"<<std::endl
             <<"   Method not yet implemented. Continue run."<<std::endl;
}

void Input_Output_Handler::AddInputMode(const iotype::code c1) 
{
  if (m_on==false) return;
  msg.Error()<<"Error in Input_Output_Handler::AddInputMode("<<c1<<")"<<std::endl
             <<"   Method not yet implemented. Continue run."<<std::endl;
}

void Input_Output_Handler::PrintEvent(ATOOLS::Blob_List *const blobs) {
  if (m_on==false || !(m_io&1) || !msg.LevelIsEvents()) return;
  switch (m_screenout) {
  case iotype::Sherpa:
    if (!blobs->empty()) {
      msg.Out()<<"  -------------------------------------------------  "<<std::endl;
      for (Blob_List::iterator blit=blobs->begin();blit!=blobs->end();++blit) 
        msg.Out()<<*(*blit)<<std::endl;
      msg.Out()<<"  -------------------------------------------------  "<<std::endl;
    }
    else msg.Out()<<"  ******** Empty event ********  "<<std::endl;
    break;
  case iotype::HepMC:
#ifdef USING__CLHEP
    p_hepmc->Sherpa2HepMC(blobs);
    p_hepmc->PrintEvent(1,msg.Out());
    break;
#else
    msg.Error()<<"Error in "<<METHOD<<": HepMC format can only be created when Sherpa "
                <<"was linked with CLHEP, please read our Howto to fix this."<<std::endl;
    abort();
#endif
  case iotype::HepEvt: 
    p_hepevt->Sherpa2HepEvt(blobs);
    p_hepevt->PrintEvent(3,msg.Out(),p_hepevt->Nhep());
    break;
  default:
    msg.Error()<<"Error in "<<METHOD<<std::endl
               <<"   Unknown Output format : "<<m_outtype<<std::endl
               <<"   No output, continue run ..."<<std::endl;
    break;
  }
}

void Input_Output_Handler::ResetInterfaces() {
#ifdef USING__CLHEP
  if (p_hepmc) p_hepmc->Reset();
#endif 
  if (p_hepevt) p_hepevt->Reset();
}


bool Input_Output_Handler::OutputToFormat(ATOOLS::Blob_List *const blobs,const double weight) 
{
  if (m_on==false)          return false;
  if (!(m_io&1)&&!(m_io&2)) return false;
  
  ResetInterfaces();

  if (!m_outmap.empty()) {
    for (std::map<iotype::code,NameStream *>::iterator oit=m_outmap.begin();
         oit!=m_outmap.end();oit++) {
      if (oit->second!=NULL) {
        switch (oit->first) {
        case iotype::Sherpa:
          SherpaOutput(oit->second->outstream,blobs);
          break;
        case iotype::HepMC:
#ifdef USING__CLHEP
          p_hepmc->Sherpa2HepMC(blobs);
          p_hepmc->PrintEvent(1,oit->second->outstream);
          break;
#else
          msg.Error()<<"Error in "<<METHOD<<": HepMC format can only be created when Sherpa "
                      <<"was linked with CLHEP, please read our Howto to fix this."<<std::endl;
          abort();
#endif
        case iotype::OldHepMC:
#ifdef USING__CLHEP
          p_hepmc->Sherpa2HepMC(blobs);
          p_hepmc->PrintEvent(2,oit->second->outstream);
          break;
#else
          msg.Error()<<"Error in "<<METHOD<<": HepMC format can only be created when Sherpa "
                      <<"was linked with CLHEP, please read our Howto to fix this."<<std::endl;
          abort();
#endif
        case iotype::HepEvt:
          p_hepevt->Sherpa2HepEvt(blobs);
          p_hepevt->PrintEvent(1,oit->second->outstream,p_hepevt->Nhep());
          break;
        case iotype::D0HepEvt:
          p_hepevt->Sherpa2HepEvt(blobs);
          p_hepevt->PrintEvent(2,oit->second->outstream,p_hepevt->Nhep());
          break;
       
        default:
          msg.Error()<<"Error in "<<METHOD<<std::endl
                     <<"   Unknown Output format : "<<oit->first<<std::endl
                   <<"   No output, continue run ..."<<std::endl;
          break;
        }
      }
    }
  }
  
  m_evtnumber++;
  m_evtcount++;
  
  if (m_evtcount%m_filesize==0) {
    m_evtcount = 0;
    std::string number(ToString(int(m_evtnumber/m_filesize)));
    if (!m_outmap.empty()) {
      for (std::map<iotype::code,NameStream *>::iterator oit=m_outmap.begin();
           oit!=m_outmap.end();oit++) {
        if (oit->second!=NULL) {
          std::string newfilename=oit->second->basicfilename+"."+number+oit->second->fileextension;
          if(oit->first==iotype::Sherpa) oit->second->outstream<<newfilename<<std::endl;
          oit->second->outstream.close();
          oit->second->outstream.open(newfilename.c_str());
          if (!oit->second->outstream.good()) { 
            msg.Error()<<"ERROR in Input_Output_Handler."<<std::endl
                       <<"   Could not open event file "
                       <<(oit->second->basicfilename+"."+number+oit->second->fileextension)
                       <<"."<<std::endl
                       <<"   Will abort the run."<<std::endl;
            abort();
          }
        }
      }
    }
  }
  return true;
}


bool Input_Output_Handler::InputFromFormat(ATOOLS::Blob_List *const blobs) 
{
  if (m_on==false) return false;
  if (!(m_io&4)) return false;
  switch (m_intype) {
  case iotype::Sherpa: return SherpaInput(blobs); 
#ifdef USING__CLHEP
  case iotype::HepMC: 
    THROW(not_implemented,"Reading input from HepMC is not yet possible.");
#endif
  case iotype::HepEvt: return p_hepevt->HepEvt2Sherpa(blobs); 
  default:
    msg.Error()<<"Error in Input_Output_Handler::InputFromFormat."<<std::endl
	       <<"   Unknown Input format : "<<m_intype<<std::endl
	       <<"   No input, continue run ... ."<<std::endl;
    break;
  }
  return false;
}

/*------------------------------------------------------------------
  Sherpa-specific I/O methods : Output is ASCII-format
  ------------------------------------------------------------------*/

void Input_Output_Handler::SherpaOutput(std::ostream & outstream,
					ATOOLS::Blob_List *const blobs,const double weight) 
{ 
  if (m_on==false) return;
  ATOOLS::Particle_Int_Map           P2I;
  ATOOLS::Particle_Int_Map::iterator P2Iiter;
  
  for (Blob_List::iterator blit=blobs->begin();blit!=blobs->end();++blit) {
    for (int i=0;i<(*blit)->NInP();i++) {
      if (P2I.find((*blit)->InParticle(i))==P2I.end()) 
	P2I.insert(std::make_pair((*blit)->InParticle(i),P2I.size()+1));
    }
    for (int i=0;i<(*blit)->NOutP();i++) {
      if (P2I.find((*blit)->OutParticle(i))==P2I.end()) 
	P2I.insert(std::make_pair((*blit)->OutParticle(i),P2I.size()+1));
    }
  }
  outstream<<m_evtnumber<<" "<<P2I.size()<<" "<<blobs->size()<<" "<<weight<<std::endl;
  Particle * part;
  int kfc;
  for (P2Iiter=P2I.begin();P2Iiter!=P2I.end();P2Iiter++) {
    part = P2Iiter->first;
    kfc  = part->Flav().Kfcode(); if (part->Flav().IsAnti()) kfc=-kfc;
    outstream<<P2Iiter->second<<" "<<part->Status()<<" "<<part->Info()<<" "<<kfc<<" "
	       <<" "<<part->Momentum()[0]<<" "<<part->Momentum()[1]
	       <<" "<<part->Momentum()[2]<<" "<<part->Momentum()[3]<<" \n";
  }
  for (Blob_List::iterator blit=blobs->begin();blit!=blobs->end();++blit) {
    outstream<<(*blit)->Id()<<" "<<(*blit)->Status()<<" "<<(int)(*blit)->Type()<<" "<<(*blit)->TypeSpec()
	       <<" "<<(*blit)->NInP()<<" "<<(*blit)->NOutP()<<" \n"
	       <<" "<<(*blit)->Position()[0]<<" "<<(*blit)->Position()[1]
	       <<" "<<(*blit)->Position()[2]<<" "<<(*blit)->Position()[3]<<" \n";
    for (int i=0;i<(*blit)->NInP();i++)  outstream<<P2I.find((*blit)->InParticle(i))->second<<" ";
    for (int i=0;i<(*blit)->NOutP();i++) outstream<<P2I.find((*blit)->OutParticle(i))->second<<" ";
    outstream<<" \n";
  }
}

bool Input_Output_Handler::SherpaInput(ATOOLS::Blob_List *const blobs) 
{ 
  return false;

//   if (m_on==false) return false;
//   blobs->clear();

//   m_evtcount++;
//   int panumber, blnumber, weight;
//   (*p_instream)>>m_evtnumber>>panumber>>blnumber>>weight;

//   ATOOLS::Int_Particle_Map           I2P;
//   ATOOLS::Int_Particle_Map::iterator I2Piter;

//   int        paid, status, kfc;
//   double     mom[4];
//   Vec4D      momentum;
//   char       info;
//   Flavour    flav;
//   Particle * part;
  
//   for (int i=0;i<panumber;i++) {
//     (*p_instream)>>paid>>status>>info>>kfc>>mom[0]>>mom[1]>>mom[2]>>mom[3];
//     flav     = Flavour(kf::code(abs(kfc))); if (kfc<0) flav=flav.Bar();
//     momentum = Vec4D(mom[0],mom[1],mom[2],mom[3]);
//     part     = new Particle(paid,flav,momentum);
//     part->SetStatus(part_status::code(status));
//     part->SetFinalMass();
//     part->SetInfo(info);
//     I2P.insert(std::make_pair(paid,part));
//   }

//   int         type, ninp, noutp, blid;
//   std::string typespec;
//   double      pos[4];
//   Vec4D       position;
//   Blob      * blob;

//   for (int i=0;i<blnumber;i++) {
//     (*p_instream)>>blid>>status>>type>>typespec>>ninp>>noutp;
//     (*p_instream)>>pos[0]>>pos[1]>>pos[2]>>pos[3];

//     position = Vec4D(pos[0],pos[1],pos[2],pos[3]);
//     blob     = new Blob(position,blid);
//     blob->SetStatus(blob_status::code(status));
//     blob->SetType((btp::code)type);
//     blob->SetTypeSpec(typespec);
//     for (int i=0;i<ninp;i++) {
//       (*p_instream)>>paid;
//       I2Piter = I2P.find(paid);
//       if (I2Piter==I2P.end()) {
// 	msg.Error()<<"Error in Input_Output_Handler::SherpaInput."<<std::endl
// 		   <<"   Particle with number "<<paid<<" not found in event."<<std::endl
// 		   <<"   Will return false, continue & hope for the best."<<std::endl;
// 	blobs->clear();
// 	return false;
//       }
//       blob->AddToInParticles(I2Piter->second);
//     }
//     for (int i=0;i<noutp;i++) {
//       (*p_instream)>>paid;
//       I2Piter = I2P.find(paid);
//       if (I2Piter==I2P.end()) {
// 	msg.Error()<<"Error in Input_Output_Handler::SherpaInput."<<std::endl
// 		   <<"   Particle with number "<<paid<<" not found in event."<<std::endl
// 		   <<"   Will return false, continue & hope for the best."<<std::endl;
// 	blobs->clear();
// 	return false;
//       }
//       blob->AddToOutParticles(I2Piter->second);
//     }

//     blobs->push_back(blob);
//   }

//   if (m_evtcount%m_filesize==0) {
//     std::string file, filename;
//     (*p_instream)>>file;
//     p_instream->close();
//     filename =  m_path+std::string("/")+file+std::string(".evts"); 
//     delete p_instream;
//     p_instream = new std::ifstream(filename.c_str()); 
//     if (!p_instream->good()) {
//       msg.Error()<<"ERROR in HepEvt_Interface."<<std::endl
// 		 <<"   Event file "<<filename<<" not found."<<std::endl
// 		 <<"   Will abort the run."<<std::endl;
//       abort();
//     }
//     std::string gentype;
//     (*p_instream)>>gentype>>m_filesize;
//     if (gentype!=std::string("Sherpa")) {
//       msg.Error()<<"ERROR in HepEvt_Interface."<<std::endl
// 		 <<"   Generator type:"<<gentype<<" not SHERPA."<<std::endl
// 		 <<"   Cannot guarantee i/o operations. Will abort the run."<<std::endl;
//       abort();
//     }
//     m_evtcount=0;
//   }
//   return true;
}
