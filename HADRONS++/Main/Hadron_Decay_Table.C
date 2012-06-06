#include "HADRONS++/Main/Hadron_Decay_Table.H"
#include "HADRONS++/Main/Hadron_Decay_Channel.H"
#include "HADRONS++/Main/Mixing_Handler.H"
#include "ATOOLS/Phys/Blob.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Shell_Tools.H"

using namespace HADRONS;
using namespace ATOOLS;
using namespace PHASIC;
using namespace std;

Hadron_Decay_Table::Hadron_Decay_Table(Flavour decayer, const Mass_Selector* ms,
                                       Mixing_Handler* mh) :
  Decay_Table(decayer, ms), p_mixinghandler(mh)
{
}

Hadron_Decay_Table::~Hadron_Decay_Table()
{
}


void Hadron_Decay_Table::Read(std::string path, std::string file)
{
  Data_Reader reader = Data_Reader("|",";","!");
  reader.SetAddCommandLine(false);
  reader.AddComment("#");
  reader.AddComment("//");
  reader.SetInputPath(path);
  reader.SetInputFile(file);
  // read decay table file
  vector<vector<string> > helpsvv;
  if(!reader.MatrixFromFile(helpsvv)) {
    msg_Error()<<"ERROR in "<<METHOD<<endl
      <<"   Read in failure "<<path<<file<<", will abort."<<endl;
    abort();
  }

  int nchannels(0), rewrite(0);
  vector<int>      helpkfc;
  double           BR, dBR;
  string           origin;
  Flavour          flav;
  Hadron_Decay_Channel* hdc;
  
  for (size_t i=0;i<helpsvv.size();i++) {
    if ( helpsvv[i][0] == string("NO_ANTI") )
      continue;
    if (ExtractFlavours(helpkfc,helpsvv[i][0])) {
      ExtractBRInfo( helpsvv[i][1], BR, dBR, origin );
      hdc = new Hadron_Decay_Channel(Flav(),p_ms,path);
      int charge = Flav().IntCharge();
      double mass = Flav().HadMass();
      for (size_t j=0;j<helpkfc.size();++j) {
        flav = Flavour(abs(helpkfc[j]));
        if (helpkfc[j]<0) flav = flav.Bar();
        hdc->AddDecayProduct(flav);
	charge-=flav.IntCharge();
	mass-=flav.HadMass();
      }
      hdc->SetWidth(BR*Flav().Width());
      hdc->SetDeltaWidth(dBR*Flav().Width());
      hdc->SetOrigin(origin);
      if(helpsvv[i].size()==3) hdc->SetFileName(helpsvv[i][2]);
      else {
        hdc->SetFileName();
        rewrite=1;
      }
      if(charge!=0)
	THROW(fatal_error,"Charge not conserved for "+hdc->FileName());
      if(mass<-Accu())
	THROW(fatal_error,"Decaying mass "+ToString(mass)+" too low in "+
              hdc->FileName());
      AddDecayChannel(hdc);
      nchannels++;
    }
  }

  if(rewrite) {
    Move(path+file, path+"."+file+".old");
    ofstream ostr( (path + file).c_str() );
    Write(ostr);
    ostr.close();
  }
  ScaleToWidth();
}


void Hadron_Decay_Table::Initialise(GeneralModel& startmd)
{
  if(size()==0) {
    msg_Error()<<"WARNING in "<<METHOD<<": "<<endl
      <<"   No decay channels found for "<<Flav()<<endl
      <<"   Will continue and hope for the best."<<endl;
  }
  else {
    msg_Tracking()<<"Initialising "<<size()
      <<" decay channels for "<<Flav()
      <<" ("<<TotalWidth()/Flav().Width()*100.0<<"%)"<<endl;
    if(msg_LevelIsDebugging()) Output();
  }

  Hadron_Decay_Channel* hdc;
  for (size_t i=0; i<size(); i++) {
    hdc = at(i);
    hdc->Initialise(startmd);
  }
}


void Hadron_Decay_Table::Write(std::ostream& ostr)
{
  ostr<<"# outgoing part. \t | BR(Delta BR) \t [Origin] \t | DC-file\n"<<endl;
  for (size_t j=0; j<size();j++) {
    Hadron_Decay_Channel* hdc = at(j);
    double dBR=hdc->DeltaWidth()/Flav().Width();
    ostr<<"{"<<int(hdc->Flavs()[0]);
    for (size_t k=0; k<hdc->Flavs().size();++k) ostr<<","<<int(hdc->Flavs()[k]);
    ostr<<"}\t | ";
    ostr<<hdc->Width()/Flav().Width();
    if(dBR>0.0)           ostr<<"("<<dBR<<")";
    if(hdc->Origin()!="") ostr<<"["<<hdc->Origin()<<"]";
    ostr<<"\t | ";
    ostr<<hdc->FileName()<<";"<<endl;
  }
}

bool Hadron_Decay_Table::ExtractFlavours(vector<int> & helpkfc,string help)
{
  helpkfc.clear();    
  size_t pos = help.find("{");
  bool             hit;
  if (pos!=string::npos) help = help.substr(pos+1);
  else {
    msg_Error()<<"WARNING in "<<METHOD<<": "<<endl
           <<"   Something wrong with final state of decay (Bracket missing) :"<<help<<endl
           <<"   Will skip it."<<endl;
    return false;
  }
  pos    = help.find("}");
  if (pos!=string::npos) help = help.substr(0,pos);
  else {
    msg_Error()<<"WARNING in "<<METHOD<<": "<<endl
           <<"   Something wrong with final state of decay (Bracket missing) :"<<help<<endl
           <<"   Will skip it."<<endl;
    return false;
  }
  hit    = true;
  while (hit) {
    pos      = help.find(",");
    if (pos!=string::npos) {
      helpkfc.push_back(atoi((help.substr(0,pos)).c_str()));
      help  = help.substr(pos+1);
    }
    else {
      helpkfc.push_back(atoi(help.c_str()));
      hit = false;
    }
  }
//  if (helpkfc.size()<2) {
//    msg_Error()<<"WARNING in Decay_Map:: : "<<endl
//           <<"   Something wrong with final state of decay (Too little particles) : ";
//    for (int j=0;j<helpkfc.size();j++) msg_Error()<<helpkfc[j]<<" ";
//    msg_Error()<<endl<<"   Will skip it."<<endl;
//    return false;
//  } 
  if (helpkfc.size()<1) {
    msg_Error()<<"WARNING in "<<METHOD<<": "<<endl
           <<"   Something wrong with final state of decay. (no particles?)"<<endl
           <<"   Will skip it and hope for the best."<<endl;
    return false;
  } 
  return true;
}
 
void Hadron_Decay_Table::ExtractBRInfo( string entry, double & br, double & dbr, string & origin )
{
  size_t posa, posb;        // start and end of things b/w brackets
  size_t posmin;            // start of first bracket

  std::string sbr, sdbr;

  // extract Delta BR
  posa = entry.find("(");
  posb = entry.find(")");
  posmin = posa;
  if(posa!=string::npos && posb!=string::npos)
    sdbr = entry.substr(posa+1,posb-posa-1);
  else sdbr = "-1.0";

  // extract Origin
  posa = entry.find("[");
  posb = entry.find("]");
  if(posmin==string::npos || (posmin!=string::npos && posmin>posa)) posmin=posa;
  if(posa!=string::npos && posb!=string::npos)
    origin = entry.substr(posa+1,posb-posa-1);
  else origin = string("");

  // extract BR
  if( posmin!=string::npos ) sbr = entry.substr(0,posmin);
  else                       sbr = entry.substr(0);

  Algebra_Interpreter ip;
  sdbr=ip.Interprete(sdbr);
  sbr=ip.Interprete(sbr);

  dbr=ToType<double>(sdbr);
  br=ToType<double>(sbr);

  if (dbr==-1.0) dbr = br;
}

void Hadron_Decay_Table::LatexOutput(std::ostream& f)
{
  f<<"\\subsection{\\texorpdfstring{Decaying Particle: $"<<Flav().TexName()<<"$"
   <<" ["<<Flav().Kfcode()<<"]}"
   <<"{"<<"["<<Flav().Kfcode()<<"] "<<Flav()<<"}}"<<endl;
  f<<"\\begin{tabular}{ll}"<<endl;
  f<<" number of decay channels:    & "<<size()<<"\\\\ "<<endl;
  f<<" total width:               & "<<TotalWidth()<<" GeV \\\\ "<<endl;
  f<<" experimental width:        & "<<Flav().Width()<<" GeV \\\\ "<<endl;
  f<<"\\end{tabular}"<<endl;
  f<<"\\begin{longtable}[l]{lll}"<<endl;
  f<<"\\multicolumn{3}{c}{\\bf Exclusive Decays}\\\\"<<endl;
  f<<"\\hline"<<endl;
  f<<"Decay Channel & Input BR [Origin]/Integrated BR [Matrix Element]\\\\"<<endl;
  f<<"\\hline\n\\hline"<<endl;
  for(size_t i=0; i<size(); ++i) {
    if(at(i)->Width()!=0.0) at(i)->LatexOutput(f, TotalWidth());
  }
  // skip inclusives for now
  f<<"\\hline"<<endl;
  f<<"\\end{longtable}"<<endl;
}

Decay_Channel * Hadron_Decay_Table::Select(Blob* blob) const
{
  Blob_Data_Base* data = (*blob)["dc"];
  if(data) {
    if(blob->Has(blob_status::internal_flag)) {
      bool partonic_finalstate(false);
      Decay_Channel* dc;
      do {
        dc = Decay_Table::Select();
        for (size_t i=0; i<dc->Flavs().size(); ++i) {
          if(dc->Flavs()[i].Strong()) {
            partonic_finalstate=true;
            break;
          }
        }
      } while (!partonic_finalstate);
      DEBUG_INFO("retrying with "<<dc->Name());
      blob->UnsetStatus(blob_status::internal_flag);
      blob->AddData("dc",new Blob_Data<Decay_Channel*>(dc));
      return dc;
    }
    else return data->Get<Decay_Channel*>();
  }
  
  Decay_Channel* dec_channel=p_mixinghandler->Select(blob->InParticle(0),*this);

  blob->AddData("dc",new Blob_Data<Decay_Channel*>(dec_channel));
  return dec_channel;
}

