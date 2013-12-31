#include "HADRONS++/Main/Hadron_Decay_Table.H"
#include "HADRONS++/Main/Hadron_Decay_Channel.H"
#include "HADRONS++/Main/Mixing_Handler.H"
#include "HADRONS++/Main/Tools.H"
#include "ATOOLS/Phys/Blob.H"
#include "ATOOLS/Org/Data_Reader.H"
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

  bool haspartonics(false);
  std::vector<int> specs;
  std::vector<double> specweights;
  for (size_t i=0;i<helpsvv.size();i++) {
    if (helpsvv[i].size()!=3 && helpsvv[i][0]==string("SPECTATORS")) {
      haspartonics = true;
      Tools::ExtractSpecs(helpsvv[i][1],specs,specweights);
      msg_Tracking()<<METHOD<<": found spectators for "<<m_flin<<": ";
      for (size_t j=0;j<specs.size();j++) msg_Tracking()<<specs[j]<<" ";
      msg_Tracking()<<"\n";
    }
  }

  int nchannels(0), rewrite(0);
  vector<int>      helpkfc;
  double           BR, dBR, totBR(0.);
  string           origin;
  Flavour          flav;
  Hadron_Decay_Channel* hdc;
  for (size_t i=0;i<helpsvv.size();i++) {
    if ( helpsvv[i][0] == string("NO_ANTI") )
      continue;
    if (helpsvv[i].size()>1 && Tools::ExtractFlavours(helpkfc,helpsvv[i][0])) {
      Tools::ExtractBRInfo(helpsvv[i][1], BR, dBR, origin);
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
      if (mass<0.) { 
	msg_Tracking()<<"Found too low mass.";
	BR = 0.; dBR = 0.; continue; 
      }
      if (haspartonics) totBR += BR;
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
  if (Flav().IsC_Hadron()) 
    msg_Tracking()<<METHOD<<" for "<<om::green<<Flav()<<om::reset<<": "
	     <<"total BR = "<<om::green<<totBR<<om::reset<<".\n";
  if (haspartonics) {
    PHASIC::Decay_Table * dectable(NULL);
    if (!Flav().IsB_Hadron() && !Flav().IsC_Hadron()) {
      msg_Error()<<"ERROR in "<<METHOD<<":\n"
		 <<"   No suitable partonic decay table found for "
		 <<Flav()<<".\n"
		 <<"   Will continue and hope for the best.\n";
      if(rewrite) {
	Move(path+file, path+"."+file+".old");
	ofstream ostr( (path + file).c_str() );
	Write(ostr);
	ostr.close();
      }
      ScaleToWidth();
      return;
    }
    double  totspec(0.);
    Flavour spec;
    for (size_t k=0;k<specs.size();k++) totspec+=specweights[k];
    for (size_t k=0;k<specs.size();k++) {
      bool isAnti(false);
      spec = Flavour(abs(specs[k]));
      if (specs[k]<0) spec = spec.Bar();
      if ((spec.IsQuark() && !spec.IsAnti()) ||
	  (spec.IsDiQuark() && spec.IsAnti())) isAnti=true;
      if (Flav()==Flavour(kf_B_c)) {
	msg_Tracking()<<METHOD<<"("<<Flav()<<"): spectator = "<<spec
		 <<" --> anti = "<<isAnti<<".\n";
	if (abs(specs[k])==5)         dectable = Tools::partonic_c;
	else if (abs(specs[k])==4)    dectable = Tools::partonic_b;
	else {
	  msg_Tracking()<<"WARNING in "<<METHOD<<" for "<<Flav()<<":\n"
		   <<"   No annihilation table yet.  Will continue.\n";
	  continue;
	}
      }
      else {
	if (Flav().IsB_Hadron())      dectable = Tools::partonic_b;
	else if (Flav().IsC_Hadron()) dectable = Tools::partonic_c;
      }
      msg_Tracking()<<"Total hadronic width for "<<Flav()<<" = "<<totBR<<".\n";
      double  partWidth((1.-totBR)*Flav().Width()/
			(dectable->TotalWidth()*totspec));
      for (size_t i=0;i<dectable->size();i++) {
	BR = ((*dectable)[i]->Width()*specweights[k]);
	int charge = Flav().IntCharge();
	double mass = Flav().HadMass();
	charge-=spec.IntCharge();
	mass-=spec.HadMass();
	for (size_t j=0;j<(*dectable)[i]->NOut();j++) {
	  mass -= ((*dectable)[i]->GetDecayProduct(j)).HadMass();
	}
	if (mass<0.) {
	  msg_Tracking()<<"Mass too low for "<<Flav()<<" -->";
	  for (size_t j=0;j<(*dectable)[i]->NOut();j++) {
	    msg_Tracking()<<" "<<(*dectable)[i]->GetDecayProduct(j);
	  }
	  msg_Tracking()<<".\n";
	  continue;
	}
	hdc = new Hadron_Decay_Channel(Flav(),p_ms,path);
	hdc->AddDecayProduct(spec,false);
	std::string filename=m_flin.IDName();
	if (filename.find("_{s}")!=string::npos) 
	  filename = StringReplace(filename, "_{s}", "s");
	if (filename.find("-_{b}")!=string::npos) 
	  filename = StringReplace(filename, "-_{b}", "b-");
	if (filename.find("_{b}")!=string::npos) 
	  filename = StringReplace(filename, "_{b}", "b");
	if (filename.find("++_{c}")!=string::npos) 
	  filename = StringReplace(filename, "++_{c}", "c++");
	if (filename.find("+_{c}")!=string::npos) 
	  filename = StringReplace(filename, "+_{c}", "c+");
	if (filename.find("_{c}")!=string::npos) 
	  filename = StringReplace(filename, "_{c}", "c");
	if (filename.find("_fict")!=string::npos) 
	  filename = StringReplace(filename, "_fict", "");
	filename += "_"+spec.IDName();
	msg_Tracking()<<"   Add partonic decay: "<<Flav()<<" --> ";
	for (size_t j=0;j<(*dectable)[i]->NOut();j++) {
	  flav = (*dectable)[i]->GetDecayProduct(j);
	  if (isAnti) flav=flav.Bar();
	  msg_Tracking()<<flav<<" ";
	  hdc->AddDecayProduct(flav,false);
	  charge   -= flav.IntCharge();
	  mass     -= flav.HadMass();
	  filename += flav.IDName();
	}
	hdc->SetWidth(BR*partWidth);
	hdc->SetDeltaWidth(0.);
	hdc->SetOrigin("");
	filename += ".dat";
	msg_Tracking()<<"  ---> "<<filename<<".\n";
	hdc->SetFileName(filename);
	AddDecayChannel(hdc);
	nchannels++;
      }
    }
  }


  if(rewrite) {
    My_Out_File ostr(path + file);
    ostr.Open();
    Write(*ostr);
    ostr.Close();
  }
  ScaleToWidth();
  if (Flav().IsC_Hadron()) 
    msg_Tracking()<<"   --> after rescaling: "
		  <<om::green<<totBR<<om::reset<<".\n";
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

void Hadron_Decay_Table::ScaleToWidth() {
  if(m_flin.Width()/m_totalwidth!=1.0) {
    double delta_tot(0.0);
    for (size_t i=0;i<size();i++)
      if (at(i)->Active()>=0)
        delta_tot+=at(i)->DeltaWidth();
    if (delta_tot>0.0) {
      for (size_t i=0;i<size();i++) {
        if (at(i)->Active()>=0) {
          double scale_fac=at(i)->DeltaWidth()/delta_tot;
          at(i)->SetWidth(at(i)->Width()+
                          scale_fac*(m_flin.Width()-m_totalwidth));
        }
      }
      UpdateWidth();
    }
  }
}

Decay_Channel * Hadron_Decay_Table::Select(Blob* blob) const
{
  Blob_Data_Base* data = (*blob)["dc"];
  //msg_Tracking()<<METHOD<<" for "<<data<<" and flag "<<blob->Status()<<"\n"
  //	   <<(*blob)<<"\n";
  if (data) {
    if (blob->Has(blob_status::internal_flag)) {
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
      //msg_Tracking()<<METHOD<<": erasing "
      //	       <<data->Get<Decay_Channel*>()->Name()<<",\n"
      //	       <<"   retrying with "<<dc->Name()<<".\n";
      DEBUG_INFO("retrying with "<<dc->Name());
      blob->UnsetStatus(blob_status::internal_flag);
      blob->AddData("dc",new Blob_Data<Decay_Channel*>(dc));
      return dc;
    }
    return data->Get<Decay_Channel*>();
  }
  
  Decay_Channel* dec_channel=p_mixinghandler->Select(blob->InParticle(0),*this);

  blob->AddData("dc",new Blob_Data<Decay_Channel*>(dec_channel));
  return dec_channel;
}

