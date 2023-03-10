#include "HADRONS++/Main/Hadron_Decay_Table.H"
#include "HADRONS++/Main/Hadron_Decay_Channel.H"
#include "HADRONS++/Main/Mixing_Handler.H"
#include "HADRONS++/Main/Tools.H"
#include "ATOOLS/Phys/Blob.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/Shell_Tools.H"
#include "ATOOLS/Org/My_MPI.H"
#include <algorithm>

using namespace HADRONS;
using namespace ATOOLS;
using namespace PHASIC;
using namespace std;

Hadron_Decay_Table::Hadron_Decay_Table(Flavour decayer, const Mass_Selector* ms,
                                       Mixing_Handler* mh) :
  Decay_Table(decayer, ms), p_mixinghandler(mh)
{
  m_flavwidth=Flav().Width();
  if (Flav().Kfcode()==kf_tau && m_flavwidth==0.0) {
    m_flavwidth = 2.26735e-12;
  }
}

Hadron_Decay_Table::~Hadron_Decay_Table()
{
}

void Hadron_Decay_Table::Read(Scoped_Settings s, GeneralModel& startmd)
{
  DEBUG_FUNC((int) Flav());

  // read "Partonics" block
  vector<int> specs;
  vector<double> specweights;
  for (auto spec: s["Spectators"].GetItems()) {
    specs.push_back(ToType<int>(spec.GetKeys()[0]));
    specweights.push_back(spec[spec.GetKeys()[0]]["Weight"].SetDefault(1).Get<double>());
  }
  for (size_t j=0;j<specs.size();j++) DEBUG_INFO("found spectator for "<<m_flin<<": "<<specs[j]);

  // read "Forced" block
  vector<string> forced_channels(s["Forced"].SetDefault(vector<string>()).GetVector<string>());

  // read all other entries (= real channels)
  double totwidth(0.);
  for (const auto& channel: s.GetKeys()) {
    if (channel=="Spectators" || channel=="Forced") continue;

    DEBUG_VAR(channel);
    vector<int> helpkfc;
    Tools::ExtractFlavours(helpkfc,channel);
    Hadron_Decay_Channel* hdc = new Hadron_Decay_Channel(Flav(),p_ms);
    int charge = Flav().IntCharge();
    double mass = Flav().HadMass();
    for (size_t j=0;j<helpkfc.size();++j) {
      Flavour flav = Flavour(abs(helpkfc[j]));
      if (helpkfc[j]<0) flav = flav.Bar();
      hdc->AddDecayProduct(flav);
      charge-=flav.IntCharge();
      mass-=flav.HadMass();
    }
    if (mass<0.) THROW(fatal_error, "Mass too low.");
    if (charge!=0) THROW(fatal_error,"Charge not conserved for "+hdc->IDCode());
    if (s[channel]["BR"].GetItemsCount()==0) {
      // some partonic channels will be constructed manually below.
      // these will not have a BR provided in the decay data.
      // (they are nevertheless stored in the decay table for the intresults)
      delete hdc;
      continue;
    }
    hdc->SetWidth(s[channel]["BR"][0].SetDefault(-1.0).Get<double>()*m_flavwidth);
    hdc->SetDeltaWidth(s[channel]["BR"][1].SetDefault(-1.0).Get<double>()*m_flavwidth);
    hdc->SetOrigin(s[channel]["Origin"].SetDefault("").Get<string>());
    hdc->SetActive(s[channel]["Status"].SetDefault(hdc->Active()).GetVector<int>());
    totwidth += hdc->Width();
    hdc->Initialise(s[channel], startmd);
    AddDecayChannel(hdc);
  }

  DEBUG_VAR(totwidth);
  if (specs.size()>0) {
    PHASIC::Decay_Table * dectable(NULL);
    if (!Flav().IsB_Hadron() && !Flav().IsC_Hadron()) {
      msg_Error()<<"ERROR in "<<METHOD<<":\n"
		 <<"   No suitable partonic decay table found for "
		 <<Flav()<<".\n"
		 <<"   Will continue and hope for the best.\n";
      return;
    }
    double  totspec(0.);
    for (size_t k=0;k<specs.size();k++) totspec+=specweights[k];
    for (size_t k=0;k<specs.size();k++) {
      bool isAnti(false);
      Flavour spec = Flavour(abs(specs[k]));
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
      msg_Tracking()<<"Total hadronic width for "<<Flav()<<" = "<<totwidth<<".\n";
      double  partWidth((m_flavwidth-totwidth)/
			(dectable->TotalWidth()*totspec));
      for (size_t i=0;i<dectable->size();i++) {
	double BR = ((*dectable)[i]->Width()*specweights[k]);
	double mass = Flav().HadMass();
        string fsidcode;
	mass-=spec.HadMass();
        fsidcode=ToString((int)spec);
	for (size_t j=0;j<(*dectable)[i]->NOut();j++) {
          Flavour flav = (*dectable)[i]->GetDecayProduct(j);
	  if (isAnti) flav=flav.Bar();
	  mass -= flav.HadMass();
          fsidcode=fsidcode+","+ToString((int)flav);
	}
	if (mass<0.) {
	  msg_Tracking()<<"Mass too low for "<<Flav()<<" -->";
	  for (size_t j=0;j<(*dectable)[i]->NOut();j++) {
	    msg_Tracking()<<" "<<(*dectable)[i]->GetDecayProduct(j);
	  }
	  msg_Tracking()<<".\n";
	  continue;
	}
	Hadron_Decay_Channel* hdc = new Hadron_Decay_Channel(Flav(),p_ms);
	hdc->AddDecayProduct(spec,false);
	msg_Tracking()<<"   Add partonic decay: "<<Flav()<<" --> ";
	for (size_t j=0;j<(*dectable)[i]->NOut();j++) {
	  Flavour flav = (*dectable)[i]->GetDecayProduct(j);
	  if (isAnti) flav=flav.Bar();
	  msg_Tracking()<<flav<<" ";
	  hdc->AddDecayProduct(flav,false);
	}
	hdc->SetWidth(BR*partWidth);
	hdc->SetDeltaWidth(0.);
	hdc->SetOrigin("");
        hdc->Initialise(s[fsidcode],startmd);
	AddDecayChannel(hdc);
      }
    }
  }

  // set forced channels
  UpdateChannelStatuses();
}


void Hadron_Decay_Table::LatexOutput(std::ostream& f)
{
  f<<"\\subsection{\\texorpdfstring{Decaying Particle: $"<<Flav().TexName()<<"$"
   <<" ["<<Flav().Kfcode()<<"]}"
   <<"{"<<"["<<Flav().Kfcode()<<"] "<<Flav()<<"}}"<<endl;
  f<<"\\begin{tabular}{ll}"<<endl;
  f<<" number of decay channels:    & "<<size()<<"\\\\ "<<endl;
  f<<" total width:               & "<<TotalWidth()<<" GeV \\\\ "<<endl;
  f<<" experimental width:        & "<<m_flavwidth<<" GeV \\\\ "<<endl;
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

Decay_Channel * Hadron_Decay_Table::Select(Blob* blob)
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

