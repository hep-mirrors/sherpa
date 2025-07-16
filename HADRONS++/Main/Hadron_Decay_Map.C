#include "HADRONS++/Main/Tools.H"
#include "HADRONS++/Main/Hadron_Decay_Map.H"
#include "HADRONS++/Main/Hadron_Decay_Table.H"
#include "HADRONS++/Main/Hadron_Decay_Channel.H"
#include "HADRONS++/PS_Library/HD_PS_Base.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "HADRONS++/ME_Library/HD_ME_Base.H"
#include "HADRONS++/ME_Library/Current_ME.H"
#include "HADRONS++/Current_Library/Current_Base.H"
#include "ATOOLS/Org/Getter_Function.H"
#include "HADRONS++/Main/Mixing_Handler.H"
#include "ATOOLS/Org/My_MPI.H"
#include "ATOOLS/Phys/KF_Table.H"
#include "ATOOLS/Org/Scoped_Settings.H"

using namespace HADRONS;
using namespace ATOOLS;
using namespace PHASIC;
using namespace std;

Hadron_Decay_Map::Hadron_Decay_Map(const Mass_Selector* ms) :
  Decay_Map(ms),
  m_fixed_next_tables(0), p_mixinghandler(NULL)
{
}

Hadron_Decay_Map::~Hadron_Decay_Map()
{
  for (map<string, Hadron_Decay_Table*>::iterator it=m_fixed_tables.begin();
       it!=m_fixed_tables.end(); ++it) {
    delete it->second;
  }
}

void Hadron_Decay_Map::Read(Scoped_Settings& s)
{
  // constants
  m_startmd.clear();
  for (const auto& key: s["Constants"].GetKeys()) {
    m_startmd[key] = s["Constants"][key].SetDefault(-1.0).Get<double>();
  }

  // aliases
  for (const auto& key: s["Aliases"].GetKeys()) {
    kf_code alias = ToType<kf_code>(key);
    kf_code real = s["Aliases"][key].SetDefault(0).Get<kf_code>();
    Tools::aliases[alias]=real;
    Particle_Info* aliasinfo = new Particle_Info(*s_kftable[real]);
    aliasinfo->m_kfc=alias;
    s_kftable[alias]=aliasinfo;
    msg_Info()<<METHOD<<" created alias "<<alias<<" for "<<Flavour(alias)<<endl;
  }

  // partonic decay tables for b and c quark
  ReadInPartonicDecays(Flavour(kf_b),s);
  ReadInPartonicDecays(Flavour(kf_c),s);

  // hadron decay tables
  Flavour fl;
  for (const auto& decayer: s["Channels"].GetKeys()) {
    Flavour decayerflav( (kf_code) abs(atoi(decayer.c_str())), atoi(decayer.c_str())<0);
    Hadron_Decay_Table * dt = new Hadron_Decay_Table(decayerflav, p_ms,
                                                     p_mixinghandler);
    dt->Read(s["Channels"][decayer], m_startmd);
    // add decayer to decaymap
    Decay_Map::iterator it = find(decayerflav);
    if (it==end()) {
      insert(make_pair(decayerflav, dt));
    }
    else {
      THROW(fatal_error, "Duplicate decay table for "+decayerflav.IDName());
    }
  }

  // create booklet
  if (s["CreateBooklet"].SetDefault(false).Get<bool>()) {
    CreateBooklet();
    THROW(normal_exit, string("Created HADRONS++ booklet. ")
          +"Run 'latex hadrons.tex' for compilation.");
  }
}


void Hadron_Decay_Map::ReadInPartonicDecays(const ATOOLS::Flavour & decflav, Scoped_Settings& ss)
{
  Flavour flav;
  std::string origin;
  Decay_Table * dt=NULL;
  if (decflav==Flavour(kf_b))
    dt = Tools::partonic_b;
  else if (decflav==Flavour(kf_c))
    dt = Tools::partonic_c;
  else
    THROW(fatal_error, "Internal error.");
  
  auto s = ss["Partonics"][decflav.IDName()];
  double width = s["Width"].SetDefault(-1.0).Get<double>();
  for (const auto& channel: s["Channels"].GetKeys()) {
    vector<int> helpkfc;
    Tools::ExtractFlavours(helpkfc,channel);
    Decay_Channel * dc(new Decay_Channel(decflav,p_ms));
    for (size_t j=0;j<helpkfc.size();++j) {
      flav = Flavour(abs(helpkfc[j]));
      if (helpkfc[j]<0) flav = flav.Bar();
      dc->AddDecayProduct(flav,false);
    }
    dc->SetWidth(s["Channels"][channel]["BR"][0].SetDefault(-1.0).Get<double>()*width);
    dc->SetDeltaWidth(s["Channels"][channel]["BR"][1].SetDefault(-1.0).Get<double>()*width);
    dt->AddDecayChannel(dc);
    msg_Tracking()<<METHOD<<" adds "<<(*dc)<<"\n";
  }
  dt->UpdateWidth();
  msg_Tracking()<<om::red<<"Read in partonic "<<decflav<<"-decays. Found "<<dt->size()
		<<" channels.\n"<<om::reset;
}


void Hadron_Decay_Map::ReadFixedTables()
{
  /* TODO
  Data_Reader reader = Data_Reader(" ",";","!","->");
  reader.AddWordSeparator("\t");
  reader.AddComment("#");
  reader.AddComment("//");
  reader.SetInputPath(path);
  reader.SetInputFile(file);
  reader.AddLineSeparator("\n");
  
  vector<vector<string> > Decayers;
  if(!reader.MatrixFromFile(Decayers,"")) {
    return;
  }

  Flavour fl;
  for (size_t i=0;i<Decayers.size();++i) {
    vector<string> line = Decayers[i];
    if (line.size()==4) {
      std::string table_id = line[0];
      int decayerkf = atoi((line[1]).c_str());
      Flavour decayerflav = Flavour( (kf_code) abs(decayerkf), decayerkf<0);
      Hadron_Decay_Table * dt = new Hadron_Decay_Table(decayerflav, p_ms,
                                                       p_mixinghandler);
      dt->Read(path+line[2], line[3]);
      pair<SDtMMapIt, SDtMMapIt> found=m_fixed_tables.equal_range(table_id);
      for (SDtMMapIt it=found.first; it!=found.second; ++it) {
        if (it->second->Flav()==decayerflav) {
          THROW(fatal_error, "Duplicate decayer "+ToString((long int)decayerflav)
                +" for fixed decay table ID="+table_id);
        }
      }
      m_fixed_tables.insert(make_pair(table_id, dt));
    }
    else {
      msg_Error()<<METHOD<<" Invalid line in FixedDecays.dat:"<<endl
                 <<"  "<<line<<endl<<"Ignoring it."<<endl;
      
    }
  }
  for (map<string, Hadron_Decay_Table*>::iterator it=m_fixed_tables.begin();
       it!=m_fixed_tables.end(); ++it) {
    it->second->Initialise(m_startmd);
  }
  */
}


void Hadron_Decay_Map::FixDecayTables(std::string table_id)
{
  pair<SDtMMapIt, SDtMMapIt> found=m_fixed_tables.equal_range(table_id);
  for (SDtMMapIt it=found.first; it!=found.second; ++it) {
    m_fixed_next_tables.push_back(it->second);
  }
}


void Hadron_Decay_Map::ClearFixedDecayTables()
{
  m_fixed_next_tables.clear();
}


void Hadron_Decay_Map::CreateBooklet()
{
  ofstream f("hadrons.tex");
  // header
  f<<"\\documentclass[a4paper]{scrartcl}\n"
   <<"\\usepackage{latexsym,amssymb,amsmath,amsxtra,longtable,fullpage}\n"
   <<"\\usepackage[ps2pdf,colorlinks,bookmarks=true,bookmarksnumbered=true]{hyperref}\n\n"
   <<"\\begin{document}\n"<<endl; 
  f<<"\\newcommand{\\m}{-}"<<endl;
  f<<"\\setlength{\\parindent}{0pt}"<<endl;
  f<<"\\newcommand{\\p}{+}"<<endl; 
  f<<"\\newcommand{\\mytarget}[1]{\\hypertarget{#1}{#1}}"<<endl;
  f<<"\\newcommand{\\mylink}[1]{\\hyperlink{#1}{#1}}"<<endl;
  f<<"\\title{Available Matrix Elements and Decay Channels of the "
   <<"{\\tt HADRONS++} Module}\n\\maketitle"<<endl;
  f<<"\\tableofcontents"<<endl<<endl;

  // MEs
  std::string indent="  \\subsubsection{ ";
  std::string separator=" } \n";
  std::string lineend=" \n";
  std::string replacefrom="_";
  std::string replaceto="\\_";
  f<<"\\section{Available Decay Matrix Elements}"<<endl;
  f<<"\\subsection{Complete Matrix Elements}"<<endl;
  Getter_Function<HD_ME_Base,ME_Parameters>::PrintGetterInfo(
    f,30,indent, separator, lineend, replacefrom, replaceto);
  f<<"\\subsection{Weak Currents}"<<endl;
  Getter_Function<Current_Base,ME_Parameters>::PrintGetterInfo(
    f,30,indent, separator, lineend, replacefrom, replaceto);

  // text 
  f<<"\\section{Decay Channels}"<<endl;
  std::vector<HD_ME_Base*> mes;
  for ( Decay_Map::iterator pos = begin(); pos != end(); ++pos) {
    Hadron_Decay_Table* dt=(Hadron_Decay_Table*) pos->second;
    if(dt==NULL) continue;
    dt->LatexOutput(f);
  }
  // end 
  f<<"\\end{document}"<<endl;
  f.close();
}

Decay_Table* Hadron_Decay_Map::FindDecay(const ATOOLS::Flavour & decayer)
{
  // first check, whether a fixed decaytable has been requested for this decayer
  for (size_t i=0; i<m_fixed_next_tables.size(); ++i) {
    if (m_fixed_next_tables[i]->Flav().Kfcode()==decayer.Kfcode()) {
      return m_fixed_next_tables[i];
    }
  }

  return Decay_Map::FindDecay(decayer);
}
