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

using namespace HADRONS;
using namespace ATOOLS;
using namespace PHASIC;
using namespace std;

Hadron_Decay_Map::Hadron_Decay_Map(const Mass_Selector* ms) :
  Decay_Map(ms),
  m_fixed_next_tables(0), p_mixinghandler(NULL)
{}

Hadron_Decay_Map::~Hadron_Decay_Map()
{
  for (map<string, Hadron_Decay_Table*>::iterator it=m_fixed_tables.begin();
       it!=m_fixed_tables.end(); ++it) {
    delete it->second;
  }
}

void Hadron_Decay_Map::ReadInConstants(const string& path, const string& file)
{
  m_startmd.clear();                            // clear model
  Data_Reader reader = Data_Reader(" ",";","!","|");
  reader.AddWordSeparator("\t");
  reader.SetAddCommandLine(false);
  reader.AddComment("#");
  reader.AddComment("//");
  reader.SetInputPath(path);
  reader.SetInputFile(file);

  vector<vector<string> > constants;
  if(!reader.MatrixFromFile(constants)) {
    msg_Error()<<"Warning! The file "<<path<<file<<" does not exist"<<endl
             <<"     or has some syntax error."<<endl;
    msg_Error()<<"     Will ignore it and hope for the best."<<endl;
    return;
  }

  for (size_t i=0;i<constants.size();++i) {
    if( constants[i][1] == "=" ) {              // <name> = <value>
      m_startmd[constants[i][0]] = ToType<double> (
          reader.Interpreter()->Interprete(constants[i][2]) );
    }
  }
  
  set<Flavour> neutral_mesons;
  neutral_mesons.insert(Flavour(kf_K));
  neutral_mesons.insert(Flavour(kf_D));
  neutral_mesons.insert(Flavour(kf_B));
  neutral_mesons.insert(Flavour(kf_B_s));
  set<Flavour>::const_iterator flavit;
  for(flavit = neutral_mesons.begin(); flavit!=neutral_mesons.end(); flavit++) {
    GeneralModel::iterator yit(m_startmd.find("y_"+flavit->IDName()));
    GeneralModel::iterator xit(m_startmd.find("x_"+flavit->IDName()));
    GeneralModel::iterator qit(m_startmd.find("qoverp2_"+flavit->IDName()));
    if(yit!=m_startmd.end()) flavit->SetDeltaGamma(2.0*flavit->Width()*yit->second);
    if(xit!=m_startmd.end()) flavit->SetDeltaM(flavit->Width()*xit->second);
    if(qit!=m_startmd.end()) flavit->SetQOverP2(qit->second);
  }
}

void Hadron_Decay_Map::ReadInPartonicDecays(const ATOOLS::Flavour & decflav,
					    const std::string& path, 
					    const std::string& file) 
{
  if (decflav!=Flavour(kf_b) && decflav!=Flavour(kf_c)) {
    msg_Error()<<"Error in "<<METHOD<<":\n"
	       <<"   No structure to read in partonic decays of "
	       <<decflav<<".\n"
	       <<"   Will continue and hope for the best.\n";
    return;
  }
  Data_Reader reader = Data_Reader(" ",";","!","|");
  reader.AddWordSeparator("\t");
  reader.SetAddCommandLine(false);
  reader.AddComment("#");
  reader.AddComment("//");
  reader.SetInputPath(path);
  reader.SetInputFile(file);
  vector<vector<string> > helpsvv;
  vector<int>             helpkfc;
  if(!reader.MatrixFromFile(helpsvv)) {
    msg_Error()<<"ERROR in "<<METHOD<<endl
      <<"   Read in failure "<<path<<file<<", will abort."<<endl;
    abort();
  }
  Flavour flav;
  double  width(0.),BR,dBR;
  std::string origin;
  Decay_Table * dt;
  if (decflav==Flavour(kf_b))
    dt = Tools::partonic_b;
  else if (decflav==Flavour(kf_c))
    dt = Tools::partonic_c;
  double decmass = decflav.HadMass(); 
  
  for (size_t i=0;i<helpsvv.size();i++) {
    if (helpsvv[i].size()==2) {
      if (helpsvv[i][0]=="Width:") width = ToType<double>(helpsvv[i][1]);
    }
    else if (helpsvv[i].size()==3) {
      if (Tools::ExtractFlavours(helpkfc,helpsvv[i][0])) {
	double mass(decmass);
	for (size_t j=0;j<helpkfc.size();++j) 
	  mass -= Flavour(abs(helpkfc[j])).HadMass();
	if (mass<0.) {
	  msg_Tracking()<<METHOD<<": decay "<<decflav<<" --> ";
	  for (size_t j=0;j<helpkfc.size();++j) 
	    msg_Tracking()<<(helpkfc[j]>0?Flavour(helpkfc[j]):
			Flavour(helpkfc[j]).Bar())<<" ";
	  msg_Tracking()<<"cannot be initialised - not enough mass.\n";
	  continue;
	}
	Decay_Channel * dc(new Decay_Channel(decflav,p_ms));
	for (size_t j=0;j<helpkfc.size();++j) {
	  flav = Flavour(abs(helpkfc[j]));
	  if (helpkfc[j]<0) flav = flav.Bar();
	  dc->AddDecayProduct(flav,false);
	}
	Tools::ExtractBRInfo(helpsvv[i][1], BR, dBR, origin);
	dc->SetWidth(BR*width);
	dt->AddDecayChannel(dc);
	msg_Tracking()<<METHOD<<" adds "<<(*dc)<<"\n";
      }
    }
  }
  dt->UpdateWidth();
  msg_Tracking()<<om::red<<"Read in partonic "<<decflav<<"-decays "
		<<"from "<<file<<".  Found "<<dt->size()
		<<" channels.\n"<<om::reset;
}


void Hadron_Decay_Map::ReadHadronAliases(const string& path, const string& file)
{
  Data_Reader reader = Data_Reader("->", ";", "#", "");
  reader.AddWordSeparator("\t");
  reader.SetAddCommandLine(false);
  reader.AddComment("#");
  reader.AddComment("//");
  reader.SetInputPath(path);
  reader.SetInputFile(file);
  
  vector<vector<string> > aliases;
  reader.MatrixFromFile(aliases);
  
  for (size_t i=0;i<aliases.size();++i) {
    if (aliases[i].size()!=2) {
      msg_Error()<<METHOD<<": Wrong syntax in hadron alias file."<<endl
          <<"  "<<aliases[i]<<endl;
    }
    kf_code alias = ToType<kf_code>(aliases[i][0]);
    kf_code real = ToType<kf_code>(aliases[i][1]);
    m_startmd.m_aliases[alias]=real;
    Particle_Info* aliasinfo = new Particle_Info(*s_kftable[real]);
    aliasinfo->m_kfc=alias;
    s_kftable[alias]=aliasinfo;
    msg_Info()<<METHOD<<" created alias "<<alias<<" for "<<Flavour(alias)<<endl;
  }
}

void Hadron_Decay_Map::Read(const string& path, const string& file, bool verify)
{
  Data_Reader reader = Data_Reader(" ",";","!","->");
  reader.AddWordSeparator("\t");
  reader.SetAddCommandLine(false);
  reader.AddComment("#");
  reader.AddComment("//");
  reader.SetInputPath(path);
  reader.SetInputFile(file);
  
  vector<vector<string> > Decayers;
  if(reader.MatrixFromFile(Decayers)) {
    msg_Info()<<METHOD<<":"
              <<"   Initializing "<<file<<". This may take some time."
              <<endl;
  }
  else {
    if (verify) {
      THROW(fatal_error, "Could not read from DECAYFILE="+file);
    }
  }
  
  Flavour fl;
  std::string createbooklet;
  for (size_t i=0;i<Decayers.size();++i) {
    vector<string> line = Decayers[i];
    if( line.size()>0 && line[0] == string("CREATE_BOOKLET") ) {
      createbooklet = line.size()>1?line[1]:"hadrons.tex";
    }
    else if (line.size()==3) {
      int decayerkf = atoi((line[0]).c_str());
      Flavour decayerflav = Flavour( (kf_code) abs(decayerkf), decayerkf<0);
      msg_Tracking()<<"New hadron decay table for "<<decayerflav<<", "
		    <<"with mass = "<<decayerflav.HadMass()<<"\n";
      Hadron_Decay_Table * dt = new Hadron_Decay_Table(decayerflav, p_ms,
                                                       p_mixinghandler);
      dt->Read(path+line[1], line[2]);
      // add decayer to decaymap
      Decay_Map::iterator it = find(decayerflav);
      if (it==end()) {
        insert(make_pair(decayerflav, vector<Decay_Table*>(1,dt)));
      }
      else {
        it->second.push_back(dt);
        m_counters.insert(make_pair(decayerflav,0));
      }
    }
    else {
      PRINT_INFO("Warning: Decay table in "<<path<<" / "<<file
                 <<" contains incorrect line: "<<endl<<line);
    }
  }
  //createbooklet = "hadrons.tex";
  if(!createbooklet.empty()) {
    Initialise();
    CreateBooklet(createbooklet);
    THROW(normal_exit, string("Created HADRONS++ booklet. ")
          +"Run 'latex hadrons.tex' for compilation.");
  }
}


void Hadron_Decay_Map::ReadFixedTables(const string& path, const string& file)
{
  Data_Reader reader = Data_Reader(" ",";","!","->");
  reader.AddWordSeparator("\t");
  reader.SetAddCommandLine(false);
  reader.AddComment("#");
  reader.AddComment("//");
  reader.SetInputPath(path);
  reader.SetInputFile(file);
  
  vector<vector<string> > Decayers;
  if(!reader.MatrixFromFile(Decayers)) {
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


void Hadron_Decay_Map::Initialise()
{
  for (Decay_Map::iterator pos = this->begin(); pos != this->end(); ++pos) {
    for(size_t i=0; i<pos->second.size(); i++) {
      Hadron_Decay_Table* dt = (Hadron_Decay_Table*) pos->second[i];
      dt->Initialise(m_startmd);
    }
  }
}


void Hadron_Decay_Map::CreateBooklet(std::string & name)
{
  ofstream f(name.c_str());
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
    Hadron_Decay_Table* dt=(Hadron_Decay_Table*) pos->second[0];
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
