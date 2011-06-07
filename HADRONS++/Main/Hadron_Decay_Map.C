#include "HADRONS++/Main/Hadron_Decay_Map.H"
#include "HADRONS++/Main/Hadron_Decay_Table.H"
#include "HADRONS++/Main/Hadron_Decay_Channel.H"
#include "HADRONS++/PS_Library/HD_PS_Base.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "HADRONS++/ME_Library/HD_ME_Base.H"
#include "HADRONS++/Current_Library/Current_Base.H"
#include "ATOOLS/Org/Getter_Function.H"

using namespace HADRONS;
using namespace ATOOLS;
using namespace std;

Hadron_Decay_Map::Hadron_Decay_Map(string path,string file,string constfile) :
  m_decaypath(path), m_decayfile(file), m_constfile(constfile),
  m_fixed_next_tables(0)
{
}

Hadron_Decay_Map::~Hadron_Decay_Map()
{
  for (map<string, Hadron_Decay_Table*>::iterator it=m_fixed_tables.begin();
       it!=m_fixed_tables.end(); ++it) {
    delete it->second;
  }
}

void Hadron_Decay_Map::ReadInConstants()
{
  m_startmd.clear();                            // clear model
  Data_Reader reader = Data_Reader(" ",";","!","|");
  reader.AddWordSeparator("\t");
  reader.SetAddCommandLine(false);
  reader.AddComment("#");
  reader.AddComment("//");
  reader.SetInputPath(m_decaypath);
  reader.SetInputFile(m_constfile);

  vector<vector<string> > constants;
  if(!reader.MatrixFromFile(constants)) {
    msg_Error()<<"Warning! The file "<<m_decaypath<<m_constfile<<" does not exist"<<endl
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
}

void Hadron_Decay_Map::SetHadronProperties()
{
  FlavourSet neutral_mesons;
  neutral_mesons.insert(Flavour(kf_K));
  neutral_mesons.insert(Flavour(kf_D));
  neutral_mesons.insert(Flavour(kf_B));
  neutral_mesons.insert(Flavour(kf_B_s));
  FlavourSet::const_iterator flavit;
  for(flavit = neutral_mesons.begin(); flavit!=neutral_mesons.end(); flavit++) {
    GeneralModel::iterator yit(m_startmd.find("y_"+flavit->IDName()));
    GeneralModel::iterator xit(m_startmd.find("x_"+flavit->IDName()));
    GeneralModel::iterator qit(m_startmd.find("qoverp2_"+flavit->IDName()));
    if(yit!=m_startmd.end()) flavit->SetDeltaGamma(2.0*flavit->Width()*yit->second);
    if(xit!=m_startmd.end()) flavit->SetDeltaM(flavit->Width()*xit->second);
    if(qit!=m_startmd.end()) flavit->SetQOverP2(qit->second);
  }
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
              <<"   Initializing from "<<file<<", this may take some time."
              <<endl;
  }
  else {
    if (verify) {
      THROW(fatal_error, "Could not read from DECAYFILE="+file);
    }
  }
  
  Flavour fl;
  bool createbooklet=false;
  for (size_t i=0;i<Decayers.size();++i) {
    vector<string> line = Decayers[i];
    if( line[0] == string("CREATE_BOOKLET") ) {
      createbooklet = true;
    }
    else {
      int decayerkf = atoi((line[0]).c_str());
      Flavour decayerflav = Flavour( (kf_code) abs(decayerkf), decayerkf<0);
      Hadron_Decay_Table * dt = new Hadron_Decay_Table(decayerflav);
      dt->Read(m_decaypath+line[1], line[2]);
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
  }
  if(createbooklet) {
    Initialise();
    CreateBooklet();
  }
}


void Hadron_Decay_Map::ReadFixedTables()
{
  Data_Reader reader = Data_Reader(" ",";","!","->");
  reader.AddWordSeparator("\t");
  reader.SetAddCommandLine(false);
  reader.AddComment("#");
  reader.AddComment("//");
  reader.SetInputPath(m_decaypath);
  reader.SetInputFile("FixedDecays.dat");
  
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
      Hadron_Decay_Table * dt = new Hadron_Decay_Table(decayerflav);
      dt->Read(m_decaypath+line[2], line[3]);
      pair<SDtMMapIt, SDtMMapIt> found=m_fixed_tables.equal_range(table_id);
      for (SDtMMapIt it=found.first; it!=found.second; ++it) {
        if (it->second->Flav()==decayerflav) {
          THROW(fatal_error, "Duplicate decayer "+ToString(decayerflav.HepEvt())
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


Hadron_Decay_Table* Hadron_Decay_Map::FindDecay(const ATOOLS::Flavour & decayer)
{
  // first check, whether a fixed decaytable has been requested for this decayer
  for (size_t i=0; i<m_fixed_next_tables.size(); ++i) {
    if (m_fixed_next_tables[i]->Flav().Kfcode()==decayer.Kfcode()) {
      return m_fixed_next_tables[i];
    }
  }

  Flavour tempdecayer=decayer;
  Decay_Map::iterator it = find(decayer);
  if(it==end()) {
    it = find(decayer.Bar());
    tempdecayer=decayer.Bar();
  }
  if(it==end()) return NULL;

  // there may be multiple decay tables for one flavour, so find the right one
  int count = it->second.size()-1; // default to last decay table available
  map<ATOOLS::Flavour,int>::iterator counterit = m_counters.find(tempdecayer);
  if(counterit!=m_counters.end() &&
     counterit->second < int(it->second.size()-1) )
  {
    count = counterit->second;
    counterit->second++;
  }
  return (Hadron_Decay_Table*) it->second[count];
}


void Hadron_Decay_Map::ResetCounters()
{
  map<ATOOLS::Flavour,int>::iterator it;
  for(it=m_counters.begin(); it!=m_counters.end(); it++) {
    it->second=0;
  }
}


void Hadron_Decay_Map::CreateBooklet()
{
  string fn = string("hadrons.tex");
  ofstream f(fn.c_str());
  
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
  Getter_Function<HD_ME_Base,Flavour_Info>::PrintGetterInfo(
    f,30,indent, separator, lineend, replacefrom, replaceto);
  f<<"\\subsection{Weak Currents}"<<endl;
  Getter_Function<Current_Base,Flavour_Info>::PrintGetterInfo(
    f,30,indent, separator, lineend, replacefrom, replaceto);

  // text 
  f<<"\\section{Decay Channels}"<<endl;
  Hadron_Decay_Channel * dc (NULL);
  std::vector<MEPair> mes; std::vector<CurrentsPair> currents;
  FlavourSet outs;
  char helpstr[100];
  int total_decaychannels(0),total_mes(0);
  for ( Decay_Map::iterator pos = begin(); pos != end(); ++pos) {
    Hadron_Decay_Table* dt=(Hadron_Decay_Table*) pos->second[0];
    if(dt==NULL) continue;
    f<<"\\subsection{\\texorpdfstring{Decaying Particle: $"<<dt->Flav().TexName()<<"$"
     <<" ["<<dt->Flav().Kfcode()<<"]}"
     <<"{"<<"["<<dt->Flav().Kfcode()<<"] "<<dt->Flav()<<"}}"<<endl;
    f<<"\\begin{tabular}{ll}"<<endl;
    f<<" number of decay channels:    & "<<dt->size()<<"\\\\ "<<endl;
    f<<" total width:               & "<<dt->TotalWidth()<<" GeV \\\\ "<<endl;
    f<<" experimental width:        & "<<dt->Flav().Width()<<" GeV \\\\ "<<endl;
    f<<"\\end{tabular}"<<endl;
    f<<"\\begin{longtable}[l]{lll}"<<endl;
    f<<"\\multicolumn{3}{c}{\\bf Exclusive Decays}\\\\"<<endl;
    f<<"\\hline"<<endl;
    f<<"Decay Channel & Input BR [Origin]/Integrated BR [Matrix Element]\\\\"<<endl;
    f<<"\\hline\n\\hline"<<endl;
    for(size_t i=0; i<dt->size(); ++i) {
      dc = dt->at(i);
      if(dc->Width()==0.0) continue;
      total_decaychannels++;
      outs = dc->GetDecayProducts();
      f<<"$"<<dc->GetDecaying().TexName()<<"$ $\\to$ ";
      for (FlSetConstIter fl=outs.begin();fl!=outs.end();++fl) f<<"$"<<fl->TexName()<<"$ ";
      f<<" & ";
      sprintf( helpstr, "%.4f", dc->Width()/dt->Flav().Width()*100. );
      f<<helpstr;
      if( dc->DeltaWidth() > 0. ) {
        sprintf( helpstr, "%.4f", dc->DeltaWidth()/dt->Flav().Width()*100. );
        f<<" $\\pm$ "<<helpstr;
      }
      f<<" \\% ";
      if(dc->Origin()!="") {
        f<<"[\\verb;"<<dc->Origin()<<";]";
      }
      f<<"\\\\"<<endl;
      if((dc->GetMEs().size()>0 && dc->GetMEs()[0].second->Name()!="Generic") ||
          dc->GetCurrents().size()>0 )
      {
        sprintf( helpstr, "%.4f", dc->GetPS()->Result()/dt->Flav().Width()*100. );
        f<<" & "<<helpstr;
        if( dc->GetPS()->Error() > 0. ) {
          sprintf( helpstr, "%.4f", dc->GetPS()->Error()/dt->TotalWidth()*100. );
          f<<" $\\pm$ "<<helpstr;
        }
        f<<" \\% ";
      }
      mes=dc->GetMEs();
      for(size_t i=0;i<mes.size();i++) {
        if(mes[i].second->Name()!="Generic") {
          f<<"\\verb;"<<mes[i].second->Name()<<"; & \\\\"<<endl;
          total_mes++;
        }
      }
      currents=dc->GetCurrents();
      for(size_t i=0;i<currents.size();i++) {
        f<<"\\verb;"<<currents[i].second.first->Name()
         <<"; $\\otimes$ \\verb;"<<currents[i].second.second->Name()<<"; & \\\\"<<endl;
        total_mes++;
      }
    }
    // skip inclusives for now
    f<<"\\hline"<<endl;
    f<<"\\end{longtable}"<<endl;
    continue;

    f<<"\\hline"<<endl;
    msg_Out()<<"Creating booklet for "<<dt->Flav().TexName()<<endl;
    map<string,pair<double,double> > brmap;
    GetInclusives(dt,brmap);
    f<<"\\multicolumn{3}{c}{\\hfill}\\\\"<<endl;
    f<<"\\multicolumn{3}{c}{\\bf Inclusive Decays with $BR>10^{-6}$}\\\\"<<endl;
    f<<"\\hline"<<endl;
    f<<"Decay Channel & Branching Ratio \\\\"<<endl;
    f<<"\\hline\n\\hline"<<endl;
    for( BRMapString::iterator brit=brmap.begin(); brit != brmap.end(); ++brit ) {
      double br  = (brit->second.first+brit->second.second)/2.;
      double dbr = (brit->second.first-brit->second.second)/2.;
      if( br>1.e-6 ) {
        f<<"$"<<dc->GetDecaying().TexName()<<"$ $\\to$ $"<<brit->first<<"$ & ";
        sprintf( helpstr, "%.4f", br*100. );
        f<<helpstr;
        if( dbr > 0. ) {
          sprintf( helpstr, "%.4f", dbr*100. );
          f<<" $\\pm$ "<<helpstr;
        }
        f<<" \\% \\\\"<<endl;
      }
    }
    f<<"\\hline"<<endl;
    f<<"\\end{longtable}"<<endl;
  }
   
  // end 
  f<<"\\end{document}"<<endl;
  THROW(normal_exit,"Created HADRONS++ booklet for "+ToString<int>(total_decaychannels)
        +" decaychannels with "+ToString<int>(total_mes)+" MEs."
        +" Run 'latex hadrons.tex' for compilation.");
  abort();
}

std::vector<BRPairFlavourSet> Hadron_Decay_Map::GetInclusives(
    Hadron_Decay_Table * dt,                        // decay table to be considered
    BRMapString & brmap,                            // map with flavourset <-> (upper, lower)
    FlavourSet flset,                               // flavour set
    DoublePair br,                                  // upper, lower br
    bool eoi )                                      // end of iteration
{
  Decay_Channel * dc;
  FlavourSet outs;
  vector<BRPairFlavourSet> new_flset, flvec;
  // for each decay channel
  for(size_t I=0; I<dt->size(); ++I) {
    dc = dt->at(I);                            // decay channel
    outs = dc->GetDecayProducts();                          // FlavourSet of products
    // for each decay product
    new_flset.clear();
    double dbr = (dc->DeltaWidth()>0.)? dc->DeltaWidth()/dt->TotalWidth() : 0.;
    new_flset.push_back(BRPairFlavourSet(
          flset,
          DoublePair(
            br.first*(dc->Width()/dt->TotalWidth()+dbr),
            br.second*(dc->Width()/dt->TotalWidth()-dbr))
          ));
    for (FlSetConstIter fl=outs.begin();fl!=outs.end();++fl) {
      if (this->find((*fl))==this->end() ||
          fl->IsStable()) {                                  // if daughter is stable
        for( size_t i=0; i<new_flset.size(); ++i ) {
          new_flset[i].first.insert(*fl);
        }
      }
      else {                                                 // if daughter has DT 
        vector< vector<BRPairFlavourSet> > dauchans; 
        for( size_t i=0; i<new_flset.size(); ++i ) {
          dauchans.push_back( GetInclusives( (Hadron_Decay_Table*) (*this)[(*fl)][0], brmap, new_flset[i].first, new_flset[i].second, 0 ) );
        }
        new_flset.clear();
        for( size_t i=0; i<dauchans.size(); ++i ) 
          for( size_t j=0; j<dauchans[i].size(); ++j )
            new_flset.push_back( dauchans[i][j] );
      }
    }
    // end of iteration?
    if( eoi ) {
      // add to map
      for( size_t i=0; i<new_flset.size(); ++i ) {
        string channel = string("");
        for (FlSetConstIter fl=new_flset[i].first.begin();fl!=new_flset[i].first.end();++fl) 
          channel += fl->TexName()+string(" ");
        brmap[channel].first += new_flset[i].second.first;          // upper limit
        brmap[channel].second += new_flset[i].second.second;        // lower limit
      }
    }
    else {
      for( size_t i=0; i<new_flset.size(); ++i ) 
        flvec.push_back(new_flset[i]);
    }
  }
  if( flvec.size() ) {
  }
  return flvec;
}
