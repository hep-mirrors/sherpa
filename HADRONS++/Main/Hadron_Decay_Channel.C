#include "HADRONS++/Main/Hadron_Decay_Channel.H"
#include "HADRONS++/ME_Library/HD_ME_Base.H"
#include "HADRONS++/ME_Library/Generic.H"
#include "HADRONS++/ME_Library/Current_ME.H"
#include "HADRONS++/Current_Library/Current_Base.H"
#include "HADRONS++/PS_Library/HD_PS_Base.H"
#include "PHASIC++/Decays/Decay_Table.H"
#include "PHASIC++/Channels/Multi_Channel.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Math/Vector.H"
#include "ATOOLS/Org/Shell_Tools.H"
#include "ATOOLS/Phys/Blob.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/My_MPI.H"

using namespace HADRONS;
using namespace ATOOLS;
using namespace METOOLS;
using namespace std;

Hadron_Decay_Channel::Hadron_Decay_Channel(Flavour fl, const Mass_Selector* ms,
                                           string _path) :
  Decay_Channel(fl, ms),
  m_path(_path), m_always_integrate(false),
  m_cp_asymmetry_C(0.0), m_cp_asymmetry_S(0.0)
{
}

Hadron_Decay_Channel::~Hadron_Decay_Channel()
{
}

void Hadron_Decay_Channel::SetFileName(std::string filename)
{
  if(filename=="") {
    filename += GetDecaying().ShellName();
    if (filename=="B_{s}") filename = "Bs";
    filename += string("_");
    for ( int i=0; i<NOut(); i++ ) {
      filename += GetDecayProduct(i).ShellName();
    }
    filename += string(".dat");
  }
  m_filename = filename;
}

vector<vector<string> > Hadron_Decay_Channel::Process(const string& content,
                                                      const string& begin,
                                                      const string& end)
{
  Data_Reader reader(" ",";","!","=");
  reader.AddComment("#");
  reader.AddComment("//");
  reader.SetMatrixType(mtc::transposed);
  reader.AddLineSeparator("\n");
  reader.RescanFileContent(SubString(content,begin,end),true);
  vector<vector<string> > helpsvv;
  reader.MatrixFromString(helpsvv,"");
  return helpsvv;
}

bool Hadron_Decay_Channel::Initialise(GeneralModel startmd)
{
  m_physicalflavours=m_flavours;
  for (size_t i=0; i<m_flavours.size(); ++i) {
    map<kf_code,kf_code>::const_iterator it = 
        startmd.m_aliases.find(m_flavours[i].Kfcode());
    if (it!=startmd.m_aliases.end())
      m_physicalflavours[i] = Flavour(it->second, m_flavours[i].IsAnti());
  }
  
  double totalmass=0.0;
  for (size_t i=1; i<m_flavours.size(); ++i) {
    totalmass+=m_flavours[i].HadMass();
  }
  if(totalmass>m_flavours[0].HadMass()) {
    msg_Error()<<"Error in "<<METHOD<<" for "<<Name()<<"\n"
	       <<"    Total outgoing mass heavier than incoming particle.\n"
	       <<"    Will return and hope for the best.\n";
    return false;
  }
  SetChannels(new PHASIC::Multi_Channel(""));
  Channels()->SetNin(1);
  Channels()->SetNout(NOut());
  m_startmd=startmd;

  // check if dc file exists
  if (FileExists(m_path+m_filename)) {
    msg_Tracking()<<METHOD<<": read "<<m_path<<m_filename<<endl;
    My_In_File infile(m_path+m_filename);
    infile.Open();
    string content = infile.FileContent();
    DEBUG_VAR(content);

    GeneralModel model_for_ps = m_startmd;
    ProcessOptions(content);
    ProcessME(content, model_for_ps);
    ProcessPhasespace(content, model_for_ps);
    ProcessResult(content);
  }
  else { // if DC file does not exist yet
    PRINT_INFO("Decay channel file "<<m_filename<<" in "<<m_path<<" does not exist yet. Will use Isotropic decay.");
    msg_Tracking()<<"No DC file yet in :"<<m_path<<"/"<<m_filename<<".\n";
    int n=NOut()+1;
    vector<int> decayindices(n);
    for(int i=0;i<n;i++) decayindices[i]=i;
    HD_ME_Base* me=new Generic(m_physicalflavours,decayindices,"Generic");
    AddDiagram(me);
    AddPSChannel( string("Isotropic"), 1., m_startmd);
    msg_Tracking()<<"Calculating width for "<<Name()<<":\n";
    CalculateWidth();
    msg_Tracking()<<"   yields "<<m_iwidth<<".\n";
    WriteOut(true,m_path,m_filename);
  }
  return true;
}

void Hadron_Decay_Channel::ProcessOptions(const string& content)
{
  vector<vector<string> > helpsvv = Process(content,"<Options>","</Options>");
  
  for (size_t i=0;i<helpsvv.size();i++) {
    if (helpsvv[i][0]==string("AlwaysIntegrate")) {
      m_always_integrate=atoi(helpsvv[i][1].c_str());
    }
    else if (helpsvv[i][0]==string("CPAsymmetryS")) {
      m_cp_asymmetry_S = ToType<double>(helpsvv[i][1]);
    }
    else if (helpsvv[i][0]==string("CPAsymmetryC")) {
      m_cp_asymmetry_C = ToType<double>(helpsvv[i][1]);
    }
  }
  // convert C and S to lambda, assuming DeltaGamma=0 for the determination 
  // of C and S.
  // this allows DeltaGamma dependent terms in the asymmetry
  double Abs2 = -1.0*(m_cp_asymmetry_C-1.0)/(m_cp_asymmetry_C+1.0);
  double Im = m_cp_asymmetry_S/(m_cp_asymmetry_C+1.0);
  double Re = sqrt(Abs2-sqr(Im));
  m_cp_asymmetry_lambda = Complex(Re, Im);
}

void Hadron_Decay_Channel::ProcessPhasespace(const string& content,
                                             const GeneralModel& model_for_ps)
{
  vector<vector<string> > ps_svv = Process(content,"<Phasespace>","</Phasespace>");
  
  int nr_of_channels=0;
  for (size_t i=0;i<ps_svv.size();i++) {
    double weight = ToType<double>(ps_svv[i][0]);
    if(AddPSChannel( ps_svv[i][1], weight, model_for_ps ) ) nr_of_channels++;
    else {
      msg_Error()<<METHOD<<":  Warning\n"
	  	 <<"   "<<ps_svv[i][1]<<" in "<<m_path<<m_filename
	  	 <<" is not a valid phase space channel.\n"
	  	 <<"   Will ignore it and hope for the best.\n";
    }
  }
  if(nr_of_channels == 0) {
    msg_Error()<<METHOD<<": Warning. No valid phase space channels found in "
	       <<m_path<<m_filename<<". Using Isotropic."<<endl;
    AddPSChannel( string("Isotropic"), 1., m_startmd );
  }
}

void Hadron_Decay_Channel::ProcessME(const string& content,
                                     GeneralModel& model_for_ps )
{
  vector<vector<string> > me_svv = Process(content,"<ME>","</ME>");

  int nr_of_mes=0;
  Algebra_Interpreter ip;
  ip.AddTag("GF", "8.24748e-6");

  for (size_t i=0;i<me_svv.size();i++) {
    if(me_svv[i].size()==3) {
      msg_Tracking()<<"Selecting ME for "<<Name()<<endl;
      HD_ME_Base* me = SelectME( me_svv[i][2] );
      me->SetPath(m_path);
      msg_Tracking()<<"  "<<me->Name()<<endl;
      GeneralModel me_model = m_startmd;
      me_model.AddParameters(SubString(content,"<"+me_svv[i][2]+">","</"+me_svv[i][2]+">"));
      model_for_ps.AddParameters(SubString(content,"<"+me_svv[i][2]+">","</"+me_svv[i][2]+">"));
      me->SetModelParameters( me_model );
      Complex factor = Complex(ToType<double>(ip.Interprete(me_svv[i][0])),
                               ToType<double>(ip.Interprete(me_svv[i][1])));
      me->SetFactor(factor);
      AddDiagram(me);
      nr_of_mes++;
    }
    if(me_svv[i].size()==4) {
      msg_Tracking()<<"Selecting currents for "<<Name()<<endl;
      Current_Base* current1 = SelectCurrent(me_svv[i][2]);
      current1->SetPath(m_path);
      GeneralModel current1_model = m_startmd;
      current1_model.AddParameters(SubString(content,"<"+me_svv[i][2]+">","</"+me_svv[i][2]+">"));
      model_for_ps.AddParameters(SubString(content,"<"+me_svv[i][2]+">","</"+me_svv[i][2]+">"));
      current1->SetModelParameters( current1_model );

      Current_Base* current2 = SelectCurrent(me_svv[i][3]);
      current2->SetPath(m_path);
      GeneralModel current2_model = m_startmd;
      current2_model.AddParameters(SubString(content,"<"+me_svv[i][3]+">","</"+me_svv[i][3]+">"));
      model_for_ps.AddParameters(SubString(content,"<"+me_svv[i][3]+">","</"+me_svv[i][3]+">"));
      current2->SetModelParameters( current2_model );

      msg_Tracking()<<"  "<<current1->Name()<<endl;
      msg_Tracking()<<"  "<<current2->Name()<<endl;

      // Sanity checks for current selection
      if(size_t(1+NOut()) != current1->DecayIndices().size()+
         current2->DecayIndices().size()) {
        msg_Error()<<"Error in "<<METHOD<<": Current selection does not look sane "
                   <<"for "<<Name()<<". Check decaychannelfile."<<std::endl;
        Abort();
      }

      Complex factor = Complex(ToType<double>(ip.Interprete(me_svv[i][0])),
                               ToType<double>(ip.Interprete(me_svv[i][1])));

      vector<int> indices (NOut()+1);
      for(int i=0; i<NOut()+1; i++) indices[i] = i;

      Current_ME* me=
        new Current_ME(m_physicalflavours, indices, "Current_ME");
      me->SetCurrent1(current1);
      me->SetCurrent2(current2);
      me->SetFactor(factor);
      AddDiagram(me);
      nr_of_mes++;
    }
  }
  if(nr_of_mes == 0) {
    msg_Error()<<METHOD<<": Warning. No valid matrix element found in "
               <<m_path<<m_filename<<". Using Generic."<<endl;
    int n=NOut()+1;
    vector<int> decayindices(n);
    for(int i=0;i<n;i++) decayindices[i]=i;
    HD_ME_Base* me=new Generic(m_physicalflavours,decayindices,"Generic");
    AddDiagram(me);
  }
}

void Hadron_Decay_Channel::ProcessResult(const string& content)
{
  vector<vector<string> > result_svv = Process(content,"<Result>","</Result>");
  
  if(result_svv.size()!=1) {
    msg_Info()<<"Calculating width (PR1) for "<<Name()<<endl;
    CalculateWidth();
    msg_Info()<<"   yields "<<m_iwidth<<".\n";
    WriteOut(false,m_path,m_filename);
  }
  else if(result_svv[0].size()==3 && m_always_integrate) {
    msg_Info()<<"Calculating width (PR2) for "<<Name()<<endl;
    CalculateWidth();
    msg_Info()<<"   yields "<<m_iwidth<<".\n";
    // check whether result is different from before and write out if it is
    double oldwidth=ToType<double>(result_svv[0][0]);
    double oldmax=ToType<double>(result_svv[0][2]);
    if(oldwidth!=m_iwidth || oldmax!=m_max) WriteOut(false,m_path,m_filename);
  }
  else if(result_svv[0].size()==3) {
    if (result_svv[0][0].find("nan")!=string::npos) {
      PRINT_INFO("Found nan in "<<Name()<<". Ignoring and continuing.");
      return;
    }
    m_iwidth=ToType<double>(result_svv[0][0]);
    m_ideltawidth=ToType<double>(result_svv[0][1]);
    m_max=ToType<double>(result_svv[0][2]);
  }
  else {
    THROW(fatal_error, "Result section of "+m_path+"/"+m_filename+" did not "+
          "contain three entries. Aborting.");
  }
}

void Hadron_Decay_Channel::WriteOut(bool newfile, string path, string file)
{
  string content, entry;  
  My_Out_File to(path+file);
  to.Open();
  to->precision(4);

  if ( newfile ) {                // if DC file doesn't exist yet
    // write header
    *to<<"# Decay: "<<Name()<<endl;
    *to<<"#        "<<setw(m_flavours[0].IDName().length())<<left<<"0"<<" --> ";
    int i=0;
    for (size_t i=1; i<m_flavours.size(); ++i) {
      *to<<setw(m_flavours[i].IDName().length()+1)<<left<<i;
    }
    *to<<endl<<endl;

    // write out options
    *to<<"<Options>"<<endl;
    *to<<"  AlwaysIntegrate = "<<m_always_integrate
       <<"    # 0...read results and skip integration"<<endl;
    *to<<"                         # 1...don't read results and integrate"<<endl;
    *to<<"</Options>"<<endl<<endl;

    // write out phasespace settings
    *to<<"<Phasespace>"<<endl;
    *to<<"  1.0 Isotropic"<<endl;
    *to<<"</Phasespace>"<<endl<<endl;

    // write out ME settings
    *to<<"<ME>"<<endl;
    *to<<"  1.0 0.0 Generic"<<endl;
    *to<<"</ME>"<<endl<<endl;

    // write out result
    *to<<"<Result>"<<endl;
    *to<<"  "<<m_iwidth<<" "<<m_ideltawidth<<" "<<m_max<<";"<<endl; 
    *to<<"</Result>"<<endl; 

  } 
  // if (read DC file)
  else {                                
    // if DC file exists
    PRINT_INFO("Edited "<<m_path<<" "<<m_filename<<" in Decaydata.zip :");
    cout<<"<Result>"<<endl;
    int oldprec=cout.precision(4);
    cout<<"  "<<m_iwidth<<" "<<m_ideltawidth<<" "<<m_max<<";"<<endl;
    cout.precision(oldprec);
    cout<<"</Result>"<<endl;

    My_In_File infile(path+file);
    infile.Open();
    content = infile.FileContent();

    // change integration results and overwrite result section
    string resultline;
    size_t found_begin = content.find("<Result>");
    size_t found_end = content.find("</Result>");
    resultline="<Result>\n  "+ToString(m_iwidth)+" "+ToString(m_ideltawidth)+" "+ToString(m_max)+";\n</Result>\n";
    if (found_begin!=string::npos) content.replace(found_begin,found_end,resultline);
    else content.append(resultline);
    msg_IODebugging()<<content<<endl;
    *to<<content;
  }
}

bool Hadron_Decay_Channel::SetColorFlow(ATOOLS::Blob* blob)
{
  int n_q(0), n_g(0);
  for(int i=0;i<blob->NOutP();i++) {
    if(blob->OutParticle(i)->Flav().IsQuark())      n_q++;
    else if(blob->OutParticle(i)->Flav().IsGluon()) n_g++;
  }
  if(n_q>0 || n_g>0) {
    blob->SetStatus(blob_status::needs_showers);
    Particle_Vector outparts=blob->GetOutParticles();
    if(m_diagrams.size()>0) {
      // try if the matrix element knows how to set the color flow
      HD_ME_Base* firstme=(HD_ME_Base*) m_diagrams[0];
      bool anti=blob->InParticle(0)->Flav().IsAnti();
      if(firstme->SetColorFlow(outparts,n_q,n_g,anti)) return true;
    }
    // otherwise try some common situations
    int n=outparts.size();
    if(n_q==2 && n_g==0 && n==2) {
      if(outparts[0]->Flav().IsAnti()) {
        outparts[0]->SetFlow(2,-1);
        outparts[1]->SetFlow(1,outparts[0]->GetFlow(2));
      }
      else {
        outparts[0]->SetFlow(1,-1);
        outparts[1]->SetFlow(2,outparts[0]->GetFlow(1));
      }
      return true;
    }
    else if(n_q==0 && n_g==2) {
      int inflow(-1), outflow(-1);
      Particle_Vector::iterator pit;
      for(pit=outparts.begin(); pit!=outparts.end(); pit++) {
        if((*pit)->Flav().IsGluon()) {
          (*pit)->SetFlow(2,inflow);
          (*pit)->SetFlow(1,outflow);
          inflow=(*pit)->GetFlow(1);
          outflow=(*pit)->GetFlow(2);
        }
      }
      return true;
    }
    else if(n_q==0 && n_g==n) {
      outparts[0]->SetFlow(2,-1);
      outparts[0]->SetFlow(1,-1);
      for(int i=1;i<n-1;++i) {
        unsigned int c=Flow::Counter();
        outparts[i]->SetFlow(2,c-1);
        outparts[i]->SetFlow(1,c);
      }
      outparts[n-1]->SetFlow(2,outparts[n-2]->GetFlow(1));
      outparts[n-1]->SetFlow(1,outparts[0]->GetFlow(2));
      return true;
    }
    else {
      msg_Error()<<METHOD<<" wasn't able to set the color flow for"<<endl
                 <<*blob<<endl;
      return false;
    }
  }
  else return true;
}

Current_Base* Hadron_Decay_Channel::SelectCurrent(string current_string)
{
  Data_Reader reader(",",";","#","]");
  reader.AddWordSeparator("[");
  vector<string> resultstrings;
  reader.SetString(current_string);
  reader.VectorFromString(resultstrings);
  
  int n=resultstrings.size()-1;
  vector<int> indices(n);
  for(int i=0; i<n; i++) indices[i] = ToType<int>(resultstrings[i+1]);
  ME_Parameters fi(m_physicalflavours, indices);

  Current_Base* current=Current_Getter_Function::GetObject(resultstrings[0],fi);
  if(current==NULL) {
    msg_Error()<<METHOD<<": Current '"<<resultstrings[0]<<"' specified in "
               <<m_path<<m_filename<<" was not recognized as a valid current. "
               <<"Will abort."<<endl;
    Abort();
  }
  return current;
}


HD_ME_Base * Hadron_Decay_Channel::SelectME(string me_string)
{
  Data_Reader reader(",",";","#","]");
  reader.AddWordSeparator("[");
  vector<string> resultstrings;
  reader.SetString(me_string);
  reader.VectorFromString(resultstrings);
  if(resultstrings.size()==1 && resultstrings[0]=="Generic") {
    for(int i=0;i<NOut()+1;i++) 
      resultstrings.push_back( ToString<size_t>(i) );
  }
  if(int(resultstrings.size())!=NOut()+2) {
    msg_Error()<<METHOD<<" Error: Number of indices in \""<<me_string<<"\" ("
      <<int(resultstrings.size())-1<<") in "<<m_path<<m_filename<<" doesn't "
      <<"equal number of particles ("<<NOut()+1<<"). Will abort."<<endl;
    Abort();
  }

  int n=NOut()+1;
  vector<int> indices(n);
  for(int i=0; i<n; i++) indices[i] = ToType<int>(resultstrings[i+1]);
  ME_Parameters fi(m_physicalflavours, indices);

  HD_ME_Base* me = HD_ME_Getter_Function::GetObject(resultstrings[0],fi);
  if(me==NULL) {
    msg_Error()<<METHOD<<": Error. Matrix element \""<<me_string<<"\" specified in "
      <<m_path<<m_filename<<" was not recognized as a valid ME. Will abort."<<endl;
    Abort();
  }
  return me;
}

void Hadron_Decay_Channel::LatexOutput(std::ostream& f, double totalwidth)
{
  f<<"$"<<GetDecaying().TexName()<<"$ $\\to$ ";
  for (size_t i=1; i<m_flavours.size(); ++i)
    f<<"$"<<m_flavours[i].TexName()<<"$ ";
  f<<" & ";
  char helpstr[100];
  sprintf( helpstr, "%.4f", Width()/totalwidth*100. );
  f<<helpstr;
  if( DeltaWidth() > 0. ) {
    sprintf( helpstr, "%.4f", DeltaWidth()/totalwidth*100. );
    f<<" $\\pm$ "<<helpstr;
  }
  f<<" \\% ";
  if(Origin()!="") {
    f<<"[\\verb;"<<Origin()<<";]";
  }
  f<<"\\\\"<<endl;
  if((m_diagrams.size()>0 &&
      ((HD_ME_Base*) m_diagrams[0])->Name()!="Generic")) {
    sprintf( helpstr, "%.4f", IWidth()/totalwidth*100. );
    f<<" & "<<helpstr;
    if( IDeltaWidth() > 0. ) {
      sprintf( helpstr, "%.4f", IDeltaWidth()/totalwidth*100. );
      f<<" $\\pm$ "<<helpstr;
    }
    f<<" \\% ";
  }
  for(size_t i=0;i<m_diagrams.size();i++) {
    HD_ME_Base* me=(HD_ME_Base*) m_diagrams[i];
    if(me->Name()=="Current_ME") {
      Current_ME* cme=(Current_ME*) me;
      f<<"\\verb;"<<cme->GetCurrent1()->Name()
       <<";$\\otimes$\\verb;"<<cme->GetCurrent2()->Name()<<"; & \\\\"<<endl;
    }
    else if (me->Name()=="Generic") {
      // do nothing
    }
    else {
      f<<"\\verb;"<<me->Name()<<"; & \\\\"<<endl;
    }
  }
}

bool Hadron_Decay_Channel::AddPSChannel(string name,double weight,
                                        GeneralModel const & md)
{
  PHASIC::Single_Channel * sc=
    HD_Channel_Selector::GetChannel(1, NOut(),&m_flavours.front(),name,md,p_ms);
  if (sc!=NULL) {
    sc->SetAlpha(weight);
    Channels()->Add(sc);
    return true;
  }
  else return false;
}

