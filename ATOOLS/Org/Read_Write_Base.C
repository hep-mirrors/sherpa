#include "ATOOLS/Org/Read_Write_Base.H"

#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/My_Limits.H"

using namespace ATOOLS;

String_Vector Read_Write_Base::s_commandline;
Buffer_Map Read_Write_Base::s_buffermap;
String_Map Read_Write_Base::s_globaltags;

Read_Write_Base::Read_Write_Base(const unsigned int infiles,
				 const unsigned int outfiles):
  File_IO_Base(infiles,outfiles),
  m_comment(1,defaultcom), m_wordsep(1,defaultwsep), m_linesep(1,defaultlsep),
  m_filecontent(infiles)
{
  Init();
}

Read_Write_Base::Read_Write_Base
(const unsigned int infiles,const unsigned int outfiles,
 const std::string &wordsep,const std::string &linesep,
 const std::string &comment,const std::string &ignore):
  File_IO_Base(infiles,outfiles),
  m_comment(1,comment), m_wordsep(1,wordsep), m_linesep(1,linesep),
  m_filecontent(infiles)
{
  if (ignore!=nullstring) m_ignore.push_back(ignore);
  Init();
}

Read_Write_Base::~Read_Write_Base() 
{
  delete p_interpreter;
}

void Read_Write_Base::Init()
{
  p_interpreter = new Algebra_Interpreter();
  m_blank.push_back(defaultblank);
  m_blank.push_back(defaulttab);
  m_vectortype=vtc::horizontal;
  m_matrixtype=mtc::transposed; 
  m_allownans=false;
  m_addcommandline=true;
  m_useglobaltags=true;
  m_ignorecase=false;
  m_ignoreblanks=false;
  m_exactmatch=true;
  m_interprete=true;
  m_cmode=false;
  m_occurrence=std::string::npos;
  m_escape='\\';
  m_namesplit='|';
}

void Read_Write_Base::SplitInFileName(const size_t &i)
{
  std::string name(InputFile(i));
  size_t sb(name.find(m_namesplit));
  if (sb==std::string::npos) return;
  size_t se(name.find(m_namesplit,sb+1));
  if (se==std::string::npos) return;
  SetInputFile(name.substr(0,sb),i);
  m_filebegin.push_back(name.substr(sb+1,se-sb-1));
  m_fileend.push_back(name.substr(se+1));
  msg_IODebugging()<<METHOD<<"(): Set '"<<m_filebegin.back()
		 <<"'->'"<<m_fileend.back()<<"'.\n"; 
}

bool Read_Write_Base::IsBlank(const char &ch) const
{
  for (Char_Vector::const_iterator bit(m_blank.begin());
       bit!=m_blank.end();++bit) if (ch==*bit) return true;
  return false;
}

size_t Read_Write_Base::IsComment(const std::string &ch) const
{
  for (String_Vector::const_iterator bit(m_comment.begin());
       bit!=m_comment.end();++bit) 
    if (ch.substr(0,bit->length())==*bit) return bit->length();
  return false;
}

size_t Read_Write_Base::IsWordSeparator(const std::string &ch) const
{
  for (String_Vector::const_iterator bit(m_wordsep.begin());
       bit!=m_wordsep.end();++bit) 
    if (ch.substr(0,bit->length())==*bit) return bit->length();
  return false;
}

size_t Read_Write_Base::IsLineSeparator(const std::string &ch) const
{
  for (String_Vector::const_iterator bit(m_linesep.begin());
       bit!=m_linesep.end();++bit) 
    if (ch.substr(0,bit->length())==*bit) return bit->length();
  return false;
}

size_t Read_Write_Base::Find(std::string input,std::string parameter,
			     size_t &length) const
{
  if (m_ignorecase) {
    for (size_t i=0;i<input.length();++i) input[i]=toupper(input[i]);
    for (size_t i=0;i<parameter.length();++i) 
      parameter[i]=toupper(parameter[i]);
  }
  size_t cutinputblanks(0), plength(parameter.length());
  if (m_ignoreblanks) {
    for (size_t i(0);i<input.length();++i) {
      if (IsBlank(input[i])) {
	input.erase(i,1);
	++cutinputblanks;
      }
    }
    for (size_t i(0);i<plength;++i) 
      if (IsBlank(parameter[i])) parameter.erase(i,1);
  }
  length=plength+cutinputblanks;
  size_t pos(input.find(parameter));
  if (pos!=std::string::npos && m_exactmatch) {
    if (pos>0) {
      size_t i(0); 
      for (;i<Blank().size();++i) if (input[pos-1]==Blank()[i]) break;
      if (i==Blank().size()) {
	for (i=0;i<WordSeparator().size();++i) 
	  if (input.rfind(WordSeparator()[i],pos)==pos-1) break;
	if (i==WordSeparator().size()) {
	  for (i=0;i<LineSeparator().size();++i) 
	    if (input.rfind(LineSeparator()[i],pos)==pos-1) break;
	  if (i==LineSeparator().size()) pos=std::string::npos;
	}
      }
    }
    if (pos+plength<input.length()) {
      size_t i(0);
      for (;i<Blank().size();++i) if (input[pos+plength]==Blank()[i]) break;
      if (i==Blank().size()) {
	for (i=0;i<WordSeparator().size();++i) 
	  if (input.find(WordSeparator()[i],pos+plength)==pos+plength) break;
	if (i==WordSeparator().size()) {
	  for (i=0;i<LineSeparator().size();++i) 
	    if (input.find(LineSeparator()[i],pos+plength)==pos+plength) break;
	  if (i==LineSeparator().size()) pos=std::string::npos;
	}
      }
    }
  }
  if (pos==std::string::npos) length=0;
  return pos;
}

char Read_Write_Base::PrevChar(String_Vector &buffer,
			       int &line,int &pos) const
{
  if (pos>0) return buffer[line][--pos];
  while (line>0) 
    if (buffer[--line].length()>0) 
      return buffer[line][pos=buffer[line].length()-1];
  return (char)0;
}

char Read_Write_Base::NextChar(String_Vector &buffer,
			       int &line,int &pos) const
{
  if (++pos<(int)buffer[line].length()) return buffer[line][pos];
  while (++line<(int)buffer.size())
    if (buffer[line].length()>0) return buffer[line][pos=0];
  return (char)0;
}

void Read_Write_Base::InterpreteBuffer(String_Vector &buffer,
				       int &line,int &pos,
				       const int level,const bool keep) 
{
  if (buffer.empty() || line>(int)buffer.size() ||
      (line==(int)buffer.size() && pos>=(int)buffer.back().length())) return;
  int oldline(line), brackets(0);
  char cur, blank(m_blank.size()?m_blank.front():defaultblank);
  bool ifcandidate(true), elsecandidate(false), res(true);
  bool first(true), bracketed(false);
  while ((cur=NextChar(buffer,line,pos))!=(char)0) {
    if (line>oldline) {
      if (brackets==0 && level>0) {
	PrevChar(buffer,line,pos);
	return;
      }
      ifcandidate=true;
      oldline=line;
    }
    std::string rem(buffer[line].substr(pos));
    size_t bl;
    if ((bl=IsLineSeparator(rem))>0) {
      if (!keep) for (size_t i(0);i<bl;++i) buffer[line][pos++]=blank;
      else pos+=bl;
      --pos;
      first=false;
      if (brackets==0 && level>0) return;
      ifcandidate=true;
      continue;
    }
    if (cur=='{') {
      if (brackets==0 && level>0 && first) bracketed=true;
      if ((bracketed && brackets==0) || !keep) buffer[line][pos]=blank;
      ++brackets;
      ifcandidate=true;
      continue;
    }
    if (cur=='}') {
      if (brackets>0) --brackets;
      if ((bracketed &&brackets==0) || !keep) buffer[line][pos]=blank;
      if (bracketed && brackets==0 && level>0) {
	PrevChar(buffer,line,pos);
	return;
      }
      ifcandidate=true;
      continue;
    }
    bool noblank(!IsBlank(cur));
    if (noblank) first=false;
    if (!keep) {
      buffer[line][pos]=blank;
      continue;
    }
    if (ifcandidate && 
	buffer[line].substr(pos,2)=="if") {
      int sline(line), spos(pos++);
      std::string arg;
      while ((cur=NextChar(buffer,line,pos))!=(char)0) 
	if (!IsBlank(cur)) break;
      if (cur!='(') continue;
      while ((cur=NextChar(buffer,line,pos))!=(char)0) {
	if (cur==')') break;
	else arg+=cur;
      }
      if (cur!=')') continue;
      p_interpreter->SetTagReplacer(this);
      res=ToType<int>(p_interpreter->Interprete(arg));
      while (sline<line || spos<=pos) {
	buffer[sline][spos]=blank;
	NextChar(buffer,sline,spos);
      }
      InterpreteBuffer(buffer,line,pos,level+1,res);
      elsecandidate=true;
      continue;
    }
    else if (elsecandidate &&
	     buffer[line].substr(pos,4)=="else") {
      for (size_t i(0);i<4;++i) buffer[line][pos++]=blank;
      InterpreteBuffer(buffer,line,pos,level+1,!res);
      elsecandidate=false;
      continue;
    }
    if ((ifcandidate || elsecandidate) && noblank) {
      ifcandidate=false;
      elsecandidate=false;
    }
  }
}

void Read_Write_Base::InterpreteBuffer(String_Vector &buffer) 
{ 
  int line(0), pos(-1); 
  InterpreteBuffer(buffer,line,pos,0,true); 
}

std::string Read_Write_Base::StripEscapes(const std::string &buffer) const
{
  if (buffer.length()==0) return buffer;
  std::string input=buffer;
  size_t pos, next=0;
  while ((pos=input.find(Escape(),next))!=std::string::npos) {
    input.erase(pos,1);
    if (input.length()>pos && input[pos]==Escape()) next=pos+1;
  }
  return input;
}

std::string Read_Write_Base::ReplaceTags(std::string &expr) const
{ 
  std::string tag=expr;
  bool success=false;
  for (std::map<std::string,std::string>::const_iterator 
	 tit=m_tags.begin();tit!=m_tags.end();++tit) {
    size_t pos=tag.find(tit->first);
    if (pos!=std::string::npos) {
      tag.replace(pos,tit->first.length(),tit->second);
      success=true;
    }
  }
  if (m_useglobaltags)
    for (std::map<std::string,std::string>::const_iterator 
	   tit=s_globaltags.begin();tit!=s_globaltags.end();++tit) {
      size_t pos=tag.find(tit->first);
      if (pos!=std::string::npos) {
	tag.replace(pos,tit->first.length(),tit->second);
	success=true;
      }
    }
  if (success && tag!=expr) return ReplaceTags(tag);
  return tag;
}

std::string &Read_Write_Base::KillBlanks(std::string &buffer) const
{
  if (buffer==nullstring) return buffer;
  bool hit=true;
  while (hit && buffer.length()>0) { 
    hit=false;
    if (IsBlank(buffer[0])) {
      buffer.erase(0,1); 
      hit=true;
      break;
    }
  }
  hit=true;
  while (hit && buffer.length()>0) { 
    if (buffer.length()>1 && buffer[buffer.length()-1]==Escape()) break;
    hit=false;
    if (IsBlank(buffer[buffer.length()-1])) {
      buffer.erase(buffer.length()-1,1);
      hit=true;
      break;
    }
  }
  return buffer;
}

std::string &Read_Write_Base::KillComments(std::string &buffer) const
{
  size_t pos;
  for (unsigned int i=0;i<Comment().size();++i) {
    size_t next=0;
    while ((pos=buffer.find(m_comment[i],next))!=std::string::npos) {
      if (pos>0 && buffer[pos-1]==Escape()) next=pos+m_comment[i].length();
      else buffer=buffer.substr(0,pos);
    }
  }
  return KillBlanks(buffer);
}

std::string &Read_Write_Base::KillIgnored(std::string &buffer) const
{
  size_t pos;
  char blank=Blank().size()?Blank().front():defaultblank;
  for (unsigned int i=0; i<Ignore().size(); ++i) {
    size_t next=0;
    while ((pos=buffer.find(m_ignore[i],next))!=std::string::npos) {
      if (pos>0 && buffer[pos-1]==Escape()) next=pos+m_ignore[i].length();
      else {
	buffer=buffer.substr(0,pos)+blank+
	  buffer.substr(pos+m_ignore[i].length());
      }
    }
  }
  return KillBlanks(buffer);
}

void Read_Write_Base::AddFileContent(std::string line,const unsigned int i)
{
  KillComments(line);
  KillIgnored(line);
  size_t length;
  bool lastword(false);
  for (int j(0);j<(int)line.length();++j) {
    std::string rem(line.substr(j));
    if ((length=IsWordSeparator(rem))>0) {
      if (j>0) {
	if (lastword) m_filecontent[i].back().push_back(line.substr(0,j));
	else m_filecontent[i].push_back(String_Vector(1,line.substr(0,j)));
	lastword=true;
      }
      KillBlanks(line=rem.substr(length));
      j=-1;
    }
    else if ((length=IsLineSeparator(rem))>0) {
      if (j>0) {
	if (lastword) m_filecontent[i].back().push_back(line.substr(0,j));
	else m_filecontent[i].push_back(String_Vector(1,line.substr(0,j)));
      }
      lastword=false;
      KillBlanks(line=rem.substr(length));
      j=-1;
    }
  }
  if (line.length()>0) {
    if (lastword) m_filecontent[i].back().push_back(line);
    else m_filecontent[i].push_back(String_Vector(1,line));
  }
}

bool Read_Write_Base::OpenInFile(const unsigned int i,const int mode)
{  
  if (InputPath(i)+InputFile(i)==nullstring) {
    if (m_addcommandline && CommandLine().size()>0) {
      for (size_t j(0);j<CommandLine().size();++j)
	AddFileContent(CommandLine()[j],i);
      return true;
    }
    return false;
  }
  if (InFileMode(i)==fom::unknown) SetInFileMode(fom::permanent);
  if (!m_filecontent[i].empty()&&mode==0) return true;
  SplitInFileName(i);
  std::string lastline, file(InputPath(i)+InputFile(i));
  String_Vector &buffer(s_buffermap[file]);
  msg_IODebugging()<<METHOD<<"(): ("<<this<<") checks buffer '"
		 <<file<<"' -> ("<<&buffer<<")\n";
  My_In_File &infile(InFile(i));
  if (mode&2 || buffer.empty()) {
    msg_IODebugging()<<METHOD<<"(): ("<<this<<") reads '"<<file
		   <<"', mode = "<<infile.Mode()<<"\n";
    infile.Open();	
    buffer.clear();
    if (*infile) {
      getline(*infile,lastline);
      do {
	if (lastline.length()>0) buffer.push_back(lastline);
	getline(*infile,lastline);
      } while (*infile);
    }
  }
  m_filecontent[i].clear();
  bool checkbegin=(bool)(m_filebegin.size()>0);
  bool checkend=(bool)(m_fileend.size()>0);
  int filebegin=0;
  unsigned int occurrence=0;
  if (!buffer.empty()) {
    String_Vector cbuffer;
    for (size_t ln(0);ln<buffer.size();++ln) {
      lastline=buffer[ln];
      if (lastline.length()==0) continue;
      if (checkbegin) {
	for (size_t length=0,j=0;j<m_filebegin.size();++j) {
	  size_t pos=Find(lastline,m_filebegin[j],length);
	  if (pos!=std::string::npos) {
	    if (filebegin==0) lastline=lastline.substr(pos+length);
	    if (occurrence==m_occurrence ||
		m_occurrence==std::string::npos) ++filebegin;
	    if (filebegin==0) ++occurrence;
	    break;
	  }
	}
	if (filebegin==0) {
	  lastline=nullstring;
	}
	else if (checkend) {
	  for (size_t length=0,j=0;j<m_fileend.size();++j) {
	    size_t pos=Find(lastline,m_fileend[j],length);
	    if (pos!=std::string::npos) {
	      if (occurrence==m_occurrence ||
		  m_occurrence==std::string::npos) --filebegin;
	      if (filebegin==0) {
		lastline=lastline.substr(0,pos);
		++occurrence;
	      }
	      break;
	    }
	  }
	}
	if (lastline.length()>0) cbuffer.push_back(lastline);
      }
      else if (lastline.length()>0) {
	if (lastline.length()>0) cbuffer.push_back(lastline);
      }
    }
    InterpreteBuffer(cbuffer);
    for (size_t j(0);j<cbuffer.size();++j)
      AddFileContent(cbuffer[j],i);
  }
  if (m_addcommandline && CommandLine().size()>0) {
    for (size_t j(0);j<CommandLine().size();++j)
      AddFileContent(CommandLine()[j],i);
  }
  msg_IODebugging()<<METHOD<<"(): Read file content '"<<InputPath()
		 <<InputFile()<<"' {\n";
  for (size_t j(0);j<m_filecontent[i].size();++j) {
    msg_IODebugging()<<"  ";
    for (size_t k(0);k<m_filecontent[i][j].size();++k)
      msg_IODebugging()<<"'"<<m_filecontent[i][j][k]<<"' ";
    msg_IODebugging()<<"\n";
  }
  msg_IODebugging()<<"}\n";
  if (buffer.empty()) infile.SetMode(fom::error);
  return m_filecontent[i].size();
}

void Read_Write_Base::CloseInFile(const unsigned int i,const int mode)
{ 
  msg_IODebugging()<<METHOD<<"(): ("<<this<<") closes file '"
		 <<InputPath(i)+InputFile(i)<<"', mode = "
		 <<InFileMode(i)<<"\n";
  My_In_File &infile(InFile(i));
  if (infile()==NULL) return;
  m_filecontent[i].clear();
  if (infile.Mode()==fom::permanent && !mode) return;
  std::string file(InputPath(i)+InputFile(i));
  if (s_buffermap.find(file)!=s_buffermap.end()) {
    msg_IODebugging()<<METHOD<<"(): ("<<this<<") clears buffer '"
                   <<file<<"' -> ("<<&s_buffermap[file]<<")\n";
    s_buffermap.erase(s_buffermap.find(file));
  }
  infile.Close(); 
}

bool Read_Write_Base::OpenOutFile(const unsigned int i)
{  
  if (OutputFile(i)==nullstring) return false;
  if (OutFileMode(i)==fom::unknown) SetOutFileMode(fom::permanent);
  My_Out_File &outfile(OutFile(i));
  if (outfile()==NULL) {
    outfile.Open();	
    if (m_filebegin.size()>0 && !outfile->bad()) {
      (*outfile)<<m_filebegin[0]<<std::endl;
    }
  }
  return !outfile->bad();
}

void Read_Write_Base::CloseOutFile(const unsigned int i,const int mode)
{ 
  My_Out_File &outfile(OutFile(i));
  if (outfile()==NULL) return;
  if (outfile.Mode()==fom::permanent && !mode) return;
  if (m_fileend.size()>0 && !outfile->bad()) {
    (*outfile)<<m_fileend[0]<<std::endl;
  }
  outfile.Close(); 
}

namespace ATOOLS {

  template <class Type> Type Read_Write_Base::Default() const
  { 
    return std::numeric_limits<Type>::max(); 
  }

  template int Read_Write_Base::Default<int>() const;
  template unsigned int Read_Write_Base::Default<unsigned int>() const;
  template long int Read_Write_Base::Default<long int>() const;
  template unsigned long int Read_Write_Base::Default<unsigned long int>() const;
  template float Read_Write_Base::Default<float>() const;
  template double Read_Write_Base::Default<double>() const;

  template <> std::string Read_Write_Base::Default<std::string>() const
  {
    return "";
  }
	
}

