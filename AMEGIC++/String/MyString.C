/*
#ifndef __GNUC__
#include "AMEGIC++/String/MyString.H"
#include "ATOOLS/Math/MathTools.H"

using namespace AMEGIC;

MyString::MyString(const char* str)
{
  if (!str) {
    len = 0;
    _string = new char[1];
    _string[0] = 0;
  }
  else {
    len = strlen(str);
    _string = new char[len+1];
    strcpy(_string,str);
  }
}

MyString::MyString(const MyString& str)
{
  len = str.len;
  if (!str._string) {
    len = 0;
    _string = new char[1];
    _string[0] = 0;
  }
  else {
    _string = new char[len+1];
    strcpy(_string,str._string);
  }
}

MyString::~MyString()
{
  if (len==0) delete _string;
         else delete[] _string;
}

MyString& MyString::operator=(const MyString& str)
{
  if (this!=&str) {
    if (len>0)  delete[] _string;
    if (len==0) delete   _string;
    len = str.length();
    _string = new char[len+1];
    strcpy(_string,str.c_str());
  }
  return *this;
}

MyString& MyString::operator+=(const MyString& str)
{ 
  if (str._string) {
    if (_string) {
      MyString tmp(*this);
      if (len==0) delete _string;
             else delete[] _string;
      len += str.len;
      _string = new char[len+1];
      strcpy(_string,tmp._string);
      strcpy(_string+tmp.len,str._string);
    }
    else {
      len = str.len;
      _string = new char[len+1];
      strcpy(_string,str._string);
    }
  }
  return *this;
}

MyString& MyString::operator+=(const char* s)
{ 
  if (_string) {
    if (s) {
      MyString tmp(*this);
      if (len==0) delete _string;
             else delete[] _string;
      len += strlen(s);
      _string = new char[len+1];
      strcpy(_string,tmp._string);
      strcpy(_string+tmp.len,s);
    }
  }
  else {
    len = strlen(s);
    _string = new char[len+1];
    strcpy(_string,s);
  }
  return *this;
}

ostream& AMEGIC::operator<<(ostream& s, const MyString& str)
{
  return s<<str.c_str();
}

char* MyString::c_str() const
{
  return _string;
}

long int MyString::length() const
{
  return len;
}

char MyString::operator[](const long int i) const
{
  if (i>len-1) {
    msg_Error()<<"MyString::Out of Range: "<<_string<<":"<<len<<"/"<<i<<endl;

    return _string[0];
  }
  return _string[i];
}

void MyString::erase(const long int& pos,const long int& _len)
{
  MyString s;
  s += substr(0,pos);
  s += substr(pos+_len);
  if (len==0) delete _string;
         else delete[] _string;
  len = s.len;
  _string = new char[len+1];
  strcpy(_string,s._string);
}

MyString MyString::substr(long int pos,long int _len) const 
{
  if (_len>0) {
    char* help;
    long int i;
    help = new char[_len+1];
    for (i=pos;i<pos+_len;i++) help[i-pos] = _string[i];
    help[_len] = '\0';
    MyString s(help);
    delete[] help;
    return s;
  }
  MyString s;
  return s;
}

MyString MyString::substr(long int pos) const
{
  return substr(pos,len-pos);
}

long int MyString::find(const MyString& str)
{
  long int i,j;

  for (i=0;i<len-str.len+1;i++) {
    if (_string[i]==(str.c_str())[0]) {
      long int sw1 = 1;
      for (j=1;j<str.len;j++) {
	if (_string[i+j]!=(str.c_str())[j]) {
	  sw1 = 0;
	  break;
	}
      }
      if (sw1) return i;
    }
  }
  return -1;
}

void MyString::insert(long int pos,const MyString& str)
{
  MyString s;
  s += substr(0,pos);
  s += str;
  s += substr(pos);
  if (len==0) delete _string;
         else delete[] _string;
  len = s.len;
  _string = new char[len+1];
  strcpy(_string,s._string);
}

int MyString::Convert(int)
{
  MyString s(*this);
  int min = 1;
  if (s[0]=='-') {
    min = -1;
    s   = s.substr(1);
  }
  int numb = 0;
  for (int i=0;i<s.length();i++) 
    numb += (s[i]-48)*int(pow(10.,double(s.length()-i-1)));
  
  return min*numb;
}

double MyString::Convert(double)
{
  // ./,
  long int i;
  for (i=0;i<len;i++) { 
    if (_string[i]=='.' || _string[i]==',') break;
  }
  
  int int_tag;

  if (i==len) return double(Convert(int_tag));

  double c;

  if (i!=0) {
    MyString s1 = substr(0,i);  
    c = double(s1.Convert(int_tag));
  }
  else 
    c = 0.;

  if (i==len-1) return c;

  MyString s2 = substr(i+1);
  c += double(s2.Convert(int_tag))*pow(10,-s2.length());

  return c;
}




//#endif
*/
