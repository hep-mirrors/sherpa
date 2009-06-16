#ifndef AMEGIC_String_MyString_H
#define AMEGIC_String_MyString_H
//String header


//#ifdef __GNUC__
#include <string> 
/*
#else
#include <iostream>
#include <string>

namespace AMEGIC {
  class MyString {
  public:
    MyString(const char* = 0);
    MyString(const MyString&);
    ~MyString();
    MyString& operator=(const MyString&);
    void erase(const long int& pos,const long int& length);
    MyString& operator+=(const MyString& str);
    MyString& operator+=(const char*);
    char operator[](const long int) const;
    char* c_str() const;
    long int length() const;
    MyString substr(long int) const;
    MyString substr(long int,long int) const;
    long int find(const MyString&);
    void insert(long int,const MyString&);
    char* _string;
    double Convert(double);
    int    Convert(int);
  private:
    long int   len;  
  };
  
  inline MyString operator+(const MyString& str1,const MyString& str2)
  {
    MyString s(str1);
    s += str2;
    return s;
  }
  
  inline bool operator==(const MyString& str1,const MyString& str2) 
  {
    if (str1.length()!=str2.length()) return false;
    return strcmp(str1.c_str(),str2.c_str()) ? false : true;
  }
  
  inline bool operator!=(const MyString& str1,const MyString& str2)
  {
    return !(str1==str2);
  }
  
  ostream& operator<<(ostream& s, const MyString&);
  
  typedef MyString string;
}
#endif
*/
#endif

