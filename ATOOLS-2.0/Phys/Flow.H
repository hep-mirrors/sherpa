#ifndef Flow_h
#define Flow_h

#include <iostream>
#include <map>
#include <set>

namespace APHYTOOLS {
  class Parton;
  class Flow {  
    friend std::ostream& operator<<( std::ostream& ostr, const Flow& f );
  private:
    Parton *            m_owner;
    std::map<int,int>   m_code;
    static long qcd_counter;
  public:
    Flow( Parton * parton_owner = 0 );
    Flow( const Flow & );
    virtual ~Flow();
    Flow &  operator=( const Flow& );
    bool    operator==( const Flow& a ) const; //compares only flow
    bool    operator!=( const Flow& a ) const; //patterns not owner

    const Parton *     Owner() const;
    int                Code( int _index = 1 ) const;
    void               SetCode( int code_index, int code );
    void               SetUniqueCode( int code_index = 1 );


    //////////////////////
    // container access //
    //////////////////////
    bool               empty() const;
    int                size() const;
    void               clear();
    bool               erase( int code_index );
    typedef std::map<int,int>::iterator       iterator;
    typedef std::map<int,int>::const_iterator const_iterator;
    iterator            begin();
    iterator            end();
    const_iterator      begin() const;
    const_iterator      end() const;

    static void ResetCounter()       { qcd_counter = 600; }
  };  

  inline const Parton * Flow::Owner() const { return m_owner; }
  inline int Flow::Code( int _index ) const {
    std::map<int,int>::const_iterator a = m_code.find(_index);
    return a==m_code.end() ? 0 : (*a).second;
  }
  inline void Flow::SetCode( int _index, int _code )  { m_code[_index] = _code; }
  inline void Flow::SetUniqueCode( int _index )       { m_code[_index] = ++qcd_counter; }

  inline bool Flow::empty() const                     { return (bool)m_code.empty(); }
  inline int Flow::size() const                       { return (int)m_code.size(); }
  inline void Flow::clear()                           { m_code.clear(); }
  inline bool Flow::erase( int _index )               { return (bool)m_code.erase( _index ); }
  inline Flow::iterator Flow::begin()                 { return m_code.begin(); }
  inline Flow::iterator Flow::end()                   { return m_code.end(); }
  inline Flow::const_iterator Flow::begin() const     { return m_code.begin(); }
  inline Flow::const_iterator Flow::end() const       { return m_code.end(); }

  inline bool Flow::operator==( const Flow& a ) const { return (m_code == a.m_code); }
  inline bool Flow::operator!=( const Flow& a ) const { return !( *this == a ); }
  inline Flow& Flow::operator=( const Flow& inflow ) {
    m_code = inflow.m_code;
    return *this;
  }
}

#endif
