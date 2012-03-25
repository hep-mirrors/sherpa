// -*- C++ -*-
///////////////////////////////////////////////////////////////////////////////
// File: circulator.h                                                        //
// Description: header file for circulator (circulator class)                //
// This file is part of the SISCone project.                                 //
// For more details, see http://projects.hepforge.org/siscone                //
//                                                                           //
// Copyright (c) 2006 Gavin Salam and Gregory Soyez                          //
//                                                                           //
// This program is free software; you can redistribute it and/or modify      //
// it under the terms of the GNU General Public License as published by      //
// the Free Software Foundation; either version 2 of the License, or         //
// (at your option) any later version.                                       //
//                                                                           //
// This program is distributed in the hope that it will be useful,           //
// but WITHOUT ANY WARRANTY; without even the implied warranty of            //
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             //
// GNU General Public License for more details.                              //
//                                                                           //
// You should have received a copy of the GNU General Public License         //
// along with this program; if not, write to the Free Software               //
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA //
//                                                                           //
// $Revision:: 103                                                          $//
// $Date:: 2007-02-18 17:07:34 +0100 (Sun, 18 Feb 2007)                     $//
///////////////////////////////////////////////////////////////////////////////

#ifndef __CIRCULATOR_H__
#define __CIRCULATOR_H__

namespace siscone{

/// a circulator that is foreseen to take as template member either a
/// pointer or an iterator;
template<class T> class circulator { 

public:
  inline circulator(T here, T begin, T end) : m_here(here), m_begin(begin), m_end(end) {}

  inline circulator(const circulator<T> & other) : m_here(other.m_here), m_begin(other.m_begin), m_end(other.m_end) {}

  /// set just the position without resetting the begin and end elements
  void set_position(const circulator<T> & other) {m_here = other.m_here;}
  void set_position(T pointer) {m_here = pointer;}

  T operator()() {return m_here;}

  inline circulator<T> & operator++() { 
    ++m_here; 
    if (m_here == m_end) m_here = m_begin;
    return *this;
  }

  inline circulator<T> & operator--() { 
    if (m_here == m_begin) m_here = m_end;
    --m_here; 
    return *this;
  }
  
  /// NB: for efficiency, this checks only the here element
  bool operator==(const circulator & other) const {return m_here == other.m_here;}
  /// NB: for efficiency, this checks only the here element
  bool operator!=(const circulator & other) const {return m_here != other.m_here;}

private:
  T  m_here, m_begin, m_end;
};

}

#endif // __CIRCULATOR_H__
