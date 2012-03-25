// -*- C++ -*-
///////////////////////////////////////////////////////////////////////////////
// File: siscone_error.h                                                     //
// Description: header file for SISCone error messages (Csiscone_error)      //
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
// $Revision:: 123                                                          $//
// $Date:: 2007-03-01 02:52:16 +0100 (Thu, 01 Mar 2007)                     $//
///////////////////////////////////////////////////////////////////////////////

#ifndef __SISCONE_ERROR_H__
#define __SISCONE_ERROR_H__

#include<iostream>
#include<string>

namespace siscone{

/// class corresponding to errors that will be thrown by siscone
class Csiscone_error {
public:
  // constructors
  Csiscone_error() {;};
  Csiscone_error(const std::string & message) {
    m_message = message; 
    if (m_print_errors) std::cerr << "siscone::Csiscone_error: "<<message << std::endl;
  };

  std::string message() const {return m_message;};

  static void setm_print_errors(bool print_errors) {
    m_print_errors = print_errors;};

private:
  std::string m_message;
  static bool m_print_errors;
};

}
#endif
