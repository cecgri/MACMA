// Copyright (C) 2011 ENIB
//   Ecole Nationale d'Ingenieurs de Brest (ENIB)
//   CS 73862 - 29238 BREST Cedex 3 - France
//   Tel: +33(0)298 05 89 89, Fax: +33(0)298 05 89 79, e-mail: combes@enib.fr
//
// This file is part of MACMA.
//
//   MACMA is free software: you can redistribute it and/or modify
//   it under the terms of the GNU General Public License as published by
//   the Free Software Foundation, either version 3 of the License, or
//   (at your option) any later version.
//
//   MACMA is distributed in the hope that it will be useful,
//   but WITHOUT ANY WARRANTY; without even the implied warranty of
//   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//   GNU General Public License for more details.
//
//   You should have received a copy of the GNU General Public License
//   along with MACMA.  If not, see <http://www.gnu.org/licenses/>.
//
//   author: leyaouanq@cervval.com 

/* august 2015: creation of two inherited classes 
   from Macma, with and without GUI.
   Macma is the parent class.
   author: cecile.grigne@univ-brest.fr (CG) */

#ifndef MACMA_H_
#define MACMA_H_

#include <iostream>
#include <string>
#include <map>

// --- BOOST
#include "boost/thread/mutex.hpp"
#include "boost/thread/thread.hpp"

// --- MACMA
#include "MACMA/physics/earth.h"
#include "MACMA/utils/utils.h"


// --- XML
#include <locale.h>
#include "MACMA/utils/RapidXML/rapidxml.hpp"
#include "MACMA/utils/RapidXML/rapidxml_print.hpp"
using namespace rapidxml;

using std::string;
using std::map;

// ==================================================
class Macma 
{
 public:
  Macma();
  ~Macma();
  
  enum MACMAState
  {
    EMPTY,
    PARAMETERS,
    INTERFACES,
    PLATES,
    READY,
    RUNNING
  };
  
  virtual bool isA(string) = 0; 

  // State
  void setMACMAState(MACMAState);
  inline MACMAState getMACMAState() { return _state; }
  string showMACMAState();

 protected:
  virtual void openXmlFile(string);
  virtual void saveXmlFile(string);
  virtual void writeXmlFile(string);
  Earth* _earth;

  boost::mutex _mutex;
  
  virtual void _earthThread() = 0; // pure virtual

  double _nextTimeWriteDiagnosis;
  double _nextTimeWriteXml;

  boost::thread _thread;

  MACMAState _state;
};
// ==================================================

#endif

