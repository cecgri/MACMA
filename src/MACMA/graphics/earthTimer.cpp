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
//   author: cecile.grigne@univ-brest.fr

#include "MACMA/graphics/earthTimer.h"
#include <iostream>
#include <sstream>

using std::endl;
using std::ostringstream;

// static (counter for number of timers created)
unsigned int EarthTimer::_timerCounter = 0;



// constructor - - - - - - 
EarthTimer::EarthTimer(Earth* earth) : _earth(earth)
{
  ostringstream oss;
  oss << "earthTimer." << _timerCounter;
  _id = oss.str();

  _timerCounter ++;
  _nextTimeWrite = 0.0;
}

EarthTimer::~EarthTimer()
{
}

// signal emit - - - - - - - 
void
EarthTimer::isTimeReached(double dWrite)
{
  if(_earth->getTimeMyr() >= _nextTimeWrite)
    {      
     _nextTimeWrite += dWrite;
      emit timeReached();
    }
}

