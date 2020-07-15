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

#include "MACMA/physics/earth.h"
#include "MACMA/physics/cell.h"

using namespace std;

// ==================================================
Cell::Cell(Earth* earth) : _earth(earth), _selected(false)
{
  _id = Earth::getId("Cell");
  
  _leftElement = NULL;
  _rightElement = NULL;
  _plate = NULL;
}

Cell::~Cell()
{
  if(_leftElement)
    _leftElement = NULL;
  if(_rightElement)
    _rightElement = NULL;
}

/*
 ============================
 Computing
 ============================
*/

/*
 ========================================
 Cell::updateOrientation()
 -------------------------
 ========================================
*/
void
Cell::updateOrientation()
{
  _orientation = 1.0;
 
  if(_leftElement && (_leftElement->isA("Ridge") || _leftElement->isA("Staple")))
    {
      if((_rightElement && _rightElement->isA("Subduction")) || !_rightElement)
	_orientation = -1.0;
    }
  else if(_leftElement && _leftElement->isA("Subduction"))
    {
      if((_rightElement && _rightElement->isA("Ridge")) || !_rightElement)
	_orientation = 1.0;
    }
  else if(_rightElement && (_rightElement->isA("Ridge") || _rightElement->isA("Staple")))
    {
      if((_leftElement && _leftElement->isA("Subduction")) || !_leftElement)
	_orientation = 1.0;
    }
  else if(_rightElement && _rightElement->isA("Subduction"))
    {
      if((_leftElement && _leftElement->isA("Ridge")) || !_leftElement)
	_orientation = -1.0;
    }
}

/*
 ============================
 Graphics
 ============================
*/
double
Cell::getStartAngle()
{
  return _rightPosition;
}

double
Cell::getStopAngle()
{
  return _leftPosition;
}

double
Cell::getSpeedFactor()
{
  double speedfactor = 4.0;
  return speedfactor;
}
