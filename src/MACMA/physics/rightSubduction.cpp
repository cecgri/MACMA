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
#include "MACMA/physics/rightSubduction.h"

using namespace std;
// ==================================================
RightSubduction::RightSubduction(Earth* earth, double position) : Subduction(earth, position)
{
  _id = Earth::getId("RightSubduction");
}

RightSubduction::~RightSubduction()
{
  if(_file.is_open())
    _file.close();
}

string
RightSubduction::getClass()
{
  return "RightSubduction";
}

void
RightSubduction::updateVelocity()
{
  /* trench for right subduction has the velocity 
     of the right(upper) plate */
  _U = _rightPlate->getU(); // in cm/yr
}

void
RightSubduction::updatePosition(double dt)
{
  if(!_earth->fixedConfiguration())
    _position = getNextPosition(dt);

  double Vup  = _rightPlate->getU(); // in cm/yr
  double Vsub = _leftPlate->getU();

  double deep = 0.01 * (Vsub - Vup) * dt; /* dt is in yrs, V in cm/yr
					     and depth in m */
  if(_depth < 2800.0e3)
    _depth -= deep;

  updateThickness(-deep);
}
