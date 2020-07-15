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
#include "MACMA/physics/plate.h"
#include "MACMA/physics/ridge.h"
#include "MACMA/physics/enum_struct.h"

using namespace std;
// ==================================================
Ridge::Ridge(Earth* earth, double position) : Interface(earth, position)
{
  _id = Earth::getId("Ridge");
  _active = true;
}

Ridge::~Ridge()
{
  if(_file.is_open())
    _file.close();
}

string
Ridge::getClass()
{
  return "Ridge";
}

bool 
Ridge::isA(string className) 
{  /* ridge is an interface only if it is active
      (useful in fillStructure()) */
  return className == "Ridge" || 
    (className == "Interface" && _active) || 
    className == "GeoElement"; 
}

// --------------------------------------------------
void
Ridge::updateVelocity()
{ // in cm/yr
  if(_active)
    _U = 0.5 * (_leftPlate->getU()  + _rightPlate->getU());
  else
    { /* could be right- of leftSection. 
	 Ridge is in the middle of a plate */
      _U = _rightSection->getPlate()->getU();
    }
}
// --------------------------------------------------
double 
Ridge::computeRidgePush(Direction direction)
{
    double ridgePush = 0.0;
    
    double T_m = _earth->getT();
    double T_surf = _earth->get_T_surf(); //surface temperature
    
    double rho = _earth->get_rho_pl_p();
    if(_earth->ridgePush_rho_depends_on_T())
        rho = _earth->get_rho_pl();
    
    // maximum age of the section (effectiveAge takes SSC into account)
    double age = 0.0; // age in Myr
    switch(direction)
    {
        case LEFT :
            age = _leftSection->getLeftElement()->getRightAge().getEffValue();
            break;
        case RIGHT :
            age = _rightSection->getRightElement()->getLeftAge().getEffValue();
            break;
        default :
            return ridgePush;
    }
    
    ridgePush = _earth->get_preFactor_ridgePush() *
    rho * (T_m - T_surf) * ::Myr_to_sec(age);
    
    return ridgePush;
}
// --------------------------------------------------
