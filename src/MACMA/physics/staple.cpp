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
#include "MACMA/physics/staple.h"

using namespace std;
// ==================================================
Staple::Staple(Earth* earth, double position) : GeoElement(earth, position)
{
	_id = Earth::getId("Staple");
	
	_leftSection = NULL;
	_rightSection = NULL;
	
	_leftNeighborExtremity = NULL;
	_rightNeighborExtremity = NULL;
}

Staple::~Staple()
{
}

bool
Staple::isA(string className)
{
	return className == "Staple" || className == "GeoElement";
}

string
Staple::getClass()
{
	return "Staple";
}

void
Staple::updateVelocity()
{ // in cm/yr 
	_U = _leftSection->getPlate()->getU(); 
	// _leftSection or _rightSection are on the same plate for a staple
}

void
Staple::updatePosition(double dt)
{ // to have uniform writing for all geoElements
	if(!_earth->fixedConfiguration())
		_position = getNextPosition(dt);
}


double 
Staple::getNextPosition(double dt)
{ // updateVelocity() must have been done before
	return ::addPosition(_position, _earth->cm_to_deg(_U) * dt);
}
