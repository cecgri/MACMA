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
#include "MACMA/physics/interface.h"
#include "MACMA/physics/force.h"
#include "MACMA/physics/enum_struct.h"

using namespace std;
// ==================================================
Interface::Interface(Earth* earth, double position) : GeoElement(earth, position)
{
  _leftPlate = NULL;
  _rightPlate = NULL;
}

Interface::~Interface()
{
  
}

bool
Interface::isA(string className)
{
  return className == "Interface" || className == "GeoElement";
}

string
Interface::getClass()
{
  return "GeoElement";
}

// -------------------------------------------------------
Direction
Interface::getDirectionDrivingForces(Direction direction)
{ 
  /* direction is for the plate that is on the "direction" 
   side of this interface. */
  Direction directionDrivingForces;

  Force* forces = new Force(_earth);
  Interface* nextInterface = NULL;

  switch(direction)
    {
    case LEFT :
      {
	nextInterface = findNextInterface(LEFT);
	if(nextInterface)
	  {	    
	    forces->addDrivingForces(this, RIGHT);
	    forces->addDrivingForces(nextInterface, LEFT);
	  }
	break;
      }
    case RIGHT :
      {
	nextInterface = findNextInterface(RIGHT);
	if(nextInterface)
	  {
	    forces->addDrivingForces(nextInterface, RIGHT);
	    forces->addDrivingForces(this, LEFT);
	  }
	break;
      }
    default :
      {
	cerr << "ERROR: Interface::getDirectionDrivingForces(direction)\n is used without given direction \n";
	exit(EXIT_FAILURE);
      }
    }

  double limit = 1.0E-6;
  if(nextInterface)
    directionDrivingForces =
      ::doubleToDirection(forces->getDrivingForces(), limit);
  else
    // couldn't find any other interface: this is the only interface on _earth
    {
      if (this->isA("LeftSubduction"))
	directionDrivingForces = LEFT;
      else if (this->isA("RightSubduction"))
	directionDrivingForces = RIGHT;
    }
	
	delete forces;

  return directionDrivingForces;
}
// -------------------------------------------------------
Interface*
Interface::findNextInterface(Direction direction)
{ /* find the closest interface on the left or the right.
     Useful for interfaces that are not yet included in earth _elements.
   */
  Interface* nextInterface = NULL;
  
  double minDistance = 360.0;
  vector<GeoElement*> elements = _earth->accessElements();

  for(unsigned int i=0; i<elements.size(); i++)
    {
      if(elements[i]->isA("Interface") && elements[i] != this)
	{
	  double distance = 360.0;
	  switch(direction)
	    {
	    case LEFT :
	      {	      
		distance =
		  ::computeRightLeftDistance(_position,
					     elements[i]->getPosition());
		break;
	      }
	    case RIGHT :
	      {
		distance =
		  ::computeRightLeftDistance(elements[i]->getPosition(),
					     _position);
		break;
	      }
	    default :
	      { 
		cerr << "ERROR: Interface::findNextInterface(direction) is used without direction given" << endl;
		exit(EXIT_FAILURE);
	      }
	    } // switch

	  minDistance = min(distance, minDistance);
	  if(distance == minDistance)
	    nextInterface = (Interface*)elements[i]; 
	} // elements[i]
    }
    
  return nextInterface;
}
// -------------------------------------------------------
bool
Interface::nearInterface(double limitDistance)
{
  bool tooClose = false;

  Interface* leftInterface = findNextInterface(LEFT);
  if(leftInterface)
      tooClose = 
	::computeRightLeftDistance(_position, 
				   leftInterface->getPosition()) < limitDistance;

  if(!tooClose)
    {
      Interface* rightInterface = findNextInterface(RIGHT);
      if(rightInterface)
	tooClose = 
	  ::computeRightLeftDistance(rightInterface->getPosition(), 
				     _position) < limitDistance;
    }

  return tooClose;
}
// -----------------------------------------------------
