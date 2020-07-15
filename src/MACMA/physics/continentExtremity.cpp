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
#include "MACMA/physics/continentExtremity.h"

using namespace std;

// ==================================================
// static for continental growth
double ContinentExtremity::referenceGrowthRate = 1.625E-8; // in kg/(m*yr)
double ContinentExtremity::minDepthSubduction = 150000.0; // in m

// ==================================================
ContinentExtremity::ContinentExtremity(Earth* earth, Continent* continent) : GeoElement(earth, 0.0), 
_continent(continent)
{
	_leftSection = NULL;
	_rightSection = NULL;
	_neighbor = NULL;
	
	_growthRate = 0.0;
	_growing = false;
    
    _activeMargin = NULL;
}
// ==================================================
ContinentExtremity::~ContinentExtremity()
{
}
// ==================================================
string
ContinentExtremity::getClass()
{
  return "ContinentExtremity";
}
// ==================================================
bool
ContinentExtremity::isA(string className)
{
  return className == "ContinentExtremity" || className == "GeoElement";
}
// ==================================================
void
ContinentExtremity::updateVelocity()
{ // in cm/yr
  _U = _continent->getU();

  updateGrowthRate();
  if(_growing)
    {
      modifyVelocity(); /* this includes changing the velocity
			   of the neighboring subduction */
    }

}
// ==================================================
void
ContinentExtremity::updateGrowthRate()
{ // in cm/yr
  _growthRate = 0.0;
  _growing = false;
  if(!_earth->growingContinents())
    return;

  if(_neighbor)
    {
      if( (this->isA("RightContinentExtremity") && _neighbor->isA("LeftSubduction")) ||
	  (this->isA("LeftContinentExtremity") && _neighbor->isA("RightSubduction")) )
	{
	  Subduction* subduction = (Subduction*)_neighbor;
	  double depth = subduction->getDepth();
	  if(depth < ContinentExtremity::minDepthSubduction)
	    return;

	  if(subduction->getRightPlate() && subduction->getLeftPlate())
	    {

	      double Vconv = 0.0;

	      // Right border of continent
	      if(this->isA("RightContinentExtremity") && _neighbor->isA("LeftSubduction"))
		{
		  // velocity of the subducting plate:
		  double Vsub = subduction->getRightPlate()->getU(); // in cm/yr
		  // velocity of the upper plate:
		  double Vup = subduction->getLeftPlate()->getU(); // in cm/yr
		  
		  // convergence velocity:
		  Vconv = Vsub - Vup;
		}

	      // Left border of continent
	      if(this->isA("LeftContinentExtremity") && _neighbor->isA("RightSubduction"))
		{
		  // velocity of the subducting plate:
		  double Vsub = subduction->getLeftPlate()->getU(); // in cm/yr
		  // velocity of the upper plate:
		  double Vup = subduction->getRightPlate()->getU(); // in cm/yr

		  // convergence velocity:
		  Vconv = Vup - Vsub;
		}
	   

	      if(Vconv > 0.0)
		{
		  double T = _earth->getT();

		  double rhoCont = _earth->get_rhoCont(); // in kg/m^2
		  double thickness = _earth->get_cont_thickness(); // in m

		  /* referenceGrowthRate is in kg/(m*yr), 
		     contGrowthCoeff is dimensionless,
		     and thickness is in m.
		     Multiplied by 100.0 to be in cm/yr */
		  double coeff = 100.0 * _earth->get_contGrowthCoeff() * 
		    ContinentExtremity::referenceGrowthRate / (rhoCont * thickness); // in cm/yr
		  
		  double T0 = 1477.0 
		    + 320.0 * pow(Vconv, -0.9) 
		    - 270.0 * pow(Vconv, -1.41);
		  
		  
		  double deltaT = max( T - T0, 0.0 );
	      
		  double a = 2.15
		    - 2.57 * pow(Vconv, -0.9) 
		    + 2.93 * pow(Vconv, -1.31);
		  
		  double b = 18.9 
		    + 6.32 * pow(Vconv, -0.5)
		    - 6.46 * pow(Vconv, -1.31);
		  
		  _growthRate = coeff * pow(deltaT, a) * exp(b);
		}
	    }
	  else
	    {
	      cerr << _earth->getTimeMyr() 
		   << "\t Problem in continent::updateGrowth()\n";
	      cerr << "\t" << subduction->getShortId()
		   << " at position " << subduction->getPosition()
		   << " does not have defined left and right plates" << endl;
	    }
	} // if _neighbor is well oriented subduction

      if(_growthRate > 0.0)
	_growing = true;

    } // if(_neighbor)

}
// ==================================================
double
ContinentExtremity::growthRate()
{ // in cm/yr
  updateGrowthRate();

  return _growthRate;
}
// ==================================================
