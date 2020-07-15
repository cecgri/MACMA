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

#include <sstream>
#include "MACMA/physics/earth.h"
#include "MACMA/physics/continent.h"
#include "MACMA/physics/leftContinentExtremity.h"
#include "MACMA/physics/rightContinentExtremity.h"
#include "MACMA/physics/warmingZone.h"
#include "MACMA/physics/plate.h"
#include "MACMA/physics/age.h"
#include "MACMA/physics/geoElement.h"

using namespace std;
// ==================================================
Continent::Continent(Earth* earth, double position, 
					 double length, string origin) :  
_earth(earth),
_position(position),
_length(length),
_origin(origin)
{
	/* Constructor for continents. 
     "origin" is the reason for creation ("initial", 
     "breakup" or "collision")*/
	_id = Earth::getId("Continent");
	
	if(_earth->isWritingContinents())
    {
		string filename = _workspace_ + "elements/" + _id + ".log";
		_file.open(filename.c_str());
		if(_file.is_open())
			_file << "# Age_creation= " << _earth->getAgeMa() 
			<< "\t Length= " << _length
			<< "\t " << _origin
			<< endl;
    }
	
	_rightExtremity = NULL;
	_leftExtremity = NULL;
	
	_plateSection = NULL;
	
	_warmingZone = new WarmingZone(_earth, this);
	_selected = false;
	
	_birthDate = Earth::getTimeMyr();
	_breakable = false;
}
// ========================================
Continent::~Continent()
{
	delete _warmingZone;

	if(_file.is_open())
		_file.close();
}
// ========================================
void
Continent::setPosition(double position)
{
	_position = position;
	_rightExtremity->setPosition(_position);
	_leftExtremity->setPosition(::addPosition(_position, _length));
}
// ========================================
void
Continent::setLength(double length)
{
	_length = length;
	_leftExtremity->setPosition(::addPosition(_position, _length));
}
// ========================================
double
Continent::getLm()
{ // length in m (_length is in degrees)
	return _earth->deg_to_m(_length);
}
// ========================================
double
Continent::getU()
{ // in cm/yr
	return _plateSection->getPlate()->getU();
}


///////////////////////////////////////////////////////
/*
 Positions
 */
double
Continent::getLeftPosition()
{
	return ::addPosition(_position, _length);
}

double
Continent::getRightPosition()
{
	return _position;
}

double
Continent::getMiddlePosition()
{
	return ::addPosition(_position, 0.5 * _length);
}

void
Continent::updatePosition(double dt)
{
	if(!_earth->fixedConfiguration())
    {
		_leftExtremity->updatePosition(dt);
		_rightExtremity->updatePosition(dt);
		
		_position = _rightExtremity->getPosition();
		
		if(_earth->growingContinents())
			_length = ::computeRightLeftDistance(_rightExtremity->getPosition(),
												 _leftExtremity->getPosition());
		
		_warmingZone->updateLifetime(dt);
		_warmingZone->updateWarming();
    }
}
// ==================================================
/* AGES */
void
Continent::shareNeighborAges()
{ /* make sure the right and left continent extremities
   have the same sectionAges as their neighbor (staple
   or subduction). */
	if(_rightExtremity->getNeighbor())
    {
		GeoElement* neighbor = _rightExtremity->getNeighbor();
		_rightExtremity->setRightAges(neighbor->getRightAges());
		_rightExtremity->setRightAge(neighbor->getRightAge());
		_rightExtremity->setLeftAges(neighbor->getLeftAges());
		_rightExtremity->setLeftAge(neighbor->getLeftAge());		
    }
	
	if(_leftExtremity->getNeighbor())
    {
		GeoElement* neighbor = _leftExtremity->getNeighbor();
		_leftExtremity->setRightAges(neighbor->getRightAges());
		_leftExtremity->setRightAge(neighbor->getRightAge());
		_leftExtremity->setLeftAges(neighbor->getLeftAges());
		_leftExtremity->setLeftAge(neighbor->getLeftAge());
    }
}
// --------------------------------------------------
void
Continent::gatherAges(vector<Age> rightAges, vector<Age> leftAges)
{
    vector<Age> newAges;
    
    /* For collisions between two continents:
     average the position and age in between
     the two continents (they should be touching) */
    
    // get the order left/right with checking the smallest distance right-left
    double dRL =
    ::computeRightLeftDistance(rightAges[rightAges.size()-1].getPosition(),
                               leftAges[0].getPosition());
    double dLR =
    ::computeRightLeftDistance(leftAges[0].getPosition(),
                               rightAges[rightAges.size()-1].getPosition());
    
    Age leftAge = Age(_earth, 0.0, 0.0);
    Age rightAge = Age(_earth, 0.0, 0.0);
	
    if (::isStrictlyLess(dLR, dRL))
    { // switch order right / left
        rightAge = leftAges[0];
        leftAge  = rightAges[rightAges.size()-1];
    }
    else
    { //keep normal order
        rightAge = rightAges[rightAges.size()-1];
        leftAge = leftAges[0];
    }
    
    double position = ::getMiddle(rightAge.getPosition(),
                                  leftAge.getPosition());
    double value = _earth->linearAge(position,
                                     rightAge,
                                     leftAge);
    
    Age midAge = Age(_earth, value, position);
    midAge.setOceanic(false);
    
    // push right ages without the points past midAge->getPosition
    unsigned int i = 0;
    int iEnd = -1;
    while(true)
    {
        if( rightAges.size() > 1 &&
           ::isStrictlyWithin(midAge.getPosition(),
                              rightAges[i].getPosition(),
                              rightAges[i+1].getPosition() ) )
        {
            iEnd = i+1;
            break;
        }
        
        i ++;
        if( i > rightAges.size()-1 )
        { // all points are right of midAge
            iEnd = rightAges.size();
            break;
        }
    }
    if(iEnd == -1)
    {
        cerr << "ERROR in gatherAges for continental collision " << endl;
        cerr << " did not find limit position of the right continent" << endl;
        exit(EXIT_FAILURE);
    }
    for(unsigned int i=0; i<iEnd; i++)
        newAges.push_back(rightAges[i]);
    
    // push mid age
    newAges.push_back(midAge);

    // push left ages from first point after midAge->getPosition()
    i = leftAges.size()-1;
    int iStart = -1;
    while(true)
    { // find first i
        if( leftAges.size() > 1 &&
           ::isStrictlyWithin(midAge.getPosition(),
                              leftAges[i-1].getPosition(),
                              leftAges[i].getPosition()) )
        {
            iStart = i;
            break;
        }
        
        i--;
        if(i < 1)
        {
            iStart = 0;
            break;
        }
    }
    if(iStart == -1)
    {
        cerr << "ERROR in gatherAges for continental collision " << endl;
        cerr << " did not find limit position of the left continent" << endl;
        exit(EXIT_FAILURE);
    }
    
    for(unsigned int i=iStart; i<leftAges.size(); i++)
        newAges.push_back(leftAges[i]);
    
    _rightExtremity->setLeftAges(newAges);
	_rightExtremity->setLeftAge(newAges[0]);
	
    _leftExtremity->setRightAges(newAges);
	_leftExtremity->setRightAge(newAges[newAges.size()-1]);
	
    
    if(_rightExtremity->getNeighbor())
	{
		_rightExtremity->getNeighbor()->setLeftAges(newAges);
		_rightExtremity->getNeighbor()->setLeftAge(newAges[0]);
	}
        
    
    if(_leftExtremity->getNeighbor())
	{
		_leftExtremity->getNeighbor()->setRightAges(newAges);
		_leftExtremity->getNeighbor()->setRightAge(newAges[newAges.size()-1]);
	}
    
}
// --------------------------------------------------
void
Continent::setAgeGroups()
{ /* _ageGroups contains the max age in a slice of ages,
   and the right and left positions of a portion of continent
   within the considered age slice. */
	_ageGroups.clear();
	vector<Age> ages = _plateSection->getRightElement()->getLeftAges();
	
	if(ages.size() < 1)
		return;
	
	vector<pair<double, double> > slices = 
    _earth->getContinentalAgeSlices();
	
	unsigned int i = 0;
	while(true)
    { // find the slice
		unsigned int iSlice = 0;
		for(unsigned int j=0; j<slices.size(); j++)
		{
			if( ages[i].getValue() < slices[j].first &&
			   ages[i].getValue() >= slices[j].second )
			{
				iSlice = j;
				break;
			}
		}
		
		// iSlice is the index of the group considered now
		AgeGroup group;
		group.ageMax = slices[iSlice].first + Earth::contAgeGroupStep;
		group.rightPosition = ages[i].getPosition();
		while(true)
		{
			i++;
			
			if(i == ages.size())
				break;
			
			if(ages[i].getValue() >= slices[iSlice].first ||
			   ages[i].getValue() < slices[iSlice].second )
				break;
		}
		
		if(i == ages.size())
			group.leftPosition = ages[i-1].getPosition();
		else
			group.leftPosition = ages[i].getPosition();
		
		_ageGroups.push_back(group);
		
		if(i == ages.size())
			break;
    }
}
// ==================================================
void
Continent::updateBreakable()
{
	_warmingZone->updateF(); 
	_breakable = (_warmingZone->getF() >= _earth->get_F_lim());
}
// ==================================================
bool
Continent::checkOpening()
{
	bool opening = false;
	
	if(_breakable)
    {
		Plate* plate = _plateSection->getPlate();
		if( plate->getRightElement()->isA("Subduction") &&
		   plate->getLeftElement()->isA("Subduction") )
			opening = true;
    }
	
	return opening;
}

// ==================================================
Plate*
Continent::getPlate()
{
	return _plateSection->getPlate();
}
// ==================================================
void
Continent::updateFile()
{
	if(_file.is_open())
		_file << _earth->getAgeMa() 
		<< "\t" << _position 
		<< "\t" << getU() << endl;
}
// ==================================================
unsigned int
Continent::getCounter()
{
	string strbase = "Continent.";
	string sscount = _id.substr(strbase.size(), _id.size()-strbase.size());
	istringstream iss(sscount);
	unsigned int count;
	iss >> count;
	
	return count + 1;  // +1 to avoid having a counter at 0, in order to use directly 1, 2, 3... in outputs
}
// ==================================================
