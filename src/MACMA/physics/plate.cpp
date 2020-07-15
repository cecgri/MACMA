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
//   author: leyaouanq@cervval.com  & cecile.grigne@univ-brest.fr

#include "MACMA/physics/cell.h"
#include "MACMA/physics/plate.h"
#include "MACMA/physics/continent.h"
#include "MACMA/physics/earth.h"
#include "MACMA/physics/force.h"

using namespace std;
// ==================================================
Plate::Plate(Earth* earth, vector<PlateSection*> sections) : _earth(earth), _sections(sections), _selected(false)
{
	_id = Earth::getId("Plate");
	
	_U = 0.0;
	
	_forces = new Force(_earth);
	
	for(unsigned int i=0; i<_sections.size(); i++)
		_sections[i]->setPlate(this);
	
	updateType();
}

Plate::~Plate()
{
	delete _forces;
	
	if(getLeftElement()->isA("Interface"))
		((Interface*)getLeftElement())->setRightPlate(NULL);
	
	if(getRightElement()->isA("Interface"))
		((Interface*)getRightElement())->setLeftPlate(NULL);
	
	for(unsigned int i=0; i<_sections.size(); i++)
		delete _sections[i];
	_sections.clear();
}

/*
 ======================
 Computing
 ======================
 
 ========================================
 Plate::updateForces()
 -------------------------
 - Updates forces right and left 
 + equivalent viscosities
 ========================================
 */
void
Plate::updateForces()
{
	_forces->zeroAll(); /* needed because forces are added
						 onto a plate, from left and right elements */
	
	if(!_earth->noSubduction())
    {
		GeoElement* leftElement = getLeftElement();
		GeoElement* rightElement = getRightElement();
		_forces->addDrivingForces(leftElement, LEFT);
		_forces->addDrivingForces(rightElement, RIGHT);
		_forces->computeViscosities(leftElement, LEFT);
		_forces->computeViscosities(rightElement, RIGHT);
		_forces->setEtaH(computeHorizontalViscosity());
    }
	
}
// --------------------------------------------------
void
Plate::completeForces()
{
	/* in _forces: all is already computed that's needed 
	 to get the velocity of the plate, but total bending, 
	 mantle drag, viscous shear are not yet computed. */
	GeoElement* leftElement = getLeftElement();
	GeoElement* rightElement = getRightElement();
	
	// mantle drag
	double MD = - _forces->getEtaH() * _U;
	_forces->setMantleDrag(MD);
	
	double ULeft = 0.0;
	double URight = 0.0;
	//  bending
	double B = 0.0;
	if(leftElement->isA("Subduction"))
    {
		B -= _forces->getEtaBLeft() * _U;
		if(leftElement->isA("LeftSubduction"))
		{
			ULeft = getLeftPlate()->getU();
			B += _forces->getEtaBLeft() * ULeft;
		}
    }
	if(rightElement->isA("Subduction"))
    {
		B -= _forces->getEtaBRight() * _U;
		if(rightElement->isA("RightSubduction"))
		{
			URight = getRightPlate()->getU();
			B += _forces->getEtaBRight() * URight;
		}
    }
	_forces->setBending(B);
	
	// viscous shear
	double VS = 0.0;
	if(leftElement->isA("Subduction"))
    {
		VS -= _forces->getEtaVLeft() * _U;
		if(leftElement->isA("LeftSubduction"))
		{
			ULeft = getLeftPlate()->getU();
			VS += _forces->getEtaVLeft() * ULeft;
		}
    }
	if(rightElement->isA("Subduction"))
    {
		VS -= _forces->getEtaVRight() * _U;
		if(rightElement->isA("RightSubduction"))
		{
			URight = getRightPlate()->getU();
			VS += _forces->getEtaVRight() * URight;
		}
    }
	_forces->setViscousShear(VS);
}
// --------------------------------------------------
double
Plate::computeHorizontalViscosity()
{
	double etaH = 0.0;
	
	double continentL = 0.0;
	double oceanL = 0.0;
	
	for(unsigned int i=0; i<_sections.size(); i++)
    {
		if(!_sections[i]->isContinental())
			oceanL += 
			::computeRightLeftDistance( _sections[i]->getRightElement()->getPosition(),
									   _sections[i]->getLeftElement()->getPosition() );
		else
			continentL += 
			::computeRightLeftDistance( _sections[i]->getRightElement()->getPosition(),
									   _sections[i]->getLeftElement()->getPosition() );
    }
	/* lengths are in degrees so far */
	Earth::Thickness thick = _earth->getThicknesses();
	
	etaH += _earth->deg_to_m(oceanL) / 
    ( thick.oceanAst / _earth->getEtaAst() + 
	 thick.oceanUM / _earth->getEtaUm() + 
	 thick.LM / _earth->getEtaM() ); 
	
	if(continentL > 0.0)
		etaH += _earth->deg_to_m(continentL) /
		( thick.subcont / _earth->getEtaSubCont() + 
		 thick.contAst / _earth->getEtaAst() +
		 thick.contUM / _earth->getEtaUm() + 
		 thick.LM / _earth->getEtaM() ); 
	
    etaH *= 2.0 * _earth->getCoeffMantleDrag();  // factor two: sigma_v = eta * V /(d/2) = 2*eta*V/d
	
	return etaH;
}
// --------------------------------------------------
bool
Plate::isTensive()
{ // not used anywhere
	// Not strong enough condition for continental breakup 
	// (new ridge can re-close depending on forces)
	bool tensive = false;
	
	Force* leftForces = new Force(_earth);
	Force* rightForces = new Force(_earth);
	
	leftForces->addDrivingForces(getLeftElement(), LEFT);
	rightForces->addDrivingForces(getRightElement(), RIGHT);
	
	double difference = leftForces->getDrivingForces() - 
    rightForces->getDrivingForces();
	
	tensive = ( difference > 0.0 );
	
	delete leftForces;
	delete rightForces;
	
	return tensive;
}
// --------------------------------------------------


/*
 ========================================
 - Computes the plate's length
 ========================================
 */
double
Plate::computeLength()
{ // in degrees
	double leftPos  = getLeftElement()->getPosition();
	double rightPos = getRightElement()->getPosition();
	
	double L = ::computeRightLeftDistance(rightPos, leftPos); // in degrees
	if(L == 0.0)
    { // check if there's only one plate
		vector<Plate*> plates = _earth->accessPlates();
		if(plates.size() == 1)
			L = 360.0;
    }
	
	return L;
}

double
Plate::computeLm()
{ // in m
	return _earth->deg_to_m(computeLength()); 
}

double
Plate::computeContinentalLength()
{ // in degrees
	double L = 0.0;
	
	vector<Continent*> continents = this->getContinents();
	for(unsigned int i=0; i<continents.size(); i++)
		L += continents[i]->getLength();
	
	return L;
}

double
Plate::getMiddlePosition()
{ // in degrees
	double leftPos  = getLeftElement()->getPosition();
	double rightPos = getRightElement()->getPosition();
	
	double middle = ::getMiddle(rightPos, leftPos);
	
	if(leftPos == rightPos && computeLength() > 180.0)
		middle = fmod(leftPos + 180.0, 360.0);//only one plate
	
	return middle;
}

/*
 ======================
 Cells
 ======================
 */
void
Plate::setCells()
{ // CG
	for(unsigned int i=0; i<_cells.size(); i++)
		delete _cells[i];
	_cells.clear();
	
	if((getLeftElement()->isA("Ridge") && getRightElement()->isA("Subduction")) ||
	   (getLeftElement()->isA("Subduction") && getRightElement()->isA("Ridge")))
    { //simple cell
		Cell* cell = new Cell(_earth);
		cell->setLeftElement(getLeftElement());
		cell->setRightElement(getRightElement());
		cell->setLeftPosition(getLeftElement()->getPosition());
		cell->setRightPosition(getRightElement()->getPosition());
		cell->setPlate(this);
		cell->updateOrientation();
		_cells.push_back(cell);
    }
	else
    { 
		if((getLeftElement()->isA("Subduction") && getRightElement()->isA("Subduction")) ||
		   (getLeftElement()->isA("Ridge") && getRightElement()->isA("Ridge")))
		{ // the borders are the same sort of interface 
			bool ridge = getLeftElement()->isA("Ridge") && getRightElement()->isA("Ridge");
			vector<Continent*> continents = getContinents(); /* continents on this plate; 
															  order: position right to left */
			if(continents.size() > 0 && ridge)
			{
				Cell* firstCell = new Cell(_earth);
				firstCell->setRightElement(getRightElement());
				firstCell->setRightPosition(getRightElement()->getPosition());
				firstCell->setLeftElement(continents[0]->getRightExtremity());
				firstCell->setLeftPosition(continents[0]->getRightExtremity()->getPosition());
				firstCell->setPlate(this);
				firstCell->setOrientation(1.0);
				_cells.push_back(firstCell);		      
				
				for(unsigned int i=0; i<continents.size(); i++)
				{
					Cell* rightCell = new Cell(_earth);
					rightCell->setRightElement(continents[i]->getRightExtremity());
					rightCell->setRightPosition(continents[i]->getRightExtremity()->getPosition());
					rightCell->setLeftElement(NULL);
					rightCell->setLeftPosition(continents[i]->getMiddlePosition());
					rightCell->setPlate(this);
					rightCell->setOrientation(-1.0);
					_cells.push_back(rightCell);
					
					Cell* leftCell = new Cell(_earth);
					leftCell->setRightElement(NULL);
					leftCell->setRightPosition(continents[i]->getMiddlePosition());
					leftCell->setLeftElement(continents[i]->getLeftExtremity());
					leftCell->setLeftPosition(continents[i]->getLeftExtremity()->getPosition());
					leftCell->setPlate(this);
					leftCell->setOrientation(1.0);
					_cells.push_back(leftCell);
					
					if(i < continents.size()-1) // cells between two continents.
					{
						double rightPos = continents[i]->getLeftExtremity()->getPosition();
						double leftPos = continents[i+1]->getRightExtremity()->getPosition();
						
						double middle = ::getMiddle(rightPos, leftPos);
						/* check if there's a ridge (inactive) 
						 in this section between rightPos and leftPos */
						vector<GeoElement*> elements = _earth->accessElements();
						for(unsigned int j=0; j<elements.size(); j++)
						{
							if( elements[j]->isA("Ridge") && 
							   elements[j]->isStrictlyWithin(rightPos, leftPos) )
							{
								middle = elements[j]->getPosition();
								break;
							}
						}
						
						Cell* rightCell2 = new Cell(_earth);
						rightCell2->setRightElement(continents[i]->getLeftExtremity());
						rightCell2->setRightPosition(rightPos);
						rightCell2->setLeftElement(NULL);
						rightCell2->setLeftPosition(middle);
						rightCell2->setPlate(this);
						rightCell2->setOrientation(-1.0);
						_cells.push_back(rightCell2);
						
						Cell* leftCell2 = new Cell(_earth);
						leftCell2->setRightElement(NULL);
						leftCell2->setRightPosition(middle);
						leftCell2->setLeftElement(continents[i+1]->getRightExtremity());
						leftCell2->setLeftPosition(leftPos);
						leftCell2->setPlate(this);
						leftCell2->setOrientation(1.0);
						_cells.push_back(leftCell2);		      
					}
				} // contine
				
				Cell* lastCell = new Cell(_earth);
				lastCell->setRightElement(continents[continents.size()-1]->getLeftExtremity());
				lastCell->setRightPosition(continents[continents.size()-1]->getLeftExtremity()->getPosition());
				lastCell->setLeftElement(getLeftElement());
				lastCell->setLeftPosition(getLeftElement()->getPosition());
				lastCell->setPlate(this);
				lastCell->setOrientation(-1.0);
				_cells.push_back(lastCell);		      
			}
			else
			{/* R-R or S-S with no continent: making two cells, 
			  separated at the middle of the plate */
				Cell* rightCell = new Cell(_earth);
				rightCell->setRightElement(getRightElement());
				rightCell->setLeftElement(NULL);
				rightCell->setLeftPosition(getMiddlePosition());
				rightCell->setRightPosition(getRightElement()->getPosition());
				rightCell->setPlate(this);
				rightCell->updateOrientation();
				_cells.push_back(rightCell);
				
				Cell* leftCell = new Cell(_earth);
				leftCell->setLeftElement(getLeftElement());
				leftCell->setRightElement(NULL);
				leftCell->setLeftPosition(getLeftElement()->getPosition());
				leftCell->setRightPosition(getMiddlePosition());
				leftCell->setPlate(this);
				leftCell->updateOrientation();
				_cells.push_back(leftCell);
			}
		}
		else
		{ // staple - staple
			if(!_earth->onlyOneSection()) // more than one section
			{ /* just make cells with downwelling at staples and upwelling 
			   in the middle between two staples  */
				for(unsigned int i=0; i<_sections.size(); i++)
				{
					double posRight = _sections[i]->getRightElement()->getPosition();
					double posLeft = _sections[i]->getLeftElement()->getPosition();
					double middle = ::getMiddle(posRight, posLeft);
					
					Cell* rightCell = new Cell(_earth);
					rightCell->setRightElement(_sections[i]->getRightElement());
					rightCell->setRightPosition(_sections[i]->getRightElement()->getPosition());
					rightCell->setLeftElement(NULL);
					rightCell->setLeftPosition(middle);
					rightCell->setPlate(this);
					rightCell->setOrientation(-1.0);
					_cells.push_back(rightCell);
					
					Cell* leftCell = new Cell(_earth);
					leftCell->setRightElement(NULL);
					leftCell->setRightPosition(middle);
					leftCell->setLeftElement(_sections[i]->getLeftElement());
					leftCell->setLeftPosition(_sections[i]->getLeftElement()->getPosition());
					leftCell->setPlate(this);
					leftCell->setOrientation(1.0);
					_cells.push_back(leftCell);		  
				}
			} 
			else
			{ // only one plate and one section -> make two cells
				double posStaple = getLeftElement()->getPosition();
				double middle = fmod(posStaple + 180.0, 360.0);
				Cell* cell1 = new Cell(_earth);
				cell1->setPlate(this);
				cell1->setOrientation(-1.0);
				cell1->setRightPosition(posStaple);
				cell1->setRightElement(getLeftElement());
				cell1->setLeftPosition(middle);
				
				Cell* cell2 = new Cell(_earth);
				cell2->setPlate(this);
				cell2->setOrientation(1.0);
				cell2->setRightPosition(middle);
				cell2->setLeftElement(getLeftElement());
				cell2->setLeftPosition(posStaple);
				
				_cells.push_back(cell1);
				_cells.push_back(cell2);
			} // !onlyOneSection()
		} // subduction-subduction, ridge-ridge or staple-staple
    } // ridge - subduction
	
}

/*
 ======================
 Sections
 ======================
 */
PlateSection*
Plate::findSection(const char* name)
{
	for(unsigned int i=0; i<_sections.size(); i++)
		if(!strcmp(_sections[i]->getId().c_str(), name))
			return _sections[i];
	cerr << "Error : section " << name << " not found" << endl;
	return NULL;
}

int
Plate::findSection(PlateSection* section)
{
	if(section)
    {
		for(unsigned int i=0; i<_sections.size(); i++)
			if(_sections[i] == section)
				return i;
		cerr << "Error : plate " + section->getId() + " not found" << endl;
    }
	return -1;
}

void
Plate::addSection(PlateSection* section)
{
	_sections.push_back(section);
}


bool
Plate::isContinental()
{ // true if plate is entirely a continent
	bool continental = true;
	for(unsigned int i=0; i<_sections.size(); i++)
    {
		if(!_sections[i]->isContinental())
		{
			continental = false;
			break;
		}
    }
	return continental;
}

bool
Plate::isOceanic()
{ // true if plate is entirely oceanic
	bool oceanic = true;
	for(unsigned int i=0; i<_sections.size(); i++)
    {
		if(_sections[i]->isContinental())
		{
			oceanic = false;
			break;
		}
    }
	return oceanic;
}

void
Plate::updateType()
{ // R: ridge; S: subducting slab; U: upper plate
	
	_subducting = false;
	
	if(getRightElement()->isA("Ridge"))
    {
		if(getLeftElement()->isA("Ridge"))
			_type = "RR";
		else if(getLeftElement()->isA("RightSubduction"))
			_type = "RU";
		else if(getLeftElement()->isA("LeftSubduction"))
		{
			_type = "RS";
			_subducting = true;
		}
    }
	else if(getRightElement()->isA("RightSubduction"))
    {
		_subducting = true;
		if(getLeftElement()->isA("Ridge"))
			_type = "RS";
		else if(getLeftElement()->isA("RightSubduction"))
			_type = "SU";
		else if(getLeftElement()->isA("LeftSubduction"))
			_type = "SS";
    }
	else if(getRightElement()->isA("LeftSubduction"))
    {
		if(getLeftElement()->isA("Ridge"))
			_type = "RU";
		else if(getLeftElement()->isA("RightSubduction"))
			_type = "UU";
		else if(getLeftElement()->isA("LeftSubduction"))
		{
			_type = "SU";
			_subducting = true;
		}
    }
	else if(_earth->noSubduction())
    {
		_type = "STUCK";
		_subducting = false;
    }
}

/*
 ======================
 Interfaces
 ======================
 */

GeoElement*
Plate::getLeftElement()
{
	return (Interface*)_sections[_sections.size()-1]->getLeftElement();
}

GeoElement*
Plate::getRightElement()
{
	return (Interface*)_sections[0]->getRightElement();
}

/*
 ==============================
 Neighbour plates
 ==============================
 */
Plate*
Plate::getLeftPlate()
{
	return getLeftElement()->getLeftSection()->getPlate();
}

Plate*
Plate::getRightPlate()
{
	return getRightElement()->getRightSection()->getPlate();
}

/*
 ==============================
 Continents on this plate
 ==============================
 */
vector<Continent*> 
Plate::getContinents()
{
	vector<Continent*> continents; // continents on this plate
	vector<Continent*> allContinents = _earth->accessContinents();
	for(unsigned int i=0; i<allContinents.size(); i++)
    {
		if(allContinents[i]->getPlate() == this)
			continents.push_back(allContinents[i]);
    }
	
	if(continents.size() > 1)
    {
		/*sort continents by position 
		 (of their right extremity compared to the right limit of the plate) */
		vector<Continent*> backup = continents;
		unsigned int size = continents.size();
		continents.clear();
		
		double rightPos = getRightElement()->getPosition();
		
		while(continents.size() < size)
		{
			double distance = numeric_limits<double>::infinity();
			Continent* continent = NULL;
			vector<Continent*>::iterator it;
			vector<Continent*>::iterator itToDel;
			for(it=backup.begin(); it!=backup.end(); it++)
			{
				if(::computeRightLeftDistance(rightPos,(*it)->getRightPosition()) < distance)
				{
					distance = ::computeRightLeftDistance(rightPos, (*it)->getRightPosition());
					continent = *it;
					itToDel = it;
				}
			}
			continents.push_back(continent);
			backup.erase(itToDel);	 
		}
    }
	
	return continents;
}

