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

#include "MACMA/physics/plateSection.h"
#include "MACMA/physics/earth.h"
#include "MACMA/physics/plate.h"
#include "MACMA/physics/continent.h"
#include "MACMA/physics/geoElement.h" 

using namespace std;
// ==================================================
PlateSection::PlateSection(Earth* earth,
                           GeoElement* rightElement,
                           GeoElement* leftElement) :
_earth(earth),
_rightElement(rightElement),
_leftElement(leftElement)
{
    _id = Earth::getId("PlateSection");
    
    _leftElement->setRightSection(this);
    _rightElement->setLeftSection(this);
    
    _selected = false;
    _length = computeLength();
    
	_plate = NULL;
    _continent = NULL;
}

PlateSection::~PlateSection()
{
    _leftElement->setRightSection(NULL);
    _rightElement->setLeftSection(NULL);
}

/*
 ========================================
 - Computes the section's length
 ========================================
 */
double
PlateSection::computeLength()
{ // in degrees 
    double L = 0.0;
    
    if(!_leftElement || !_rightElement)
    {
        cerr << "PlateSection::computeLength(): \
        Error, a plate section must have left and right elements defined" << endl;
    }
    else
    {
        if(_leftElement == _rightElement)
            L = 360.0;
        else
        {
            L = ::computeRightLeftDistance(_rightElement->getPosition(),
                                           _leftElement->getPosition()); // in deg.
        }
    }
    
    return L;
}
// --------------------------------------------------
double
PlateSection::computeLm()
{ // in m
    return _earth->deg_to_m(computeLength());
}
// --------------------------------------------------
void
PlateSection::initAges()
{ /* initiates ages (at beginning of run and for newly created sections,
   after continental breakup or new ridges for instance),
   and computes effective ages, useful
   to compute plate thicknesses and heat flux */
    
    vector<Age>& ages = _rightElement->getLeftAges();
    
    if(!ages.empty())
        cerr << "Warning: leftAges for "
        << _rightElement->getShortId()
        << " is not initially empty" << endl;
    
    // first element is the _rightElement
    ages.push_back(_rightElement->getLeftAge());
    
    while(true)
    {
		double position = ::addPosition(ages[ages.size()-1].getPosition(),
										_earth->get_resolution());
        if( _leftElement->isWithinOrEqual(ages[ages.size()-1].getPosition(),
                                          position,
                                          LEFT) )
        { // stops when the leftElement has been reached
            break;
        }
        
        double value =  _earth->linearAge( position,
                                          _rightElement->getLeftAge(),
                                          _leftElement->getRightAge() );
        Age newAge = Age(_earth, value, position);
        ages.push_back(newAge);
        
    }
    ages.push_back(_leftElement->getRightAge());
	
    /* continents */
    if(isContinental())
    {
        for(unsigned int i=0; i<ages.size(); i++)
		{
			ages[i].setOceanic(false);
			if (_continent->getOrigin().compare("initial") == 0) 
			{
				ages[i].setOriginContCounter(_continent->getCounter());
			}
		}
    }
    
    /*compute ages on borders and effective
     ages for heat flux and thicknesses. */
    effectiveAges();
    
    /* give exact same age to borders */
    _rightElement->setLeftAge(ages[0]);
    
    /* give the same ages to the leftelement: */
    _leftElement->setRightAges(ages);
    _leftElement->setRightAge(ages[ages.size()-1]);
    
}
// --------------------------------------------------
bool
PlateSection::isContinental()
{
    if(_continent)
        return true;
    
    return false;
}
// --------------------------------------------------
void
PlateSection::updateQtot()
{ // total heat flow of a plate section
    _Qtot = 0.0;
    
    vector<Age> ages = _rightElement->getLeftAges();

    if(!isContinental())
    { // ocean
        double deltaT = _earth->getT() - _earth->get_T_surf();
        double coeffHeatFlow = _earth->get_preFactor_heatFlow() * deltaT;
        
        for(unsigned int i=0; i<ages.size()-1; i++)
        {
            double surf =
            ::computeRightLeftDistance(ages[i].getPosition(),
                                       ages[i+1].getPosition());
            surf = _earth->deg_to_squareM(surf);
            
            double age1 = ::Myr_to_sec(ages[i].getEffValue());
            if(!::areEqual(ages[i].getEffValue(),
                           ages[i+1].getEffValue(),
                           1.0E-4)) // 1.0E-4Myr = 100yr
            { // average between age[i] and age[i+1] (trapeze)
                double age2 = ::Myr_to_sec(ages[i+1].getEffValue());
                
                // mean heat flux multiplied by surface
                _Qtot += 2.0 * coeffHeatFlow *
                ( sqrt(age2) - sqrt(age1) ) /
                ( age2 - age1 ) * surf;
            }
            else
            { // simple
                _Qtot += coeffHeatFlow / sqrt(max(age1, 3.0E10)) * surf;
                /* little min (~1000yr = 3E10s) to avoid NaN */
            }
        } // ages[i]
    }
    else
    { // continent
        double surf = 2.0 * computeLm() * _earth->get_radius(); /* surf=2*alpha*R^2=2*L*R
                                                                 (L=alpha*R) */
        _Qtot += _earth->computeHeatFlux(-1.0) * surf; // any negative age for continent
    }
}
// --------------------------------------------------
void
PlateSection::updateBathymetry()
{
    ////////////////////////////////////////
    // TODO: recheck this
    ////////////////////////////////////////
    
    _VbelowRidge_HSCM = 0.0; //Halfspace cooling model
    _VbelowRidge_PM95 = 0.0; //Plate model for zm=95km
    _VbelowRidge_PM125 = 0.0; //Plate model for zm=125km
    
    if(!isContinental())
    { // ocean
        
        vector<Age> ages = _rightElement->getLeftAges();
        
        /* Constants for Halfspace cooling model HSCM */
        double deltaT = _earth->getT() - _earth->get_T_surf();
        double coeffBathy =
        _earth->get_preFactor_bathyHSCM() * deltaT *
        _earth->get_rho_um() /
        (_earth->get_rho_um() - _earth->get_rho_seawater());
        
        for(unsigned int i=0; i<ages.size()-1; i++)
        {
            double surf =
            ::computeRightLeftDistance(ages[i].getPosition(),
                                       ages[i+1].getPosition());
            
            surf = _earth->deg_to_squareM(surf);
            
            double age1 = ::Myr_to_sec(ages[i].getEffValue());
            if(!::areEqual(ages[i].getEffValue(),
                           ages[i+1].getEffValue(),
                           1.0E-4)) // 1.0E-4Myr = 100yr
            {
                /* HSCM: integration between age1 and age2 to
                 get mean value */
                double age2 = ::Myr_to_sec(ages[i+1].getEffValue());
                
                // mean depth multiplied by surface
                _VbelowRidge_HSCM += 2.0 / 3.0 * coeffBathy *
                ( pow(age2, 1.5) - pow(age1, 1.5) ) /
                ( age2 - age1 ) * surf;
                
                
                /* Plate Model PM: trapeze */
                double b1 =
                _earth->computeDepthBelowRidgePM(ages[i].getEffValue(),
                                                 _earth->getZM95());
                double b2 =
                _earth->computeDepthBelowRidgePM(ages[i+1].getEffValue(),
                                                 _earth->getZM95());
                
                _VbelowRidge_PM95 += 0.5 * (b1 + b2) * surf;
                
                
                b1 =
                _earth->computeDepthBelowRidgePM(ages[i].getEffValue(),
                                                 _earth->getZM125());
                b2 =
                _earth->computeDepthBelowRidgePM(ages[i+1].getEffValue(),
                                                 _earth->getZM125());
                
                _VbelowRidge_PM125 += 0.5 * (b1 + b2) * surf;
            }
            else
            {
                _VbelowRidge_HSCM += coeffBathy * sqrt(age1) * surf;
                
                _VbelowRidge_PM95 +=
                _earth->computeDepthBelowRidgePM(ages[i].getEffValue(),
                                                 _earth->getZM95()) * surf;
                
                _VbelowRidge_PM125 +=
                _earth->computeDepthBelowRidgePM(ages[i].getEffValue(),
                                                 _earth->getZM125()) * surf;
            }
        }
        
    } //!isContinental()
}
// --------------------------------------------------
void
PlateSection::removeExtraAges()
{
    /* remove points that were added for the computation of heat flux,
     or that have extra=true because they are outside the section */
    vector<Age>& ages = _rightElement->getLeftAges();
    vector<Age>::iterator it = ages.begin();
    while(true)
    {
        if((*it).isExtra())
        {
            ages.erase(it);
            it --; // in case two consecutive points are "extra"
        }
        it ++;
        if(it == ages.end())
            break;
    }
    _leftElement->setRightAges(ages);
}
// --------------------------------------------------
void
PlateSection::moveAges()
{
    /* remove points that were added by effectiveAges()
     for the computation of heatflux */
    if(_earth->insertAges())
        removeExtraAges();
    
    /* remove points that are exactly at ridges before moving,
     in order not to create many many close points. */
    vector<Age>& ages = _rightElement->getLeftAges();
	
	if (ages.empty()) {
		cerr << "WARNING: section " << _id << " has empty vector of ages " << endl;
	}
	
	if(!_earth->noSubduction())
    {
        if(_rightElement->isA("Ridge"))
        {
            Ridge* ridge = (Ridge*)_rightElement;
            if(ridge->isActive())
                ages.erase(ages.begin());
        }
        if(_leftElement->isA("Ridge"))
        {
            Ridge* ridge = (Ridge*)_leftElement;
            if(ridge->isActive())
                ages.pop_back();
        }
    }
    
    /* getU() is in cm/yr and dt is in yrs. displacement is in degrees */
    double displacement = 0.0;
    if(!_earth->fixedConfiguration() || _plate->hasType("RS"))
    {
        displacement =
        _earth->cm_to_deg(_plate->getU()) * _earth->getDt();
    }
    
    for(unsigned int i=0; i<ages.size(); i++)
    {
        ages[i].move(displacement);
        if(ages[i].isOceanic())
            ages[i].addValue(_earth->getDt() * 1.0E-6);  // _dt is in years and age in Myr
    }
    
    /* create new ages at ridges */
    createNewOceanAges();
    
    /* create new continental ages
     if growing continents */
    createNewContinentAges();
    
    /* remove points that are subducted */
    removeSubductedAges();
    
    /* effectiveAges for heat flux and plate thickness */
    effectiveAges();
    
    /* give section ages to elements */
    _rightElement->setLeftAge(ages[0]);
    _leftElement->setRightAge(ages[ages.size()-1]);
    _leftElement->setRightAges(ages);
    
    /* if this is a continental section and RCE and/or LCE have neighbors,
     then RCE and LCE are not up to date for ages
     (sections are defined on the neighbor then, not on the extremity).
     Put ages on RCE / LCE also: */
    if (isContinental())
    {
        _continent->getRightExtremity()->setLeftAge(ages[0]);
        _continent->getRightExtremity()->setLeftAges(ages);
        _continent->getLeftExtremity()->setRightAge(ages[ages.size()-1]);
        _continent->getLeftExtremity()->setRightAges(ages);
    }    
}
// --------------------------------------------------
void
PlateSection::removeSubductedAges()
{
    if(isContinental())
        return;
    
    /* remove ages points that were subducted.
     Using the bool extra to mark the points
     that have to be removed. */
    int iNotFound = -1;
    int iRight = iNotFound;
    int iLeft  = iNotFound;
    
    bool addRightAge = false;
    bool addLeftAge  = false;
    
    vector<Age>& ages = _rightElement->getLeftAges();
    
    // on the right:
    if(_rightElement->isA("RightSubduction") || _rightElement->isA("Staple"))
    {
        for(unsigned int i=1; i<ages.size(); i++)
        {
            if( _rightElement->isWithinOrEqual(ages[i-1].getPosition(),
                                               ages[i].getPosition(),
                                               LEFT ) )
            {
                iRight = i;
                break;
            }
            /* stop the search for iRight when far enough in the section
             (for very large sections that might loop around) */
            if(ages.size() > 20 && i > 20 &&
               contains(ages[i].getPosition()) &&
               contains(ages[i-1].getPosition()))
                break;
        }
        
        if(iRight != iNotFound)
        { // iRight is the first position inside the section
            addRightAge = !::areEqual(ages[iRight].getPosition(),
                                      _rightElement->getPosition());
            for(unsigned int i=0; i<iRight; i++)
                ages[i].setExtra(true);
        }
        
        /* find age exactly at the position of the right element */
		Age rightElementAge = Age(_earth, 0.0, 0.0);
        if(addRightAge)
        {
            rightElementAge.setPosition(_rightElement->getPosition());
			rightElementAge.setValue(_earth->linearAge(rightElementAge.getPosition(),
													   ages[iRight-1],
													   ages[iRight]));
        }        
        /* remove points outside of the section now */
        removeExtraAges();
		
		if (addRightAge)
		{ // put the age exactly at the position of rightElement;
			ages.insert(ages.begin(), rightElementAge);
		}
		
    } // rightElement is a subduction or staple
    
    // on the left:
    if(_leftElement->isA("LeftSubduction") || _leftElement->isA("Staple"))
    {
        for(unsigned int i=ages.size()-1; i>0; i--)
        {
            if( _leftElement->isWithinOrEqual(ages[i-1].getPosition(),
                                              ages[i].getPosition(),
                                              LEFT ) )
            {
                iLeft = i-1;
                break;
            }
            
            /* stop the search for iLeft when far enough in the section
             (for very large sections that might loop around) */
            if( ages.size() > 20 && i < ages.size()-20 &&
               contains(ages[i].getPosition()) &&
               contains(ages[i-1].getPosition()) )
                break;
        }
        
        if(iLeft != iNotFound)
        { //iLeft is the last position inside the section
            addLeftAge = !::areEqual(ages[iLeft].getPosition(),
                                     _leftElement->getPosition());
            for(unsigned int i=ages.size()-1; i>iLeft; i--)
                ages[i].setExtra(true);
        }
        
        /* find age exactly at the position of the left element */
		Age leftElementAge = Age(_earth, 0.0, 0.0);
        if(addLeftAge)
        {
			leftElementAge.setPosition(_leftElement->getPosition());
			leftElementAge.setValue(_earth->linearAge(leftElementAge.getPosition(),
													  ages[iLeft],
													  ages[iLeft+1]));
        }
        /* remove points outside of the section now */
        removeExtraAges();
		
		if (addLeftAge) 
		{ // put the age exactly at the position of the element;
			ages.push_back(leftElementAge);
		}
		
    } // leftElement is a subduction or a staple
    
    /* redefine ages at the left and right elements: */
    if(addRightAge)
        _rightElement->setLeftAge(ages[0]);
    if(addLeftAge)
        _leftElement->setRightAge(ages[ages.size()-1]);
    
    /* gives these "section ages" to left element also */
    _leftElement->setRightAges(ages);
}
// --------------------------------------------------
void
PlateSection::createNewOceanAges()
{
    if(isContinental())
        return;
    
    if(!_earth->noSubduction())
    {
        vector<Age>& ages = _rightElement->getLeftAges();
        
        if( _rightElement->isA("Ridge") &&
           ((Ridge*)_rightElement)->isActive())
        {
            Age ridgeAge = Age(_earth, 0.0, _rightElement->getPosition());
            
            while(true)
            {
                double position = ::addPosition(ages[0].getPosition(),
												-_earth->get_resolution());
                
                if( !contains(position) ||
                   _rightElement->isWithinOrEqual(position,
                                                  ages[0].getPosition(),
                                                  RIGHT) )
                {
                    ages.insert(ages.begin(), ridgeAge);
                    break;
                }
                else
                {
                    double value =
                    _earth->linearAge( position,
                                      ridgeAge,
                                      ages[0] );
                    Age newAge = Age(_earth, value, position);
                    ages.insert(ages.begin(), newAge);
                }
                
            } // while(true)
        } // _rightElement->isA("Ridge")
        
        if( _leftElement->isA("Ridge") &&
           ((Ridge*)_leftElement)->isActive())
        {
            Age ridgeAge = Age(_earth, 0.0, _leftElement->getPosition());
            while(true)
            {
                double position = ::addPosition(ages[ages.size()-1].getPosition(),
												_earth->get_resolution());
                
                if( !contains(position) ||
                   _leftElement->isWithinOrEqual(ages[ages.size()-1].getPosition(),
                                                 position,
                                                 LEFT) )
                {
                    ages.push_back(ridgeAge);
                    break;
                }
                else
                {
                    double value = _earth->linearAge(position,
                                                     ages[ages.size()-1],
                                                     ridgeAge);
                    Age newAge = Age(_earth, value, position);
                    ages.push_back(newAge);
                }
            }
        } // _leftElement->isA("Ridge")
        
		_rightElement->setLeftAge(ages[0]);
		_leftElement->setRightAge(ages[ages.size()-1]);
		
        _leftElement->setRightAges(ages);
    } // !_earth->noSubduction()
}
// --------------------------------------------------
void
PlateSection::createNewContinentAges()
{
    if( _earth->noSubduction() ||
       !_earth->growingContinents() ||
       !isContinental() )
        return;
    
    Continent* continent = getContinent();
    double resolution = _earth->get_resolution();
    RightContinentExtremity* rightExtremity =
    continent->getRightExtremity();
    LeftContinentExtremity* leftExtremity =
    continent->getLeftExtremity();
    
    vector<Age>& ages = _rightElement->getLeftAges();
    
    if(leftExtremity->growing())
    {
        double distance =
        ::computeRightLeftDistance( ages[ages.size()-2].getPosition(),
                                   ages[ages.size()-1].getPosition() );
        if(::isStrictlyLess(distance, resolution))
        { // remove last point to avoid having many close points
            ages.pop_back();
        }
        
        Age leftAge =  Age(_earth, _earth->getAgeMa(),
                           leftExtremity->getPosition());
        leftAge.setOceanic(false);
        while(true)
        {
            double lastPos = ages[ages.size()-1].getPosition();
            double newPos = ::addPosition(lastPos, resolution);
            if( !contains(newPos) ||
               leftExtremity->isWithinOrEqual( lastPos, newPos, LEFT ) )
            { /* leftAge is less than "resolution" away from
               last point: add only one point */
                ages.push_back(leftAge);
                break;
            }
            else
            { /* leftAge is more than resolution away from last
               point: add several points */
                double value =
                _earth->linearAge( newPos,
                                  ages[ages.size()-1],
                                  leftAge );
                Age newAge = Age(_earth, value, newPos);
                newAge.setOceanic(false);
                ages.push_back(newAge);
            }
        } // while(true)
    } // leftExtremity->growing()
    
    if(rightExtremity->growing())
    {
        double distance =
        ::computeRightLeftDistance(ages[0].getPosition(),
                                   ages[1].getPosition());
        if(::isStrictlyLess(distance, resolution))
        { // remove first point to avoid many close points
            ages.erase(ages.begin());
        }
        
        Age rightAge = Age(_earth, _earth->getAgeMa(),
                           rightExtremity->getPosition());
        rightAge.setOceanic(false);
        while(true)
        {
            double firstPos = ages[0].getPosition();
			double newPos = ::addPosition(firstPos, -resolution);

            if( !contains(newPos) ||
               rightExtremity->isWithinOrEqual( newPos, firstPos, RIGHT ) )
            { /* rightAge is less than resolution away from
               first point: add only one point */
                ages.insert(ages.begin(), rightAge);
                break;
            }
            else
            { /* rightAge is more than resolution away from
               first point: add several points */
                double value =
                _earth->linearAge( newPos,
                                  rightAge,
                                  ages[0] );
                Age newAge = Age(_earth, value, newPos);
                newAge.setOceanic(false);
                ages.insert(ages.begin(), newAge);
            }
        } //while(true)
    } // rightExtremity->growing()
    
	_rightElement->setLeftAge(ages[0]);
	_leftElement->setRightAge(ages[ages.size()-1]);
	
    _leftElement->setRightAges(ages);
}
// --------------------------------------------------
void
PlateSection::correctAges()
{ 
    /* after a collision or continental breakup:
     because of numeric error or slight continental motion,
     the positions of ages[0] and ages[size-1] are
     not exactly at the elements' positions.
     */
    vector<Age>& ages = _rightElement->getLeftAges();
    
    if(ages.size() > 2)
        removeSubductedAges();
    
    ages[0].setPosition(_rightElement->getPosition());
    ages[ages.size()-1].setPosition(_leftElement->getPosition());
    
    /* can also happen that a value = -0.000 -> problem with sqrt then */
    for (unsigned int i=0; i<ages.size(); i++)
    {
        if (ages[i].getValue() < 0.0)
        {
            if (::isStrictlyLess(ages[i].getValue(), -Earth::agePrecision))
            { // make a warning only if significative negative age
                cout << "WARNING: negative age = "
                << ages[i].getValue() << " at position " << ages[i].getPosition() << endl;
            }
            ages[i].setValue(abs(ages[i].getValue()));
        }
    }
    
    /* effectiveAges for heat flux and plate thickness */
    effectiveAges();
    
    /* give section ages to elements */
    _rightElement->setLeftAge(ages[0]);
    _leftElement->setRightAge(ages[ages.size()-1]);
    _leftElement->setRightAges(ages);
    
}
// --------------------------------------------------
void
PlateSection::effectiveAges()
{
    vector<Age>& ages = _rightElement->getLeftAges();
    
    if(isContinental())
    {// put all effValue = -1.0
        for(unsigned int i=0; i<ages.size(); i++)
        {
            ages[i].setOceanic(false);
            ages[i].setEffValue(-1.0);
        }
        return;
    }
    
    _earth->effectiveAges(ages);
    
}
// --------------------------------------------------
void
PlateSection::update_U_for_Ages()
{
    vector<Age>& ages = _rightElement->getLeftAges();
    for(unsigned int i=0; i<ages.size(); i++)
        ages[i].setU(_plate->getU());
    _leftElement->setRightAges(ages);
}
// --------------------------------------------------
bool
PlateSection::contains(double position)
{ /* position is strictly inside the section */
    if(_rightElement == _leftElement)
        return true;
    
    return ::isStrictlyWithin(position,
                              _rightElement->getPosition(),
                              _leftElement->getPosition());
}
// --------------------------------------------------
vector<Age>
PlateSection::splitAges(double breakupPos, double offset, Direction direction)
{
    /* for continental break up: split ages in the left or
     right part of the section, and move all the ages' points
     by +/- offset. */
    
    vector<Age> splitAges;
    
    if(!isContinental())
    {
        cerr << "WARNING: PlateSection::splitAges() is used on a non-continental section\n"
        << "-> nothing done"
        << endl;
        return splitAges; // nothing defined
    }
    
    if( !::isStrictlyWithin(breakupPos,
                            _rightElement->getPosition(),
                            _leftElement->getPosition()) )
    {
        cerr << "ERROR in PlateSection::splitAges()\n"
        << "the breakup position is not within the section"
        << endl;
        exit(EXIT_FAILURE);
    }
    
    vector<Age> ages = _rightElement->getLeftAges();
    
    int iLimit = -1;
    for(unsigned int i=0; i<ages.size()-1; i++)
    {
        if ( ::isWithinOrEqual( breakupPos,
                               ages[i].getPosition(),
                               ages[i+1].getPosition(),
                               RIGHT ) )
        { // breakupPos is between i and i+1 position
            iLimit = i; /* last position within the new right continent or
                         first position out of the section for new
                         left continent */
            break;
        }
        
        if(i == ages.size()-1)
        {
            cerr << "ERROR in PlateSection::splitAges(): breakupPos was not found within the section\n";
            cerr << "\t breakupPos: " << breakupPos << endl;
            cerr << "\t rightElement: " << _rightElement->getShortId()
            << "\t" << _rightElement->getPosition()
            << endl;
            cerr << "\t leftElement: " << _leftElement->getShortId()
            << "\t" << _leftElement->getPosition()
            << endl;
            exit(EXIT_FAILURE);
        }
    }
    
    switch(direction)
    {
        case RIGHT :
        {
            for(unsigned int i=0; i<=iLimit; i++)
            {
                splitAges.push_back(ages[i]);
            }
            
            if(!::areEqual(breakupPos, splitAges[splitAges.size()-1].getPosition()))
            { /* add one point for the border,
               with linear interpolation for the age */
                double value =
                _earth->linearAge( breakupPos,
                                  ages[iLimit],
                                  ages[iLimit+1] );
                Age limitAge = Age(_earth, value, breakupPos);
                limitAge.setOceanic(false);
                splitAges.push_back(limitAge);
            }
            splitAges[splitAges.size()-1].setPosition(breakupPos);
            
            // move the positions by -offset
            for(unsigned int i=0; i<splitAges.size(); i++)
                splitAges[i].move(-offset);
            
            break;
        }
        case LEFT :
        {
            for(unsigned int i=iLimit+1; i<ages.size(); i++)
            {
                splitAges.push_back(ages[i]);
            }
            
            if(!::areEqual(breakupPos, splitAges[0].getPosition()))
            {
                double value =
                _earth->linearAge( breakupPos,
                                  ages[iLimit],
                                  ages[iLimit+1] );
                Age limitAge = Age(_earth, value, breakupPos);
                limitAge.setOceanic(false);
                splitAges.insert(splitAges.begin(), limitAge);
            }
            splitAges[0].setPosition(breakupPos);
            
            // move the positions by +offset
            for(unsigned int i=0; i<splitAges.size(); i++)
                splitAges[i].move(offset);
            
            break;
        }
        default :
        {
            cerr << "ERROR : PlateSection::splitAges() is used without direction" << endl;
            exit(EXIT_FAILURE);
        }
    }
    
    return splitAges;
}
// --------------------------------------------------
void
PlateSection::setAgeGroups()
{ /* _ageGroups contains the max age in a slice of ages,
   and the right and left positions of a portion of continent
   within the considered age slice. */
    if(isContinental())
        return;
    
    _ageGroups.clear();
    
    vector<Age> ages = _rightElement->getLeftAges();
    
    if(ages.size() < 1)
        return;
    
    vector<pair<double, double> > slices =
    _earth->getOceanAgeSlices();
    
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
        group.ageMax = slices[iSlice].first;
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
// --------------------------------------------------
double
PlateSection::computeMeanThickness()
{ // results in m
    double meanH = 0.0;
    
    if(!isContinental())
    {
        vector<Age> ages = _rightElement->getLeftAges();
        
        double coeff = Earth::coeff_plate_thickness * sqrt(_earth->get_kappa());
        double twothird = 2.0 / 3.0;
        
        for(unsigned int i=0; i<ages.size()-1; i++)
        {
            double L =
            ::computeRightLeftDistance(ages[i].getPosition(),
                                       ages[i+1].getPosition());
            
            double age1 = ::Myr_to_sec(ages[i].getEffValue());
            if(!::areEqual(ages[i].getEffValue(),
                           ages[i+1].getEffValue(),
                           1.0E-4)) // 1.0E-4Myr = 100yr
            { // integral between age[i] and age[i+1]
                double age2 = ::Myr_to_sec(ages[i+1].getEffValue());
                meanH += twothird * coeff *
                ( pow(age2, 1.5) - pow(age1, 1.5) )
                / ( age2 - age1 ) * L;
            }
            else
            { // simple
                meanH += coeff * sqrt(age1) * L;
            }
        } // ages[i]
    }
    
    
    meanH /= computeLength();
    
    return meanH;
}
// ==================================================
