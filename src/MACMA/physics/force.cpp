// This file is part of MACMA
// 
// author: cecile.grigne@univ-brest.fr
//

#include "MACMA/physics/earth.h"
#include "MACMA/physics/force.h"
#include "MACMA/physics/geoElement.h"
#include "MACMA/physics/subduction.h"
#include "MACMA/physics/rightSubduction.h"
#include "MACMA/physics/leftSubduction.h"
#include "MACMA/physics/ridge.h"
#include "MACMA/physics/enum_struct.h"

using namespace std;

// ==================================================
Force::Force(Earth* earth) : _earth(earth)
{
    _ridgePush = 0.0;
    _slabPull = 0.0;
    _slabSuction = 0.0;
    _etaH = 0.0;
    _etaVRight = 0.0;
    _etaVLeft  = 0.0;
    _etaBRight = 0.0;
    _etaBLeft  = 0.0;
    _mantleDrag = 0.0;
    _viscousShear = 0.0;
    _bending = 0.0;
    _URight = 0.0;
    _ULeft  = 0.0;
}

Force::~Force()
{
}

// ------------------------------
void
Force::zeroAll()
{
    _ridgePush = 0.0;
    _slabPull = 0.0;
    _slabSuction = 0.0;
    _etaH = 0.0;
    _etaVRight = 0.0;
    _etaVLeft  = 0.0;
    _etaBRight = 0.0;
    _etaBLeft  = 0.0;
    _mantleDrag = 0.0;
    _viscousShear = 0.0;
    _bending = 0.0;
    _URight = 0.0;
    _ULeft  = 0.0;
}
// -------------------------------------------------------
void 
Force::addDrivingForces(GeoElement* element, Direction direction)
{ /* updates forces due to element, taking into consideration
   whether this element is on the LEFT or the RIGHT of the plate
   */
    switch(direction)
    {
        case LEFT :
        {
            if(element->isA("Subduction"))
            {
                Subduction* subduction = (Subduction*)element;
                _slabSuction += subduction->computeSlabSuction();
                if(subduction->isA("LeftSubduction"))
                {
                    LeftSubduction* leftSubduction = (LeftSubduction*)subduction;
                    _slabPull += leftSubduction->computeSlabPull();
                }
            }
            else if(element->isA("Ridge"))
            {
                Ridge* ridge = (Ridge*)element;
                _ridgePush -= ridge->computeRidgePush(RIGHT);
            }
            break;
        }
        case RIGHT :
        {
            if(element->isA("Subduction"))
            {
                Subduction* subduction = (Subduction*)element;
                _slabSuction -= subduction->computeSlabSuction();
                if(subduction->isA("RightSubduction"))
                {
                    RightSubduction* rightSubduction = (RightSubduction*)subduction;
                    _slabPull -= rightSubduction->computeSlabPull();
                }
            }
            else if(element->isA("Ridge"))
            {
                Ridge* ridge = (Ridge*)element;
                _ridgePush += ridge->computeRidgePush(LEFT);
            }
            break;
        }
        default :
        {
            cerr << "ERROR: Force::addDrivingForces(direction) is used without direction given\n";
            exit(EXIT_FAILURE);
        }
    } // switch(direction)
}
// -------------------------------------------------------
void
Force::computeViscosities(GeoElement* element, Direction direction)
{
    if(element->isA("Subduction"))
    {
        Subduction* subduction = (Subduction*)element;
        switch(direction)
        {
            case LEFT :
            {
                _etaVLeft = subduction->computeVerticalViscosity();
                _etaBLeft = subduction->computeBending();
                break;
            }
            case RIGHT :
            {
                _etaVRight = subduction->computeVerticalViscosity();
                _etaBRight = subduction->computeBending();
                break;
            }
            default :
            {
                cerr << "ERROR: Force::computeViscosities(direction) is used without direction given\n";
                exit(EXIT_FAILURE);
            }
        } // switch(direction)
    }
}
// -------------------------------------------------------
