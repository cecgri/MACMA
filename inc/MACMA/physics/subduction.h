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

#ifndef SUBDUCTION_H_
#define SUBDUCTION_H_

#include <string>
#include <vector>
#include "MACMA/physics/interface.h"
#include "MACMA/physics/enum_struct.h"
using std::string;
using std::vector;
using std::pair;

class RightContinentExtremity;
class LeftContinentExtremity;
class ActiveMargin;

/*
 ==================================================
 Class Subduction : extends Interface, abstract class
 -------------------------------------------------
 - Unify left and right subductions
 - Has a depth
 ==================================================
 */

class Subduction : public Interface
{
public:
    virtual ~Subduction();
    
    // --- Interface
    virtual bool isA(string className);
    virtual string getClass();
    virtual void updatePosition(double dt) = 0;
    virtual double getNextPosition(double dt);
    virtual void updateVelocity() = 0;
    
    // --- Subduction
    inline double getDepth() { return _depth; }
    inline void setDepth(double depth) { _depth = depth; }
    
    // --- Thickness of slab at different depths
    inline vector<pair<double, double> > getThickness() { return _thickness; }
    inline void setThickness(vector<pair<double, double> > thick) { _thickness = thick; }
    double getMeanThickness();
    void initThickness(double, double);
    void updateThickness(double);
    inline double getThicknessAtDepth() { return _thicknessAtDepth; }
    inline void setThicknessAtDepth(double H) { _thicknessAtDepth = H; }
    
    // --- Continents
    inline void setLeftNeighborExtremity(RightContinentExtremity* extremity) { _leftNeighborExtremity = extremity; }
    inline RightContinentExtremity* getLeftNeighborExtremity() { return _leftNeighborExtremity; }
    inline void setRightNeighborExtremity(LeftContinentExtremity* extremity) { _rightNeighborExtremity = extremity; }
    inline LeftContinentExtremity* getRightNeighborExtremity() { return _rightNeighborExtremity; }
    
    // --- markers
    inline void setActiveMargin(ActiveMargin* margin) { _activeMargin = margin; }
    inline ActiveMargin* getActiveMargin() { return _activeMargin; }
    
    // compute -----
    double computeSlabSuction();
    double computeSlabPull();
    double computeBending();
    double computeVerticalViscosity();
    
protected:
    Subduction(Earth* earth, double position);
    double _depth;
    
    RightContinentExtremity* _leftNeighborExtremity;
    LeftContinentExtremity* _rightNeighborExtremity;
    
    vector<pair<double, double> > _thickness;
    double _thicknessAtDepth;
    ActiveMargin* _activeMargin;
    
};

#endif
