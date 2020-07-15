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

#ifndef CONTINENT_EXTREMITY_H_
#define CONTINENT_EXTREMITY_H_

#include <string>
#include "MACMA/physics/geoElement.h"

using std::string;

class Earth;
class Continent;

/*
 ==================================================
 Class ContinentExtremity : extends GeoElement
 -------------------------------------------------
 ==================================================
 */

class ContinentExtremity : public GeoElement
{
public:
    ~ContinentExtremity();
    
    // static for continental growth
    static double referenceGrowthRate;
    static double minDepthSubduction;
    
    // --- GeoElement
    virtual void updateVelocity();
    virtual bool isA(string className);
    virtual string getClass();
    
    inline void setContinent(Continent* continent) { _continent = continent; }
    inline Continent* getContinent() { return _continent; }
    
    // --- ContinentExtremity
    inline void setNeighbor(GeoElement* element) { _neighbor = element; }
    inline GeoElement* getNeighbor() { return _neighbor; }
    
    // --- Continental growth
    virtual double growthRate(); // in degrees
    virtual void updateGrowthRate();
    inline double getGrowthRate() { return _growthRate; };
    inline bool growing() { return _growing; }
    virtual void modifyVelocity() = 0;
    
    // --- Marker
    inline void setActiveMargin(ActiveMargin* margin) { _activeMargin = margin; }
    inline ActiveMargin* getActiveMargin() { return _activeMargin; }
    
protected:
    ContinentExtremity(Earth* earth, Continent* continent);
    Continent* _continent;
    
    GeoElement* _neighbor;
    
    ActiveMargin* _activeMargin;  // just a tracer
    
    double _growthRate;
    bool _growing;
};

#endif
