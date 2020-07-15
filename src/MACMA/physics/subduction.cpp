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
#include "MACMA/physics/subduction.h"
#include "MACMA/physics/continent.h"
#include "MACMA/physics/enum_struct.h"

using namespace std;
// ==================================================
Subduction::Subduction(Earth* earth, double position) : Interface(earth, position), _depth(0.0)
{
    _id = Earth::getId("Subduction");
    
    _leftNeighborExtremity = NULL;
    _rightNeighborExtremity = NULL;
    
    _depth = 0.0;
    _thicknessAtDepth = 0.0;
    
    _activeMargin = NULL;
}

Subduction::~Subduction()
{
}

bool
Subduction::isA(string className)
{
    return className == "Subduction" || className == "Interface" || className == "GeoElement";
}

string
Subduction::getClass()
{
    return "Subduction";
}

// --------------------------------------------------
double 
Subduction::getNextPosition(double dt)
{
    return ::addPosition(_position, _earth->cm_to_deg(_U) * dt);
}
// --------------------------------------------------
double
Subduction::computeSlabPull()
{
    double slabPull = 0.0;
    
    double Z = 0.0;
    if(_earth->getSlabPullDepth() == Earth::UPPER_MANTLE)
        Z = min(_depth, _earth->get_D());
    else
        Z = _depth;
    
    double H = getMeanThickness();
    
    // TEST ///////////////
    H = max(H, 10.0E3); // at least 10km for the plate here.
    ///////////////////////
    
    return _earth->computeSlabPull(H, Z); // is with dimensions
}
// --------------------------------------------------
double
Subduction::computeSlabSuction()
{ /* computed using eta_m_p * V_sink_p (present reference values),
   which is equivalent to considering that DeltaRho * H
   (DeltaRho: density contrast; H: thickness of the slab) does not
   vary over time (Stokes velocity so that
   DeltaRho * H * g * Z ~ eta_m * V_sink / (d * Z )
   */
    
    double slabSuction = 0.0;
    
    double eta_m_p = _earth->get_eta_m_p();
    double Vsink_p   = _earth->get_V_sink_p();
    double width = _earth->get_d() - _earth->get_D(); // thickness of the lower mantle
    double Z = max(_depth - _earth->get_D(), 0.0); // depth in lower mantle only
    
    slabSuction = _earth->getCoeffSlabSuction() * 2.0 / 9.0 *
    eta_m_p * Vsink_p * Z / width;  // 2/9: not exact coeff for this shape (Stokes' velocity for a sphere)
    
    slabSuction *= Earth::force0; // eta and V are adimensional;
    
    return slabSuction;
}
// --------------------------------------------------
double 
Subduction::computeVerticalViscosity()
{
    double etaV = 0.0;
    
    double eta_um  = _earth->getEtaUm();
    double eta_ast = _earth->getEtaAst();
    double eta_m   = _earth->getEtaM();
    double width   = _earth->get_D(); // to divide by 2?
    double asth    = _earth->get_thick_ast();
    double D       = _earth->get_D();
    
    double Z = 0.0;
    double d1 = 0.0;
    double d2 = 0.0;
    double d3 = 0.0;
    
    if(_earth->getViscousShearDepth() == Earth::UPPER_MANTLE)
        Z = min(_depth, D);
    else
        Z = _depth;
    
    d1 = min(Z, asth); // thickness in asthenosphere if it exists.
    d2 = min(D - asth, max(0.0, Z - asth));  // thickness in upper mantle without asthenosphere
    
    if(_earth->getViscousShearDepth() == Earth::WHOLE_MANTLE)
        d3 = max(0.0, Z - D); // thickness in lower mantle
    
    etaV = _earth->getCoeffViscousShear() * 2.0 *
    (eta_ast * d1 + eta_um * d2 + eta_m * d3) / width;
    
    return etaV;
}
// --------------------------------------------------
double
Subduction::computeBending()
{
    double bending = 0.0;
    double H = 0.0;
    if(this->isA("LeftSubduction"))
        H = computePlateThickness(RIGHT);
    else
        H = computePlateThickness(LEFT);
    
    bending = (2.0/3.0) * _earth->getEtaPl() * pow(H / _earth->get_R_min(), 3.0);
    bending *= _earth->getCoeffBending();
    
    return bending;
}
// --------------------------------------------------

/* -----------------------------------
 | Thickness of the subducted slabs |
 ----------------------------------- */

void
Subduction::initThickness(double depth, double thickness)
{ /* CG: _thickness: vector of pairs, with .first: depth, 
   and .second: thickness at this depth */
    _thickness.push_back(make_pair(depth, thickness));
    if(depth > 0.0)
    {
        double thickness0 = 0.0;
        if(this->isA("LeftSubduction"))
            thickness0 = computePlateThickness(RIGHT);
        else
            thickness0 = computePlateThickness(LEFT);
        
        // linear
        double depthResolution = _earth->deg_to_m(_earth->get_resolution());
        double dth = (thickness - thickness0) / depth * depthResolution;
        
        double d = depth;
        double thick = thickness;
        while(true)
        {
            d -= depthResolution;
            thick -= dth;
            if(d < depthResolution)
                break;
            _thickness.insert(_thickness.begin(), make_pair(d, thick));
        }
    }
}

void
Subduction::updateThickness(double deep)
{ // deep: amount of sinking in the mantle
    for(unsigned int i=0; i<_thickness.size(); i++)
        _thickness[i].first += deep;
    
    /* remove points below 2800.0E3 */
    while(true)
    {
        if(_thickness[_thickness.size()-1].first > 2800.0E3)
            _thickness.pop_back();
        else
            break;
    }
    
    /* add new points */
    double depthResolution = _earth->deg_to_m(_earth->get_resolution());
    if(_thickness[0].first > depthResolution)
    {
        double thickness0 = 0.0;
        if(this->isA("LeftSubduction"))
            thickness0 = computePlateThickness(RIGHT);
        else
            thickness0 = computePlateThickness(LEFT);
        
        double d = _thickness[0].first;
        double thick = _thickness[0].second;
        double dth = (thick - thickness0) / d * depthResolution;
        while(true)
        {
            d -= depthResolution;
            thick -= dth;
            if(d < depthResolution)
                break;
            _thickness.insert(_thickness.begin(), make_pair(d, thick));
        }
    }
}
// --------------------------------------------------
double
Subduction::getMeanThickness()
{
    double thickness = 0.0;
    double addDepth = 0.0;
    
    double totalDepth = 0.0;
    if(_earth->getSlabPullDepth() == Earth::WHOLE_MANTLE)
        totalDepth = _depth;
    else // upper mantle only
        totalDepth = min(_depth, _earth->get_D()); // max upper mantle thickness
    
    double totalThickness = 0.0;
    /* trapeze method */
    // from 0 to first recorded depth/thickness:
    double thickness0 = 0.0;
    if(this->isA("LeftSubduction"))
        thickness0 = computePlateThickness(RIGHT);
    else // right sub
        thickness0 = computePlateThickness(LEFT);
    
    totalThickness += (thickness0 + _thickness[0].second)*0.5 * _thickness[0].first;
    addDepth += _thickness[0].first;
    
    
    // middle trapezes
    if(_thickness.size() > 1)
    {
        unsigned int i = 1;
        while(true)
        {
            if(_thickness[i].first > totalDepth || i == _thickness.size()-1)
                break;
            
            double width = _thickness[i].first - _thickness[i-1].first;
            totalThickness += (_thickness[i-1].second + _thickness[i].second)*0.5 * width;
            addDepth += width;
            
            i ++;
        }
        _thicknessAtDepth = _thickness[i-1].second;
    }
    else
        _thicknessAtDepth = _thickness[0].second;
    
    return totalThickness / addDepth;
}
// --------------------------------------------------
