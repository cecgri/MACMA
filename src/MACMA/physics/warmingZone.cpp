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
#include "MACMA/physics/warmingZone.h"

using namespace std;
// ==================================================
WarmingZone::WarmingZone(Earth* earth, Continent* continent): 
_earth(earth), 
_continent(continent)
{
	_id = Earth::getId("WarmingZone");
    
	_lifetime = 0.0;
	_warming = 0.0;
	
	_thermalPressure = computeThermalPressure();
}

WarmingZone::~WarmingZone()
{
}

/*
 ============================
 Computing
 ============================
 */
double
WarmingZone::computeThermalPressure()
{ /* term used in the expression of warming.
   Unit is Pa/K */
	double h = _earth->get_subcont_warming_H(); // in m
	double a = _continent->getLm();  // in m
	double a_h = 0.5 * a / h;
	
	double rho = _earth->get_rho_um_p();
	
	double Pth = 2.0 * h * _earth->get_alpha_um() *
	rho * _earth->get_g() /
	( pow(a_h, 2.0) + pow(a_h, -2.0) );
	
	return Pth;
}

double 
WarmingZone::computeTmax()
{
	double eta = _earth->getEtaUm() * Earth::eta0;
	return sqrt(_earth->radioHeatWperKg() * eta /
				(_earth->get_C_p() * _thermalPressure) );
}

double
WarmingZone::computeWarmingTimeMyr()
{
	double eta = _earth->getEtaUm() * Earth::eta0;
	return ::sec_to_Myr( sqrt( _earth->get_C_p() * eta / 
							  (_earth->radioHeatWperKg() * _thermalPressure) ) ); 
}

void
WarmingZone::updateWarming()
{
	// maximum warming temperature: 
	double Tmax = computeTmax();
	
	// characteristic time of warming:
	double tauMyr = computeWarmingTimeMyr();
	
	_warming = Tmax * tanh( _lifetime / tauMyr);
}

void
WarmingZone::updateF()
{ /* */
	_F = 0.25 * pow( _continent->getLm(), 2.0 ) / _earth->get_subcont_warming_H() *
    _thermalPressure * _warming;
}


/*
 ========================================
 WarmingZone:: lifetime
 -------------------------
 ========================================
 */

void
WarmingZone::updateLifetime(double dt)
{ // dt is in years
	_lifetime += dt / 1.0E6; // in Myr
}

void
WarmingZone::computeLifetime0(double warming)
{
	// Time since the beginning of warming
	double Tmax = computeTmax();
	double tauMyr = computeWarmingTimeMyr();
	
	_lifetime = 0.5 * tauMyr * 
    log( (Tmax + warming) / (Tmax - warming) );
}
