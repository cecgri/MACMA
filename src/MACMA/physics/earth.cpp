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
//   initial author: leyaouanq@cervval.com
//   modifications: cecile.grigne@univ-brest.fr

#include "MACMA/physics/earth.h"
#include <limits>

using namespace std;

/* =====================================
 constants
 ======================================
 */
const double seconds_per_yrs = 365.0 * 24.0 * 3600.0;
const double seconds_per_Myrs = 1.0E6 * seconds_per_yrs;
const double deg_per_rad = 180.0 / M_PI;

/* statics for time/age */
double Earth::zero = 0.0;
double Earth::loop = 360.0;
double Earth::tinyDouble = numeric_limits<double>::epsilon();
double Earth::ageEarth = 4500.0;
double Earth::contAgeGroupStep = 250.0;
double Earth::contAgeGroupMax = 5000.0;
double Earth::contAgeGroupMin = -2000.0;
double Earth::oceanAgeGroupStep = 25.0;
double Earth::oceanAgeGroupMax = 250.0;
double Earth::oceanAgeGroupMin = 0.0;

/* statics for dimensionless equations */
double Earth::eta0 = 1.0E22; // in Pa.s
double Earth::u0 = 0.01 / seconds_per_yrs; // (1cm/yr given in m/s)
double Earth::force0 = eta0 * u0; // in N/m (3.171TN/m with values above)

/* other statics */
double Earth::gasConstant = 8.314;
double Earth::creationDistance = 0.01; // in degrees
double Earth::positionPrecision = 1.0E-10; // precision for equality on positions, in degrees
double Earth::agePrecision = 1.0; // precision for equality of age, in year
double Earth::coeff_plate_thickness = 2.0 * sqrt(log(3.0));
double Earth::coeff_T_plate = 0.486065; /* erf(1) + (exp(-1)-1) / sqrt(pi) 
                                         coming from integral of
                                         erf(z/(2(sqrt(kappa t)) thru plate */

/*
 ==============================
 global (free) functions - CG
 ==============================
 */
bool areEqual(double a, double b)
{
    return fabs(a-b) < Earth::positionPrecision;
}
// ---------------------------------------------------
bool areEqual(double a, double b, double epsilon)
{ // CG checks if a = b more or less epsilon
    return fabs(a-b) < epsilon;
}
// --------------------------------------------------
bool isStrictlyLess(double a, double b)
{ // checks if a < b
    return a-b < -Earth::positionPrecision;
}
// --------------------------------------------------
bool isStrictlyLess(double a, double b, double epsilon)
{ // checks if a < b
    double  eps = fabs(epsilon);
    return a-b < -eps;
}
// --------------------------------------------------
bool isStrictlyMore(double a, double b)
{ // checks if a > b
    return a-b > Earth::positionPrecision;
}
// --------------------------------------------------
bool isStrictlyMore(double a, double b, double epsilon)
{ // checks if a > b
    double eps = fabs(epsilon);
    return a-b > eps;
}
// --------------------------------------------------
bool isStrictlyWithin(double pos, double right, double left)
{ // pos is strictly comprised between right and left
    if(::areEqual(right, left))
    {
        if(::areEqual(right, pos))
            return true;
        else
            return false;
    }
    
    if(::isStrictlyLess(left, right))
    {
        return (::isStrictlyLess(right, pos) || ::isStrictlyLess(pos, left));
    }
    else
    {
        return (::isStrictlyLess(right, pos) && ::isStrictlyLess(pos, left));
        
    }
}
// --------------------------------------------------
bool isWithinOrEqual(double pos, double right, double left)
{ // right <= pos <= left
    if(::isStrictlyWithin(pos, right, left))
        return true;
    
    if(::areEqual(pos, right))
        return true;
    
    if(::areEqual(pos, left))
        return true;
    
    return false;
}
// --------------------------------------------------
bool isWithinOrEqual(double pos, double right, double left, Direction dir)
{ // right < pos <= left (if dir==left) or right <= pos < left (if dir==right)
    if(::isStrictlyWithin(pos, right, left))
        return true;
    
    switch(dir)
    {
        case LEFT :
            return ::areEqual(pos, left);
        case RIGHT :
            return ::areEqual(pos, right);
        default :
            return false;
    }
}
// --------------------------------------------------
double computeAbsDistance(double posA, double posB)
{ // CG: absolute distance between two points ( 0 < distance < 180.0)
    double distance = min(360.0 - fabs(posA-posB), fabs(posA-posB));
    return distance;
}
// --------------------------------------------------
double computeDistance(double posA, double posB)
{ // CG: posA - posB, so that (-180 < distance < 180)
    /* if A is left of B: positive
     if B is left of A: negative */
    double distance = posA - posB;
    if(distance < -180.0)
        distance += 360.0;
    if(distance > 180.0)
        distance -= 360.0;
    
    return distance;
}
// --------------------------------------------------
double computeRightLeftDistance(double right, double left)
{ // distance from right to left position (0 < distance < 360)
    if(::isStrictlyLess(left, right))
        left += Earth::loop;
    return left - right;
}
// --------------------------------------------------
double getMiddle(double rightPos, double leftPos)
{
    if(::isStrictlyLess(leftPos, rightPos))
        leftPos += Earth::loop;
    
    double middle = (leftPos + rightPos) / 2.0;
    middle = fmod(middle, Earth::loop);
    
    return middle;
}
// --------------------------------------------------
void changePosition(double &position, double add)
{ // takes into account periodicity (0 <= p < 360)
    position = fmod(position + add, Earth::loop);
    if(::isStrictlyLess(position, Earth::zero))
        position += Earth::loop;
}
// --------------------------------------------------
double addPosition(double initPos, double add)
{
    double newPos = fmod(initPos + add, Earth::loop);
    if(::isStrictlyLess(newPos, Earth::zero))
        newPos += Earth::loop;
    return newPos;
}
// --------------------------------------------------
pair<double, double> getMeanSD(vector<double> vec)
{ /* computes mean and SD of a vec.
   getMeanSD.fist : mean
   getMeanSD.second : SD */
    pair<double, double> meanSD;
    if(vec.size() > 0)
    {
        meanSD.first = 0.0;
        meanSD.second = 0.0;
        
        for(unsigned int i=0; i<vec.size(); i++)
            meanSD.first += vec[i];
        meanSD.first = meanSD.first / vec.size();
        
        for(unsigned int i=0; i<vec.size(); i++)
            meanSD.second += pow(vec[i] - meanSD.first, 2);
        meanSD.second = sqrt(meanSD.second / vec.size());
    }
    return meanSD;
}
// --------------------------------------------------
vector<double> getWhiskers(vector<double> vec)
{ /* computes the decile 0.1 and 0.9, 
   and the quartile .25, .5 (median) and .75,
   e.g. for plotting a boxplot (whisker box).
   vec[1]: D1 (10%)
   vec[2]: Q1 (25%)
   vec[3]: Q2 (50%)
   vec[4]: Q3 (75%)
   vec[5]: D9 (90%)
   WARNING: vec must be in ascending order!
   */
    vector<double> whiskers;
    
    double percArray[] = {.1, .25, .5, .75, .9};
    vector<double> perc( percArray,
                        percArray + sizeof(percArray) / sizeof(percArray[0]) );
    
    for(unsigned int i=0; i<perc.size(); i++)
    {
        unsigned int iInf = (unsigned int)((vec.size()-1) * perc[i]);
        double p = (double)(vec.size()-1) * perc[i] - (double)(iInf);
        
        double value = vec[iInf]*(1.0-p) + vec[iInf+1]*p;
        whiskers.push_back(value);
    }
    
    return whiskers;
}
// --------------------------------------------------
vector<pair<double, double> > getDistribution(vector<double> vec, double slice)
{ /* used to plot histogram of ages. 
   In pair: first is age of each slice,
   second is number of times this age appears.
   */
    vector<pair<double, double> > distribution;
    int cpt = 0;
    int nbSlices = 0;
    unsigned int i=0;
    while(i < vec.size())
    {
        if(vec[i] < nbSlices * slice)
        {
            cpt ++;
            i ++;
        }
        else
        {
            distribution.push_back(make_pair(nbSlices * slice, (double)cpt));
            nbSlices ++;
            cpt = 0;
        }
    }
    
    return distribution;
}
// --------------------------------------------------
double lin_Pos_to_Y(double p, double posR, double posL, double yR, double yL)
{ /* linear interpolation for p between posR (on the right)
   and posL (on the left). Because of periodicity (0<pos<360),
   the point with the jump is chosen as the antipode of the middle
   between posR and posL (oppMiddle). */
    double halfloop = Earth::loop / 2.0;
    
    double oppMiddle = halfloop - getMiddle(posR, posL);
    changePosition(p, oppMiddle);
    changePosition(posR, oppMiddle);
    changePosition(posL, oppMiddle);
    
    if(::areEqual(posR, posL))
        posR += Earth::loop;
    
    return yR + (yL - yR) * (p - posR) / (posL - posR);
}
// --------------------------------------------------
double lin_Y_to_Pos(double y, double posR, double posL, double yR, double yL)
{ /* find position that gives y for the linear interpolation *
   between posR and posL (y going from yR to yL) */
    double oppMiddle = 180.0 - getMiddle(posR, posL);
    changePosition(posR, oppMiddle);
    changePosition(posL, oppMiddle);
    
    if(yR == yL)
    {
        cerr << "ERROR in ::lin_Y_to_Pos: yR = yL" << endl;
        exit(EXIT_FAILURE);
    }
    
    double p = posR + (posL - posR) * (y - yR) / (yL - yR);
    
    changePosition(p, -oppMiddle);
    
    return p;
}
// --------------------------------------------------
// For writing outputs: 
// --------------------
string doubleToString(double value, unsigned int precision)
{
    ostringstream oss;
    oss << setprecision(precision) << fixed << value;
    
    return oss.str();
}
// --------------------
string doubleToString(double value, unsigned int precision,
                      string prefix, string suffix)
{
    ostringstream oss;
    oss << prefix
    << setprecision(precision) << fixed << value
    << suffix;
    
    return oss.str();
}
// --------------------
Direction doubleToDirection(double value, double epsilon)
{
    Direction direction = NONE;
    if(value > epsilon)
        direction = LEFT;
    else if(value < -epsilon)
        direction = RIGHT;
    
    return direction;
}
// --------------------
string idSuffixe(string id)
{
    string number = "";
    unsigned int pos = id.find(".");
    if(pos < 1000)
        number = id.erase(0, pos+1);
    return number;
}

// conversions ---------------------------
double yr_to_sec(double tyr)
{
    return tyr * seconds_per_yrs;
}

double Myr_to_sec(double tMyr)
{
    return tMyr * seconds_per_Myrs;
}

double sec_to_Myr(double tsec)
{
    return tsec / seconds_per_Myrs;
}
// ----------------------------------------

/*
 =======================================
 Initialization of _earth static members
 =======================================
 */
double Earth::_time = 0.0;
unsigned int Earth::_timestep = 0;
unsigned int Earth::_instanceCounter = 0;
unsigned int Earth::_ridgeCounter = 0;
unsigned int Earth::_subductionCounter = 0;
unsigned int Earth::_stapleCounter = 0;
unsigned int Earth::_plateCounter = 0;
unsigned int Earth::_plateSectionCounter = 0;
unsigned int Earth::_cellCounter = 0;
unsigned int Earth::_continentCounter = 0;
unsigned int Earth::_continentExtremityCounter = 0;
unsigned int Earth::_warmingZoneCounter = 0;
unsigned int Earth::_collisionCounter = 0;
unsigned int Earth::_activeMarginCounter = 0;

/*
 =========================================
 Methods for class Earth
 =========================================
 */

Earth::Earth(double dt) : _dt(dt)
{
    cout << "... Starting a new Earth" << endl;
    
    _random = NULL;
    
    // --- Properties
    _ready = false;
    _running = false;
    _runStarted = false;
    
    _ageSlice = 5.0; // in Myr, for distribution of ages.
    
    // --- Graphics & outputs
    _needRefresh = false;
    _drawCells = true;
    _showEvents = true;
    _drawAgeOceans = false;
    
    _nextTimeWriteLogs = 0.0;
    _nextTimeWriteAges = 0.0;
    
    // --- Initialization
    prepareFiles();
    
    defaultParameters();
    initFixedParameters();
    
    initParameters(); // time, physical and model
    
    // -------------------------------------------------------
    // some booleans for special testing and/or writings
    _computeEta = true;
    _write_continents = true; // aside from continents, all should be removed
    _write_plate_forces = false; // true for debugging mainly
    // should be let to false to have plates files ready to draw
    // -------------------------------------------------------
    
    // OTHERS:
    _depthFactor = 0.1;
    _insertAges = true; /* gives better precision on heat flux
                         and avoids jumps */
    
    _noSubduction = false;
    _timespan_noSubduction = 0.0;
    _inactiveRidges = false;
    
    _T_cond_reached = false;
    _age_cond_reached = false;
    
    saveState();
}
// ----------------------------------------
Earth::~Earth()
{
    clear();
    
    if(_logFile.is_open())
        _logFile.close();
    if(_eventsFile.is_open())
        _eventsFile.close();
    if(_velocityFile.is_open())
        _velocityFile.close();
    if(_geometryFile.is_open())
        _geometryFile.close();
    if(_outputFile.is_open())
        _outputFile.close();
    if(_continentsFile.is_open())
        _continentsFile.close();
    if(_ageFile.is_open())
        _ageFile.close();
    if(_forcesFile.is_open())
        _forcesFile.close();
    if(_radioContFile.is_open())
        _radioContFile.close();
    if(_testFile.is_open())
        _testFile.close();
    if(_bathymetryFile.is_open())
        _bathymetryFile.close();
    if(_platesFile.is_open())
        _platesFile.close();
    if(_elementsFile.is_open())
        _elementsFile.close();
    if (_collisionsFile.is_open())
        _collisionsFile.close();
}

// ----------------------------------------
void
Earth::clear()
{
    lock();
    
    // --- Plates
    for(unsigned int i=0; i<_plates.size(); i++)
        delete _plates[i];
    _plates.clear();
    
    // --- Continents
    for(unsigned int i=0; i<_continents.size(); i++)
        delete _continents[i];
    _continents.clear();
    
    // --- Ages
    _ages.clear();
    
    // --- GeoElements
    _elements.clear();
    
    // --- Static
    _instanceCounter = 0;
    _ridgeCounter = 0;
    _subductionCounter = 0;
    _stapleCounter = 0;
    _plateCounter = 0;
    _plateSectionCounter = 0;
    _cellCounter = 0;
    _continentCounter = 0;
    _continentExtremityCounter = 0;
    _warmingZoneCounter = 0;
    _collisionCounter = 0;
    _activeMarginCounter = 0;
    
    _time = 0.0;
    _timestep = 0;
    
    // --- Properties
    _running = false;
    _ready = false;
    
    _T = 1600.0;
    
    unlock();
}
// ----------------------------------------
void
Earth::clearPlates()
{
    _needRefresh = true;
    
    _plateCounter = 0;
    _plateSectionCounter = 0;
    _cellCounter = 0;
    
    for(unsigned int i=0; i<_elements.size(); i++)
    {
        _elements[i]->setLeftSection(NULL);
        _elements[i]->setRightSection(NULL);
        
        if(_elements[i]->isA("Interface"))
        {
            ((Interface*)_elements[i])->setLeftPlate(NULL);
            ((Interface*)_elements[i])->setRightPlate(NULL);
        }
    }
    
    for(unsigned int i=0; i<_continents.size(); i++)
        _continents[i]->setPlateSection(NULL);
    
    for(unsigned int i=0; i<_plates.size(); i++)
        delete _plates[i];
    _plates.clear();
}

/* ==============================
 INITIALIZATION
 ============================== */
void 
Earth::prepareFiles()
{
    string filename = _workspace_ + "logs/earth.log";
    _logFile.open(filename.c_str(), ios::out | ios::trunc);
    
    if(!_logFile.is_open())
        cerr << "Error while opening " << filename << endl;
    else
    {
        _logFile << "# 1-Age(Ma) \t"
        << "2-T(K) \t"
        << "3-Qtot(TW) \t"
        << "4-Qoc(mW/m2) \t"
        << "5-Qcont(mW/m2) \t"
        << "6-H(TW)\t"
        << "7-Urey\t"
        << "8-Urms(cm/yr)\t"
        << "9-UrmsL \t"
        << "10-Nb_plates \t"
        << "11-Nb_continents \t"
        << "12-Age_max(Myr) \t"
        << "13-Age_mean \t"
        << "14-Age_SD \n";
    }
    
    // --- Events
    filename = _workspace_ + "logs/events.log";
    _eventsFile.open(filename.c_str());
    if(!_eventsFile.is_open())
        cerr << "Error while opening " << filename << endl;
    else
    {
        _eventsFile << "# 1-Age(Ma) \t"
        << "2-T(K) \t"
        << "3-Type \n";
    }
    
    // --- velocities output
    filename = _workspace_ + "logs/velocity.log";
    _velocityFile.open(filename.c_str(), ios::out | ios::trunc);
    if(!_velocityFile.is_open())
        cerr << "Error while opening" << filename << endl;
    else
    {
        _velocityFile << "# 1-Age(Ma) \t"
        << "2-T(K) \t"
        << "3-Uabs(cm/yr) \t"
        << "4-UabsL \t"
        << "5-Umean \t"
        << "6-UmeanL \t"
        << "7-Umin \t"
        << "8-Umax \t"
        << "9-UElements \t"
        << "10-URidges \t"
        << "11-USubductions \t"
        << "12-UrmsCont \t"
        << "13-UrmsContL \t"
        << "14-UrmsContLpl \t"
        << "15-UabsCont \t"
        << "16-UabsContL \t"
        << "17-UabsContLpl \t"
        << "18-UOceans \t"
        << "19-UOceansL \t"
        << "20-Usubd \t"
        << "21-Uupper \t"
        << "22-UsubdL \t"
        << "23-UupperL \n" ;
    }
    
    // --- ages output
    filename = _workspace_ + "logs/age.log";
    _ageFile.open(filename.c_str(), ios::out | ios::trunc);
    if(!_ageFile.is_open())
        cerr << "Error while opening " << filename << endl;
    else
    {
        _ageFile << "# 1-Age(Ma) \t"
        << "2-T(K) \t"
        << "3-Age_min(Myr) \t"
        << "4-D1 \t"
        << "5-Q1 \t"
        << "6-Q2 \t"
        << "7-Q3 \t"
        << "8-D9 \t"
        << "9-Age_max \t"
        << "10-maxAgeSubduction \t"
        << "11-maxAgeUpperPlate \n";
    }
    
    // --- forces output
    filename = _workspace_ + "logs/forces.log";
    _forcesFile.open(filename.c_str(), ios::out | ios::trunc);
    if(!_forcesFile.is_open())
        cerr << "Error while opening " << filename << endl;
    else
    {
        _forcesFile << "# 1-Age(Ma) \t"
        << "2-T(K) \t"
        << "3-nb_ridges \t"
        << "4-nb_subductions \t"
        << "5-SP(TN/m) \t"
        << "6-abs(SP) \t"
        << "7-RP \t"
        << "8-abs(RP) \t"
        << "9-SS \t"
        << "10-abs(SS) \t"
        << "11-MD \t"
        << "12-abs(MD) \t"
        << "13-VS \t"
        << "14-abs(VS) \t"
        << "15-B \t"
        << "16-abs(B) \t"
        << "17-meanPlateThickness(km)\n" ;
    }
    
    // plates and continents length
    filename = _workspace_ + "logs/geometry.log";
    _geometryFile.open(filename.c_str(), ios::out | ios::trunc);
    if(!_geometryFile.is_open())
        cerr << "Error while opening " << filename << endl;
    else
    {
        _geometryFile << "# 1-Age(Ma) \t"
        << "2-T(K) \t"
        << "3-Rms_PlateLength (km)\t"
        << "4-Min_PlateLength \t"
        << "5-Max_PlateLength \t"
        << "6-Rms_ContinentLength \t"
        << "7-Min_ContinentLength \t"
        << "8-Max_ContinentLength \t"
        << "9-seafloor_prod(km^2/yr) \t"
        << "10-slab_flux(km^2/yr) \t"
        << "11-slab_fluxL(km^2/yr) \t"
        << "12-NContPieces \n";
    }
    
    // checking radioactive heating
    filename = _workspace_ + "logs/radioHeat_continent.log";
    _radioContFile.open(filename.c_str());
    if(!_radioContFile.is_open())
        cerr << "Error while opening " << filename << endl;
    else
    {
        _radioContFile << "# 1-Age(Myr) \t"
        << "2-T(K) \t"
        << "3-Ocean(km^2) \t"
        << "4-Continent(km^2) \t"
        << "5-surf_cont_ratio \t"
        << "6-H(TW) \t"
        << "7-H_primitive(TW) \t"
        << "8-H_depleted(TW) \n";
    }
    
    // TEST on density computations ////////
    filename = _workspace_ + "logs/rho.log";
    _testFile.open(filename.c_str(), ios::out | ios::trunc);
    if(_testFile.is_open())
    {
        _testFile << "# 1-Age(Ma) \t"
        << "2-T(K) \t"
        << "3-T_plate \t"
        << "4-rho_um \t"
        << "5-rho_pl \t"
        << endl;
    }
    // END TEST //////////////
    
    // sea level
    filename = _workspace_ + "logs/bathymetry.log";
    _bathymetryFile.open(filename.c_str(), ios::out | ios::trunc);
    if(_bathymetryFile.is_open())
    {
        _bathymetryFile << "# 1-Age(Ma) \t"
        << "2-T(K) \t"
        << "3-VbelowRidge_HSCM(m^3) \t"
        << "4-ridgeDepth_HSCM(m) \t"
        << "5-sealevel_HSCM(m)\t"
        << "6-VbelowRidge_PM95(m^3) \t"
        << "7-ridgeDepth_PM95(m) \t"
        << "8-sealevel_PM95(m)\t"
        << "9-VbelowRidge_PM125(m^3) \t"
        << "10-ridgeDepth_PM125(m) \t"
        << "11-sealevel_PM125(m)"
        << endl;
    }
    
    filename = _workspace_ + "logs/continents.log";
    _continentsFile.open(filename.c_str(), ios::out | ios::trunc);
    if(_continentsFile.is_open())
    {
        _continentsFile << "# 1-Age(Ma) \t"
        << "2-NCont\t"
        << "3-Id\t"
        << "4-Length(deg)\t"
        << "5-Pos(deg)\t"
        << "6-V(cm/yr)\t"
        << "7-Age(Ma)\t"
        << "8-Twarm(K)\t"
        << "9-F(TN)\t ..."
        << endl;
    }
    
    filename = _workspace_ + "logs/plates.log";
    _platesFile.open(filename.c_str(), ios::out | ios::trunc);
    if(_platesFile.is_open())
    {
        _platesFile << "# 1-Age(Ma) \t"
        << "2-NPlates\t"
        << "3-RightElement\t"
        << "4-LeftElement\t"
        << "5-V(cm/yr)\t"
        << "6-Length(deg)\t"
        << "7-LengthCont\t"
        << "8-Ncont\t...\n";
    }
    
    filename = _workspace_ + "logs/elements.log";
    _elementsFile.open(filename.c_str(), ios::out | ios::trunc);
    if(_elementsFile.is_open())
    {
        _elementsFile << "# 1-Age(Ma) \t"
        << "2-NElements\t"
        << "3-Id\t"
        << "4-position(deg)\t...\n";
    }
    
    filename = _workspace_	 + "logs/collisions.log";
    _collisionsFile.open(filename.c_str(), ios::out | ios::trunc);
    if (_collisionsFile.is_open())
    {
        _collisionsFile << "# 1-Age(Ma) \t"
        << "2-NCollisions\t"
        << "3-Id\t"
        << "4-position(deg)\t"
        << "5-age(Ma)\t...\n";
    }

    filename = _workspace_     + "logs/activeMargins.log";
    _activeMarginsFile.open(filename.c_str(), ios::out | ios::trunc);
    if (_activeMarginsFile.is_open())
    {
        _activeMarginsFile << "# 1-Age(Ma) \t"
        << "2-NActiveMargins\t"
        << "3-Id\t"
		<< "4-Status\t"
        << "5-position(deg)\t"
        << "6-age(Ma)\t...\n";
    }
    
    filename = _workspace_ + "logs/run.output";
    _outputFile.open(filename.c_str());
    
    if(!_outputFile.is_open())
        cerr << "Error while opening " << filename << endl;
    
} // prepareFiles

// ----------------------------------------
void
Earth::defaultParameters()
{
    // --- Fixed parameters: not to change through the graphic interface
    _fixedParameters["earthRadius"]     =      6371.0e3; // m
    _fixedParameters["g"]    	 =      10.0;
    
    _fixedParameters["d"]     	 = 	2900e3; // m
    _fixedParameters["D"]     	 =	670e3;
    
    _fixedParameters["T_surf"] 	 =	300.0;
    
    _fixedParameters["kappa"] 	 =	8e-7;
    _fixedParameters["C_p"]   	 =	1200.0;
    
    _fixedParameters["earthMass"]  =      5.9736e24;
    _fixedParameters["mantleMass"] =      4.079e24;
    
    // for radioactive decay and heat
    _fixedParameters["halfLife235U"] = 704.0; // Myr
    _fixedParameters["halfLife238U"] = 4470.0; // Myr
    _fixedParameters["halfLife232Th"] = 14.0E3; // Myr
    _fixedParameters["halfLife40K"] = 1.25E3; // Myr
    _fixedParameters["heat235U"] = 5.69E-4; // W/kg
    _fixedParameters["heat238U"] = 9.46E-5; // W/kg
    _fixedParameters["heat232Th"] = 2.64E-5; // W/kg
    _fixedParameters["heat40K"] = 2.92E-5; // W/kg
    _fixedParameters["235U_over_U"] = 0.0072; // ratio 235U over total U
    _fixedParameters["238U_over_U"] = 0.9927; // ratio 238U over total U
    _fixedParameters["40K_over_K"] = 1.28E-4; // ratio 40K over total K
    _fixedParameters["present_cont_surfRatio"] = 0.40;
    _fixedParameters["rho_cont"] = 2700;
    _fixedParameters["cont_thickness"] = 40.0E3; // in mA
    
    // computation of rho differences
    _fixedParameters["slabPull_rho_depends_on_T"]  = 1.0; // 1.0 for true
    _fixedParameters["ridgePush_rho_depends_on_T"] = 1.0; // 1.0 for true
    
    // for sea level
    _fixedParameters["rho_seawater"] = 1030.0;
    _fixedParameters["present_oceanVolume"] = 1.31750E18; // m^3
    _fixedParameters["present_ridgeDepth"] = 2446.0; // m
    _fixedParameters["zm125"] = 125.0E3; // in m
    _fixedParameters["zm95"] = 95.0E3; // in m
    
    // -- other parameters: appear in graphic interface
    
    // Time parameters
    _timeParameters["startAge"]  = 3000.0;
    _timeParameters["endAge"]    = -100.0;
    _timeParameters["startingTime"] = -1.0; // to fit former versions
    _timeParameters["runningTime"] = -1.0; // to fit former versions
    _timeParameters["T_final"]         = 1450.0;
    _timeParameters["writeAges"]       = 5.0; // in Myr
    _timeParameters["writeLogs"]       = 0.5;
    _timeParameters["writeXml"]        = 100.0;
    _timeParameters["writeImages"]     = 5.0; //
    _timeParameters["resolution"]      = 0.5; // in degrees
    _timeParameters["courantNumber"]   = 0.50;
    _timeParameters["timeStep"]        = 10000.0; // in years
    _timeParameters["fixed_timestep"]  = 0; // non zero if fixed
    _timeParameters["age_T_and_condition"] = 0.0; // non zero if AND, 0.0 if OR condition
    
    // Physical parameters
    _physicalParameters["T_m_init"]      =	1850.0;
    _physicalParameters["T_p"]  	          =	1625.0;
    _physicalParameters["Tau_sub_p"]     =        180.0;
    _physicalParameters["Tau_ssc_p"]     =        80.0;
    _physicalParameters["V_sink_p"]      =        1.0; // in cm/yr
    _physicalParameters["R_min"] 	       =	390; // in km
    _physicalParameters["F_lim"] 	       =	3.0E12; // in N/m
    _physicalParameters["Qmax"]          =        0.000;  // en W/m^2
    _physicalParameters["min_plate_thick"]   =    0.000; // in km
    _physicalParameters["subcont_warming_H"] =    350; // in km
    _physicalParameters["E_m"]    	  =	300; // in kJ/mol
    _physicalParameters["E_um"] 	          =	300;
    _physicalParameters["E_pl"]    	  =	300;
    _physicalParameters["eta_m_p"]     	  =     1e22;
    _physicalParameters["eta_um_p"] 	  =     1e21;
    _physicalParameters["eta_ast_p"] 	  =     1e21;
    _physicalParameters["eta_pl_p"] 	  =     1e23;
    _physicalParameters["eta_subcont_p"]    =     1e21;
    _physicalParameters["thick_ast"]        =     0.0; // in km
    _physicalParameters["thick_subcont"]    =     0.0;
    _physicalParameters["thick_continent"]  =     200.0;
    _physicalParameters["k_continent"]      =     0.0;
    _physicalParameters["k_ocean"]          =     3.0;
    _physicalParameters["rho_um_p"]        =     3300.0; // in kg.m^-3
    _physicalParameters["DeltaRho_p"]      =     65.0; // in kg.m^-3
    _physicalParameters["alpha_um"]          =     2.0e-5; // in K^-1
    _physicalParameters["alpha_pl"]          =     3.25e-5; // in K^-1
    
    // Model parameters
    _modelParameters["coeffSlabPull"]     = 1.0; // all booleans below
    _modelParameters["coeffRidgePush"]    = 1.0;
    _modelParameters["coeffSlabSuction"]  = 1.0;
    _modelParameters["coeffMantleDrag"]   = 1.0;
    _modelParameters["coeffViscousShear"] = 1.0;
    _modelParameters["coeffBending"]      = 1.0;
    _modelParameters["slabPullWholeMantle"] = 0.0;
    _modelParameters["viscousShearWholeMantle"] = 0.0;
    _modelParameters["initMode_brittle"] = 1.0;
    _modelParameters["initMode_constant"] = 0.0;
    _modelParameters["initMode_convective"] = 0.0;
    _modelParameters["initPlace_continents"] = 1.0;
    _modelParameters["initPlace_staples"] = 1.0;
    _modelParameters["initPlace_upperPlates"] = 0.0;
    _modelParameters["upperPlates_ageRatioCriterion"] = 0.0; // all booleans above
    _modelParameters["upperPlates_ratioAgeValue"] = 3.0; // dimensionless value
    _modelParameters["upperPlates_minimumAge"] = 50.0; // Myr
    _modelParameters["reverseSubduction"] = 0.0; // bool
    _modelParameters["tauSub_randomNoise"] = 0.0; // in Myr
    _modelParameters["alwaysSubduct"] = 0.0; //bool
    _modelParameters["small_scale_convection"] = 0.0; // bool
    _modelParameters["sscMode_constant"] = 0.0; // bool
    _modelParameters["sscMode_convective"] = 1.0; // bool
    _modelParameters["middle_breakup"] = 1.0; // bool
    _modelParameters["breakup_position"] = 0.333;
    _modelParameters["fixed_configuration"] = 0.0; // bool
    _modelParameters["U_BSE"] = 0.020; // ppm from McDonough and Sun (1995)
    _modelParameters["Th_BSE"] = 0.079; // ppm from McDonough and Sun (1995)
    _modelParameters["K_BSE"] = 240.0; // ppm from McDonough and Sun (1995)
    _modelParameters["U_cont"] = 1.30; // ppm from Rudnick and Gao (2003)
    _modelParameters["Th_cont"] = 5.60; // ppm from Rudnick and Gao (2003)
    _modelParameters["K_cont"] = 1.50E4; // ppm from Rudnick and Gao (2003)
    _modelParameters["depletion_always"] = 1.0; // 1.0 for true
    _modelParameters["depletion_never"] = 0.0; // 1.0 for true
    _modelParameters["depletion_adaptative"] = 0.0; // 1.0 for true.
    _modelParameters["continental_growth"] = 0.0; // 0.0 for false, 1.0 for true
    _modelParameters["contGrowthCoeff"]    =  1.0; // no dimension
    _modelParameters["randomSeed"] = 1;
    _modelParameters["constant_internal_heating"] = 0;
    _modelParameters["constant_heat_value"] = 40.0;// in TW
    
    
} // default parameter
// ----------------------------------------
void
Earth::initParameters()
{ /* initParameters is first done in Earth constructor.
   initFixedParameters() is done before,
   and called only once in the constructor
   */
    
    initTimeParameters();
    initPhysicalParameters();
    initModelParameters();
    
    // some computations that can be done initially:
    _T_plate_p = compute_T_plate(_T_p);
    _rho_pl_p = _rho_um_p + _DeltaRho_p;
    
    _preFactor_ridgePush = _coeffRidgePush * _alpha_pl * _g * _kappa;
    _preFactor_heatFlow = _k_ocean / sqrt( M_PI * _kappa);
    _preFactor_bathyHSCM = 2.0 * _alpha_pl * sqrt(_kappa / M_PI);
    _preFactor_bathyPM2 = _kappa * seconds_per_Myrs * pow(M_PI, 2);
}
// --------------------------------------------------
void
Earth::initFixedParameters()
{
    _T_surf = getFixedParameter("T_surf");
    
    _d = getFixedParameter("d");
    _D = getFixedParameter("D");
    
    _g = getFixedParameter("g");
    
    _kappa = getFixedParameter("kappa");
    _C_p = getFixedParameter("C_p");
    
    _earthRadius = getFixedParameter("earthRadius");
    _earthMass = getFixedParameter("earthMass");
    _present_mantleMass = getFixedParameter("mantleMass");
    
    _heat235U = getFixedParameter("heat235U");
    _heat238U = getFixedParameter("heat238U");
    _heat232Th = getFixedParameter("heat232Th");
    _heat40K = getFixedParameter("heat40K");
    _rho_cont = getFixedParameter("rho_cont");
    _cont_thickness = getFixedParameter("cont_thickness");
    _present_cont_surfRatio = getFixedParameter("present_cont_surfRatio");
    
    _slabPull_rho_depends_on_T  =
    (getFixedParameter("slabPull_rho_depends_on_T") == 1.0);
    _ridgePush_rho_depends_on_T =
    (getFixedParameter("ridgePush_rho_depends_on_T") == 1.0);
    
    _rho_seawater = getFixedParameter("rho_seawater");
    _present_oceanVolume = getFixedParameter("present_oceanVolume");
    _present_ridgeDepth = getFixedParameter("present_ridgeDepth");
    _zm95 = getFixedParameter("zm95");
    _zm125 = getFixedParameter("zm125");
    
    // computations that can be done already with fixed parameters:
    _surfTot = 4.0 * M_PI * pow(_earthRadius, 2);
    
    _present_continentMass =  _present_cont_surfRatio *
    _surfTot * _rho_cont * _cont_thickness;
    
    _primitive_mantleMass =
    _present_mantleMass + _present_continentMass;
    
    _present_VaboveRidge = _present_ridgeDepth * 4.0 * M_PI *
    (1.0 - _present_cont_surfRatio) * pow( _earthRadius, 2.0 );
    
}
// -----------------------------------------
void
Earth::formerTimeSetting()
{
    /* check if former version of time are given
     in config file.
     */
    if(getTimeParameter("startingTime") > 0.0)
    {
        double startAge = Earth::ageEarth - getTimeParameter("startingTime");
        setTimeParameter("startAge", startAge);
        
        cout << "** Warning: startingTime is given in input configuration file\n";
        cout << "\t will not be supported in future versions -> please switch to startAge/endAge\n";
        cout << "\t startingTime = " << getTimeParameter("startingTime")
        << "\t --> startAge = " << getTimeParameter("startAge") << endl;
        
        setTimeParameter("startingTime", -1.0);
    }
    
    
    if(getTimeParameter("runningTime") > 0.0)
    {
        double endAge = getTimeParameter("startAge") - getTimeParameter("runningTime");
        setTimeParameter("endAge", endAge);
        
        cout << "** Warning: runningTime is given in input configuration file\n";
        cout << "\t will not be supported in future versions -> please switch to startAge/endAge\n";
        cout << "\t runningTime = " << getTimeParameter("runningTime")
        << "\t --> endAge = " << getTimeParameter("endAge") << endl;
        
        setTimeParameter("runningTime", -1.0);
    }
    
}
// --------------------------------------------------
void
Earth::initTimeParameters()
{
    
    _startAge = getTimeParameter("startAge");
    _endAge = getTimeParameter("endAge");
    
    if(!_runStarted)
        _ageMa = _startAge;
    
    _T_final = getTimeParameter("T_final");
    _timeStep = getTimeParameter("timeStep");
    _writeLogs = getTimeParameter("writeLogs");
    _writeAges = getTimeParameter("writeAges");
    _writeImages = getTimeParameter("writeImages");
    _writeXml = getTimeParameter("writeXml");
    _resolution = getTimeParameter("resolution");
    _courant = getTimeParameter("courantNumber");
    _fixed_timestep = (getTimeParameter("fixed_timestep") == 0 ? false : true);
    _age_T_and_condition = (getTimeParameter("age_T_and_condition") == 1.0 ? true : false);
    if(_fixed_timestep)
        _dt = _timeStep;
    else
        _dt = 200.0;
    
    
    // make sure the age resolution is smaller than creationDistance
    if(_resolution <= 2.0 * Earth::creationDistance)
    {
        cout << "WARNING: _resolution " << _resolution
        << " is smaller than 2*creationDistance " << 2.0*Earth::creationDistance
        << "\n \t -> lowering resolution to " << Earth::creationDistance / 2.0
        << endl;
        _resolution = Earth::creationDistance / 2.0;
    }
    
}
// --------------------------------------------------
void
Earth::initPhysicalParameters()
{
    _T_m_init = getPhysicalParameter("T_m_init");
    _T_p =  getPhysicalParameter("T_p");
    _tau_sub_p = getPhysicalParameter("Tau_sub_p"); // in Myr
    _tau_ssc_p = getPhysicalParameter("Tau_ssc_p");
    _V_sink_p = getPhysicalParameter("V_sink_p");
    _F_lim = getPhysicalParameter("F_lim");
    _subcont_warming_H = getPhysicalParameter("subcont_warming_H") * 1.0E3; // from km to m
    _R_min =  getPhysicalParameter("R_min") * 1.0E3; // from km to m
    _Qmax =  getPhysicalParameter("Qmax");
    _min_plate_thick = getPhysicalParameter("min_plate_thick") * 1.0E3; // from km to m
    _eta_m_p =  getPhysicalParameter("eta_m_p") / Earth::eta0; //non-dimensional eta
    _eta_um_p =  getPhysicalParameter("eta_um_p") / Earth::eta0;
    _eta_pl_p =  getPhysicalParameter("eta_pl_p") / Earth::eta0;
    _eta_ast_p =  getPhysicalParameter("eta_ast_p") / Earth::eta0;
    _eta_subcont_p =  getPhysicalParameter("eta_subcont_p") / Earth::eta0;
    _E_m = getPhysicalParameter("E_m") * 1.0E3; // from kJ/mol to J/mol
    _E_um = getPhysicalParameter("E_um") * 1.0E3;
    _E_pl = getPhysicalParameter("E_pl") * 1.0E3;
    _thick_ast = getPhysicalParameter("thick_ast") * 1.0E3; // from km to m
    _thick_subcont = getPhysicalParameter("thick_subcont") * 1.0E3;
    _thick_continent = getPhysicalParameter("thick_continent") * 1.0E3;
    _k_ocean = getPhysicalParameter("k_ocean");
    _k_continent = getPhysicalParameter("k_continent");
    _rho_um_p = getPhysicalParameter("rho_um_p");
    _DeltaRho_p = getPhysicalParameter("DeltaRho_p");
    _alpha_um = getPhysicalParameter("alpha_um");
    _alpha_pl = getPhysicalParameter("alpha_pl");
    
    if(!_runStarted)
        _T = _T_m_init;
    
}
// --------------------------------------------------
void
Earth::initModelParameters()
{
    _coeffSlabPull = getModelParameter("coeffSlabPull");
    _coeffRidgePush = getModelParameter("coeffRidgePush");
    _coeffSlabSuction = getModelParameter("coeffSlabSuction");
    _coeffMantleDrag = getModelParameter("coeffMantleDrag");
    _coeffViscousShear = getModelParameter("coeffViscousShear");
    _coeffBending = getModelParameter("coeffBending");
    
    _slabPullDepth =
    (getModelParameter("slabPullWholeMantle") == 0 ? UPPER_MANTLE : WHOLE_MANTLE);
    _viscousShearDepth =
    (getModelParameter("viscousShearWholeMantle") == 0 ? UPPER_MANTLE : WHOLE_MANTLE);
    
    // subduction
    setSubductionMode(); // using _modelParameters("initMode_XXX")
    setSubductionPlace(); // using _modelParameters("initPlace_XXX")
    
    // SSC or not SSC
    setLimitedThickening();
    setSSCMode();
    
    // special case
    setFixedConfiguration();
    
    // continental breakup
    _middleBreakup = (getModelParameter("middle_breakup") != 0.0);
    _breakupPosition = getModelParameter("breakup_position");
    
    // radio heating
    setDepletionMode();
    _constant_heating = (getModelParameter("constant_internal_heating") == 0.0 ? false : true);
    _constant_heat_value = getModelParameter("constant_heat_value");
    
    // continental growth
    _continental_growth = (getModelParameter("continental_growth") == 0.0 ? false : true);
    _contGrowthCoeff = getModelParameter("contGrowthCoeff");
    
    // randomSeed used for continental breakup and age for subduction
    _randomSeed = (unsigned int)(getModelParameter("randomSeed"));
    if(_random)
        _random = NULL;
    
    _random = new Random(_randomSeed);
}
// --------------------------------------------------
void
Earth::makeReady()
{ 
    /* make sure the structure is set,
     with ages and slabs initialized,
     as well as all values depending on initial parameters.
     Also write initial parameters in _outputFile. */
    
    //  checkHeatFlux();
    checkFixedConfiguration();
    
    _T = _T_m_init;
    initTData();
    
    initAges();
    initSlabs();
    
    updateCells();
    
    _age_min_thickness = thickness_to_age(_min_plate_thick);
    
    ostringstream oss;
    
    // - - - - - - - - - - - - - -
    oss << "\n#########################" << endl;
    oss << "age_min_thickness = " << _age_min_thickness
    << " Myr" << endl;
    
    // - - - - - - - - - - - - - -
    _earthPerimeter = 2.0 * M_PI * _earthRadius; // in m
    
    updateSurfaces(); // ocean and continental surfaces
    updateMasses();
    
    
    oss << "#########################" << endl;
    oss << "Oceanic total surface = " << _surfTotOcean << " m^2" << endl;
    oss << "Oceanic percent = " << _surfTotOcean / _surfTot * 100.0 << " %" << endl;
    oss << "Continental total surface = " << _surfTotContinent << " m^2" << endl;
    oss << "Continental percent = " << _cont_surfRatio * 100.0 << " %" << endl;
    oss << "Present day ocean volume above ridge = "
    << _present_VaboveRidge << " m^3" << endl;
    // - - - - - - - - - - - - - -
    computeThicknessesMantle(); /* thicknesses used in
                                 Plate::computeHorizontalViscosity() */
    oss << "#########################" << endl;
    oss << "LAYERS (in km):" << endl;
    oss << "\t\t SUBCONT  |    OCEAN " << endl;
    oss << " subcont. layer =  " << setw(5) << _thicknesses.subcont*1.0E-3
    << "\t | " << setw(5) << 0 << endl;
    oss << " asthenosphere  =  " << setw(5) << _thicknesses.contAst*1.0E-3
    << "\t | "  << setw(5) << _thicknesses.oceanAst*1.0E-3 << endl;
    oss << " upper mantle   =  " << setw(5) << _thicknesses.contUM*1.0E-3
    << "\t | "  << setw(5) << _thicknesses.oceanUM*1.0E-3 << endl;
    oss << " lower mantle   =       " << setw(5) << _thicknesses.LM*1.0E-3 << endl;
    
    // for color of continents and oceans sections in GUI - - - - - - -
    setContinentalAgeSlices();
    if(_drawAgeOceans)
        setOceanAgeSlices();
    
    // PARAMETERS ----------------------
    oss << "\n#########################\n";
    
    map<string,double>::iterator it;
    // fixed
    oss << "FIXED PARAMETERS:" << endl;
    for(it=_fixedParameters.begin(); it!=_fixedParameters.end(); it++)
        oss << "\t" << it->first << " = " << it->second << endl;
    // time
    oss << "\nTIME PARAMETERS:" << endl;
    for(it=_timeParameters.begin(); it!=_timeParameters.end(); it++)
        oss << "\t" << it->first << " = " << it->second << endl;
    // physics
    oss << "\nPHYSICAL PARAMETERS:" << endl;
    for(it=_physicalParameters.begin(); it!=_physicalParameters.end(); it++)
        oss << "\t" << it->first << " = " << it->second << endl;
    // model
    oss << "\nMODEL PARAMETERS:" << endl;
    for(it=_modelParameters.begin(); it!=_modelParameters.end(); it++)
        oss << "\t" << it->first << " = " << it->second << endl;
    
    oss << "#########################\n\n";
    
    
    // ELEMENTS -------------------------
    for(unsigned int i=0; i<_elements.size(); i++)
    {
        oss << setw(18) << _elements[i]->getId()
        << setw(14) << " Position = " << setw(6) << _elements[i]->getPosition()
        << setw(15) << " Left age  = " << setw(6) << _elements[i]->getLeftAge().getValue() << setw(3) << " Ma"
        << setw(15) << " Right age = " << setw(6) << _elements[i]->getRightAge().getValue() << setw(3) << " Ma" ;
        if(_elements[i]->isA("Subduction"))
        {
            Subduction* subduction = (Subduction*)_elements[i];
            oss << setw(12) << " Depth = "
            << setw(6) << subduction->getDepth()/1.0E3 << setw(3) << " km";
        }
        oss << endl;
    }
    
    // CONTINENTS -----------------------------------
    if(_continents.size() > 0)
    {
        for(unsigned int i=0; i<_continents.size(); i++)
        {
            double leftage = _continents[i]->getLeftExtremity()->getLeftAge().getValue();
            double rightage = _continents[i]->getRightExtremity()->getRightAge().getValue();
            oss << "\n" << setw(18) << _continents[i]->getId()
            << setw(14) << " Position = "  << setw(6) << _continents[i]->getRightPosition()
            << setw(15) << " Left age  = " << setw(6) << leftage << setw(3) << " Ma"
            << setw(15) << " Right age = " << setw(6) << rightage << setw(3) << " Ma"
            << setw(15) << " Length = " << setw(4) << _continents[i]->getLength() << setw(4) << " deg";
        }
    }
    oss << "\n#########################\n";
    
    // radioactive heating - - - - -
    initConcentrations();
    
    oss << "fixed depletion ratio: " << boolalpha << _fixed_depletion << endl;
    oss << "depletion ratio = " << _depletionRatio << endl;
    oss << "#########################\n";
    
    // Warnings about some tests
    if( _fixed_configuration || !_computeEta )
    {
        oss << "                 WARNING                " << endl;
        oss << "                 -------                " << endl;
        if(_fixed_configuration)
            oss << " _fixed_configuration is true " << endl;
        if(!_computeEta)
            oss << " _computeEta is false: FIXED VISCOSITIES!" << endl;
        oss << "###########################################" << endl;
    }
    
    cout << oss.str();
    if(_outputFile.is_open())
        _outputFile << oss.str() << flush;
    
    _ready = true;
    
} // makeReady()


/*
 =======================
 Static
 =======================
 */

double
Earth::getTime()
{
    return _time; // in years
}

double
Earth::getTimeMyr()
{
    return _time / 1.0e6;
}

int
Earth::getTimestep()
{
    return _timestep;
}

string
Earth::getId(string className)
{
    string id = className + ".";
    ostringstream oss;
    
    if(className == "Ridge")
    {
        oss << _ridgeCounter;
        _ridgeCounter++;
    }
    else if(className == "RightSubduction" || className == "LeftSubduction")
    {
        oss << _subductionCounter;
        _subductionCounter++;
    }
    else if(className == "Staple")
    {
        oss << _stapleCounter;
        _stapleCounter++;
    }
    else if(className == "Plate")
    {
        oss << _plateCounter;
        _plateCounter++;
    }
    else if(className == "PlateSection")
    {
        oss << _plateSectionCounter;
        _plateSectionCounter++;
    }
    else if(className == "Cell")
    {
        oss << _cellCounter;
        _cellCounter++;
    }
    else if(className == "Continent")
    {
        oss << _continentCounter;
        _continentCounter++;
    }
    else if(className == "LeftContinentExtremity"
            || className == "RightContinentExtremity")
    {
        oss << _continentExtremityCounter;
        _continentExtremityCounter++;
    }
    else if(className == "WarmingZone")
    {
        oss << _warmingZoneCounter;
        _warmingZoneCounter++;
    }
    else if(className == "Collision")
    {
        oss << _collisionCounter;
        _collisionCounter++;
    }
    else if(className == "ActiveMargin")
    {
        oss << _activeMarginCounter;
        _activeMarginCounter++;
    }
    else
    {
        oss << _instanceCounter;
    }
    
    id += oss.str();
    _instanceCounter++;
    return id;
}

/*
 =======================
 Structure
 =======================
 */

void
Earth::checkNeighbors()
{
    // is used at initialization of the run, when no plate velocity is known yet
    // and neighbor needs to be imposed
    if(_continents.size() == 0)
        return;
    
    for(unsigned int i=0; i<_continents.size(); i++)
    {
        // Left
        for(unsigned int j=0; j<_elements.size(); j++)
        {
            double distance =
            ::computeAbsDistance(_continents[i]->getLeftPosition(),
                                 _elements[j]->getPosition());
            
            if(distance < Earth::creationDistance)
            {
                // position
                _elements[j]->setPosition(_continents[i]->getLeftPosition());
                
                // neighbor
                _continents[i]->getLeftExtremity()->setNeighbor(_elements[j]);
                if(_elements[j]->isA("Subduction"))
                    ((Subduction*)_elements[j])->setRightNeighborExtremity(_continents[i]->getLeftExtremity());
                else if(_elements[j]->isA("Staple"))
                    ((Staple*)_elements[j])->setRightNeighborExtremity(_continents[i]->getLeftExtremity());
                
                // age
                _elements[j]->setRightAge(_continents[i]->getLeftExtremity()->getRightAge());
            }
        }
        
        // Right
        for(unsigned int j=0; j<_elements.size(); j++)
        {
            double distance = ::computeAbsDistance(_continents[i]->getRightPosition(),
                                                   _elements[j]->getPosition());
            if(distance < Earth::creationDistance)
            {
                // position
                _elements[j]->setPosition(_continents[i]->getRightPosition());
                
                // neighbor
                _continents[i]->getRightExtremity()->setNeighbor(_elements[j]);
                if(_elements[j]->isA("Subduction"))
                    ((Subduction*)_elements[j])->setLeftNeighborExtremity(_continents[i]->getRightExtremity());
                else if(_elements[j]->isA("Staple"))
                    ((Staple*)_elements[j])->setLeftNeighborExtremity(_continents[i]->getRightExtremity());
                
                // age
                _elements[j]->setLeftAge(_continents[i]->getRightExtremity()->getLeftAge());
            }
        }
    }
}

// --------------------------------------------------
void
Earth::fillStructure(bool putAges)
{ /* Creates plates and sections from elements. 
   putAges is optional (default = true) */
    
    sortElements();
	
    // --- Elements + continent extremities
    vector<GeoElement*> elements;
    for(unsigned int i=0; i<_elements.size(); i++)
        elements.push_back(_elements[i]);
    
    for(unsigned int i=0; i<_continents.size(); i++)
    { /* add the border of continents as new elements if
       there are not already subductions or staple */
        if(!_continents[i]->getRightExtremity()->getNeighbor())
            elements.push_back(_continents[i]->getRightExtremity());
        if(!_continents[i]->getLeftExtremity()->getNeighbor())
            elements.push_back(_continents[i]->getLeftExtremity());
        // Test : double checking continents
        if (!_continents[i]->getRightExtremity())
        {
            cerr << "ERROR: " << _continents[i]->getId() << " does not have a right extremity" << endl;
        }
        if (!_continents[i]->getLeftExtremity())
        {
            cerr << "ERROR: " << _continents[i]->getId() << " does not have a left extremity" << endl;
        }
    }
    
    sortElements(elements);
    
    // --- create all plate sections
    vector<PlateSection*> sections;
    
    vector<GeoElement*>::iterator it = elements.begin();
    vector<GeoElement*>::iterator itp = it + 1;
    while(true)
    {
        if(itp == elements.end())
        {
            sections.push_back(new PlateSection(this, (*it), elements[0]));
            break;
        }
        
        sections.push_back(new PlateSection(this, (*it), (*itp)));
        it ++;
        itp ++;
    }
    
    // --- create vector of only interfaces (ridge & subduction)
    vector<Interface*> interfaces;
    for(unsigned int i=0; i<_elements.size(); i++)
    {
        if(_elements[i]->isA("Interface"))
            interfaces.push_back((Interface*)_elements[i]);
    }
    
    // --- check if there is at least one subduction
    _noSubduction = true;
    it = elements.begin();
    while(true)
    {
        if((*it)->isA("Subduction"))
        {
            _noSubduction = false;
            break;
        }
        it ++;
        if(it == elements.end())
            break;
    }
    
    // --- vector of sections that will form a plate:
    vector<PlateSection*> localSections;
    
    // --- create plates
    if(_noSubduction)
    { /* create only one plate if there is no subduction
       ( = no mobile plate, stuck plate tectonics ) */
        vector<GeoElement*>::iterator itStart = elements.begin();
        if(!interfaces.empty())
        { /* starts the plate at the first found interface (ridge)
           if there is one (otherwise will start at elements.begin() */
            while(true)
            {
                if((*itStart) == interfaces[0]) //correct to write like that?
                    break;
                itStart ++;
            }
        }
        
        it = itStart;
        while(true)
        {
            localSections.push_back((*it)->getLeftSection());
            it++;
            if(it == elements.end())
                it = elements.begin();
            if(it == itStart)
            { // did the full loop:
                break;
            }
        }
        
        createPlate(localSections);
        localSections.clear();
        
        /* in that case: some ridges may not have yet _leftPlate
         and _rightPlate defined */
        for(unsigned int i=0; i<interfaces.size(); i++)
        {
            if(!interfaces[i]->getRightPlate())
                interfaces[i]->setRightPlate(_plates[0]);
            if(!interfaces[i]->getLeftPlate())
                interfaces[i]->setLeftPlate(_plates[0]);
        }
        
    } // _noSubduction
    else  // !_noSubduction
    { /* there is at least one subduction: at least 2 plates */
        vector<GeoElement*>::iterator itStart = elements.begin();
        for(unsigned int i=0; i<interfaces.size(); i++)
        { // put itStart at interfaces[i]
            while(true)
            {
                if((*itStart) == interfaces[i])
                    break;
                itStart++;
            }
            // push sections between itStart and next interface:
            it = itStart;
            while(true)
            {
                localSections.push_back((*it)->getLeftSection());
                it++;
                if(it == elements.end())
                    it = elements.begin();
                if((*it)->isA("Interface"))
                    break;
            }
            // create plate:
            createPlate(localSections);
            localSections.clear();
        }
    } // !_noSubduction
    sections.clear();
    
    if(_continents.size() > 0)
    {
        resetContinents(); // give a section to each continent
        /* give left/right section to extremity that have neighbor elements
         (allows shorter writing when needing length of right or left section*/
        for(unsigned int i=0; i<_continents.size(); i++)
        {
            if(_continents[i]->getRightExtremity()->getNeighbor())
            {
                PlateSection* rightSection =
                _continents[i]->getRightExtremity()->getNeighbor()->getRightSection();
                _continents[i]->getRightExtremity()->setRightSection(rightSection);
            }
            if(_continents[i]->getLeftExtremity()->getNeighbor())
            {
                PlateSection* leftSection =
                _continents[i]->getLeftExtremity()->getNeighbor()->getLeftSection();
                _continents[i]->getLeftExtremity()->setLeftSection(leftSection);
            }
        }
    }
    
    // - - - - - - - - - - - - - - - - - - - - -
    if(putAges)
    { /* each old element 'knows' the ages on its right and left
       -> use them to re-put ages for each section.
       */
        for(unsigned int i=0; i<_plates.size(); i++)
        {
            vector<PlateSection*>& sections = _plates[i]->accessSections();
            for(unsigned int j=0; j<sections.size(); j++)
            {
                GeoElement* rightElement = sections[j]->getRightElement();
                
                if(rightElement->sectionAgesIsEmpty(LEFT))
                { // only for newly created sections (e.g. around new ridges)
                    sections[j]->initAges();
                }
            } // sections
        } // _plates
    } //putAges
    
    
    // updateCells() is needed to avoid jumps in the graphic interface
    updateCells();
}
// --------------------------------------------------
bool
Earth::onlyOneSection()
{
    unsigned int nSections = 0;
    for(unsigned int i=0; i<_plates.size(); i++)
    {
        vector<PlateSection*> sections = _plates[i]->accessSections();
        nSections += sections.size();
    }
    return (nSections < 2);
}
// --------------------------------------------------
void
Earth::updateCells()
{ // CG
    if(_drawCells)
    {
        for(unsigned int i=0; i<_plates.size(); i++)
            _plates[i]->setCells();
    }
}
// --------------------------------------------------
vector<pair<double, double> >
Earth::setAgeSlices(double step, double aMax, double aMin)
{
    vector<pair<double, double> > slices;
    
    double agemax = aMax;
    while(true)
    {
        double agemin = agemax - step;
        if(agemin < aMin)
            break;
        
        pair<double, double> slice =
        make_pair(agemax, agemin);
        
        slices.push_back(slice);
        
        agemax -= step;
    }
    
    return slices;
}
// --------------------------------------------------
void
Earth::setContinentalAgeSlices()
{ /* creates the vector of pair<double,double>
   _continentalAgeSlices, which contains age slices
   that will be used for different colors of continents */
    double step = Earth::contAgeGroupStep;
    double agemax = Earth::contAgeGroupMax;
    double agemin = Earth::contAgeGroupMin;
    
    _continentalAgeSlices = setAgeSlices(step, agemax, agemin);
}
// --------------------------------------------------
void
Earth::setOceanAgeSlices()
{ /* creates the vector of pair<double,double>
   _oceanAgeSlices, which contains age slices
   that will be used for different colors of oceans */
    double step = Earth::oceanAgeGroupStep;
    double agemax = Earth::oceanAgeGroupMax;
    double agemin = Earth::oceanAgeGroupMin;
    
    _oceanAgeSlices = setAgeSlices(step, agemax, agemin);
}

/*
 =======================
 Computing
 =======================
 */
void
Earth::compute()
{
    if(_running)
    {
        lock();
        
        _timeLoop(_time, _dt);
        
        _time += _dt;
        _ageMa -= _dt / 1.0E6;
        
        _timestep ++;
        
        unlock();
    }
}
// --------------------------------------------------
/* actual computing steps are here */
void
Earth::_timeLoop(double time, double dt)
{
    if(_timestep != 0)
    {
        if(!_fixed_timestep)
            dt = _adaptDt();
        
        _noSubductionPrec = _noSubduction;
        
        /* move elements, continents, ages, and take
         care of collisions */
        updatePositions();
        
        // --- updating cells (CG March 2013)
        updateCells();
        
        _checkChangePlateTectonics(); // just for inline output
    }
    
    // writing log, ages etc. if time is reached
    _writeOutputs();
    
    _noSubductionPrec = _noSubduction;
    if(_inactiveRidges)
        _activateAllRidges();
    // =================================
    // --- Old oceanic seafloor dives
    _checkDiving();
    _checkAlwaysSubduct();
    
    _checkChangePlateTectonics(); // just for inline output
	
    // ==============================
    // --- Continents break up
    _checkOceanOpening();
    
    // ==============================
    // --- Extension in upper plate above subduction: new ridge
    _checkRidgeCreation();
    
	// ==============================
    // compute plate velocities
    _computeVelocities();
    // ==============================
    // check that all subductions are really diving    	
	_verifySubductions();
	
    // check that all the ridges are really extensive
    _verifyRidges();
    
    // ==============================
    // CG - Updating temperature, viscosities and critical ages
    if(_continental_growth)
    {
        updateSurfaces(); //oc. and cont. surfaces needed for heat flux
        updateMasses();
    }
    // ==============================
    // -- Update heat flux, using age distribution
    updateQtot();
    
    // -- radiogenic heating
    updateConcentrations();
    updateRadioHeatTW();
    
    // ==============================
    // -- Update temperature and its dependencies
    updateT();
    updateTData(); // data that depend on global T (eta, tau_ssc, rho etc.)

}
// --------------------------------------------------
void
Earth::updateQtot()
{ 
    _QtotTW = 0.0;
    
    for(unsigned int i=0; i<_plates.size(); i++)
    {
        vector<PlateSection*>& sections = _plates[i]->accessSections();
        for(unsigned int j=0; j<sections.size(); j++)
        {
            sections[j]->updateQtot();
            _QtotTW += sections[j]->getQtot() / 1.0E12;            
        }
    }
}
// --------------------------------------------------
void
Earth::updateT()
{ // dt is in years
    _T += computedTdt() * ::yr_to_sec(_dt);
}
// --------------------------------------------------
void
Earth::initTData()
{ 
    updateViscosities();
    updateDensities();
    
    if (_SSCMode == SSCCONSTANT)
    {
        _tau_ssc = _tau_ssc_p;
    }
    else
    {
        updateTauC();
    }
    
    updateTauSub();
}
// --------------------------------------------------
void
Earth::updateTData()
{ // CG
    updateViscosities();
    updateDensities();
    updateTauC();
    updateTauSub();
}
// --------------------------------------------------
void 
Earth::updateViscosities()
{ // CG
    _eta_m = computeEta(_eta_m_p, _E_m);
    _eta_um = computeEta(_eta_um_p, _E_um);
    _eta_pl = computeEta(_eta_pl_p, _E_pl);
    
    _eta_ast = _eta_um;
    if(_thick_ast > 0.0)
        _eta_ast = computeEta(_eta_ast_p, _E_um);
    
    _eta_subcont = _eta_um;
    if(_thick_subcont > 0.0)
        _eta_subcont = computeEta(_eta_subcont_p, _E_um);
}
// --------------------------------------------------
void 
Earth::updateDensities()
{
    _rho_um = computeRho(_rho_um_p, _alpha_um, _T, _T_p);
    
    updateTPlate();
    _rho_pl = computeRho(_rho_pl_p, _alpha_pl, _T_plate, _T_plate_p);
}
// --------------------------------------------------
void
Earth::updateTPlate()
{
    _T_plate = compute_T_plate(_T);
}
// --------------------------------------------------
double
Earth::computeQmeanOcean()
{ // in W/m^2
    double QmeanOcean = 0.0;
    double QtotOcean = 0.0;
    for(unsigned int i=0; i<_plates.size(); i++)
    {
        vector<PlateSection*> sections = _plates[i]->accessSections();
        for(unsigned int j=0; j<sections.size(); j++)
        {
            if(!sections[j]->isContinental())
            {
                QtotOcean += sections[j]->getQtot(); // W
            }
        }
    }
    QmeanOcean = QtotOcean / _surfTotOcean;
    return QmeanOcean;
}
// --------------------------------------------------
double 
Earth::computeQmeanCont()
{ // in W/m^2
    double QmeanCont = 0.0;
    double QtotCont = 0.0;
    if(_continents.size() > 0)
    {
        for(unsigned int i=0; i<_continents.size(); i++)
        {
            QtotCont += _continents[i]->getPlateSection()->getQtot();
        }
        QmeanCont = QtotCont / _surfTotContinent;
    }
    
    return QmeanCont;
}
// --------------------------------------------------
double
Earth::computeMeanPlateThickness()
{ // results in m
    double meanH = 0.0;
    double totalLength = 0.0;
    
    for(unsigned int i=0; i<_plates.size(); i++)
    {
        vector<PlateSection*> sections = _plates[i]->accessSections();
        for(unsigned int j=0; j<sections.size(); j++)
        {
            if(!sections[j]->isContinental())
            {
                double L = sections[j]->computeLength(); // length in deg.
                totalLength += L;
                meanH += L * sections[j]->computeMeanThickness();
            }
        }
    }
    
    meanH /= totalLength;
    return meanH;
}
// --------------------------------------------------
void
Earth::computeSeaLevel()
{
    _VbelowRidge_HSCM = 0.0;
    _VbelowRidge_PM95 = 0.0;
    _VbelowRidge_PM125 = 0.0;
    
    for(unsigned int i=0; i<_plates.size(); i++)
    {
        vector<PlateSection*>& sections = _plates[i]->accessSections();
        for(unsigned int j=0; j<sections.size(); j++)
        {
            sections[j]->updateBathymetry();
            _VbelowRidge_HSCM += sections[j]->getVbelowRidge_HSCM();
            _VbelowRidge_PM95 += sections[j]->getVbelowRidge_PM95();
            _VbelowRidge_PM125 += sections[j]->getVbelowRidge_PM125();
        }
    }
    
    _VaboveRidge_HSCM = _present_oceanVolume - _VbelowRidge_HSCM;
    _VaboveRidge_PM95 = _present_oceanVolume - _VbelowRidge_PM95;
    _VaboveRidge_PM125 = _present_oceanVolume - _VbelowRidge_PM125;
    
    //////////////////////////////////
    // warning: computation not valid if _VaboveRidge < 0.0
    //////////////////////////////////
    //double ratio = (1.0 - _cont_surfRatio) / (1.0 - _present_cont_surfRatio);
    double ratio = (1.0 - _present_cont_surfRatio) / (1.0 - _cont_surfRatio); // correction on April 17 2018
    _ridgeDepth_HSCM = _present_ridgeDepth *
    _VaboveRidge_HSCM / _present_VaboveRidge * ratio;
    
    _ridgeDepth_PM95 = _present_ridgeDepth *
    _VaboveRidge_PM95 / _present_VaboveRidge * ratio;
    
    _ridgeDepth_PM125 = _present_ridgeDepth *
    _VaboveRidge_PM125 / _present_VaboveRidge * ratio;
    
    _sealevel_HSCM = (1.0 - _rho_seawater / _rho_um) *
    (_ridgeDepth_HSCM - _present_ridgeDepth);
    
    _sealevel_PM95 = (1.0 - _rho_seawater / _rho_um) *
    (_ridgeDepth_PM95 - _present_ridgeDepth);
    
    _sealevel_PM125 = (1.0 - _rho_seawater / _rho_um) *
    (_ridgeDepth_PM125 - _present_ridgeDepth);
    
    
}
// --------------------------------------------------
double
Earth::computeSeafloorProduction()
{ // in km^2/yr
    double prod = 0.0;
    for (unsigned int i=0; i<_elements.size(); i++)
    {
        if (_elements[i]->isA("Ridge"))
        {
            Ridge* ridge = (Ridge*)_elements[i];
            double URight = 0.0;
            double ULeft = 0.0;
            if(ridge->getRightPlate())
                URight = cm_to_deg(ridge->getRightPlate()->getU());
            if(ridge->getLeftPlate())
                ULeft = cm_to_deg(ridge->getLeftPlate()->getU());
            
            if (::isStrictlyLess(ULeft, URight))
            {
                /* Do nothing. send error message only if significant
                 velocities */
                if(URight - ULeft > 1.0E-3)
                {
                    cerr << "Error in computeSeafloorProduction: ULeft < URight "
                    << "  (" << ULeft << ", " << URight << ")"
                    << " for " << ridge->getShortId()
                    << " at pos = " <<  ridge->getPosition()
                    << endl;
                }
            }
            else
                prod += deg_to_squareM((ULeft - URight)) * 1.0E-6; // in km2/yr
        }
    }
    
    return prod;
}
// --------------------------------------------------
double
Earth::computeSlabFlux(bool weight)
{ // in km^2/yr
    double flux = 0.0;
	double sum_w = 0.0;
	unsigned int nsubd = 0;
    for(unsigned int i=0; i<_elements.size(); i++)
    {
        if(_elements[i]->isA("Subduction"))
        {
            Subduction* subduction = (Subduction*)_elements[i];
            double URight = 0.0;
            double ULeft = 0.0;
            if(subduction->getRightPlate())
                URight = cm_to_deg(subduction->getRightPlate()->getU());
            if(subduction->getLeftPlate())
                ULeft = cm_to_deg(subduction->getLeftPlate()->getU());
            
            if(::isStrictlyLess(URight, ULeft)) // going in the wrong direction
            {
                // do nothing. Send message only if significant neg.
                if(URight - ULeft < -1.0E-3)
                {
                    cerr << "Error in computeSlabFlux():"
                    << " subduction not converging."
                    << endl;
                }
            }
            else 
			{		
				double w = 1.0;
				if (weight)
				{
					if(subduction->isA("RightSubduction"))
						w = subduction->getLeftPlate()->computeLength();
					else
						w = subduction->getRightPlate()->computeLength();
					sum_w += w;
					nsubd ++;
				}						
				flux += w * deg_to_squareM((URight - ULeft)) * 1.0E-6; // in km2/yr
			}
        } // isA("Subduction")
    } // elements[i]

	if (weight && sum_w > 0.0)
		flux *= double(nsubd) / sum_w;
    
    return flux;
}
// --------------------------------------------------
double
Earth::computedTdt()
{
    double dT_dt = 0.0;
    
    double radio = _radioHeatTW * 1.0E12;
    double Q     = _QtotTW * 1.0E12;
    
    dT_dt = ( -Q + radio) / ( _earthMass * _C_p );
    
    return dT_dt;
}
// ---------------------------------------------------
void
Earth::initConcentrations()
{
    /* order for _concentrations is
     i=0-> U235, 1->U238, 2->Th232, 3->K40 */
    
    // U235
    Concentration U235;
    U235.id = "235U";
    U235.heat = getFixedParameter("heat235U"); // W/kg
    U235.lambda = log(2) / getFixedParameter("halfLife235U"); // 1/Myr
    U235.presentPrimitive = getModelParameter("U_BSE") * getFixedParameter("235U_over_U");
    U235.continent = getModelParameter("U_cont") * getFixedParameter("235U_over_U");
    _concentrations.push_back(U235);
    
    //U238
    Concentration U238;
    U238.id = "238U";
    U238.heat = getFixedParameter("heat238U"); // W/kg
    U238.lambda = log(2) / getFixedParameter("halfLife238U"); // 1/Myr
    U238.presentPrimitive = getModelParameter("U_BSE") * getFixedParameter("238U_over_U");
    U238.continent = getModelParameter("U_cont") * getFixedParameter("238U_over_U");
    _concentrations.push_back(U238);
    
    //Th232
    Concentration Th232;
    Th232.id = "232Th";
    Th232.heat = getFixedParameter("heat232Th"); // W/kg
    Th232.lambda = log(2) / getFixedParameter("halfLife232Th"); // 1/Myr
    Th232.presentPrimitive = getModelParameter("Th_BSE");
    Th232.continent = getModelParameter("Th_cont");
    _concentrations.push_back(Th232);
    
    //K40
    Concentration K40;
    K40.id = "40K";
    K40.heat = getFixedParameter("heat40K"); // W/kg
    K40.lambda = log(2) / getFixedParameter("halfLife40K"); // 1/Myr
    K40.presentPrimitive = getModelParameter("K_BSE") * getFixedParameter("40K_over_K");
    K40.continent = getModelParameter("K_cont") * getFixedParameter("40K_over_K");
    _concentrations.push_back(K40);
    
    
    cout << "present-day concentrations(ppm)" << endl;
    cout << "U235 \t U238 \t Th232 \t K40" << endl;
    
    for(unsigned int i=0; i<_concentrations.size(); i++)
    {
        _concentrations[i].presentDepleted =
        ( _concentrations[i].presentPrimitive * _primitive_mantleMass  -
         _concentrations[i].continent * _present_continentMass ) /
        _present_mantleMass;
        
        cout << "\t" <<  _concentrations[i].presentDepleted;
    }
    cout << endl;
    
    /* depletion ratio fixed or not:
     fixed if depletionMode is DEPLETED or PRIMITIVE.
     Also fixed if continents are not growing. */
    switch(_depletionMode)
    {
        case ADAPTATIVE :
        {
            _fixed_depletion = false;
            if( !_continental_growth )
            {
                _fixed_depletion = true;
                _depletionRatio =
                _cont_surfRatio / _present_cont_surfRatio;
            }
            break;
        }
        case DEPLETED :
        {
            _fixed_depletion = true;
            _depletionRatio = 1.0;
            break;
        }
        case PRIMITIVE :
        {
            _fixed_depletion = true;
            _depletionRatio = 0.0;
            break;
        }
        default :
        {
            _fixed_depletion = false;
            break;
        }
    }
    
    
    // for comparison: prepare always depleted and always primitive cases
    for(unsigned int i=0; i<_concentrations.size(); i++)
    {
        _depleted.push_back(_concentrations[i]);
        _primitive.push_back(_concentrations[i]);
    }
    
}
// --------------------------------------------------
void
Earth::updateConcentrations()
{
    if(_constant_heating)
        return;
    
    if(!_fixed_depletion)
        _depletionRatio = _cont_surfRatio / _present_cont_surfRatio;
    
    for(unsigned int i=0; i<_concentrations.size(); i++)
    {
        _concentrations[i].value =
        ( (1.0 - _depletionRatio) * _present_continentMass * _concentrations[i].continent  +
         _present_mantleMass * _concentrations[i].presentDepleted
         ) * exp(_ageMa * _concentrations[i].lambda) /
        ( (1.0 - _depletionRatio) * _present_continentMass + _present_mantleMass );
    }
}
// --------------------------------------------------
void
Earth::updateReferenceConcentrations()
{
    for(unsigned int i=0; i<_primitive.size(); i++)
    {
        _primitive[i].value =
        _primitive[i].presentPrimitive * exp(_ageMa * _primitive[i].lambda);
    }
    
    for(unsigned int i=0; i<_depleted.size(); i++)
    {
        _depleted[i].value =
        _depleted[i].presentDepleted * exp(_ageMa * _depleted[i].lambda);
    }
    
}
// --------------------------------------------------
void
Earth::updateRadioHeatTW()
{
    _radioHeatTW = computeRadioHeatTW(_concentrations);
}
// --------------------------------------------------
double
Earth::computeRadioHeatTW(vector<Concentration> conc)
{
    if(_constant_heating)
        return _constant_heat_value;
    
    
    double radio = 0.0;
    
    double coeff = _mantleMass * 1.0E-18;// 1E-18 = 1E-6 (ppm) * 1E-12 (TW)
    for(unsigned int i=0; i<conc.size(); i++)
    {// .heat: in W/kg; .value: in ppm; coeff: in kg*1e-18
        radio += coeff *
        conc[i].heat * conc[i].value;
    }
    
    return radio;
}
// --------------------------------------------------
double
Earth::radioHeatWperKg()
{ // output is in W/kg
    return _radioHeatTW * 1.0E12 / _mantleMass;
}
// --------------------------------------------------
double
Earth::computeEta(double eta_p, double E)
{
    double eta = eta_p;
    if(_computeEta)
        eta = eta_p * exp( E * ( (1.0/_T) - (1.0/_T_p)) / Earth::gasConstant);
    
    return eta;
}
// --------------------------------------------------
double
Earth::computeRho(double rho_p, double alpha, double T, double Tp)
{
    return rho_p * (1.0 - alpha * (T - Tp));
}
// --------------------------------------------------
double
Earth::compute_T_plate(double T)
{
    return _T_surf + (T - _T_surf) * Earth::coeff_T_plate;
}
// --------------------------------------------------
void
Earth::updateTauC()
{ //CG
    if (_SSCMode == SSCCONVECTIVE)
    {
        double twothird = 2.0 / 3.0;
        _tau_ssc = _tau_ssc_p * pow( (_T_p - _T_surf) / (_T - _T_surf), twothird ) *
        exp( twothird * _E_um * ((1.0/_T) - (1.0/_T_p)) / Earth::gasConstant);
    }
}
// --------------------------------------------------
void
Earth::updateTauSub()
{ //CG
    // Brittle criterion
    if(_initSubduction.mode == BRITTLE)
    {
        _tau_sub = _tau_sub_p * pow( (_T_p - _T_surf) / (_T - _T_surf), 2.0 );
    }
    // Convective criterion
    else if(_initSubduction.mode == CONVECTIVE)
    {
        double twothird = 2.0 / 3.0;
        _tau_sub = _tau_sub_p * pow( (_T_p - _T_surf) / (_T - _T_surf), twothird )
        * exp( twothird * _E_um * ((1.0/_T) - (1.0/ _T_p)) / Earth::gasConstant);
    }
    // Constant critical age
    else if(_initSubduction.mode == CONSTANT)
        _tau_sub = _tau_sub_p;
    
    
    if(false)
    {
        // TEST //////////////////////////////////////////
        if(_noSubduction)
            _timespan_noSubduction += _dt * 1.0E-6; // in Myr
        else if(_timespan_noSubduction > 0.0)
            _timespan_noSubduction = 0.0; // back to zero when restarts
        // END TEST ////////////////////////////////////////////////
        
        if(_noSubduction)
        {
            cout << "no subd: "
            << _timespan_noSubduction
            << "\t tausub = " << _tau_sub
            << endl;;
            
            if(_timespan_noSubduction > _tau_sub*0.333)
            {
                _tau_sub *= 0.333;
                cout << "\t -> taueff = " << _tau_sub << endl;
            }
        }
    }
    
}
// --------------------------------------------------
void
Earth::updateAgesDistribution()
{
    vector<double> ages = getOceanAges();
    sort(ages.begin(), ages.end());
    
    _agesDistribution = ::getDistribution(ages, _ageSlice);
}
// --------------------------------------------------
vector<pair<double, double> > 
Earth::getAgesDistribution()
{
    updateAgesDistribution();
    return _agesDistribution;
}
// --------------------------------------------------
void
Earth::updateSurfaces()
{ /* _surfTot must be computed before.
   This also updates masses of continents and
   mantle (if true). */
    
    _surfTotOcean = 0.0;
    _surfTotContinent = 0.0;
    
    for(unsigned int i=0; i<_plates.size(); i++)
    {
        vector<PlateSection*> sections = _plates[i]->accessSections();
        for(unsigned int j=0; j<sections.size(); j++)
        {
            if(!sections[j]->isContinental())
            {
                _surfTotOcean +=
                ::computeRightLeftDistance( sections[j]->getRightElement()->getPosition(),
                                           sections[j]->getLeftElement()->getPosition() );
            }
            else
                _surfTotContinent +=
                ::computeRightLeftDistance( sections[j]->getRightElement()->getPosition(),
                                           sections[j]->getLeftElement()->getPosition() );
        }
    }
    /* _surfTotOcean and ...Continent are in degrees so far */
    _surfTotOcean = deg_to_squareM(_surfTotOcean);
    _surfTotContinent = deg_to_squareM(_surfTotContinent);
    
    /* ratio between cont and total surface: used in radioactive heating */
    _cont_surfRatio = _surfTotContinent / _surfTot;
}
// --------------------------------------------------
void
Earth::updateMasses()
{
    _continentMass = _surfTotContinent * _rho_cont * _cont_thickness;
    _mantleMass = _primitive_mantleMass - _continentMass;
}
// --------------------------------------------------
double
Earth::thickness_to_age(double h)
{ /* returns age (in Myr) corresponding to a given 
   thickness h (in m) of the lithosphere */
    double age = - numeric_limits<double>::infinity();
    if(h > 0.0)
        age = ::sec_to_Myr( pow( (h / coeff_plate_thickness), 2.0 )
                           / _kappa );
    
    return age;
}
// --------------------------------------------------
void
Earth::updateElementsVelocity()
{
    if(_fixed_configuration)
        return;
    
    for(unsigned int i=0; i<_elements.size(); i++)
    {
        _elements[i]->updateVelocity();
    }
    
    if(_continents.size() > 0)
    {
        for(unsigned int i=0; i<_continents.size(); i++)
        { /* includes continental growth */
            _continents[i]->getRightExtremity()->updateVelocity();
            _continents[i]->getLeftExtremity()->updateVelocity();
        } // _continents[i]
        
    } // _continents.size > 0
}
// --------------------------------------------------
void
Earth::updatePositions()
{ /* --------------------------------------------------
   update positions of elements (RCE and LCE included),
   with taking possible collisions into account, and
   adapting _dt to not have sections overlapping
   -------------------------------------------------- */
    
    updateElementsVelocity(); // give the correct velocity to each element to be moved
    
    double dtMin = _dt;
    
    /* Collision = struct{double dt, GeoElement* right, GeoElement* left}*/
    vector<Contact> collisions = findCollisions(dtMin); /* dtMin is an output:
                                                         dt for the 'closest' collision */
    _dt = dtMin;
    
    /* move elements */
    for(unsigned int i=0; i<_elements.size(); i++)
    {
        _elements[i]->updatePosition(_dt);
    }
    
    for(unsigned int i=0; i<_continents.size(); i++)
    {
        _continents[i]->updatePosition(_dt);
    }
    
    for(unsigned int i=0; i<_markers.size(); i++)
    { /* Markers are moving like their continents */
        _markers[i]->updatePosition(_dt);
    }
    
    /* move ages of each sections + creating new age points at ridges
     + removing subducted points. */
    moveAges();
    
    /* deal with collisions */
    if(!collisions.empty())
    {
        treatCollisions(collisions);
        correctAges();
    }
    
    /* make sure continents' extremities have the same
     ages as neighbors */
    for(unsigned int i=0; i<_continents.size(); i++)
        _continents[i]->shareNeighborAges();
    
}
// =======================================================
double
Earth::_adaptDt()
{ //CG: adapt timestep to resolution and max velocity of plates
    double dtMax = _timeStep; // in years
    double dtMin = 1.0;
    
    double vmax = 0.0;
    for(unsigned int i=0; i<_plates.size(); i++)
        vmax = max(vmax, fabs(_plates[i]->getU())); // cm/yr
    
    vmax = cm_to_deg(vmax); // degrees/yr
    
    if(vmax != 0.0)
        _dt = _courant * _resolution / vmax; // in yr
    else
        _dt = dtMax;
    
    _dt = max(_dt, dtMin);
    _dt = min(_dt, dtMax);
    return _dt;
}
// --------------------------------------------------
void
Earth::moveAges()
{
    for(unsigned int i=0; i<_plates.size(); i++)
    {
        vector<PlateSection*>& sections = _plates[i]->accessSections();
        for(unsigned int j=0; j<sections.size(); j++)
        { // continental sections are done here too
            sections[j]->moveAges();
        }
    }
}
// =======================================================
string
Earth::getTimeString()
{
    ostringstream oss;
    oss << setw(7)
    << setprecision(1)
    << fixed << setfill('0') << getTimeMyr()
    << "Myr";
    
    return oss.str();
}
// --------------------------------------------------
void
Earth::logToFile(bool all)
{ /*  write average data (temperature, heat flow, velocity... ) 
   to earth.log only if !all,
   and to earth.log, ages.log and forces.log if all.
   */
    
    if(_logFile.is_open())
    {
        double QmeanOcean = computeQmeanOcean() * 1.0E3; //mW.m^-2
        double QmeanCont = computeQmeanCont() * 1.0E3; //mW.m^-2
        double Urey = _radioHeatTW / _QtotTW;
        
        // for checking total plate length
        double Lsum = 0.0;
        
        // Plates velocity and size + continents - - - - - - - - - -
        //    pair: first is unweighted, second is weighted by plates' length
        pair<double, double> Uabs = make_pair(0.0, 0.0);
        pair<double, double> Umean = make_pair(0.0, 0.0);
        pair<double, double> Urms = make_pair(0.0, 0.0);
        double UabsMax = - numeric_limits<double>::infinity();
        double UabsMin = + numeric_limits<double>::infinity();
        for(unsigned int i=0; i<_plates.size(); i++)
        {
            double L = _plates[i]->computeLm(); // in m
            Lsum += L;
            
            double U = _plates[i]->getU(); // in cm/yr
            Urms.first += pow(U, 2.0);
            Urms.second += pow(U, 2.0) * L;
            if(all)
            {
                Uabs.first += fabs(U);
                Uabs.second += fabs(U) * L;
                Umean.first += U;
                Umean.second += U * L;
                UabsMax = max(UabsMax, fabs(U));
                UabsMin = min(UabsMin, fabs(U));
            }
        } // _plates.size()
        
        Urms.first = sqrt(Urms.first / (double)_plates.size());
        Urms.second = sqrt(Urms.second / _earthPerimeter);
        if(all)
        {
            Uabs.first  /= ((double)_plates.size());
            Uabs.second  /= _earthPerimeter;
            
            Umean.first /= ((double)_plates.size());
            Umean.second /= _earthPerimeter;
        }
        
        
        if(!areEqual(Lsum, _earthPerimeter, 1.0))
        {
            cerr << "WARNING: total length of plates is not equal to Earth perimeter\t"
            << "t=" << _ageMa
            << "\t L, perim = " << Lsum << "\t" << _earthPerimeter << endl;
        } // check Lsum and earth total perimeter
        
        
        // Ages - - - - - - - - - - -
        vector<double> ages = getOceanAges(); /* all ages of all plates
                                               (except continental sections)*/
        sort(ages.begin(), ages.end()); // put ages in ascending order.
        pair<double,double> meanAge = ::getMeanSD(ages); /* first is mean,
                                                          second is SD */
        
        /* - - - - - - - - - - -
         writing in _logFile
         - - - - - - - - - - - */
        _logFile << setprecision(3) << fixed << _ageMa << "\t"
        << _T << "\t"
        << _QtotTW << "\t"
        << QmeanOcean << "\t"
        << QmeanCont << "\t"
        << _radioHeatTW << "\t"
        << Urey << "\t"
        << Urms.first << "\t"
        << Urms.second << "\t"
        << _plates.size() << "\t"
        << _continents.size() << "\t"
        << ages[ages.size()-1] << "\t"
        << meanAge.first << "\t"
        << meanAge.second << "\n";
        
        if(all)
        {
            unsigned int Nb_ridges = 0.0;
            unsigned int Nb_subductions = 0.0;
            
            for(unsigned int i=0; i<_plates.size(); i++)
                _plates[i]->updateType();
            
            /* - - - - -
             VELOCITIES
             - - - - - */
            if(_velocityFile.is_open())
            {
                /* ELEMENTS */
                double UrmsElements = 0.0;
                double UrmsRidges = 0.0;
                double UrmsSubductions = 0.0;
                
                for(unsigned int i=0; i<_elements.size(); i++)
                {
                    double U2 = pow(_elements[i]->getU(), 2.0);
                    UrmsElements += U2;
                    if(_elements[i]->isA("Ridge"))
                    {
                        UrmsRidges += U2;
                        Nb_ridges ++;
                    }
                    if(_elements[i]->isA("Subduction"))
                    {
                        UrmsSubductions += U2;
                        Nb_subductions ++;
                    }
                }
                UrmsElements = sqrt(UrmsElements / (double)_elements.size());
                if(Nb_ridges != 0)
                    UrmsRidges = sqrt(UrmsRidges / (double)Nb_ridges);
                if(Nb_subductions != 0)
                    UrmsSubductions = sqrt(UrmsSubductions / (double)Nb_subductions);
                
                /* CONTINENTS */
                pair<double, pair<double, double> > UrmsContinents =
                make_pair(0.0, make_pair(0.0, 0.0));
                pair<double, pair<double, double> > UabsContinents =
                make_pair(0.0, make_pair(0.0, 0.0));
                double lengthContinents = 0.0;
                double lengthPlates = 0.0;
                if(_continents.size() > 0)
                {
                    for(unsigned int i=0; i<_continents.size(); i++)
                    {
                        double U2 = pow(_continents[i]->getU(), 2.0);
                        UrmsContinents.first += U2;
                        double Uabs = fabs(_continents[i]->getU());
                        UabsContinents.first += Uabs;
                        double L = _continents[i]->getLength();
                        double Lplate = _continents[i]->getPlate()->computeLength();
                        lengthContinents += L;
                        lengthPlates += Lplate;
                        UrmsContinents.second.first += U2 * L;
                        UrmsContinents.second.second += U2 * Lplate;
                        UabsContinents.second.first += Uabs * L;
                        UabsContinents.second.second += Uabs * Lplate;
                    }
                    // UrmsContinents:
                    UrmsContinents.first =
                    sqrt(UrmsContinents.first / (double)_continents.size());
                    UrmsContinents.second.first =
                    sqrt(UrmsContinents.second.first / lengthContinents);
                    UrmsContinents.second.second =
                    sqrt(UrmsContinents.second.second / lengthPlates);
                    // UabsContinents:
                    UabsContinents.first /= (double)_continents.size();
                    UabsContinents.second.first /= lengthContinents;
                    UabsContinents.second.second /= lengthPlates;
                }
                
                /* PURELY OCEANIC PLATES */
                pair<double, double> UrmsOceans =
                make_pair(0.0, 0.0);
                double lengthOceans = 0.0;
                unsigned int nOceans = 0;
                for(unsigned int i=0; i<_plates.size(); i++)
                {
                    if(_plates[i]->getContinents().size() == 0)
                    {
                        nOceans += 1;
                        double U2 = pow(_plates[i]->getU(), 2.0);
                        double L = _plates[i]->computeLength();
                        lengthOceans += L;
                        UrmsOceans.first += U2;
                        UrmsOceans.second += U2 * L;
                    }
                }
                if(nOceans > 0)
                {
                    UrmsOceans.first = sqrt(UrmsOceans.first / (double)nOceans);
                    UrmsOceans.second = sqrt(UrmsOceans.second / lengthOceans);
                }

                /* SUBDUCTING AND UPPER PLATES */
                // first: simple rms, second: weight of plate length
                pair<double,double> UrmsUpperPlates =
                make_pair(0.0, 0.0);
                pair<double,double> UrmsSubdPlates =
                make_pair(0.0, 0.0);
                for(unsigned int i=0; i<_plates.size(); i++)
                {
                    if(_plates[i]->isSubducting())
                    {
                        double U2 = pow(_plates[i]->getU(), 2.0);
                        double L = _plates[i]->computeLm();
                        UrmsSubdPlates.first += U2;
                        UrmsSubdPlates.second += U2 * L;
                    }
                    else
                    {
                        double U2 = pow(_plates[i]->getU(), 2.0);
                        double L = _plates[i]->computeLm();
                        UrmsUpperPlates.first += U2;
                        UrmsUpperPlates.second += U2 * L;
                    }
                }
                UrmsUpperPlates.first =
                sqrt(UrmsUpperPlates.first / (double)_plates.size());
                UrmsUpperPlates.second =
                sqrt(UrmsUpperPlates.second / _earthPerimeter);
                UrmsSubdPlates.first =
                sqrt(UrmsSubdPlates.first / (double)_plates.size());
                UrmsSubdPlates.second =
                sqrt(UrmsSubdPlates.second / _earthPerimeter);
                
                _velocityFile << setprecision(3) << fixed << _ageMa << "\t"
                << _T << "\t"
                << Uabs.first << "\t"
                << Uabs.second << "\t"
                << Umean.first << "\t"
                << Umean.second << "\t"
                << UabsMin << "\t"
                << UabsMax << "\t"
                << UrmsElements << "\t"
                << UrmsRidges << "\t"
                << UrmsSubductions << "\t"
                << UrmsContinents.first << "\t"
                << UrmsContinents.second.first << "\t"
                << UrmsContinents.second.second << "\t"
                << UabsContinents.first << "\t"
                << UabsContinents.second.first << "\t"
                << UabsContinents.second.second << "\t"
                << UrmsOceans.first << "\t"
                << UrmsOceans.second << "\t"
                << UrmsSubdPlates.first << "\t"
                << UrmsUpperPlates.first << "\t"
                << UrmsSubdPlates.second << "\t"
                << UrmsUpperPlates.second << "\n";
                
            } // _velocityFile.is_open()
            
            
            /*  - - - - -
             AGES
             - - - - - */
            if(_ageFile.is_open())
            {
                vector<double> whiskers = ::getWhiskers(ages);
                
                // max age for subducting plate and upper plate
                double maxAgeSubduction = 0.0;
                double maxAgeUpperPlate = 0.0;
                for(unsigned int i=0; i<_elements.size(); i++)
                {
                    if(_elements[i]->isA("LeftSubduction"))
                    {
                        Subduction* subduction = (Subduction*)_elements[i];
                        maxAgeSubduction = max(maxAgeSubduction,
                                               subduction->getRightAge().getValue());
                        if(!subduction->getLeftNeighborExtremity())
                            maxAgeUpperPlate = max(maxAgeUpperPlate,
                                                   subduction->getLeftAge().getValue());
                    }
                    if(_elements[i]->isA("RightSubduction"))
                    {
                        Subduction* subduction = (Subduction*)_elements[i];
                        maxAgeSubduction = max(maxAgeSubduction,
                                               subduction->getLeftAge().getValue());
                        if(!subduction->getRightNeighborExtremity())
                            maxAgeUpperPlate = max(maxAgeUpperPlate,
                                                   subduction->getRightAge().getValue());
                    }
                } // _elements[i]
                
                // Writing in _ageFile
                _ageFile << setprecision(3) << fixed << _ageMa << "\t"
                << _T << "\t"
                << ages[0] << "\t";
                for(unsigned int i=0; i<whiskers.size(); i++)
                    _ageFile << whiskers[i] << "\t";
                _ageFile << ages[ages.size()-1] << "\t"
                << maxAgeSubduction << "\t"
                << maxAgeUpperPlate
                << endl;
            } // _ageFile.is_open()
            
            
            
            /* - - - - -
             FORCES
             - - - - -*/
            if(_forcesFile.is_open())
            {
                // get forces on each plate
                completeForces(); /* visco. are computed at each time step,
                                   but not total forces MD, VS and B.
                                   */
                
                // pair: first is normal sum, second is abs.
                pair<double, double> SP = make_pair(0.0, 0.0);
                pair<double, double> RP = make_pair(0.0, 0.0);
                pair<double, double> SS = make_pair(0.0, 0.0);
                pair<double, double> MD = make_pair(0.0, 0.0);
                pair<double, double> VS = make_pair(0.0, 0.0);
                pair<double, double> BD = make_pair(0.0, 0.0);
                
                for(unsigned int i=0; i<_plates.size(); i++)
                {
                    Force* forces = _plates[i]->accessForces();
                    // Driving forces
                    SP.first += forces->getSlabPull();
                    SP.second += fabs(forces->getSlabPull());
                    RP.first += forces->getRidgePush();
                    RP.second += fabs(forces->getRidgePush());
                    SS.first += forces->getSlabSuction();
                    SS.second += fabs(forces->getSlabSuction());
                    
                    // Resistive forces
                    MD.first += forces->getMantleDrag();
                    MD.second += fabs(forces->getMantleDrag());
                    VS.first += forces->getViscousShear();
                    VS.second += fabs(forces->getViscousShear());
                    BD.first += forces->getBending();
                    BD.second += fabs(forces->getBending());
                } // _plates[i]
                
                // RP, SP, and SS are with dimensions.
                SP.first /= 1.0E12;
                SP.second /= 1.0E12;
                RP.first /= 1.0E12;
                RP.second /= 1.0E12;
                SS.first /= 1.0E12;
                SS.second /= 1.0E12;
                
                // MD, VS, and BD are dimensionless before this.
                MD.first *= Earth::force0 / 1.E12;
                MD.second *= Earth::force0 / 1.E12;
                VS.first *= Earth::force0 / 1.E12;
                VS.second *= Earth::force0 / 1.E12;
                BD.first *= Earth::force0 / 1.E12;
                BD.second *= Earth::force0 / 1.E12;
                
                // mean thickness of plates
                double meanPlateThickness =
                computeMeanPlateThickness() * 1.0E-3; // from m to km
                
                // writing in _forcesFile
                _forcesFile << setprecision(3) << fixed << _ageMa << "\t"
                << _T << "\t"
                << Nb_ridges << "\t"
                << Nb_subductions << "\t"
                << scientific
                << SP.first << "\t"
                << SP.second << "\t"
                << RP.first << "\t"
                << RP.second << "\t"
                << SS.first << "\t"
                << SS.second << "\t"
                << MD.first << "\t"
                << MD.second << "\t"
                << VS.first << "\t"
                << VS.second << "\t"
                << BD.first << "\t"
                << BD.second << "\t"
                << meanPlateThickness << "\n";
            } // _forcesFile.is_open()
            
            
            /* - - - - -
             GEOMETRY
             - - - - -*/
            if(_geometryFile.is_open())
            {
                double Lmin = numeric_limits<double>::infinity();
                double Lmax = 0.0;
                double Lrms = 0.0;
                for(unsigned int i=0; i<_plates.size(); i++)
                {
                    double lengthKm = _plates[i]->computeLm() * 1.0E-3;
                    Lmin = min(Lmin, lengthKm);
                    Lmax = max(Lmax, lengthKm);
                    Lrms += pow(lengthKm, 2);
                }
                Lrms = sqrt(Lrms / _plates.size());
                
                double ContLmin = 0.0;
                double ContLmax = 0.0;
                double ContLrms = 0.0;
                unsigned int nContPieces = 0;

                if(_continents.size() > 0)
                {
                    ContLmin = numeric_limits<double>::infinity();
                    for(unsigned int i=0; i<_continents.size(); i++)
                    {
                        double ContLengthKm = _continents[i]->getLm() * 1.0E-3;
                        ContLmin = min(ContLmin, ContLengthKm);
                        ContLmax = max(ContLmax, ContLengthKm);
                        ContLrms += pow(ContLengthKm, 2);
                        
                        unsigned int nCollisions = 0;
                        for (unsigned int j=0; j<_markers.size(); j++)
                        {
                            if (_markers[j]->isA("Collision") && _markers[j]->getContinent() == _continents[i])
                            {
                                nCollisions ++;
                            }
                        }
                        nContPieces += nCollisions + 1;
                    }
                    ContLrms = sqrt(ContLrms / _continents.size());
                } // _continents.size() > 0
                
                double seafloorProduction = computeSeafloorProduction();
                double slabFlux = computeSlabFlux();
                double slabFluxL = computeSlabFlux(true); 
				
                _geometryFile << setprecision(3) << fixed << _ageMa << "\t"
                << _T << "\t"
                << Lrms << "\t"
                << Lmin << "\t"
                << Lmax << "\t"
                << ContLrms << "\t"
                << ContLmin << "\t"
                << ContLmax << "\t"
                << seafloorProduction << "\t"
                << slabFlux << "\t"
                << slabFluxL << "\t"
                << nContPieces << "\n";
            } // _geometryFile.is_open()
            
            
            /* - - - - - - - - - - - - - - - - - - -
             radioactive heat + continental growth
             - - - - - - - - - - - - - - - - - - - */
            if(_radioContFile.is_open())
            {
                updateReferenceConcentrations();
                double Hprim = computeRadioHeatTW(_primitive);
                double Hdepl = computeRadioHeatTW(_depleted);
                _radioContFile << setprecision(3) << fixed  << _ageMa << "\t"
                << _T << "\t"
                << setprecision(6) << scientific << _surfTotOcean * 1.0E-6 << "\t"
                << _surfTotContinent * 1.0E-6 << "\t"
                << setprecision(6) << fixed << _cont_surfRatio << "\t"
                << setprecision(5) << _radioHeatTW << "\t"
                << Hprim << "\t"
                << Hdepl << "\n";
            } // _radioContFile.is_open()
            
            /* - - - - - -
             TEST : rho
             - - - - - - */
            if(_testFile.is_open())
            {
                _testFile << setprecision(3) << fixed << _ageMa << "\t"
                << _T << "\t"
                << _T_plate << "\t"
                << _rho_um << "\t"
                << _rho_pl << "\n";
            }
            
            
            
            /* - - - - - - - - - - - -
             sea level - bathymetry
             - - - - - - - - - - - - -*/
            if(_bathymetryFile.is_open())
            {
                computeSeaLevel();
                
                _bathymetryFile << setprecision(3) << fixed  << _ageMa << "\t"
                << _T << "\t"
                << setprecision(4) << scientific << _VbelowRidge_HSCM << "\t"
                << setprecision(3) << fixed << _ridgeDepth_HSCM << "\t"
                << _sealevel_HSCM << "\t"
                << setprecision(4) << scientific << _VbelowRidge_PM95 << "\t"
                << setprecision(3) << fixed << _ridgeDepth_PM95 << "\t"
                << _sealevel_PM95 << "\t"
                << setprecision(4) << scientific << _VbelowRidge_PM125 << "\t"
                << setprecision(3) << fixed << _ridgeDepth_PM125 << "\t"
                << _sealevel_PM125
                << "\n";
            }
            
            
            /* - - - - - - - - - - - -
             continents
             - - - - - - - - - - - - */
            if(_write_continents)
            {
                _continentsFile << setprecision(3) << fixed << _ageMa
                << "\t" << _continents.size() << "\t";
                for(unsigned int i=0; i<_continents.size(); i++)
                {
                    _continentsFile << ::idSuffixe(_continents[i]->getId()) << "\t"
                    << setprecision(2) << fixed << _continents[i]->getLength() << "\t"
                    << setprecision(5) << fixed << _continents[i]->getPosition() << "\t"
                    << _continents[i]->getU() << "\t"
                    << setprecision(2) << fixed
                    << _continents[i]->getWarmingZone()->getLifetime() << "\t"
                    << _continents[i]->getWarmingZone()->getWarming() << "\t"
                    << setprecision(4) << fixed
                    << _continents[i]->getWarmingZone()->getF()*1.0E-12 << "\t";
                }
                _continentsFile << endl;
                
                
                for(unsigned int i=0; i<_continents.size(); i++)
                {
                    _continents[i]->updateFile();
                }
            } // _write_continents
            
            
            /* - - - - - - - - -
             plates
             - - - - - - - - - - */
            if(_platesFile.is_open())
            {
                _platesFile << setprecision(3) << fixed << _ageMa
                << "\t" << _plates.size() << "\t";
                for(unsigned int i=0; i<_plates.size(); i++)
                {
                    _platesFile << _plates[i]->getRightElement()->getShortId() << "\t"
                    << _plates[i]->getLeftElement()->getShortId() << "\t"
                    << setprecision(5) << fixed << _plates[i]->getU() << "\t"
                    << setprecision(2) << fixed << _plates[i]->computeLength() << "\t"
                    << _plates[i]->computeContinentalLength() << "\t"
                    << _plates[i]->getContinents().size() << "\t";
                }
                _platesFile << endl;
            }
            
            /* - - - - - - -
             geoElements
             - - - - - - - */
            if(_elementsFile.is_open())
            {
                vector<GeoElement*> elements = findInterfaces();
                _elementsFile << setprecision(3) << fixed << _ageMa
                << "\t" << elements.size() << "\t";
                for(unsigned int i=0; i<elements.size(); i++)
                {
                    _elementsFile << elements[i]->getShortId() << "\t"
                    << setprecision(3) << fixed << elements[i]->getPosition() << "\t";
                }
                _elementsFile << endl;
            }
            
            /* - - - - - - - -
             markers
             - - - - - - - - - */
            vector<Collision*> collisions;
            vector<ActiveMargin*> margins;
            for (unsigned int i=0; i<_markers.size(); i++)
            {
                if (_markers[i]->isA("Collision"))
                    collisions.push_back((Collision*)_markers[i]);
                if (_markers[i]->isA("ActiveMargin"))
                    margins.push_back((ActiveMargin*)_markers[i]);
            }
            
            if (_collisionsFile.is_open())
            {
                _collisionsFile << setprecision(3) << fixed << _ageMa
                << "\t" << collisions.size() << "\t";
                for (unsigned int i=0; i<collisions.size(); i++)
                {
                    _collisionsFile << collisions[i]->getShortId() << "\t"
                    << setprecision(3) << fixed << collisions[i]->getPosition()
                    << "\t" << collisions[i]->getAgeMa() << "\t";
                }
                _collisionsFile << endl;
            }
            
            if (_activeMarginsFile.is_open())
            {
                _activeMarginsFile << setprecision(3) << fixed << _ageMa
                << "\t" << margins.size() << "\t";
                for (unsigned int i=0; i<margins.size(); i++)
                {
                    _activeMarginsFile << margins[i]->getShortId() << "\t"
					<< margins[i]->getStatus() << "\t"
                    << setprecision(3) << fixed << margins[i]->getPosition()
                    << "\t" << margins[i]->getAgeMa() << "\t";
                }
                _activeMarginsFile << endl;
            }
        } // all
        
    } // if(_logFile.is_open())
    
}

/*
 =======================
 Events
 =======================
 */
// --------------------------------------------------
vector<Earth::Contact>
Earth::findCollisions(double& dtMin)
{
    vector<Contact> contacts;
    
    for(unsigned int i=0; i<_plates.size(); i++)
    {
        vector<PlateSection*> sections = _plates[i]->accessSections();
        for(unsigned int j=0; j<sections.size(); j++)
        {
            
            // timestep dtColl to have right- and left-Element colliding:
            double URight = cm_to_deg(sections[j]->getRightElement()->getU());
            double ULeft = cm_to_deg(sections[j]->getLeftElement()->getU());
    
            if(::isStrictlyLess(ULeft, URight, Earth::tinyDouble))
            {
                double dtColl = ::computeRightLeftDistance(sections[j]->getRightElement()->getPosition(),
                                                           sections[j]->getLeftElement()->getPosition())
                / ( URight - ULeft );
                
                // collision:
                if (::isStrictlyLess(dtColl, Earth::zero))
                {
                    cerr << "Warning: found dtColl < 0.0 for rightElement = "
                    << sections[j]->getRightElement()->getShortId()
                    << " at pos = " << sections[j]->getRightElement()->getPosition()
                    << "\t and leftElement = " << sections[j]->getLeftElement()->getShortId()
                    << " at pos = " << sections[j]->getLeftElement()->getPosition()
                    << endl;
                }
                
                if(dtColl <= _dt)
                {
                    Contact coll;
                    coll.dt = dtColl;
                    coll.rightElement = sections[j]->getRightElement();
                    coll.leftElement = sections[j]->getLeftElement();
                    contacts.push_back(coll);
                }
            } // URight != ULeft
        } // sections[j]
    } // _plates[i]
    
    if(!contacts.empty())
    { /* keep only collisions happening for the smallest dt */
        vector<Contact>::iterator it = contacts.begin();
        dtMin = (*it).dt;
        while(true)
        { // put the collisions with smallest dt first
            it ++;
            if(it == contacts.end())
                break;
            if((*it).dt <= dtMin)
            {
                dtMin = (*it).dt;
                Contact collToKeep = (*it);
                contacts.erase(it);
                contacts.insert(contacts.begin(), collToKeep);
            }
        }
        // remove collisions with dt > dtMin
        it = contacts.begin();
        while(true)
        {
            it++;
            if(it == contacts.end())
                break;
            if((*it).dt > dtMin)
            {
                contacts.erase(it);
                it --;
            }
        }
        if (contacts.empty())
        {
            cerr << "Warning: some possible collisions were found but vector contacts is empty.\n";
        }
    }
    
    return contacts;
}
// --------------------------------------------------
void
Earth::treatCollisions(vector<Contact> collisions)
{
    for(unsigned int i=0; i<collisions.size(); i++)
    {
        GeoElement* rightElement = collisions[i].rightElement;
        GeoElement* leftElement  = collisions[i].leftElement;
        
        ostringstream oss;
        oss << "\t -> " << setw(8) << setprecision(2) << fixed << getTimeMyr()
        << " Myr \t" << setw(8) << _ageMa
        << " Ma \t" << setw(20) << "-collision: "
        << setw(7) << rightElement->getShortId() << "\t"
        << setw(7) << leftElement->getShortId() << "\t"
        << setw(10) << rightElement->getPosition() << "\t"
        << setw(10) << leftElement->getPosition() << "\n";
        if(_outputFile.is_open())
            _outputFile << oss.str() << flush;
        
        if(_showEvents)
            cout << oss.str();
        
        if(rightElement->isA("RightSubduction"))
        {
            if(leftElement->isA("Staple"))
            {
                if(((Staple*)leftElement)->getLeftNeighborExtremity())
                {
                    _collisionContinentSubduction(rightElement, leftElement, RIGHT);
                }
                else
                { // no continent with the staple: staple is simply subducted
                    rightElement->setLeftAges(leftElement->getLeftAges());
                    rightElement->setLeftAge(leftElement->getLeftAge());
                    removeElement(leftElement);
                    clearPlates();
                    fillStructure();
                }
            } //leftElement->isA("Staple")
            else if(leftElement->isA("Ridge"))
            { // remove subduction and ridge, and create a staple
                Staple* staple = createStaple(rightElement->getPosition());
                staple->setRightAges(rightElement->getRightAges());
                staple->setRightAge(rightElement->getRightAge());
                staple->setLeftAges(leftElement->getLeftAges());
                staple->setLeftAge(leftElement->getLeftAge());
                
                staple->setRightNeighborExtremity(((Subduction*)rightElement)->getRightNeighborExtremity());
                if(staple->getRightNeighborExtremity())
                    staple->getRightNeighborExtremity()->setNeighbor(staple);
                
                // marker:
				string status = "other";
                changeActiveMargin((Subduction*)rightElement, status);
                
                // remove subduction and ridge:
                removeElement(leftElement);
                removeElement(rightElement);
                
                clearPlates();
                fillStructure();
            } //leftElement->isA("Ridge")
            else if(leftElement->isA("ContinentExtremity"))
            {
                _collisionContinentSubduction(rightElement, leftElement, RIGHT);
            } //leftElement->isA("ContinentExtremity")
            else if(leftElement->isA("LeftSubduction"))
            { // remove the two subductions and make a staple
                Staple* staple = createStaple(rightElement->getPosition());
                staple->setRightAges(rightElement->getRightAges());
                staple->setRightAge(rightElement->getRightAge());
                staple->setLeftAges(rightElement->getLeftAges());
                staple->setLeftAge(rightElement->getLeftAge());
                
                staple->setRightNeighborExtremity(((Subduction*)rightElement)->getRightNeighborExtremity());
                if(staple->getRightNeighborExtremity())
                    staple->getRightNeighborExtremity()->setNeighbor(staple);
                staple->setLeftNeighborExtremity(((Subduction*)leftElement)->getLeftNeighborExtremity());
                if(staple->getLeftNeighborExtremity())
                    staple->getLeftNeighborExtremity()->setNeighbor(staple);
                
                // marker:
				string status = "other";
                changeActiveMargin((Subduction*)rightElement, status);
                changeActiveMargin((Subduction*)leftElement, status);
            
                // remove the two subductions:
                removeElement(leftElement);
                removeElement(rightElement);
                clearPlates();
                fillStructure();
            }
            else
            {
                cerr << "Collision between rightSubduction and "
                << leftElement->getClass() << endl;
                cerr << "CASE NOT IMPLEMENTED" << endl;
                exit(EXIT_FAILURE);
            }
        } // rightElement->isA("RightSubduction")
        else if(leftElement->isA("LeftSubduction"))
        {
            if(rightElement->isA("Staple"))
            {
                if(((Staple*)rightElement)->getRightNeighborExtremity())
                {
                    _collisionContinentSubduction(rightElement, leftElement, LEFT);
                }
                else
                {
                    leftElement->setRightAges(rightElement->getRightAges());
                    leftElement->setRightAge(rightElement->getRightAge());
                    
                    removeElement(rightElement);
                    clearPlates();
                    fillStructure();
                }
            } // rightElement->isA("Staple")
            else if(rightElement->isA("Ridge"))
            {
                Staple* staple = createStaple(leftElement->getPosition());
                staple->setRightAges(rightElement->getRightAges());
                staple->setRightAge(rightElement->getRightAge());
                staple->setLeftAges(leftElement->getLeftAges());
                staple->setLeftAge(leftElement->getLeftAge());
                
                staple->setLeftNeighborExtremity(((Subduction*)leftElement)->getLeftNeighborExtremity());
                if(staple->getLeftNeighborExtremity())
                    staple->getLeftNeighborExtremity()->setNeighbor(staple);
                
                // marker:
				string status = "other";
                changeActiveMargin((Subduction*)leftElement, status);
                
                // remove ridge and subduction:
                removeElement(rightElement);
                removeElement(leftElement);
                clearPlates();
                fillStructure();
            } // rightElement->isA("Ridge")
            else if(rightElement->isA("ContinentExtremity"))
            {
                _collisionContinentSubduction(rightElement, leftElement, LEFT);
            } // rightElement->isA("ContinentExtremity")
            else
            {
                cerr << "Collision between leftSubduction and "
                << rightElement->getClass() << endl;
                cerr << "CASE NOT IMPLEMENTED" << endl;
                exit(EXIT_FAILURE);
            }
        } // leftElement->isA("LeftSubduction")
        else
        { // no active subduction on border of section:
            cerr << "Collision without a subduction" << endl;
            cerr << "CASE NOT IMPLEMENTED" << endl;
            
            // TEST ////////////////////////////////////////
            cout << endl;
            
            cout << "Right Element:"
            << "\t" << rightElement->getId()
            << "\t" << rightElement->getPosition()
            << "\t" << rightElement->getU()
            << endl;
            
            cout << "Left Element:"
            << "\t" << leftElement->getId()
            << "\t" << leftElement->getPosition()
            << "\t" << leftElement->getU()
            << endl;
            ////////////////////////////////////////////////
            
            exit(EXIT_FAILURE);
        }
    } // collisions[i]
}
// --------------------------------------------------
void
Earth::correctAges()
{ /* after a collision, because of numeric error, or after 
   opening a continent, there can be
   a tiny offset between the positions of ages and elements.
   Here: make sure all is exactly equal. */
    
    for(unsigned int i=0; i<_plates.size(); i++)
    {
        vector<PlateSection*>& sections = _plates[i]->accessSections();
        for(unsigned int j=0; j<sections.size(); j++)
        {
            sections[j]->correctAges();
        }
    }
}
// --------------------------------------------------
void
Earth::_collisionContinentSubduction(GeoElement* rightElement,
                                     GeoElement* leftElement,
                                     Direction direction)
{
    switch(direction)
    {
        case RIGHT :
        {
            /* rightElement is a rightSubduction and leftElement is
             either a staple touching a continent or
             a RightContinentExtremity */
            Continent* leftContinent = NULL;
            if(leftElement->isA("Staple"))
                leftContinent =
                ((Staple*)leftElement)->getLeftNeighborExtremity()->getContinent();
            else if(leftElement->isA("RightContinentExtremity"))
                leftContinent =
                ((ContinentExtremity*)leftElement)->getContinent();
            // just for checking
            if(!leftContinent)
            {
                cerr << "ERROR in _collisionContinentSubduction\n"
                << "\t leftContinent not found"
                << endl;
                exit(EXIT_FAILURE);
            }
            
            if(((Subduction*)rightElement)->getRightNeighborExtremity())
            {
                /* ========================================
                 collision continent/continent
                 ======================================== */
                Continent* rightContinent =
                ((Subduction*)rightElement)->getRightNeighborExtremity()->getContinent();
                                
                /* create a new collision (before removing leftContinent, to know the position) */
				Collision* collision = new Collision(this,
													 leftContinent->getPosition());
				changeActiveMargin((Subduction*)rightElement, collision->getShortId());				
                changeActiveMargin(rightContinent->getLeftExtremity(), collision->getShortId());
                changeActiveMargin(leftContinent->getRightExtremity(), collision->getShortId());
                
				/* remove staple and subduction in between */
                removeElement(rightElement);
                if(leftElement->isA("Staple"))
                    removeElement(leftElement);
								
                /* create new continent and delete the two old ones */
                Continent* newContinent = _collisionContinentContinent(rightContinent, leftContinent);
				
				/* markers*/
				collision->setContinent(newContinent);
				_markers.push_back(collision);
                
            }
            else
            {
                /* ==========================================
                 continent arrives at the trench -> staple
                 ========================================== */
                if(leftElement->isA("Staple"))
                {
                    Staple* staple = (Staple*)leftElement;
                    staple->setRightAges(rightElement->getRightAges());
                    staple->setRightAge(rightElement->getRightAge());
                }
                else if(leftElement->isA("RightContinentExtremity"))
                { // create a staple
                    Staple* staple = createStaple(leftElement->getPosition());
                    staple->setRightAges(rightElement->getRightAges());
                    staple->setRightAge(rightElement->getRightAge()); // ocean section
                    staple->setLeftAges(leftElement->getLeftAges());
                    staple->setLeftAge(leftElement->getLeftAge()); // continent section
                    
                    staple->setLeftNeighborExtremity((RightContinentExtremity*)leftElement);
                    ((RightContinentExtremity*)leftElement)->setNeighbor(staple);
                }
                else
                {
                    cerr << "ERROR in _collisionContinentSubduction\n"
                    << "\t leftElement is neither a staple nor a RCE"
                    << endl;
                    exit(EXIT_FAILURE);
                }
                
                // remove the subduction
                removeElement(rightElement);
                clearPlates();
                fillStructure();
            }
            
            break; // case RIGHT
        }
            
        case LEFT :
        {
            /* leftElement is a leftSubduction and rightElement is
             either a staple touching a continent or
             a leftContinentExtremity */
            Continent* rightContinent = NULL;
            if(rightElement->isA("Staple"))
                rightContinent =
                ((Staple*)rightElement)->getRightNeighborExtremity()->getContinent();
            else if(rightElement->isA("LeftContinentExtremity"))
                rightContinent =
                ((ContinentExtremity*)rightElement)->getContinent();
            
            // just for checking
            if(!rightContinent)
            {
                cerr << "ERROR in _collisionContinentSubduction\n"
                << "\t rightContinent not found"
                << endl;
                exit(EXIT_FAILURE);
            }
            
            if(((Subduction*)leftElement)->getLeftNeighborExtremity())
            {
                /* ========================================
                 collision continent/continent
                 ======================================== */
                Continent* leftContinent =
                ((Subduction*)leftElement)->getLeftNeighborExtremity()->getContinent();
                
				/* create a new collision (before removing leftContinent, to know the position) */
				Collision* collision = new Collision(this,
													 leftContinent->getPosition());
                changeActiveMargin((Subduction*)leftElement, collision->getShortId());
                changeActiveMargin(leftContinent->getRightExtremity(), collision->getShortId());
                changeActiveMargin(rightContinent->getLeftExtremity(), collision->getShortId());
                
                /*remove staple and subduction in between */
                removeElement(leftElement);
                if(rightElement->isA("Staple"))
                    removeElement(rightElement);
                
                /* create new continent and delete the two old ones */
				Continent* newContinent = _collisionContinentContinent(rightContinent, leftContinent);
				
				/* markers*/
				collision->setContinent(newContinent);
				_markers.push_back(collision);
                
            }
            else
            {
                /* =========================================
                 continent arrives at the trench -> staple
                 ========================================= */
                if(rightElement->isA("Staple"))
                {
                    Staple* staple = (Staple*)rightElement;
                    staple->setLeftAges(leftElement->getLeftAges());
                    staple->setLeftAge(leftElement->getLeftAge());
                }
                else if(rightElement->isA("LeftContinentExtremity"))
                { // create a staple
                    Staple* staple = createStaple(rightElement->getPosition());
                    staple->setRightAges(rightElement->getRightAges());
                    staple->setRightAge(rightElement->getRightAge());
                    staple->setLeftAges(leftElement->getLeftAges());
                    staple->setLeftAge(leftElement->getLeftAge());
                    
                    staple->setRightNeighborExtremity((LeftContinentExtremity*)rightElement);
                    ((LeftContinentExtremity*)rightElement)->setNeighbor(staple);
                }
                else
                {
                    cerr << "ERROR in _collisionContinentSubduction\n"
                    << "\t rightElement is neither a staple nor a LCE"
                    << endl;
                    exit(EXIT_FAILURE);
                }
                
                // remove the subduction
                removeElement(leftElement);
                clearPlates();
                fillStructure();
            }
            
            break; // case LEFT
        }
        default :
            cerr << "ERROR: _collisionContinentSubduction() used without appropriate DIRECTION"
            << endl;
            exit(EXIT_FAILURE);
            
    }
}
// --------------------------------------------------
Continent*
Earth::_collisionContinentContinent(Continent* rightContinent,
                                    Continent* leftContinent)
{
    string strColl = "collision";
    Continent* newContinent =
    createContinent(rightContinent->getPosition(),
                    rightContinent->getLength() +
                    leftContinent->getLength(),
                    strColl);
    
    /* switch continents for all markers belonging
     to leftContinent and rightContinent */
    _markersChangeContinents(rightContinent, newContinent);
    _markersChangeContinents(leftContinent, newContinent);
        
    /* events file */
    if(_eventsFile.is_open())
    {
        _eventsFile << _ageMa << "\t"
        << _T << "\t"
        << "A" << endl;
    }
    /* exchange of extremities (it includes sections' ages) */
    newContinent->setRightExtremity(rightContinent->getRightExtremity());
    newContinent->getRightExtremity()->setContinent(newContinent);
    
    newContinent->setLeftExtremity(leftContinent->getLeftExtremity());
    newContinent->getLeftExtremity()->setPosition(newContinent->getLeftPosition());
    newContinent->getLeftExtremity()->setContinent(newContinent);
    
    /* get the ages of the collided continents */
    vector<Age> rightAges = rightContinent->getRightExtremity()->getLeftAges();
    vector<Age> leftAges = leftContinent->getRightExtremity()->getLeftAges();
    
    newContinent->gatherAges(rightAges, leftAges);
    
    /* if leftContinent has a staple/subd on the border:
     there can be a slight change of position */
    if (newContinent->getLeftExtremity()->getNeighbor())
    {
        newContinent->getLeftExtremity()->getNeighbor()->setPosition(newContinent->getLeftPosition());
        newContinent->getLeftExtremity()->getNeighbor()->setRightAges(newContinent->getLeftExtremity()->getRightAges());
        newContinent->getLeftExtremity()->getNeighbor()->setRightAge(newContinent->getLeftExtremity()->getRightAge());
    }
    if (newContinent->getRightExtremity()->getNeighbor())
    {
        newContinent->getRightExtremity()->getNeighbor()->setPosition(newContinent->getRightPosition());
        newContinent->getRightExtremity()->getNeighbor()->setLeftAges(newContinent->getRightExtremity()->getLeftAges());
        newContinent->getRightExtremity()->getNeighbor()->setLeftAge(newContinent->getRightExtremity()->getLeftAge());
    }
    
    
    /* delete internal extremities */
    delete leftContinent->getRightExtremity();
    delete rightContinent->getLeftExtremity();
    
    /* subcontinental warming */
    double warming = ( leftContinent->getWarmingZone()->getWarming() * leftContinent->getLength() +
                      rightContinent->getWarmingZone()->getWarming() * rightContinent->getLength()
                      ) / newContinent->getLength();
    
    newContinent->getWarmingZone()->computeLifetime0(warming); // equivalent lifetime of new cont.
    
    // --- Remove continents
    removeContinent(leftContinent);
    removeContinent(rightContinent);
    
    clearPlates();
    fillStructure();
    
	return newContinent;
}
// --------------------------------------------------
void
Earth::_checkDiving()
{  
    if(_fixed_configuration)
        return;
    
    bool diving = false;
    
    if(_initSubduction.staples)
        diving = _divingStaple();
    
    if(!diving && _initSubduction.continents)
        diving = _divingContinentExtremity();
    
    if(!diving && _initSubduction.upperPlates)
        diving = _divingUpperPlate();
    
    
    if(diving)
        _checkDiving();
}
// --------------------------------------------------
bool
Earth::_divingStaple()
{
    bool diving = false;
    double tauCritical = _tau_sub;
    
    for(unsigned int i=0; i<_elements.size(); i++)
    {
        if(_elements[i]->isA("Staple"))
        {
            Staple* staple = (Staple*)_elements[i];
            Subduction* subduction = NULL;
            
            // RIGHT SUBDUCTION - - - - -
            /* conditions: limit age */
            if(_initSubduction.randomTauSub)
                tauCritical = _randomAge(_tau_sub,
                                         _initSubduction.tauSub_randomNoise);
            
            if(staple->getLeftAge().isCritical(tauCritical))
            {
                subduction = new RightSubduction(this, staple->getPosition());
                // check that the subduction is in the good direction:
                diving = _realDivingSubduction(subduction, staple, RIGHT);
            } // limit age for right subduction
            
            // LEFT SUBDUCTION - - - - - - - -
            if(!diving)
            {
                if(_initSubduction.randomTauSub)
                    tauCritical = _randomAge(_tau_sub,
                                             _initSubduction.tauSub_randomNoise);
                
                if(staple->getRightAge().isCritical(tauCritical))
                {
                    subduction = new LeftSubduction(this, staple->getPosition());
                    // check that the subduction is in the good direction:
                    diving = _realDivingSubduction(subduction, staple, LEFT);
                } // limit age for left subduction
            } // !diving (for left subduction)
            
            // really create the subduction now - - - - - - - - -
            if(diving)
            {
                _setNewSubduction(subduction, staple);
                
                clearPlates();
                
                insertElement(subduction);
                removeElement(staple);
                
                fillStructure();
                break;
                
            } // actually diving
            else // !diving
            { /* decrement the number of subduction to not count
               the new one here that was created just for testing */
                if(subduction)
                {
                    delete subduction;
                    _subductionCounter --;
                }
            }
        } // staple
    } // elements[i]
    
    return diving;
}

// --------------------------------------------------
bool
Earth::_divingContinentExtremity()
{
    bool diving = false;
    double tauCritical = _tau_sub;
    
    for(unsigned int i=0; i<_continents.size(); i++)
    {
        Subduction* subduction = NULL;
        LeftContinentExtremity* leftExtremity = NULL;
        RightContinentExtremity* rightExtremity = NULL;
        
        // RIGHT SUBDUCTION ON LCE - - - - - - - - - -
        if(!_continents[i]->getLeftExtremity()->getNeighbor() ||
           _continents[i]->getLeftExtremity()->getNeighbor()->isA("Staple"))
        {
            leftExtremity = _continents[i]->getLeftExtremity();
            
            // condition on age
            if(_initSubduction.randomTauSub)
                tauCritical = _randomAge(_tau_sub,
                                         _initSubduction.tauSub_randomNoise);
            
            
            if(leftExtremity->getLeftAge().isCritical(tauCritical))
            {
                subduction = new RightSubduction(this,
                                                 leftExtremity->getPosition());
                // check that the subduction is in the good direction:
                diving = _realDivingSubduction(subduction,
                                               leftExtremity, RIGHT);
            } // limit age for right subduction
        } // left extremity has no neighbor
        
        if(!diving)
        { // LEFT SUBDUCTION ON RCE - - - - - - - - -
            if(!_continents[i]->getRightExtremity()->getNeighbor() ||
               _continents[i]->getRightExtremity()->getNeighbor()->isA("Staple"))
            {
                rightExtremity = _continents[i]->getRightExtremity();
                
                // condition on age
                if(_initSubduction.randomTauSub)
                    tauCritical = _randomAge(_tau_sub,
                                             _initSubduction.tauSub_randomNoise);
                
                if(rightExtremity->getRightAge().isCritical(tauCritical))
                {
                    subduction =
                    new LeftSubduction(this,
                                       rightExtremity->getPosition());
                    // check that the subduction is in the good direction:
                    diving = _realDivingSubduction(subduction,
                                                   rightExtremity, LEFT);
                } // limit age for left subduction
            } // right extremity has no neighbor
        } // !diving (for right extremity -> left subduction)
        
        // now actually create the subduction - - - - - - - -
        if(diving)
        {
            if(subduction->isA("RightSubduction"))
            {
                _setNewRightSubduction(subduction, leftExtremity);
            }
            else //leftSubduction
            {
                _setNewLeftSubduction(subduction, rightExtremity);
            }
            
            clearPlates();
            
            insertElement(subduction);
            fillStructure();
            break;
        } // actually diving
        else // !diving
        { /* decrement the number of subduction to not count
           the new one here that was created just for testing */
            if(subduction)
            {
                delete subduction;
                _subductionCounter --;
            }
        }
    } // continents[i]
    
    return diving;
}
// -------------------------------------------------------
bool
Earth::_divingUpperPlate()
{
    bool diving = false;
    
    for(unsigned int i=0; i<_elements.size(); i++)
    {
        
        Subduction* subduction = NULL;
        Subduction* newSubduction = NULL;
        
        if(_elements[i]->isA("RightSubduction"))
        {
            subduction = (RightSubduction*)_elements[i];
            
            if(subduction->getRightAge().isOceanic())
            { // upper plate is oceanic
                
                bool ageReached = false;
                double rightAge = subduction->getRightAge().getValue();
                
                if(rightAge >= _initSubduction.upperPlateMinAge)
                {
                    if(_initSubduction.ageRatioCriterion)
                    {
                        double leftAge = subduction->getLeftAge().getValue();
                        if( rightAge > _initSubduction.ratioAgeValue * leftAge
                           && leftAge > 0.0 )
                            ageReached = true;
                    }
                    else
                    { // criterion is on _tau_sub
                        double tauCritical = _tau_sub;
                        if(_initSubduction.randomTauSub)
                            tauCritical =
                            _randomAge(_tau_sub,
                                       _initSubduction.tauSub_randomNoise);
                        if( rightAge > tauCritical )
                            ageReached = true;
                    }
                }
                
                if(ageReached)
                {
                    /*  needs to have a position slightly different
                     from subduction so that the next interface on the right
                     is not the old subduction (for _realDivingSubduction() */
                    double position =
                    offsetPosition(subduction,
                                   2.0 * Earth::creationDistance, RIGHT);
                    // factor 2 is to let room for a new ridge if necessary //
                    
                    newSubduction =
                    new LeftSubduction(this,
                                       position);
                    
                    diving =
                    _realDivingSubduction(newSubduction, subduction,
                                          LEFT, false);
                    
                } // ageReached
            } // upper oceanic plate
            
        } // rightSubduction
        
        if(!diving)
        {
            if(_elements[i]->isA("LeftSubduction"))
            {
                subduction = (LeftSubduction*)_elements[i];
                
                if(subduction->getLeftAge().isOceanic())
                { // upper plate is oceanic
                    
                    bool ageReached = false;
                    double leftAge = subduction->getLeftAge().getValue();
                    
                    if(leftAge >= _initSubduction.upperPlateMinAge)
                    {
                        if(_initSubduction.ageRatioCriterion)
                        {
                            double rightAge = subduction->getRightAge().getValue();
                            if( leftAge > _initSubduction.ratioAgeValue * rightAge
                               && rightAge > 0.0 )
                                ageReached = true;
                        }
                        else
                        { // criterion is on _tau_sub
                            double tauCritical = _tau_sub;
                            if(_initSubduction.randomTauSub)
                                tauCritical =
                                _randomAge(_tau_sub,
                                           _initSubduction.tauSub_randomNoise);
                            if( leftAge > tauCritical )
                                ageReached = true;
                        }
                    }
                    
                    if(ageReached)
                    {
                        double position =
                        offsetPosition(subduction,
                                       2.0 * Earth::creationDistance, LEFT);
                        newSubduction =
                        new RightSubduction(this,
                                            position);
                        
                        diving =
                        _realDivingSubduction(newSubduction, subduction,
                                              RIGHT, false);
                        
                    } // ageReached
                    
                } // upper oceanic plate
            } // leftSubduction
        } // !diving
        
        if(diving) // now actually create the new subduction - - - - - - -
        {
            if(_initSubduction.reverse)
            { // remove the old subduction
                newSubduction->setPosition(subduction->getPosition());
                
                _setNewSubduction(newSubduction, subduction);
                
                removeElement(subduction);
            }
            else // not reverse
            { /* create a new ridge
               in between the two subductions */
                _setNewSubduction(newSubduction, subduction);
                
                if(newSubduction->isA("RightSubduction"))
                {
                    double posRidge =
                    ::getMiddle(subduction->getPosition(),
                                newSubduction->getPosition());
                    
                    createRidge(posRidge);
                    newSubduction->clearRightAges();
                    subduction->clearLeftAges();
                }
                else // LEFT newSubduction
                {
                    double posRidge =
                    ::getMiddle(newSubduction->getPosition(),
                                subduction->getPosition());
                    
                    createRidge(posRidge);
                    newSubduction->clearLeftAges();
                    subduction->clearRightAges();
                }
                
                ostringstream oss;
                oss << "\t -> " << setw(8) << setprecision(2) << fixed << getTimeMyr()
                << " Myr \t" << setw(8) << _ageMa
                << " Ma \t" << setw(20) << "+creation ridge at: "
                << setw(7) << newSubduction->getShortId()
                << "\t" << setw(7) << newSubduction->getPosition() << "\n";
                if(_outputFile.is_open())
                    _outputFile << oss.str() << flush;
                if(_showEvents)
                    cout << oss.str() << flush;
                
            } // not reverse
            
            clearPlates();
            insertElement(newSubduction);
            fillStructure();
            break;
        }
        else // !diving
        { /* decrement the number of subduction to not count
           the new one here that was created just for testing */
            if(newSubduction)
            {
                delete newSubduction;
                _subductionCounter --;
            }
        }
        
    } // elements[i]
    
    return diving;
}
// -------------------------------------------------------
// --------------------------------------------------
// void
// Earth::_checkReverse()
// { //June 15th 2017 : NOT WORKING SO FAR.........
//   if(!_reverseSubduction)
//     return;

//   bool reverse = false;
//   double tauCritical = _tau_sub;

//   for(unsigned int i=0; i<_elements.size(); i++)
//     {
//       if(_elements[i]->isA("Subduction"))
// 	{
// 	  Subduction* oldSubduction = (Subduction*)_elements[i];
// 	  Subduction* newSubduction = NULL;

// 	  // RS reversing to LS - - - - - -
// 	  if(oldSubduction->isA("RightSubduction")) 
// 	    {
// 	      if(_initSubduction.randomTauSub) 
// 		tauCritical = _randomAge(_tau_sub, _initSubduction.tauSub_randomNoise);

// 	      if(oldSubduction->getRightAge().ocean && oldSubduction->getRightAge().value > tauCritical)
// 		{
// 		  newSubduction = new LeftSubduction(this, oldSubduction->getPosition());
// 		  reverse = _realDivingSubduction(newSubduction, oldSubduction, LEFT, false);
// 		}
// 	    } // RS reversing to LS

// 	  // LS reversing to RS
// 	  if(!reverse)
// 	    {
// 	      if(oldSubduction->isA("LeftSubduction"))
// 		{
// 		  if(_initSubduction.randomTauSub) 
// 		    tauCritical = _randomAge(_tau_sub, _initSubduction.tauSub_randomNoise);

// 		  if(oldSubduction->getLeftAge().ocean && oldSubduction->getLeftAge().value > tauCritical)
// 		    {
// 		      newSubduction = new RightSubduction(this, oldSubduction->getPosition());
// 		      reverse = _realDivingSubduction(newSubduction, oldSubduction, RIGHT, false);
// 		    }
// 		} // oldSubduction->isA("LeftSubduction")
// 	    } // !reverse

// 	  // Really reverse subduction - - - - - - - - 
// 	  if(reverse)
// 	    { /* remove former subduction and create new reversed one, 
// 		 with exactly same neighbors, ages etc. */
// 	      newSubduction->setLeftSectionAges(oldSubduction->getLeftSectionAges());
// 	      newSubduction->setRightSectionAges(oldSubduction->getRightSectionAges());

// 	      ostringstream oss;
// 	      oss << "\t -> " << setw(8) << setprecision(2) << fixed << getTimeMyr()
// 		  << " Myr \t" << setw(8) << _ageMa
// 		  << " Ma \t" << setw(20) << "reverse: " 
// 		  << setw(7) << oldSubduction->getShortId()
// 		  << "\t" << setw(10) << oldSubduction->getPosition() << "\n";

// 	      if(_outputFile.is_open())
// 		_outputFile << oss.str() << flush;
// 	      if(_showEvents)
// 		cout << oss.str() << flush;

// 	      clearPlates();

// 	      insertElement(newSubduction);
// 	      removeElement(oldSubduction);

// 	      fillStructure();
// 	      break;

// 	    } // reverse
// 	  else // !reverse
// 	    { /* decrement the number of subduction to not count 
// 		 the new one here that was created just for testing */
// 	      if(newSubduction)
// 		{
// 		  delete newSubduction; 
// 		  _subductionCounter --;
// 		}
// 	    }	      

// 	} // isA("Subduction")
//     } // elements[i]

// }
// --------------------------------------------------
bool
Earth::_realDivingSubduction(Subduction* subduction, 
                             GeoElement* element,
                             Direction direction,
                             bool checkNearInterface)
{
    /* checkNearInterface: optional (default is true) */
    bool diving = false;
    double limitDistance = 2.0 * Earth::creationDistance;
    double thickness = 0.0;
    double depth = 0.0;
    
    // check that there is not a very close interface
    if(checkNearInterface)
        diving = !(subduction->nearInterface(limitDistance));
    else
        diving = true;
    
    if(!diving)
        return false;
    
    switch(direction)
    {
        case LEFT :
        {
            /* checking if a new left subduction can really start */
            thickness = element->computePlateThickness(RIGHT);
            depth = thickness * _depthFactor;
            subduction->setRightAge(element->getRightAge());
            subduction->setDepth(depth);
            subduction->initThickness(depth, thickness);
            
            /* check that the plate between rightInterface
             and this left subduction is going to go to the left. */
            Interface* rightInterface = subduction->findNextInterface(RIGHT);
            
            if(rightInterface)
                diving = (subduction->getDirectionDrivingForces(RIGHT) == LEFT);
            else // if there's only one plate:
                diving = true;
            break;
        }
            
        case RIGHT :
        {
            /* checking if a new right subduction can really start */
            thickness = element->computePlateThickness(LEFT);
            depth = thickness * _depthFactor;
            subduction->setLeftAge(element->getLeftAge());
            subduction->setDepth(depth);
            subduction->initThickness(depth, thickness);
            
            /* check that the plate between leftInterface
             and this right subduction is going to go to the right. */
            Interface* leftInterface = subduction->findNextInterface(LEFT);
            
            if(leftInterface)
                diving = (subduction->getDirectionDrivingForces(LEFT) == RIGHT);
            else // if there's only one plate:
                diving = true;
            break;
        }
            
        default :
        {
            cerr << "ERROR: Earth::_realDivingSubduction() is used with no direction"
            << endl;
            exit(EXIT_FAILURE);
        }
    }
    
    return diving;
}
// --------------------------------------------------
void
Earth::_checkAlwaysSubduct()
{
    if(!_initSubduction.alwaysOne || !_noSubduction)
        return;
    
    // if _noSubduction:
    //first search for oldest continental margins
    /* warning: fast easy solution here;
     if two points have the same age, the last that was
     found is the one where the subduction starts
     */
    
    if(_continents.size() > 0)
    {
        double oldest = 0.0;
        ContinentExtremity* extremity = NULL;
        ContinentExtremity* oldestExtremity = NULL;
        for(unsigned int i=0; i<_continents.size(); i++)
        {
            if(!_continents[i]->getLeftExtremity()->getNeighbor() ||
               _continents[i]->getLeftExtremity()->getNeighbor()->isA("Staple"))
            {
                extremity = _continents[i]->getLeftExtremity();
                oldest = max(extremity->getLeftAge().getValue(), oldest);
                if(::areEqual(oldest, extremity->getLeftAge().getValue(), tinyDouble))
                    oldestExtremity = extremity;
            }
            
            if(!_continents[i]->getRightExtremity()->getNeighbor() ||
               _continents[i]->getRightExtremity()->getNeighbor()->isA("Staple"))
            {
                extremity = _continents[i]->getRightExtremity();
                oldest = max(extremity->getRightAge().getValue(), oldest);
                if(::areEqual(oldest, extremity->getRightAge().getValue(), tinyDouble))
                    oldestExtremity = extremity;
            }
        } // continents.size()
        
        if(oldestExtremity)
        {
            bool diving = false;
            
            LeftContinentExtremity* leftExtremity = NULL;
            RightContinentExtremity* rightExtremity = NULL;
            Subduction* subduction = NULL;
            if(oldestExtremity->isA("LeftContinentExtremity"))
            {
                leftExtremity = (LeftContinentExtremity*)oldestExtremity;
                subduction = new RightSubduction(this, leftExtremity->getPosition());
                
                // mandatory to init _thickness of subduction
                diving = _realDivingSubduction(subduction,
                                               leftExtremity, RIGHT);
            }
            else
            {
                rightExtremity = (RightContinentExtremity*)oldestExtremity;
                subduction = new LeftSubduction(this,
                                                rightExtremity->getPosition());
                
                // mandatory to init _thickness of subduction
                diving = _realDivingSubduction(subduction,
                                               rightExtremity, LEFT);
            }
            
            
            if(diving)
            {
                if(subduction->isA("RightSubduction"))
                    _setNewRightSubduction(subduction, leftExtremity);
                else
                    _setNewLeftSubduction(subduction, rightExtremity);
                
                clearPlates();
                insertElement(subduction);
                fillStructure();
            }
            else
            {
                if(subduction)
                {
                    delete subduction;
                    _subductionCounter --;
                }
            }
        } // oldestExtremity
    } // continents.size() > 0
    else
    { // if no continent: search for oldest staples
        double oldest = 0.0;
        Staple* staple = NULL;
        Staple* oldestStaple = NULL;
        Direction direction;
        for(unsigned int i=0; i<_elements.size(); i++)
        {
            if(_elements[i]->isA("Staple"))
            {
                staple = (Staple*)_elements[i];
                
                if(!staple->getRightNeighborExtremity())
                {
                    oldest = max(staple->getRightAge().getValue(), oldest);
                    if(::areEqual(staple->getRightAge().getValue(), oldest, tinyDouble))
                    {
                        oldestStaple = staple;
                        direction = RIGHT;
                    }
                }
                if(!staple->getLeftNeighborExtremity())
                {
                    oldest = max(staple->getLeftAge().getValue(), oldest);
                    if(::areEqual(staple->getLeftAge().getValue(), oldest, tinyDouble))
                    {
                        oldestStaple = staple;
                        direction = LEFT;
                    }
                }
            } // elements is a staple
        } // for(elements)
        
        if(oldestStaple)
        {
            bool diving = false;
            Subduction* subduction = NULL;
            
            if(direction == LEFT)
            {
                subduction = new RightSubduction(this, staple->getPosition());
                diving = _realDivingSubduction(subduction, staple, RIGHT);
            }
            else
            {
                subduction = new LeftSubduction(this, staple->getPosition());
                diving = _realDivingSubduction(subduction, staple, LEFT);
            } // direction
            
            if(diving)
            {
                _setNewSubduction(subduction, staple);
                clearPlates();
                insertElement(subduction);
                removeElement(staple);
                
                fillStructure();
            }
            else
            {
                if(subduction)
                {
                    delete subduction;
                    _subductionCounter --;
                }
            } // diving
            
            
        } // oldestStaple
        
    } // No continent
}
// --------------------------------------------------
void
Earth::_setNewRightSubduction(Subduction* subduction, LeftContinentExtremity* leftExtremity)
{
    subduction->setLeftAges(leftExtremity->getLeftAges());
    subduction->setLeftAge(leftExtremity->getLeftAge());
    subduction->setRightAges(leftExtremity->getRightAges());
    subduction->setRightAge(leftExtremity->getRightAge());
    
    subduction->setRightNeighborExtremity(leftExtremity);
    
    if(leftExtremity->getNeighbor()) //can only be a staple so far
        removeElement(leftExtremity->getNeighbor());
    
    leftExtremity->setNeighbor(subduction);
    
    ostringstream oss;
    oss << "\t -> " << setw(8) << setprecision(2) << fixed << getTimeMyr()
    << " Myr \t" << setw(8) << _ageMa
    << " Ma \t" << setw(20) << "+creation: "
    << setw(7) << subduction->getShortId()
    << setw(5) << " at "
    << setw(7) << leftExtremity->getShortId()
    << "\t" << setw(10) << subduction->getPosition()
    << "\n";
    
    if(_outputFile.is_open())
        _outputFile << oss.str() << flush;
    if(_showEvents)
        cout << oss.str() << flush;

    /* create a new activeMargin marker */
    ActiveMargin* margin = new ActiveMargin(this,
                                            subduction->getPosition());
    margin->setSubduction(subduction);
	margin->setContinent(leftExtremity->getContinent());
    margin->setContinentExtremity(leftExtremity);
	string status = "active";
	margin->setStatus(status);
    
    subduction->setActiveMargin(margin);
    leftExtremity->setActiveMargin(margin);
    _markers.push_back(margin);
}
// --------------------------------------------------
void
Earth::_setNewLeftSubduction(Subduction* subduction, RightContinentExtremity* rightExtremity)
{
    subduction->setRightAges(rightExtremity->getRightAges());
    subduction->setRightAge(rightExtremity->getRightAge());
    subduction->setLeftAges(rightExtremity->getLeftAges());
    subduction->setLeftAge(rightExtremity->getLeftAge());
    
    subduction->setLeftNeighborExtremity(rightExtremity);
    
    if(rightExtremity->getNeighbor()) //can only be a staple so far
        removeElement(rightExtremity->getNeighbor());
    
    rightExtremity->setNeighbor(subduction);
    
    ostringstream oss;
    oss << "\t -> " << setw(8) << setprecision(2) << fixed << getTimeMyr()
    << " Myr \t" << setw(8) << _ageMa
    << " Ma \t" << setw(20) << "+creation: "
    << setw(7) << subduction->getShortId()
    << setw(5) << " at "
    << setw(7) << rightExtremity->getShortId()
    << "\t" << setw(10) << subduction->getPosition()
    << "\n";
    
    if(_outputFile.is_open())
        _outputFile << oss.str() << flush;
    if(_showEvents)
        cout << oss.str() << flush;
    
    /* create a new activeMargin marker */
    ActiveMargin* margin = new ActiveMargin(this,
                                            subduction->getPosition());
    margin->setSubduction(subduction);
	margin->setContinent(rightExtremity->getContinent());
    margin->setContinentExtremity(rightExtremity);
	string status = "active";
	margin->setStatus(status);

    subduction->setActiveMargin(margin);
    rightExtremity->setActiveMargin(margin);
    _markers.push_back(margin);
    
    
}
// -------------------------------------------------------
void
Earth::_setNewSubduction(Subduction* subduction, GeoElement* element)
{
    // ages
    subduction->setLeftAges(element->getLeftAges());
    subduction->setLeftAge(element->getLeftAge());
    subduction->setRightAges(element->getRightAges());
    subduction->setRightAge(element->getRightAge());
    if(element->isA("Staple"))
    {
        subduction->setLeftNeighborExtremity(((Staple*)element)->getLeftNeighborExtremity());
        subduction->setRightNeighborExtremity(((Staple*)element)->getRightNeighborExtremity());
    }
    else if(element->isA("Subduction"))
    {
        subduction->setLeftNeighborExtremity(((Subduction*)element)->getLeftNeighborExtremity());
        subduction->setRightNeighborExtremity(((Subduction*)element)->getRightNeighborExtremity());
    }
    else
    {
        cerr << "ERROR: Earth::_setNewSubduction(Subduction*,GeoElement*) is implemented only for GeoElement=Staple or Subduction" << endl;
        exit(EXIT_FAILURE);
    }
    
    /* give neighbor to continent, if element touches continent */
    if(subduction->getRightNeighborExtremity())
        subduction->getRightNeighborExtremity()->setNeighbor(subduction);
    if(subduction->getLeftNeighborExtremity())
        subduction->getLeftNeighborExtremity()->setNeighbor(subduction);
    
    ostringstream oss;
    oss << "\t -> " << setw(8) << setprecision(2) << fixed << getTimeMyr()
    << " Myr \t" << setw(8) << _ageMa
    << " Ma \t" << setw(20) << "+creation: "
    << setw(7) << subduction->getShortId()
    << setw(5) << " at " << setw(7) << element->getShortId();
    
    if(subduction->isA("RightSubduction") &&
       subduction->getRightNeighborExtremity())
        oss << " - " << setw(7) << subduction->getRightNeighborExtremity()->getShortId();
    else if(subduction->isA("LeftSubduction") &&
            subduction->getLeftNeighborExtremity())
        oss << " - " << setw(7) << subduction->getLeftNeighborExtremity()->getShortId();
    oss << "\t" << setw(10) << subduction->getPosition() << "\n";
    
    if(_outputFile.is_open())
        _outputFile << oss.str() << flush;
    if(_showEvents)
        cout << oss.str() << flush;
    
}
// -------------------------------------------------------
void
Earth::_checkOceanOpening()
{
    if(_continents.size() == 0 || _fixed_configuration)
        return;
    
    Continent* continent = NULL;
    bool opening = false;
    
    for(unsigned int i=0; i<_continents.size(); i++)
    {
        _continents[i]->updateBreakable();
        if(_plates.size() > 1)
        {
            if(_continents[i]->checkOpening())
            {
                continent = _continents[i];
                opening = true;
                if(opening)
                    break;
            }
        }
    } // for(unsigned int i=0; i<_continents.size(); i++)
    
    if(opening)
    {
        double breakupPos = _openOcean(continent);  /* returns the position of
                                                     continental breakup */
        
        ostringstream oss;
        oss << "\t -> " << setw(8) << setprecision(2) << fixed << getTimeMyr()
        << " Myr \t" << setw(8) << _ageMa
        << " Ma \t" << setw(20) << " breakup:"
        << setw(10) << breakupPos << "\n";
        
        if(_outputFile.is_open())
            _outputFile << oss.str() << flush;
        if(_showEvents)
            cout << oss.str() << flush;
        
        _checkOceanOpening();
    }
}
// --------------------------------------------------
double
Earth::_openOcean(Continent* continent)
{ 
    double middle = continent->getMiddlePosition(); // in degrees
    double length = continent->getLength(); // in degrees
    
    double breakupPos = middle;
    if(!_middleBreakup)
    { /* opening continents not exactly in their middle,
       but at more or less  "p*length"   */
        breakupPos = _random->nextDouble(middle - _breakupPosition*length,
                                         middle + _breakupPosition*length);
        breakupPos = fmod(breakupPos, 360.0);
        if(::isStrictlyLess(breakupPos, 0.0))
            breakupPos += 360.0;
    }
    
    double leftContLength =
    ::computeRightLeftDistance(breakupPos, continent->getLeftPosition());
    double rightContLength =
    ::computeRightLeftDistance(continent->getRightPosition(), breakupPos);
    
    createRidge(breakupPos);
    
    /* create the two new continents,
     with checking that their outer sides do not collide with
     other elements (criterion on the size of plateSections left and right)
     */
    double openLength = Earth::creationDistance;
    openLength =
    min( openLength,
        0.9*continent->getRightExtremity()->getRightSection()->getLength() );
    openLength =
    min( openLength,
        0.9*continent->getLeftExtremity()->getLeftSection()->getLength() );
    
    double leftContPosition = ::addPosition(breakupPos, openLength);
    double rightContPosition = ::addPosition(continent->getRightPosition(), -openLength);
    
    /* create the new continents */
    string strBreak = "breakup";
    
    Continent* rightContinent = createContinent(rightContPosition,
                                                rightContLength, strBreak);
    Continent* leftContinent  = createContinent(leftContPosition,
                                                leftContLength, strBreak);
    
    /* new extremities: inside */
    createRightContinentExtremity(leftContinent);
    createLeftContinentExtremity(rightContinent);
    
    /* give outside extremities to new continents */
    leftContinent->setLeftExtremity(continent->getLeftExtremity());
    leftContinent->getLeftExtremity()->setContinent(leftContinent);
    leftContinent->getLeftExtremity()->setPosition(leftContinent->getLeftPosition());
    
    rightContinent->setRightExtremity(continent->getRightExtremity());
    rightContinent->getRightExtremity()->setContinent(rightContinent);
    rightContinent->getRightExtremity()->setPosition(rightContPosition);
    
    /* set the ages on the inside borders of the new continents
     and inside the continents */
    PlateSection* continentalSection = continent->getPlateSection();
    vector<Age> rightAges = continentalSection->splitAges(breakupPos,
                                                          openLength, RIGHT);
    vector<Age> leftAges  = continentalSection->splitAges(breakupPos,
                                                          openLength, LEFT);
    
    rightContinent->getRightExtremity()->setLeftAge(rightAges[0]);
    rightContinent->getLeftExtremity()->setRightAge(rightAges[rightAges.size()-1]);
    rightContinent->getRightExtremity()->setLeftAges(rightAges);
    rightContinent->getLeftExtremity()->setRightAges(rightAges);
    
    leftContinent->getRightExtremity()->setLeftAge(leftAges[0]);
    leftContinent->getLeftExtremity()->setRightAge(leftAges[leftAges.size()-1]);
    leftContinent->getRightExtremity()->setLeftAges(leftAges);
    leftContinent->getLeftExtremity()->setRightAges(leftAges);
    
    double tinyAge = 1.0E-6;
    /* for newly created ocean: give tiny non zero age
     to continents inside borders to not have plate
     of age zero everywhere (=jump in Qtot) */
    Age insideAgeLeft = Age(this, tinyAge,
                            leftContinent->getRightExtremity()->getPosition());
    leftContinent->getRightExtremity()->setRightAge(insideAgeLeft);
    Age insideAgeRight = Age(this, tinyAge,
                             rightContinent->getLeftExtremity()->getPosition());
    rightContinent->getLeftExtremity()->setLeftAge(insideAgeRight);
    
    /* if RCE or LCE have neighbor elements:
     change their position and share section ages */
    if(rightContinent->getRightExtremity()->getNeighbor())
    {
        rightContinent->getRightExtremity()->getNeighbor()->setLeftAges(rightAges);
        rightContinent->getRightExtremity()->getNeighbor()->setLeftAge(rightAges[0]);
        rightContinent->getRightExtremity()->getNeighbor()->setPosition(rightContPosition);
    }
    if(leftContinent->getLeftExtremity()->getNeighbor())
    {
        leftContinent->getLeftExtremity()->getNeighbor()->setRightAges(leftAges);
        leftContinent->getLeftExtremity()->getNeighbor()->setRightAge(leftAges[leftAges.size()-1]);
        leftContinent->getLeftExtremity()->getNeighbor()->setPosition(leftContinent->getLeftPosition());
    }
    
    /* markers --------------------- */
    // change continents for existing tectonics markers - - - -
    _markersChangeContinents(continent, rightContinent, leftContinent,
							 breakupPos, openLength);
    /* ------------------ markers */
    
    /* events file */
    if(_eventsFile.is_open())
        _eventsFile << _ageMa << "\t"
        << _T << "\t"
        << "B" << endl;
    
    removeContinent(continent);
    
    clearPlates();
    fillStructure();
    
    /* correct the ages of sections that were
     affected by the slight offset (openLength)
     of continents*/
    correctAges();
    
    return breakupPos;
}
// --------------------------------------------------
void
Earth::_checkRidgeCreation()
{
    if(_noSubduction || _fixed_configuration)
        return;
    
    bool creation = false;
    
    Subduction* subduction = NULL;
    for(unsigned int i=0; i<_elements.size(); i++)
    {
        if(_elements[i]->isA("RightSubduction"))
        { // orientation of the driving forces for the plate on the right:
            subduction = (Subduction*)_elements[i];
            creation = (subduction->getDirectionDrivingForces(RIGHT) == RIGHT);
            if(creation)
                break;
        } // _elements[i]->isA("RightSubduction")
        
        if(_elements[i]->isA("LeftSubduction"))
        { // orientation of the driving forces for the plate on the left:
            subduction = (Subduction*)_elements[i];
            creation = (subduction->getDirectionDrivingForces(LEFT) == LEFT);
            if(creation)
                break;
        } // _elements[i]->isA("LeftSubduction")
    } // _elements[i]
    
    if(creation)
    {
        ostringstream oss;
        oss << "\t -> " << setw(8) << setprecision(2) << fixed << getTimeMyr()
        << " Myr \t" << setw(8) << _ageMa
        << " Ma \t" << setw(20) << "+creation ridge at: "
        << setw(7) << subduction->getShortId()
        << "\t" << setw(7) << subduction->getPosition() << "\n";
        if(_outputFile.is_open())
            _outputFile << oss.str() << flush;
        if(_showEvents)
            cout << oss.str() << flush;
        
        _ridgeCreation(subduction);
        _checkRidgeCreation();
    }
}
// --------------------------------------------------
void
Earth::_ridgeCreation(Subduction* subduction)
{
    double posStaple = subduction->getPosition(); /* used only if subduction
                                                   not on continent border */
    
    if(subduction->isA("RightSubduction"))
    {
        double openLength = Earth::creationDistance;
        openLength = min( openLength,
                         0.45*subduction->getLeftSection()->getLength() );
        
        double posRidge = ::addPosition(subduction->getPosition(), openLength);
        double posSubduction = ::addPosition(subduction->getPosition(), 2.0 * openLength);
        
        if(subduction->getRightNeighborExtremity())
        {
            LeftContinentExtremity* extremity =
            subduction->getRightNeighborExtremity();
            extremity->setNeighbor(NULL);
            subduction->setRightNeighborExtremity(NULL);
            // put age at zero on the left of the continent
            extremity->clearLeftAges();
            // markers:
			string status = "other";
            changeActiveMargin(subduction, status); /* change this margin from currently active 
													 (with associated subduction to a passive marker 
													 (associated to a continent) */
        }
        else
        { // create a staple
            Staple* newStaple = createStaple(posStaple);
            newStaple->setRightAges(subduction->getRightAges());
            newStaple->setRightAge(subduction->getRightAge());
        }
        
        createRidge(posRidge);
        subduction->setPosition(posSubduction);
        subduction->clearRightAges(); // includes setting _rightAge._value=0.0 and position
        
    }
    else // LEFT
    {
        double openLength = Earth::creationDistance;
        openLength = min( openLength,
                         0.45*subduction->getRightSection()->getLength() );
        
        double posRidge = ::addPosition(subduction->getPosition(), -openLength);
        double posSubduction = ::addPosition(subduction->getPosition(), -2.0 * openLength);
        
        
        if(subduction->getLeftNeighborExtremity())
        {
            RightContinentExtremity* extremity =
            subduction->getLeftNeighborExtremity();
            extremity->setNeighbor(NULL);
            subduction->setLeftNeighborExtremity(NULL);
            // put ages at zero on the right of the continent
            extremity->clearRightAges();
            // markers:
			string status = "other";
            changeActiveMargin(subduction, status); /* change this margin from currently active 
													 (with associated subduction to a passive marker 
													 (associated to a continent) */
        }
        else
        {
            Staple* newStaple = createStaple(posStaple);
            newStaple->setLeftAges(subduction->getLeftAges());
            newStaple->setLeftAge(subduction->getLeftAge());
        }
        
        createRidge(posRidge);
        subduction->setPosition(posSubduction);
        subduction->clearLeftAges();
    }
    
    clearPlates();
    fillStructure();
}
// --------------------------------------------------
void
Earth::_activateAllRidges()
{ /* check if any ridge was set to inactive, set them back
   to active and reset Plates (fillStructure)
   */
    bool change = false;
    for(unsigned int i=0; i<_elements.size(); i++)
    {
        if(_elements[i]->isA("Ridge"))
        {
            Ridge* ridge = (Ridge*)_elements[i];
            if(!ridge->isActive())
            {
                change = true;
                ridge->setActive(true);
            }
        } // isA("Ridge")
    } // elements[i]
    
    if(change)
    {
        clearPlates();
        fillStructure();
    }
    
    _inactiveRidges = false;
}
// --------------------------------------------------
void
Earth::_verifySubductions()
{   /* checks that a subduction does not go in the wrong direction
     (up instead of sinking). Normally useful only in the beginning when the initial
     configuration is not right. */
    bool sinking = true;
    for(unsigned int i=0; i<_elements.size(); i++)
    {
        
        if(_elements[i]->isA("LeftSubduction"))
        {
            LeftSubduction* leftSubduction = (LeftSubduction*)_elements[i];
            double Vsub = leftSubduction->getRightPlate()->getU();
            double Vup  = leftSubduction->getLeftPlate()->getU();
            if(Vup > Vsub)
            {
                // replace the subduction by a staple
                sinking = false;
                Staple* newStaple = createStaple(leftSubduction->getPosition());
                newStaple->setRightAges(leftSubduction->getRightAges());
                newStaple->setRightAge(leftSubduction->getRightAge());
                newStaple->setLeftAges(leftSubduction->getLeftAges());
                newStaple->setLeftAge(leftSubduction->getLeftAge());
                if (leftSubduction->getLeftNeighborExtremity())
                {
                    newStaple->setLeftNeighborExtremity(leftSubduction->getLeftNeighborExtremity());
                    newStaple->getLeftNeighborExtremity()->setNeighbor(newStaple);
                }
                
                
                removeElement(leftSubduction);
                clearPlates();
                fillStructure();
                break;
            }
        } // leftSubduction
        
        if(_elements[i]->isA("RightSubduction"))
        {
            RightSubduction* rightSubduction = (RightSubduction*)_elements[i];
            double Vsub = rightSubduction->getLeftPlate()->getU();
            double Vup  = rightSubduction->getRightPlate()->getU();
            if(Vsub > Vup)
            {
                sinking = false;
                Staple* newStaple = createStaple(rightSubduction->getPosition());
                newStaple->setRightAges(rightSubduction->getRightAges());
                newStaple->setRightAge(rightSubduction->getRightAge());
                newStaple->setLeftAges(rightSubduction->getLeftAges());
                newStaple->setLeftAge(rightSubduction->getLeftAge());
                if(rightSubduction->getRightNeighborExtremity())
                {
                    newStaple->setRightNeighborExtremity(rightSubduction->getRightNeighborExtremity());
                    newStaple->getRightNeighborExtremity()->setNeighbor(newStaple);
                }
                removeElement(rightSubduction);
                clearPlates();
                fillStructure();
                break;
            }
        } // rightSubduction
        
    }
    
    if(!sinking)
        _verifySubductions();
    
}
// --------------------------------------------------
void
Earth::_verifyRidges()
{
    bool tension = true;
    
    for(unsigned int i=0; i<_elements.size(); i++)
    {
        if(_elements[i]->isA("Ridge"))
        {
            Ridge* ridge = (Ridge*)_elements[i];
            if(ridge->isActive())
            {
                double Vright = ridge->getRightPlate()->getU();
                double Vleft = ridge->getLeftPlate()->getU();
                
                if(Vleft < Vright)
                {
                    tension = false;
                    ridge->setActive(false);
                    _inactiveRidges = true;
                    break;
                }
            } // Vleft < Vright
        } // ridge
    }
    
    if(!tension)
    {
        /* ridge that are inactive are not considered
         as interfaces -> will be in the middle of plates */
        clearPlates();
        fillStructure();
        _computeVelocities();
        
        _verifyRidges();
    }
}
// --------------------------------------------------
void
Earth::_checkChangePlateTectonics()
{
    if(_initSubduction.alwaysOne)
        return;
    
    if(_noSubductionPrec && !_noSubduction)
        cout << "+++++ Restart at "
        << _ageMa << " Ma +++++" << endl;
    if(!_noSubductionPrec && _noSubduction)
        cout << "+++++ Plate tectonics stuck at "
        << _ageMa << " Ma +++++" << endl;
}
// --------------------------------------------------
void
Earth::_markersChangeContinents(Continent* oldContinent,
								Continent* newContinent)
{ 
    /* for continental assembly: simply change from old-
     to new-continent */
    for(unsigned int i=0; i<_markers.size(); i++)
    {
		if (!_markers[i]->getContinent())
		{
			cerr << "Warning: collision " << _markers[i]->getShortId()
			<< " does not have continent" << endl;
		}
		else
		{
			if (_markers[i]->getContinent() == oldContinent)
			{
				_markers[i]->setContinent(newContinent);
			}
		}
	}
}
// --------------------------------------------------
void
Earth::_markersChangeContinents(Continent* oldContinent,
								Continent* rightContinent,
								Continent* leftContinent,
								double breakupPos,
								double openLength)
{
    /* for breakup: check where markers are,
     in order to put them in the good continent */
    for(unsigned int i=0; i<_markers.size(); i++)
    {
		if (!_markers[i]->getContinent())
		{
			cerr << "Warning: marker " << _markers[i]->getShortId()
			<< " does not have continent" << endl;
		}
		else
		{
			if (_markers[i]->getContinent() == oldContinent)
			{
				Direction dir;
				if (::isWithinOrEqual(_markers[i]->getPosition(),
									  oldContinent->getPosition(),
									  breakupPos,
									  RIGHT))
				{
					dir = RIGHT;
				}
				else
				{
					dir = LEFT;
				}
				
				double newPos = 0.0;
				if (dir == RIGHT)
				{
					_markers[i]->setContinent(rightContinent);
					newPos = ::addPosition(_markers[i]->getPosition(), -openLength);
				}
				else  // LEFT
				{
					_markers[i]->setContinent(leftContinent);
					newPos = ::addPosition(_markers[i]->getPosition(), openLength);
				}
				_markers[i]->setPosition(newPos);
			}  // if (_markers[i]->getContinent() == oldContinent)
		}
    } // for
}
// ---------------------------------------------
/*
 =======================
 Elements
 =======================
 */

Ridge*
Earth::createRidge(double position)
{
    Ridge* newRidge = new Ridge(this, position);
    // newRidge is	created with _leftAge and _rightAge at 0.0
    
    insertElement(newRidge);
    
    return newRidge;
}
// --------------------------------------------------
Subduction*
Earth::createSubduction(double position, Direction direction)
{
    Subduction* newSubduction;
    
    switch(direction)
    {
        case LEFT :
        {
            newSubduction = new LeftSubduction(this, position);
            break;
        }
            
        case RIGHT :
        {
            newSubduction = new RightSubduction(this, position);
            break;
        }
            
        default :
        {
            cerr << "ERROR: createSubduction used without proper direction" << endl;
            exit(EXIT_FAILURE);
        }
    }
    
    insertElement(newSubduction);
    
    return newSubduction;
}
// --------------------------------------------------
Staple*
Earth::createStaple(double position)
{
    Staple* newStaple = new Staple(this, position);
    insertElement(newStaple);
    
    return newStaple;
}
// --------------------------------------------------
void
Earth::insertElement(GeoElement* element)
{
    unsigned int i = 0;
    for(i=0; i<_elements.size(); i++)
    {
        if(element->getPosition() < _elements[i]->getPosition())
            break;
    }
    
    vector<GeoElement*>::iterator it = _elements.begin() + i;
    _elements.insert(it, element);
}
// --------------------------------------------------
GeoElement*
Earth::findElement(const char* name)
{
    for(unsigned int i=0; i<_elements.size(); i++)
        if(!strcmp(_elements[i]->getId().c_str(), name))
            return _elements[i];
    cerr << "Error: element " << name << " not found" << endl;
    return NULL;
}
// --------------------------------------------------
int
Earth::findElement(GeoElement* element)
{
    if(element)
    {
        for(unsigned int i=0; i<_elements.size(); i++)
            if(_elements[i] == element)
                return i;
        
        cerr << "Error: element " + element->getId() + " not found" << endl;
    }
    
    return -1;
}
// --------------------------------------------------
bool
Earth::removeElement(GeoElement* element)
{
    for(vector<GeoElement*>::iterator it=_elements.begin(); it!=_elements.end(); it++)
        if((*it) == element)
        {
            _elements.erase(it);
            return true;
        }
    
    return false;
}
// --------------------------------------------------
void
Earth::sortElements()
{
    if (_elements.empty())
        return;
    
    unsigned int size = _elements.size();
    vector<GeoElement*> backup = _elements;
    _elements.clear();
    
    while(_elements.size()!=size)
    {
        double position = numeric_limits<double>::infinity(); // "infini"
        GeoElement* element = NULL;
        vector<GeoElement*>::iterator it;
        vector<GeoElement*>::iterator itToDel;
        for(it=backup.begin(); it!=backup.end(); it++)
        {
            if((*it)->getPosition() < position)
            {
                position = (*it)->getPosition();
                element = *it;
                itToDel = it;
            }
        }
        _elements.push_back(element);
        backup.erase(itToDel);
    }
}
// --------------------------------------------------
void
Earth::sortElements(vector<GeoElement*>& elements)
{
    if (elements.empty())
        return;
    
    unsigned int size = elements.size();
    vector<GeoElement*> backup = elements;
    elements.clear();
    
    // re-classes the elements in their order counterclockwise
    while(elements.size()!=size)
    {
        double position = numeric_limits<double>::infinity();
        GeoElement* element = NULL;
        vector<GeoElement*>::iterator it;
        vector<GeoElement*>::iterator itToDel;
        for(it=backup.begin(); it!=backup.end(); it++)
        {
            if((*it)->getPosition() < position)
            {
                position = (*it)->getPosition();
                element = *it;
                itToDel = it;
            }
        }
        elements.push_back(element);
        backup.erase(itToDel);
    }
    
    if(_continents.size() == 0)
        return;
    
    /* for continent extremities: put staples or subductions
     that are at the extremities on the good side of the continent*/
    vector<GeoElement*>::iterator it = elements.begin();
    while(true)
    {
        if( (*it)->isA("RightContinentExtremity") &&
           ((ContinentExtremity*)(*it))->getNeighbor() )
        { // if element on the left shares the same RCE, then switch
            vector<GeoElement*>::iterator itLeft = it + 1;
            if(itLeft == elements.end())
                itLeft = elements.begin();
            
            if( ((ContinentExtremity*)(*it))->getNeighbor() == (*itLeft) )
            { // switch the order of RCE and subduction/staple (subd/stap should be before RCE)
                GeoElement* toKeep = (*it);
                (*it) = (*itLeft);
                (*itLeft) = toKeep;
            }
        } // *it is a RCE
        
        if( (*it)->isA("LeftContinentExtremity") &&
           ((ContinentExtremity*)(*it))->getNeighbor() )
        { // if element on the right share the same LCE: switch
            vector<GeoElement*>::iterator itRight;
            if(it == elements.begin())
                itRight = elements.end() - 1;
            else
                itRight = it - 1;
            
            if( ((ContinentExtremity*)(*it))->getNeighbor() == (*itRight) )
            { // switch the order of LCE and subduction/staple (subd/stap should be after LCE)
                GeoElement* toKeep = (*it);
                (*it) = (*itRight);
                (*itRight) = toKeep;
            }
        } // it is LCE
        
        it++;
        if(it == elements.end())
            break;
    } // while true
    
}
// --------------------------------------------------
vector<GeoElement*>
Earth::findInterfaces()
{
    // --- Interfaces
    vector<GeoElement*> elements;
    for(unsigned int i=0; i<_elements.size(); i++)
    {
        if(!_elements[i]->isA("Staple"))
            elements.push_back(_elements[i]);
    }
    sortElements(elements);
    
    return elements;
}

/*
 =======================
 Plates
 =======================
 */

Plate*
Earth::createPlate(vector<PlateSection*>& sections)
{
    Plate* newPlate = new Plate(this, sections);
    _plates.push_back(newPlate);
    
    if(newPlate->getLeftElement() && newPlate->getLeftElement()->isA("Interface"))
        ((Interface*)newPlate->getLeftElement())->setRightPlate(newPlate);
    if(newPlate->getRightElement() && newPlate->getRightElement()->isA("Interface"))
        ((Interface*)newPlate->getRightElement())->setLeftPlate(newPlate);
    
    return newPlate;
}
// --------------------------------------------------
Plate*
Earth::findPlate(const char* plateName)
{
    for(unsigned int i=0; i<_plates.size(); i++)
        if(!strcmp(_plates[i]->getId().c_str(), plateName))
            return _plates[i];
    cerr << "Error: plate " << plateName << " not found" << endl;
    return NULL;
}
// --------------------------------------------------
int
Earth::findPlate(Plate* plate)
{
    if(plate)
    {
        for(unsigned int i=0; i<_plates.size(); i++)
            if(_plates[i] == plate)
                return i;
        
        cerr << "Error: plate " + plate->getId() + " not found" << endl;
    }
    
    return -1;
}
// --------------------------------------------------
bool
Earth::removePlate(Plate* plate)
{
    for(vector<Plate*>::iterator it=_plates.begin(); it!=_plates.end(); it++)
        if((*it) == plate)
        {
            _plates.erase(it);
            delete plate;
            return true;
        }
    return false;
}
// -----------------------------------------------------
bool
Earth::removeMarker(Marker* marker)
{
    for (vector<Marker*>::iterator it=_markers.begin(); it!=_markers.end(); it++)
    {
        if ((*it) == marker)
        {
            _markers.erase(it);
            delete marker;
            return true;
        }
    }
    return false;
}
// ------------------------------------------------------
bool
Earth::changeActiveMargin(Subduction* subduction, string status)
{
    if (subduction->getActiveMargin())
    {
		ActiveMargin* margin = subduction->getActiveMargin();
		margin->setStatus(status);
		margin->setSubduction(NULL);
		
        subduction->setActiveMargin(NULL);
        return true;
    }
    return false;
}
// ------------------------------------------------------
bool
Earth::changeActiveMargin(ContinentExtremity* extremity, string status)
{
    if (extremity->getActiveMargin())
    {
        ActiveMargin* margin = extremity->getActiveMargin();
        if (margin->getStatus().compare("other") == 0)
        {
            margin->setStatus(status);
            margin->setContinentExtremity(NULL);
            
            extremity->setActiveMargin(NULL);
            return true;
        }
    }
    return false;
}

/*========================================
 various computations
 =======================================*/
void
Earth::computeThicknessesMantle()
{ /* Thicknesses that are used to compute the horizontal drag on the mantle.
   Lower mantle: half the layer only.
   */
    
    // subcontinental layer
    _thicknesses.subcont = _thick_subcont;
    
    // asthenosphere
    _thicknesses.oceanAst = _thick_ast;
    _thicknesses.contAst =      // = 0 if subcontinental layer thicker than asthenosphere
    max(0.0, _thick_ast - _thick_subcont);
    
    // upper mantle
    _thicknesses.oceanUM = _D - _thick_ast;
    _thicknesses.contUM = _D - max(_thick_ast, _thick_subcont);
    
    // lower mantle
    _thicknesses.LM = 0.5 * _d - _D;
    
    if(_thicknesses.contUM < 0.0 || _thicknesses.oceanUM < 0.0)
    {
        cerr << "ERROR: upper mantle must be thicker than asthenosphere and/or subcontinental layer" << endl;
        exit(EXIT_FAILURE);
    }
}

// --------------------------------------------------
double 
Earth::computeHeatFlux(double age)
{ // input age in Myr, must already take into account SSC, Qmax etc.
    double Q = 0.0;
    
    if(age >= 0.0) // ocean
    {
        Q = _preFactor_heatFlow * (_T - _T_surf) /
        sqrt( max(3.0E10, ::Myr_to_sec(age)) ); // limit min (~1000yr) to avoid NaN
    }
    else // continent (imposed age=-1 when using computeHeatFlux())
    {
        if(_thick_continent > 0.0)
            Q = _k_continent * (_T - _T_surf) / _thick_continent;
    }
    return Q;
}
// --------------------------------------------------
double
Earth::linearAge(double position, Age age1, Age age2)
{ /* finds the age for position comprised between
   age1 and age2 (linear interpolation) */
    return ::lin_Pos_to_Y( position,
                          age1.getPosition(),
                          age2.getPosition(),
                          age1.getValue(),
                          age2.getValue() );
}
// --------------------------------------------------
double
Earth::linearPosition(double age, Age age1, Age age2)
{ /* find the value of the age so that it's at the good position
   for the linear interpolation between age1 and age2*/
    return ::lin_Y_to_Pos( age,
                          age1.getPosition(),
                          age2.getPosition(),
                          age1.getValue(),
                          age2.getValue() );
}
// --------------------------------------------------
double
Earth::computeBathyHSCM(double age)
{ /* input age in Myr, must already take into account SSC, Qmax etc.
   Continents are recognized because of negative age */
    
    double bathy = 1000.0; // continental elevation
    if(age >= 0.0) // ocean
    {
        // b: depth below ridge
        double b = _preFactor_bathyHSCM * _rho_um /
        ( _rho_um - _rho_seawater ) * ( _T - _T_surf) *
        sqrt( ::Myr_to_sec(age) );
        
        bathy = _sealevel_HSCM - _ridgeDepth_HSCM - b;
    }
    
    return bathy;
}
// --------------------------------------------------
pair<double, double>
Earth::computeBathyPM(double ageMyr)
{
    double b95 = computeDepthBelowRidgePM(ageMyr, _zm95);
    double b125 = computeDepthBelowRidgePM(ageMyr, _zm125);
    
    double bathy95 = 1000.0;
    double bathy125 = 1000.0;
    if(ageMyr >= 0.0)
    {
        bathy95 = _sealevel_PM95 - _ridgeDepth_PM95 - b95;
        bathy125 = _sealevel_PM125 - _ridgeDepth_PM125 - b125;
    }
    
    return make_pair(bathy95, bathy125);
}
// --------------------------------------------------
double
Earth::computeDepthBelowRidgePM(double ageMyr, double zm)
{ /* input age in Myr; computes bathy for plate model 
   with zm=95km and zm=125km */
    double b = 0.0;
    
    if(ageMyr >= 0.0) // ocean
    {
        double c1 = _alpha_pl * (_T - _T_surf) *
        _rho_um / (_rho_um - _rho_seawater) * zm;
        
        double c2 = _preFactor_bathyPM2 / pow(zm, 2);
        
        // adapt nmax to age:
        unsigned int nmax = 200;
        if(ageMyr < 0.1)
            nmax = 99;
        else if(ageMyr < 1.0)
            nmax = 19;
        else if(ageMyr < 10.0)
            nmax = 9;
        else
            nmax = 3;
        
        //b: depth below ridge
        b = c1 * ( 0.5 - 4.0 / pow(M_PI, 2) *
                  _sumForBathyPM(c2, ageMyr, nmax) );
    }
    
    return b;
}
// --------------------------------------------------
double
Earth::_sumForBathyPM(double coeff, double ageMyr, unsigned int nmax)
{
    double som = 0.0;
    for(unsigned int n=0; n<=nmax; n++)
    {
        unsigned int m = pow(1.0 + 2.0*n, 2.0);
        som += exp(-coeff * m * ageMyr) / m;
    }
    
    return som;
}
// --------------------------------------------------
double
Earth::computeSlabPull(double thickness, double depth)
{ 
    double slabPull = 0.0;
    
    double deltaRho = _DeltaRho_p;
    if(_slabPull_rho_depends_on_T)
        deltaRho = _rho_pl - _rho_um;
    
    slabPull = _coeffSlabPull * deltaRho * _g * thickness * depth;
    
    return slabPull;
}
// --------------------------------------------------
void
Earth::completeForces()
{ /* CG: complete the computation of Mantle Drag, 
   Bending and Viscous Shear (only equivalent viscosities
   were computed so far here) */
    for(unsigned int i=0; i<_plates.size(); i++)
        _plates[i]->completeForces();
}
// --------------------------------------------------
void
Earth::_computeVelocities()
{
    if(_noSubduction)
    {
        for(unsigned int i=0; i<_plates.size(); i++)
            _plates[i]->setU(0.0);
        return;
    }
    
    _simpleVelocities();
    
}
// --------------------------------------------------
void 
Earth::_simpleVelocities()
{ /* simple force balance : 
   RP + SP + SS  = etah*U + etal*U if subduction,
   for subduction upper or driving plate.
   */
    
    for(unsigned int i=0; i<_plates.size(); i++)
    {
        _plates[i]->updateForces();
        Force* forces = _plates[i]->accessForces();
        
        double U = (forces->getDrivingForces() / Earth::force0) /
        (forces->getEtaH() + forces->getEtaLeft() + forces->getEtaRight());
        _plates[i]->setU(U);
    }
    
}
// --------------------------------------------------
void
Earth::_solveVelocities()
{
    /* solve "almost" tridiagonal system to get velocity
     of each plate. "almost" means there are non zero
     values in the corner of the matrix (periodic system) */
    
    // CREATE THE MATRIX ---------------------------------
    vector<Tridiag> matrix;/* actually contains upper and lower-diagonals,
                            diagonal, and RHS vector */
    /*
     [c0  l0  0   0   0  ... 0  r0]    [rhs0]
     [r1  c1  l1  0   0  ...     0]    [rhs1]
     [0   r2  c2  l2  0  ...     0]    [rhs2]
     [ ...                        ]    [... ]
     [ln   0 ...            rn  cn]    [rhsn]
     */
    
    for(unsigned int i=0; i<_plates.size(); i++)
    {
        _plates[i]->updateForces();
        Force* forces = _plates[i]->accessForces();
        
        Tridiag line;
        line.rhs = 0.0;
        line.center = 0.0;
        line.left = 0.0;
        line.right = 0.0;
        
        // rhs of the matrix (driving forces)
        line.rhs = forces->getDrivingForces() / Earth::force0;
        
        // coeff for _plates[i]
        line.center = forces->getEtaH()
        + forces->getEtaRight() + forces->getEtaLeft();
        /* upper diagonal: coeff related to _plates[i+1] (left plate) */
        if(_plates[i]->getLeftElement()->isA("LeftSubduction"))
            line.left  = - forces->getEtaLeft();
        /* lower diagonal: coeff related to _plates[i-1] (right plate) */
        if(_plates[i]->getRightElement()->isA("RightSubduction"))
            line.right = - forces->getEtaRight();
        
        matrix.push_back(line);
    } // _plates[i]
    
    if(_plates.size() == 2)
    { // the left and right plates are then the same
        matrix[0].left += matrix[0].right;
        matrix[1].right += matrix[1].left;
    }
    
    // solve M U = RHS:
    vector<double> U = solverTridiagMatrix(matrix);
    
    for(unsigned int i=0; i<_plates.size(); i++)
        _plates[i]->setU(U[i]);
    
    // TEST
    _verifyVelocities(); // To comment at a point...
    
}

// --------------------------------------------------
vector<double>
Earth::solverTridiagMatrix(vector<Tridiag> matrix)
{ // CG
    /* solve the tridiagonal system for plates' velocities */
    vector<double> U;
    unsigned int N = matrix.size();
    
    // to simplify writing hereafter:
    vector<double> L;
    vector<double> C;
    vector<double> R;
    vector<double> rhs;
    for(unsigned int i=0; i<N; i++)
    {
        L.push_back(matrix[i].left);
        C.push_back(matrix[i].center);
        R.push_back(matrix[i].right);
        rhs.push_back(matrix[i].rhs);
        U.push_back(0.0); // all velocities set at zero initially
    }
    
    if(N < 4) // 2x2 or 3x3
    { // direct inversion of the matrix
        double m[N][N]; // inverse matrix
        if(N == 2)
        { // 2x2
            double det = C[0] * C[1] - L[0] * R[1];
            // first row:
            m[0][0] = C[1] / det;
            m[0][1] = -L[0] / det;
            // second row:
            m[1][0] = -R[1] / det;
            m[1][1] = C[0] / det;
        }
        else
        { // 3x3
            double det =
            C[0] * (C[1]*C[2] - R[2]*L[1]) -
            R[1] * (L[0]*C[2] - R[2]*R[0]) +
            L[2] * (L[0]*L[1] - C[1]*R[0]) ;
            // first row:
            m[0][0] = ( C[1]*C[2] - L[1]*R[2]) / det;
            m[0][1] = (-L[0]*C[2] + R[0]*R[2]) / det;
            m[0][2] = ( L[0]*L[1] - R[0]*C[1]) / det;
            // second row:
            m[1][0] = (-R[1]*C[2] + L[1]*L[2]) / det;
            m[1][1] = ( C[0]*C[2] - R[0]*L[2]) / det;
            m[1][2] = (-C[0]*L[1] + R[0]*R[1]) / det;
            // third row
            m[2][0] = ( R[1]*R[2] - C[1]*L[2]) / det;
            m[2][1] = (-C[0]*R[2] + L[0]*L[2]) / det;
            m[2][2] = ( C[0]*C[1] - L[0]*R[1]) / det;
        }
        
        for(unsigned int i=0; i<N; i++)
        {
            U[i] = 0.0;
            for(unsigned int j=0; j<N; j++)
                U[i] += m[i][j] * rhs[j];
        }
    } // 2x2 or 3x3
    else // N > 3
    {
        // SOLVER for a matrix at least 4x4 ------------------
        
        double ell[N]; // arrays go from 0 to N-1
        double u[N];
        double lambda[N];
        double nu[N];
        double y[N];
        double sol[N];
        double som = 0.0;
        
        /* - - - - - - - - - -
         LU Decomposition:
         - - - - - - - - - -  */
        
        // diagonal of L and upper diag of U - - - -
        ell[0] = C[0];
        u[0] = L[0] / ell[0];
        for(unsigned int i = 1; i<N-2; i++)
        {
            ell[i] = C[i] - R[i]*u[i-1];
            u[i] = L[i] / ell[i];
        }
        ell[N-2] = C[N-2] -
        R[N-2] * u[N-3];
        
        // last colum of U - - - - - - - - - -
        nu[0] = R[0] / ell[0];
        for(unsigned int i = 1; i<N-2; i++)
        {
            nu[i] = - R[i]*nu[i-1] / ell[i];
        }
        nu[N-2] = (L[N-2] - R[N-2]*nu[N-3]) /
        ell[N-2];
        
        // last row of L - - - - - - - - - - -
        lambda[0] = L[N-1];
        for(unsigned int i = 1; i<N-2; i++)
        {
            lambda[i] = - lambda[i-1] * u[i-1];
        }
        lambda[N-2] = R[N-1] -
        lambda[N-3]*u[N-3];
        
        for(unsigned int i=0; i<N-1; i++)
            som += lambda[i]*nu[i];
        lambda[N-1] = C[N-1] - som;
        
        
        // solving L y = rhs - - - - - - - - -
        y[0] = rhs[0] / ell[0];
        som = lambda[0] * y[0];
        for(unsigned int i = 1; i<N-1; i++)
        {
            y[i] = (rhs[i] - R[i]*y[i-1]) / ell[i];
            som += lambda[i]*y[i];
        }
        y[N-1] = (rhs[N-1] - som) /
        lambda[N-1];
        
        // solving U sol = y - - - - - - - - -
        sol[N-1] = y[N-1];
        sol[N-2] = y[N-2] -
        nu[N-2] * sol[N-1];
        for(int i = N-3; i>-1; i--)
        {
            sol[i] = y[i] - nu[i]*sol[N-1] -
            u[i]*sol[i+1];
        }
        // the solution for each plate is in sol[i]
        
        for(unsigned int i=0; i<N; i++)
            U[i] = sol[i];
    } // N > 3
    
    return U;
}
// --------------------------------------------------
bool
Earth::_verifyVelocities()
{ // CG
    bool correct = true;
    
    if(_timestep < 1)
        return correct;
    
    for(unsigned int i=0; i<_plates.size(); i++)
    {
        Plate* plate = _plates[i];
        
        double etaH = plate->accessForces()->getEtaH();
        double etaL = plate->accessForces()->getEtaLeft();
        double etaR = plate->accessForces()->getEtaRight();
        
        double driving = plate->accessForces()->getDrivingForces() / Earth::force0; // dimensionless
        
        double U = plate->getU();
        double ULeft = 0.0;
        double URight = 0.0;
        if(plate->getLeftElement()->isA("LeftSubduction"))
            ULeft  = plate->getLeftPlate()->getU();
        if(plate->getRightElement()->isA("RightSubduction"))
            URight = plate->getRightPlate()->getU();
        
        correct = ::areEqual(driving,
                             etaH*U + etaL*(U-ULeft) + etaR*(U-URight),
                             0.001);
        if(!correct)
        {
            cerr << setprecision(3) << fixed << _ageMa
            << "\t ERROR: wrong velocity for plate " << _plates[i]->getId()
            << "\t at right pos=" << _plates[i]->getRightElement()->getPosition()
            << endl;
            cerr << "\t Driving = " << driving
            << "\t Resistive = " << (etaH*U + etaL*(U-ULeft) + etaR*(U-URight))
            << endl;
            cerr << "\t right: " << _plates[i]->getRightElement()->getId()
            << "\t" << _plates[i]->getRightElement()->getPosition() << endl;
            cerr << "\t left: " << _plates[i]->getLeftElement()->getId()
            << "\t" << _plates[i]->getLeftElement()->getPosition() << endl;
            cerr << "\t eta:" << etaH << "\t" << etaR << "\t" << etaL
            << "\t U:" << U << "\t" << URight << "\t" << ULeft << endl;
        }
    }
    
    return correct;
}



/*
 =======================
 Cells
 =======================
 */

vector<Cell*>
Earth::getCells()
{
    vector<Cell*> cells;
    
    for(unsigned int i=0; i<_plates.size(); i++)
    {
        vector<Cell*> cell = _plates[i]->accessCells();
        for(unsigned int j=0; j<cell.size(); j++)
            cells.push_back(cell[j]);
    }
    
    return cells;
}
// --------------------------------------------------
Cell*
Earth::findCell(const char* name)
{
    Cell* cell = NULL;
    
    vector<Cell*> cells = getCells();
    
    for(unsigned int i=0; i<cells.size(); i++)
    {
        if(!strcmp(cells[i]->getId().c_str(), name))
        {
            cell = cells[i];
            break;
        }
    }
    
    if(!cell)
        cerr << "Earth::findCell error : cell " << name << " not found" << endl;
    
    return cell;
}
// --------------------------------------------------
int
Earth::findCell(Cell* cell)
{
    if(cell)
    {
        vector<Cell*> cells = getCells();
        for(unsigned int i=0; i<cells.size(); i++)
            if(cells[i] == cell)
                return i;
        
        cerr << "Error: cell " + cell->getId() + " not found" << endl;
    }
    
    return -1;
}

/*
 ======================
 Continents
 ======================
 */

Continent*
Earth::createContinent(double position, double length, string origin)
{
    Continent* newContinent = new Continent(this, position, length, origin);
    _continents.push_back(newContinent);
    
    //resetContinents();
    
    return newContinent;
}
// --------------------------------------------------
Continent*
Earth::findContinent(const char* name)
{
    for(unsigned int i=0; i<_continents.size(); i++)
        if(!strcmp(_continents[i]->getId().c_str(), name))
            return _continents[i];
    cerr << "Error: continent " << name << " not found" << endl;
    return NULL;
}
// --------------------------------------------------
int
Earth::findContinent(Continent* continent)
{
    if(continent)
    {
        for(unsigned int i=0; i<_continents.size(); i++)
            if(_continents[i] == continent)
                return i;
        cerr << "Error: continent " + continent->getId() + " not found" << endl;
    }
    return -1;
}
// --------------------------------------------------
bool
Earth::removeContinent(Continent* continent)
{
    for(vector<Continent*>::iterator it=_continents.begin(); it!=_continents.end(); it++)
        if((*it) == continent)
        {
            _continents.erase(it);
            delete continent;
            return true;
        }
    return false;
}
// --------------------------------------------------
void
Earth::resetContinents()
{ // find which section is linked to each continent
    for(unsigned int i=0; i<_continents.size(); i++)
    {
        bool found = false;
        double middle = _continents[i]->getMiddlePosition();
        
        for(unsigned int j=0; j<_plates.size(); j++)
        {
            vector<PlateSection*> sections = _plates[j]->accessSections();
            for(unsigned int k=0; k<sections.size(); k++)
            {
                if(sections[k]->contains(middle))
                {
                    _continents[i]->setPlateSection(sections[k]);
                    sections[k]->setContinent(_continents[i]);
                    found = true;
                    break;
                }
            } // sections[k]
            if(found)
                break;
        } // plates[j]
    } // continents[i]
}
// -------------------------------------------------
RightContinentExtremity*
Earth::createRightContinentExtremity(Continent* continent)
{
    RightContinentExtremity* rightExtremity = new RightContinentExtremity(this, continent);
    
    if (continent->getRightExtremity())
        cout << "WARNING: rightExtremity created for a continent that already has one." << endl;
    else
        continent->setRightExtremity(rightExtremity);
    
    return rightExtremity;
}
// -------------------------------------------------
LeftContinentExtremity*
Earth::createLeftContinentExtremity(Continent* continent)
{
    LeftContinentExtremity* leftExtremity =
    new LeftContinentExtremity(this, continent);
    
    if(continent->getLeftExtremity())
        cout << "WARNING: leftExtremity created for a continent that already has one." << endl;
    else
        continent->setLeftExtremity(leftExtremity);
    
    return leftExtremity;
}

/*
 =======================
 Properties
 =======================
 */
void
Earth::initSlabs()
{ /* CG: compute thicknesses of slabs in the upper mantle initially,
   considering the seafloor age at the trench was the initial given one
   before starting the run */
    for(unsigned int i=0; i<_elements.size(); i++)
    {
        if(_elements[i]->isA("Subduction"))
        {
            Subduction* subduction = (Subduction*)_elements[i];
            double depth = subduction->getDepth();
            
            double thickness = 0.0;
            if(subduction->isA("LeftSubduction"))
                thickness = subduction->computePlateThickness(RIGHT);
            else // right subduction
                thickness = subduction->computePlateThickness(LEFT);
            
            // put a minimum depth:
            depth = max(depth, thickness * _depthFactor);
            subduction->setDepth(depth);
            
            subduction->initThickness(depth, thickness);
        } // subduction
    }
}
// --------------------------------------------------
double
Earth::offsetPosition(GeoElement* element, 
                      double length, Direction direction)
{
    /* computes a position LEFT or RIGHT of the position of
     element, located at "length" from the element but
     checking that the platesection is wide enough. */
    
    double newPosition = 0.0;
    
    switch(direction)
    {
        case RIGHT :
        {
            double offset =
            min ( length,
                 0.45*element->getRightSection()->getLength() );
            newPosition = ::addPosition(element->getPosition(), -offset);
            
            break;
        }
            
        case LEFT :
        {
            double offset =
            min ( length,
                 0.45*element->getLeftSection()->getLength() );
            newPosition = ::addPosition(element->getPosition(), offset);
            break;
        }
            
        default :
        {
            cerr << "ERROR: Earth::offsetPosition() is used with no direction"
            << endl;
            exit(EXIT_FAILURE);
        }
            
    }
    
    return newPosition;
}
// --------------------------------------------------


/*
 =======================
 Global parameters
 =======================
 */

void
Earth::setTimeParameter(string name, double value)
{
    _timeParameters[name] = value;
}
// --------------------------------------------------
void
Earth::setFixedParameter(string name, double value)
{
    _fixedParameters[name] = value;
}
// --------------------------------------------------
void
Earth::setPhysicalParameter(string name, double value)
{
    _physicalParameters[name] = value;
}
// --------------------------------------------------
void
Earth::setModelParameter(string name, double value)
{
    _modelParameters[name] = value;
}

// --------------------------------------------------
double
Earth::getTimeParameter(string name)
{
    double value = 0.0;
    map<string, double>::iterator it = _timeParameters.find(name);
    if(it!=_timeParameters.end())
        value = it->second;
    else
        cerr << "Parameter \"" + name + "\" not found for Earth --> dummy value" << endl;
    
    return value;
}
// --------------------------------------------------
double
Earth::getFixedParameter(string name)
{
    double value = 0.0;
    map<string, double>::iterator it = _fixedParameters.find(name);
    if(it!=_fixedParameters.end())
        value = it->second;
    else
        cerr << "Parameter \"" + name + "\" not found for Earth --> dummy value" << endl;
    
    return value;
}
// --------------------------------------------------
double
Earth::getPhysicalParameter(string name)
{
    double value = 0.0;
    map<string, double>::iterator it = _physicalParameters.find(name);
    if(it!=_physicalParameters.end())
        value = it->second;
    else
        cerr << "Parameter \"" + name + "\" not found for Earth --> dummy value" << endl;
    
    return value;
}
// --------------------------------------------------
double
Earth::getModelParameter(string name)
{
    double value = 0.0;
    map<string, double>::iterator it = _modelParameters.find(name);
    if(it!=_modelParameters.end())
        value = it->second;
    else
        cerr << "Parameter \"" + name + "\" not found for Earth --> dummy value" << endl;
    
    return value;
}
// --------------------------------------------------
void
Earth::changeTimeParameters()
{ /* change value of startAge
   _timeParameters[] (for Macma::writeXmlFile()) */
    setTimeParameter("startAge", _ageMa);
}
// --------------------------------------------------
void
Earth::changePhysicalParameters()
{ /* change values of T_m_init in the map
   _physicalParameters[] (for Macma::writeXmlFile()) */
    setPhysicalParameter("T_m_init", _T);
}
// --------------------------------------------------
void
Earth::setSubductionMode()
{
    bool set = false;
    if(getModelParameter("initMode_brittle") != 0.0)
    {
        set = true;
        _initSubduction.mode = BRITTLE;
        setModelParameter("initMode_constant", 0.0);
        setModelParameter("initMode_convective", 0.0);
    }
    else if(!set && getModelParameter("initMode_constant") != 0.0)
    {
        set = true;
        _initSubduction.mode = CONSTANT;
        setModelParameter("initMode_brittle", 0.0);
        setModelParameter("initMode_convective", 0.0);
    }
    else if(!set && getModelParameter("initMode_convective") != 0.0)
    {
        set = true;
        _initSubduction.mode = CONVECTIVE;
        setModelParameter("initMode_brittle", 0.0);
        setModelParameter("initMode_constant", 0.0);
    }
    else
    {
        cerr << "ERROR: subduction initiation mode is badly set" << endl;
        exit(EXIT_FAILURE);
    }
}
// --------------------------------------------------
void
Earth::setSubductionPlace()
{
    _initSubduction.continents = false;
    _initSubduction.staples = false;
    _initSubduction.upperPlates = false;
    
    _initSubduction.ratioAgeValue = numeric_limits<double>::infinity();
    _initSubduction.reverse = false;
    
    if(getModelParameter("initPlace_continents") != 0.0)
        _initSubduction.continents = true;
    if(getModelParameter("initPlace_staples") != 0.0)
        _initSubduction.staples = true;
    if(getModelParameter("initPlace_upperPlates") != 0.0)
    {
        _initSubduction.upperPlates = true;
        
        _initSubduction.upperPlateMinAge =
        getModelParameter("upperPlates_minimumAge");
        
        if(getModelParameter("upperPlates_ageRatioCriterion") != 0.0)
        {
            _initSubduction.ageRatioCriterion = true;
            _initSubduction.ratioAgeValue =
            getModelParameter("upperPlates_ratioAgeValue");
        }
        else
            _initSubduction.ageRatioCriterion = false;
        
        if(getModelParameter("reverseSubduction") != 0.0)
            _initSubduction.reverse = true;
    }
    
    _initSubduction.tauSub_randomNoise = getModelParameter("tauSub_randomNoise");
    _initSubduction.randomTauSub =
    ( _initSubduction.tauSub_randomNoise > 0.0 ? true : false );
    _initSubduction.alwaysOne =
    (getModelParameter("alwaysSubduct") != 0.0 ? true : false );
}
// --------------------------------------------------
void
Earth::setLimitedThickening()
{
    _limitedThickening = false;
    if(getModelParameter("small_scale_convection") != 0.0)
        _limitedThickening = true;
}
// --------------------------------------------------
void
Earth::setSSCMode()
{
    if (getModelParameter("sscMode_constant") != 0.0)
    {
        _SSCMode = SSCCONSTANT;
        setModelParameter("sscMode_convective", 0.0);
    }
    else
    {
        _SSCMode = SSCCONVECTIVE;
        setModelParameter("sscMode_constant", 0.0);
    }
}
// --------------------------------------------------
void
Earth::setFixedConfiguration()
{
    _fixed_configuration = false;
    if(getModelParameter("fixed_configuration") != 0.0)
        _fixed_configuration = true;
}
// --------------------------------------------------
void
Earth::setDepletionMode()
{
    bool set = false;
    if(getModelParameter("depletion_always") != 0.0)
    {
        set = true;
        _depletionMode = DEPLETED;
    }
    else if(!set && getModelParameter("depletion_never") != 0.0)
    {
        set = true;
        _depletionMode = PRIMITIVE;
    }
    else if(!set && getModelParameter("depletion_adaptative") != 0.0)
    {
        set = true;
        _depletionMode = ADAPTATIVE;
    }
    else
    { // either several or no mode are chosen...
        cerr << "ERROR: depletion mode is badly set" << endl;
        exit(EXIT_FAILURE);
    }
}
// --------------------------------------------------

/*
 =====================================
 Ages
 =====================================
 */
void
Earth::initAges()
{
    for(unsigned int i=0; i<_plates.size(); i++)
    {
        vector<PlateSection*>& sections = _plates[i]->accessSections();
        for(unsigned int j=0; j<sections.size(); j++)
        {
            sections[j]->initAges();
        }
    }
    
    /* make sure continents' extremities have the same
     ages as neighbors */
    for(unsigned int i=0; i<_continents.size(); i++)
        _continents[i]->shareNeighborAges();
    
}

// --------------------------------------------------
vector<double>
Earth::getOceanAges()
{
    vector<double> ages;
    for(unsigned int i=0; i<_plates.size(); i++)
    {
        vector<PlateSection*> sections = _plates[i]->accessSections();
        for(unsigned int j=0; j<sections.size(); j++)
        {
            if(!sections[j]->isContinental())
            {
                vector<Age> currentAges = sections[j]->getRightElement()->getLeftAges();
                for(unsigned int k=0; k<currentAges.size(); k++)
                    ages.push_back(currentAges[k].getValue());
            }
        }
    }
    return ages;
}

// --------------------------------------------------
void
Earth::getAllAges()
{
    _ages.clear();
    
    for(unsigned int i=0; i<_plates.size(); i++)
    {
        vector<PlateSection*> sections = _plates[i]->accessSections();
        for(unsigned int j=0; j<sections.size(); j++)
        {
            sections[j]->update_U_for_Ages(); // in order to set veloc for _ages[]
            vector<Age> ages = sections[j]->getRightElement()->getLeftAges();
            
            for(unsigned int k=0; k<ages.size(); k++)
            {
                _ages.push_back(ages[k]);
            }
        } // sections[j]
    } // _plates[i]
    
    if(_ages.size() > 0)
    {
        /* Put in correct order of positions.
         Works only if ages of sections do not overlap or
         are mixed */
        while(_ages[0].getPosition() >= _ages[_ages.size()-1].getPosition())
        {
            _ages.insert(_ages.begin(), _ages[_ages.size()-1]);
            _ages.pop_back();
        }
        
        // for drawing: add points for positions 0 and 360.0
        if(! ::areEqual(_ages[0].getPosition(), 0.0))
        {
            Age borderAge = _ages[0];
            borderAge.setPosition(0.0);
            _ages.insert(_ages.begin(), borderAge);
        }
        if(! ::areEqual(_ages[_ages.size()-1].getPosition(), 360.0))
        {
            Age borderAge = _ages[_ages.size()-1];
            borderAge.setPosition(360.0);
            _ages.push_back(borderAge);
        }
    }
}
// ----------------------------------------
void
Earth::effectiveAges(vector<Age>& ages)
{
    for(unsigned int i=0; i<ages.size(); i++)
        ages[i].setEffValue(ages[i].getValue());
    
    /* take SSC into account */
    if(_limitedThickening)
        maxAgesLimit(ages);
    
    /* limit in heat flux close to the ridge */
    minAgesLimit(ages);
}
// --------------------------------------------------
void
Earth::maxAgesLimit(vector<Age>& ages)
{ // modifies ages to take into account SSC
    
    if (ages.size() > 1)
    {
        vector<Age>::iterator it;
        vector<Age>::iterator itp;
        
        if(ages[0].getValue() < ages[ages.size()-1].getValue())
        { // ages in increasing order from right to left
            itp = ages.end()-1;
            it = itp - 1;
            while(true)
            {
                if((*it).getValue() > _tau_ssc)
                { // both it and itp are above max age
                    (*itp).setEffValue(_tau_ssc);
                }
                else if((*itp).getValue() > _tau_ssc)
                { // the limit is between it and itp
                    (*itp).setEffValue(_tau_ssc);
                    if(_insertAges)
                    {
                        double position =
                        linearPosition( _tau_ssc, (*it), (*itp) );
                        Age newAge = Age(this, _tau_ssc, position);
                        newAge.setExtra(true);
                        ages.insert(itp, newAge);
                    }
                    break;
                }
                else
                { // both it and itp are below max age
                    break;
                }
                
                itp --;
                it --;
                if(itp == ages.begin())
                {
                    if((*itp).getValue() > _tau_ssc)
                        (*itp).setEffValue(_tau_ssc);
                    break;
                }
            } // while(true)
        }
        else
        { // ages in decreasing order from right to left
            it = ages.begin();
            itp = it + 1;
            while(true)
            {
                if((*itp).getValue() > _tau_ssc)
                { // both itp and it are above max age
                    (*it).setEffValue(_tau_ssc);
                }
                else if((*it).getValue() > _tau_ssc)
                { // limit is between it and itp
                    (*it).setEffValue(_tau_ssc);
                    if(_insertAges)
                    {
                        double position =
                        linearPosition( _tau_ssc,  (*it), (*itp) );
                        Age newAge = Age(this, _tau_ssc, position);
                        newAge.setExtra(true);
                        ages.insert( itp, newAge );
                    }
                    break;
                }
                else
                { // both it and itp are below max age
                    break;
                }
                
                it ++;
                itp ++;
                if(itp == ages.end())
                {
                    if((*it).getValue() > _tau_ssc)
                        (*it).setEffValue(_tau_ssc);
                    break;
                }
                
            } // while(true)
        }
    } // ages.size() > 1
    else
    {
        ages[0].setEffValue(_tau_ssc);
    }
}
// ----------------------------------------
void
Earth::minAgesLimit(vector<Age>& ages)
{ //modifies _ages to take into account Qmax and min_plate_thick
    double minAge = - numeric_limits<double>::infinity();
    
    if(_Qmax > 0.0)
    {
        minAge =
        ::sec_to_Myr( pow( _k_ocean * (_T - _T_surf)/_Qmax, 2 )
                     / (M_PI * _kappa) );
    }
    
    if(_age_min_thickness > 0.0)
        minAge = max(minAge, _age_min_thickness);
    
    if(minAge < 0.0)
        return;
    
    if (ages.size() > 1)
    {
        vector<Age>::iterator it;
        vector<Age>::iterator itp;
        if(ages[0].getValue() < ages[ages.size()-1].getValue())
        { // ages are in increasing order from right to left
            it = ages.begin();
            itp = it + 1;
            while(true)
            {
                if((*itp).getValue() < minAge)
                { // both it and itp are below minAge
                    (*it).setEffValue(minAge);
                }
                else if((*it).getValue() < minAge)
                { // limit is between it and itp
                    (*it).setEffValue(minAge);
                    if(_insertAges)
                    {
                        double position =
                        linearPosition( minAge,  (*it), (*itp) );
                        Age newAge =  Age(this, minAge, position);
                        newAge.setExtra(true);
                        ages.insert( itp, newAge);
                    }
                    break;
                }
                else
                { // both ages are above minAge
                    break;
                }
                
                it ++;
                itp ++;
                if(itp == ages.end())
                {
                    if((*it).getValue() < minAge)
                        (*it).setEffValue(minAge);
                    break;
                }
            } // while(true)
        }
        else
        { // ages are in decreasing order from right to left;
            itp = ages.end()-1;
            it = itp - 1;
            while(true)
            {
                if((*it).getValue() < minAge)
                { // both it and itp are below minAge
                    (*itp).setEffValue(minAge);
                }
                else if((*itp).getValue() < minAge)
                { // limit is between it and itp
                    (*itp).setEffValue(minAge);
                    if(_insertAges)
                    {
                        double position =
                        linearPosition( minAge,  (*it), (*itp) );
                        Age newAge = Age(this, minAge, position);
                        newAge.setExtra(true);
                        ages.insert( itp, newAge);
                    }
                    break;
                }
                else
                { // both ages are above minAge
                    break;
                }
                
                it --;
                itp --;
                if(itp == ages.begin())
                {
                    if((*itp).getValue() < minAge)
                        (*itp).setEffValue(minAge);
                    break;
                }
            } // while(true)
        }
    } // ages.size() > 1
    else
    {
        ages[0].setEffValue(minAge);
    }
}
// ----------------------------------------
void
Earth::lock()
{
    _mutex.lock();
}

void
Earth::unlock()
{
    _mutex.unlock();
}
// --------------------------------------------------
void
Earth::saveState()
{
    EarthState state;
    
    lock();
    state.tMyr = Earth::getTimeMyr();
    state.ageMa = _ageMa;
    state.T = _T;
    state.Q = _QtotTW;
    state.ScontRatio = _cont_surfRatio;
    
    state.selectedInterface = -1;
    for(unsigned int i=0; i<_elements.size(); i++)
    {
        if(_elements[i]->getClass() == "Ridge")
            state.interfaces.push_back(make_pair("Ridge", make_pair(_elements[i]->getPosition(), 0.0)));
        if(_elements[i]->getClass() == "LeftSubduction" || _elements[i]->getClass() == "RightSubduction")
            state.interfaces.push_back(make_pair(_elements[i]->getClass(), make_pair(_elements[i]->getPosition(), ((Subduction*)_elements[i])->getDepth())));
        if(_elements[i]->getClass() == "Staple")
            state.interfaces.push_back(make_pair("Staple", make_pair(_elements[i]->getPosition(), 0.0)));
        
        if(_elements[i]->isSelected())
            state.selectedInterface = i;
    }
    
    state.selectedCell = -1;
    vector<Cell*> cells = getCells();
    for(unsigned int i=0; i<cells.size(); i++)
    {
        state.cells.push_back(make_pair( make_pair(cells[i]->getStartAngle(), cells[i]->getStopAngle()),
                                        make_pair(cells[i]->getOrientation(), cells[i]->getSpeedFactor()) ) );
        if(cells[i]->isSelected())
            state.selectedCell = i;
    }
    
    state.selectedPlate = -1;
    state.selectedSection = -1;
    for(unsigned int i=0; i<_plates.size(); i++)
    {
        state.plates.push_back(make_pair(_plates[i]->getRightElement()->getPosition(), _plates[i]->getLeftElement()->getPosition()));
        
        if(_plates[i]->isSelected())
            state.selectedPlate = i;
        
        if(state.selectedPlate == (int)i)
        {
            vector<PlateSection*>& sections = _plates[i]->accessSections();
            for(unsigned int j=0; j<sections.size(); j++)
            {
                state.sections.push_back(make_pair(sections[j]->getRightElement()->getPosition(),
                                                   sections[j]->getLeftElement()->getPosition()));
                if(sections[j]->isSelected())
                    state.selectedSection = j;
            } // sections[j]
        } // selected
        
        if(_drawAgeOceans)
        {
            state.drawAgeOceans = true;
            vector<PlateSection*> sections = _plates[i]->accessSections();
            for(unsigned int j=0; j<sections.size(); j++)
            {
                if(!sections[j]->isContinental())
                {
                    sections[j]->setAgeGroups();
                    vector<AgeGroup> groups = sections[j]->getAgeGroups();
                    for(unsigned int k=0; k<groups.size(); k++)
                        state.ageOceans.push_back(make_pair(groups[k].ageMax,
                                                            make_pair(groups[k].rightPosition,
                                                                      groups[k].leftPosition)
                                                            )
                                                  );
                }
            }
        }
        else
            state.drawAgeOceans = false;
    } // plates[i]
    
    state.selectedContinent = -1;
    for(unsigned int i=0; i<_continents.size(); i++)
    {
        state.continents.push_back(make_pair(_continents[i]->canBreak(),
                                             make_pair(_continents[i]->getRightPosition(),
                                                       _continents[i]->getLeftPosition())));
        if(_continents[i]->isSelected())
            state.selectedContinent = i;
        
        _continents[i]->setAgeGroups();
        vector<AgeGroup> groups = _continents[i]->getAgeGroups();
        for(unsigned int j=0; j<groups.size(); j++)
            state.ageContinents.push_back(make_pair(groups[j].ageMax,
                                                    make_pair(groups[j].rightPosition,
                                                              groups[j].leftPosition)
                                                    )
                                          );
    }
    
    
    _state = state;
    unlock();
}
// --------------------------------------------------
EarthState
Earth::getState()
{
    saveState();
    return _state;
}

/* =======================================================
 Outputs
 ======================================================= */
void
Earth::_writeOutputs()
{
    if(getTimeMyr() > _nextTimeWriteLogs)
    {
        logToFile(true);
        _nextTimeWriteLogs += _writeLogs;
    }
    
    // --- writing ages and positions of elements
    if(_timestep == 0 || getTimeMyr() > _nextTimeWriteAges)
    {
        // AGES --------------------------
        ofstream ageFile;
        ostringstream ossAge;
        ossAge << _workspace_ << "ages/ages"
        << getTimeString() << ".log";
        
        ageFile.open(ossAge.str().c_str());
        if(ageFile.is_open())
        {
            getAllAges();
            if(_timestep == 0)
                computeSeaLevel();
            
            ageFile << "# Age = "
            << setprecision(4) << fixed << _ageMa << " Ma\n";
            ageFile << "# T = " << _T << " K\n";
            ageFile << "# tauSub = " << _tau_sub << " Myr\n";
            ageFile << "# Qtot = " << _QtotTW << " TW\n";
            ageFile << "# sealevel(HSCM,PM95,PM125) = " << _sealevel_HSCM << " \t"
            << _sealevel_PM95 << " \t"
            << _sealevel_PM125 << " m\n";
            ageFile << "# ridgeDepth(HSCM,PM95,PM125) = " << _ridgeDepth_HSCM << "\t"
            << _ridgeDepth_PM95 << "\t"
            << _ridgeDepth_PM125 << " m\n";
            ageFile << "# VbelowRidge(HSCM,PM95,PM125) = "
            << scientific << _VbelowRidge_HSCM << " \t"
            << _VbelowRidge_PM95 << " \t"
            << _VbelowRidge_PM125 << " m^3\n";
            ageFile << "# ----------------------------------------\n";
            ageFile << "# 1-Position(degrees) \t"
            << "2-Age(Myr) \t"
            << "3-U(cm/yr) \t"
            << "4-heat flux(mW/m^2) \t"
            << "5-bathymetry_HSCM(m) \t"
            << "6-bathymetry_PM95(m) \t"
            << "7-bathymetry_PM125(m) \t"
            << "8-Effective age(Myr) \t\n";
            for(unsigned int i=0; i<_ages.size(); i++)
            {
                pair<double, double> bathyPM =
                computeBathyPM(_ages[i].getEffValue());
                ageFile << setw(12) << setprecision(9) << fixed << _ages[i].getPosition()
                << "\t" << setw(10) << setprecision(5) << fixed << _ages[i].getValue() + (double)_ages[i].getOriginContCounter()
                << "\t" << setw(10) << _ages[i].getU()
                << "\t" << setw(10) << computeHeatFlux(_ages[i].getEffValue()) * 1.0E3
                << "\t" << setw(10) << computeBathyHSCM(_ages[i].getEffValue())
                << "\t" << setw(10) << bathyPM.first
                << "\t" << setw(10) << bathyPM.second
                << "\t" << setw(10) << _ages[i].getEffValue()
                << "\n";
            }
            
            ageFile.close();
        }
        else
        {
            ostringstream ossErr;
            ossErr << "Error while opening " << ossAge.str().c_str() << endl;
            cerr << ossErr.str();
            if(_outputFile.is_open())
                _outputFile << ossErr.str() << flush;
        }
        
        
        // POSITIONS --------------------
        ofstream posFile;
        ostringstream ossPos;
        ossPos << _workspace_ << "ages/positions"
        << getTimeString() << ".log";
        
        posFile.open(ossPos.str().c_str());
        
        if(posFile.is_open())
        {
            posFile << "# Age = " << _ageMa
            << "\t T = " << _T << "\n";
            _writeGeoElementsPositions(posFile);
            
            posFile.close();
        }
        else
        {
            ostringstream ossErr;
            ossErr << "Error while opening " << ossPos.str().c_str() << endl;
            cerr << ossErr.str();
            if(_outputFile.is_open())
                _outputFile << ossErr.str() << flush;
        }
        
        
        // PLATES ----------
        ofstream plateFile;
        ostringstream ossPlates;
        ossPlates << _workspace_ << "ages/plates"
        << getTimeString() << ".log";
        
        plateFile.open(ossPlates.str().c_str(), ios::out | ios::trunc);
        if(plateFile.is_open())
        {
            plateFile << "# Age = " << _ageMa
            << "\t T = " << _T << "\n";
            
            if(_write_plate_forces)  // for debugging mainly
            {
                for(unsigned int i=0; i<_plates.size(); i++)
                {
                    plateFile << setw(16) << _plates[i]->getRightElement()->getShortId()
                    << setw(10) << _plates[i]->getRightElement()->getPosition()
                    << setw(10) << _plates[i]->getRightElement()->getLeftAge().getValue()
                    << setw(10) << _plates[i]->getLeftElement()->getShortId()
                    << setw(10) << _plates[i]->getLeftElement()->getPosition()
                    << setw(10) << _plates[i]->getLeftElement()->getRightAge().getValue()
                    << setprecision(3) << fixed
                    << setw(10) << _plates[i]->getU()
                    << setprecision(3) << scientific
                    << setw(11) << _plates[i]->accessForces()->getRidgePush() / 1.0E12
                    << setw(11) << _plates[i]->accessForces()->getSlabPull() / 1.0E12
                    << setw(11) << _plates[i]->accessForces()->getSlabSuction() / 1.0E12
                    << setw(11) << _plates[i]->accessForces()->getMantleDrag() * Earth::force0 / 1.0E12
                    << setw(11) << _plates[i]->accessForces()->getViscousShear() * Earth::force0 / 1.0E12
                    << setw(11) << _plates[i]->accessForces()->getBending() * Earth::force0 / 1.0E12
                    << setw(11) << _plates[i]->accessForces()->getEtaH()
                    << setw(11) << _plates[i]->accessForces()->getEtaRight()
                    << setw(11) << _plates[i]->accessForces()->getEtaLeft()
                    << fixed
                    << "\n";
                }
            } // _write_plate_forces
            else
            { // write width, type (oceanic, continental, mixed) and velocity of plates
                for(unsigned int i=0; i<_plates.size(); i++)
                {
                    plateFile << setw(10) << _plates[i]->getRightElement()->getShortId()
                    << setprecision(4) << fixed << setw(11) << _plates[i]->getRightElement()->getPosition()
                    << setw(10) << _plates[i]->getLeftElement()->getShortId()
                    << setprecision(4) << fixed << setw(11) << _plates[i]->getLeftElement()->getPosition()
                    << setw(11) << _plates[i]->getU()
                    << setw(11) << _plates[i]->computeLength();
                    // type of plate
                    if(_plates[i]->isOceanic())
                        plateFile << setw(13) << "oceanic";
                    else if(_plates[i]->isContinental())
                        plateFile << setw(13) << "continental";
                    else
                        plateFile << setw(13) << "mixed";
                    // type of right and left elements
                    plateFile << setw(8) << _plates[i]->getType()
                    << endl;
                }
            }
            
            plateFile.close();
        }
        else
        {
            ostringstream ossErr;
            ossErr << "Error while opening " << ossPlates.str().c_str() << endl;
            cerr << ossErr.str();
            if(_outputFile.is_open())
                _outputFile << ossErr.str() << flush;
        }
        
        _nextTimeWriteAges += _writeAges;
        
    } // timestep || getTime for writing Ages and Positions
}
// --------------------------------------------------
void 
Earth::_writeGeoElementsPositions(ofstream& posFile)
{ // CG
    vector<GeoElement*> elements;
    
    for(unsigned int i=0; i<_elements.size(); i++)
        elements.push_back(_elements[i]);
    
    // add continents' extremities:
    for(unsigned int i=0; i<_continents.size(); i++)
    {
        elements.push_back(_continents[i]->getRightExtremity());
        elements.push_back(_continents[i]->getLeftExtremity());
    }
    
    sortElements(elements);
    
    for(unsigned int i=0; i<elements.size(); i++)
    {
        posFile << setw(12) << setprecision(3) << fixed << elements[i]->getPosition()
        << setw(10) << elements[i]->getShortId()
        << setw(12) << elements[i]->getU();
        
        if(elements[i]->isA("Ridge"))
        {
            Ridge* ridge = (Ridge*)elements[i];
            posFile << "\t" << boolalpha << ridge->isActive();
        }
        else if(elements[i]->isA("Subduction"))
        {
            Subduction* subduction = (Subduction*)elements[i];
            //vector<pair<double, double> > thickness = subduction->getThickness();
            
            posFile << "\t" << setw(8) << subduction->getDepth() / 1.0E3;
            if(subduction->isA("LeftSubduction"))
                posFile << "\t" << setw(8) << subduction->getRightAge().getValue()
                << "\t" << setw(8) << subduction->computePlateThickness(RIGHT) / 1.0E3;
            else // right sub
                posFile << "\t" << setw(8) << subduction->getLeftAge().getValue()
                << "\t" << setw(8) << subduction->computePlateThickness(LEFT) / 1.0E3;
            
            double mean = subduction->getMeanThickness(); // also computes _thicknessAtDetph
            posFile << "\t" << setw(8) << mean / 1.0E3
            << "\t" << setw(8) << subduction->getThicknessAtDepth() / 1.0E3;
            
            if(subduction->isA("LeftSubduction"))
            {
                if(subduction->getLeftNeighborExtremity())
                {
                    ContinentExtremity* extremity = subduction->getLeftNeighborExtremity();
                    posFile << "\t"
                    << extremity->getShortId();
                }
            }
            else // right sub
            {
                if(subduction->getRightNeighborExtremity())
                {
                    ContinentExtremity* extremity = subduction->getRightNeighborExtremity();
                    posFile << "\t"
                    << extremity->getShortId();
                }
            }
            
            
        } // elements[i]->isA("Subduction")
        else if( elements[i]->isA("ContinentExtremity"))
        {
            // add the other continent extremity for the same continent
            Continent* continent = ((ContinentExtremity*)elements[i])->getContinent();
            if(elements[i]->isA("LeftContinentExtremity"))
                posFile << "\t" << setw(10)
                << continent->getRightExtremity()->getShortId();
            else
                posFile << "\t" << setw(10)
                << continent->getLeftExtremity()->getShortId();
            
            if(((ContinentExtremity*)elements[i])->getNeighbor())
            {
                posFile << "\t" << setw(10)
                << ((ContinentExtremity*)elements[i])->getNeighbor()->getShortId();
            } //elements[i]->isA("ContinentExtremity") with neighbor
        }
        posFile << "\n";
    }
}

// =======================================================
// conversions ----------------------------------------
double
Earth::deg_to_m(double ddeg)
{ // converts distance in degrees to distance in m
    return ddeg / deg_per_rad * _earthRadius;
}
double
Earth::deg_to_cm(double ddeg)
{ // converts distance in degrees to distance in cm
    return deg_to_m(ddeg) * 100.0;
}

double 
Earth::m_to_deg(double dm)
{ // converts distance in m to distance in deg
    return dm * deg_per_rad / _earthRadius;
}

double 
Earth::cm_to_deg(double dcm)
{ // converts distance in cm to distance in deg
    return m_to_deg(dcm * 1.0E-2);
}

pair<double, double> 
Earth::deg_to_m(pair<double, double> l)
{
    double coeff = _earthRadius / deg_per_rad;
    double l1 = l.first * coeff;
    double l2 = l.second * coeff;
    return make_pair(l1, l2);
}

double
Earth::deg_to_squareM(double ddeg)
{ // converts distance in deg to surface (m^2)
    return ddeg / 360.0 * _surfTot;
}
// =======================================================
// FINAL 
// -------------------------------------------------------
bool
Earth::endAge_is_reached()
{
    bool isReached = false;
    if(getTimeMyr() > 0.0 && _ageMa < _endAge)
    {
        isReached = true;
        cout << "*********************************************\n";
        cout << "Final age is reached: " << _ageMa << endl;
        cout << "*********************************************\n";
    }
    
    return isReached;
}
// --------------------------------------------------
bool
Earth::T_final_is_reached()
{
    bool isReached = false;
    if(getTimeMyr() > 0.0 && _T < _T_final)
    {
        isReached = true;
        cout << "*********************************************\n";
        cout << "Final temperature is reached: " << _T << " K" << endl;
        cout << "*********************************************\n";
    }
    
    return isReached;
}
// --------------------------------------------------
bool
Earth::final_condition_is_reached()
{
    bool isReached = false;
    if(!_T_cond_reached) {
        _T_cond_reached = T_final_is_reached();
    }
    if(!_age_cond_reached) {
        _age_cond_reached = endAge_is_reached();
    }

    if(_age_T_and_condition)
    {
        if(_T_cond_reached && _age_cond_reached)
            isReached = true;
    }
    else
    {
        if(_T_cond_reached || _age_cond_reached)
            isReached = true;
    }
    return isReached;
}
// --------------------------------------------------
void 
Earth::diagnosis()
{
    vector<double> ages = getOceanAges();
    sort(ages.begin(), ages.end());
    
    cout << setw(6) << (int)(getTimeMyr() + .5) << " Myr "
    << setw(6) << (int)(_ageMa + .5) << " Ma "
    << " T=" << setw(9) << fixed << setprecision(3) << _T
    << "\t age_max = " << ages[ages.size()-1]
    << "\t #plates: " << _plates.size();
    if(_continents.size() != 0)
        cout << "\t #continents: " << _continents.size();
    cout << endl;
}
// --------------------------------------------------
void
Earth::checkHeatFlux()
{
    _T = 1600.0;
    updateTData();
    
    while(true)
    {
        vector<Age> ages;
        double position = 0.0;
        double value = 0.0;
        while(true)
        {
            Age newAge = Age(this, value, position);
            ages.push_back(newAge);
            
            position += 1.0;
            value += 1.0;
            if(value > 200.0)
                break;
        }
        
        effectiveAges(ages);
        
        stringstream sTemp;
        sTemp << (int)_T;
        string filename =
        _workspace_ + "logs/heatflux_T" + sTemp.str() + ".log";
        ofstream heatfluxFile;
        heatfluxFile.open(filename.c_str(), ios::out | ios::trunc);
        if(heatfluxFile.is_open())
        {
            for(unsigned int i=0; i<ages.size(); i++)
                heatfluxFile << ages[i].getPosition()
                << "\t" << ages[i].getValue()
                << "\t" << ages[i].getEffValue()
                << "\t" << computeHeatFlux(ages[i].getValue()) * 1.0E3
                << "\t" << computeHeatFlux(ages[i].getEffValue()) * 1.0E3
                << endl;
            heatfluxFile.close();
        }
        
        _T += 100.0;
        if(_T > 1900.0)
            break;
        
    }
}
// --------------------------------------------------
void
Earth::checkFixedConfiguration()
{
    if(!_fixed_configuration)
        return;
    
    bool fixedPossible = true;
    
    for(unsigned int i=0; i<_plates.size(); i++)
    {
        if(_plates[i]->hasType("RS"))  // ridge and active subduction
        {
            vector<PlateSection*>& sections = _plates[i]->accessSections();
            if(sections.size() > 1)
            {
                for(unsigned int j=0; j<sections.size(); j++)
                { /* fixed is not possible if one section is continental
                   (the continent has to move with the plate) */
                    if(sections[j]->isContinental())
                    {
                        fixedPossible = false;
                        break;
                    }
                }
            }
        }
        else
        {
            /* fixed is possible only if plate is an upper plate
             for subduction on both sides */
            fixedPossible = _plates[i]->hasType("UU");
            if(!fixedPossible)
                break;
        }
    }
    
    if(!fixedPossible)
    {
        cout << "**************************************************" << endl;
        cout << "ERROR: fixedConfiguration is not possible with this configuration" << endl;
        cout << "**************************************************" << endl;
        exit(EXIT_FAILURE);
    }
}
// --------------------------------------------------
double
Earth::_randomAge(double ageIn, double noise)
{
    /* picks an age in (age +/- noise) */
    double ageOut = ageIn;
    
    if(ageIn < noise)
        ageOut = _random->nextDouble(ageIn,
                                     ageIn + noise);
    else
        ageOut = _random->nextDouble(ageIn - noise,
                                     ageIn + noise);
    
    return ageOut;
}
// --------------------------------------------------
