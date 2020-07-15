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
//          modifications: cecile.grigne@univ-brest.fr

#ifndef EARTH_H_
#define EARTH_H_

#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>
#include <utility>
#include <algorithm>
#include <stdio.h>
#include <stdlib.h>
#include <tuple>

#include "boost/thread/mutex.hpp"

// --- MACMA
#include "MACMA/physics/enum_struct.h"
#include "MACMA/physics/cell.h"
#include "MACMA/physics/ridge.h"
#include "MACMA/physics/leftSubduction.h"
#include "MACMA/physics/rightSubduction.h"
#include "MACMA/physics/staple.h"
#include "MACMA/physics/plate.h"
#include "MACMA/physics/plateSection.h"
#include "MACMA/physics/continent.h"
#include "MACMA/physics/leftContinentExtremity.h"
#include "MACMA/physics/rightContinentExtremity.h"
#include "MACMA/physics/warmingZone.h"
#include "MACMA/physics/marker.h"
#include "MACMA/physics/collision.h"
#include "MACMA/physics/activeMargin.h"
#include "MACMA/physics/force.h"

// --- Utils
#include "MACMA/utils/random.h"

using std::cout;
using std::vector;
using std::string;
using std::map;
using std::pair;
using std::ofstream;
using std::stringstream;

// --- Global
extern string _workspace_;

/* global (free) functions */
bool areEqual(double, double);
bool areEqual(double, double, double);
bool isStrictlyWithin(double, double, double);
bool isStrictlyLess(double, double);
bool isStrictlyLess(double, double, double);
bool isStrictlyMore(double, double);
bool isStrictlyMore(double, double, double);
bool isWithinOrEqual(double, double, double);
bool isWithinOrEqual(double, double, double, Direction);
double getMiddle(double, double);
double computeAbsDistance(double, double);
double computeDistance(double, double);
double computeRightLeftDistance(double, double);
pair<double, double> getMeanSD(vector<double>);
vector<double> getWhiskers(vector<double>);
double yr_to_sec(double);
double Myr_to_sec(double);
double sec_to_Myr(double);
vector<pair<double, double> > getDistribution(vector<double>, double);
double lin_Pos_to_Y(double, double, double, double, double);
double lin_Y_to_Pos(double, double, double, double, double);
void changePosition(double &pos, double);
double addPosition(double, double);
string doubleToString(double, unsigned int);
string doubleToString(double, unsigned int, string prefix, string suffix);
Direction doubleToDirection(double, double);
string idSuffixe(string);
/*
 ====================================================
 Class Earth
 ====================================================
 */

typedef struct _EarthState{
    double tMyr;
    double ageMa;
    double T;
    double Q;
    double ScontRatio;
    
    vector<pair<string, pair<double, double> > > interfaces;
    int selectedInterface;
    
    vector<pair<pair<double, double> , pair<double, double> > > cells;
    int selectedCell;
    
    vector<pair<double, double> > plates;
    int selectedPlate;
    
    vector<pair<double, double> > sections;
    int selectedSection;
    
    vector<pair<bool, pair<double, double> > > continents;
    int selectedContinent;
    
    bool drawAgeOceans;
    vector<pair<double, pair<double, double> > > ageOceans;
    
    vector<pair<double, pair<double, double> > > ageContinents;
    
} EarthState;


class Earth
{
public:
    // --- Initialization
    Earth(double dt);
    ~Earth();
    
    typedef enum _ContinentOpeningMode
    {
        TWO_SUBDUCTIONS,
        EXTENSION
    } ContinentOpeningMode;
    
    typedef enum _DepthInMantle
    {
        UPPER_MANTLE,
        WHOLE_MANTLE
    } DepthInMantle;
    
    enum SubductionInitiationMode
    {
        BRITTLE,
        CONVECTIVE,
        CONSTANT
    };
    
	enum SSCMode 
	{
		SSCCONVECTIVE,
		SSCCONSTANT
	};
	
    enum EventType
    {
        COLLISION,
        RIFT
    };
    
    struct InitSubduction
    {
        SubductionInitiationMode mode;
        bool continents;
        bool staples;
        bool upperPlates;
        bool ageRatioCriterion;
        double ratioAgeValue;
        double upperPlateMinAge;
        bool reverse;
        double tauSub_randomNoise;
        bool randomTauSub;
        bool alwaysOne;
    };
    
    enum DepletionMode
    {
        DEPLETED,
        PRIMITIVE,
        ADAPTATIVE
    };
    
    struct Contact // useful for collisions
    {
        double dt;
        GeoElement* rightElement;
        GeoElement* leftElement;
    };
    
    // initialisation - - - - -
    void prepareFiles();
    void defaultParameters();
    void initParameters();
    void initFixedParameters();
    void initTimeParameters();
    void initPhysicalParameters();
    void initModelParameters();
    void makeReady();
    void formerTimeSetting();
    void checkHeatFlux();
    void checkFixedConfiguration();
    
    // Model parameters - - - -
    inline DepthInMantle getSlabPullDepth() { return _slabPullDepth; }
    inline DepthInMantle getViscousShearDepth() { return _viscousShearDepth; }
    
    inline SubductionInitiationMode getSubductionMode() { return _initSubduction.mode; }
    void setSubductionMode();
    void setSubductionPlace();
    
	inline SSCMode getSSCMode() { return _SSCMode; }
	void setSSCMode();
    void setLimitedThickening();
    void setFixedConfiguration();
    
    inline bool middleBreakup() { return _middleBreakup; }
    inline double getBreakupPosition() { return _breakupPosition; }
    
    void setDepletionMode();
    inline DepletionMode getDepletionMode() { return _depletionMode; }
    
    // - - - - - -
    
    inline ContinentOpeningMode getContinentOpeningMode() {return _continentOpeningMode; }
    inline void setContinentOpeningMode(ContinentOpeningMode mode) {_continentOpeningMode = mode; }
    void setContinentOpeningMode(int twoSub, int extension);
    
    virtual void clear();
    virtual void clearPlates();
    
    // --- Static
    static double zero;
    static double loop;
    static double tinyDouble;
    static double ageEarth;
    static string getId(string className);	// Instance counter and id generator
    static double getTime();    // Simulated time in year
    static double getTimeMyr(); // Simulated time in Myr
    static int getTimestep();
    static double contAgeGroupStep;
    static double contAgeGroupMax;
    static double contAgeGroupMin;
    static double oceanAgeGroupStep;
    static double oceanAgeGroupMax;
    static double oceanAgeGroupMin;
    
    // static variables for dimensionless equations
    static double force0, eta0, u0;
    // other statics
    static double gasConstant;
    static double creationDistance;
    static double positionPrecision;
    static double agePrecision;
    static double coeff_plate_thickness;
    static double coeff_T_plate;
    
    
    // --- Structure
    virtual void checkNeighbors();
    virtual void fillStructure(bool optional = true);
    virtual void resetContinents();
    virtual void updateCells();
    
    // --- Computing
    inline double getAgeMa() { return _ageMa; }
    virtual void compute();
    virtual void updateT();
	void initTData();
    void updateTData();
    void updateViscosities();
    void updateDensities();
    void updateTPlate();
    void updateTauC();
    void updateTauSub();
    void initSlabs();
    inline double getT() { return _T; }
    inline void setT(double T) { _T = T; }
    inline double getEtaM() {return _eta_m; }
    inline double getEtaUm() {return _eta_um; }
    inline double getEtaPl() {return _eta_pl; }
    inline double getEtaAst() {return _eta_ast; }
    inline double getEtaSubCont() {return _eta_subcont; }
    inline double getTauSSC() {return _tau_ssc; }
    inline double getTauSub() {return _tau_sub; }
    inline double getRadioHeatTW() { return _radioHeatTW; }
    inline double getQtotTW() { return _QtotTW; }
    virtual double computedTdt();
    virtual double computeQmeanOcean();
    virtual double computeQmeanCont();
    virtual double computeMeanPlateThickness();
    virtual double computeSeafloorProduction();
    virtual double computeSlabFlux(bool optional = false);
    
    void completeForces();
    
    void computeThicknessesMantle();
    double computeSlabPull(double, double);
    
    double computeHeatFlux(double);
    virtual void updateQtot();
    
    double linearAge(double, Age, Age);
    double linearPosition(double, Age, Age);
    
    double computeBathyHSCM(double);
    pair<double, double> computeBathyPM(double);
    double computeDepthBelowRidgePM(double, double);
    double _sumForBathyPM(double, double, unsigned int);
    virtual void computeSeaLevel();
    
    virtual void initConcentrations();
    virtual void updateConcentrations();
    virtual void updateReferenceConcentrations();
    
    virtual void updateRadioHeatTW();
    virtual double computeRadioHeatTW(vector<Concentration>);
    virtual double radioHeatWperKg();
    
    virtual double computeEta(double eta_p, double E);
    virtual double computeRho(double rho_p, double alpha, double T, double Tp);
    virtual double compute_T_plate(double T);
    virtual void updateSurfaces();
    virtual void updateMasses();
    
    virtual double thickness_to_age(double);
    
    virtual double offsetPosition(GeoElement*, double, Direction);
    
    
    /* conversions */
    double deg_to_m(double);
    double deg_to_cm(double);
    double deg_to_squareM(double);
    pair<double, double> deg_to_m(pair<double, double>);
    double m_to_deg(double);
    double cm_to_deg(double);
    
    // --- Elements
    inline vector<GeoElement*>& accessElements() { return _elements; }
    virtual Ridge* createRidge(double position);
    virtual Subduction* createSubduction(double position, Direction direction);
    virtual Staple* createStaple(double position);
    virtual void insertElement(GeoElement* element);
    virtual int findElement(GeoElement* element);
    virtual GeoElement* findElement(const char* name);
    virtual bool removeElement(GeoElement* element);
    virtual void sortElements();
    virtual void sortElements(vector<GeoElement*>& elements);
    virtual vector<GeoElement*> findInterfaces();
    virtual void updateElementsVelocity();
    virtual void updatePositions();
    virtual void moveAges();
    virtual void effectiveAges(vector<Age>& ages);
    virtual void maxAgesLimit(vector<Age>& ages);
    virtual void minAgesLimit(vector<Age>& ages);
    
    // --- Plates
    inline vector<Plate*>& accessPlates(){ return _plates; }
    virtual Plate* createPlate(vector<PlateSection*>& sections);
    virtual Plate* findPlate(const char* plateName);
    virtual int findPlate(Plate* plate);
    virtual bool removePlate(Plate* plate);
    
    // bool for special configurations ----------------
    inline bool noSubduction() {return _noSubduction; }
    bool onlyOneSection();
        
    // --- Cells and other drawings
    virtual vector<Cell*> getCells();
    virtual Cell* findCell(const char* name);
    virtual int findCell(Cell* cell);
    inline bool drawCells() { return _drawCells; }
    inline void setDrawCells(bool draw) { _drawCells = draw; }
    inline bool drawAgeOceans() { return _drawAgeOceans; }
    
    
    // --- Continents
    virtual Continent* createContinent(double position,
                                       double length, string origin);
    virtual Continent* findContinent(const char* name);
    virtual int findContinent(Continent* continent);
    virtual bool removeContinent(Continent* continent);
    inline vector<Continent*>& accessContinents() { return _continents; }
    
	virtual RightContinentExtremity* createRightContinentExtremity(Continent* continent);
	virtual LeftContinentExtremity* createLeftContinentExtremity(Continent* continent);
	
	
    // --- Properties
    inline bool isReady() { return _ready; }
    inline void setReady(bool ready) { _ready = ready; }
    inline bool isRunning() { return _running; }
    inline void setRunning(bool running) { _running = running; }
    inline void setRunStarted(bool run) { _runStarted = run; }
    inline bool runStarted() { return _runStarted; }
    
    inline void setDt(double dt) { _dt = dt; }
    inline double getDt() { return _dt; }
    inline void setComputeEta(bool computeEta) { _computeEta = computeEta; }
    
    virtual bool T_final_is_reached();
    virtual bool endAge_is_reached();
    virtual bool final_condition_is_reached();
    
    //  -- special properties for testing or writing (CG)
    inline bool fixedConfiguration() {return _fixed_configuration; }
    inline bool isWritingContinents() { return _write_continents; }
    // - - -
    inline bool slabPull_rho_depends_on_T() {return _slabPull_rho_depends_on_T; }
    inline bool ridgePush_rho_depends_on_T() {return _ridgePush_rho_depends_on_T; }
    
    inline double get_T_m_init() {return _T_m_init; }
    inline double get_tau_sub_p() {return _tau_sub_p;}
    inline double get_tau_ssc_p() {return _tau_ssc_p; }
    inline double get_V_sink_p() {return _V_sink_p; }
    inline double get_F_lim() {return _F_lim; }
    inline double get_subcont_warming_H() {return _subcont_warming_H; }
    inline double get_R_min() {return _R_min; }
    inline double get_Qmax() {return _Qmax; }
    inline double get_min_plate_thick() {return _min_plate_thick; }
    inline double get_age_min_thickness() { return _age_min_thickness; }
    inline double get_T_p() {return _T_p; }
    inline double get_eta_m_p() {return _eta_m_p; }
    inline double get_eta_um_p() {return _eta_um_p; }
    inline double get_eta_pl_p() {return _eta_pl_p; }
    inline double get_eta_ast_p() {return _eta_ast_p; }
    inline double get_eta_subcont_p() {return _eta_subcont_p; }
    inline double get_E_m() {return _E_m; }
    inline double get_E_um() {return _E_um; }
    inline double get_E_pl() {return _E_pl; }
    inline double get_thick_ast() {return _thick_ast; }
    inline double get_thick_subcont() {return _thick_subcont; }
    inline double get_thick_continent() {return _thick_continent; }
    inline double get_k_continent() { return _k_continent; }
    inline double get_k_ocean() { return _k_ocean; }
    inline bool growingContinents() { return _continental_growth; }
    inline double get_contGrowthCoeff() { return _contGrowthCoeff; }
    inline double get_rhoCont() { return _rho_cont; }
    inline double get_cont_thickness() { return _cont_thickness; }
    
    // fixed parameters not appearing in graphic interface
    inline double get_T_surf() {return _T_surf; }
    inline double get_d() {return _d; }
    inline double get_D() {return _D; }
    inline double get_g() {return _g; }
    inline double get_kappa() {return _kappa; }
    inline double get_C_p() {return _C_p; }
    inline double get_rho_um() {return _rho_um; }
    inline double get_rho_um_p() {return _rho_um_p; }
    inline double get_rho_pl() {return _rho_pl; }
    inline double get_rho_pl_p() {return _rho_pl_p; }
    inline double get_rho_seawater() {return _rho_seawater; }
    inline double get_alpha_um() {return _alpha_um; }
    inline double get_alpha_pl() {return _alpha_pl; }
    inline double get_radius() {return _earthRadius; }
    
    // writing and time parameters
    inline double getWriteImages() {return _writeImages; }
    inline double getWriteXml() {return _writeXml; }
    inline double getStartAge() {return _startAge; }
    inline double getEndAge() {return _endAge; }
    inline double getStartTime() {return _startTime; }
    inline double getEndTime() {return _endTime; }
    inline double get_T_final() { return _T_final; }
    
    inline double get_resolution() {return _resolution; }
    inline double get_courant() {return _courant; }
    
    // -- others
    inline double getDepthFactor() {return _depthFactor; }
    inline void setDepthFactor(double d) {_depthFactor = d; }
    
    
    // --- importance of the different forces
    inline void setCoeffRidgePush(double coeff) { _coeffRidgePush = coeff; }
    inline void setCoeffSlabPull(double coeff) { _coeffSlabPull = coeff; }
    inline void setCoeffSlabSuction(double coeff) { _coeffSlabSuction = coeff; }
    inline void setCoeffMantleDrag(double coeff) { _coeffMantleDrag = coeff; }
    inline void setCoeffViscousShear(double coeff) { _coeffViscousShear = coeff; }
    inline void setCoeffBending(double coeff) { _coeffBending = coeff; }
    
    inline double getCoeffRidgePush() { return _coeffRidgePush; }
    inline double getCoeffSlabPull() { return _coeffSlabPull; }
    inline double getCoeffSlabSuction() { return _coeffSlabSuction; }
    inline double getCoeffMantleDrag() { return _coeffMantleDrag; }
    inline double getCoeffViscousShear() { return _coeffViscousShear; }
    inline double getCoeffBending() { return _coeffBending; }
    
    // to avoid extra computations at each time step
    inline double get_preFactor_ridgePush() { return _preFactor_ridgePush; }
    inline double get_preFactor_heatFlow() { return _preFactor_heatFlow; }
    inline double get_preFactor_bathyHSCM() { return _preFactor_bathyHSCM; }
    
    inline double getZM95() { return _zm95; }
    inline double getZM125() { return _zm125; }
    // -- computation (plates' velocities)
    vector<double> solverTridiagMatrix(vector<Tridiag>);
    
    // --- Graphics
    inline bool needRefresh() { return _needRefresh; }
    inline void setNeedRefresh(bool need) { _needRefresh = need; }
    
    // --- Parameters
    virtual void setTimeParameter(string name, double value);
    virtual double getTimeParameter(string name);
    virtual void setFixedParameter(string name, double value);
    virtual double getFixedParameter(string name);
    virtual void setPhysicalParameter(string name, double value);
    virtual double getPhysicalParameter(string name);
    virtual void setModelParameter(string name, double value);
    virtual double getModelParameter(string name);
    
    inline map<string, double> accessTimeParameters() { return _timeParameters; }
    inline map<string, double> accessFixedParameters() { return _fixedParameters; }
    inline map<string, double> accessPhysicalParameters() { return _physicalParameters; }
    inline map<string, double> accessModelParameters() { return _modelParameters; }
    virtual void changeTimeParameters();
    virtual void changePhysicalParameters();
    
    inline void setFixedTimestep(bool fixed) {_fixed_timestep = fixed; }
    inline bool isFixedTimestep() {return _fixed_timestep; }
    
    // --- Logs
    virtual void logToFile(bool);
    virtual string getTimeString();
    
    // --- Events
    vector<Contact> findCollisions(double& dt);
    virtual void treatCollisions(vector<Contact>);
    virtual void correctAges();
    
    // --- Ages
    virtual void initAges();
    vector<double> getOceanAges();
    void getAllAges();
    void updateAgesDistribution();
    vector<pair<double, double> > getAgesDistribution();
    
    vector<pair<double, double> > setAgeSlices(double, double, double);
    void setContinentalAgeSlices();
    inline vector<pair<double, double> > getContinentalAgeSlices() { return _continentalAgeSlices; }
    void setOceanAgeSlices();
    inline vector<pair<double, double> > getOceanAgeSlices() { return _oceanAgeSlices; }
    
    inline bool insertAges() { return _insertAges; }
    
    double _randomAge(double, double);
    
    // --- State
    virtual void lock();
    virtual void unlock();
    virtual void saveState();
    virtual EarthState getState();
    
    // --- Diagnosis
    void diagnosis();
    inline void setShowEvents(bool shw) { _showEvents = shw; }
    inline bool showEvents() { return _showEvents; }
    
    // --- thicknesses of the different layers in the mantle
    struct Thickness {
        double subcont;
        double contAst;
        double oceanAst;
        double contUM;
        double oceanUM;
        double LM;
    };
    inline Thickness getThicknesses() { return _thicknesses; }
    
    // --- markers
    bool removeMarker(Marker*);
    bool changeActiveMargin(Subduction*, string status);
    bool changeActiveMargin(ContinentExtremity*, string status);
    
    // --- miscellaneous
    inline Random* getRandom() { return _random; }
    
protected:
    // --- Static
    static double _time;
    static unsigned int _timestep;
    static unsigned int _instanceCounter;
    
    static unsigned int _ridgeCounter;
    static unsigned int _subductionCounter;
    static unsigned int _stapleCounter;
    static unsigned int _plateCounter;
    static unsigned int _plateSectionCounter;
    static unsigned int _cellCounter;
    static unsigned int _continentCounter;
    static unsigned int _continentExtremityCounter;
    static unsigned int _warmingZoneCounter;
    static unsigned int _collisionCounter;
    static unsigned int _activeMarginCounter;
    
    // --- Computing
    virtual void _timeLoop(double t, double dt);
    double _adaptDt();
    double _ageMa;
    
    /* ----------
     EVENTS
     ---------- */
    // --- collision
    virtual void _collisionContinentSubduction(GeoElement*, GeoElement*, Direction);
    virtual Continent* _collisionContinentContinent(Continent*, Continent*);
    
    // --- subduction
    virtual void _checkDiving();
    virtual void _checkAlwaysSubduct();
    virtual void _setNewRightSubduction(Subduction*, LeftContinentExtremity*);
    virtual void _setNewLeftSubduction(Subduction*, RightContinentExtremity*);
    virtual void _setNewSubduction(Subduction*, GeoElement*);
    
    
    virtual bool _divingStaple();
    virtual bool _divingContinentExtremity();
    virtual bool _divingUpperPlate();
    
    virtual bool _realDivingSubduction(Subduction*, GeoElement*, Direction, bool checkNearInterface = true);
    
    virtual void _verifySubductions();
    
    // --- ridge
    virtual void _checkOceanOpening();
    virtual double _openOcean(Continent* continent);
    
    virtual void _checkRidgeCreation();
    virtual void _ridgeCreation(Subduction* subduction);
    
    virtual void _verifyRidges();
    virtual void _activateAllRidges();
    
    // --- others
    virtual void _checkChangePlateTectonics();
    virtual void _markersChangeContinents(Continent*, Continent*);
    virtual void _markersChangeContinents(Continent*, Continent*, Continent*,
										  double, double);
    // --- computation
    void _computeVelocities();
    void _solveVelocities();
    void _simpleVelocities();
    bool _verifyVelocities();
    
    // --- Outputs
    void _writeGeoElementsPositions(ofstream&);
    void _writeOutputs();
    
    // --- Properties
    ContinentOpeningMode _continentOpeningMode;
    double _dt;
    bool _ready;
    bool _running;
    bool _runStarted; // set to true only once
    
    bool _fixed_timestep;
    bool _age_T_and_condition;
    bool _T_cond_reached;
    bool _age_cond_reached;
    
    // Temperature, viscosities and densities ----
    bool _computeEta;
    double _T;
    double _QtotTW;
    double _eta_m;
    double _eta_um;
    double _eta_pl;
    double _eta_ast;
    double _eta_subcont;
    double _rho_um;
    double _rho_pl;
    double _T_plate;
    double _T_plate_p;
    double _rho_pl_p;
    
    // concentration and heat radioactive elements
    vector<Concentration> _concentrations;  // order is U235, U238, Th232, K40
    vector<Concentration> _depleted;
    vector<Concentration> _primitive;
    
    double _cont_surfRatio;
    double _radioHeatTW;
    
    // thicknesses of different mantle layers ----
    Thickness _thicknesses;
    
    // surfaces of oceans and continents ---------
    double _surfTotOcean;
    double _surfTotContinent;
    double _surfTot;
    
    // --- limit ages ----------
    double _tau_sub;
    double _tau_ssc;
    
    // --- ages -----------
    vector<Age> _ages;
    bool _insertAges;
    vector<pair<double, double> > _continentalAgeSlices;
    vector<pair<double, double> > _oceanAgeSlices;
    
    // --- Physical parameters
    double _T_m_init;
    double _tau_sub_p;
    double _tau_ssc_p;
    double _V_sink_p;
    double _F_lim;
    double _subcont_warming_H;
    double _R_min;
    double _Qmax;
    double _min_plate_thick;
    double _T_p;
    double _eta_m_p;
    double _eta_um_p;
    double _eta_pl_p;
    double _eta_ast_p;
    double _eta_subcont_p;
    double _E_m;
    double _E_um;
    double _E_pl;
    double _thick_ast;
    double _thick_subcont;
    double _thick_continent;
    double _k_continent;
    double _k_ocean;
    double _rho_um_p;
    double _DeltaRho_p;
    double _alpha_um;
    double _alpha_pl;
    bool _continental_growth;
    double _contGrowthCoeff;
    
    // -- fixed parameters not appearing in graphic interface:
    double _T_surf;
    
    double _d;
    double _D;
    
    double _g;
    
    double _kappa;
    double _C_p;
    
    double _earthRadius;
    double _earthPerimeter;
    double _earthMass;
    
    // for radioactive heat
    double _heat235U;
    double _heat238U;
    double _heat232Th;
    double _heat40K;
    double _heatTW_235U;
    double _heatTW_238U;
    double _heatTW_232Th;
    double _heatTW_40K;
    double _present_cont_surfRatio;
    double _cont_thickness;
    double _rho_cont;
    double _present_mantleMass;
    double _mantleMass;
    double _primitive_mantleMass;
    double _present_continentMass;
    double _continentMass;
    double _depletionRatio;
    bool _fixed_depletion;
    
    // for sea level variations
    double _rho_seawater;
    double _present_oceanVolume;
    double _present_ridgeDepth;
    double _present_VaboveRidge;
    
    double _ridgeDepth_HSCM;
    double _VbelowRidge_HSCM;
    double _VaboveRidge_HSCM;
    double _meanDepthBelowRidge_HSCM;
    double _sealevel_HSCM;
    
    double _ridgeDepth_PM95;
    double _VbelowRidge_PM95;
    double _VaboveRidge_PM95;
    double _meanDepthBelowRidge_PM95;
    double _sealevel_PM95;
    
    double _ridgeDepth_PM125;
    double _VbelowRidge_PM125;
    double _VaboveRidge_PM125;
    double _meanDepthBelowRidge_PM125;
    double _sealevel_PM125;
    
    // --- time + writing
    double _startAge;
    double _endAge;
    double _startTime;
    double _endTime;
    double _writeLogs;
    double _writeAges;
    double _writeImages;
    double _writeXml;
    double _timeStep;
    double _T_final;
    
    double _nextTimeWriteLogs;
    double _nextTimeWriteAges;
    
    // --- resolution
    double _resolution;
    double _courant;
    
    // --- others (graphics)
    double _depthFactor;
    
    // --- Booleans for special testing of the code
    bool _fixed_configuration;
    bool _write_continents;
    bool _write_plate_forces; // mainly for debugging
    
    // --- importance of the different forces (in force balance)
    double _coeffRidgePush;
    double _coeffSlabPull;
    double _coeffSlabSuction;
    double _coeffMantleDrag;
    double _coeffViscousShear;
    double _coeffBending;
    
    // --- pre-factors for forces, heat flux and bathymetry
    double _preFactor_ridgePush;
    double _preFactor_heatFlow;
    double _preFactor_bathyHSCM;
    double _preFactor_bathyPM2;
    double _zm95;
    double _zm125;
    
    // --- densities dependence on temperature
    bool _slabPull_rho_depends_on_T;
    bool _ridgePush_rho_depends_on_T;
    
    // --- radioactive isotopes concentrations
    DepletionMode _depletionMode;
    double _constant_heat_value;
    bool _constant_heating;
    
    // --- subductions
    InitSubduction _initSubduction;
    
    // --- depth for forces related to slabs
    DepthInMantle _slabPullDepth, _viscousShearDepth;
    
    // --- small scale convection
	SSCMode _SSCMode;
    bool _limitedThickening;
    
    // --- continental breakup
    bool _middleBreakup;
    double _breakupPosition;
    
    // --- minimum thickness of the lithosphere
    double _age_min_thickness; // computed when initializing
    
    // --- Parameters maps
    map<string, double> _timeParameters;
    map<string, double> _fixedParameters;
    map<string, double> _physicalParameters;
    map<string, double> _modelParameters;
    
    // -- Interfaces
    vector<GeoElement*> _elements;
    
    // --- Plates
    vector<Plate*> _plates;
    
    // --- Special positions' markers
    vector<Marker*> _markers;
    
    // special configurations :
    bool _noSubduction;
    bool _noSubductionPrec;
    double _timespan_noSubduction;
    bool _inactiveRidges;
    
    // --- Continents
    vector<Continent*> _continents;
    
    // --- Graphics and outputs
    bool _needRefresh;
    bool _drawCells;
    bool _showEvents;
    bool _drawAgeOceans;
    
    // --- Logs
    ofstream _logFile;
    ofstream _eventsFile;
    ofstream _ageFile;
    ofstream _velocityFile;
    ofstream _forcesFile;
    ofstream _outputFile;
    ofstream _geometryFile;
    ofstream _continentsFile;
    ofstream _testFile;
    ofstream _radioContFile;
    ofstream _bathymetryFile;
    ofstream _platesFile;
    ofstream _elementsFile;
    ofstream _collisionsFile;
    ofstream _activeMarginsFile;
	
    // --- State
    EarthState _state;
    
    // age distribution
    double _ageSlice;
    vector<pair<double, double> > _agesDistribution;
    
    
    // -- miscellaneous
    Random* _random;
    unsigned int _randomSeed;
    
    boost::mutex _mutex;
};

#endif 
