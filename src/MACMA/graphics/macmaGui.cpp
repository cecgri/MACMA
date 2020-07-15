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

/* august 2015: creation of two inherited classes 
 from Macma, with and without GUI.
 Macma is the parent class.
 author: cecile.grigne@univ-brest.fr (CG) */

#include "MACMA/graphics/macmaGui.h"
#include "MACMA/graphics/earthTimer.h"

#define MAX_FRAME 5
#define GRAPHICS_DT 10
#define RECORD_DT 500

using namespace std;
// ==================================================
MacmaGui::MacmaGui(QWidget* parent) : QMainWindow(parent), Macma()
{
	
	askForWorkspace();
	
	_startT = 0.0;
	
	_wait_usec = 50000.0; // in microsec
	
	_earth = new Earth(200.0);
	
	_frameCounter = 0;
	
	_selectedElement = NULL;
	_selectedPlate = NULL;
	_selectedSection = NULL;
	_selectedContinent = NULL;
	
	_thread = boost::thread(&MacmaGui::_earthThread, this);
	
	_earth->setDrawCells(true);
	_earth->setShowEvents(true); // TODO: put this as a choice.
	//  _earth->setShowEvents(false);  
	
	_loadWidgets();
	_connectSignals();
	
	putParametersInUI();
	setMACMAState(PARAMETERS);
	
}

MacmaGui::~MacmaGui()
{
}
// ==================================================
// --- Widgets
void
MacmaGui::_loadWidgets()
{
	_ui.setupUi(this);
	
	_ui.actionRun->setEnabled(false);
	
	_ui.startTime->setEnabled(false);
	_ui.endTime->setEnabled(false);
	_ui.interfaceParametersGroup->setEnabled(false);
	_ui.nextInterfaceButton->setEnabled(false);
	_ui.interfacePosition->setEnabled(false);
	_ui.interfaceLeftAge->setEnabled(false);
	_ui.interfaceRightAge->setEnabled(false);
	_ui.slabDepth->setEnabled(false);
	_ui.saveInterfaceParameters->setEnabled(false);
	_ui.removeInterfaceButton->setEnabled(false);
	
	_ui.plateGroup->setEnabled(false);
	_ui.continentGroup->setEnabled(false);
	_ui.continentParametersGroup->setEnabled(false);
	
	_graphics = new Graphics(this, _earth);
	_ui.splitter->addWidget(_graphics);
	_ui.splitter->setStretchFactor(0, 2);
	_ui.splitter->setStretchFactor(1, 20);
}
// --------------------------------------------------
void
MacmaGui::_connectSignals()
{
	// Toolbar
	connect(_ui.actionOpen, SIGNAL(triggered()), this, SLOT(_slotOpen()));
	connect(_ui.actionSave, SIGNAL(triggered()), this, SLOT(_slotSave()));
	connect(_ui.actionRun, SIGNAL(triggered()), this, SLOT(_slotRun()));
	
	// ALL EARTH - - - - - - - - - - - -
	connect(_ui.earthButton, SIGNAL(clicked()), this, SLOT(_slotEarthButton()));
	// EARTH : specific for time:
	connect(_ui.startAge, SIGNAL(valueChanged(double)), this, SLOT(_slotChangeStartAge(double)));
	connect(_ui.endAge, SIGNAL(valueChanged(double)), this, SLOT(_slotChangeEndAge(double)));
	
	// wait time in ms to avoid jumpls
	connect(_ui.waitms, SIGNAL(valueChanged(int)), this, SLOT(_slotChangeWaitms(int)));
	
	
	// MODEL - - - - - - - - - - - - 
	connect(_ui.modelButton, SIGNAL(clicked()), this, SLOT(_slotModelButton()));
	
	// Interfaces
	connect(_ui.interfaceList, SIGNAL(currentRowChanged(int)), this, SLOT(_slotInterfaceList(int)));
	connect(_ui.addInterfaceButton, SIGNAL(clicked()), this, SLOT(_slotAddInterface()));
	connect(_ui.removeInterfaceButton, SIGNAL(clicked()), this, SLOT(_slotRemoveInterface()));
	connect(_ui.interfacePosition, SIGNAL(valueChanged(double)), this, SLOT(_slotInterfacePosition(double)));
	connect(_ui.slabDepth, SIGNAL(valueChanged(int)), this, SLOT(_slotSlabDepth(int)));
	connect(_ui.saveInterfaceParameters, SIGNAL(clicked()), this, SLOT(_slotSaveInterfaceParameters()));
	connect(_ui.nextInterfaceButton, SIGNAL(clicked()), this, SLOT(_slotNextInterface()));
	
	// Plates
	connect(_ui.plateList, SIGNAL(currentRowChanged(int)), this, SLOT(_slotPlateList(int)));
	connect(_ui.sectionList, SIGNAL(currentRowChanged(int)), this, SLOT(_slotSectionList(int)));
	connect(_ui.continentList, SIGNAL(currentRowChanged(int)), this, SLOT(_slotContinentList(int)));
	connect(_ui.addContinentButton, SIGNAL(clicked()), this, SLOT(_slotAddContinent()));
	connect(_ui.removeContinentButton, SIGNAL(clicked()), this, SLOT(_slotRemoveContinent()));
	connect(_ui.continentPosition, SIGNAL(valueChanged(double)), this, SLOT(_slotContinentPosition(double)));
	connect(_ui.continentLength, SIGNAL(valueChanged(double)), this, SLOT(_slotContinentLength(double)));
	connect(_ui.saveContinentParameters, SIGNAL(clicked()), this, SLOT(_slotSaveContinentParameters()));
	
	// Timers
	_graphicsTimer = new QTimer(this);
	connect(_graphicsTimer, SIGNAL(timeout()), this, SLOT(_slotUpdateGraphics()));
	_graphicsTimer->start(GRAPHICS_DT);
	
	// _recordTimer = new QTimer(this);
	// connect(_recordTimer, SIGNAL(timeout()), this, SLOT(_slotRecordTimeout()));
	
	_imagesTimer = new EarthTimer(_earth);
	connect(_imagesTimer, SIGNAL(timeReached()), this, SLOT(_slotWriteImages()));
	
	_xmlTimer = new EarthTimer(_earth);
	connect(_xmlTimer, SIGNAL(timeReached()), this, SLOT(_slotWriteXml()));
	
}
// --------------------------------------------------

/* ----------
 |  SLOTS |
 ---------- */

void
MacmaGui::_slotOpen()
{
	QString filename = QFileDialog::getOpenFileName(this, "Open ...", ".","MACMA configuration files (*.macma)");
	if(!filename.isEmpty())
		open(filename); /* makes _state = PLATES */
}

void
MacmaGui::_slotSave()
{
	QString filename = QFileDialog::getSaveFileName(this, "Save ...", ".", "MACMA configuration file (*.macma)");
	if(!filename.isEmpty()) 
		save(filename);
}

void
MacmaGui::_slotRun()
{
	toggleSuspended();
}

// void
// MacmaGui::_slotRecord()
// {
//   if(_recordTimer->isActive())
//     _recordTimer->stop();
//   else
//     _recordTimer->start(RECORD_DT);
// }

// void
// MacmaGui::_slotGraphs()
// {
//   _showGraphs = !_showGraphs;
//   if(_showGraphs)
//     _graphs->show();
//   else
//     _graphs->hide();
// }

void
MacmaGui::_slotEarthButton()
{ /* all parameters in the Earth tab 
   are saved and set when
   the button "save" is clicked */
	
	_earth->lock();
	
	// Time/Resolution/Writing - - - 
	_earth->setTimeParameter("startAge", _ui.startAge->value());
	_earth->setTimeParameter("endAge", _ui.endAge->value());
	_earth->setTimeParameter("T_final", _ui.T_final->text().toDouble());
	_earth->setTimeParameter("writeLogs", _ui.writeLogs->text().toDouble());
	_earth->setTimeParameter("writeAges", _ui.writeAges->text().toDouble());
	_earth->setTimeParameter("writeImages", _ui.writeImages->text().toDouble());
	_earth->setTimeParameter("writeXml", _ui.writeXml->text().toDouble());
	_earth->setTimeParameter("resolution", _ui.resolution->text().toDouble());
	_earth->setTimeParameter("courantNumber", _ui.courantNumber->text().toDouble());
	_earth->setTimeParameter("timeStep", _ui.timeStep->text().toDouble());
    _earth->setTimeParameter("age_T_and_condition",
                             (_ui.andConditionCheck->isChecked() ? 1.0 : 0.0));
	_earth->setTimeParameter("fixed_timestep",
							 (_ui.fixedTimestepCheck->isChecked() ? 1.0 : 0.0)); /* 0.0 if is unchecked,
																				  non zero otherwise */
	
	// Physical parameters - - - - - 
	_earth->setPhysicalParameter("T_m_init", _ui.T_m_init->text().toDouble());
	_earth->setPhysicalParameter("T_p", _ui.T_p->text().toDouble());
	_earth->setPhysicalParameter("Tau_sub_p", _ui.Tau_sub_p->text().toDouble());
	_earth->setPhysicalParameter("Tau_ssc_p", _ui.Tau_ssc_p->text().toDouble());
	_earth->setPhysicalParameter("F_lim", _ui.F_lim->text().toDouble());
	_earth->setPhysicalParameter("R_min", _ui.R_min->text().toDouble());
	_earth->setPhysicalParameter("subcont_warming_H", _ui.subcont_warming_thickness->text().toDouble());
	_earth->setPhysicalParameter("E_m", _ui.E_m->text().toDouble());
	_earth->setPhysicalParameter("E_um", _ui.E_um->text().toDouble());
	_earth->setPhysicalParameter("E_pl", _ui.E_pl->text().toDouble());
	_earth->setPhysicalParameter("V_sink_p", _ui.V_sink_p->text().toDouble());
	_earth->setPhysicalParameter("Qmax", _ui.Qmax->text().toDouble());
	_earth->setPhysicalParameter("min_plate_thick", _ui.min_plate_thickness->text().toDouble());
	_earth->setPhysicalParameter("eta_pl_p", _ui.eta_pl_p->text().toDouble());
	_earth->setPhysicalParameter("eta_m_p", _ui.eta_m_p->text().toDouble());
	_earth->setPhysicalParameter("eta_um_p", _ui.eta_um_p->text().toDouble());
	_earth->setPhysicalParameter("eta_ast_p", _ui.eta_ast_p->text().toDouble());
	_earth->setPhysicalParameter("eta_subcont_p", _ui.eta_subcont_p->text().toDouble());
	_earth->setPhysicalParameter("thick_ast", _ui.thick_ast->text().toDouble());
	_earth->setPhysicalParameter("thick_subcont", _ui.thick_subcont->text().toDouble());
	_earth->setPhysicalParameter("thick_continent", _ui.thick_continent->text().toDouble());
	_earth->setPhysicalParameter("k_ocean", _ui.k_ocean->text().toDouble());
	_earth->setPhysicalParameter("k_continent", _ui.k_continent->text().toDouble());
	_earth->setPhysicalParameter("rho_um_p", _ui.rho_um_p->text().toDouble());
	_earth->setPhysicalParameter("DeltaRho_p", _ui.DeltaRho_p->text().toDouble());
	_earth->setPhysicalParameter("alpha_um", _ui.alpha_um->text().toDouble());
	_earth->setPhysicalParameter("alpha_pl", _ui.alpha_pl->text().toDouble());
	
	_earth->unlock();
	
	_earth->initTimeParameters(); // TODO: really needed with new slots?
	_earth->initPhysicalParameters();
	
	if(_state < PARAMETERS)
		setMACMAState(PARAMETERS);
	
	_ui.interfaceGroup->setEnabled(true);
	_ui.interfaceList->setEnabled(true);
	_ui.interfaceCombo->setEnabled(true);
	_ui.addInterfaceButton->setEnabled(true);
	
	_ui.interfaceParametersGroup->setEnabled(false);
	//_ui.nextInterfaceButton->setEnabled(false);
	
}

void
MacmaGui::_slotChangeStartAge(double value)
{
	_earth->setTimeParameter("startAge", value);
	
	// endAge cannot be bigger than this startAge
	_ui.endAge->setMaximum(value);
	
	double start = Earth::ageEarth - value;
	_ui.startTime->setText(QString::number(start));
}

void
MacmaGui::_slotChangeEndAge(double value)
{
	_earth->setTimeParameter("endAge", value);
	
	// startAge cannot be smaller than this startAge
	_ui.startAge->setMinimum(value);
	
	double endtime = 
    Earth::ageEarth - _earth->getTimeParameter("endAge");
	_ui.endTime->setText(QString::number(endtime));
}

void
MacmaGui::_slotChangeWaitms(int i)
{
	_wait_usec = 1000.0 * (double)i;
}



void
MacmaGui::_slotModelButton()
{ /* all parameters in the Model tab 
   are saved and set when
   the button "save" is clicked */
	
	_earth->lock();
	
	// Coefficients for forces
	_earth->setModelParameter("coeffSlabPull", _ui.coeffSlabPull->text().toDouble());
	_earth->setModelParameter("coeffRidgePush", _ui.coeffRidgePush->text().toDouble());
	_earth->setModelParameter("coeffSlabSuction", _ui.coeffSlabSuction->text().toDouble());
	_earth->setModelParameter("coeffViscousShear", _ui.coeffViscousShear->text().toDouble());
	_earth->setModelParameter("coeffMantleDrag", _ui.coeffMantleDrag->text().toDouble());
	_earth->setModelParameter("coeffBending", _ui.coeffBending->text().toDouble());
	
	// Depth for SP and VS
	_earth->setModelParameter("slabPullWholeMantle", 
							  (_ui.SPWholeMantleCheck->isChecked() ? 1.0 : 0.0));
	_earth->setModelParameter("viscousShearWholeMantle", 
							  (_ui.VSWholeMantleCheck->isChecked() ? 1.0 : 0.0));
	
	
	// Subduction initiation mode
	_earth->setModelParameter("initMode_brittle",
							  (_ui.brittleCheck->isChecked() ? 1.0 : 0.0));
	_earth->setModelParameter("initMode_constant",
							  (_ui.constantCheck->isChecked() ? 1.0 : 0.0));
	_earth->setModelParameter("initMode_convective",
							  (_ui.convectiveCheck->isChecked() ? 1.0 : 0.0));
	
	
	// Subduction initiation place and critical age
	_earth->setModelParameter("initPlace_continents",
							  (_ui.initSubdContinentsCheck->isChecked() ? 1.0 : 0.0));
	_earth->setModelParameter("initPlace_staples",
							  (_ui.initSubdStaplesCheck->isChecked() ? 1.0 : 0.0));
	_earth->setModelParameter("initPlace_upperPlates",
							  (_ui.initSubdUpperPlatesCheck->isChecked() ? 1.0 : 0.0));
	
	_earth->setModelParameter("upperPlates_ageRatioCriterion",
							  (_ui.upperPlatesAgeRatioCheck->isChecked() ? 1.0 : 0.0));
	_earth->setModelParameter("upperPlates_ratioAgeValue",
							  _ui.upperPlatesAgeRatio->text().toDouble());
	
	_earth->setModelParameter("upperPlates_minimumAge",
							  _ui.minUpperPlatesAge->text().toDouble());
	
	_earth->setModelParameter("reverseSubduction",
							  (_ui.reverseSubductionCheck->isChecked() ? 1.0 : 0.0));
	
	_earth->setModelParameter("tauSub_randomNoise",
							  _ui.tauSubNoise->text().toDouble());
	
	_earth->setModelParameter("alwaysSubduct",
							  (_ui.alwaysSubCheck->isChecked() ? 1.0 : 0.0));
	
	// SSC
	_earth->setModelParameter("small_scale_convection",
							  (_ui.SSCCheck->isChecked() ? 1.0 : 0.0));
	_earth->setModelParameter("sscMode_constant",
							  (_ui.SSCConstantCheck->isChecked() ? 1.0 : 0.0));
	_earth->setModelParameter("sscMode_convective",
							  (_ui.SSCConvectiveCheck->isChecked() ? 1.0 : 0.0));
	
	// fixed plate boundaries
	_earth->setModelParameter("fixed_configuration",
							  (_ui.fixedConfigurationCheck->isChecked() ? 1.0 : 0.0));
	
	// radioactive heating
	_earth->setModelParameter("U_BSE", _ui.U_BSE->text().toDouble());
	_earth->setModelParameter("Th_BSE", _ui.Th_BSE->text().toDouble());
	_earth->setModelParameter("K_BSE", _ui.K_BSE->text().toDouble());
	_earth->setModelParameter("U_cont", _ui.U_cont->text().toDouble());
	_earth->setModelParameter("Th_cont", _ui.Th_cont->text().toDouble());
	_earth->setModelParameter("K_cont", _ui.K_cont->text().toDouble());
	
	_earth->setModelParameter("depletion_always",
							  (_ui.depletedMantleCheck->isChecked() ? 1.0 : 0.0));
	_earth->setModelParameter("depletion_never",
							  (_ui.primitiveMantleCheck->isChecked() ? 1.0 : 0.0));
	_earth->setModelParameter("depletion_adaptative",
							  (_ui.adaptativeDepletionCheck->isChecked() ? 1.0 : 0.0));
	
	_earth->setModelParameter("constant_internal_heating",
							  (_ui.constantHCheck->isChecked() ? 1.0 : 0.0));
	_earth->setModelParameter("constant_heat_value", _ui.constantH->text().toDouble());
	
	
	// continental breakup position
	_earth->setModelParameter("middle_breakup",
							  (_ui.middleBreakupCheck->isChecked() ? 1.0 : 0.0));
	_earth->setModelParameter("breakup_position", _ui.breakupPosWidth->value());
	
	// continental growth
	_earth->setModelParameter("continental_growth",
							  _ui.contGrowthCheck->isChecked() ? 1.0 : 0.0);
	_earth->setModelParameter("contGrowthCoeff", 
							  _ui.contGrowthCoeff->text().toDouble());
	
	
	
	// miscellaneous
	_earth->setModelParameter("randomSeed", _ui.randomSeed->value());
	
	
	_earth->unlock();
	
	_earth->initModelParameters();
}


void
MacmaGui::_slotInterfaceList(int row)
{
	_ui.interfacePosition->setEnabled(false);
	_ui.interfaceLeftAge->setEnabled(false);
	_ui.interfaceRightAge->setEnabled(false);
	_ui.slabDepth->setEnabled(false);
	_ui.saveInterfaceParameters->setEnabled(false);
	_ui.removeInterfaceButton->setEnabled(false);
	
	if(row != -1)
    {
		if(_selectedElement)
			_selectedElement->setSelected(false);
		_selectedElement = NULL;
		
		string interfaceName = _ui.interfaceList->currentItem()->text().toStdString();
		
		double position = 0.0;
		double slabDepth = 0.0;
		
		_earth->lock();
		_selectedElement = _earth->findElement(interfaceName.c_str());
		if(_selectedElement)
		{
			_selectedElement->setSelected(true);
			position = _selectedElement->getPosition();
			
			if(_selectedElement->isA("Subduction"))
				slabDepth = ((Subduction*)_selectedElement)->getDepth() / 1e3;
		}
		_earth->unlock();
		
		if(_selectedElement->isA("Ridge"))
		{
			_ui.interfacePosition->setEnabled(true);
			_ui.interfaceLeftAge->setEnabled(false);
			_ui.interfaceRightAge->setEnabled(false);
			_ui.slabDepth->setEnabled(false);
			_ui.saveInterfaceParameters->setEnabled(false);
			_ui.removeInterfaceButton->setEnabled(true);
		}
		else if(_selectedElement->isA("Staple"))
		{
			_ui.interfacePosition->setEnabled(true);
			_ui.interfaceLeftAge->setEnabled(true);
			_ui.interfaceRightAge->setEnabled(true);
			_ui.slabDepth->setEnabled(false);
			_ui.saveInterfaceParameters->setEnabled(true);
			_ui.removeInterfaceButton->setEnabled(true);
		}
		else if(_selectedElement->isA("Subduction"))
		{
			_ui.interfacePosition->setEnabled(true);
			_ui.interfaceLeftAge->setEnabled(true);
			_ui.interfaceRightAge->setEnabled(true);
			_ui.slabDepth->setEnabled(true);
			_ui.saveInterfaceParameters->setEnabled(true);
			_ui.removeInterfaceButton->setEnabled(true);
		}
		
		_ui.interfacePosition->setValue(position);
		_ui.slabDepth->setValue(slabDepth);
		_ui.interfaceLeftAge->setText(QString::number(_selectedElement->getLeftAge().getValue()));
		_ui.interfaceRightAge->setText(QString::number(_selectedElement->getRightAge().getValue()));
    }
	else
    {
		if(_selectedElement)
			_selectedElement->setSelected(false);
		_selectedElement = NULL;
		_ui.interfacePosition->setDecimals(2);
		_ui.interfacePosition->setValue(0.0);
		_ui.slabDepth->setValue(0);
		_ui.interfaceLeftAge->setText("0");
		_ui.interfaceRightAge->setText("0");
    }
}

void
MacmaGui::_slotAddInterface()
{
	int interfaceType = _ui.interfaceCombo->currentIndex();
	
	GeoElement* interface;
	_earth->lock();
	switch(interfaceType)
    {
		case 0:
		{
			interface = _earth->createRidge(0.0);
			break;
		}
			
		case 1:
		{
			interface = _earth->createSubduction(0.0, LEFT);
			break;
		}
			
		case 2:
		{
			interface = _earth->createSubduction(0.0, RIGHT);
			break;
		}
			
		case 3:
		{
			interface = _earth->createStaple(0.0);
			break;
		}
			
		default:
			interface = NULL;
			break;
    }
	_earth->unlock();
	
	_ui.interfaceList->addItem(QString(interface->getId().c_str()));
	_ui.interfaceList->setCurrentRow(_ui.interfaceList->count()-1);
	
	elementActivation(true);
	_ui.removeInterfaceButton->setEnabled(true);
	
	if(_state > INTERFACES)
		setMACMAState(INTERFACES);
	
	_ui.nextInterfaceButton->setEnabled(true);
}

void
MacmaGui::_slotRemoveInterface()
{
	string elementName = _ui.interfaceList->currentItem()->text().toStdString();
	_ui.interfaceList->takeItem(_ui.interfaceList->currentRow());
	unselectAll();
	
	_earth->lock();
	_earth->removeElement(_earth->findElement(elementName.c_str()));
	_earth->unlock();
	
	_ui.removeInterfaceButton->setEnabled(false);
	_ui.interfaceParametersGroup->setEnabled(false);
	
	_ui.interfacePosition->setValue(0);
	_ui.interfaceLeftAge->setText("0");
	_ui.interfaceRightAge->setText("0");
	_ui.slabDepth->setValue(0);
	
	if(_state > INTERFACES)
		setMACMAState(INTERFACES);
	
	if(_ui.interfaceList->count() < 1)
		_ui.nextInterfaceButton->setEnabled(false);
}

void
MacmaGui::_slotInterfacePosition(double position)
{
	if(_selectedElement)
    {
		string elementName = _ui.interfaceList->currentItem()->text().toStdString();
		
		_earth->lock();
		GeoElement* element = _earth->findElement(elementName.c_str());
		if(element)
		{
			element->setPosition(position);
			Age leftAge = Age(_earth, element->getLeftAge().getValue(), position);
			Age rightAge = Age(_earth, element->getRightAge().getValue(), position);
			element->setLeftAge(leftAge);
			element->setRightAge(rightAge);
		}
		_earth->sortElements();
		_earth->unlock();
		
		if(_state > INTERFACES)
			setMACMAState(INTERFACES); // set one level lower
    }
}

void
MacmaGui::_slotSlabDepth(int value)
{
	if(_ui.interfaceList->currentItem())
    {
		string elementName = _ui.interfaceList->currentItem()->text().toStdString();
		
		_earth->lock();
		GeoElement* element = _earth->findElement(elementName.c_str());
		if(element && element->isA("Subduction"))
			((Subduction*)element)->setDepth(value*1e3);
		_earth->unlock();
    }
}

void
MacmaGui::_slotSaveInterfaceParameters()
{
	if(_selectedElement)
    {
		Age leftAge = Age(_earth,
						  _ui.interfaceLeftAge->text().toDouble(),
						  _ui.interfacePosition->text().toDouble());
		Age rightAge = Age(_earth,
						   _ui.interfaceRightAge->text().toDouble(),
						   _ui.interfacePosition->text().toDouble());
		_earth->lock();
		_selectedElement->setLeftAge(leftAge);
		_selectedElement->setRightAge(rightAge);
		_earth->unlock();
		
		if(_state > INTERFACES)
			setMACMAState(INTERFACES);
		
    }
}

void
MacmaGui::_slotNextInterface()
{
	_ui.plateList->clear();
	
	_earth->lock();
	_earth->clearPlates();
	_earth->checkNeighbors();
	_earth->fillStructure(false);
	_earth->updateCells();
	_earth->unlock();
	
	setMACMAState(PLATES);
	
	_earth->lock();
	vector<Plate*> & plates = _earth->accessPlates();
	for(unsigned int i=0; i<plates.size(); i++)
    {
		_ui.plateList->addItem(QString(plates[i]->getId().c_str()));
    }
	_earth->unlock();
	
	elementActivation(false);
	_ui.interfaceList->setEnabled(true);
	
	plateActivation(true);
	
	continentActivation(true);
	_ui.addContinentButton->setEnabled(false);
	_ui.removeContinentButton->setEnabled(false);
	_ui.continentPosition->setEnabled(false);
	_ui.continentLength->setEnabled(false);
	_ui.continentLeftAge->setEnabled(false);
	_ui.continentRightAge->setEnabled(false);
	_ui.saveContinentParameters->setEnabled(false);
	
	_ui.actionRun->setEnabled(true);
	_ui.actionRecord->setEnabled(true);
	
	unselectAll();
}

void
MacmaGui::_slotPlateList(int row)
{
	if(row != -1)
    {
		_earth->lock();
		string name = _ui.plateList->currentItem()->text().toStdString();
		
		if(_selectedPlate)
			_selectedPlate->setSelected(false);
		_selectedPlate = NULL;
		
		_selectedPlate = _earth->findPlate(name.c_str());
		if(_selectedPlate)
		{
			_selectedPlate->setSelected(true);
			
			if(_selectedSection)
				_selectedSection->setSelected(false);
			_selectedSection = NULL;
			
			_ui.sectionList->clear();
			
			vector<PlateSection*> & sections = _selectedPlate->accessSections();
			for(unsigned int i=0; i<sections.size(); i++)
				_ui.sectionList->addItem(QString(sections[i]->getId().c_str()));
			
			if(_selectedContinent)
				_selectedContinent->setSelected(false);
			_selectedContinent = NULL;
			
			_ui.continentList->clear();
			
			continentActivation(true);
			_ui.sectionList->setEnabled(true);
			_ui.continentList->setEnabled(true);
			
			if(_state < READY)
				plateActivation(true);
		}
		_earth->unlock();
    }
	else
    {
		if(_selectedPlate)
			_selectedPlate->setSelected(false);
		_selectedPlate = NULL;
		
		_ui.sectionList->clear();
		_ui.continentList->clear();
    }
}

void
MacmaGui::_slotSectionList(int row)
{
	if(row != -1)
    {
		_earth->lock();
		string name = _ui.sectionList->currentItem()->text().toStdString();
		
		_selectedSection = _selectedPlate->findSection(name.c_str());
		if(_selectedSection)
		{
			_selectedSection->setSelected(true);
			
			_ui.continentList->clear();
			
			vector<Continent*> & continents = _earth->accessContinents();
			for(unsigned int i=0; i<continents.size(); i++)
			{
				if(continents[i]->getPlateSection() == _selectedSection)
					_ui.continentList->addItem(QString(continents[i]->getId().c_str()));
			}
			
			_ui.addContinentButton->setEnabled(true);
			
			double minPos = _selectedSection->getRightElement()->getPosition();
			double maxPos = _selectedSection->getLeftElement()->getPosition();
			
			if(maxPos <= minPos)
				maxPos += 360.0;
			
			_ui.continentPosition->setDecimals(2);
			_ui.continentPosition->setMinimum(minPos);
			_ui.continentPosition->setMaximum(maxPos);
			
			_ui.continentLength->setDecimals(2);
			_ui.continentLength->setMinimum(1.0);
			_ui.continentLength->setMaximum(maxPos-minPos);
		}
		_earth->unlock();
    }
	else
    {
		if(_selectedSection)
			_selectedSection->setSelected(false);
		_selectedSection = NULL;
		
		_ui.addContinentButton->setEnabled(false);
    }
}

void
MacmaGui::_slotContinentList(int row)
{
	if(row != -1)
    {
		string name = _ui.continentList->currentItem()->text().toStdString();
		
		_selectedContinent = _earth->findContinent(name.c_str());
		if(_selectedContinent)
		{
			_selectedContinent->setSelected(true);
			
			_ui.removeContinentButton->setEnabled(true);
			_ui.continentPosition->setEnabled(true);
			_ui.continentLength->setEnabled(true);
			_ui.continentLeftAge->setEnabled(true);
			_ui.continentRightAge->setEnabled(true);
			_ui.saveContinentParameters->setEnabled(true);
			
			_ui.continentPosition->setValue(_selectedContinent->getPosition());
			_ui.continentLength->setValue(_selectedContinent->getLength());
			_ui.continentLeftAge->setText(QString::number(_selectedContinent->getLeftExtremity()->getLeftAge().getValue()));
			_ui.continentRightAge->setText(QString::number(_selectedContinent->getRightExtremity()->getRightAge().getValue()));
		}
    }
	else
    {
		if(_selectedContinent)
			_selectedContinent->setSelected(false);
		_selectedContinent = NULL;
		
		_ui.removeContinentButton->setEnabled(false);
		_ui.continentPosition->setEnabled(false);
		_ui.continentLength->setEnabled(false);
		_ui.continentLeftAge->setEnabled(false);
		_ui.continentRightAge->setEnabled(false);
		_ui.saveContinentParameters->setEnabled(false);
		
		_ui.continentPosition->setDecimals(2);
		_ui.continentPosition->setValue(0.0);
		
		_ui.continentLength->setDecimals(2);
		_ui.continentLength->setValue(0.0);
		
		_ui.continentLeftAge->setText("0");
		_ui.continentRightAge->setText("0");
    }
}

void
MacmaGui::_slotAddContinent()
{
	_earth->lock();
	
	string strInit = "initial";
	Continent* continent =
    _earth->createContinent(_selectedSection->getRightElement()->getPosition(),
							0.25, strInit);
	_earth->createRightContinentExtremity(continent);
	_earth->createLeftContinentExtremity(continent);
	
	_ui.continentList->addItem(QString(continent->getId().c_str()));
	
	_ui.removeContinentButton->setEnabled(true);
	_earth->unlock();
	
	setMACMAState(INTERFACES);
}

void
MacmaGui::_slotRemoveContinent()
{
	if(_selectedContinent)
		_selectedContinent->setSelected(false);
	_selectedContinent = NULL;
	
	string name = _ui.continentList->currentItem()->text().toStdString();
	_ui.continentList->takeItem(_ui.continentList->currentRow());
	unselectAll();
	
	_ui.continentPosition->setEnabled(false);
	_ui.continentLength->setEnabled(false);
	_ui.continentLeftAge->setEnabled(false);
	_ui.continentRightAge->setEnabled(false);
	_ui.saveContinentParameters->setEnabled(false);
	
	_ui.continentPosition->setValue(0.0);
	_ui.continentLength->setValue(0.0);
	_ui.continentLeftAge->setText("0");
	_ui.continentRightAge->setText("0");
	
	_earth->lock();
	_earth->removeContinent(_earth->findContinent(name.c_str()));
	_earth->unlock();
	
	setMACMAState(INTERFACES);
}

void
MacmaGui::_slotContinentPosition(double value)
{
	if(_selectedContinent)
    {
		_earth->lock();
		_selectedContinent->setPosition(value);
		_earth->unlock();
		setMACMAState(INTERFACES);
		
		double maxPos = _selectedSection->getLeftElement()->getPosition();
		if(value > maxPos)
			maxPos += 360.0;
		
		double maxLength = maxPos - value;
		double currentLength = _selectedContinent->getLength();
		
		_ui.continentLength->setDecimals(2);
		_ui.continentLength->setMinimum(1.0);
		_ui.continentLength->setMaximum(maxLength);
		
		if(currentLength > maxLength)
			_ui.continentLength->setValue(_ui.continentLength->maximum());
    }
}

void
MacmaGui::_slotContinentLength(double value)
{
	if(_selectedContinent)
    {
		_earth->lock();
		_selectedContinent->setLength(value);
		
		setMACMAState(INTERFACES);
		_earth->unlock();
    }
}

void
MacmaGui::_slotSaveContinentParameters()
{
	_earth->lock();
	if(_selectedContinent)
    {
		_selectedContinent->setPlateSection(_selectedSection);
		
		Age leftAge = Age(_earth,
						  _ui.continentLeftAge->text().toDouble(),
						  fmod(_ui.continentPosition->text().toDouble() + 
							   _ui.continentLength->text().toDouble(),
							   360.0));  
		_selectedContinent->getLeftExtremity()->setLeftAge(leftAge);	
		
		Age rightAge = Age(_earth,
						   _ui.continentRightAge->text().toDouble(),
						   _ui.continentPosition->text().toDouble());
		_selectedContinent->getRightExtremity()->setRightAge(rightAge);
    }
	
	_earth->clearPlates();
	_earth->checkNeighbors();
	_earth->fillStructure(false);
	_earth->updateCells();
	_earth->unlock();
	
	setMACMAState(PLATES);
}

void
MacmaGui::_slotUpdateGraphics()
{
	_graphics->update();
}

void
MacmaGui::_slotRecordTimeout()
{
	record();
}

void
MacmaGui::_slotWriteImages()
{
	ostringstream oss;
	oss << _workspace_ << "images/MACMA_"
	<< _earth->getTimeString() << ".png";
	
	_graphics->saveToPNG(oss.str().c_str());
}

void
MacmaGui::_slotWriteXml()
{
	ostringstream oss;
	oss << _workspace_ << "logs/param_"
	<< _earth->getTimeString() << ".macma";
	
	writeXmlFile(oss.str());
}
// --------------------------------------------------
/* ----------
 | Utils  |
 ---------- */

void
MacmaGui::open(QString filename)
{
	openXmlFile(filename.toStdString()); /* makes _state = PLATES */
	
	earthActivation(true);
	elementActivation(true);
	plateActivation(true);
	continentActivation(true);
	
	_ui.addContinentButton->setEnabled(false);
	_ui.removeContinentButton->setEnabled(false);
	_ui.continentPosition->setEnabled(false);
	_ui.continentLength->setEnabled(false);
	
	_ui.actionRun->setEnabled(true);
	_ui.actionRecord->setEnabled(true);
	
	return;
}

void
MacmaGui::save(QString filename)
{
	saveXmlFile(filename.toStdString());
	return;
}

void
MacmaGui::record()
{
	string frameNumber = "";
	char buf[10];
	sprintf(buf, "%d", _frameCounter);
	int length = strlen(buf);
	
	for(int i=0; i<MAX_FRAME-length; i++)
		frameNumber+="0";
	frameNumber += string(buf);
	
	ostringstream oss;
	oss << _workspace_ + _recordFolder + "/" + "MACMA_" << frameNumber << ".png";
	
	_graphics->saveToPNG(oss.str().c_str());
	
	_frameCounter++;
}


void
MacmaGui::unselectAll()
{
	QItemSelectionModel* selectionModel;
	
	if(_selectedElement)
		_selectedElement->setSelected(false);
	_selectedElement = NULL;
	
	selectionModel = _ui.interfaceList->selectionModel();
	selectionModel->clearSelection();
	
	if(_selectedPlate)
		_selectedPlate->setSelected(false);
	_selectedPlate = NULL;
	
	selectionModel = _ui.plateList->selectionModel();
	selectionModel->clearSelection();
	
	if(_selectedSection)
		_selectedSection->setSelected(false);
	_selectedSection = NULL;
	
	selectionModel = _ui.sectionList->selectionModel();
	selectionModel->clearSelection();
	
	if(_selectedContinent)
		_selectedContinent->setSelected(false);
	_selectedContinent = NULL;
	
	selectionModel = _ui.continentList->selectionModel();
	selectionModel->clearSelection();
}

void
MacmaGui::toggleSuspended()
{
	if(_state < RUNNING)
    {
		if(_state < PLATES)
		{
			cerr << "ERROR: Run was launched before the structure was set"
			<< endl;
			exit(EXIT_FAILURE);
		}
		
		unselectAll();
		_ui.actionRun->setIcon(QIcon(":/icons/icons/pause.png"));
		
		earthActivation(false);
		elementActivation(false);
		plateActivation(false);
		continentActivation(false);
		
		if(_state < READY)
		{
			_earth->makeReady();
			setMACMAState(READY);
			cout << "STARTING RUN... " << endl;
		}
		else
			cout << "...RUNNING" << endl;
		
		setMACMAState(RUNNING);
		_earth->setRunning(true);
		
    }
	else
    {
		_ui.actionRun->setIcon(QIcon(":/icons/icons/play.png"));
		// if(_recordTimer->isActive())
		// 	_recordTimer->stop();
		
		_ui.timeGroup->setEnabled(true);
		
		setMACMAState(READY);
		_earth->setRunning(false);
		cout << "SUSPENDED..." << endl;
    }
}

double
MacmaGui::getTime()
{
	struct timeval tp;
	gettimeofday(&tp,NULL);
	return double(tp.tv_sec)*1e6 + double(tp.tv_usec);
}

EarthState
MacmaGui::getEarthState()
{
	EarthState state;
	
	_earth->lock();
	state.tMyr = _earth->getTimeMyr();
	state.ageMa = _earth->getAgeMa();
	state.T = _earth->getT();
	
	vector<GeoElement*>& elements = _earth->accessElements();
	for(unsigned int i=0; i<elements.size(); i++)
    {
		if(elements[i]->getClass() == "Ridge")
			state.interfaces.push_back(make_pair("Ridge", make_pair(elements[i]->getPosition(), 0.0)));
		if(elements[i]->getClass() == "LeftSubduction" || elements[i]->getClass() == "RightSubduction")
			state.interfaces.push_back(make_pair(elements[i]->getClass(), make_pair(elements[i]->getPosition(), ((Subduction*)elements[i])->getDepth())));
		if(elements[i]->getClass() == "Staple")
			state.interfaces.push_back(make_pair("Staple", make_pair(elements[i]->getPosition(), 0.0)));
    }
	state.selectedInterface = _earth->findElement(_selectedElement);
	
	state.selectedSection = -1;
	vector<Plate*>& plates = _earth->accessPlates();
	for(unsigned int i=0; i<plates.size(); i++)
    {
		state.plates.push_back(make_pair(plates[i]->getRightElement()->getPosition(), plates[i]->getLeftElement()->getPosition()));
		
		if(_selectedPlate == plates[i])
		{
			vector<PlateSection*>& sections = plates[i]->accessSections();
			for(unsigned int j=0; j<sections.size(); j++)
				state.sections.push_back(make_pair(sections[j]->getRightElement()->getPosition(), sections[j]->getLeftElement()->getPosition()));
			if(state.selectedSection == -1)
				state.selectedSection = plates[i]->findSection(_selectedSection);
		}
    }
	state.selectedPlate = -1;
	if(_selectedPlate)
		state.selectedPlate = _earth->findPlate(_selectedPlate);
	
	vector<Continent*>& continents = _earth->accessContinents();
	for(unsigned int i=0; i<continents.size(); i++)
		state.continents.push_back(make_pair(continents[i]->canBreak(), make_pair(continents[i]->getRightPosition(), continents[i]->getLeftPosition())));
	state.selectedContinent = _earth->findContinent(_selectedContinent);
	
	_earth->unlock();
	
	return state;
}

void
MacmaGui::lock()
{
	_mutex.lock();
}

void
MacmaGui::unlock()
{
	_mutex.unlock();
}

void
MacmaGui::earthActivation(bool activated)
{
	// in Time group
	_ui.timeGroup->setEnabled(activated);
	_ui.startTime->setEnabled(false);
	_ui.endTime->setEnabled(false); // always deactivated
	
	// Model
	_ui.forceBalanceGroup->setEnabled(activated);
	_ui.modelGroup->setEnabled(activated);
	
	// Physical parameters
	_ui.earthParametersGroup->setEnabled(activated);
}

void
MacmaGui::elementActivation(bool activated)
{
	_ui.interfaceGroup->setEnabled(activated);
	_ui.interfaceParametersGroup->setEnabled(activated);
	_ui.nextInterfaceButton->setEnabled(activated);
}

void
MacmaGui::plateActivation(bool activated)
{
	_ui.plateGroup->setEnabled(activated);
}

void
MacmaGui::continentActivation(bool activated)
{
	_ui.continentGroup->setEnabled(activated);
	_ui.continentParametersGroup->setEnabled(activated);
}

void
MacmaGui::askForWorkspace()
{
	QString path = QFileDialog::getExistingDirectory(this, tr("Select a workspace ..."),
													 QDir::homePath(),
													 QFileDialog::ShowDirsOnly
													 | QFileDialog::DontResolveSymlinks);
	if(!path.isEmpty())
    {
		string slog = path.toStdString() + "/logs";
		string sages = path.toStdString() + "/ages";
		string selements = path.toStdString() + "/elements";
		string simages = path.toStdString() + "/images"; 
		QDir dirlog(QString(slog.c_str()));
		QDir dirages(QString(sages.c_str()));
		QDir direlements(QString(selements.c_str()));
		QDir dirimages(QString(simages.c_str()));
		
		
		if(dirlog.exists())
			cerr << "Files in logs/ : " << dirlog.count()-2 << endl;
		else if(dirages.exists())
			cerr << "Files in ages/ : " << dirages.count()-2 << endl;
		else if(direlements.exists())
			cerr << "Files in elements/ : " << direlements.count()-2 << endl;
		else if(dirimages.exists())
			cerr << "Files in images/ : " << dirimages.count()-2 << endl;
		else
			cerr << "Empty workspace" << endl;
		if((dirlog.exists() && dirlog.count()>2) 
		   || (dirages.exists() && dirages.count()>2)
		   || (direlements.exists() && direlements.count()>2)
		   || (dirimages.exists() && dirimages.count()>2))
		{
			int answer2 = QMessageBox::warning(this, 
											   tr("MACMA"),
											   tr("You have selected : ") + path + tr("\n\n"
																					  "This workspace is not empty.\n"
																					  " All data in logs/, ages/, elements/ and images/\n"
																					  " will be erased.\n\n"
																					  "Are you sure you want to continue with this workspace ?"),
											   QMessageBox::Ok | QMessageBox::Cancel);
			if(answer2 == QMessageBox::Cancel)
				askForWorkspace();
			else
				prepareWorkspace(path.toStdString());
		}
		else
			prepareWorkspace(path.toStdString());
    }
	else
    {
		int answer = QMessageBox::critical(this, 
										   tr("MACMA"),
										   tr("No workspace selected, /tmp/MACMA will be used by default"),
										   QMessageBox::Ok | QMessageBox::Cancel);
		if(answer == QMessageBox::Cancel)
			askForWorkspace();
		else
		{
			QDir::root().QDir::mkpath("/tmp/MACMA");
			prepareWorkspace("/tmp/MACMA");
		}
    }
}

void
MacmaGui::prepareWorkspace(string path)
{
	_workspace_ = path + "/";
	
	string pathlog = path + "/logs";
	string pathages = path + "/ages";
	string pathelements = path + "/elements";
	string pathimages = path + "/images"; 
	cleanDirectory(pathlog);
	cleanDirectory(pathages);
	cleanDirectory(pathelements);
	cleanDirectory(pathimages);
	
	QDir dir(QString(path.c_str()));
	dir.mkdir("logs");
	dir.mkdir("ages");
	dir.mkdir("elements");
	dir.mkdir("images");
}

void
MacmaGui::cleanDirectory(string path)
{
	QDir dir(QString(path.c_str()));
	foreach(QFileInfo fileInfo, dir.entryInfoList()) 
    {
		if(fileInfo.fileName() != "." && fileInfo.fileName() != "..")
		{
			if(fileInfo.isDir())
			{
				cleanDirectory(fileInfo.absoluteFilePath().toStdString());
				dir.rmdir(fileInfo.absoluteFilePath());
			}
			else
				dir.remove(fileInfo.absoluteFilePath());
		}
    }
}
// --------------------------------------------------
void
MacmaGui::putParametersInUI()
{ 
	// PUT READ PARAMETERS IN GRAPHIC INTERFACE :
	_ui.startAge->setValue(_earth->getTimeParameter("startAge"));
	_ui.endAge->setValue(_earth->getTimeParameter("endAge"));
	_ui.T_final->setText(QString::number(_earth->getTimeParameter("T_final")));
	_ui.writeAges->setText(QString::number(_earth->getTimeParameter("writeAges")));
	_ui.writeLogs->setText(QString::number(_earth->getTimeParameter("writeLogs")));
	_ui.writeImages->setText(QString::number(_earth->getTimeParameter("writeImages")));
	_ui.writeXml->setText(QString::number(_earth->getTimeParameter("writeXml")));
	_ui.resolution->setText(QString::number(_earth->getTimeParameter("resolution")));
	_ui.courantNumber->setText(QString::number(_earth->getTimeParameter("courantNumber")));
	_ui.timeStep->setText(QString::number(_earth->getTimeParameter("timeStep")));
	_ui.fixedTimestepCheck->setChecked((_earth->getTimeParameter("fixed_timestep") != 0.0));
    _ui.andConditionCheck->setChecked((_earth->getTimeParameter("age_T_and_condition") != 0.0));
	
	// PARAMETERS ---------------
	_ui.T_m_init->setText(QString::number(_earth->getPhysicalParameter("T_m_init")));
	_ui.T_p->setText(QString::number(_earth->getPhysicalParameter("T_p")));
	_ui.Tau_sub_p->setText(QString::number(_earth->getPhysicalParameter("Tau_sub_p")));
	_ui.Tau_ssc_p->setText(QString::number(_earth->getPhysicalParameter("Tau_ssc_p")));
	_ui.F_lim->setText(QString::number(_earth->getPhysicalParameter("F_lim")));
	_ui.subcont_warming_thickness->setText(QString::number(_earth->getPhysicalParameter("subcont_warming_H")));
	_ui.V_sink_p->setText(QString::number(_earth->getPhysicalParameter("V_sink_p")));
	_ui.eta_m_p->setText(QString::number(_earth->getPhysicalParameter("eta_m_p")));
	_ui.eta_um_p->setText(QString::number(_earth->getPhysicalParameter("eta_um_p")));
	_ui.eta_pl_p->setText(QString::number(_earth->getPhysicalParameter("eta_pl_p")));
	_ui.eta_ast_p->setText(QString::number(_earth->getPhysicalParameter("eta_ast_p")));
	_ui.eta_subcont_p->setText(QString::number(_earth->getPhysicalParameter("eta_subcont_p")));
	_ui.E_m->setText(QString::number(_earth->getPhysicalParameter("E_m")));
	_ui.E_um->setText(QString::number(_earth->getPhysicalParameter("E_um")));
	_ui.E_pl->setText(QString::number(_earth->getPhysicalParameter("E_pl")));
	_ui.R_min->setText(QString::number(_earth->getPhysicalParameter("R_min")));
	_ui.Qmax->setText(QString::number(_earth->getPhysicalParameter("Qmax")));
	_ui.min_plate_thickness->setText(QString::number(_earth->getPhysicalParameter("min_plate_thick")));
	_ui.thick_ast->setText(QString::number(_earth->getPhysicalParameter("thick_ast")));
	_ui.thick_subcont->setText(QString::number(_earth->getPhysicalParameter("thick_subcont")));
	_ui.thick_continent->setText(QString::number(_earth->getPhysicalParameter("thick_continent")));
	_ui.k_ocean->setText(QString::number(_earth->getPhysicalParameter("k_ocean")));
	_ui.k_continent->setText(QString::number(_earth->getPhysicalParameter("k_continent")));
	_ui.rho_um_p->setText(QString::number(_earth->getPhysicalParameter("rho_um_p")));
	_ui.DeltaRho_p->setText(QString::number(_earth->getPhysicalParameter("DeltaRho_p")));
	_ui.alpha_um->setText(QString::number(_earth->getPhysicalParameter("alpha_um")));
	_ui.alpha_pl->setText(QString::number(_earth->getPhysicalParameter("alpha_pl")));
	
	// MODEL ---------------
	_ui.coeffSlabPull->setText(QString::number(_earth->getModelParameter("coeffSlabPull")));
	_ui.coeffRidgePush->setText(QString::number(_earth->getModelParameter("coeffRidgePush")));
	_ui.coeffSlabSuction->setText(QString::number(_earth->getModelParameter("coeffSlabSuction")));
	_ui.coeffMantleDrag->setText(QString::number(_earth->getModelParameter("coeffMantleDrag")));
	_ui.coeffViscousShear->setText(QString::number(_earth->getModelParameter("coeffViscousShear")));
	_ui.coeffBending->setText(QString::number(_earth->getModelParameter("coeffBending")));
	
	_ui.SPUpperMantleCheck->setChecked((_earth->getModelParameter("slabPullWholeMantle") == 0.0));
	_ui.SPWholeMantleCheck->setChecked((_earth->getModelParameter("slabPullWholeMantle") != 0.0));
	_ui.VSUpperMantleCheck->setChecked((_earth->getModelParameter("viscousShearWholeMantle") == 0.0));
	_ui.VSWholeMantleCheck->setChecked((_earth->getModelParameter("viscousShearWholeMantle") != 0.0));
	
	_ui.brittleCheck->setChecked((_earth->getModelParameter("initMode_brittle") != 0.0));
	_ui.convectiveCheck->setChecked((_earth->getModelParameter("initMode_convective") != 0.0));
	_ui.constantCheck->setChecked((_earth->getModelParameter("initMode_constant") != 0.0));
	
	_ui.initSubdContinentsCheck->setChecked((_earth->getModelParameter("initPlace_continents") != 0.0));
	_ui.initSubdStaplesCheck->setChecked((_earth->getModelParameter("initPlace_staples") != 0.0));
	_ui.initSubdUpperPlatesCheck->setChecked((_earth->getModelParameter("initPlace_upperPlates") != 0.0));
	_ui.upperPlatesAgeRatioCheck->setChecked((_earth->getModelParameter("upperPlates_ageRatioCriterion") != 0.0));
	_ui.upperPlatesTauSubCheck->setChecked((_earth->getModelParameter("upperPlates_ageRatioCriterion") == 0.0));
	_ui.minUpperPlatesAge->setText(QString::number(_earth->getModelParameter("upperPlates_minimumAge")));
	_ui.reverseSubductionCheck->setChecked((_earth->getModelParameter("reverseSubduction") != 0.0));
	_ui.upperPlatesAgeRatio->setText(QString::number(_earth->getModelParameter("upperPlates_ratioAgeValue")));
	
	_ui.tauSubNoise->setText(QString::number(_earth->getModelParameter("tauSub_randomNoise")));
	
	_ui.alwaysSubCheck->setChecked((_earth->getModelParameter("alwaysSubduct") != 0.0));
	
	_ui.SSCCheck->setChecked((_earth->getModelParameter("small_scale_convection") != 0.0));
	_ui.SSCConstantCheck->setChecked((_earth->getModelParameter("sscMode_constant") != 0.0));
	_ui.SSCConvectiveCheck->setChecked((_earth->getModelParameter("sscMode_convective") != 0.0));
	
	_ui.fixedConfigurationCheck->setChecked((_earth->getModelParameter("fixed_configuration") != 0.0));  
	
	_ui.U_BSE->setText(QString::number(_earth->getModelParameter("U_BSE")));
	_ui.Th_BSE->setText(QString::number(_earth->getModelParameter("Th_BSE")));
	_ui.K_BSE->setText(QString::number(_earth->getModelParameter("K_BSE")));
	_ui.U_cont->setText(QString::number(_earth->getModelParameter("U_cont")));
	_ui.Th_cont->setText(QString::number(_earth->getModelParameter("Th_cont")));
	_ui.K_cont->setText(QString::number(_earth->getModelParameter("K_cont")));
	_ui.depletedMantleCheck->setChecked((_earth->getModelParameter("depletion_always") != 0.0));
	_ui.primitiveMantleCheck->setChecked((_earth->getModelParameter("depletion_never") != 0.0));
	_ui.adaptativeDepletionCheck->setChecked((_earth->getModelParameter("depletion_adaptative") != 0.0));
	
	_ui.constantHCheck->setChecked((_earth->getModelParameter("constant_internal_heating") != 0.0));
	_ui.constantH->setText(QString::number(_earth->getModelParameter("constant_heat_value")));
	
	_ui.middleBreakupCheck->setChecked((_earth->getModelParameter("middle_breakup") != 0.0));
	_ui.randomBreakupCheck->setChecked((_earth->getModelParameter("middle_breakup") == 0.0));
	_ui.breakupPosWidth->setValue(_earth->getModelParameter("breakup_position"));
	_ui.contGrowthCheck->setChecked(_earth->getModelParameter("continental_growth") != 0.0);
	_ui.contGrowthCoeff->setText(QString::number(_earth->getModelParameter("contGrowthCoeff")));
	
	_ui.randomSeed->setValue(_earth->getModelParameter("randomSeed"));
}
// --------------------------------------------------
void
MacmaGui::putInterfacesInUI()
{
	// PUT INTERFACES IN UI
	vector<GeoElement*> & elements = _earth->accessElements();
	for(unsigned int i=0; i<elements.size(); i++)
		_ui.interfaceList->addItem(QString(elements[i]->getId().c_str()));
}
// --------------------------------------------------
void
MacmaGui::putPlatesInUI()
{
	// PUT PLATES IN UI
	vector<Plate*> & plates = _earth->accessPlates();
	for(unsigned int i=0; i<plates.size(); i++)
		_ui.plateList->addItem(QString(plates[i]->getId().c_str()));
}
// ==================================================
/*
 ============================================
 Earth
 ============================================
 */

void
MacmaGui::_earthThread()
{
	double t0, t1;
	double time0, time1;
	while(true)
    {
		if(_state == RUNNING && _earth->isRunning())
		{
			
			_earth->setRunStarted(true);
			
			time0 = _earth->getTimeMyr();
			time1 = _earth->getTimeMyr();	  
			while(time1 - time0 < _earth->getWriteImages())
			{
				_earth->compute();
				_xmlTimer->isTimeReached(_earth->getWriteXml());	  
				_imagesTimer->isTimeReached(_earth->getWriteImages());
				
				time1 = _earth->getTimeMyr();	  
			}
			struct timeval toWait;
			toWait.tv_sec = 0.0;
			toWait.tv_usec = _wait_usec;
			
			select(0,NULL,NULL,NULL,&toWait);
			
			if(_earth->final_condition_is_reached())
			{
				_startT = _earth->getTimeMyr();
				cout << "Final condition is reached: "
                << _earth->getEndAge() << "\t" << _earth->getT() << endl;
				toggleSuspended();
			}
		}
		
    }
}
// ==================================================
