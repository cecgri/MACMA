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


#ifndef MACMAGUI_H_
#define MACMAGUI_H_

// --- STD
#include <iostream>

// --- QT
#include <QMainWindow>
#include <QtGui>

// --- MACMA
#include "MACMA/graphics/macma.h"

#include "MACMA/graphics/graphics.h"
/* #include "MACMA/graphics/macmaGraphs.h" */
#include "MACMA/graphics/earthTimer.h"

// --- UI
#include "ui_MACMA.h"

class Gnuplot;

class MacmaGui : public QMainWindow, public Macma
{
  Q_OBJECT
    
 public:
  MacmaGui(QWidget* parent = 0);
  ~MacmaGui();
  
  virtual void open(QString filename);
  virtual void save(QString filename);
  virtual void record();
  
  // --- Utils
  virtual void toggleSuspended();
  virtual void unselectAll();
  virtual EarthState getEarthState();
  
  virtual double getTime();
  virtual void lock();
  virtual void unlock();
  
  void earthActivation(bool);
  void elementActivation(bool);
  void plateActivation(bool);
  void continentActivation(bool);
  
  // --- Workspace
  void askForWorkspace();
  void prepareWorkspace(string path);
  void cleanDirectory(string path);
  
  inline bool isA(string what) { return what == "MacmaGui"; }

  virtual void putParametersInUI();
  virtual void putInterfacesInUI();
  virtual void putPlatesInUI();

  private slots:

    void _slotOpen();
    void _slotSave();
    void _slotRun();
    /* void _slotRecord(); */
    /* void _slotGraphs(); */
    
    // EARTH ---------------------------
    void _slotEarthButton();
    // time:
    void _slotChangeStartAge(double);
    void _slotChangeEndAge(double);

    void _slotChangeWaitms(int);

    // MODEL ---------------------------
    void _slotModelButton();

    // interface ----------------------
    void _slotInterfaceList(int);
    void _slotAddInterface();
    void _slotRemoveInterface();
    void _slotInterfacePosition(double);
    void _slotSlabDepth(int);
    void _slotSaveInterfaceParameters();
    void _slotNextInterface();
  
    void _slotPlateList(int);
    void _slotSectionList(int);
    void _slotContinentList(int);
    void _slotAddContinent();
    void _slotRemoveContinent();
    void _slotContinentPosition(double);
    void _slotContinentLength(double);
    void _slotSaveContinentParameters();
  
    // WRITINGS (IMAGES AND FILES) ------
    void _slotUpdateGraphics();
    void _slotRecordTimeout();
    void _slotWriteXml();
    void _slotWriteImages();

 private:

    void _loadWidgets();
    void _connectSignals();

    void _earthThread();
  
    Ui_Macma        _ui;
  
    double _startT;
    int         _frameCounter;
    /* map<string, QLineEdit*> _physicalQLineEdit; */
    /* map<string, QLineEdit*> _timeQLineEdit; */
    /* map<string, QDoubleSpinBox*> _timeQDoubleSpinBox; */
    /* map<string, QCheckBox*> _timeQCheckBox; */
    /* map<string, QLineEdit*> _modelQLineEdit; */
    /* map<string, QRadioButton*> _modelQRadioButton; */
    /* map<string, QCheckBox*>  _modelQCheckBox; */

    string      _recordFolder;
  
    GeoElement*     _selectedElement;
    Plate*          _selectedPlate;
    PlateSection*   _selectedSection;
    Continent*      _selectedContinent;
  
    Graphics*       _graphics;
  
    QTimer*         _graphicsTimer;
    /* QTimer*         _recordTimer; */

    EarthTimer* _imagesTimer;
    EarthTimer* _xmlTimer;
    double _wait_usec;
  
    // --- Graphs
    /* bool _showGraphs; */
    /* MacmaGraphs* _graphs;     */

};

#endif


