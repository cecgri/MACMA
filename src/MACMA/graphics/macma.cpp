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

#include "MACMA/graphics/macma.h"
#include "MACMA/graphics/macmaGui.h"

using namespace std;


// --- Global
string _workspace_ = "";


// ==================================================

Macma::Macma()
{
  setMACMAState(EMPTY);

  _nextTimeWriteXml = 0.0;
  _nextTimeWriteDiagnosis = 0.0;
}

Macma::~Macma()
{
}

// ==================================================
void
Macma::openXmlFile(string filename)
{
	setlocale(LC_NUMERIC, "C");
	
	FILE* fp = fopen(filename.c_str(), "r");
	if(!fp) 
    {
		cerr << "Error while loading parameters from " + filename << endl;
		return;
    }
	
	string buffer;
	char line[255];
	while(fgets(line, 255, fp) != NULL) 
		buffer.append(line);
	fclose(fp);
	
	xml_document<> doc;
	doc.parse<0>((char*)buffer.c_str());
	
	_earth->lock();
	
	// TIME ---------------- 
	for(xml_node<>* node = doc.first_node("MACMA")->first_node("Time")->first_node("Parameter"); 
		node;
		node = node->next_sibling()
		)
    {
		char* name = node->first_attribute("name")->value();
		char* svalue = node->first_attribute("value")->value();
		
		double value = 0.0;
		sscanf(svalue, "%lg", &value);
		_earth->setTimeParameter(string(name), value);
    }
	_earth->formerTimeSetting();
	
	// PARAMETERS -----------------
	for(xml_node<>* node = doc.first_node("MACMA")->first_node("Parameters")->first_node("Parameter"); 
		node;
		node = node->next_sibling()
		)
    {
		char* name = node->first_attribute("name")->value();
		char* svalue = node->first_attribute("value")->value();
		
		double value = 0.0;
		sscanf(svalue, "%lg", &value);
		_earth->setPhysicalParameter(string(name), value);
    }
	
	// MODEL -----------------
	for(xml_node<>* node = doc.first_node("MACMA")->first_node("Model")->first_node("Parameter"); 
		node;
		node = node->next_sibling())
    {
		char* name = node->first_attribute("name")->value();
		char* svalue = node->first_attribute("value")->value();
		
		double value = 0.0;
		sscanf(svalue, "%lg", &value);
		_earth->setModelParameter(string(name), value);
    }
	
	_earth->unlock();
	
	/* convert double to bool or structure and 
     put correct dimensions */
	_earth->initParameters();
	
	/* if GUI mode: put read parameters in UI */
	if(this->isA("MacmaGui"))
    {
		((MacmaGui*)this)->putParametersInUI();
    }
	setMACMAState(PARAMETERS);
	
	// INTERFACES --------------------
	for(xml_node<>* node = 
		doc.first_node("MACMA")->first_node("Interfaces")->first_node("Interface"); 
		node;
		node = node->next_sibling()
		)
    {
		char* className = node->first_attribute("className")->value();
		char* sposition = node->first_attribute("position")->value();
		char* sleftAge = node->first_attribute("leftAge")->value();
		char* srightAge = node->first_attribute("rightAge")->value();
		char* sslabDepth = node->first_attribute("slabDepth")->value();
		
		double position = 0.0;
		sscanf(sposition, "%lg", &position);
		double leftAge = 0.0;
		sscanf(sleftAge, "%lg", &leftAge);
		double rightAge = 0.0;
		sscanf(srightAge, "%lg", &rightAge);
		double slabDepth = 0.0;
		sscanf(sslabDepth, "%lg", &slabDepth);
		
		GeoElement* element = NULL;
		
		if(!strcmp(className, "Ridge"))
			element = _earth->createRidge(position);
		
		else if(!strcmp(className, "LeftSubduction"))
			element = _earth->createSubduction(position, LEFT);
		
		else if(!strcmp(className, "RightSubduction"))
			element = _earth->createSubduction(position, RIGHT);
		
		else if(!strcmp(className, "Staple"))
			element = _earth->createStaple(position);
		
		if(element)
		{ /* leftAge and rightAge are double here.
		   Make class Age for elements.
		   */
			element->setRightAge(_earth, rightAge, element->getPosition());
			element->setLeftAge(_earth, leftAge, element->getPosition());
			
			if(element->isA("Subduction"))
			{
				double depth = slabDepth * 1.0E3;
				((Subduction*)element)->setDepth(depth);
			}
		}
    }
	
	_earth->sortElements();
	
	if(this->isA("MacmaGui"))
    {
		((MacmaGui*)this)->putInterfacesInUI();
    }
	setMACMAState(INTERFACES);
	
	
	// CONTINENTS ---------------------------
	for(xml_node<>* node = doc.first_node("MACMA")->first_node("Continents")->first_node("Continent"); 
		node;
		node = node->next_sibling())
    {
		char* id = node->first_attribute("id")->value();
		char* sposition = node->first_attribute("position")->value();
		char* slength = node->first_attribute("length")->value();
		char* sleftAge = node->first_attribute("leftAge")->value();
		char* srightAge = node->first_attribute("rightAge")->value();
		
		double position = 0.0;
		sscanf(sposition, "%lg", &position);
		double length = 0.0;
		sscanf(slength, "%lg", &length);
		double leftAge = 0.0;
		sscanf(sleftAge, "%lg", &leftAge);
		double rightAge = 0.0;
		sscanf(srightAge, "%lg", &rightAge);
		
		string strInit = "initial";
		Continent* continent = _earth->createContinent(position, length, strInit);
		_earth->createRightContinentExtremity(continent);
		_earth->createLeftContinentExtremity(continent);
		
		/* rightExtremity */		
		double rightPos = continent->getRightExtremity()->getPosition();
		continent->getRightExtremity()->setLeftAge(_earth, Earth::ageEarth, rightPos);
		continent->getRightExtremity()->setRightAge(_earth, rightAge, rightPos);
		continent->getRightExtremity()->getLeftAge().setOceanic(false);
		
		/* leftExtremity */
		double leftPos = continent->getLeftExtremity()->getPosition();
		continent->getLeftExtremity()->setRightAge(_earth, Earth::ageEarth, leftPos);
		continent->getLeftExtremity()->setLeftAge(_earth, leftAge, leftPos);
		continent->getLeftExtremity()->getRightAge().setOceanic(false);
    }
	
	_earth->clearPlates();
	_earth->checkNeighbors();
	_earth->fillStructure(false);
	
	if(this->isA("MacmaGui"))
    {
		((MacmaGui*)this)->putPlatesInUI();
    }
	
	setMACMAState(PLATES);
	
}
// ==================================================
void
Macma::writeXmlFile(string filename)
{
  // changing time-dependent parameters
  _earth->changeTimeParameters();
  _earth->changePhysicalParameters();

  saveXmlFile(filename);
}
// ==================================================
void
Macma::saveXmlFile(string filename)
{
  map<string, double> timeParameters = _earth->accessTimeParameters();
  map<string, double> physicalParameters = _earth->accessPhysicalParameters();
  map<string, double> modelParameters = _earth->accessModelParameters();

  // create xml doc
  xml_document<> doc;

  xml_node<>* root = doc.allocate_node(node_element, "MACMA");
  doc.append_node(root);
        

  // TIME --------------------------------
  xml_node<> *nodeTime = doc.allocate_node(node_element, "Time");
  root->append_node(nodeTime);

  map<string,double>::iterator it;
  for(it=timeParameters.begin();it!=timeParameters.end();it++)
    {
      xml_node<> *node = doc.allocate_node(node_element, "Parameter");
      nodeTime->append_node(node);

      xml_attribute<> *name = doc.allocate_attribute("name", it->first.c_str());
      node->append_attribute(name);

      ostringstream oss;
      oss << it->second << flush;
      char *node_value = doc.allocate_string(oss.str().c_str());
      xml_attribute<> *value = doc.allocate_attribute("value", node_value);
      node->append_attribute(value);
    }

  // PHYSICAL PARAMETERS --------------------
  xml_node<> *nodePhysical = doc.allocate_node(node_element, "Parameters");
  root->append_node(nodePhysical);

  for(it=physicalParameters.begin();it!=physicalParameters.end();it++)
    {
      xml_node<> *node = doc.allocate_node(node_element, "Parameter");
      nodePhysical->append_node(node);

      xml_attribute<> *name = doc.allocate_attribute("name", it->first.c_str());
      node->append_attribute(name);

      ostringstream oss;
      oss << it->second << flush;
      char *node_value = doc.allocate_string(oss.str().c_str());
      xml_attribute<> *value = doc.allocate_attribute("value", node_value);
      node->append_attribute(value);
    }

  // MODEL ---------------------------------
  xml_node<> *nodeModel = doc.allocate_node(node_element, "Model");
  root->append_node(nodeModel);

  for(it=modelParameters.begin();it!=modelParameters.end();it++)
    {
      xml_node<> *node = doc.allocate_node(node_element, "Parameter");
      nodeModel->append_node(node);

      xml_attribute<> *name = doc.allocate_attribute("name", it->first.c_str());
      node->append_attribute(name);

      ostringstream oss;
      oss << it->second << flush;
      char *node_value = doc.allocate_string(oss.str().c_str());
      xml_attribute<> *value = doc.allocate_attribute("value", node_value);
      node->append_attribute(value);
    }
  
  // INTERFACES ----------------------------
  xml_node<>* nodeInterfaces = doc.allocate_node(node_element, "Interfaces");
  root->append_node(nodeInterfaces);

  vector<GeoElement*> elements = _earth->accessElements();
  for(unsigned int i=0; i<elements.size(); i++)
    {
      if(elements[i]->isA("Interface") || elements[i]->isA("Staple"))
	{
	  xml_node<>* node = doc.allocate_node(node_element, "Interface");
	  nodeInterfaces->append_node(node);
	  
	  // className
	  char* node_className = doc.allocate_string(elements[i]->getClass().c_str());
	  xml_attribute<>* className = doc.allocate_attribute("className", node_className);
	  node->append_attribute(className);
	  
	  ostringstream oss;
	  // position
	  oss << elements[i]->getPosition() << flush;
	  char* node_position = doc.allocate_string(oss.str().c_str());
	  xml_attribute<>* position = doc.allocate_attribute("position", node_position);
	  node->append_attribute(position);
	  oss.str("");
	  
	  // left age
	  oss << elements[i]->getLeftAge().getValue() << flush;
	  char* node_leftAge = doc.allocate_string(oss.str().c_str());
	  xml_attribute<>* leftAge = doc.allocate_attribute("leftAge", node_leftAge);
	  node->append_attribute(leftAge);
	  oss.str("");
	  
	  // right age
	  oss << elements[i]->getRightAge().getValue() << flush;
	  char* node_rightAge = doc.allocate_string(oss.str().c_str());
	  xml_attribute<>* rightAge = doc.allocate_attribute("rightAge", node_rightAge);
	  node->append_attribute(rightAge);
	  oss.str("");
	  
	  // slab depth
	  double slabDepthValue = 0.0;
	  if(elements[i]->isA("Subduction"))
	    slabDepthValue = ((Subduction*)elements[i])->getDepth() / 1.0e3;

	  oss << slabDepthValue;
	  char *node_slabDepth = doc.allocate_string(oss.str().c_str());
	  xml_attribute<> *slabDepth = doc.allocate_attribute("slabDepth", node_slabDepth);
	  node->append_attribute(slabDepth);
	  oss.str("");
	    
	} // elements[i] isA interface or staple
    } // elements[i]
  
  // Continents - - - - - - - - - - 
  xml_node<>* nodeContinents = doc.allocate_node(node_element, "Continents");
  root->append_node(nodeContinents);
  
  vector<Continent*> continents = _earth->accessContinents();
  for(unsigned int i=0; i<continents.size(); i++)
    {
      xml_node<> *node = doc.allocate_node(node_element, "Continent");
      nodeContinents->append_node(node);
      
      // id
      char* node_continent = doc.allocate_string(continents[i]->getId().c_str());
      xml_attribute<> *id = doc.allocate_attribute("id", node_continent);
      node->append_attribute(id);
      
      ostringstream oss;
      // position
      oss << continents[i]->getPosition() << flush;
      char *node_position = doc.allocate_string(oss.str().c_str());
      xml_attribute<> *position = doc.allocate_attribute("position", node_position);
      node->append_attribute(position);
      oss.str("");
      
      // length
      oss << continents[i]->getLength() << flush;
      char *node_length = doc.allocate_string(oss.str().c_str());
      xml_attribute<> *length = doc.allocate_attribute("length", node_length);
      node->append_attribute(length);
      oss.str("");
      
      // leftAge
      oss << continents[i]->getLeftExtremity()->getLeftAge().getValue() << flush;
      char *node_leftAge = doc.allocate_string(oss.str().c_str());
      xml_attribute<> *leftAge = doc.allocate_attribute("leftAge", node_leftAge);
      node->append_attribute(leftAge);
      oss.str("");
      
      // rightAge
      oss << continents[i]->getRightExtremity()->getRightAge().getValue() << flush;
      char *node_rightAge = doc.allocate_string(oss.str().c_str());
      xml_attribute<> *rightAge = doc.allocate_attribute("rightAge", node_rightAge);
      node->append_attribute(rightAge);
      oss.str("");
    } // continents[i]
  
  // writing into file
  ofstream xmlFile(filename.c_str());
  if(xmlFile.is_open())
    {
      xmlFile << doc;
      xmlFile.close();
    }
  else
    cerr << "ERROR while saving " << filename << endl;
}
// ==================================================
void
Macma::setMACMAState(MACMAState state)
{
  _state = state; 
  //  cout << state << " - " << showMACMAState() << endl;
}
// ==================================================
string
Macma::showMACMAState()
{
  string sstate;

  if(_state == EMPTY)
    sstate = "NOT INITIALIZED";
  else if(_state == PARAMETERS)
    sstate = "PARAMETERS ARE SET";
  else if(_state == INTERFACES)
    sstate = "CHANGING INTERFACES";
  else if(_state == PLATES)
    sstate = "PLATES ARE SET";
  else if(_state == READY)
    sstate = "ALL IS READY";
  else if(_state == RUNNING)
    sstate = "...RUNNING";
  else
    cerr << "MACMAState not recognized" << endl;

  return sstate;
}
// ==================================================
