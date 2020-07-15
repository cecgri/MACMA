// this file is part of MACMA
// author: cecile.grigne@univ-brest.fr

#ifndef ACTIVEMARGIN_H_
#define ACTIVEMARGIN_H_

#include <string>
#include "MACMA/physics/marker.h"

using std::string;

class Subduction;
class Continent;
class ContinentExtremity;

/*
 class ActiveMargin : extends marker
 Has possibly a subduction or a continentExtremity
 Can be currently active or not (just a tracer then)
 */

class ActiveMargin : public Marker
{
public:
    ActiveMargin(Earth* earth,
                 double position);
    
    ~ActiveMargin();
    
    bool isA(string);
    void updatePosition(double);
    
    inline string getClass() { return "ActiveMargin"; }
    
    inline void setSubduction(Subduction* subduction) { _subduction = subduction; }
    inline Subduction* getSubduction() { return _subduction; }

    inline void setContinentExtremity(ContinentExtremity* extremity) { _extremity = extremity; }
    inline ContinentExtremity* getContinentExtremity() { return _extremity; }
    
	inline string getStatus() { return _status; }
	inline void setStatus(string status) { _status = status; }
	
protected:
    Subduction* _subduction;
    ContinentExtremity* _extremity;
	string _status;   /* can be "active", "CO.XX" (shortId of the collision 
					   that ended this active margin) or "other"
					   */
};

#endif


