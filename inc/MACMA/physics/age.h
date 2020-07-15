// this file is part of MACMA
// author: cecile.grigne@univ-brest.fr

#ifndef AGE_H_
#define AGE_H_

#include <string>
using std::string;

class Earth;

// =============================================
class Age
{
public:
	Age(Earth* earth, double val, double pos);
	~Age();
	
	// set/get
	inline void setPosition(double pos) { _position = pos; }
	inline void setValue(double val) { _value = val; }
	inline void setEffValue(double val) { _effValue = val; }
	inline void setOceanic(bool oc) { _ocean = oc; }
	inline void setExtra(bool extra) { _extra = extra; }
	inline void setU(double U) { _U = U; }
	inline void setOriginContCounter(unsigned int count) { _originContCounter = count; }
	
	inline double getPosition() { return _position; }
	inline double getValue() { return _value; }
	inline double getEffValue() { return _effValue; }
	inline bool isOceanic() { return _ocean; }
	inline bool isExtra() { return _extra; }
	inline double getU() { return _U; }
	inline unsigned int getOriginContCounter() { return _originContCounter; }
	
	inline bool isCritical(double val) { return _ocean && _value > val; }
	
	void reset(double);
	void addValue(double);
	void move(double);
	
protected:
	Earth* _earth;
	double _position;
	double _value;
	double _effValue;
	bool _ocean;
	bool _extra;
	double _U;
	unsigned int _originContCounter;
};

#endif
