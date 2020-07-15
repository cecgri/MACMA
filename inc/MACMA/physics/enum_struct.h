// this file is part of MACMA
// author: cecile.grigne@univ-brest.fr

#ifndef ENUM_STRUCT_H_
#define ENUM_STRUCT_H_

#include <vector>
#include <string>

/* ==========
   enums
   ========== */ 
enum Direction {
  LEFT,
  RIGHT,
  NONE
};
 
/* ==========
   structures 
   ========== */
struct Tridiag {
  double center, right, left, rhs;
};
 
struct AgeGroup {
  double ageMax;
  double rightPosition;
  double leftPosition;
};

struct Concentration {
  std::string id;
  double value; //value recomputed all the time
  double continent; // all the rest are constants
  double presentPrimitive;
  double presentDepleted;
  double lambda;
  double heat;
};

#endif

