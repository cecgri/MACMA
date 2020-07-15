/* April 2015.
   This file is part of MACMA
   Class MacmaNonGui: Command line version of MACMA (author: leyaouanq@cervval.com)
  
   author: cecile.grigne@univ-brest.fr
*/

#ifndef MACMANONGUI_H_
#define MACMANONGUI_H_

// --- Parent class Macma
#include "MACMA/graphics/macma.h"

// --- Boost
#include <boost/filesystem.hpp>

class Earth;

// --------------------------------------------------
class MacmaNonGui : public Macma
{
 public:
  MacmaNonGui();
  ~MacmaNonGui();
  
  bool checkFolder(); 
  void run(string);
  
  inline bool isA(string what) { return what == "MacmaNonGui"; }

 protected:
    bool prepareWorkspace();

 private:
    void _earthThread();
};

// --------------------------------------------------
#endif

