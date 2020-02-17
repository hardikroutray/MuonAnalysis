// -*- C++ -*-
//Add includes for your classes here
#include <vector>
#include "DataFormats/Common/interface/Wrapper.h"

namespace {
   struct dictionary {
      std::vector<std::vector<float> > vf2d;
      edm::Wrapper<std::vector<std::vector<float> > > wvf2d;

      std::vector<std::vector<bool> > vb2d;
      edm::Wrapper<std::vector<std::vector<bool> > > wvb2d;

      std::vector<std::vector<int> > vi2d;
      edm::Wrapper<std::vector<std::vector<int> > > wvi2d;

      
      
   };
}
