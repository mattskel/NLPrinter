//
//  nlprettyprinter.h
//  libminizinc
//
//  Created by Matthew Skelley on 27/01/2016.
//
//

#ifndef __MINIZINC_NLPRETTYPRINTER_HH__
#define __MINIZINC_NLPRETTYPRINTER_HH__

#include <iostream>
#include <vector>

#include <minizinc/ast.hh>

namespace MiniZinc {
  
  class NLPrinter {
  private:
    std::ostream& _os;
    
  public:
    NLPrinter(std::ostream& os) : _os(os) {}
    
    std::vector<Expression*> outputVariables;
    std::vector<Expression*> outputArrays;

    void print(const Model* m);
    
    void getOutputVars(std::vector<int> v);
    
  };
}

#endif