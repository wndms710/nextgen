#include "SimpleFunc.hpp"

SimpleFunc * processSimpleFunc(Param * param) {
  
  if (param == NULL) {
    CERR("processSimpleFunc passed NULL param");
  }
  
  SimpleFunc * ptr = NULL;

  double x0_check,dx_check,x1_check;
  bool b_check = false;
  
  int iarg = 0;
  while (iarg < param->size()) {
    string token = param->getString(iarg++);
    if (token == "CONSTANT") {
      assert(ptr == NULL);
      const double value = param->getDouble(iarg++);
      ptr = new ConstantFunc(value);
    }
    else if (token == "TABLE") {
      assert(ptr == NULL);
      vector<pair<double,double> > xyVec;
      while (iarg < param->size()) {
        if (param->getString(iarg) == "CHECK")
          break;
        const double x = param->getDouble(iarg++);
        const double y = param->getDouble(iarg++);
        xyVec.push_back(pair<double,double>(x,y));
      }
      ptr = new TableFunc(xyVec);
    }
    else if (token == "CLIPPED_CUBIC") {
      assert(ptr == NULL);
      // expect:
      // COEFF c0 c1 c2 c3
      // RANGE <min> <max>
      bool b_range = false;
      double range[2];
      bool b_coeff = false;
      double coeff[4];
      while (iarg < param->size()) {
        token = param->getString(iarg++);
        if (token == "COEFF") {
          coeff[0] = param->getDouble(iarg++);
          coeff[1] = param->getDouble(iarg++);
          coeff[2] = param->getDouble(iarg++);
          coeff[3] = param->getDouble(iarg++);
          b_coeff = true;
        }
        else if (token == "RANGE") {
          range[0] = param->getDouble(iarg++);
          range[1] = param->getDouble(iarg++);
          b_range = true;
        }
        else if (token == "CHECK") {
          if (iarg+3 > param->size()) {
            CERR("CHECK for SimpleFunc " << param->getName() << " requires <x0> <dx> <x1>");
          }
          x0_check = param->getDouble(iarg++); 
          dx_check = param->getDouble(iarg++); 
          x1_check = param->getDouble(iarg++);
          b_check = true;
        }
        else {
          CERR("unrecognized CLIPPED_CUBIC token: " << token);
        }
      }
      if (!b_coeff && !b_range) {
        CERR("CLIPPED_CUBIC missing COEFF <c0> <c1> <c2> <c3> and RANGE <xmin> <xmax>");
      }
      if (!b_coeff) {
        CERR("CLIPPED_CUBIC missing COEFF <c0> <c1> <c2> <c3>");
      }
      if (!b_range) {
        CERR("CLIPPED_CUBIC missing RANGE <xmin> <xmax>");
      }
      ptr = new ClippedCubicFunc(coeff,range);
    }
    else if (token == "CLIPPED_QUARTIC") {
      assert(ptr == NULL);
      // expect:
      // COEFF c0 c1 c2 c3 c4
      // RANGE <min> <max>
      bool b_range = false;
      double range[2];
      bool b_coeff = false;
      double coeff[5];
      while (iarg < param->size()) {
        token = param->getString(iarg++);
        if (token == "COEFF") {
          coeff[0] = param->getDouble(iarg++);
          coeff[1] = param->getDouble(iarg++);
          coeff[2] = param->getDouble(iarg++);
          coeff[3] = param->getDouble(iarg++);
          coeff[4] = param->getDouble(iarg++);
          b_coeff = true;
        }
        else if (token == "RANGE") {
          range[0] = param->getDouble(iarg++);
          range[1] = param->getDouble(iarg++);
          b_range = true;
        }
        else if (token == "CHECK") {
          if (iarg+3 > param->size()) {
            CERR("CHECK for SimpleFunc " << param->getName() << " requires <x0> <dx> <x1>");
          }
          x0_check = param->getDouble(iarg++); 
          dx_check = param->getDouble(iarg++); 
          x1_check = param->getDouble(iarg++);
          b_check = true;
        }
        else {
          CERR("unrecognized CLIPPED_QUARTIC token: " << token);
        }
      }
      if (!b_coeff && !b_range) {
        CERR("CLIPPED_QUARTIC missing COEFF <c0> <c1> <c2> <c3> and RANGE <xmin> <xmax>");
      }
      if (!b_coeff) {
        CERR("CLIPPED_QUARTIC missing COEFF <c0> <c1> <c2> <c3>");
      }
      if (!b_range) {
        CERR("CLIPPED_QUARTIC missing RANGE <xmin> <xmax>");
      }
      ptr = new ClippedQuarticFunc(coeff,range);
    }
    else if (token == "CHECK") {
      if (iarg+3 > param->size()) {
        CERR("CHECK for SimpleFunc " << param->getName() << " requires <x0> <dx> <x1>");
      }
      x0_check = param->getDouble(iarg++); 
      dx_check = param->getDouble(iarg++); 
      x1_check = param->getDouble(iarg++);
      b_check = true;
    }
    else {
      CERR(param->getName() << " has unrecognized SimpleFunc param: " << token);
    }
  }
  
  assert(ptr);

  // process CHECK if requested. This uses the eval machinery to 
  // allow the user to confirm their SimpleFunc's are as expected...
  
  if (b_check) {
    if (mpi_rank == 0) {
      cout << "CHECK for " << param->getName() << endl;
      for (double x = x0_check; x <= x1_check+1.0E-6*dx_check; x += dx_check) {
        const double y = ptr->eval(x);
        cout << "   " << x << " " << y << endl;
      }
    }
  }

  return ptr;
  
}

