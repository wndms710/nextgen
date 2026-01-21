#ifndef AVERAGEOPERATOR_HPP
#define AVERAGEOPERATOR_HPP

#include "DistributedDataExchanger.hpp"

#define FOR_IGR for (int igr = 0; igr < ngr; ++igr)
#define FOR_ACTIVE_IGR for (int igr = 0; igr < ngr_a; ++igr)

class AverageOperator : public DistributedDataExchanger {
  
public:
  int n;
  int * group;
  const double * weight;
  
  int ngr;
  int ngr_a;
  double * inv_weight_sum;
  
public:
  AverageOperator(int * my_ida_global,const int my_nda,const int daora[]);
  AverageOperator(const AverageOperator& copy);
  ~AverageOperator();
  
  int size() const { return n; }
  void max(double * phi);
  void min(double * phi);
  void apply(double * phi);
  void apply(double (*u)[3]);
};

#endif
