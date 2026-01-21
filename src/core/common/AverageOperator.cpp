#include "AverageOperator.hpp"

AverageOperator::AverageOperator(int * my_ida_global,const int my_nda,const int daora[]) :
DistributedDataExchanger(my_ida_global,my_nda,daora) {
  
  group = NULL;
  weight = NULL;
  inv_weight_sum = NULL;
}

AverageOperator::AverageOperator(const AverageOperator& copy) : DistributedDataExchanger(copy) {
  // do not allow a copy constructor (unless you want to write one!)...
  throw(0);
}

AverageOperator::~AverageOperator() {
  
  DELETE(group);
  DELETE(inv_weight_sum);
}

void AverageOperator::max(double * phi) {
  
  double * phi_max = new double[ngr];
  FOR_IGR phi_max[igr] = -HUGE_VAL;
  
  for (int i = 0; i < n; ++i) {
    const int igr = group[i];
    phi_max[igr] = std::max(phi_max[igr],phi[i]);
  }
  
  push(phi_max,phi_max+ngr_a,MAX_DATA);
  pull(phi_max+ngr_a,phi_max);

  for (int i = 0; i < n; ++i) {
    const int igr = group[i];
    phi[i] = phi_max[igr];
  }
  
  delete[] phi_max;
}

void AverageOperator::min(double * phi) {
  
  double * phi_min = new double[ngr];
  FOR_IGR phi_min[igr] = HUGE_VAL;
  
  for (int i = 0; i < n; ++i) {
    const int igr = group[i];
    phi_min[igr] = std::min(phi_min[igr],phi[i]);
  }
  
  push(phi_min,phi_min+ngr_a,MIN_DATA);
  pull(phi_min+ngr_a,phi_min);
  
  for (int i = 0; i < n; ++i) {
    const int igr = group[i];
    phi[i] = phi_min[igr];
  }
  
  delete[] phi_min;
}

void AverageOperator::apply(double * phi) {
  
  assert(group != NULL);
  assert(weight != NULL);
  assert(inv_weight_sum != NULL);
  
  double * phi_avg = new double[ngr];
  FOR_IGR phi_avg[igr] = 0.0;

  for (int i = 0; i < n; ++i) {
    const int igr = group[i];
    phi_avg[igr] += weight[i]*phi[i];
  }
  
  push(phi_avg,phi_avg+ngr_a,ADD_DATA);
  
  FOR_ACTIVE_IGR {
    phi_avg[igr] *= inv_weight_sum[igr];
  }
  
  pull(phi_avg+ngr_a,phi_avg);
  
  for (int i = 0; i < n; ++i) {
    const int igr = group[i];
    phi[i] = phi_avg[igr];
  }
  
  delete[] phi_avg;
}

void AverageOperator::apply(double (*u)[3]) {
  
  assert(group != NULL);
  assert(weight != NULL);
  assert(inv_weight_sum != NULL);
  
  double (*u_avg)[3] = new double[ngr][3];
  FOR_IGR FOR_J3 u_avg[igr][j] = 0.0;

  for (int i = 0; i < n; ++i) {
    const int igr = group[i];
    FOR_J3 u_avg[igr][j] += weight[i]*u[i][j];
  }
  
  push(u_avg,u_avg+ngr_a,ADD_DATA);

  FOR_ACTIVE_IGR {
    FOR_J3 u_avg[igr][j] *= inv_weight_sum[igr];
  }
  
  pull(u_avg+ngr_a,u_avg);
  
  for (int i = 0; i < n; ++i) {
    const int igr = group[i];
    FOR_J3 u[i][j] = u_avg[igr][j];
  }
  
  delete[] u_avg;
}
