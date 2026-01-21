#ifndef TRANSPOSEOPERATOR_HPP
#define TRANSPOSEOPERATOR_HPP

#include "Common.hpp"
#include "MiscUtils.hpp"
using namespace MiscUtils;

#include "DistributedDataExchanger.hpp"

class TransposeOperator : public DistributedDataExchanger {

public:
  
  TransposeOperator(const int iora[],const int jora[]) {
    
    //cout << "iora[mpi_size]: " << iora[mpi_size] << endl;
    //cout << "jora[mpi_size]: " << jora[mpi_size] << endl;

    int8 * daora = new int8[mpi_size+1];
    for (int rank = 0; rank <= mpi_size; ++rank) 
      daora[rank] = int8(jora[rank])*int8(iora[mpi_size]);

    int my_nda = (iora[mpi_rank+1]-iora[mpi_rank])*jora[mpi_size];
    int8 * my_ida_global = new int8[my_nda];
    int ida = 0;
    for (int i_global = iora[mpi_rank]; i_global != iora[mpi_rank+1]; ++i_global) {
      const int i_local = i_global - iora[mpi_rank];
      for (int j_global = 0; j_global < jora[mpi_size]; ++j_global) {
	my_ida_global[i_local*jora[mpi_size] + j_global] = int8(j_global)*int8(iora[mpi_size]) + i_global;
	assert(ida == i_local*jora[mpi_size] + j_global);
	++ida;
      }
    }
    
    DistributedDataExchanger::init(my_ida_global,my_nda,daora);
    
    delete[] my_ida_global;
    delete[] daora;
        
  }

  void apply(double *pji[],const double * const pij[]) {

    // here we have been passed a 2d array pij where the i-index is distributed across the
    // processors according to iora and the full j index is present. We want to convert to 
    // a transposed representation where the j-index is distributed across the processors
    // according to jora and the full i-index is present.
    
    // note that we have to allow for zero-size array2d's here, so thus the
    // logic and NULL's...

    if ((pji == NULL)&&(pij == NULL))
      push((double*)NULL,(double*)NULL);
    else if (pij == NULL)
      push(pji[0],(double*)NULL);
    else if (pji == NULL)
      push((double*)NULL,pij[0]);
    else 
      push(pji[0],pij[0]);
  
  }

};

#endif

