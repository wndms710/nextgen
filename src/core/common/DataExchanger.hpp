
using namespace CtiRegister;

class DataExchanger { 
public: 

  StaticSolver* solver;

  int nx_pack;
  int nx_unpack;

  int* send_count;
  int* recv_count;
  int* send_disp;
  int* recv_disp;

  double (*x_pack)[3];
  int *cvopt;
  int* idopt_recv;

  DataExchanger(StaticSolver* _solver) { 

    solver     = _solver;

    nx_pack     = 0;
    nx_unpack   = 0;

    send_count = NULL;
    recv_count = NULL;
    send_disp  = NULL;
    recv_disp  = NULL;

    x_pack     = NULL;
    cvopt      = NULL;
    idopt_recv = NULL;

  }


  ~DataExchanger() { 

    DELETE(send_count);
    DELETE(recv_count);
    DELETE(send_disp);
    DELETE(recv_disp);

    DELETE(x_pack);
    DELETE(cvopt);
    DELETE(idopt_recv);

  }


  void init(double (*x)[3], const int nx) { 

    if ( solver->cvAdt == NULL) 
      solver->buildCvAdt();

    send_count = new int[mpi_size];
    for (int rank = 0; rank < mpi_size; ++rank) 
      send_count[rank] = 0;

    vector<int> int_vec;

    double * send_buf = NULL;
    int8* send_buf_i8 = NULL;

    int send_count_sum;

    for (int iter = 0; iter < 2; ++iter) { 

      for (int ix = 0; ix < nx; ++ix) { 
        
        assert( int_vec.empty());
        solver->cvBboxAdt->buildListForPoint(int_vec,x[ix]);
        
        for (int ii =0, ii_end = int_vec.size(); ii < ii_end; ++ii) { 
          const int rank = int_vec[ii];  assert( (rank >=0) && (rank < mpi_size));

          if ( iter == 0) { 

            send_count[rank] += 3;

          } else { 

            send_buf[send_disp[rank]+0] = x[ix][0];
            send_buf[send_disp[rank]+1] = x[ix][1];
            send_buf[send_disp[rank]+2] = x[ix][2];

            send_buf_i8[send_disp[rank]/3] = BitUtils::packRankBitsIndex(mpi_rank,0,ix);
            send_disp[rank]            += 3;

          }
        }
        
        int_vec.clear();
      }

      if ( iter == 0) { 

        assert( send_disp == NULL); send_disp = new int[mpi_size];
        send_disp[0] = 0;
        for (int rank =1 ; rank < mpi_size; ++rank) 
          send_disp[rank] = send_count[rank-1] + send_disp[rank-1];

        send_count_sum = send_disp[mpi_size-1] + send_count[mpi_size-1];
        assert( send_buf == NULL); send_buf = new double[send_count_sum];
        assert( send_buf_i8 == NULL); send_buf_i8 = new int8[send_count_sum/3];

      }

    }

    // rewind the disp vec .. 

    send_disp[0] = 0;
    for (int rank = 1; rank < mpi_size; ++rank) 
      send_disp[rank] = send_count[rank-1] + send_disp[rank-1];

    assert( recv_count == NULL); recv_count = new int[mpi_size];
    MPI_Alltoall(send_count,1,MPI_INT,recv_count,1,MPI_INT,mpi_comm);

    assert( recv_disp == NULL); recv_disp = new int[mpi_size];
    recv_disp[0] = 0;
    for (int rank = 1; rank < mpi_size; ++rank) 
      recv_disp[rank] = recv_count[rank-1] + recv_disp[rank-1];

    const int recv_count_sum = recv_disp[mpi_size-1] + recv_count[mpi_size-1];

    double * recv_buf = new double[recv_count_sum];
    MPI_Alltoallv(send_buf,send_count,send_disp,MPI_DOUBLE, 
                  recv_buf,recv_count,recv_disp,MPI_DOUBLE,mpi_comm);


    int8* recv_buf_i8 = new int8[recv_count_sum/3];
    { 

      // send over the index as well.. 

      for (int rank = 0; rank < mpi_size; ++rank) { 

        send_count[rank] /= 3;
        send_disp[rank]  /= 3;

        recv_count[rank] /= 3;
        recv_disp[rank]  /= 3;

      }

      MPI_Alltoallv(send_buf_i8,send_count,send_disp,MPI_INT8, 
                    recv_buf_i8,recv_count,recv_disp,MPI_INT8, mpi_comm);


      for (int rank = 0; rank < mpi_size; ++rank) { 

        send_count[rank] *= 3;
        send_disp[rank]  *= 3;

        recv_count[rank] *= 3;
        recv_disp[rank]  *= 3;

      }

    }

    delete[] send_buf_i8;

    // now determine if the point is owned by this rank and which cv 
    // owns the value .. 

    int * x_flag = new int[recv_count_sum/3];
    for (int ix = 0; ix < recv_count_sum/3; ++ix) 
      x_flag[ix] = -1;

    nx_pack = 0; // number of x data members owned by this rank .. 

    for (int ir = 0; ir < recv_count_sum; ir += 3) { 

      double xp[3]; 
      for (int i = 0; i < 3; ++i) 
        xp[i] = recv_buf[ir+i];

      int_vec.clear();
      solver->cvAdt->buildListForPoint(int_vec,xp);

      const int ix = ir/3;

      double d2 = HUGE_VAL;
      for (int ii = 0, ii_end = int_vec.size(); ii < ii_end; ++ii) { 

        const int icv        = int_vec[ii]; 
        
        if ( icv < 0 ) { 
          cout << " icv invalid : " << icv << endl;
          assert(0);
        } else if ( icv >= solver->ncv) { 
          if ( icv < solver->ncv_g) { 
            cout << " icv is ghost? : " << icv << endl;
            assert(0);
          } else { 
            cout  << " icv is bigger than ncv : " << icv << endl;
            assert(0);
          }
        }
        
        //assert( (icv >=0) && (icv < solver->ncv));
        const double this_d2 = DIST2(xp,solver->x_vv[icv]);

        if ( (this_d2 <= d2) && (this_d2 <= (1.0+1.0e-10)*solver->r_vv[icv]*solver->r_vv[icv])) { 

          x_flag[ix] = icv;
          d2         = this_d2;

        }
      }

      if ( x_flag[ix] >= 0) { 

        // check against ghost data as well .. 

        const int icv = x_flag[ix]; assert( (icv >=0) && (icv < solver->ncv));

        if ( icv >= solver->ncv_i ) { 

          for (int coc = solver->cvocv_i[icv]; coc != solver->cvocv_i[icv+1]; ++coc) { 

            const int icv_nbr = solver->cvocv_v[coc];
            if ( icv_nbr >= solver->ncv) { 

              const double this_d2 = DIST2(xp,solver->x_vv[icv_nbr]);
              if ( this_d2 < (1.0-1.0E-10)*d2) { 

                //cout << " nbr closer: " << COUT_VEC(xp) << "   " << COUT_VEC(solver->x_vv[icv]) << "    " << COUT_VEC(solver->x_vv[icv_nbr]) << endl;
                x_flag[ix] = -1;
              
              } else if ( this_d2 < (1.0+1.0E-10)*d2) { 

                int rank,bits,index;
                BitUtils::unpackRankBitsIndex(rank,bits,index,solver->rbi_g[icv_nbr-solver->ncv]);

                if ( solver->icv_global[icv] > solver->icv_global[icv_nbr]) { 

                  //cout << " matched dist, nbr wins: " << COUT_VEC(xp) << "    " << COUT_VEC(solver->x_vv[icv]) << "    " << COUT_VEC(solver->x_vv[icv_nbr]) << endl;

                  if (bits == 0)
                    x_flag[ix] = -1;      // other proc owns it
                  else  
                    x_flag[ix] = icv_nbr; // we will get from our ranks ghost data
                  break;
                }
                
              }
            }
          }
        }
      }

      if ( x_flag[ix] >= 0) 
        ++nx_pack; 
    }

    
    { 
      int ntotal; 
      int nx_tmp = nx;
      MPI_Reduce(&nx_tmp,&ntotal,1,MPI_INT,MPI_SUM,0,mpi_comm);
      
      int npack_total;
      MPI_Reduce(&nx_pack,&npack_total,1,MPI_INT,MPI_SUM,0,mpi_comm);

      if ( mpi_rank == 0) { 

        cout << " > ntotal, npack_total : " << ntotal << "   " << npack_total << endl;
        cout.flush();

      }

    }
    
    //throw(0);

    DELETE(send_buf);

    // now build the x_pack vectors .. we're going to swap the 
    // meaning of the send and recv sides as well now... 

    for (int rank = 0; rank < mpi_size; ++rank) { 

      send_count[rank] = 0;
      send_disp[rank]  = 0;
    
    }

    assert( x_pack == NULL); x_pack = new double[nx_pack][3];
    assert( cvopt  == NULL); cvopt  = new int[nx_pack];
    int* idopt_send  = new int[nx_pack];

    for (int iter = 0; iter < 2; ++iter) { 

      for (int ir = 0; ir < recv_count_sum; ir += 3) { 

        const int ix = ir/3;
        if ( x_flag[ix] >= 0) {

          int rank,bits,index; 
          BitUtils::unpackRankBitsIndex(rank,bits,index,recv_buf_i8[ix]);
                          
          if ( iter == 0) { 
            
            send_count[rank]++;

          } else { 

            // note that the x_pack is reordered so that you can 
            // the send size is ready to go ... 

            for (int i = 0; i < 3; ++i) 
              x_pack[send_disp[rank]][i] = recv_buf[ir+i];
            cvopt[send_disp[rank]] = x_flag[ix];
            idopt_send[send_disp[rank]] = index;

            send_disp[rank]++;

          }

        }
      }

      // rebuild the send_disp on each iter ... 

      send_disp[0] = 0;
      for (int rank = 1; rank < mpi_size; ++rank) 
        send_disp[rank] = send_disp[rank-1] + send_count[rank-1];

      if ( iter == 1) { // check the packing 

        assert( nx_pack == send_disp[mpi_size-1] + send_count[mpi_size-1]);

      }

    }

    delete[] x_flag;
    delete[] recv_buf;
    delete[] recv_buf_i8;
    recv_buf = NULL;

    // now send the counts back to the recv side .. 

    MPI_Alltoall(send_count,1,MPI_INT,recv_count,1,MPI_INT,mpi_comm);

    // check to make sure that we are getting back exactly the number 
    // of requests that we sent out, ie that the requests are being 
    // singly serviced .. 

    recv_disp[0] = 0;
    for (int rank = 1; rank < mpi_size; ++rank) 
      recv_disp[rank] = recv_disp[rank-1] + recv_count[rank-1];

    nx_unpack = recv_count[mpi_size-1] + recv_disp[mpi_size-1];

    // we also need to know the order in which data is coming back. 

    assert( idopt_recv == NULL); idopt_recv = new int[nx_unpack];

    MPI_Alltoallv(idopt_send,send_count,send_disp,MPI_INT, 
                  idopt_recv,recv_count,recv_disp,MPI_INT,mpi_comm);


    delete[] idopt_send;

    
    // lastly let's check the x_transfer .. 

    if ( mpi_rank == 0 ) 
      cout << " > DataExchanger::about to check the coordinates ... " << endl;

    { 

      double *x_check = new double[3*nx_unpack];
      
      for (int rank = 0; rank < mpi_size; ++rank) { 

        send_count[rank] *= 3;
        send_disp[rank]  *= 3;

        recv_count[rank] *= 3;
        recv_disp[rank]  *= 3;

      }

      MPI_Alltoallv((double*)x_pack,send_count,send_disp,MPI_DOUBLE,
                    x_check,recv_count,recv_disp,MPI_DOUBLE, mpi_comm);

      for (int ii = 0; ii < nx_unpack; ++ii) { 

        const int ix = idopt_recv[ii];
        for(int i = 0; i < 3; ++i) 
          assert( abs(x_check[3*ii+i] - x[ix][i]) < 1.0e-14);

      }
    
      for (int rank = 0; rank < mpi_size; ++rank) { 

        send_count[rank] /= 3;
        send_disp[rank]  /= 3;

        recv_count[rank] /= 3;
        recv_disp[rank]  /= 3;

      }

      delete[] x_check;

    }

    if ( mpi_rank == 0) 
      cout << " > DataExchanger:: coordinate check passed. " << endl;

  }

  void checkMaxDist() { 

    double my_max_d2 = 0.0;

    for (int ix = 0; ix < nx_pack; ++ix) { 

      const int icv    = cvopt[ix];
      const double dx[3] = DIFF(x_pack[ix],solver->x_cv[icv]);
      my_max_d2          = max(my_max_d2, DOT_PRODUCT(dx,dx)/(solver->r_vv[icv]*solver->r_vv[icv]));

    }

    double max_d2;
    MPI_Reduce(&my_max_d2,&max_d2,1,MPI_DOUBLE,MPI_MAX,0,mpi_comm);
    
    if ( mpi_rank == 0) 
      cout << " > max dist : " << max_d2 << endl;

  }

  void exchangeCvDataWithGhosts(double (*r2)[3], const double (*u)[3], const double (*dudx)[3][3]) {
    
    // u and dudx MUST be updated prior to exchange

    double (*send_buf)[3] = new double[nx_pack][3];

    // for now we'll assume that we can grab a ptr to the data .. 
    // and that the ghosts are populated for this data array (level 1..) 
    // achtung, please be careful. 

    for (int ix = 0; ix < nx_pack; ++ix) { 
      
      const int icv = cvopt[ix]; assert((icv >= 0)&&(icv < solver->ncv_g));
      const double dx[3] = DIFF(x_pack[ix],solver->x_cv[icv]);
      for (int i = 0; i < 3; ++i) 
        send_buf[ix][i] = u[icv][i] + DOT_PRODUCT(dudx[icv][i],dx);

    }

    // exchange ...

    for (int rank = 0; rank < mpi_size; ++rank) { 

      send_count[rank] *= 3;
      send_disp[rank]  *= 3;

      recv_count[rank] *= 3;
      recv_disp[rank]  *= 3;

    }

    double *recv_buf = new double[3*nx_unpack];

    MPI_Alltoallv((double*)send_buf,send_count,send_disp,MPI_DOUBLE,
                  recv_buf,recv_count,recv_disp,MPI_DOUBLE,mpi_comm);


    delete[] send_buf;

    for (int ii = 0; ii < nx_unpack; ++ii) { 

      const int ix = idopt_recv[ii];
      for (int i = 0; i < 3; ++i) 
        r2[ix][i] = recv_buf[3*ii+i];

    }

    delete[] recv_buf;

    for (int rank = 0; rank < mpi_size; ++rank) { 

      send_count[rank] /= 3;
      send_disp[rank]  /= 3;

      recv_count[rank] /= 3;
      recv_disp[rank]  /= 3;

    }

  } 

  void exchangeCvDataWithGhosts(double *r1, const double *T, const double (*dTdx)[3]) {
    
    // T and dTdx MUST be updated prior to exchange

    double *send_buf = new double[nx_pack];

    // for now we'll assume that we can grab a ptr to the data .. 
    // and that the ghosts are populated for this data array (level 1..) 
    // achtung, please be careful. 

    for (int ix = 0; ix < nx_pack; ++ix) { 
      
      const int icv = cvopt[ix]; assert((icv >= 0)&&(icv < solver->ncv_g));
      const double dx[3] = DIFF(x_pack[ix],solver->x_cv[icv]);
      send_buf[ix] = T[icv] + DOT_PRODUCT(dTdx[icv],dx);

    }

    double *recv_buf = new double[nx_unpack];

    MPI_Alltoallv(send_buf,send_count,send_disp,MPI_DOUBLE,
                  recv_buf,recv_count,recv_disp,MPI_DOUBLE,mpi_comm);


    delete[] send_buf;

    for (int ii = 0; ii < nx_unpack; ++ii) { 

      const int ix = idopt_recv[ii];
      r1[ix] = recv_buf[ii];

    }

    delete[] recv_buf;
    
  } 

};
