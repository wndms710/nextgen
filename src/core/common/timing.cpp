
#include <iostream> 
#include <mpi.h> 
#include <cstdlib> 
#include <assert.h>
#include <fstream>
#include <math.h> 
#include <string.h>
#include <stdio.h>

using namespace std ; 

// 
// testing point to point latencies, 
// bandwidth, etc from the mpi 
// 
// 
class MpiTimer { 
public: 
  static const int timer_size = 7 ; 
  double buf[timer_size] ; 
  int flag ; 
  
  MpiTimer() { 
    flag = 0; 
  }
  
  void start() { 
    flag   = 0 ;
      buf[flag++] = MPI_Wtime() ; 
  } 
  
  void split() { 
    buf[flag++] = MPI_Wtime() ; 
  } 
  
  void stop() { 
    split() ; 
  } 
  
  // report the elapsed time from the start 
  // to split j
  
  double elapsed(const int j1, const int j0) const{ 
    return buf[j1] - buf[j0]; 
  }
  
  double elapsed(const int j) const { 
    return elapsed(j,0); 
  }
  
  double elapsed() const { 
    return elapsed(1,0) ; 
  }
  
}; 


class MpiProbe { 

public:  

  MPI_Comm mpi_comm_shared ; 
  int mpi_rank_shared ; 
  int mpi_size_shared ; 

  MPI_Comm intranode_comm ; 
  int mpi_node_rank ; 
  int mpi_node_size ; 

  MPI_Comm mpi_comm ; 
  int mpi_rank ; 
  int mpi_size ; 

  char name[MPI_MAX_PROCESSOR_NAME];

  void initSplitComms() { 



    MPI_Comm_dup(MPI_COMM_WORLD, &mpi_comm); 
    MPI_Comm_rank(mpi_comm,&mpi_rank); 
    MPI_Comm_size(mpi_comm,&mpi_size); 
    int len;
    MPI_Get_processor_name(name, &len);


    MPI_Comm_split_type(mpi_comm, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &mpi_comm_shared);
    if ( mpi_comm_shared == MPI_COMM_NULL) { 
      if ( mpi_rank == 0 ) 
        cerr << "> Failed to allocated mpi shared communicator "  << endl ;
      throw(0); 
    } 
    
    MPI_Comm_rank(mpi_comm_shared, &mpi_rank_shared) ; 
    MPI_Comm_size(mpi_comm_shared, &mpi_size_shared) ;

    MPI_Comm_split(mpi_comm, mpi_rank_shared, mpi_rank, &intranode_comm); 
    MPI_Comm_rank(intranode_comm, &mpi_node_rank ) ; 
    MPI_Comm_size(intranode_comm, &mpi_node_size ) ; 

    if ( mpi_rank == 0 ) 
      cout << " > finished init split comms " << endl; 
  } 


  void p2p() { 

    // 
    // there are potential a larger number of ranks, but 
    // assume that we can perform an exhaustive search on the nodes.
    // This is O(N^2) for the nodes.  first, identify a rank 
    // on each of the nodes. 
    //

    if ( mpi_rank == 0 ) 
      cout << " > p2p tests .. " << endl ; 

    double my_latencies[2] = {-1.0e+20, 0.0}; // my minimum and max latencies .. 
    double my_lat_stats[2] = {0.0,0.0};
    double my_bw[2]        = {1.0e+20,0.0}; 
    MpiTimer timer ; 

    const int BUF_SIZE = 1024 * 1024 / 4 ; // or int_size 
    int * my_buf = new int[BUF_SIZE]; 

    double *my_lat_node = NULL;
    if (mpi_rank_shared == 0) {
      my_lat_node = new double[mpi_node_size]; // latency of node
      for (int rank = 0 ; rank < mpi_node_size ; ++rank ) {
        my_lat_node[rank] = 0.0;
      }
    }

    for (int rank = 0 ; rank < mpi_node_size ; ++rank ) { 
      
      MPI_Barrier(mpi_comm) ; // make sure everyone reports in ... 
     
      if ( mpi_rank_shared == 0 ) {

        if ( rank == mpi_node_rank ) {
          cout << " Node " << mpi_node_rank << " : " << name << endl;
          cout.flush();
          for (int rank2 = 0 ; rank2 < mpi_node_size ; ++rank2 ) { 
          
            if ( rank2 == rank) continue ; 

            int my_int = rand() % 1000 ; 
            int my_int_check ; 
            MPI_Request req ;  
            MPI_Status status ; 
            
            timer.start() ; 
            MPI_Isend(&my_int,1,MPI_INT,rank2,12345,intranode_comm, &req); 
                
            // wait for it to come back ... 
            MPI_Recv(&my_int_check,1,MPI_INT,rank2,23456,intranode_comm,&status); 
            timer.stop() ; 


            double this_lat = timer.elapsed() ; 
            assert ( my_int == my_int_check) ; 

            my_latencies[0] = max(my_latencies[0], -timer.elapsed()) ; 
            my_latencies[1] = max(my_latencies[1],  timer.elapsed()) ; 

            my_lat_stats[0] += timer.elapsed() ; 
            my_lat_stats[1] += timer.elapsed() * timer.elapsed() ; 

            my_lat_node[rank] += this_lat;
            my_lat_node[rank2] += this_lat;

            // check a prototype 1Mb msg 
            for (int i =0 ;i < BUF_SIZE ; ++i) 
              my_buf[i] = rand()% 1000 ; 

            timer.start() ; 
            MPI_Isend(my_buf,BUF_SIZE,MPI_INT,rank2,12345,intranode_comm,&req); 
            MPI_Recv(my_buf,BUF_SIZE,MPI_INT,rank2,23456,intranode_comm,&status); 
            timer.stop() ; 

            double this_bw = 2.0/double(timer.elapsed()) ; // Mb/sec
            my_bw[0] = min(my_bw[0],this_bw); 
            my_bw[1] = max(my_bw[1],this_bw) ; 

            if ( rank2 > rank ) {
              cout << rank << " ---> " << rank2 << "  lat, bw " << this_lat << "   " << this_bw << endl ; 
              cout.flush();
            }

          }//rank2 ... 

          const int npeer = mpi_node_size-1 ; // skip yourself 
          my_lat_stats[0] /= double(npeer) ; 
          my_lat_stats[1] /= double(npeer) ; 
          my_lat_stats[1] -= my_lat_stats[0]*my_lat_stats[0] ;
          my_lat_stats[1] = sqrt(my_lat_stats[1]) ;

        }
        else { 

          int my_int ; 
          MPI_Request req ; 
          MPI_Status status ; 


          // get an int and send it right back... 
          MPI_Recv(&my_int,1,MPI_INT,rank,12345,intranode_comm,&status); 
          MPI_Send(&my_int,1,MPI_INT,rank,23456,intranode_comm) ;


          // get a 1Mb buffer and send it back .. 
          MPI_Recv(my_buf,BUF_SIZE,MPI_INT,rank,12345,intranode_comm,&status); 
          MPI_Send(my_buf,BUF_SIZE,MPI_INT,rank,23456,intranode_comm) ; 


        }
      }
    }//rank .. 

    // assemble the statistics on the latencies ...  
    double ex_latencies[2]; 
    double stat_lat[2] ; 

    MPI_Reduce(my_latencies,ex_latencies,2,MPI_DOUBLE,MPI_MAX,0,intranode_comm); 
    MPI_Reduce(my_lat_stats,stat_lat,2,MPI_DOUBLE,MPI_MAX,0,intranode_comm) ; 

    if ( (mpi_rank_shared == 0)&&(mpi_node_rank==0)) { 
      cout << " Min, max of p2p latencies (round-trip) [secs]: " << -ex_latencies[0] << "    " << ex_latencies[1] << endl ; 
      cout << " Worst case mean p2p (round-trip) [secs]  : " <<  stat_lat[0] << endl ; 
    } 

    // sum the rows and columns...
    double my_avg_lat,avg_lat;
    if (mpi_rank_shared == 0) {
      double *lat_node = NULL;
      if (mpi_rank == 0)
        lat_node = new double[mpi_node_size];
      MPI_Reduce(my_lat_node,lat_node,mpi_node_size,MPI_DOUBLE,MPI_SUM,0,intranode_comm) ; 
      delete[] my_lat_node;
      if (mpi_rank == 0) {
        avg_lat = 0.0;
        for (int rank = 0; rank < mpi_node_size; ++rank) {
          lat_node[rank] /= double(2*(mpi_node_size-1)); // row and column, skipping diagonal
          avg_lat += lat_node[rank];
        }
        avg_lat /= double(mpi_node_size);
      }
      MPI_Bcast(&avg_lat,1,MPI_DOUBLE,0,intranode_comm);
      MPI_Scatter(lat_node,1,MPI_DOUBLE,&my_avg_lat,1,MPI_DOUBLE,0,intranode_comm);
      if (mpi_rank == 0)
        delete[] lat_node;
    }
    for (int rank = 0; rank < mpi_node_size; ++rank) {
      if ((mpi_rank_shared == 0) && (mpi_node_rank == rank)) {
        if (my_avg_lat > 10.0*avg_lat) 
          cout << " Node "  << name << " p2p latency: " << my_avg_lat << " >> system p2p latency: " << avg_lat << endl; 
      }
      MPI_Barrier(mpi_comm);
    }

    // assemble statistics on the mean bandwidth .. 
    double min_bw, max_bw; 
    MPI_Reduce(&my_bw[0],&min_bw,1,MPI_DOUBLE,MPI_MIN,0,intranode_comm); 
    MPI_Reduce(&my_bw[1],&max_bw,1,MPI_DOUBLE,MPI_MAX,0,intranode_comm); 

    if ( (mpi_rank_shared ==0)&&(mpi_node_rank==0)) { 
      cout << " Min bandwidth [MB/sec] (1MB msg) = " << min_bw << endl ; 
      cout << " Max bandwidth [MB/sec] (1MB msg) = " << max_bw << endl ; 
    } 


    delete[] my_buf ; 


  }//p2p  


  void collective() { 

    // sizes of 1 double, 32 doubles, 1Mb, 10Mb
    const int sizes[4] = {8,256,1024,10*1024} ; 

    for (int iter = 0 ; iter < 2 ; ++iter ) { 
      bool isAllreduce = (iter != 0) ; 
      
      for (int i = 0; i < 4 ; ++i) {
        double time = reduce(sizes[i],isAllreduce); 
        double mb   = double(sizes[i])/double(1024*1024) ; 
        if ( (iter==0)&&(mpi_rank==0)) 
          cout << " > Reduce [mpi_size, time, Mb, Mb/s] = " << mpi_size << "    " << time << "    " << mb 
               << "    " << mb/time << endl ; 
        else if ( (iter==1)&&(mpi_rank==0) ) 
          cout << " > Allreduce [mpi_size, time, Mb, Mb/s] = " << mpi_size << "    " << time << "    " << mb 
               << "    " << mb/time << endl ; 
       
      }//i
    }//iter



    // alltoall communications of 1k,4k... 
    a2a(1024) ; 
    a2a(4096) ; 

  }//collective 


  double reduce(const int _size,bool isAllreduce=false) { 


    // 
    // dumpRanges, solvers and dumping information all 
    // extensively make use of reduction... 
    // NOTE size is in bytes ... 
    // 

    const int double_size = 8 ; 
    const int n           = _size/double_size ;

    double * my_buf       = new double[n]; 
    double * buf          = new double[n]; 

    MpiTimer timer ; 

    for (int i = 0 ; i < n ; ++i) 
      my_buf[i] = double(rand())/double(RAND_MAX) ; 

    timer.start() ; 
    if ( isAllreduce )  
      MPI_Allreduce(my_buf,buf,n,MPI_DOUBLE,MPI_SUM,mpi_comm); 
    else  
      MPI_Reduce(my_buf,buf,n,MPI_DOUBLE,MPI_SUM,0,mpi_comm); 
    timer.stop(); 


    delete[] my_buf; 
    delete[] buf ; 
    return timer.elapsed(); 
  } 


  void a2a(const int _size) { 

    // 
    // flood an all-to-all communication through the network 
    // of a given size ... 
    // 

    const int int_size = 4 ; 
    const int n        = _size/int_size ; 

    MpiTimer timer ; 

    //
    // although we're going to send this to everyone in a 
    // true O(p^2) operation (which could be done with a 
    // simple alltoall call), we're going to mimic the 
    // alltoall ,alltoallv pattern that appears in the code.
    // 

    int * send_count = new int[mpi_size] ; 
    int * recv_count = new int[mpi_size] ; 
    int * send_disp  = new int[mpi_size] ; 
    int * recv_disp  = new int[mpi_size] ; 
    int send_count_sum,recv_count_sum ; 

    for (int rank = 0 ; rank < mpi_size ; ++rank) 
      send_count[rank]  = n ; 

    timer.start() ; 
    MPI_Alltoall(send_count,1,MPI_INT,recv_count,1,MPI_INT,mpi_comm); 
    timer.split() ;

    send_disp[0] = 0 ; 
    recv_disp[0] = 0 ; 

    for (int i=1; i < mpi_size ; ++i) { 
      send_disp[i] = send_count[i-1] + send_disp[i-1]; 
      recv_disp[i] = recv_count[i-1] + recv_disp[i-1]; 
    } 

    send_count_sum = send_disp[mpi_size-1] + send_count[mpi_size-1]; 
    recv_count_sum = recv_disp[mpi_size-1] + recv_count[mpi_size-1]; 

    int * send_buffer = new int[send_count_sum] ; 
    int * recv_buffer = new int[recv_count_sum] ; 

    for (int i = 0 ; i < send_count_sum ; ++i) 
      send_buffer[i] = rand() ; 


    timer.split() ; 
    MPI_Alltoallv(send_buffer,send_count,send_disp,MPI_INT, 
                  recv_buffer,recv_count,recv_disp,MPI_INT,mpi_comm);  

    timer.split(); 


    const double t_a2a  = timer.elapsed(1,0); 
    const double t_a2av = timer.elapsed(3,2); 
    const double mb     = double(_size)/double(1024*1024); 

    if ( mpi_rank ==0 ) 
      cout << " mb, alltoall (single int) time, a2v time, bw [MB/s] = " 
           << mb << "   " << t_a2a << "    " << t_a2av << "    " << t_a2av/mb << endl ; 
  
    delete[] send_count ; 
    delete[] recv_count ; 
    delete[] send_disp ; 
    delete[] recv_disp ; 
    delete[] send_buffer ; 
    delete[] recv_buffer ; 
  } 

  
}; 


class CommBurst : public MpiProbe { 

public:
  int checkAndRemoveKillFile(const string& filename) {

    int found =-1; 
    if ( mpi_rank == 0 ) { 

      ifstream ifile ; 
      ifile.open(filename.c_str()) ; 
      if ( ifile.is_open())  { 
        found = 1;
        cout << " > Found killburst... going down." << endl; 
        ifile.close(); 
        remove(filename.c_str()) ; 
      }
      else 
        found = 0; 
    }

    MPI_Bcast(&found,1,MPI_INT,0,mpi_comm); 
    assert( (found ==0)||(found==1)) ; 
    return found ; 
  } 
      



  void burst() { 


    if ( mpi_rank == 0 )  
      cout << " > entering a2a burst mode... " << endl ;

    const int size = 4096 ; // 4k msgs from everyone .. 
    const int nmsg = 2048 ; // 2048 repeated msgs ..


    for (int n = 0 ; n < nmsg ; ++n) { 

      a2a(size); 

      // 
      // use killburst to check if we're done.. 
      // 
      if ( checkAndRemoveKillFile("killburst") == 1) 
        break; 
    }
  }


}; 


int main(int argc, char * argv[]) {
  
  try {
  
    // initialize the environment: basically 
    // initializes MPI and parameters...
   
    MPI_Init(&argc,&argv); 


    int mode = 1;
    if ( argc > 1 ) { 
      if ( strcmp(argv[1], "--p2p") == 0 )
        mode = 0; 
      else if ( strcmp(argv[1], "--alltoall") == 0 ) 
        mode = 1; 
    } 
    
    if ( mode == 0 ) {
      
      MpiProbe probe ;
      probe.initSplitComms() ; 
      probe.p2p() ; 
      probe.collective() ; 

    }
    else if ( mode == 1) { 
 
      // run in the comm burst mode.  this should 
      // flood the network with some communication to 
      // see if this get any worse.
      CommBurst cb ; 
      cb.initSplitComms();  
      cb.burst() ; 
    } 
    


    // finalize the environment: reports parameter usage, shuts down MPI...
    
    MPI_Finalize();
    
    // rudimentary error management...

  }
  catch(...) {
    MPI_Finalize();
  }
    
  return(0);

}


