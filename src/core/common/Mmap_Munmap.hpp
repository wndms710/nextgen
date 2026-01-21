// in the including code, be sure to define T appropriately
int CTI_Mmap(T* &ptr,const size_t size) {

  assert(ptr == NULL);

#ifdef WITH_SHM

  if (mpi_size_shared == 1) {

    ptr = new T[size];

  }
  else {

    assert(mpi_size_shared > 1);

    int fd_shm = -1;
    // the first rank in the shared comm should create the shared memory...
    if (mpi_rank_shared == 0) {
      fd_shm = shm_open(SHM_NAME, O_CREAT | O_EXCL | O_RDWR, 0600);
      if (fd_shm == -1) {
	cerr << " 1 mpi_rank " << mpi_rank << " had shm_open error: " << strerror(errno) << endl;
	throw(-1);
      }
      // set file size..
      if (ftruncate(fd_shm, sizeof(T)*size) == -1) {
	cerr << " 2 mpi_rank " << mpi_rank << " had ftrucate error: " << strerror(errno) << endl;
	throw(-1);
      }
    }

    MPI_Barrier(mpi_comm_shared);

    if (mpi_rank_shared != 0) {
      fd_shm = shm_open(SHM_NAME, O_RDWR, 0600);
      if (fd_shm == -1) {
	cerr << " 3 mpi_rank " << mpi_rank << " had shm_open error: " << strerror(errno) << endl;
	throw(-1);
      }
    }

    if (mpi_rank_shared == 0) {
      // only the first rank in the shared comm is able to write...
      ptr = (T *)mmap(NULL, sizeof(T)*size, PROT_READ | PROT_WRITE, MAP_SHARED, fd_shm, 0);
    }
    else {
      ptr = (T *)mmap(NULL, sizeof(T)*size, PROT_READ, MAP_SHARED, fd_shm, 0);
    }
    assert(ptr != NULL);

    MPI_Barrier(mpi_comm_shared);

    // don't need dev/shm file anymore, unlink removes it...
    if (mpi_rank_shared == 0) {
      if (shm_unlink(SHM_NAME) == -1) {
	cerr << " 4 mpi_rank " << mpi_rank << " had shm_unlink error: " << strerror(errno) << endl;
	throw(-1);
      }
    }

    // also don't need file descriptor..
    if (close(fd_shm) == -1) {
      cerr << " 5 mpi_rank " << mpi_rank << " had close error: " << strerror(errno) << endl;
      throw(-1);
    }

  }

#else

  ptr = new T[size];

#endif

  return(0);

}

int CTI_Mmap(T (*&ptr)[2],const size_t size) {

  assert(ptr == NULL);

#ifdef WITH_SHM

  if (mpi_size_shared == 1) {

    ptr = new T[size][2];

  }
  else {

    assert(mpi_size_shared > 1);

    int fd_shm = -1;
    // the first rank in the shared comm should create the shared memory...
    if (mpi_rank_shared == 0) {
      fd_shm = shm_open(SHM_NAME, O_CREAT | O_EXCL | O_RDWR, 0600);
      if (fd_shm == -1) {
	cerr << " 6 mpi_rank " << mpi_rank << " had shm_open error" << endl;
	throw(-1);
      }
      // set file size..
      if (ftruncate(fd_shm, sizeof(T)*2*size) == -1) {
	cerr << " 7 mpi_rank " << mpi_rank << " had ftrucate error" << endl;
	throw(-1);
      }
    }

    MPI_Barrier(mpi_comm_shared);

    if (mpi_rank_shared != 0) {
      fd_shm = shm_open(SHM_NAME, O_RDWR, 0600);
      if (fd_shm == -1) {
	cerr << " 8 mpi_rank " << mpi_rank << " had shm_open error" << endl;
	throw(-1);
      }
    }

    if (mpi_rank_shared == 0) {
      // only the first rank in the shared comm is able to write...
      ptr = (T (*)[2])mmap(NULL, sizeof(T)*2*size, PROT_READ | PROT_WRITE, MAP_SHARED, fd_shm, 0);
    }
    else {
      ptr = (T (*)[2])mmap(NULL, sizeof(T)*2*size, PROT_READ, MAP_SHARED, fd_shm, 0);
    }
    assert(ptr != NULL);

    MPI_Barrier(mpi_comm_shared);

    // don't need dev/shm file anymore, unlink removes it...
    if (mpi_rank_shared == 0) {
      if (shm_unlink(SHM_NAME) == -1) {
	cerr << " 9 mpi_rank " << mpi_rank << " had shm_unlink error" << endl;
	throw(-1);
      }
    }

    // also don't need file descriptor..
    if (close(fd_shm) == -1) {
      cerr << " 10 mpi_rank " << mpi_rank << " had close error" << endl;
      throw(-1);
    }

  }

#else

  ptr = new T[size][2];

#endif

  return(0);

}

int CTI_Mmap(T (*&ptr)[3],const size_t size) {

  assert(ptr == NULL);

#ifdef WITH_SHM

  if (mpi_size_shared == 1) {

    ptr = new T[size][3];

  }
  else {

    assert(mpi_size_shared > 1);

    int fd_shm = -1;
    // the first rank in the shared comm should create the shared memory...
    if (mpi_rank_shared == 0) {
      fd_shm = shm_open(SHM_NAME, O_CREAT | O_EXCL | O_RDWR, 0600);
      if (fd_shm == -1) {
	cerr << " 6 mpi_rank " << mpi_rank << " had shm_open error" << endl;
	throw(-1);
      }
      // set file size..
      if (ftruncate(fd_shm, sizeof(T)*3*size) == -1) {
	cerr << " 7 mpi_rank " << mpi_rank << " had ftrucate error" << endl;
	throw(-1);
      }
    }

    MPI_Barrier(mpi_comm_shared);

    if (mpi_rank_shared != 0) {
      fd_shm = shm_open(SHM_NAME, O_RDWR, 0600);
      if (fd_shm == -1) {
	cerr << " 8 mpi_rank " << mpi_rank << " had shm_open error" << endl;
	throw(-1);
      }
    }

    if (mpi_rank_shared == 0) {
      // only the first rank in the shared comm is able to write...
      ptr = (T (*)[3])mmap(NULL, sizeof(T)*3*size, PROT_READ | PROT_WRITE, MAP_SHARED, fd_shm, 0);
    }
    else {
      ptr = (T (*)[3])mmap(NULL, sizeof(T)*3*size, PROT_READ, MAP_SHARED, fd_shm, 0);
    }
    assert(ptr != NULL);

    MPI_Barrier(mpi_comm_shared);

    // don't need dev/shm file anymore, unlink removes it...
    if (mpi_rank_shared == 0) {
      if (shm_unlink(SHM_NAME) == -1) {
	cerr << " 9 mpi_rank " << mpi_rank << " had shm_unlink error" << endl;
	throw(-1);
      }
    }

    // also don't need file descriptor..
    if (close(fd_shm) == -1) {
      cerr << " 10 mpi_rank " << mpi_rank << " had close error" << endl;
      throw(-1);
    }

  }

#else

  ptr = new T[size][3];

#endif

  return(0);

}

int CTI_Munmap(T* &ptr,const size_t size) {

#ifdef WITH_SHM

  if (mpi_size_shared == 1) {

    delete[] ptr;

  }
  else {

    assert(mpi_size_shared > 1);

    if (munmap(ptr, sizeof(T)*size) == -1) {
      cerr << " 11 rank " << mpi_rank << " had munmap() error: " << strerror(errno) << endl;
      throw(-1);
    }

  }

#else

  delete[] ptr;

#endif

  ptr = NULL;

  return(0);

}

int CTI_Munmap(T (*&ptr)[2],const size_t size) {

#ifdef WITH_SHM

  if (mpi_size_shared == 1) {

    delete[] ptr;

  }
  else {

    assert(mpi_size_shared > 1);

    if (munmap(ptr, sizeof(T)*2*size) == -1) {
      cerr << " 12 rank " << mpi_rank << " had munmap() error: " << strerror(errno) << endl;
      throw(-1);
    }

  }

#else

  delete[] ptr;

#endif

  ptr = NULL;

  return(0);

}

int CTI_Munmap(T (*&ptr)[3],const size_t size) {

#ifdef WITH_SHM

  if (mpi_size_shared == 1) {

    delete[] ptr;

  }
  else {

    assert(mpi_size_shared > 1);

    if (munmap(ptr, sizeof(T)*3*size) == -1) {
      cerr << " 12 rank " << mpi_rank << " had munmap() error: " << strerror(errno) << endl;
      throw(-1);
    }

  }

#else

  delete[] ptr;

#endif

  ptr = NULL;

  return(0);

}

// in the including code, be sure to define T appropriately
int CTI_Mmap_rw(T* &ptr,const size_t size) {

  assert(ptr == NULL);

#ifdef WITH_SHM

  if (mpi_size_shared == 1) {

    ptr = new T[size];

  }
  else {

    assert(mpi_size_shared > 1);

    int fd_shm = -1;
    // the first rank in the shared comm should create the shared memory...
    if (mpi_rank_shared == 0) {
      fd_shm = shm_open(SHM_NAME, O_CREAT | O_EXCL | O_RDWR, 0600);
      if (fd_shm == -1) {
	cerr << " 13 mpi_rank " << mpi_rank << " had shm_open error: " << strerror(errno) << endl;
	throw(-1);
      }
      // set file size..
      if (ftruncate(fd_shm, sizeof(T)*size) == -1) {
	cerr << " 14 mpi_rank " << mpi_rank << " had ftruncate error: " << strerror(errno) << endl;
	throw(-1);
      }
    }

    MPI_Barrier(mpi_comm_shared);

    if (mpi_rank_shared != 0) {
      fd_shm = shm_open(SHM_NAME, O_RDWR, 0600);
      if (fd_shm == -1) {
	cerr << " 15 mpi_rank " << mpi_rank << " had shm_open error: " << strerror(errno) << endl;
	throw(-1);
      }
    }

    ptr = (T *)mmap(NULL, sizeof(T)*size, PROT_READ | PROT_WRITE, MAP_SHARED, fd_shm, 0);
    assert(ptr != NULL);

    MPI_Barrier(mpi_comm_shared);

    // don't need dev/shm file anymore, unlink removes it...
    if (mpi_rank_shared == 0) {
      if (shm_unlink(SHM_NAME) == -1) {
	cerr << " 16 mpi_rank " << mpi_rank << " had shm_unlink error: " << strerror(errno) << endl;
	throw(-1);
      }
    }

    // also don't need file descriptor..
    if (close(fd_shm) == -1) {
      cerr << " 17 mpi_rank " << mpi_rank << " had close error: " << strerror(errno) << endl;
      throw(-1);
    }

  }

#else

  ptr = new T[size];

#endif

  return(0);

}

int CTI_Mmap_rw(T (*&ptr)[2],const size_t size) {

  assert(ptr == NULL);

#ifdef WITH_SHM

  if (mpi_size_shared == 1) {

    ptr = new T[size][2];

  }
  else {

    assert(mpi_size_shared > 1);

    int fd_shm = -1;
    // the first rank in the shared comm should create the shared memory...
    if (mpi_rank_shared == 0) {
      fd_shm = shm_open(SHM_NAME, O_CREAT | O_EXCL | O_RDWR, S_IRUSR | S_IWUSR | S_IWGRP | S_IRGRP | S_IWOTH | S_IROTH );
      if (fd_shm == -1) {
	cerr << " 6 mpi_rank " << mpi_rank << " had shm_open error" << endl;
	throw(-1);
      }
      // set file size..
      if (ftruncate(fd_shm, sizeof(T)*2*size) == -1) {
	cerr << " 7 mpi_rank " << mpi_rank << " had ftruncate error" << endl;
	throw(-1);
      }
    }

    MPI_Barrier(mpi_comm_shared);

    if (mpi_rank_shared != 0) {
      fd_shm = shm_open(SHM_NAME, O_RDWR, S_IRUSR | S_IWUSR | S_IWGRP | S_IRGRP | S_IWOTH | S_IROTH );
      if (fd_shm == -1) {
	cerr << " 8 mpi_rank " << mpi_rank << " had shm_open error" << endl;
	throw(-1);
      }
      // set file size..
      // this broke for mac clang (we don't think we need it)
      //if (ftruncate(fd_shm, sizeof(T)*2*size) == -1) {
	//cerr << " 7 mpi_rank " << mpi_rank << " had ftruncate error" << endl;
	//throw(-1);
      //}
    }

    ptr = (T (*)[2])mmap(NULL, sizeof(T)*2*size, PROT_READ | PROT_WRITE, MAP_SHARED, fd_shm, 0);
    assert(ptr != NULL);

    MPI_Barrier(mpi_comm_shared);

    // don't need dev/shm file anymore, unlink removes it...
    if (mpi_rank_shared == 0) {
      if (shm_unlink(SHM_NAME) == -1) {
	cerr << " 9 mpi_rank " << mpi_rank << " had shm_unlink error" << endl;
	throw(-1);
      }
    }

    // also don't need file descriptor..
    if (close(fd_shm) == -1) {
      cerr << " 10 mpi_rank " << mpi_rank << " had close error" << endl;
      throw(-1);
    }

  }

#else

  ptr = new T[size][2];

#endif

  return(0);

}

int CTI_Mmap_rw(T (*&ptr)[3],const size_t size) {

  assert(ptr == NULL);

#ifdef WITH_SHM

  if (mpi_size_shared == 1) {

    ptr = new T[size][3];

  }
  else {

    assert(mpi_size_shared > 1);

    int fd_shm = -1;
    // the first rank in the shared comm should create the shared memory...
    if (mpi_rank_shared == 0) {
      fd_shm = shm_open(SHM_NAME, O_CREAT | O_EXCL | O_RDWR, S_IRUSR | S_IWUSR | S_IWGRP | S_IRGRP | S_IWOTH | S_IROTH );
      if (fd_shm == -1) {
	cerr << " 6 mpi_rank " << mpi_rank << " had shm_open error" << endl;
	throw(-1);
      }
      // set file size..
      if (ftruncate(fd_shm, sizeof(T)*3*size) == -1) {
	cerr << " 7 mpi_rank " << mpi_rank << " had ftruncate error" << endl;
	throw(-1);
      }
    }

    MPI_Barrier(mpi_comm_shared);

    if (mpi_rank_shared != 0) {
      fd_shm = shm_open(SHM_NAME, O_RDWR, S_IRUSR | S_IWUSR | S_IWGRP | S_IRGRP | S_IWOTH | S_IROTH );
      if (fd_shm == -1) {
	cerr << " 8 mpi_rank " << mpi_rank << " had shm_open error" << endl;
	throw(-1);
      }
      // set file size..
      // this broke for mac clang (we don't think we need it)
      //if (ftruncate(fd_shm, sizeof(T)*3*size) == -1) {
	//cerr << " 7 mpi_rank " << mpi_rank << " had ftruncate error" << endl;
	//throw(-1);
      //}
    }

    ptr = (T (*)[3])mmap(NULL, sizeof(T)*3*size, PROT_READ | PROT_WRITE, MAP_SHARED, fd_shm, 0);
    assert(ptr != NULL);

    MPI_Barrier(mpi_comm_shared);

    // don't need dev/shm file anymore, unlink removes it...
    if (mpi_rank_shared == 0) {
      if (shm_unlink(SHM_NAME) == -1) {
	cerr << " 9 mpi_rank " << mpi_rank << " had shm_unlink error" << endl;
	throw(-1);
      }
    }

    // also don't need file descriptor..
    if (close(fd_shm) == -1) {
      cerr << " 10 mpi_rank " << mpi_rank << " had close error" << endl;
      throw(-1);
    }

  }

#else

  ptr = new T[size][3];

#endif

  return(0);

}
