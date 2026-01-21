// in the including code, be sure to define T appropriately
int CTI_Mmap(T* &ptr,const size_t size);
int CTI_Mmap(T (*&ptr)[2],const size_t size);
int CTI_Mmap(T (*&ptr)[3],const size_t size);
int CTI_Munmap(T* &ptr,const size_t size);
int CTI_Munmap(T (*&ptr)[2],const size_t size);
int CTI_Munmap(T (*&ptr)[3],const size_t size);
int CTI_Mmap_rw(T* &ptr,const size_t size);
int CTI_Mmap_rw(T (*&ptr)[2],const size_t size);
int CTI_Mmap_rw(T (*&ptr)[3],const size_t size);
