#ifndef CTIMEMORY_HPP
#define CTIMEMORY_HPP

#include <cstdlib>

#ifndef NO_CTI_PREFETCH
#define PFD_DISTANCE 8
#define cti_prefetch(a,b,c) __builtin_prefetch(a,b,c);
#else 
#define PFD_DISTANCE 0
#define cti_prefetch(a,b,c)
#endif

// some architectures can seg-fault when asked to 
// prefetch beyond memory bound...

template<typename T> 
T* align_new_(const size_t n) { 

  const size_t align_size = __alignof__(T);
  const size_t s  = sizeof(T) + sizeof(T)%align_size; // round up..
  assert( s%align_size == 0);
  char* pool      = (char*) malloc(s*n);
  //char* pool        = (char*) aligned_alloc(align_size,s*n);

  T* head = NULL;
  for (size_t i = 0; i < n; ++i) { 
    T* tmp = new (pool) T();
    if ( head == NULL)
      head = tmp;
    pool += align_size;
  }
  return head;
}
#endif
