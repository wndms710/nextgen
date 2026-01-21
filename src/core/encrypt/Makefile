
include ../../../Makefile.in

OBJS=	tomcrypt.o tommath.o
default: libtomlight.a

libtomlight.a: $(OBJS)
	$(AR) $@ $^
	$(RANLIB) $@

tomcrypt.o: 	tomcrypt.cpp tomcrypt.hpp 
	$(CXX) -O0 $< -c -o $@

tommath.o: 	tommath.cpp tommath_superclass.hpp tommath_class.hpp 
	$(CXX) -O0 $< -c -o $@

##########################################################
#Serial Files
OBJS_S= tomcrypt_s.o tommath_s.o

libtomlight_s.a: $(OBJS_S)
	$(AR) $@ $^
	$(RANLIB) $@

%_s.o : %.cpp
	$(CXX_NO_MPI) -fPIC -O0 $< -c -o $@

tomcrypt_s.o: 	tomcrypt.hpp tomcrypt.cpp
tommath_s.o: 	tommath_superclass.hpp tommath_class.hpp tommath.cpp

clean: 
	rm -f *.o *.a 
