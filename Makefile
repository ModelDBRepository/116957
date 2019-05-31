export C++ = g++
export RM= rm -f
export MAKE= make

all: testNT

testNT: testNT.cc
	$(C++) $(FLAGS) $(INCLUDE)  -o testNT testNT.cc

clean:
	$(RM) *.o testNT
