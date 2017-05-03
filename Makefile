
VERSION=3.4

INCLUDES=-I/usr/include/python$(VERSION)m/ 
CPPFLAGS=-shared -lpython3 -fPIC -O3 -std=c++11 
LDFLAGS=-lboost_system 

all: sandbox/optgrouping.so

sandbox/optgrouping.so: sandbox/optgrouping_incremental.cpp
	g++ $(INCLUDES) $(CPPFLAGS) -o $@ $< $(LDFLAGS)

clean:
	rm -f sandbox/optgrouping.so
