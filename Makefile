VERSION=3.3

INCLUDES=-I/usr/include/python$(VERSION)m/ 
CPPFLAGS=-shared -lpython3 -fPIC -O3 -std=c++11 
LDFLAGS=-lboost_system 

all: sandbox/dpfusion.so

sandbox/dpfusion.so: sandbox/dpfusion/dpfusion.cpp
	g++ $(INCLUDES) $(CPPFLAGS) -o $@ $< $(LDFLAGS)

clean:
	rm -f sandbox/dpfusion.so
