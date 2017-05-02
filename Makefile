all: 
	g++ -shared -I/usr/include/python3.3m/ -lpython3 -o sandbox/optgrouping.so sandbox/optgrouping_incremental.cpp -fPIC -O3 -std=c++11 -lboost_system

clean:
	rm sandbox/optgrouping.so
