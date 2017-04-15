# 
# To use GNU compilers instead of Intel's C++ compiler, comment the two lines 
# below and uncomment the subsequent ones
#

#CXX=icpc
#CXX_FLAGS=-openmp -O3 -xhost -fPIC -shared

CXX=g++
CXX_FLAGS=-fopenmp -O3 -march=native -fPIC -shared

LIB_SRC=../simple_pool_allocator.cpp

all: $(APP)

$(APP): polymage naive

polymage: $(APP).so

naive: $(APP)_naive.so

$(APP).so: $(APP)_polymage.cpp
	$(CXX) $(CXX_FLAGS) $(LIB_SRC) $< -o $@

$(APP)_naive.so: $(APP)_naive.cpp
	$(CXX) $(CXX_FLAGS) $(LIB_SRC) $< -o $@

clean:
	rm -rf *.pyc *.so __pycache__
