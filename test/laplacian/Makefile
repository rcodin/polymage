CXX=icpc

CFLAGS=-O3 -xhost -ipo -qopenmp

all: autovec handvec handvecg

autovec:  laplacian_unvectorized.cpp main.cpp simple_pool_allocator.cpp
	$(CXX) $(CFLAGS) simple_pool_allocator.cpp main.cpp $< -o $@

handvec:  laplacian_vectorized.cpp main.cpp simple_pool_allocator.cpp
	$(CXX) $(CFLAGS) simple_pool_allocator.cpp main.cpp $< -o $@

handvecg:  laplacian_vectorized_gather.cpp main.cpp simple_pool_allocator.cpp
	$(CXX) $(CFLAGS) simple_pool_allocator.cpp main.cpp $< -o $@

clean:
	rm -f autovec handvec handvecg
