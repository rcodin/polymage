
CFLAGS=-O3 -xhost -ipo -qopenmp `pkg-config --cflags --libs opencv`

all: autovec handvec handvecg

autovec:  laplacian_unvectorized.cpp
	icpc $(CFLAGS) simple_pool_allocator.cpp main.cpp $< -o $@

handvec:  laplacian_vectorized.cpp
	icpc $(CFLAGS) simple_pool_allocator.cpp main.cpp $< -o $@

handvecg:  laplacian_vectorized_gather.cpp
	icpc $(CFLAGS) simple_pool_allocator.cpp main.cpp $< -o $@

clean:
	rm autovec
	rm handvec
	rm handvecg
