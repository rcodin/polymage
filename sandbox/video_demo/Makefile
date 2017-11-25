DIRS = bilateral_grid harris_corner local_laplacian unsharp_mask darknet

all: naive polymage

naive:
	    @-for d in $(DIRS); do \
			        make -C $$d naive; \
					        done

polymage:
	    @-for d in $(DIRS); do \
			        make -C $$d polymage; \
					        done

clean:
	    @-for d in $(DIRS); do \
			        make -C $$d clean; \
					        done
		rm -rf *.pyc __pycache__
