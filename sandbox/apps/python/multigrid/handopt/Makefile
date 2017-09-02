APPS = 2d \
	   3d 

all: 
	@-for d in $(APPS); do \
			make -C $$d $@; \
			done

clean:
	@-for d in $(APPS); do \
			make -C $$d clean; \
			done

