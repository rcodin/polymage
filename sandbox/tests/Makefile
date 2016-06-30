PYTEST=py.test-3.3

TESTS=`ls -1u | grep -i .py$ | tr '\n' ' '`

all: test

test:
	echo $(TESTS)
	@-for test in $(TESTS); do \
		$(PYTEST) -vv $$test; \
			done

clean:
	rm -rf __pycache__ *.pyc *.dot *.png *.cpp
