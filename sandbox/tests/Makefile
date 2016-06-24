PYTEST=py.test-3.3

TESTS=`ls -1u | grep -i .py$ | tr '\n' ' '`
# TESTS=test_t_stencil.py test_stencil.py test_sampling.py test_math.py test_liveness_computation.py test_hist.py test_harris.py test_grouping.py test_expression_simplify.py test_expression.py test_constructs.py test_blur.py test_alignment.py

all: test

test:
	echo $(TESTS)
	@-for test in $(TESTS); do \
		$(PYTEST) -vv $$test; \
			done

clean:
	rm -rf __pycache__ *.pyc *.dot *.png *.cpp
