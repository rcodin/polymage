from __future__ import absolute_import, division, print_function

from fractions import Fraction
import sys
import subprocess
sys.path.insert(0, '../')

from compiler import *
from constructs import *

# The purpose of this test file is to try and provide various kinds of
# pipeline DAGS and ensure, through stress testing, that PolyMage functions
# smoothly with practical/impractical corner cases.

R1 = Parameter(Int, "R1")
C1 = Parameter(Int, "C1")
R2 = Parameter(Int, "R2")
C2 = Parameter(Int, "C2")

x = Variable(Int, "x")
y = Variable(Int, "y")
z = Variable(Int, "z")

row1 = Interval(Int, 0, R1-1)
col1 = Interval(Int, 0, C1-1)
row2 = Interval(Int, 0, R2-1)
col2 = Interval(Int, 0, C2-1)

def code_and_graph_gen(pipeline, filename):
    dot_file = filename+".dot"
    png_file = filename+".png"
    g = pipeline.pipeline_graph
    g.write(filename+".dot")
    dotty_str = "dot -Tpng "+dot_file+" -o "+png_file
    subprocess.check_output(dotty_str, shell=True)

    c_file_name = filename+".cpp"
    c_file = open(c_file_name, 'w')
    c_file.write(pipeline.generate_code().__str__())
    c_file.close()

    return

def test_dag1():
    img1 = Image(Float, "input", [R1, C1])

    pipe1 = Function(([x, y], [row1, col1]), Float, "pipe1")
    pipe1.defn = [img1(x, y)]

    # build the pipeline
    pipeline = buildPipeline([pipe1],
                             pipe_name = "dag1")

    filename = "test_dag1"
    code_and_graph_gen(pipeline, filename)

    return

def test_dag2():
    img1 = Image(Float, "input", [R1, C1])
    img2 = Image(Float, "input", [R2, C2])

    cond = Condition(C1, "==", R2)

    pipe2 = Reduction(
        ([x, y], [row1, col2]), 
        ([x, z, y], [row1, col1, col2]), 
        Float, "pipe2")
    pipe2.defn = [ Case(cond, Reduce(pipe2(x, y), img1(x, z) + img2(z, y), Op.Sum)) ]

    # build the pipeline
    pipeline = buildPipeline([pipe2],
                             pipe_name = "dag2")

    filename = "test_dag2"
    code_and_graph_gen(pipeline, filename)

    return
