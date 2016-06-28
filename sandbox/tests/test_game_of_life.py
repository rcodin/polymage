from __future__ import absolute_import, division, print_function

from fractions import Fraction
import sys
import subprocess
sys.path.insert(0, '../')

from compiler import *
from constructs import *

def test_game_of_life():

    R = Parameter(Int, "R")
    C = Parameter(Int, "C")
    x = Variable(Int, "x")
    y = Variable(Int, "y")

    row = Interval(Int, 0, R - 1)
    col = Interval(Int, 0, C - 1)

    img = Image(Int, "img", [R, C])

    life_bounds_cond = Condition(x, '>=', 1) & Condition(x, '<=', R - 2) & \
        Condition(y, '<=', C - 2) & Condition(y, '>=', 1)

    accum = Function(([x, y], [row, col]), Int, "accum")
    accum.defn = [Case(life_bounds_cond, (img(x, y + 1) +
                                     img(x + 1, y + 1) +
                                     img(x + 1, y) +
                                     img(x + 1, y - 1) +
                                     img(x, y - 1) +
                                     img(x - 1, y - 1) +
                                     img(x - 1, y) +
                                     img(x - 1, y + 1)))]

    live_cond = Condition(img(x, y), '==', 1)
    dead_cond = Condition(img(x, y), '==', 0)

    lonely_cond = live_cond & Condition(accum(x, y), '<', 2)
    propogate_cond = live_cond & \
        Condition(accum(x, y), '>=', 2) & \
        Condition(accum(x, y), '<=', 3)
    suffocate_cond = live_cond & Condition(accum(x, y), '>', 2)
    progeny_cond = dead_cond & Condition(accum(x, y), '==', 3)

    life = Function(([x, y], [row, col]), Int, "life")
    life.defn = [Case(life_bounds_cond,
                 Select(suffocate_cond, 0,
                        Select(lonely_cond, 0,
                               Select(propogate_cond, 1,
                                      Select(progeny_cond, 1, img(x, y))))))]

    groups = [(accum, life)]
    p_est = [(R, 5000), (C, 5000)]

    pipeline = buildPipeline([life],
                             grouping=groups,
                             param_estimates=p_est,
                             pipe_name="gameoflife")

    filename = "gameoflife_graph"
    dot_file = filename + ".dot"
    png_file = filename + ".png"
    g = pipeline.pipeline_graph
    g.write(filename + ".dot")
    dotty_str = "dot -Tpng " + dot_file + " -o " + png_file
    subprocess.check_output(dotty_str, shell=True)

    filename = 'gameoflife_naive.cpp'
    c_file = open(filename, 'w')
    c_file.write(pipeline.generate_code().__str__())
    c_file.close()
