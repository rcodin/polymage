from __future__ import absolute_import, division, print_function

from fractions import Fraction
import sys
import subprocess
sys.path.insert(0, '../')

import libpluto
import islpy as isl

def test_pluto_schedule():
    pluto = libpluto.PlutoFFI()

    print(pluto)
    ctx = isl.Context()
    schedule_str = "[p_0, p_1, p_2, p_3, p_4, p_5, p_7] -> { S_1[i0, i1] : i0 >= 0 and i0 <= p_0 and i1 >= 0 and i1 <= p_3 and p_2 >= 0; S_0[i0] : i0 >= 0 and i0 <= p_0}";
    deps_str = "[p_0, p_1, p_2, p_3, p_4, p_5, p_7] -> { S_0[i0] -> S_1[o0, o1] : (exists (e0 = [(p_7)/8]: 8o1 = -p_5 + p_7 + 8192i0 - 8192o0 and 8e0 = p_7 and i0 >= 0 and o0 <= p_0 and 8192o0 >= -8p_3 - p_5 + p_7 + 8192i0 and 8192o0 <= -p_5 + p_7 + 8192i0 and p_2 >= 0 and o0 >= 1 + i0)); S_1[i0, i1] -> S_0[o0] : (exists (e0 = [(p_1)/8], e1 = [(p_4)/8], e2 = [(-p_1 + p_7)/8184]: 8192o0 = p_5 - p_7 + 8192i0 + 8i1 and 8e0 = p_1 and 8e1 = p_4 and 8184e2 = -p_1 + p_7 and i1 >= 0 and 8i1 <= 8192p_0 - p_5 + p_7 - 8192i0 and 8184i1 >= 1024 + 1024p_1 - 1023p_5 - p_7 - 8380416i0 and p_2 >= 0 and p_7 <= -1 + p_5 and 8i1 >= 1 + 8p_3 + p_4 - p_5 - 8192i0 and i1 <= p_3 and i0 >= 0 and 8i1 >= 8192 - p_5 + p_7))}"
    schedule = isl.UnionSet.read_from_str(ctx, schedule_str)
    deps = isl.UnionMap.read_from_str(ctx, deps_str)

    opts = pluto.options_alloc()
    opts.tile = 0;
    opts.parallel = 0;
    opts.efup = 0;
    opts.debug = 0;
    opts.islsolve = 1;
    opts.glpksolve = 1;

    out_schedule = pluto.schedule(schedule, deps, opts)
    print(out_schedule)