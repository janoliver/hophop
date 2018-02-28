Parallel execution
==================

*hophop* can use OpenMP to execute multiple realizations of a simulation in
parallel. This can be specified with the :code:`-P, --parallel` and the
:code:`-t, --nthreads` CLI parameters, in combination with :code:`-i, --nruns`.

Typically, one averages over many realizations of the system, so :code:`-i, --nruns`
is greater than 1. For best performance, the value of :code:`-i, --nruns` can be
evenly distributed between the number of threads.

