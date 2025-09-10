CPT beauty 2Dmass and time fitter
=================================

Analysis code for CPT studies in semileptonic B0 decays

Install
-------

Requires recent GCC and ROOT. Tested with GCC 13.1 and ROOT 6.24.

### CIŚ cluster

```
srun -p INTEL_CASCADE -c 80 --pty bash
source multithreaded/setup.sh
make -j 10
```

2D Mass Fit
-----------

Read the [2D README](multithreaded/README.md).
