CPT beauty 2Dmass and time fitter
=================================

Analysis code for CPT studies in semileptonic B0 decays

Install
-------

Requires recent GCC and ROOT. Tested with GCC 13.1 and ROOT 6.24.
For automatic style formating use clang-format. Tested with clang-format 19.

### CIŚ cluster

Compile the code at CIŚ:

```
srun -p INTEL_CASCADE -c 80 --pty bash
source setup.sh
make -j 10
```

Install clang-format at CIŚ:

```
srun --pty bash
source setup.sh
python3 -m venv venv
source venv/bin/activate
pip install -r requirements.txt
```

To apply automatic style formating run:

```
srun --pty bash
source setup.sh
make format
```

1D and 2D Mass Fits
-------------------

Read the [2D README](multithreaded/README.md).
