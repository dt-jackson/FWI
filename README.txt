Seismic Full Waveform Inversion with SPECFEM2D



To compile:
    Update CFLAGS, INCDIR, and LIBDIR in makefile if needed
    run "make"
    (There will likely be several warnings about unused variables/parameters.
     This is nothing to worry about, these are mostly from functions that are
     not fully implemented yet.)
    The executable bin/xfwi will be created

dependencies:
    qhull       (libqhull_r)
    fftw3       (libfftw3)
    lapacke     (liblapacke)
    openssl     (libcrypto)


Running: 
    -Adjust parameter, src_loc, and rec_loc files as needed
     (see comments in files for instructions)
    -Make sure SPECFEM2D files and directories are setup properly
     (see below)
    -run ./xfwi parameters





SPECFEM2D setup:
    One modification to the SPECFEM2D code must be made for the FWI
    program to run correctly. See the instruction in the specfem2d_modification
    directory.

    SPECFEM2D must be setup to run multiple simulations at once through the 
    parameter NUMBER_OF_SIMULTANEOUS_RUNS. This should be set to the number 
    of source locations and the number of processors should be adjusted 
    accordingly. Some example scripts are provided to help setup the
    directory structure and call SPECFEM2D.

    To setup the initial model, generate *.bin files using SPECFEM2D with the
    option 
        SAVE_MODEL = binary
    After the bins have been generated, set
        SAVE_MODEL = default
    and 
        MODEL = binary



