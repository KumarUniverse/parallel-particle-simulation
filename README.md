# parallel-particle-simulation
A parallel particle simulation using OpenMP and Open MPI.

## Modified Files
- common.cpp
- common.h
- serial.cpp
- openmp.cpp
- mpi.cpp

## Instructions for Making files
1. Run `source modules.sh` so that MPI works.
2. Run `make` on the directory containing the _Makefile_.

## Instructions for Testing
1. Run the `projects` command on your Bridges terminal and note down your default charging ID.
2. Edit line 2 of the bash script (ex: _auto-bridges-serial_) to include your default charging ID instead of the already present charging ID which does not work.
3. Run `sbatch auto-bridges-serial`
4. Run `squeue -u <bridges-username>`
4. Wait 5 - 10 minutes for your submission to get graded. Your submission will go from PD status (pending) to CF status (resources allocated) to then R status (running). Sometimes you may have to wait longer than 10 minutes, so maybe grab a cup of coffee or do something else.
5. Once that's finished, check the output file (_auto-particle-serial.stdout_) in the directory you ran the `sbatch` command from to see your grade.

## Useful Resources
- [OpenMP API 4.5 C/C++](https://www.openmp.org/wp-content/uploads/OpenMP-4.5-1115-CPP-web.pdf)
- [Open MPI v4.0.3 Documentation](https://www.open-mpi.org/doc/current/)
- [Intro to MPI by Dr. Ngo](https://wcupa-my.sharepoint.com/:p:/r/personal/lngo_wcupa_edu/_layouts/15/guestaccess.aspx?e=g2SbSs&share=EeJ37JWSfoZOvyPx2LOe1DEBHTqqpmfu02ESjC-UoKbwCw)
- [MPI_COMM_WORLD](https://www.codingame.com/playgrounds/349/introduction-to-mpi/mpi_comm_world-size-and-ranks)
- [MPI Scatter & Gather](https://mpitutorial.com/tutorials/mpi-scatter-gather-and-allgather/)