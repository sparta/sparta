add_test(NAME adapt_mpi_1rank
        COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 $<TARGET_FILE:${TARGET_SPARTA}> -in in.adapt.rotate -log ${SPARTA_MACHINE}.adapt.rotate.log)
