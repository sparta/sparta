add_test(NAME ablation_mpi_1rank
        COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 $<TARGET_FILE:${TARGET_SPARTA}> -in in.ablation.2d -log ${SPARTA_MACHINE}.ablation.2d.log)
set_tests_properties(ablation_mpi_1rank PROPERTIES
                    PASS_REGULAR_EXPRESSION ""
                    FAIL_REGULAR_EXPRESSION "error;Error;ERROR;exited on signal"
)
