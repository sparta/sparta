#!/bin/bash

# To run the ParaView tests, uncomment and set the following shell
# variables to your installed ParaView and Sparta
#
# SPARTA_EXECUTABLE=<PATH TO SPARTA EXECUTABLE>
# SPARTA_MPI_EXEC=<PATH TO MPI LAUNCH EXECUTABLE FOR SPARTA>
#
# PARAVIEW_PVPYTHON=<PATH TO PARAVIEW PVPYTHON EXECUTABLE>
# PARAVIEW_PVBATCH=<PATH TO PARAVIEW PVBATCH EXECUTABLE>
# PARAVIEW_MPI_EXEC=<PATH TO MPI LAUNCH EXECUTABLE FOR PARAVIEW>

usage() {
  echo ""
  echo "Edit the top section of runIntegrationTests.sh to set installed Paraview"
  echo "and Sparta executable paths for running tests."
  echo ""
}

get_abs_filename() {
    echo "$(cd "$(dirname "$1")" && pwd)/$(basename "$1")"
}

check_variable() {
    if [ $# = 1 ]; then
        echo ""
        echo $1 "is not set. Unable to run tests."
        usage
        exit 1
    fi
    
    if [ ! -x "$1" ]; then
        echo ""
        echo $2 "=" $1 "is not executable. Unable to run tests."
        usage
        exit 1
    fi
}

check_file() {
    if [ ! -s "$1" ]; then
        echo "File $1 does not exist. Test $CURRENT_TEST failed."
        echo "Exiting."
        exit 1
    fi
}

test_circle() {
    if [ ! -d "/circle" ]; then
        mkdir circle
    fi

    cd circle
    cp $CIRCLE_EXAMPLE_DIR/air.species .
    cp $CIRCLE_EXAMPLE_DIR/air.vss .
    cp $CIRCLE_EXAMPLE_DIR/data.circle .
    cp $INPUT_FILE_DIR/in.circle .
    cp $INPUT_FILE_DIR/circle.txt .
    cp $INPUT_FILE_DIR/circle_slice.txt .
    cp $INPUT_FILE_DIR/circle_read_grid.txt .
    cp $INPUT_FILE_DIR/circle_catalyst_script.py .
    cp $INPUT_FILE_DIR/in.adapt.static .

    CURRENT_TEST="surf2paraview circle"
    echo ""
    echo "Running Sparta circle example"
    echo ""
    $SPARTA_MPI_EXEC -np 1 $SPARTA_EXECUTABLE < in.circle
    echo ""
    echo "Finished"
    echo ""
    echo "Running surf2paraview"
    $PARAVIEW_PVPYTHON $SURF2PARAVIEW data.circle circle_surf -r tmp_surf.*
    echo ""
    echo "Checking $CURRENT_TEST output"
    check_file "circle_surf.pvd" 
    check_file "circle_surf/circle_surf_0.vtu" 
    check_file "circle_surf/circle_surf_100.vtu" 
    check_file "circle_surf/circle_surf_200.vtu" 
    check_file "circle_surf/circle_surf_300.vtu" 
    check_file "circle_surf/circle_surf_400.vtu" 
    check_file "circle_surf/circle_surf_500.vtu" 
    check_file "circle_surf/circle_surf_600.vtu" 
    check_file "circle_surf/circle_surf_700.vtu" 
    check_file "circle_surf/circle_surf_800.vtu" 
    check_file "circle_surf/circle_surf_900.vtu" 
    check_file "circle_surf/circle_surf_1000.vtu" 
    echo ""
    echo "$CURRENT_TEST test passed"

    CURRENT_TEST="surf2paraview circle exodus"
    echo ""
    echo "Running surf2paraview"
    $PARAVIEW_PVPYTHON $SURF2PARAVIEW data.circle circle_surf -r tmp_surf.* --exodus
    echo ""
    echo "Checking $CURRENT_TEST output"
    check_file "circle_surf.ex2" 
    echo ""
    echo "$CURRENT_TEST test passed"

    CURRENT_TEST="grid2paraview pvpython circle"
    echo ""
    echo "Checking $CURRENT_TEST output"
    $PARAVIEW_PVPYTHON $GRID2PARAVIEW circle.txt circle_grid_pvpython \
        -xc 10 -yc 10 -r tmp_flow.*
    echo ""
    echo "Checking $CURRENT_TEST output"
    check_file "circle_grid_pvpython.pvd" 
    check_file "circle_grid_pvpython/circle_grid_pvpython_0.pvtu" 
    check_file "circle_grid_pvpython/circle_grid_pvpython_1_0.vtu" 
    check_file "circle_grid_pvpython/circle_grid_pvpython_2_0.vtu" 
    check_file "circle_grid_pvpython/circle_grid_pvpython_500.pvtu" 
    check_file "circle_grid_pvpython/circle_grid_pvpython_1_500.vtu" 
    check_file "circle_grid_pvpython/circle_grid_pvpython_2_500.vtu" 
    check_file "circle_grid_pvpython/circle_grid_pvpython_1000.pvtu" 
    check_file "circle_grid_pvpython/circle_grid_pvpython_1_1000.vtu" 
    check_file "circle_grid_pvpython/circle_grid_pvpython_2_1000.vtu" 
    echo ""
    echo "$CURRENT_TEST test passed"

    CURRENT_TEST="grid2paraview pvpython slice circle"
    echo ""
    echo "Checking $CURRENT_TEST output"
    $PARAVIEW_PVPYTHON $GRID2PARAVIEW circle_slice.txt \
        circle_slice_grid_pvpython -xc 10 -yc 10 -r tmp_flow.*
    echo ""
    echo "Checking $CURRENT_TEST output"
    check_file "circle_slice_grid_pvpython.pvd" 
    check_file "circle_slice_grid_pvpython/circle_slice_grid_pvpython_0.pvtu"
    check_file "circle_slice_grid_pvpython/circle_slice_grid_pvpython_500.pvtu"
    check_file "circle_slice_grid_pvpython/circle_slice_grid_pvpython_0_400.vtu"
    check_file "circle_slice_grid_pvpython/circle_slice_grid_pvpython_1_400.vtu"
    check_file "circle_slice_grid_pvpython/circle_slice_grid_pvpython_2_400.vtu"
    check_file "circle_slice_grid_pvpython/circle_slice_grid_pvpython_3_400.vtu"
    check_file "circle_slice_grid_pvpython/circle_slice_grid_pvpython_1000.pvtu"
    echo ""
    echo "$CURRENT_TEST test passed"

    CURRENT_TEST="grid2paraview pvbatch slice circle"
    $PARAVIEW_MPI_EXEC -np 4 $PARAVIEW_PVBATCH -sym $GRID2PARAVIEW circle_slice.txt \
        circle_slice_grid_pvbatch -xc 10 -yc 10 -r tmp_flow.*
    echo ""
    echo "Checking $CURRENT_TEST output"
    check_file "circle_slice_grid_pvbatch.pvd" 
    check_file "circle_slice_grid_pvbatch/circle_slice_grid_pvbatch_0.pvtu"
    check_file "circle_slice_grid_pvbatch/circle_slice_grid_pvbatch_500.pvtu"
    check_file "circle_slice_grid_pvbatch/circle_slice_grid_pvbatch_0_600.vtu"
    check_file "circle_slice_grid_pvbatch/circle_slice_grid_pvbatch_1_600.vtu"
    check_file "circle_slice_grid_pvbatch/circle_slice_grid_pvbatch_2_600.vtu"
    check_file "circle_slice_grid_pvbatch/circle_slice_grid_pvbatch_3_600.vtu"
    check_file "circle_slice_grid_pvbatch/circle_slice_grid_pvbatch_1000.pvtu"
    echo ""
    echo "$CURRENT_TEST test passed"

    CURRENT_TEST="grid2paraview pvbatch circle"
    $PARAVIEW_MPI_EXEC -np 4 $PARAVIEW_PVBATCH -sym $GRID2PARAVIEW circle.txt \
        circle_grid_pvbatch -xc 10 -yc 10 -r tmp_flow.*
    echo ""
    echo "Checking $CURRENT_TEST output"
    check_file "circle_grid_pvbatch.pvd"
    check_file "circle_grid_pvbatch/circle_grid_pvbatch_0.pvtu"
    check_file "circle_grid_pvbatch/circle_grid_pvbatch_1_0.vtu"
    check_file "circle_grid_pvbatch/circle_grid_pvbatch_2_0.vtu"
    check_file "circle_grid_pvbatch/circle_grid_pvbatch_500.pvtu"
    check_file "circle_grid_pvbatch/circle_grid_pvbatch_1_500.vtu"
    check_file "circle_grid_pvbatch/circle_grid_pvbatch_2_500.vtu"
    check_file "circle_grid_pvbatch/circle_grid_pvbatch_1000.pvtu"
    check_file "circle_grid_pvbatch/circle_grid_pvbatch_1_1000.vtu"
    check_file "circle_grid_pvbatch/circle_grid_pvbatch_2_1000.vtu"
    echo ""
    echo "$CURRENT_TEST test passed"

    CURRENT_TEST="grid2paraview catalyst circle"
    $PARAVIEW_MPI_EXEC -np 4 $PARAVIEW_PVBATCH -sym $GRID2PARAVIEW \
        circle_read_grid.txt circle_grid_pvbatch -xc 10 -yc 10 \
        -c circle_catalyst_script.py -r tmp_flow.*
    echo ""
    echo "Checking $CURRENT_TEST output"
    check_file "RenderView1_0.png"
    check_file "RenderView1_1.png"
    check_file "RenderView1_2.png"
    check_file "RenderView1_3.png"
    check_file "RenderView1_4.png"
    check_file "RenderView1_5.png"
    check_file "RenderView1_6.png"
    check_file "RenderView1_7.png"
    check_file "RenderView1_8.png"
    check_file "RenderView1_9.png"
    check_file "RenderView1_10.png"
    echo ""
    echo "$CURRENT_TEST test passed"

    echo ""
    echo "Running Sparta circle adapt static example"
    echo ""
    $SPARTA_MPI_EXEC -np 4 $SPARTA_EXECUTABLE < in.adapt.static
    echo ""
    echo "Finished"
    echo ""

    CURRENT_TEST="grid2paraview pvpython circle adapt static"
    echo ""
    echo "Checking $CURRENT_TEST output"
    $PARAVIEW_PVPYTHON $GRID2PARAVIEW circle_read_grid.txt \
        circle_adapt_static_grid_pvpython -xc 10 -yc 10 -r tmp_flow.*
    echo ""
    echo "Checking $CURRENT_TEST output"
    check_file "circle_adapt_static_grid_pvpython.pvd"
    check_file "circle_adapt_static_grid_pvpython/circle_adapt_static_grid_pvpython_0_1100.vtu"
    check_file "circle_adapt_static_grid_pvpython/circle_adapt_static_grid_pvpython_4500.pvtu"
    echo ""
    echo "$CURRENT_TEST test passed"

    cd ../
}

test_sphere() {
    if [ ! -d "/sphere" ]; then
        mkdir sphere
    fi

    cd sphere
    cp $SPHERE_EXAMPLE_DIR/air.species .
    cp $SPHERE_EXAMPLE_DIR/air.vss .
    cp $SPHERE_EXAMPLE_DIR/data.sphere .
    cp $INPUT_FILE_DIR/in.sphere .
    cp $INPUT_FILE_DIR/sphere.txt .
    cp $INPUT_FILE_DIR/sphere_read_grid.txt .
    cp $INPUT_FILE_DIR/sphere_slice.txt .
    cp $INPUT_FILE_DIR/sphere_catalyst_script.py .

    CURRENT_TEST="surf2paraview sphere"
    echo ""
    echo "Running Sparta sphere example"
    echo ""
    $SPARTA_MPI_EXEC -np 1 $SPARTA_EXECUTABLE < in.sphere
    echo ""
    echo "Finished"
    echo ""
    echo "Running surf2paraview"
    $PARAVIEW_PVPYTHON $SURF2PARAVIEW data.sphere sphere_surf -r tmp_surf.*
    echo ""
    echo "Checking $CURRENT_TEST output"
    check_file "sphere_surf.pvd" 
    check_file "sphere_surf/sphere_surf_0.vtu" 
    check_file "sphere_surf/sphere_surf_100.vtu" 
    check_file "sphere_surf/sphere_surf_200.vtu" 
    check_file "sphere_surf/sphere_surf_300.vtu" 
    check_file "sphere_surf/sphere_surf_400.vtu" 
    check_file "sphere_surf/sphere_surf_500.vtu" 
    check_file "sphere_surf/sphere_surf_600.vtu" 
    check_file "sphere_surf/sphere_surf_700.vtu" 
    check_file "sphere_surf/sphere_surf_800.vtu" 
    check_file "sphere_surf/sphere_surf_900.vtu" 
    check_file "sphere_surf/sphere_surf_1000.vtu" 
    echo ""
    echo "$CURRENT_TEST test passed"

    CURRENT_TEST="surf2paraview sphere exodus"
    echo ""
    echo "Running surf2paraview"
    $PARAVIEW_PVPYTHON $SURF2PARAVIEW data.sphere sphere_surf -r tmp_surf.* --exodus
    echo ""
    echo "Checking $CURRENT_TEST output"
    check_file "sphere_surf.ex2"
    echo ""
    echo "$CURRENT_TEST test passed"

    CURRENT_TEST="grid2paraview pvpython sphere"
    echo ""
    echo "Checking $CURRENT_TEST output"
    $PARAVIEW_PVPYTHON $GRID2PARAVIEW sphere.txt sphere_grid_pvpython \
        -xc 10 -yc 10 -zc 10 -r tmp_flow.*
    echo ""
    echo "Checking $CURRENT_TEST output"
    check_file "sphere_grid_pvpython.pvd" 
    check_file "sphere_grid_pvpython/sphere_grid_pvpython_0.pvtu" 
    check_file "sphere_grid_pvpython/sphere_grid_pvpython_500.pvtu" 
    check_file "sphere_grid_pvpython/sphere_grid_pvpython_0_500.vtu" 
    check_file "sphere_grid_pvpython/sphere_grid_pvpython_1_500.vtu" 
    check_file "sphere_grid_pvpython/sphere_grid_pvpython_2_500.vtu" 
    check_file "sphere_grid_pvpython/sphere_grid_pvpython_3_500.vtu" 
    check_file "sphere_grid_pvpython/sphere_grid_pvpython_4_500.vtu" 
    check_file "sphere_grid_pvpython/sphere_grid_pvpython_5_500.vtu" 
    check_file "sphere_grid_pvpython/sphere_grid_pvpython_6_500.vtu" 
    check_file "sphere_grid_pvpython/sphere_grid_pvpython_7_500.vtu" 
    check_file "sphere_grid_pvpython/sphere_grid_pvpython_1000.pvtu" 
    echo ""
    echo "$CURRENT_TEST test passed"

    CURRENT_TEST="grid2paraview pvpython slice sphere"
    echo ""
    echo "Checking $CURRENT_TEST output"
    $PARAVIEW_PVPYTHON $GRID2PARAVIEW sphere_slice.txt \
        sphere_slice_grid_pvpython -xc 5 -yc 5 -zc 5 -r tmp_flow.*
    echo ""
    echo "Checking $CURRENT_TEST output"
    check_file "sphere_slice_grid_pvpython.pvd" 
    check_file "sphere_slice_grid_pvpython/sphere_slice_grid_pvpython_0.pvtu"
    check_file "sphere_slice_grid_pvpython/sphere_slice_grid_pvpython_500.pvtu"
    check_file "sphere_slice_grid_pvpython/sphere_slice_grid_pvpython_33_100.vtu"
    check_file "sphere_slice_grid_pvpython/sphere_slice_grid_pvpython_50_100.vtu"
    check_file "sphere_slice_grid_pvpython/sphere_slice_grid_pvpython_27_100.vtu"
    check_file "sphere_slice_grid_pvpython/sphere_slice_grid_pvpython_1000.pvtu"
    echo ""
    echo "$CURRENT_TEST test passed"

    CURRENT_TEST="grid2paraview pvbatch slice sphere"
    $PARAVIEW_MPI_EXEC -np 4 $PARAVIEW_PVBATCH -sym $GRID2PARAVIEW \
        sphere_slice.txt sphere_slice_grid_pvbatch -xc 10 -yc 10 -zc 10 \
        -r tmp_flow.*
    echo ""
    echo "Checking $CURRENT_TEST output"
    check_file "sphere_slice_grid_pvbatch.pvd"
    check_file "sphere_slice_grid_pvbatch/sphere_slice_grid_pvbatch_0.pvtu"
    check_file "sphere_slice_grid_pvbatch/sphere_slice_grid_pvbatch_500.pvtu"
    check_file "sphere_slice_grid_pvbatch/sphere_slice_grid_pvbatch_0_200.vtu"
    check_file "sphere_slice_grid_pvbatch/sphere_slice_grid_pvbatch_1_200.vtu"
    check_file "sphere_slice_grid_pvbatch/sphere_slice_grid_pvbatch_2_200.vtu"
    check_file "sphere_slice_grid_pvbatch/sphere_slice_grid_pvbatch_3_200.vtu"
    check_file "sphere_slice_grid_pvbatch/sphere_slice_grid_pvbatch_1000.pvtu"
    echo ""
    echo "Checking $CURRENT_TEST output"
    echo ""
    echo "$CURRENT_TEST test passed"

    CURRENT_TEST="grid2paraview pvbatch sphere"
    $PARAVIEW_MPI_EXEC -np 4 $PARAVIEW_PVBATCH -sym $GRID2PARAVIEW sphere.txt sphere_grid_pvbatch \
        -xc 10 -yc 10 -zc 10 -r tmp_flow.*
    echo ""
    echo "Checking $CURRENT_TEST output"
    check_file "sphere_grid_pvbatch.pvd"
    check_file "sphere_grid_pvbatch/sphere_grid_pvbatch_0.pvtu"
    check_file "sphere_grid_pvbatch/sphere_grid_pvbatch_500.pvtu"
    check_file "sphere_grid_pvbatch/sphere_grid_pvbatch_0_500.vtu"
    check_file "sphere_grid_pvbatch/sphere_grid_pvbatch_1_500.vtu"
    check_file "sphere_grid_pvbatch/sphere_grid_pvbatch_2_500.vtu"
    check_file "sphere_grid_pvbatch/sphere_grid_pvbatch_3_500.vtu"
    check_file "sphere_grid_pvbatch/sphere_grid_pvbatch_1000.pvtu"
    echo ""
    echo "$CURRENT_TEST test passed"

    CURRENT_TEST="grid2paraview catalyst sphere"
    $PARAVIEW_MPI_EXEC -np 4 $PARAVIEW_PVBATCH -sym $GRID2PARAVIEW sphere_read_grid.txt sphere_grid_pvbatch \
        -xc 10 -yc 10 -zc 10 -c sphere_catalyst_script.py -r tmp_flow.*
    echo ""
    echo "Checking $CURRENT_TEST output"
    check_file "RenderView1_0.png"
    check_file "RenderView1_1.png"
    check_file "RenderView1_2.png"
    check_file "RenderView1_3.png"
    check_file "RenderView1_4.png"
    check_file "RenderView1_5.png"
    check_file "RenderView1_6.png"
    check_file "RenderView1_7.png"
    check_file "RenderView1_8.png"
    check_file "RenderView1_9.png"
    check_file "RenderView1_10.png"
    echo ""
    echo "$CURRENT_TEST test passed"

    cd ../
}

check_variable $SPARTA_EXECUTABLE "SPARTA_EXECUTABLE"
check_variable $SPARTA_MPI_EXEC "SPARTA_MPI_EXEC"
check_variable $PARAVIEW_PVPYTHON "PARAVIEW_PVPYTHON"
check_variable $PARAVIEW_PVBATCH "PARAVIEW_PVBATCH"
check_variable $PARAVIEW_MPI_EXEC "PARAVIEW_MPI_EXEC"

GRID2PARAVIEW=$(get_abs_filename "grid2paraview.py")
SURF2PARAVIEW=$(get_abs_filename "surf2paraview.py")
CIRCLE_EXAMPLE_DIR=$(get_abs_filename "../../examples/circle/")
SPHERE_EXAMPLE_DIR=$(get_abs_filename "../../examples/sphere/")
INPUT_FILE_DIR=$(get_abs_filename "tests/input_files/")
TEST_DIR="integration_test_scratch_directory"
if [ ! -d "${TEST_DIR}" ]; then
    mkdir ${TEST_DIR}
fi
cd ${TEST_DIR}

test_circle
test_sphere

echo ""
echo "All Sparta/Paraview integration tests passed"
