import argparse
import os
import sys
import json
import pathlib
import paraview.simple
import inspect

def main():
    args = get_command_line()
    test_file_name = get_test_file_name(args)
    if not does_file_exist(test_file_name):
        print_test_failed_message()
        sys.exit("Unable to open test file: " + test_file_name)
    test_suite = read_test_file(test_file_name)
    if not is_valid_test_suite(test_suite):
        print_test_failed_message()
        sys.exit("Invalid test file contents: " + test_file_name)
    run_sparta(args)
    run_tests(args, test_suite)
    print_test_passed_message()

def get_command_line():
    parser = argparse.ArgumentParser()
    parser.add_argument("num_procs", help="Number of MPI ranks", type=int)
    parser.add_argument("sparta_command", help="Command line to run Sparta")
    parser.add_argument("test_input_file", help="Sparta test input file")
    parser.add_argument("paraview_mpi_exec", help="ParaView MPI executable")
    parser.add_argument("pvpython_exe", help="pvpython executable path")
    parser.add_argument("pvbatch_exe", help="pvbatch executable path")
    parser.add_argument("grid2paraview_module", help="grid2paraview.py module path")
    parser.add_argument("surf2paraview_module", help="surf2paraview.py module path")
    parser.add_argument("grid2paraviewcells_module", help="grid2paraview_cells.py module path")
    return parser.parse_args()

def run_sparta(args):
    os.system("rm -fr tmp_surf.* tmp_flow.* log.sparta data.grid")
    os.system(args.sparta_command + " -in " + args.test_input_file)

def run_surf2_paraview(args):
    print(args.pvpython_exe + " " + args.surf2paraview_module + " data.circle circle_surf -r tmp_surf.*")
    os.system(args.pvpython_exe + " " + args.surf2paraview_module + " data.circle circle_surf -r tmp_surf.*")

def get_test_file_name(args):
    TEST_PREFIX = "in."
    TEST_SUFFIX = ".json"
    TEST_MPI = ".test.mpi_"
    test_file_name = args.test_input_file
    if test_file_name.startswith(TEST_PREFIX):
        test_file_name = test_file_name[len(TEST_PREFIX):]
    return test_file_name + TEST_MPI + str(args.num_procs) + TEST_SUFFIX

def does_file_exist(test_file_name):
    return os.path.exists(test_file_name)

def print_test_passed_message():
    print("passed;no failures")

def print_test_failed_message():
    print("FAILED")

def read_test_file(test_file_name):
    try:
        with open(test_file_name) as f:
            data = json.load(f)
        return data
    except Exception as e:
        print_test_failed_message()
        print(str(e))
        sys.exit("Unable to parse json in test file: " + test_file_name)

def is_valid_test_suite(test_suite):
    ret_val = True
    if "tests" not in test_suite:
        print('No "tests" key found with list of tests')
        ret_val = False
    else:
        for test in test_suite["tests"]:
            if not is_valid_test(test):
                ret_val = False
                break
    return ret_val

def is_valid_test(test):
    ret_val = True
    if "name" not in test:
        print('No "name" key found in test')
        ret_val = False
    elif "clean" not in test:
        print('No "clean" key found in test')
        ret_val = False
    elif "run" not in test:
        print('No "run" key found in test')
        ret_val = False
    elif "assert_files" in test and \
        not is_list_of_strings(test["assert_files"]):
            print('"assert_files" is not list of str')
            ret_val = False
    elif "assert_read" in test and \
        not is_valid_assert(test, "assert_read"):
        ret_val = False
    elif "assert_fetch" in test and \
        not is_valid_assert(test, "assert_fetch"):
        ret_val = False
    else:
        if not is_list_of_strings(test["run"]):
            print('"run" is not list of str')
            ret_val = False
    return ret_val

def is_valid_assert(test, assert_name):
    ret_val = True
    if "file" not in test[assert_name]:
        print('No "file" key in ' + assert_name)
        ret_val = False
    if "attributes" not in test[assert_name]:
        print('No "attributes" key in ' + assert_name)
        ret_val = False
    return ret_val

def is_list_of_strings(inp_list):
    return bool(inp_list) and all(isinstance(elem, str) for elem in inp_list)

def run_tests(args, test_suite):
    for test in test_suite["tests"]:
        try:
            print("Running ParaView test: " + test["name"])
            os.system(test["clean"])
            preprocess(test)
            os.system(get_run_command_line(args, test["run"]))
            assert_files(test)
            print("Assert files passed.")
            assert_read(test)
            print("Assert read passed.")
            assert_fetch(test)
            print("Assert fetch passed.")
        except Exception as e:
            print_test_failed_message()
            print(str(e))
            sys.exit()

def preprocess(test):
    if "preprocess" in test:
        os.system(test["preprocess"])

def assert_files(test):
    if "assert_files" in test:
        for f in test["assert_files"]:
            if not does_file_exist(f):
                raise Exception('assert_files file not found: ' + f)

def assert_read(test):
    if "assert_read" in test:
        pv = open_data_file(test["assert_fetch"]["file"])
        check_attributes(pv, test["assert_read"]["attributes"])
        
def assert_fetch(test):
    if "assert_fetch" in test:
        pv = open_data_file(test["assert_fetch"]["file"])
        d = paraview.servermanager.Fetch(pv)
        check_attributes(d, test["assert_fetch"]["attributes"])

def check_attributes(reader, attributes):
    for key, value in attributes.items():
        if not hasattr(reader, key):
            raise Exception("reader has no attribute: " + key)
        found_value = call_attribute(reader, key)
        if found_value != value:
            raise Exception(key + " mismatch, expected: " +\
                str(value) + ", found: " + str(found_value))

def call_attribute(reader, key):
    f_attr = getattr(reader, key)
    if isinstance(getattr(type(reader), key, None), property):
        return f_attr
    else:
        return f_attr()

def open_data_file(data_file_name):
    f_ext = pathlib.Path(data_file_name).suffix
    if f_ext == ".pvd":
        return paraview.simple.PVDReader(FileName=data_file_name)
    else:
        return paraview.simple.ExodusIIReader(FileName=data_file_name)

def get_run_command_line(args, run_command_list):
    command_line = ""
    for token in run_command_list:
        if token == "PARAVIEW_PVPYTHON": 
            command_line += " " + args.pvpython_exe
        elif token == "SURF2PARAVIEW":
            command_line += " " + args.surf2paraview_module
        elif token == "GRID2PARAVIEW":
            command_line += " " + args.grid2paraview_module
        elif token == "GRID2PARAVIEWCELLS":
            command_line += " " + args.grid2paraviewcells_module
        elif token == "PARAVIEW_MPI_EXEC":
            command_line += " " + args.paraview_mpi_exec + " -np " +\
                str(args.num_procs)
        elif token == "PARAVIEW_PVPBATCH":
            command_line += " " + args.pvbatch_exe + " --sym "
        else:
            command_line += " " + token
    return command_line
    
if __name__ == '__main__':
    main()
