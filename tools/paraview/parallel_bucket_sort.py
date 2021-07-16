from __future__ import print_function
from mpi4py import MPI
from random import sample
import math
import os
from functools import cmp_to_key

def parallel_sort(data, compare=None, use_file_buckets=False):
    COMM = MPI.COMM_WORLD
    _check_inputs(data, compare)
    pivots = _create_bucket_pivots(data, compare)
    if COMM.Get_rank() == 0:
        print("Found " + str(len(pivots)) + " pivot(s) in sort data")
    if not pivots:
        return _gather_to_proc_zero_and_sort(data, compare)
    
    if use_file_buckets:
        return _sort_with_file_buckets(pivots, data, compare)
    else:
        return _sort_with_memory_buckets(pivots, data, compare)

def _sort_with_file_buckets(pivots, data, compare):
    _write_file_buckets(pivots, data, compare)
    result = _read_file_bucket()
    sort_list(result, compare)
    _remove_file_buckets()
    return result

def _sort_with_memory_buckets(pivots, data, compare):
    COMM = MPI.COMM_WORLD
    buckets = _create_buckets(pivots, data, compare)
    if COMM.Get_rank() == 0:
        print("Created " + str(len(buckets)) + " bucket(s) from sort data")
    return _distribute_buckets(buckets, compare)

def _gather_to_proc_zero_and_sort(data, compare):
    COMM = MPI.COMM_WORLD
    global_data = COMM.gather(data, root = 0)
    if COMM.Get_rank() == 0:
        result = flatten_list(global_data)
        sort_list(result, compare)
        return result
    else:
        return []

def _write_file_buckets(pivots, data, compare):
    COMM = MPI.COMM_WORLD
    _remove_file_buckets()
    buckets = _create_empty_bucket_list()
    buckets_size = 0
    if COMM.Get_rank() == 0:
        print("Writing "  + str(len(buckets)) + " file bucket(s)")
    count = 1
    for element in data:
        _put_element_in_bucket(element, pivots, buckets, compare)
        buckets_size += 1
        if buckets_size > 1000000:
            _empty_buckets(buckets)
            buckets_size = 0
        if COMM.Get_rank() == 0 and count % 100000 == 0:
            print("Bucketed " + str(count) + " cell(s) of " + str(len(data)))
        count += 1
    _empty_buckets(buckets)

def _empty_buckets(buckets):
    COMM = MPI.COMM_WORLD
    if COMM.Get_rank() == 0:
        print("Emptying bucket(s)")
    for idx, b in enumerate(buckets):
        filename = _get_bucket_file_name_for_rank(idx)
        with open(filename, 'a') as f:
            for element in b:
                f.write("%s\n" % element)
        del b[:]
    if COMM.Get_rank() == 0:
        print("Finished emptying bucket(s)")

def _read_file_bucket():
    COMM = MPI.COMM_WORLD
    COMM.Barrier()
    filename = _get_bucket_file_name_for_rank(COMM.Get_rank())
    result = []
    if COMM.Get_rank() == 0:
        print("Reading file bucket(s)")
    with open(filename, 'r') as f:
        for line in f:
            result.append(line.strip())
    return result

def _remove_file_buckets():
    COMM = MPI.COMM_WORLD
    COMM.Barrier()
    if COMM.Get_rank() == 0:
        for rank in range(COMM.Get_size()):
            filename = _get_bucket_file_name_for_rank(rank)
            if os.path.exists(filename):
                os.remove(filename)
    COMM.Barrier()

def _get_bucket_file_name_for_rank(rank):
    return "sort_bucket_rank_" + str(rank) + ".txt"

def _distribute_buckets(buckets, compare):
    COMM = MPI.COMM_WORLD
    result = []
    for idx, b in enumerate(buckets):
        bucket_data = COMM.gather(b, root = idx)
        if COMM.Get_rank() == idx:
            result = flatten_list(bucket_data)
            sort_list(result, compare)
        if COMM.Get_rank() == 0:
            print("Distributed buckets to rank " + str(idx))
    return result

def _create_buckets(pivots, data, compare):
    buckets = _create_empty_bucket_list()
    for element in data:
        _put_element_in_bucket(element, pivots, buckets, compare)
    return buckets

def _create_empty_bucket_list():
    COMM = MPI.COMM_WORLD
    return [[] for x in range(COMM.Get_size())]

def _put_element_in_bucket(element, pivots, buckets, compare):
    found_bucket = False
    for idx, p in enumerate(pivots):
        less_than = element < p
        if compare is not None:
            less_than = compare(element, p) == -1
        if less_than:
            buckets[idx].append(element)
            found_bucket = True
            break
    if not found_bucket:
        buckets[-1].append(element)

def _create_bucket_pivots(data, compare):
    COMM = MPI.COMM_WORLD
    oversample_size = _get_oversample_size(data)
    local_samples = []
    if oversample_size <= len(data):
        local_samples = sample(data, oversample_size)
    global_samples = flatten_list(COMM.allgather(local_samples))
    gs_length = len(global_samples)
    if COMM.Get_size() == 1 or gs_length == 0 or \
        gs_length != oversample_size*COMM.Get_size():
        return []
    else:
        sort_list(global_samples, compare)
        increment = int(gs_length/COMM.Get_size())
        return list(global_samples[i] for i in range(increment, gs_length, increment))

def _get_oversample_size(data):
    return_value = 0
    total_data_size = _get_total_data_size(data)
    if total_data_size > 1:
        return_value = int(12*math.log(total_data_size))
    return return_value

def _get_total_data_size(data):
    COMM = MPI.COMM_WORLD
    return COMM.allreduce(len(data), op=MPI.SUM)

def _check_inputs(data, compare):
    if type(data) is not list:
        raise ValueError("data input must be type list")
    if compare is not None and not callable(compare):
        raise ValueError("compare input must be a function")

def flatten_list(input_list):
    return [item for sublist in input_list for item in sublist]

def sort_list(input_list, compare):
    if compare is not None:
        input_list.sort(key=cmp_to_key(compare))
    else:
        input_list.sort()
