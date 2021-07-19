from __future__ import print_function
from mpi4py import MPI
from random import sample
import math
import os
from functools import cmp_to_key

def parallel_sort(data, compare=None, use_file_buckets=False):
    pivots = _parallel_sort_init(data, compare)
    if not pivots:
        return _gather_to_proc_zero_and_sort(data, compare)
    
    if use_file_buckets:
        return _sort_with_file_buckets(pivots, data, compare)
    else:
        return _sort_with_memory_buckets(pivots, data, compare)

def parallel_sort_to_file_buckets(data, compare=None, prefix=None):
    pivots = _parallel_sort_init(data, compare)
    _write_file_buckets(pivots, data, compare, prefix)

def _parallel_sort_init(data, compare):
    _check_inputs(data, compare)
    pivots = _create_bucket_pivots(data, compare)
    if is_rank_zero():
        print("Found " + str(len(pivots)) + " pivot(s) in sort data")
    return pivots

def _sort_with_file_buckets(pivots, data, compare):
    _write_file_buckets(pivots, data, compare)
    result = _read_file_bucket()
    sort_list(result, compare)
    remove_file_buckets()
    return result

def _sort_with_memory_buckets(pivots, data, compare):
    buckets = _create_buckets(pivots, data, compare)
    if is_rank_zero():
        print("Created " + str(len(buckets)) + " bucket(s) from sort data")
    return _distribute_buckets(buckets, compare)

def _gather_to_proc_zero_and_sort(data, compare):
    global_data = get_comm_world().gather(data, root = 0)
    if is_rank_zero():
        result = flatten_list(global_data)
        sort_list(result, compare)
        return result
    else:
        return []

def _write_file_buckets(pivots, data, compare, prefix=None):
    remove_file_buckets(prefix)
    buckets = _create_empty_bucket_list()
    buckets_size = 0
    if is_rank_zero():
        print("Writing "  + str(len(buckets)) + " file bucket(s)")
    count = 1
    for element in data:
        _put_element_in_bucket(element, pivots, buckets, compare)
        buckets_size += 1
        if buckets_size > 1000000:
            _empty_buckets(buckets, prefix)
            buckets_size = 0
        if is_rank_zero() and count % 100000 == 0:
            print("Bucketed " + str(count) + " cell(s) of " + str(len(data)))
        count += 1
    _empty_buckets(buckets, prefix)
    barrier()

def _empty_buckets(buckets, prefix=None):
    if is_rank_zero():
        print("Emptying bucket(s)")
    for idx, b in enumerate(buckets):
        filename = get_bucket_file_name_for_rank(idx, prefix)
        with open(filename, 'a') as f:
            for element in b:
                f.write("%s\n" % element)
        del b[:]
    if is_rank_zero():
        print("Finished emptying bucket(s)")

def _read_file_bucket():
    filename = get_bucket_file_name_for_rank(get_rank())
    result = []
    if is_rank_zero():
        print("Reading file bucket(s)")
    with open(filename, 'r') as f:
        for line in f:
            result.append(line.strip())
    return result

def remove_file_buckets(prefix=None):
    barrier()
    if is_rank_zero():
        for rank in range(get_size()):
            filename = get_bucket_file_name_for_rank(rank, prefix)
            if os.path.exists(filename):
                os.remove(filename)
    barrier()

def get_bucket_file_name_for_rank(rank, prefix=None):
    if prefix is not None:
        return prefix + "_sort_bucket_rank_" + str(rank) + ".txt"
    else:
        return "sort_bucket_rank_" + str(rank) + ".txt"

def _distribute_buckets(buckets, compare):
    result = []
    for idx, b in enumerate(buckets):
        bucket_data = get_comm_world().gather(b, root = idx)
        if get_rank() == idx:
            result = flatten_list(bucket_data)
            sort_list(result, compare)
        if is_rank_zero():
            print("Distributed buckets to rank " + str(idx))
    return result

def _create_buckets(pivots, data, compare):
    buckets = _create_empty_bucket_list()
    for element in data:
        _put_element_in_bucket(element, pivots, buckets, compare)
    return buckets

def _create_empty_bucket_list():
    return [[] for x in range(get_size())]

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
    oversample_size = _get_oversample_size(data)
    local_samples = []
    if oversample_size <= len(data):
        local_samples = sample(data, oversample_size)
    global_samples = flatten_list(get_comm_world().allgather(local_samples))
    gs_length = len(global_samples)
    if get_size() == 1 or gs_length == 0 or \
        gs_length != oversample_size*get_size():
        return []
    else:
        sort_list(global_samples, compare)
        increment = int(gs_length/get_size())
        return list(global_samples[i] for i in range(increment, gs_length, increment))

def _get_oversample_size(data):
    return_value = 0
    total_data_size = _get_total_data_size(data)
    if total_data_size > 1:
        return_value = int(12*math.log(total_data_size))
    return return_value

def _get_total_data_size(data):
    return get_comm_world().allreduce(len(data), op=MPI.SUM)

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

def is_rank_zero():
    return get_comm_world().Get_rank() == 0

def get_rank():
    return get_comm_world().Get_rank()

def get_size():
    return get_comm_world().Get_size() 

def barrier():
    get_comm_world().Barrier()

def error_found_on_rank_zero(error_flag):
    return get_comm_world().bcast(error_flag, root = 0)

def get_comm_world():
    return MPI.COMM_WORLD
