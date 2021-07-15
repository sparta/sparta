
from mpi4py import MPI
from random import sample
import math

def parallel_sort(data, compare=None):
    _check_inputs(data, compare)
    pivots = _create_bucket_pivots(data, compare)
    if not pivots:
        return _gather_to_proc_zero_and_sort(data, compare)    
    buckets = _create_buckets(pivots, data, compare)
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

def _distribute_buckets(buckets, compare):
    COMM = MPI.COMM_WORLD
    result = []
    for idx, b in enumerate(buckets):
        bucket_data = COMM.gather(b, root = idx)
        if COMM.Get_rank() == idx:
            result = flatten_list(bucket_data)
            sort_list(result, compare)
    return result

def _create_buckets(pivots, data, compare):
    COMM = MPI.COMM_WORLD
    buckets = [[] for x in range(COMM.Get_size())]
    for d in data:
        found_bucket = False
        for idx, p in enumerate(pivots):
            less_than = d < p
            if compare is not None:
                less_than = compare(d, p) == -1
            if less_than:
                buckets[idx].append(d)
                found_bucket = True
                break
        if not found_bucket:
            buckets[-1].append(d)
    return buckets

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
