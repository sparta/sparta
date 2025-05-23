//@HEADER
// ************************************************************************
//
//                        Kokkos v. 4.0
//       Copyright (2022) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Part of Kokkos, under the Apache License v2.0 with LLVM Exceptions.
// See https://kokkos.org/LICENSE for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//@HEADER

#ifndef KOKKOS_SIMPLETASKSCHEDULER_HPP
#define KOKKOS_SIMPLETASKSCHEDULER_HPP

//----------------------------------------------------------------------------

#include <Kokkos_Macros.hpp>
#if defined(KOKKOS_ENABLE_TASKDAG)

#include <Kokkos_Core_fwd.hpp>
#include <Kokkos_TaskScheduler_fwd.hpp>
//----------------------------------------------------------------------------

#include <Kokkos_MemoryPool.hpp>

#include <Kokkos_Future.hpp>
#include <impl/Kokkos_TaskQueue.hpp>
#include <impl/Kokkos_SingleTaskQueue.hpp>
#include <impl/Kokkos_MultipleTaskQueue.hpp>
#include <impl/Kokkos_TaskQueueMultiple.hpp>
#include <impl/Kokkos_TaskPolicyData.hpp>
#include <impl/Kokkos_TaskTeamMember.hpp>
#include <impl/Kokkos_EBO.hpp>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#ifdef KOKKOS_ENABLE_DEPRECATION_WARNINGS
// We allow using deprecated classes in this file
KOKKOS_IMPL_DISABLE_DEPRECATED_WARNINGS_PUSH()
#endif

namespace Kokkos {

namespace Impl {

// TODO @tasking @cleanup move this
template <class T>
struct DefaultDestroy {
  T* managed_object;
  KOKKOS_FUNCTION
  void destroy_shared_allocation() { managed_object->~T(); }
};

}  // namespace Impl

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

template <class ExecSpace, class QueueType>
// requires ExecutionSpace<ExecSpace> && TaskQueue<QueueType>
class SimpleTaskScheduler
    : public Impl::TaskSchedulerBase,
      private Impl::ExecutionSpaceInstanceStorage<ExecSpace>,
      private Impl::MemorySpaceInstanceStorage<
          typename QueueType::memory_space>,
      private Impl::NoUniqueAddressMemberEmulation<
          typename QueueType::team_scheduler_info_type> {
 public:
  // TODO @tasking @generalization (maybe?) don't force QueueType to be complete
  // here

  using scheduler_type  = SimpleTaskScheduler;  // tag as scheduler concept
  using execution_space = ExecSpace;
  using task_queue_type = QueueType;
  using memory_space    = typename task_queue_type::memory_space;
  using memory_pool     = typename task_queue_type::memory_pool;

  using team_scheduler_info_type =
      typename task_queue_type::team_scheduler_info_type;
  using task_scheduling_info_type =
      typename task_queue_type::task_scheduling_info_type;
  using specialization = Impl::TaskQueueSpecialization<SimpleTaskScheduler>;
  using member_type    = typename specialization::member_type;

  template <class Functor>
  using runnable_task_type =
      typename QueueType::template runnable_task_type<Functor,
                                                      SimpleTaskScheduler>;

  using task_base_type = typename task_queue_type::task_base_type;
  using runnable_task_base_type =
      typename task_queue_type::runnable_task_base_type;

  using task_queue_traits = typename QueueType::task_queue_traits;

  template <class ValueType>
  using future_type = Kokkos::BasicFuture<ValueType, SimpleTaskScheduler>;
  template <class FunctorType>
  using future_type_for_functor = future_type<typename FunctorType::value_type>;

 private:
  template <typename, typename>
  friend class BasicFuture;

  using track_type = Kokkos::Impl::SharedAllocationTracker;
  using execution_space_storage =
      Impl::ExecutionSpaceInstanceStorage<execution_space>;
  using memory_space_storage = Impl::MemorySpaceInstanceStorage<memory_space>;
  using team_scheduler_info_storage =
      Impl::NoUniqueAddressMemberEmulation<team_scheduler_info_type>;

  track_type m_track;
  task_queue_type* m_queue = nullptr;

  KOKKOS_INLINE_FUNCTION
  static constexpr task_base_type* _get_task_ptr(std::nullptr_t) {
    return nullptr;
  }

  template <class ValueType>
  KOKKOS_INLINE_FUNCTION static constexpr task_base_type* _get_task_ptr(
      future_type<ValueType>&& f) {
    return f.m_task;
  }

  template <int TaskEnum, class DepTaskType, class FunctorType>
  KOKKOS_FUNCTION future_type_for_functor<std::decay_t<FunctorType>>
  _spawn_impl(
      DepTaskType arg_predecessor_task, TaskPriority arg_priority,
      typename runnable_task_base_type::function_type apply_function_ptr,
      typename runnable_task_base_type::destroy_type /*destroy_function_ptr*/,
      FunctorType&& functor) {
    KOKKOS_EXPECTS(m_queue != nullptr);

    using functor_future_type =
        future_type_for_functor<std::decay_t<FunctorType>>;
    using task_type =
        typename task_queue_type::template runnable_task_type<FunctorType,
                                                              scheduler_type>;

    // Reference count starts at two:
    //   +1 for the matching decrement when task is complete
    //   +1 for the future
    auto& runnable_task = *m_queue->template allocate_and_construct<task_type>(
        /* functor = */ std::forward<FunctorType>(functor),
        /* apply_function_ptr = */ apply_function_ptr,
        /* task_type = */ static_cast<Impl::TaskType>(TaskEnum),
        /* priority = */ arg_priority,
        /* queue_base = */ m_queue,
        /* initial_reference_count = */ 2);

    if (arg_predecessor_task != nullptr) {
      m_queue->initialize_scheduling_info_from_predecessor(
          runnable_task, *arg_predecessor_task);
      runnable_task.set_predecessor(*arg_predecessor_task);
      arg_predecessor_task->decrement_and_check_reference_count();
    } else {
      m_queue->initialize_scheduling_info_from_team_scheduler_info(
          runnable_task, team_scheduler_info());
    }

    auto rv = functor_future_type(&runnable_task);

    Kokkos::memory_fence();  // fence to ensure dependent stores are visible

    m_queue->schedule_runnable(std::move(runnable_task), team_scheduler_info());
    // note that task may be already completed even here, so don't touch it
    // again

    return rv;
  }

 public:
  //----------------------------------------------------------------------------
  // <editor-fold desc="Constructors, destructor, and assignment"> {{{2

  SimpleTaskScheduler() = default;

  explicit SimpleTaskScheduler(execution_space const& arg_execution_space,
                               memory_space const& arg_memory_space,
                               memory_pool const& arg_memory_pool)
      : execution_space_storage(arg_execution_space),
        memory_space_storage(arg_memory_space) {
    // Ask the task queue how much space it needs (usually will just be
    // sizeof(task_queue_type), but some queues may need additional storage
    // dependent on runtime conditions or properties of the execution space)
    auto const allocation_size = task_queue_type::task_queue_allocation_size(
        arg_execution_space, arg_memory_space, arg_memory_pool);

    // TODO @tasking @generalization DSH better encapsulation of the
    // SharedAllocationRecord pattern
    using record_type =
        Impl::SharedAllocationRecord<memory_space,
                                     Impl::DefaultDestroy<task_queue_type>>;

    // Allocate space for the task queue
    auto* record = record_type::allocate(memory_space(), "Kokkos::TaskQueue",
                                         allocation_size);
    m_queue      = new (record->data())
        task_queue_type(arg_execution_space, arg_memory_space, arg_memory_pool);
    record->m_destroy.managed_object = m_queue;
    m_track.assign_allocated_record_to_uninitialized(record);
  }

  explicit SimpleTaskScheduler(execution_space const& arg_execution_space,
                               memory_pool const& pool)
      : SimpleTaskScheduler(arg_execution_space, memory_space{},
                            pool) { /* forwarding ctor, must be empty */
  }

  explicit SimpleTaskScheduler(memory_pool const& pool)
      : SimpleTaskScheduler(execution_space{}, memory_space{},
                            pool) { /* forwarding ctor, must be empty */
  }

  SimpleTaskScheduler(memory_space const& arg_memory_space,
                      size_t const mempool_capacity,
                      unsigned const mempool_min_block_size,  // = 1u << 6
                      unsigned const mempool_max_block_size,  // = 1u << 10
                      unsigned const mempool_superblock_size  // = 1u << 12
                      )
      : SimpleTaskScheduler(
            execution_space{}, arg_memory_space,
            memory_pool(
                arg_memory_space, mempool_capacity, mempool_min_block_size,
                mempool_max_block_size,
                mempool_superblock_size)) { /* forwarding ctor, must be empty */
  }

  // </editor-fold> end Constructors, destructor, and assignment }}}2
  //----------------------------------------------------------------------------

  // Note that this is an expression of shallow constness
  KOKKOS_INLINE_FUNCTION
  task_queue_type& queue() const {
    KOKKOS_EXPECTS(m_queue != nullptr);
    return *m_queue;
  }

  KOKKOS_INLINE_FUNCTION
  SimpleTaskScheduler get_team_scheduler(int rank_in_league) const noexcept {
    KOKKOS_EXPECTS(m_queue != nullptr);
    auto rv = SimpleTaskScheduler{*this};
    rv.team_scheduler_info() =
        m_queue->initial_team_scheduler_info(rank_in_league);
    return rv;
  }

  KOKKOS_INLINE_FUNCTION
  execution_space const& get_execution_space() const {
    return this->execution_space_instance();
  }

  KOKKOS_INLINE_FUNCTION
  team_scheduler_info_type& team_scheduler_info() & {
    return this->team_scheduler_info_storage::no_unique_address_data_member();
  }

  KOKKOS_INLINE_FUNCTION
  team_scheduler_info_type const& team_scheduler_info() const& {
    return this->team_scheduler_info_storage::no_unique_address_data_member();
  }

  //----------------------------------------------------------------------------

  template <int TaskEnum, typename DepFutureType, typename FunctorType>
  KOKKOS_FUNCTION static Kokkos::BasicFuture<typename FunctorType::value_type,
                                             scheduler_type>
  spawn(Impl::TaskPolicyWithScheduler<TaskEnum, scheduler_type, DepFutureType>&&
            arg_policy,
        typename runnable_task_base_type::function_type arg_function,
        typename runnable_task_base_type::destroy_type arg_destroy,
        FunctorType&& arg_functor) {
    return std::move(arg_policy.scheduler())
        .template _spawn_impl<TaskEnum>(
            _get_task_ptr(std::move(arg_policy.predecessor())),
            arg_policy.priority(), arg_function, arg_destroy,
            std::forward<FunctorType>(arg_functor));
  }

  template <int TaskEnum, typename DepFutureType, typename FunctorType>
  KOKKOS_FUNCTION Kokkos::BasicFuture<typename FunctorType::value_type,
                                      scheduler_type>
  spawn(Impl::TaskPolicyWithPredecessor<TaskEnum, DepFutureType>&& arg_policy,
        FunctorType&& arg_functor) {
    static_assert(
        std::is_same_v<typename DepFutureType::scheduler_type, scheduler_type>,
        "Can't create a task policy from a scheduler and a future "
        "from a different scheduler");

    using task_type = runnable_task_type<FunctorType>;
    typename task_type::function_type const ptr = task_type::apply;
    typename task_type::destroy_type const dtor = task_type::destroy;

    auto const priority = arg_policy.priority();
    return _spawn_impl<TaskEnum>(std::move(arg_policy).predecessor().m_task,
                                 priority, ptr, dtor,
                                 std::forward<FunctorType>(arg_functor));
  }

  template <class FunctorType, class ValueType, class Scheduler>
  KOKKOS_FUNCTION static void respawn(
      FunctorType* functor,
      BasicFuture<ValueType, Scheduler> const& predecessor,
      TaskPriority priority = TaskPriority::Regular) {
    using task_type =
        typename task_queue_type::template runnable_task_type<FunctorType,
                                                              scheduler_type>;

    auto& task = *static_cast<task_type*>(functor);

    KOKKOS_EXPECTS(!task.get_respawn_flag());

    task.set_priority(priority);
    task.set_predecessor(*predecessor.m_task);
    task.set_respawn_flag(true);
  }

  template <class FunctorType>
  KOKKOS_FUNCTION static void respawn(
      FunctorType* functor, scheduler_type const&,
      TaskPriority priority = TaskPriority::Regular) {
    using task_type =
        typename task_queue_type::template runnable_task_type<FunctorType,
                                                              scheduler_type>;

    auto& task = *static_cast<task_type*>(functor);

    KOKKOS_EXPECTS(!task.get_respawn_flag());

    task.set_priority(priority);
    KOKKOS_ASSERT(!task.has_predecessor());
    task.set_respawn_flag(true);
  }

  template <class ValueType>
  KOKKOS_FUNCTION future_type<void> when_all(
      BasicFuture<ValueType, scheduler_type> const predecessors[],
      int n_predecessors) {
    // TODO @tasking @generalization DSH propagate scheduling info

    using task_type = typename task_queue_type::aggregate_task_type;

    future_type<void> rv;

    if (n_predecessors > 0) {
      task_queue_type* queue_ptr = nullptr;

      // Loop over the predecessors to find the queue and increment the
      // reference counts
      for (int i_pred = 0; i_pred < n_predecessors; ++i_pred) {
        auto* predecessor_task_ptr = predecessors[i_pred].m_task;

        if (predecessor_task_ptr != nullptr) {
          // TODO @tasking @cleanup DSH figure out when this is allowed to be
          // nullptr (if at all anymore)

          // Increment reference count to track subsequent assignment.
          // TODO @tasking @optimization DSH figure out if this reference count
          // increment is necessary
          predecessor_task_ptr->increment_reference_count();

          // TODO @tasking @cleanup DSH we should just set a boolean here
          // instead to make this more readable
          queue_ptr = m_queue;
        }

      }  // end loop over predecessors

      // This only represents a non-ready future if at least one of the
      // predecessors has a task (and thus, a queue)
      if (queue_ptr != nullptr) {
        auto& q = *queue_ptr;

        auto* aggregate_task_ptr =
            q.template allocate_and_construct_with_vla_emulation<
                task_type, task_base_type*>(
                /* n_vla_entries = */ n_predecessors,
                /* aggregate_predecessor_count = */ n_predecessors,
                /* queue_base = */ &q,
                /* initial_reference_count = */ 2);

        rv = future_type<void>(aggregate_task_ptr);

        for (int i_pred = 0; i_pred < n_predecessors; ++i_pred) {
          aggregate_task_ptr->vla_value_at(i_pred) =
              predecessors[i_pred].m_task;
        }

        Kokkos::memory_fence();  // we're touching very questionable memory, so
                                 // be sure to fence

        q.schedule_aggregate(std::move(*aggregate_task_ptr),
                             team_scheduler_info());
        // the aggregate may be processed at any time, so don't touch it after
        // this
      }
    }

    return rv;
  }

  template <class F>
  KOKKOS_FUNCTION future_type<void> when_all(int n_calls, F&& func) {
    // TODO @tasking @generalization DSH propagate scheduling info?

    // later this should be std::invoke_result_t
    using generated_type = decltype(func(0));
    using task_type      = typename task_queue_type::aggregate_task_type;

    static_assert(is_future<generated_type>::value,
                  "when_all function must return a Kokkos future (an instance "
                  "of Kokkos::BasicFuture)");

    // see #7779
    // There are issues with the implementation of std::is_base_of in NVCC
    // <= 12.5 for C++ 20
#if defined(KOKKOS_ENABLE_CXX17) || \
    (!defined(KOKKOS_COMPILER_NVCC) || KOKKOS_COMPILER_NVCC >= 1250)
    static_assert(
        std::is_base_of_v<scheduler_type,
                          typename generated_type::scheduler_type>,
        "when_all function must return a Kokkos::BasicFuture of a compatible "
        "scheduler type");
#endif

    auto* aggregate_task =
        m_queue->template allocate_and_construct_with_vla_emulation<
            task_type, task_base_type*>(
            /* n_vla_entries = */ n_calls,
            /* aggregate_predecessor_count = */ n_calls,
            /* queue_base = */ m_queue,
            /* initial_reference_count = */ 2);

    auto rv = future_type<void>(aggregate_task);

    for (int i_call = 0; i_call < n_calls; ++i_call) {
      auto generated_future = func(i_call);

      if (generated_future.m_task != nullptr) {
        generated_future.m_task->increment_reference_count();
        aggregate_task->vla_value_at(i_call) = generated_future.m_task;

        KOKKOS_ASSERT(m_queue ==
                          generated_future.m_task->ready_queue_base_ptr() &&
                      "Queue mismatch in when_all");
      }
    }

    Kokkos::memory_fence();

    m_queue->schedule_aggregate(std::move(*aggregate_task),
                                team_scheduler_info());
    // This could complete at any moment, so don't touch anything after this

    return rv;
  }
};

template <class ExecSpace, class QueueType>
inline void wait(SimpleTaskScheduler<ExecSpace, QueueType> const& scheduler) {
  using scheduler_type = SimpleTaskScheduler<ExecSpace, QueueType>;
  scheduler_type::specialization::execute(scheduler);
}

}  // namespace Kokkos

#ifdef KOKKOS_ENABLE_DEPRECATION_WARNINGS
KOKKOS_IMPL_DISABLE_DEPRECATED_WARNINGS_POP()
#endif

//----------------------------------------------------------------------------
//---------------------------------------------------------------------------#endif
///* #if defined( KOKKOS_ENABLE_TASKDAG ) */

#endif /* #if defined( KOKKOS_ENABLE_TASKDAG ) */
#endif /* #ifndef KOKKOS_SIMPLETASKSCHEDULER_HPP */
