# Generate spartagitversion.h with the git commit and branch of the SPARTA
# sources, as LAMMPS does for lmpgitversion.h.  This runs as a script
# (cmake -P) at every build, but the header is only replaced when its contents
# change, and only sparta.cpp includes it, so a new commit recompiles just that
# one file instead of every source.
#
# Input variables: SPARTA_DIR, SPARTA_GIT_HEADER, GIT_FOUND, GIT_EXECUTABLE

set(temp_git_commit "")
set(temp_git_branch "")

if(GIT_FOUND AND EXISTS ${SPARTA_DIR}/.git)
  execute_process(COMMAND ${GIT_EXECUTABLE} describe --always --dirty --abbrev=10
    OUTPUT_VARIABLE temp_git_commit
    ERROR_QUIET
    RESULT_VARIABLE temp_git_rc
    WORKING_DIRECTORY ${SPARTA_DIR}
    OUTPUT_STRIP_TRAILING_WHITESPACE)
  if(temp_git_rc EQUAL 0)
    execute_process(COMMAND ${GIT_EXECUTABLE} rev-parse --abbrev-ref HEAD
      OUTPUT_VARIABLE temp_git_branch
      ERROR_QUIET
      WORKING_DIRECTORY ${SPARTA_DIR}
      OUTPUT_STRIP_TRAILING_WHITESPACE)
  else()
    set(temp_git_commit "")
  endif()
endif()

# escape backslashes and double quotes (allowed in branch names) for C strings
foreach(var temp_git_commit temp_git_branch)
  string(REPLACE "\\" "\\\\" ${var} "${${var}}")
  string(REPLACE "\"" "\\\"" ${var} "${${var}}")
endforeach()

set(temp "#ifndef SPARTA_GIT_VERSION_H\n#define SPARTA_GIT_VERSION_H\n")
set(temp "${temp}const char *SPARTA_NS::SPARTA::git_commit() { return \"${temp_git_commit}\"; }\n")
set(temp "${temp}const char *SPARTA_NS::SPARTA::git_branch() { return \"${temp_git_branch}\"; }\n")
set(temp "${temp}#endif\n")

file(WRITE "${SPARTA_GIT_HEADER}.tmp" "${temp}")
execute_process(COMMAND ${CMAKE_COMMAND} -E copy_if_different "${SPARTA_GIT_HEADER}.tmp" "${SPARTA_GIT_HEADER}")
file(REMOVE "${SPARTA_GIT_HEADER}.tmp")
