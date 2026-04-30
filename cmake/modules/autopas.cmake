# autopas library
message(STATUS "Adding AutoPas.")

# Enable ExternalProject CMake module
include(FetchContent)

# Select https (default) or ssh path.
set(autopasRepoPath https://github.com/AutoPas/AutoPas.git)
if (GIT_SUBMODULES_SSH)
    set(autopasRepoPath git@github.com:AutoPas/AutoPas.git)
endif ()

FetchContent_Declare(
        autopasfetch
        GIT_REPOSITORY ${autopasRepoPath}
        GIT_TAG feature/3xa/atm-soa # Name of the branch with functor to test
)
# Populate dependency
FetchContent_MakeAvailable(autopasfetch)

# Disable warnings from the library target
target_compile_options(autopas PRIVATE -w)
# Disable warnings from included headers
get_target_property(propval autopas INTERFACE_INCLUDE_DIRECTORIES)
target_include_directories(autopas SYSTEM PUBLIC "${propval}")

# Get the current branch name
execute_process(
        COMMAND git rev-parse --abbrev-ref HEAD
        WORKING_DIRECTORY "${autopasfetch_SOURCE_DIR}"
        OUTPUT_VARIABLE AUTOPAS_BRANCH
        OUTPUT_STRIP_TRAILING_WHITESPACE
)
# Get the hash for the new version
execute_process(
        COMMAND git rev-parse HEAD
        WORKING_DIRECTORY "${autopasfetch_SOURCE_DIR}"
        OUTPUT_VARIABLE AUTOPAS_COMMIT_HASH
        OUTPUT_STRIP_TRAILING_WHITESPACE
)
