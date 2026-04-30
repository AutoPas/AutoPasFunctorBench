include(FetchContent)
FetchContent_Declare(
        CLI11
        GIT_REPOSITORY https://github.com/CLIUtils/CLI11.git
        GIT_TAG        v2.6.2 # You can update this to the latest release if needed
)
FetchContent_MakeAvailable(CLI11)