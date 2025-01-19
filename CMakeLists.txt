cmake_minimum_required(VERSION 3.10)
set(CMAKE_CXX_COMPILER /usr/bin/g++)
project(fast-sssp)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wshadow -Wno-sign-conversion -Werror=uninitialized \
-Werror=return-type -Werror=return-local-addr")

set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O1 -g -ggdb3")
set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g -ggdb3 -fsanitize=address,undefined -fno-omit-frame-pointer \
-fno-sanitize-recover=all -fmax-errors=2 -D_GLIBCXX_DEBUG -D_LIBCPP_DEBUG -D_GLIBCXX_ASSERTIONS")

add_executable(sssp main.cpp)
add_executable(bellmanford bellman-ford.cpp)