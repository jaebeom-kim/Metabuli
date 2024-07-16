# kaiju/CMakeLists.txt
cmake_minimum_required(VERSION 3.10)
project(kaiju)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Compiler flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -pthread -DNDEBUG")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -DNDEBUG")

# Include directories
include_directories(${PROJECT_SOURCE_DIR}/include)
include_directories(${PROJECT_SOURCE_DIR}/include/ncbi-blast+)

# Source files
set(BLAST_SOURCES
    include/ncbi-blast+/algo/blast/core/pattern.c
    include/ncbi-blast+/algo/blast/core/blast_posit.c
    include/ncbi-blast+/algo/blast/composition_adjustment/matrix_frequency_data.c
    include/ncbi-blast+/algo/blast/core/blast_dynarray.c
    include/ncbi-blast+/algo/blast/core/matrix_freq_ratios.c
    include/ncbi-blast+/algo/blast/core/blast_encoding.c
    include/ncbi-blast+/algo/blast/core/blast_stat.c
    include/ncbi-blast+/algo/blast/core/blast_filter.c
    include/ncbi-blast+/algo/blast/core/blast_util.c
    include/ncbi-blast+/algo/blast/core/blast_message.c
    include/ncbi-blast+/algo/blast/core/ncbi_erf.c
    include/ncbi-blast+/algo/blast/core/blast_options.c
    include/ncbi-blast+/algo/blast/core/ncbi_math.c
    include/ncbi-blast+/algo/blast/core/blast_program.c
    include/ncbi-blast+/algo/blast/core/ncbi_std.c
    include/ncbi-blast+/algo/blast/core/blast_psi_priv.c
    include/ncbi-blast+/util/tables/raw_scoremat.c
    include/ncbi-blast+/algo/blast/core/blast_query_info.c
    include/ncbi-blast+/algo/blast/core/blast_seg.c
)

set(BLAST_HEADERS
    include/ncbi-blast+/algo/blast/core/pattern.h
    include/ncbi-blast+/algo/blast/core/blast_posit.h
    include/ncbi-blast+/algo/blast/composition_adjustment/matrix_frequency_data.h
    include/ncbi-blast+/algo/blast/core/blast_dynarray.h
    include/ncbi-blast+/algo/blast/core/matrix_freq_ratios.h
    include/ncbi-blast+/algo/blast/core/blast_encoding.h
    include/ncbi-blast+/algo/blast/core/blast_stat.h
    include/ncbi-blast+/algo/blast/core/blast_filter.h
    include/ncbi-blast+/algo/blast/core/blast_util.h
    include/ncbi-blast+/algo/blast/core/blast_message.h
    include/ncbi-blast+/algo/blast/core/blast_options.h
    include/ncbi-blast+/algo/blast/core/ncbi_math.h
    include/ncbi-blast+/algo/blast/core/blast_program.h
    include/ncbi-blast+/algo/blast/core/ncbi_std.h
    include/ncbi-blast+/algo/blast/core/blast_psi_priv.h
    include/ncbi-blast+/util/tables/raw_scoremat.h
    include/ncbi-blast+/algo/blast/core/blast_query_info.h
    include/ncbi-blast+/algo/blast/core/blast_seg.h
)

set(BWT_SOURCES
    bwt/bwt.c
    bwt/compactfmi.c
    bwt/sequence.c
    bwt/suffixArray.c
)

set(BWT_HEADERS
    bwt/bwt.h
    bwt/compactfmi.h
    bwt/sequence.h
    bwt/suffixArray.h)

set(KAIJU_SOURCES
    ReadItem.cpp
    Config.cpp
    util.cpp
    ConsumerThread.cpp
    ConsumerThreadp.cpp
    ConsumerThreadx.cpp
)

set(KAIJU_HEADERS
    ReadItem.hpp
    Config.hpp
    util.hpp
    ConsumerThread.hpp
    ConsumerThreadp.hpp
    ConsumerThreadx.hpp
)


# Create the library
add_library(kaiju_lib STATIC
    ${BLAST_SOURCES}
    ${BLAST_HEADERS}
    ${BWT_SOURCES}
    ${BWT_HEADERS}
    ${KAIJU_SOURCES}    
    ${KAIJU_HEADERS}
)

# Link libraries
target_link_libraries(kaiju_lib pthread z)

set_target_properties(kaiju_lib PROPERTIES COMPILE_FLAGS "${MMSEQS_CXX_FLAGS} -w" LINK_FLAGS "${MMSEQS_CXX_FLAGS} -w")