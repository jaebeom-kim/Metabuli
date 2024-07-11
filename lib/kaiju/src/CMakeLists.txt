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

set(BWT_SOURCES
    bwt/bwt.c
    bwt/compactfmi.c
    bwt/sequence.c
    bwt/suffixArray.c
)

set(KAIJU_SOURCES
    kaiju.cpp
    ReadItem.cpp
    Config.cpp
    ConsumerThread.cpp
    util.cpp
)

# Create the library
add_library(kaiju_lib STATIC
    ${BLAST_SOURCES}
    ${BWT_SOURCES}
    ${KAIJU_SOURCES}    
)

# Link libraries
target_link_libraries(kaiju_lib pthread z)

set_target_properties(kaiju_lib PROPERTIES COMPILE_FLAGS "${MMSEQS_CXX_FLAGS} -w" LINK_FLAGS "${MMSEQS_CXX_FLAGS} -w")