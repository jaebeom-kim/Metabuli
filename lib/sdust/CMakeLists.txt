add_library(sdust
        sdust.c
        sdust.h
        kdq.h
        kvec.h
        kalloc.h)

# sdust.c is C, but it is -Wc++-compat clean and sdust.h wraps its API in
# extern "C". Compile it as C++ so the mmseqs C++ flags below apply to a single
# language (the metabuli target also builds everything as C++17).
set_source_files_properties(sdust.c PROPERTIES LANGUAGE CXX)

set_target_properties(sdust PROPERTIES COMPILE_FLAGS "${MMSEQS_CXX_FLAGS} -w" LINK_FLAGS "${MMSEQS_CXX_FLAGS} -w")
