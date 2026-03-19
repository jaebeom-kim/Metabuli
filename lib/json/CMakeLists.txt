add_library(jsonlib
    json.h
    json.cpp
)
set_target_properties(jsonlib PROPERTIES COMPILE_FLAGS "${MMSEQS_CXX_FLAGS} -w" LINK_FLAGS "${MMSEQS_CXX_FLAGS} -w")