include_directories(commons)
include_directories(workflow)
include_directories(functional)
include_directories(../lib/prodigal)
include_directories(../lib/kaiju/src)
include_directories(../lib/kaiju/src/include/ncbi-blast+)
include_directories(../lib/kaiju/src/include/)
add_subdirectory(commons)
add_subdirectory(util)
add_subdirectory(version)
add_subdirectory(workflow)
add_subdirectory(functional)

add_executable(metabuli
        ${commons_source_files}
        ${util_source_files}
        ${workflow_source_files}
        ${functional_source_files}
        metabuli.cpp
        LocalCommandDeclarations.h util/filter_by_genus.cpp)

mmseqs_setup_derived_target(metabuli)
target_link_libraries(metabuli version)
target_link_libraries(metabuli prodigal)
target_link_libraries(metabuli kaiju_lib)
add_dependencies(metabuli local-generated)

install(TARGETS metabuli DESTINATION bin)
