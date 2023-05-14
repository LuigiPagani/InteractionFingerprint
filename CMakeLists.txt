cmake_minimum_required (VERSION 3.11)
project(PP)

add_executable(main main.cpp)

set(RDKIT_INSTALL_PREFIX "/opt/luigi/rdkit/2022_09")
target_include_directories(main PRIVATE "${RDKIT_INSTALL_PREFIX}/include/rdkit")
target_link_directories(main PRIVATE "${RDKIT_INSTALL_PREFIX}/lib")
target_link_directories(main PRIVATE "${RDKIT_INSTALL_PREFIX}/lib64")
target_link_libraries(main
        RDKitChemReactions
        RDKitFileParsers
        RDKitSmilesParse
        RDKitDepictor
        RDKitRDGeometryLib
        RDKitRDGeneral
        RDKitGraphMol
        RDKitDistGeometry
        RDKitDistGeomHelpers
        RDKitMolAlign
        RDKitOptimizer
        RDKitForceField
        RDKitForceFieldHelpers
        RDKitAlignment
        RDKitForceField
        RDKitMolTransforms
        RDKitEigenSolvers
        RDKitDataStructs
        RDKitSubstructMatch
        )