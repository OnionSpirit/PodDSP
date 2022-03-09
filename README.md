# PodDSP
### Installation
Clone to project repository

Add to cmake

add_subdirectory(poddsp)/
target_include_directories(${PROJECT_NAME} PRIVATE poddsp/include)/
target_link_libraries(${PROJECT_NAME} PRIVATE poddsp)/
