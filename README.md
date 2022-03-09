# PodDSP
## Installation
#### 1)Clone to project repository

#### 2)Add to cmake:

  add_subdirectory(poddsp)  
  target_include_directories(${PROJECT_NAME} PRIVATE poddsp/include)  
  target_link_libraries(${PROJECT_NAME} PRIVATE poddsp)  
