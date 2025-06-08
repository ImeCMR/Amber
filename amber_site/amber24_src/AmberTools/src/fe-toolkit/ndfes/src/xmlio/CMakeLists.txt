if(NOT TARGET xmlio)
  file(GLOB SOURCES *.cpp *.hpp)
  add_library(xmlio STATIC ${SOURCES})
endif()
