project(edlib CXX)

file(GLOB EDLIB_HEADERS *.h)

add_custom_target(edlib SOURCES ${EDLIB_HEADERS})
