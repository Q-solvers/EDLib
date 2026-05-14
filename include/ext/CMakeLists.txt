project(extlib CXX)

file(GLOB EXT_HEADERS *.h)

add_custom_target(extlib SOURCES ${EXT_HEADERS})
