add_library(common-lib CRSStorage.h CRSStorage.cpp)

target_link_libraries(common-lib ${extlibs})