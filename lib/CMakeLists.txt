add_library(common-lib CRSStorage.h CRSStorage.cpp fortranbinding.h darnoldi.f90)

target_link_libraries(common-lib ${extlibs})