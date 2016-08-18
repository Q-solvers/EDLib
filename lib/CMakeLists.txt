add_library(common-lib CRSStorage.h CRSStorage.cpp fortranbinding.h Storage.cpp Storage.h dummy.f90 SOCRSStorage.h EDParams.h)

target_link_libraries(common-lib ${extlibs})