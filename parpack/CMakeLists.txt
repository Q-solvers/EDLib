

include_directories(${MPI_Fortran_INCLUDE_PATH})

add_library(parpack pivout.f psmout.f psvout.f pdmout.f pdvout.f
            pssaitr.f pssapps.f pssaup2.f pssaupd.f psseigt.f psseupd.f pssgets.f psgetv0.f pslarnv.f psnorm2.f
            pdsaitr.f pdsapps.f pdsaup2.f pdsaupd.f pdseigt.f pdseupd.f pdsgets.f pdgetv0.f pdlarnv.f pdnorm2.f debug.h pdlamch10.f pslamch.f stat.h)

target_link_libraries(parpack ${ARPACK_LIB})