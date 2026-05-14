#ifndef EDLIB_H
#define EDLIB_H

// New ALPSCore-free API in namespace `edlib::`. Always available.
#include "edlib/ChiLoc.h"
#include "edlib/CommonUtils.h"
#include "edlib/Dyson.h"
#include "edlib/ExecutionStatistic.h"
#include "edlib/GreensFunction.h"
#include "edlib/Hamiltonian.h"
#include "edlib/HubbardModel.h"
#include "edlib/Lanczos.h"
#include "edlib/Mesh.h"
#include "edlib/MeshFactory.h"
#include "edlib/Parameters.h"
#include "edlib/SingleImpurityAndersonModel.h"
#include "edlib/StaticObservables.h"

// Legacy ALPSCore-based API in namespace `EDLib::`. Only when
// EDLIB_WITH_ALPSCORE is defined (set by CMake when the option is ON).
#ifdef EDLIB_WITH_ALPSCORE
#include "edlib/alpscore/EDParams.h"
#include "edlib/alpscore/Hamiltonian.h"
#include "edlib/alpscore/SzSymmetry.h"
#include "edlib/alpscore/SOCRSStorage.h"
#include "edlib/alpscore/CRSStorage.h"
#include "edlib/alpscore/HubbardModel.h"
#include "edlib/alpscore/GreensFunction.h"
#include "edlib/alpscore/ChiLoc.h"
#include "edlib/alpscore/HDF5Utils.h"
#include "edlib/alpscore/SpinResolvedStorage.h"
#include "edlib/alpscore/StaticObservables.h"
#include "edlib/alpscore/MeshFactory.h"
#include "edlib/alpscore/ExecutionStatistic.h"
#endif

#endif //EDLIB_H
