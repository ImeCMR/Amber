set(LIBTORCH_VERSION "1.12.1")
set(LIBTORCH_DEVICE_TYPE "cpu")
set(LIBTORCH_DEVICE_MAJOR_VERSION "")
set(LIBTORCH_DEVICE_MINOR_VERSION "")
set(LIBTORCH_CUDA_IS_AVALIABLE FALSE)

if(CUDA AND CUDNN)
    set(LIBTORCH_DEVICE_TYPE "cu")
    set(LIBTORCH_DEVICE_MAJOR_VERSION ${CUDA_VERSION_MAJOR})
    set(LIBTORCH_CUDA_IS_AVALIABLE TRUE)

    if("${CUDA_VERSION_MAJOR}" STREQUAL "11")
        set(LIBTORCH_DEVICE_MINOR_VERSION 3) # use CUDA 11.3
    elseif("${CUDA_VERSION_MAJOR}" STREQUAL "10")
        set(LIBTORCH_DEVICE_MINOR_VERSION 2) # use CUDA 10.2, win doesn't support cuda-10
        if(TARGET_WINDOWS)
            message(STATUS "Libtorch with CUDA-10 are no longer available for Windows, please use CUDA-11 to enable GPU computing.")
            set(LIBTORCH_CUDA_IS_AVALIABLE FALSE)
        endif()
    else()
        message(STATUS "CUDA version ${CUDA_VERSION_MAJOR} is too old to use libtorch on GPU.")
        set(LIBTORCH_CUDA_IS_AVALIABLE FALSE)
    endif()
endif()

if(NOT LIBTORCH_CUDA_IS_AVALIABLE)
    set(LIBTORCH_DEVICE_TYPE "cpu")
    set(LIBTORCH_DEVICE_MAJOR_VERSION "")
    set(LIBTORCH_DEVICE_MINOR_VERSION "")
endif()

# Version identifier
set(LIBTORCH_DEVICE_ID ${LIBTORCH_DEVICE_TYPE}${LIBTORCH_DEVICE_MAJOR_VERSION}${LIBTORCH_DEVICE_MINOR_VERSION})
set(LIBTORCH_FILE_ID "libtorch_${LIBTORCH_VERSION}${LIBTORCH_DEVICE_ID}")

# Summary Libtorch information
message(STATUS "Libtorch installation information: " "AMBER trys to install libtorch library built for ${LIBTORCH_DEVICE_ID}")

# Set download url https://pytorch.org/get-started/locally/
if(TARGET_WINDOWS)
    set(LIBTORCH_DOWNLOAD_URL "https://download.pytorch.org/libtorch/${LIBTORCH_DEVICE_ID}/libtorch-win-shared-with-deps-${LIBTORCH_VERSION}%2B${LIBTORCH_DEVICE_ID}.zip")
elseif(TARGET_LINUX)
    set(LIBTORCH_DOWNLOAD_URL "https://download.pytorch.org/libtorch/${LIBTORCH_DEVICE_ID}/libtorch-shared-with-deps-${LIBTORCH_VERSION}%2B${LIBTORCH_DEVICE_ID}.zip")
elseif(TARGET_OSX)
    message(STATUS "For mac os, Libtorch only supports CPU computing platform.")
    set(LIBTORCH_DOWNLOAD_URL "https://download.pytorch.org/libtorch/cpu/libtorch-macos-${LIBTORCH_VERSION}.zip")
endif()

# Get Libtorch zip
if(NOT (EXISTS ${CMAKE_CURRENT_BINARY_DIR}/${LIBTORCH_FILE_ID} AND IS_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${LIBTORCH_FILE_ID}))
    if(NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/${LIBTORCH_FILE_ID}.zip)
        message(STATUS "Start downloading LibTorch from " ${LIBTORCH_DOWNLOAD_URL} " ...")
        execute_process(COMMAND wget ${LIBTORCH_DOWNLOAD_URL} -O ${CMAKE_CURRENT_BINARY_DIR}/${LIBTORCH_FILE_ID}.zip COMMAND_ERROR_IS_FATAL ANY)
    endif()

    message(STATUS "Unzip LibTorch at " ${CMAKE_CURRENT_BINARY_DIR}/${LIBTORCH_FILE_ID}.zip " ...")
    execute_process(COMMAND unzip -q ${CMAKE_CURRENT_BINARY_DIR}/${LIBTORCH_FILE_ID}.zip -d ${CMAKE_CURRENT_BINARY_DIR}/${LIBTORCH_FILE_ID} COMMAND_ERROR_IS_FATAL ANY)
    message(STATUS "LibTorch is downloaded in " ${CMAKE_CURRENT_BINARY_DIR}/${LIBTORCH_FILE_ID})
endif()

# Patch update
execute_process(COMMAND python ${CMAKE_CURRENT_LIST_DIR}/fix_cuda.py ${CMAKE_CURRENT_BINARY_DIR}/${LIBTORCH_FILE_ID}/libtorch/share/cmake/Caffe2/public/cuda.cmake)

# Find Torch
# stack mkl from amber to protect libtorch env
list(REMOVE_ITEM CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/hanjianwei")
set(OLD_MKL_INCLUDE_DIR ${MKL_INCLUDE_DIR})
set(OLD_MKL_LIBRARIES ${MKL_LIBRARIES})
set(MKL_INCLUDE_DIR "")
set(MKL_LIBRARIES "")

# find libtorch
set(CONTAIN_INSTALL_PREFIX FALSE)
if(CMAKE_INSTALL_PREFIX)
    execute_process(COMMAND mkdir -p ${CMAKE_INSTALL_PREFIX}/lib COMMAND_ERROR_IS_FATAL ANY)
    execute_process(COMMAND cp -r ${CMAKE_CURRENT_BINARY_DIR}/${LIBTORCH_FILE_ID}/libtorch ${CMAKE_INSTALL_PREFIX}/lib/libtorch COMMAND_ERROR_IS_FATAL ANY)
    set(CONTAIN_INSTALL_PREFIX TRUE)
endif()

if(CONTAIN_INSTALL_PREFIX)
    set(Torch_HOME ${CMAKE_INSTALL_PREFIX}/lib/libtorch)
else()
    set(Torch_HOME ${CMAKE_CURRENT_BINARY_DIR}/${LIBTORCH_FILE_ID}/libtorch)
endif()
find_package(Torch REQUIRED PATHS ${Torch_HOME} NO_DEFAULT_PATH)

if (Torch_FOUND)
    message(STATUS "Libtorch: ${Torch_HOME}")
endif()

# recover original mkl vars for further amber building
set(MKL_INCLUDE_DIR ${OLD_MKL_INCLUDE_DIR})
set(MKL_LIBRARIES ${OLD_MKL_LIBRARIES})
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/hanjianwei")

# set export variable
set(LIBTORCH_INCLUDE_DIRS ${TORCH_INCLUDE_DIRS} CACHE INTERNAL "Libtorch header paths")
set(LIBTORCH_LIBRARIES ${TORCH_LIBRARIES} CACHE INTERNAL "Libtorch library paths")

# add_library(libtorch INTERFACE IMPORTED)
# set_target_properties(libtorch PROPERTIES INTERFACE_LINK_LIBRARIES torch INTERFACE_LINK_LIBRARIES torch_library)
# set_target_properties(libtorch PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${TORCH_INCLUDE_DIRS}")
