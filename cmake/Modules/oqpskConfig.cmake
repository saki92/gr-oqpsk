INCLUDE(FindPkgConfig)
PKG_CHECK_MODULES(PC_OQPSK oqpsk)

FIND_PATH(
    OQPSK_INCLUDE_DIRS
    NAMES oqpsk/api.h
    HINTS $ENV{OQPSK_DIR}/include
        ${PC_OQPSK_INCLUDEDIR}
    PATHS ${CMAKE_INSTALL_PREFIX}/include
          /usr/local/include
          /usr/include
)

FIND_LIBRARY(
    OQPSK_LIBRARIES
    NAMES gnuradio-oqpsk
    HINTS $ENV{OQPSK_DIR}/lib
        ${PC_OQPSK_LIBDIR}
    PATHS ${CMAKE_INSTALL_PREFIX}/lib
          ${CMAKE_INSTALL_PREFIX}/lib64
          /usr/local/lib
          /usr/local/lib64
          /usr/lib
          /usr/lib64
)

INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(OQPSK DEFAULT_MSG OQPSK_LIBRARIES OQPSK_INCLUDE_DIRS)
MARK_AS_ADVANCED(OQPSK_LIBRARIES OQPSK_INCLUDE_DIRS)

