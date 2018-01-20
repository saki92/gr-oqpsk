# Copyright 2011,2012 Free Software Foundation, Inc.
#
# This file is part of GNU Radio
#
# GNU Radio is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3, or (at your option)
# any later version.
#
# GNU Radio is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with GNU Radio; see the file COPYING.  If not, write to
# the Free Software Foundation, Inc., 51 Franklin Street,
# Boston, MA 02110-1301, USA.

########################################################################
# Setup library
########################################################################
include(GrPlatform) #define LIB_SUFFIX

include_directories(${Boost_INCLUDE_DIR})
link_directories(${Boost_LIBRARY_DIRS})

list(APPEND oqpsk_sources
    pfb_oqpsk_clock_sync_cc_impl.cc
    oqpsk_phase_sync_cc_impl.cc
)

set(oqpsk_sources "${oqpsk_sources}" PARENT_SCOPE)
if(NOT oqpsk_sources)
	MESSAGE(STATUS "No C++ sources... skipping lib/")
	return()
endif(NOT oqpsk_sources)

add_library(gnuradio-oqpsk SHARED ${oqpsk_sources})
target_link_libraries(gnuradio-oqpsk ${Boost_LIBRARIES} ${GNURADIO_ALL_LIBRARIES})
set_target_properties(gnuradio-oqpsk PROPERTIES DEFINE_SYMBOL "gnuradio_oqpsk_EXPORTS")

if(APPLE)
    set_target_properties(gnuradio-oqpsk PROPERTIES
        INSTALL_NAME_DIR "${CMAKE_INSTALL_PREFIX}/lib"
    )
endif(APPLE)

########################################################################
# Install built library files
########################################################################
install(TARGETS gnuradio-oqpsk
    LIBRARY DESTINATION lib${LIB_SUFFIX} # .so/.dylib file
    ARCHIVE DESTINATION lib${LIB_SUFFIX} # .lib file
    RUNTIME DESTINATION bin              # .dll file
)

########################################################################
# Build and register unit test
########################################################################
include(GrTest)

include_directories(${CPPUNIT_INCLUDE_DIRS})

list(APPEND test_oqpsk_sources
    ${CMAKE_CURRENT_SOURCE_DIR}/test_oqpsk.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/qa_oqpsk.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/qa_oqpsk_phase_sync_cc.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/qa_pfb_oqpsk_clock_sync_cc.cc
)

add_executable(test-oqpsk ${test_oqpsk_sources})

target_link_libraries(
  test-oqpsk
  ${GNURADIO_RUNTIME_LIBRARIES}
  ${Boost_LIBRARIES}
  ${CPPUNIT_LIBRARIES}
  gnuradio-oqpsk
)

GR_ADD_TEST(test_oqpsk test-oqpsk)
