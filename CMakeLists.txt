## Copyright (C) 2011 ENIB
##   Ecole Nationale d'Ingenieurs de Brest (ENIB)
##   CS 73862 - 29238 BREST Cedex 3 - France
##   Tel: +33(0)298 05 89 89, Fax: +33(0)298 05 89 79, e-mail: combes@enib.fr
##
## This file is part of MACMA.
##
##   MACMA is free software: you can redistribute it and/or modify
##   it under the terms of the GNU General Public License as published by
##   the Free Software Foundation, either version 3 of the License, or
##   (at your option) any later version.
##
##   MACMA is distributed in the hope that it will be useful,
##   but WITHOUT ANY WARRANTY; without even the implied warranty of
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##   GNU General Public License for more details.
##
##   You should have received a copy of the GNU General Public License
##   along with MACMA.  If not, see <http://www.gnu.org/licenses/>.
##
##   author: leyaouanq@cervval.com

cmake_minimum_required(VERSION 2.6)

if(COMMAND cmake_policy)
    cmake_policy(SET CMP0003 NEW)
    cmake_policy(SET CMP0011 NEW)
endif(COMMAND cmake_policy)

project(MACMA)

if(APPLE)
    set(MACOSX_BUNDLE_SHORT_VERSION_STRING 0.1)
    set(MACOSX_BUNDLE_VERSION 0.1)
    set(MACOSX_BUNDLE_LONG_VERSION_STRING Version 0.1)
    set(CMAKE_OSX_ARCHITECTURES x86_64)
    set(MACOSX_BUNDLE_ICON_FILE icon.icns) 
    set(MAC_ICON_FILE ihm/icon.icns)
    set_source_files_properties(${MAC_ICON_FILE} PROPERTIES MACOSX_PACKAGE_LOCATION Resources)
endif(APPLE)

set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

# Dependencies
#--------------

### QT
Find_Package(Qt4 4.7.4 REQUIRED)
if ( NOT QT4_FOUND )
    message(FATAL_ERROR "QT4 not found !")
endif ( NOT QT4_FOUND )
include_directories(${QT4_INCLUDE_DIR})

include(${QT_USE_FILE})
add_definitions(${QT_DEFINITIONS})

### Boost
set(Boost_USE_MULTITHREADED      ON)
set(Boost_ADDITIONAL_VERSIONS 1.46.0)
include(FindBoost)

Find_Package(Boost 1.40.0 COMPONENTS thread iostreams system filesystem)
if(Boost_FOUND)
	include_directories(${Boost_INCLUDE_DIRS})
else(Boost_FOUND)
	message(FATAL_ERROR "Boost not found !")
endif()

### APPLE only
#if(APPLE)
#    SET(CMAKE_CXX_FLAGS "-arch i386 -arch x86_64 -sdk /Developer/SDKs/MacOSX10.6.sdk")
#endif(APPLE)

# Files
#--------------
file(
        GLOB_RECURSE
        MACMA_SOURCES
        src/*.cpp
)

file(
        GLOB_RECURSE
        MACMA_HEADERS
        inc/*.h
)

set(MACMA_FORMS ihm/MACMA.ui)
set(MACMA_RESOURCES ihm/resources.qrc)

QT4_WRAP_CPP(MACMA_HEADERS_MOC ${MACMA_HEADERS})
QT4_WRAP_UI(MACMA_FORMS_HEADERS ${MACMA_FORMS})
QT4_ADD_RESOURCES(MACMA_RESOURCES_RCC ${MACMA_RESOURCES})

INCLUDE_DIRECTORIES (${CMAKE_CURRENT_BINARY_DIR})

# Targets
#--------
include_directories(inc)

if(APPLE)
    add_executable(MACMA 
        ### Uncomment next line to build the Mac OS X package
        #MACOSX_BUNDLE
        ${MAC_ICON_FILE}
	    ${MACMA_SOURCES} 
        ${MACMA_HEADERS_MOC} 
        ${MACMA_FORMS_HEADERS} 
        ${MACMA_RESOURCES_RCC}
    )
else(APPLE)
	if(WIN32)
		add_executable(MACMA 
			WIN32
			${MACMA_SOURCES} 
			${MACMA_HEADERS_MOC} 
			${MACMA_FORMS_HEADERS} 
			${MACMA_RESOURCES_RCC}
		)
	else(WIN32)
		add_executable(MACMA 
			${MACMA_SOURCES} 
			${MACMA_HEADERS_MOC} 
			${MACMA_FORMS_HEADERS} 
			${MACMA_RESOURCES_RCC}
		)
	endif(WIN32)
endif(APPLE)

target_link_libraries(
        MACMA
	${QT_LIBRARIES}
        ${Boost_LIBRARIES}
)

# Installation
#-------------
SET(CPACK_PACKAGE_FILE_NAME "MACMA-installer")
SET(CPACK_PACKAGE_INSTALL_DIRECTORY "MACMA")
SET(CPACK_PACKAGE_INSTALL_REGISTRY_KEY "MACMA")
SET(CPACK_PACKAGE_NAME "MACMA")

if(APPLE)
	SET(CPACK_BINARY_PACKAGEMAKER ON)
	SET(CPACK_INSTALL_PREFIX "/Applications")
	SET(CPACK_PACKAGING_INSTALL_PREFIX "/Applications")
	SET(plugin_dest_dir MACMA.app/Contents/MacOS)
	SET(qtconf_dest_dir MACMA.app/Contents/Resources)
endif(APPLE)

if(UNIX AND NOT APPLE)
	SET(MAJOR_VERSION 1)
	SET(MINOR_VERSION 0)
	SET(PATCH_VERSION 0)
	
	INCLUDE(InstallRequiredSystemLibraries)
	
	SET(CPACK_BINARY_DEB ON)
	SET(CPACK_GENERATOR "DEB")

	SET(CPACK_SET_DESTDIR ON)
	SET(CPACK_INSTALL_PREFIX "/usr/local")
	SET(CMAKE_INSTALL_PREFIX "/usr/local/MACMA")
	SET(CPACK_PACKAGING_INSTALL_PREFIX "/usr/local")
	
	SET(CPACK_PACKAGE_DESCRIPTION "MACMA : Multiagent simulation of mantle convection")
	SET(CPACK_PACKAGE_DESCRIPTION_SUMMARY "MACMA : Multiagent simulation of mantle convection")
	SET(CPACK_PACKAGE_VENDOR "combes@enib.fr")
	SET(CPACK_PACKAGE_CONTACT "combes@enib.fr")
	SET(CPACK_PACKAGE_VERSION_MAJOR "${MAJOR_VERSION}")
	SET(CPACK_PACKAGE_VERSION_MINOR "${MINOR_VERSION}")
	SET(CPACK_PACKAGE_VERSION_PATCH "${PATCH_VERSION}")
	SET(CPACK_PACKAGE_FILE_NAME "${CMAKE_PROJECT_NAME}_${MAJOR_VERSION}.${MINOR_VERSION}.${CPACK_PACKAGE_VERSION_PATCH}")
	
	SET(CPACK_DEBIAN_PACKAGE_DEPENDS "libqt4-dev, libboost-thread-dev, libboost-iostreams-dev, gnuplot")
	
	SET(CPACK_DEBIAN_PACKAGE_PRIORITY "optional")
	SET(CPACK_DEBIAN_PACKAGE_SECTION "science")
	SET(CPACK_DEBIAN_ARCHITECTURE ${CMAKE_SYSTEM_PROCESSOR})
	
	SET(CPACK_COMPONENTS_ALL Libraries ApplicationData)
endif(UNIX AND NOT APPLE)

# --- Application
IF(APPLE)
    SET(APPS "\${CMAKE_INSTALL_PREFIX}/MACMA.app")
ELSE(APPLE)
    IF(WIN32)
        SET(APPS "\${CMAKE_INSTALL_PREFIX}/bin/MACMA.exe")
    ELSE(WIN32)
        SET(APPS "\${CMAKE_INSTALL_PREFIX}/MACMA")
    ENDIF(WIN32)
ENDIF(APPLE)

INSTALL(TARGETS MACMA 
    BUNDLE DESTINATION . COMPONENT Runtime
    RUNTIME DESTINATION bin COMPONENT Runtime
)

if(APPLE)
    INSTALL(DIRECTORY "${QT_PLUGINS_DIR}/imageformats" 
            DESTINATION ${plugin_dest_dir}/plugins 
            COMPONENT Runtime)

    INSTALL(CODE "
                 file(WRITE \"\${CMAKE_INSTALL_PREFIX}/${qtconf_dest_dir}/qt.conf\" \"\")
                 " 
            COMPONENT Runtime)

    SET(DIRS ${QT_LIBRARY_DIRS})

    INSTALL(CODE "
                 file(GLOB_RECURSE QTPLUGINS
                 	\"\${CMAKE_INSTALL_PREFIX}/${plugin_dest_dir}/plugins/*${CMAKE_SHARED_LIBRARY_SUFFIX}\"
		 )
                 include(BundleUtilities)
                 fixup_bundle(\"${APPS}\" \"\${QTPLUGINS}\" \"${DIRS}\")
                 " 
            COMPONENT Runtime)
endif(APPLE)

# --- Call CPack
include(CPack)
