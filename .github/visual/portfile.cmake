vcpkg_minimum_required(VERSION 2022-10-12) # for ${VERSION}

vcpkg_from_github(
	OUT_SOURCE_PATH 	SOURCE_PATH 
	REPO	P-p-H-d/mpfr
	HEAD_REF 	master
	PATCHES	dll.patch	tst.patch
)
 
vcpkg_configure_make(
    SOURCE_PATH "${SOURCE_PATH}"
    AUTOCONFIG
)

vcpkg_build_make(
	BUILD_TARGET check
)

file(REMOVE_RECURSE "${CURRENT_PACKAGES_DIR}/debug/share")
file(REMOVE_RECURSE "${CURRENT_PACKAGES_DIR}/debug/include")

file(REMOVE
    "${CURRENT_PACKAGES_DIR}/share/${PORT}/AUTHORS"
    "${CURRENT_PACKAGES_DIR}/share/${PORT}/BUGS"
    "${CURRENT_PACKAGES_DIR}/share/${PORT}/COPYING"
    "${CURRENT_PACKAGES_DIR}/share/${PORT}/COPYING.LESSER"
    "${CURRENT_PACKAGES_DIR}/share/${PORT}/NEWS"
    "${CURRENT_PACKAGES_DIR}/share/${PORT}/TODO"
)

file(INSTALL "${CMAKE_CURRENT_LIST_DIR}/usage" DESTINATION "${CURRENT_PACKAGES_DIR}/share/${PORT}")
vcpkg_install_copyright(FILE_LIST "${SOURCE_PATH}/COPYING" "${SOURCE_PATH}/COPYING.LESSER")
