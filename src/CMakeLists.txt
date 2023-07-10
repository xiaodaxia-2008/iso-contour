add_library(iso_contour)
target_sources(iso_contour PRIVATE iso_contour.cpp)
target_link_libraries(iso_contour PUBLIC CGAL::CGAL spdlog::spdlog)
target_include_directories(iso_contour PUBLIC ${CMAKE_CURRENT_LIST_DIR})

include(GNUInstallDirs)
install(TARGETS iso_contour EXPORT ${PROJECT_NAME}Targets
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

install(FILES "iso_contour.h" DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})