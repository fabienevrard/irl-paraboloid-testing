set(CURRENT_DIR "${CMAKE_CURRENT_SOURCE_DIR}")
target_sources(amr_generate_result PRIVATE "${CURRENT_DIR}/amr_main.cpp")
target_sources(amr_generate_result PRIVATE "${CURRENT_DIR}/amr_run.cpp")
target_sources(amr_generate_result PRIVATE "${CURRENT_DIR}/amr_run.h")


target_sources(irl_confirm_result PRIVATE "${CURRENT_DIR}/irl_main.cpp")
target_sources(irl_confirm_result PRIVATE "${CURRENT_DIR}/irl_run.h")
target_sources(irl_confirm_result PRIVATE "${CURRENT_DIR}/irl_run.tpp")
target_sources(irl_confirm_result PRIVATE "${CURRENT_DIR}/my_dodecahedron.h")
target_sources(irl_confirm_result PRIVATE "${CURRENT_DIR}/my_dodecahedron.tpp")


