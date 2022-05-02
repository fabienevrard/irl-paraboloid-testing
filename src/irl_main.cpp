// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2022 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <cstdio>
#include <fstream>
#include <iomanip>
#include <iostream>

#include <nlohmann/json.hpp>

#include "src/irl_run.h"

static constexpr std::size_t MAX_INPUT_STRING_SIZE = 2000;

int main(int argc, char** argv) {
  if (argc < 2) {
    std::cout
        << "Must supply name for file to write results to as first argument.\n";
    std::exit(-1);
  }
  if (argc < 3) {
    std::cout
        << "Must supply list of results files generated from AMR executable"
        << std::endl;
    std::exit(-1);
  }

  std::vector<std::string> file_names(argc - 2);
  for (int i = 2; i < argc; ++i) {
    file_names[i - 2] = std::string(argv[i]);
  }

  auto out_file = fopen(argv[1], "w");
  for (const auto& file_name : file_names) {
    std::cout << "Confirming results from " << file_name << std::endl;
    std::fstream file(file_name, file.binary | file.in);

    std::size_t input_length;
    file.read(reinterpret_cast<char*>(&input_length), sizeof(std::size_t));
    std::string input_file;
    input_file.resize(input_length);
    file.read(input_file.data(), input_length);
    auto input = nlohmann::json::parse(input_file, nullptr, true, true);

    const auto test_name = input["test_name"].get<std::string>();
    // Functions should return L1, L2, Linf, and total time
    std::size_t total_cases = 0;
    std::array<double, 10> results;
    if (test_name == "face_only") {
    } else if (test_name == "basic_cube") {
    } else if (test_name == "translating_cube") {
      const auto side_length =
          input["TranslatingCube"]["side_length"].get<double>();
      const auto range_steps =
          input["TranslatingCube"]["steps"].get<std::vector<std::size_t>>();
      total_cases = 0;
      for (auto& steps : range_steps) {
        total_cases += steps;
      }
      results = confirmTranslatingCubeResult(side_length, total_cases, file);

    } else if (test_name == "parameter_sweep") {
      auto geometry = input["Sweep"]["geometry"].get<std::string>();

      if (input["Sweep"]["random"].get<bool>()) {
        total_cases = input["Sweep"]["number_of_tests"].get<std::size_t>();
      } else {
        const auto translation_steps =
            input["Sweep"]["translation_steps"].get<std::vector<std::size_t>>();
        const auto rotation_steps =
            input["Sweep"]["rotation_steps"].get<std::vector<std::size_t>>();
        const auto coefficient_steps =
            input["Sweep"]["coefficient_steps"].get<std::vector<std::size_t>>();
        total_cases = translation_steps[0] * translation_steps[1] *
                      translation_steps[2] * rotation_steps[0] *
                      rotation_steps[1] * rotation_steps[2] *
                      coefficient_steps[0] * coefficient_steps[1];
      }
      results = confirmParameterSweepResult(geometry, total_cases, file);
    } else {
      std::cout << "Unkown test type \"" + test_name + "\"" << std::endl;
      std::exit(-1);
    }
    fprintf(out_file, "%s %10i", file_name.c_str(),
            static_cast<int>(total_cases));
    for (const auto& entry : results) {
      fprintf(out_file, " %15.8E", entry);
    }
    fprintf(out_file, "\n");
    std::cout << '\n';
  }
  fclose(out_file);

  // std::ifstream input_file(argv[1]);
  // if (!input_file.good()) {
  //   std::cout << "Problem reading input file" << std::endl;
  //   std::exit(-1);
  // }
  // auto input = nlohmann::json::parse(input_file, nullptr, true, true);

  // const auto test_name = input["test_name"].get<std::string>();
  // const auto file_name = input["out_file_name"].get<std::string>();
  // std::fstream file(file_name, file.binary | file.trunc | file.out);

  // std::string json_string = input.dump();
  // file << json_string;

  // if (test_name == "face_only") {
  // } else if (test_name == "basic_cube") {
  // } else if (test_name == "translating_cube") {
  //   const auto side_length =
  //       input["TranslatingCube"]["side_length"].get<double>();
  //   const auto range_vector = input["TranslatingCube"]["ranges"]
  //                                 .get<std::vector<std::vector<double>>>();
  //   const auto range_steps =
  //       input["TranslatingCube"]["steps"].get<std::vector<std::size_t>>();
  //   if (range_vector.size() != range_steps.size()) {
  //     std::cout << "A number of steps must be given for each range in the "
  //                  "supplied \"ranges\" vector.\n";
  //     std::exit(-1);
  //   }
  //   for (auto& range : range_vector) {
  //     if (range.size() != 2) {
  //       std::cout << "Each provided range in the \"ranges\" vector must be "
  //                    "exactly two entries.\n";
  //       std::exit(-1);
  //     }
  //   }

  //   runTranslatingCube(side_length, range_vector, range_steps, file);

  // } else if (test_name == "parameter_sweep") {
  //   const auto geometry = input["Sweep"]["geometry"].get<std::string>();

  //   const auto translation =
  //       input["Sweep"]["translation"].get<std::vector<std::vector<double>>>();
  //   if (translation.size() != 2 || translation[0].size() != 3 ||
  //       translation[1].size() != 3) {
  //     std::cout << "[\"Sweep\"][\"translation\"] must specify a starting "
  //                  "and ending point as two vectors, such as "
  //                  "[[0.0,0.0,0.0],[1.0,0.0,1.0]] \n";
  //     std::exit(-1);
  //   }

  //   const auto rotation =
  //       input["Sweep"]["rotation"].get<std::vector<std::vector<double>>>();
  //   if (rotation.size() != 2 || rotation[0].size() != 3 ||
  //       rotation[1].size() != 3) {
  //     std::cout << "[\"Sweep\"][\"rotation\"] must specify a starting "
  //                  "and ending rotation (in radians) as two vectors, such as
  //                  "
  //                  "[[-3.14,-3.14,-3.14],[3.14,3.14,3.14]] \n";
  //     std::exit(-1);
  //   }

  //   const auto coefficient =
  //       input["Sweep"]["coefficient"].get<std::vector<std::vector<double>>>();
  //   if (coefficient.size() != 2 || coefficient[0].size() != 2 ||
  //       coefficient[1].size() != 2) {
  //     std::cout << "[\"Sweep\"][\"coefficient\"] must specify a starting "
  //                  "and ending rotation for the a and b value of the "
  //                  "paraboloid, such as "
  //                  "[[-5.0,-5.0],[5.0,5.0]] \n";
  //     std::exit(-1);
  //   }

  //   const auto random = input["Sweep"]["random"].get<bool>();
  //   if (random) {
  //     const auto fix_to_paraboloid =
  //         input["Sweep"]["fix_to_paraboloid"].get<bool>();
  //     const auto number_of_tests =
  //         input["Sweep"]["number_of_tests"].get<std::size_t>();
  //     runRandomSweep(geometry, translation, rotation, coefficient,
  //                    fix_to_paraboloid, number_of_tests, file);

  //   } else {
  //     const auto translation_steps =
  //         input["Sweep"]["translation_steps"].get<std::vector<std::size_t>>();
  //     const auto rotation_steps =
  //         input["Sweep"]["rotation_steps"].get<std::vector<std::size_t>>();
  //     const auto coefficient_steps =
  //         input["Sweep"]["coefficient_steps"].get<std::vector<std::size_t>>();

  //     if (translation_steps.size() != 3 || rotation_steps.size() != 3 ||
  //         coefficient_steps.size() != 2) {
  //       std::cout << "The step counts for translation, rotation, and "
  //                    "coefficient must be 3, 3, and 2, respectively.\n";
  //       std::exit(-1);
  //     }
  //     runOrganizedSweep(geometry, translation, translation_steps, rotation,
  //                       rotation_steps, coefficient, coefficient_steps,
  //                       file);
  //   }
  // } else {
  //   std::cout << "Unkown test type \"" + test_name + "\"" << std::endl;
  //   std::exit(-1);
  // }
}
