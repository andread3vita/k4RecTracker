/*
 * Copyright (c) 2020-2024 Key4hep-Project.
 *
 * This file is part of Key4hep.
 * See https://key4hep.github.io/key4hep-doc/ for further info.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
#include "ONNXRuntime.hpp"

#include "onnxruntime_cxx_api.h"
// #include "experimental_onnxruntime_cxx_api.h"

#include <fstream>
#include <iostream>

#include <algorithm>
#include <numeric>

ONNXRuntime::ONNXRuntime(const std::string& model_path, const std::vector<std::string>& input_names)
    : env_(new Ort::Env(OrtLoggingLevel::ORT_LOGGING_LEVEL_WARNING, "onnx_runtime")), allocator(),
      input_names_(input_names) {
  if (model_path.empty())
    throw std::runtime_error("Path to ONNX model cannot be empty!");
  Ort::SessionOptions options;
  options.SetIntraOpNumThreads(1);
  std::string model{model_path}; // fixes a poor Ort experimental API
  session_ = std::make_unique<Ort::Session>(*env_, model.c_str(), options);

  // Get input names and shapes
  input_node_strings_.clear();
  input_node_dims_.clear();

  for (size_t i = 0; i < session_->GetInputCount(); ++i) {
    // get input names
    const auto input_name =
        session_->GetInputNameAllocated(i, allocator).release(); // release the ownership of the pointer
    input_node_strings_.emplace_back(input_name);

    // get input shapes
    const auto nodeInfo = session_->GetInputTypeInfo(i);
    input_node_dims_[input_name] = nodeInfo.GetTensorTypeAndShapeInfo().GetShape();
  }

  // Get output names and shapes
  output_node_strings_.clear();
  output_node_dims_.clear();
  for (size_t i = 0; i < session_->GetOutputCount(); ++i) {
    // Get output names
    const auto output_name = session_->GetOutputNameAllocated(i, allocator).release();
    output_node_strings_.emplace_back(output_name);

    // get output shapes
    const auto nodeInfo = session_->GetOutputTypeInfo(i);
    output_node_dims_[output_name] = nodeInfo.GetTensorTypeAndShapeInfo().GetShape();

    // the 0th dim depends on the batch size
    output_node_dims_[output_name].at(0) = -1;
  }
}

ONNXRuntime::~ONNXRuntime() {}

template <typename T>
ONNXRuntime::Tensor<T> ONNXRuntime::run(Tensor<T>& input, const Tensor<long>& input_shapes,
                                        unsigned long long batch_size) const {
  std::vector<Ort::Value> tensors_in;
  for (const auto& name : input_node_strings_) {
    auto input_pos = variablePos(name);
    auto value = input.begin() + input_pos;
    std::vector<int64_t> input_dims;
    if (input_shapes.empty()) {
      input_dims = input_node_dims_.at(name);
      input_dims[0] = batch_size;
    } else {
      input_dims = input_shapes[input_pos];
    }
    // rely on the given input_shapes to set the batch size
    if (input_dims[0] != static_cast<long>(batch_size)) {
      throw std::runtime_error("The first element of `input_shapes` (" + std::to_string(input_dims[0]) +
                               ") does not match the given `batch_size` (" + std::to_string(batch_size) + ")");
    }
    auto expected_len = std::accumulate(input_dims.begin(), input_dims.end(), 1, std::multiplies<int64_t>());
    if (expected_len != (int64_t)value->size())
      throw std::runtime_error("Input array '" + name + "' has a wrong size of " + std::to_string(value->size()) +
                               ", expected " + std::to_string(expected_len));

    const OrtMemoryInfo* fInfo;
    fInfo = Ort::MemoryInfo::CreateCpu(OrtArenaAllocator, OrtMemTypeDefault);
    auto input_tensor =
        Ort::Value::CreateTensor<float>(fInfo, value->data(), value->size(), input_dims.data(), input_dims.size());
    if (!input_tensor.IsTensor())
      throw std::runtime_error("Failed to create an input tensor for variable '" + name + "'.");
    tensors_in.emplace_back(std::move(input_tensor));
  }

  // convert to char*
  std::vector<const char*> input_node_names;
  for (const auto& name_i : input_node_strings_) {
    input_node_names.push_back(name_i.c_str());
  }

  std::vector<const char*> output_node_names;
  for (const auto& name_j : output_node_strings_) {
    output_node_names.push_back(name_j.c_str());
  }

  // run the inference
  auto output_tensors = session_->Run(Ort::RunOptions{nullptr}, input_node_names.data(), tensors_in.data(),
                                      tensors_in.size(), output_node_names.data(), output_node_names.size());
  // convert output tensor to values
  Tensor<T> outputs;
  size_t i = 0;
  for (auto& output_tensor : output_tensors) {
    if (!output_tensor.IsTensor())
      throw std::runtime_error("(at least) inference output " + std::to_string(i) + " is not a tensor.");
    // get output shape
    auto tensor_info = output_tensor.GetTensorTypeAndShapeInfo();
    auto length = tensor_info.GetElementCount();

    auto floatarr = output_tensor.GetTensorData<float>();
    outputs.emplace_back(floatarr, floatarr + length);
    ++i;
  }
  if (outputs.size() != session_->GetOutputCount())
    throw std::runtime_error("Number of outputs differ from the expected one: got " + std::to_string(outputs.size()) +
                             ", expected " + std::to_string(session_->GetOutputCount()));

  return outputs;
}

size_t ONNXRuntime::variablePos(const std::string& name) const {
  auto iter = std::find(input_names_.begin(), input_names_.end(), name);
  if (iter == input_names_.end())
    throw std::runtime_error("Input variable '" + name + " is not provided");
  return iter - input_names_.begin();
}

template ONNXRuntime::Tensor<float> ONNXRuntime::run(Tensor<float>&, const Tensor<long>&, unsigned long long) const;
