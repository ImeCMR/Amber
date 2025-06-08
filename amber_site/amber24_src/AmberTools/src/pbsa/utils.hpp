#pragma once

#include <fstream>
#include <vector>
#include <string>
#include <limits>
#include "torch/torch.h"
#include "torch/script.h"

namespace mlses {
    // read data from txt file
    std::vector<std::vector<float>> read_txt_file(
        std::string txt_file,
        int max_rows = std::numeric_limits<int>::max()
    );

    // save data to txt file
    void save_txt_file(std::string txt_file, const std::vector<std::vector<float>>& data);

    // batch vector to tensor
    torch::Tensor batch_vector_to_tensor(std::vector<std::vector<float>>&);

    // tensor to batch vector
    std::vector<std::vector<float>> tensor_to_batch_vector(const torch::Tensor&);

    // load jit model
    std::pair<torch::jit::script::Module, torch::Device> load_jit_model(const std::string);
}
