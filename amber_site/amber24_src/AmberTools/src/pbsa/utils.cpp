/**
 * @author Yongxian Wu @ UC Irvine, March 2023
 * @brief refinedMLSES System with the LibTorch Library
 * @version 0.2
 *
 */

#include "utils.hpp"
#include <sstream>
#include <stdexcept>


namespace mlses {
    std::vector<std::vector<float>> read_txt_file(std::string txt_file, int max_rows) {
        std::ifstream input_file(txt_file);
        if (!input_file) {
            throw std::invalid_argument(txt_file + " is not a valid file path.");
        };

        int line_idx = 0;
        std::vector<std::vector<float>> features;
        std::string line;
        while (std::getline(input_file, line) && line_idx < max_rows) {
            std::istringstream ss(line);
            features.push_back({});

            float x;
            while (ss >> x) {
                features.back().push_back(x);
            }

            ++line_idx;
        }

        return features;
    }


    void save_txt_file(std::string txt_file, const std::vector<std::vector<float>>& data) {
        std::ofstream output_file(txt_file);
        if(!output_file) {
            throw std::invalid_argument(txt_file + " is not a valid file path.");
        }

        for (auto iter = data.begin(); iter != data.end(); ++iter) {
            std::ostringstream ss;
            for (auto line_iter = iter->begin(); line_iter != iter->end(); ++line_iter) {
                if (line_iter + 1 == iter->end()) {
                    ss << *line_iter << "\n";
                } else {
                    ss << *line_iter << " ";
                }
            }
            output_file << ss.str();
        }
    }

    // batch vector to tensor
    torch::Tensor batch_vector_to_tensor(std::vector<std::vector<float>>& input) {
        // prepare input feats to input tensor
        int batch_sz = input.size();
        int feats_sz = input[0].size();

        auto options = torch::TensorOptions().dtype(torch::kFloat);
        torch::Tensor tensor = torch::zeros({batch_sz, feats_sz}, options);
        for (int i = 0; i < batch_sz; ++i) {
            tensor.slice(0, i, i + 1) = torch::from_blob(
                input[i].data(), {feats_sz}, options
            );
        }

        return tensor;
    }

    // tensor to batch vector
    std::vector<std::vector<float>> tensor_to_batch_vector(const torch::Tensor& tensor) {
        int batch_sz = tensor.size(0);
        std::vector<std::vector<float>> value(batch_sz);
        for (int i = 0; i < batch_sz; ++i) {
            torch::Tensor tensor_slice = tensor.slice(0, i, i + 1);
            value[i] = std::vector<float>(
                tensor_slice.data_ptr<float>(),
                tensor_slice.data_ptr<float>() + tensor_slice.numel()
            );
        }

        return value;
    }

    std::pair<torch::jit::script::Module, torch::Device> load_jit_model(const std::string jit_path) {
        torch::jit::script::Module model;

        // load pretrained model
        try {
            model = torch::jit::load(jit_path);
        }
        catch(const torch::Error& error) {
            std::cerr << "Could not load pretrained scriptmodule from the file path " << jit_path << " ." << "Please check the file path.\n";
        }

        // check cuda context
        auto cuda_available = torch::cuda::is_available();
        torch::Device device = cuda_available? torch::kCUDA: torch::kCPU;
        model.to(device);

        // log out device
        std::cout << (cuda_available ? "CUDA available. MLSES running on GPU." : "CUDA unavailable. MLSES running on CPU.") << std::endl;

        return {model, device};
    }
}
