#pragma once

#include <cstddef>
#include <utility>

#include "torch/torch.h"
#include "torch/script.h"
#include "data.hpp"


namespace mlses {
    // const jit model path
    const std::string amber_home = std::string(std::getenv("AMBERHOME"));
    const static std::string pretrained_mlses_model_path = amber_home + "/lib/mlses_traced.pt";

    // load jit traced model
    static bool cuda_available = torch::cuda::is_available();
    static torch::Device device = cuda_available? torch::kCUDA: torch::kCPU;
    static torch::jit::script::Module model = torch::jit::load(pretrained_mlses_model_path, device);
    static std::string device_info = cuda_available ? "CUDA available. MLSES running on GPU.": "CUDA unavailable. MLSES running on CPU.";

    // compute predicted value based on torch model
    template<typename Dataset = mlses::MLSESDataset>
    torch::Tensor compute_predict(
        Dataset&,
        int64_t = 128,
        int64_t = 4
    );

    // compute gradient based on torch model
    template<typename Dataset = mlses::MLSESDataset>
    torch::Tensor compute_gradient(
        Dataset&,
        int64_t = 128,
        int64_t = 4
    );

    // compute predicted value on one tensor
    torch::Tensor compute_predict_on_tensor(torch::jit::script::Module&, const torch::Tensor&);

    // compute gradient value on one tensor
    torch::Tensor compute_gradient_on_tensor(torch::jit::script::Module&, const torch::Tensor&);
}
