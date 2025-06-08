/**
 * @author Yongxian Wu @ UC Irvine, March 2023
 * @brief refinedMLSES System with the LibTorch Library
 * @version 0.2
 *
 */

#include <c10/core/Device.h>
#include <c10/core/DeviceType.h>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <torch/data.h>
#include <torch/data/dataloader.h>
#include <torch/data/samplers/sequential.h>
#include <utility>
#include <cstdlib>

#include "mlses.hpp"
#include "data.hpp"
#include "utils.hpp"

namespace mlses {
    template<typename Dataset>
    torch::Tensor compute_predict(
        Dataset& dataset,
        int64_t batch_size,
        int64_t dataloader_worker
    ) {
        auto jit = load_jit_model(mlses::pretrained_mlses_model_path);
        auto model = jit.first;
        auto device = jit.second;
        int64_t sample_num = dataset.size().value();

        // prepare dataset and dataloader
        int64_t inference_batch_size = batch_size > 0? batch_size: sample_num;
        auto inference_dataset = dataset.map(torch::data::transforms::Stack<>());
        auto dataloader_opt = torch::data::DataLoaderOptions(inference_batch_size).workers(dataloader_worker);
        auto inference_dataloader = torch::data::make_data_loader<torch::data::samplers::SequentialSampler>(std::move(inference_dataset), dataloader_opt);

        // forward without taping gradient
        torch::NoGradGuard no_grad;
        model.eval();
        std::vector<torch::Tensor> batch_predict;
        for (auto& batch_x : *inference_dataloader) {
            auto input_batch_x = batch_x.data;
            std::vector<torch::jit::IValue> inputs{input_batch_x};
            torch::Tensor batch_y = model.forward(inputs).toTensor() - 1.5;

            // post-processing tensor to vector
            // std::vector<std::vector<float>> batch_predict = tensor_to_batch_vector(batch_y);
            // predict.insert(predict.end(), batch_predict.begin(), batch_predict.end());
            batch_predict.push_back(batch_y.detach());
        }

        // concat final result
        torch::TensorList list(batch_predict);
        auto predict = torch::cat(list, 0);
        return predict;
    }

    template<typename Dataset>
    torch::Tensor compute_gradient(
        Dataset& dataset,
        int64_t batch_size,
        int64_t dataloader_worker
    ) {
        auto jit = load_jit_model(mlses::pretrained_mlses_model_path);
        auto model = jit.first;
        auto device = jit.second;
        int64_t sample_num = dataset.size().value();

        // prepare dataset and dataloader
        int64_t inference_batch_size = batch_size > 0? batch_size: sample_num;
        auto inference_dataset = dataset.map(torch::data::transforms::Stack<>());
        auto dataloader_opt = torch::data::DataLoaderOptions(inference_batch_size).workers(dataloader_worker);
        auto inference_dataloader = torch::data::make_data_loader<torch::data::samplers::SequentialSampler>(std::move(inference_dataset), dataloader_opt);

        // forward and backward
        std::vector<torch::Tensor> batch_gradient;
        for (auto& batch_x : *inference_dataloader) {
            torch::Tensor batch_grad_x = batch_x.data.requires_grad_();
            batch_grad_x.retain_grad();
            std::vector<torch::jit::IValue> inputs{batch_grad_x};
            torch::Tensor batch_y = model.forward(inputs).toTensor();
            torch::Tensor gradient_output = torch::ones_like(batch_y);
            batch_y.backward(gradient_output);

            // post-porcessing tensor to vector
            // std::vector<std::vector<float>> batch_gradient = tensor_to_batch_vector(batch_grad_x.grad().to(torch::kCPU));
            // gradient.insert(gradient.end(), batch_gradient.begin(), batch_gradient.end());
            batch_gradient.push_back(batch_grad_x.grad().detach().clone());

            // zero model gradient
            batch_grad_x.grad().zero_();
            for (const auto& parameter: model.named_parameters()) {
                if (parameter.value.requires_grad()) {
                    parameter.value.grad().zero_();
                }
            }
        }

        // concat final result
        torch::TensorList list(batch_gradient);
        auto gradient = torch::cat(list, 0);
        return gradient;
    }

    // instance functions
    template torch::Tensor compute_predict<mlses::MLSESDataset>(
        mlses::MLSESDataset&,
        int64_t,
        int64_t
    );

    template torch::Tensor compute_predict<mlses::MLSESArrayDataset>(
        mlses::MLSESArrayDataset&,
        int64_t,
        int64_t
    );

    template torch::Tensor compute_gradient<mlses::MLSESDataset>(
        mlses::MLSESDataset&,
        int64_t,
        int64_t
    );

    template torch::Tensor compute_gradient<mlses::MLSESArrayDataset>(
        mlses::MLSESArrayDataset&,
        int64_t,
        int64_t
    );

    torch::Tensor compute_predict_on_tensor(
        torch::jit::script::Module& model,
        const torch::Tensor& input_feature
    ) {
        torch::InferenceMode guard; // runnning on inference mode
        // compute forward
        std::vector<torch::jit::IValue> inputs{std::move(input_feature)};
        torch::Tensor output = model.forward(inputs).toTensor() - 1.5;
        return output;
    }

    torch::Tensor compute_gradient_on_tensor(
        torch::jit::script::Module& model,
        const torch::Tensor& input_feature
    ) {
        auto gradable_input_feature = input_feature.requires_grad_();
        gradable_input_feature.retain_grad();
        std::vector<torch::jit::IValue> inputs{gradable_input_feature};
        torch::Tensor output = model.forward(inputs).toTensor();
        torch::Tensor gradient_output = torch::ones_like(output);
        output.backward(gradient_output);
        return gradable_input_feature.grad();
    }
}
