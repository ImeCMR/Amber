/**
 * @author Yongxian Wu @ UC Irvine, March 2023
 * @brief refinedMLSES System with the LibTorch Library
 * @version 0.2
 *
 */

#include <algorithm>
#include <memory>
#include <numeric>
#include <tuple>
#include <vector>
#include <chrono>

#include "data.hpp"
#include "mlses.hpp"

// the input feature size of the model
constexpr unsigned input_feature_size = 96;


/**
 * @brief MLSES Forward Wrapper
 *  computes the MLSES model prediction based on custom input features
 *
 * @param isamples_nums the number of input samples
 * @param ibatch_size the number of samples for each batch
 * @param iworker_num the number of process for data loading
 * @param in the sorted input tensor where all samples are flattened into a 1-D array.
 * @param out the output allocated variable where the feature size of prediction output is 1
 */
void mlses_predict_wrapper_(int *isamples_nums, int *ibatch_size, int *iworker_num, const torch::Tensor& in, float *out){
    int64_t samples_nums = *isamples_nums;
    int64_t batch_size = *ibatch_size;
    int64_t worker_num = *iworker_num;

    // make array dataset
    mlses::MLSESArrayDataset dataset(in, input_feature_size, samples_nums);

    // make prediction
    auto prediction_tensor = mlses::compute_predict<mlses::MLSESArrayDataset>(dataset, batch_size, worker_num);

    // convert prediction to output variable
    auto output_feat_sz = prediction_tensor.size(1);
    for (int64_t i = 0; i < samples_nums; ++i) {
       for (int64_t j = 0; j < output_feat_sz; ++j) {
        // *(out + output_feat_sz * i + j) = prediction[i][j];
        *(out + output_feat_sz * i + j) = prediction_tensor.index({i, j}).item().to<float>();
       }
    }
}

/**
 * @brief MLSES Gradient Wrapper
 *  computes the gradient w.r.t the custom input features
 *
 * @param isamples_nums the number of input samples
 * @param ibatch_size the number of samples for each batch
 * @param iworker_num the number of process for data loading
 * @param in the sorted input tensor where all samples are flattened into a 1-D array.
 * @param out the output allocated variable where the feature size of prediction output is 200
 */
void mlses_gradient_wrapper_(int *isamples_nums, int *ibatch_size, int *iworker_num, const torch::Tensor& in, float *out) {
    int64_t samples_nums = *isamples_nums;
    int64_t batch_size = *ibatch_size;
    int64_t worker_num = *iworker_num;

    // make array dataset
    mlses::MLSESArrayDataset dataset(in, input_feature_size, samples_nums);

    // make gradient
    auto gradient_tensor = mlses::compute_gradient<mlses::MLSESArrayDataset>(dataset, batch_size, worker_num);

    // convert gradient to output variable
    auto output_feat_sz = gradient_tensor.size(1);
    for (int64_t i = 0; i < samples_nums; ++i) {
        for (int64_t j = 0; j < output_feat_sz; ++j) {
            *(out + output_feat_sz * i + j) = gradient_tensor.index({i, j}).item().to<float>();
        }
    }
}

/**
 * @brief MLSES Forward Wrapper
 *  computes the MLSES model prediction based on custom input features
 *
 * @param isamples_nums the number of input samples
 * @param in the sorted input tensor where all samples are flattened into a 1-D array.
 * @param out the output allocated variable where the feature size of prediction output is 1
 */
extern "C" void mlses_predict_on_tensor_wrapper_(int *isamples_nums, float *in, float *out){
    // log out device
    std::cout << mlses::device_info << std::endl;

    // prepare input
    int64_t samples_nums = *isamples_nums;
    auto tensor_opt = torch::TensorOptions().dtype(torch::kFloat32).device(torch::kCPU);
    auto input_tensor = torch::from_blob(in, {samples_nums, input_feature_size}, tensor_opt).to(mlses::device);

    // make prediction
    auto prediction_tensor = mlses::compute_predict_on_tensor(mlses::model, input_tensor);
    auto cpu_prediction_tensor = prediction_tensor.to(torch::kCPU);

    // convert prediction to output variable
    auto pred = cpu_prediction_tensor.data_ptr<float>();
    std::memcpy(out, pred, samples_nums * sizeof(float));
}

/**
 * @brief MLSES Backward Wrapper
 *  computes the MLSES model prediction based on custom input features
 *
 * @param isamples_nums the number of input samples
 * @param in the sorted input tensor where all samples are flattened into a 1-D array.
 * @param out the output allocated variable where the feature size of prediction output is 200
 */
extern "C" void mlses_gradient_on_tensor_wrapper_(int *isamples_nums, float *in, float *out){
    // log out device
    std::cout << mlses::device_info << std::endl;

    // prepare input
    int64_t samples_nums = *isamples_nums;
    auto tensor_opt = torch::TensorOptions().dtype(torch::kFloat32).device(torch::kCPU);
    auto input_tensor = torch::from_blob(in, {samples_nums, input_feature_size}, tensor_opt).to(mlses::device);

    // make prediction
    auto gradient_tensor = mlses::compute_gradient_on_tensor(mlses::model, input_tensor);
    auto cpu_gradient_tensor = gradient_tensor.to(torch::kCPU);

    // convert prediction to output variable
    auto grad = cpu_gradient_tensor.data_ptr<float>();
    std::memcpy(out, grad, samples_nums * input_feature_size * sizeof(float));
}

/**
 * @brief MLSES Preprocess Helper: sort input features based on atom distances
 *
 * @param isamples_nums the number of input samples
 * @param in the input features where all samples are flattened into a 1-D array. This is the distance vectors and
 * radi of atoms
 * @param key the input key (1-D) but 4 times as short as in() because it is the list of distances of atoms
 * @param padding_val the padding distance value indicates the atom is null
 * @return torch::Tensor the sorted input feature array in according device (CPU or GPU)
 */
torch::Tensor mlses_preprocess_helper_(int *isamples_nums, float *in, float *key, float *padding_val){
    size_t samples_nums = *isamples_nums;

    std::cout << "*****Following are debug information*****" << std::endl;
    std::cout << "\n\nOringial unsorted input features:" << std::endl;
    for (int i = 0; i < 100; ++i) {
        for (int j = 0; j < 200; ++j) {
            std::cout << *(in + 200 * i + j) << " ";
        }
        std::cout << std::endl;
    }
    std::cout << "\n\nKeys for sorting:" << std::endl;
    for (int i = 0; i < 100; ++i) {
        for (int j = 0; j < 50; ++j) {
            std::cout << *(key + 50 * i + j) << " ";
        }
        std::cout << std::endl;
    }

    // check whether cuda is avalible
    auto cuda_is_avalible = torch::cuda::is_available();

    // get sort index for each sample
    auto key_num = 50;
    int sample_num = *isamples_nums;
    auto key_tensor_opt = torch::TensorOptions().dtype(torch::kFloat32).device(torch::kCPU);
    auto key_tensor = torch::from_blob(key, {sample_num, key_num}, key_tensor_opt).to(cuda_is_avalible? torch::kCUDA: torch::kCPU);
    auto sort_result = torch::sort(key_tensor, -1);
    auto sort_value = std::get<0>(sort_result);
    auto sort_index = std::get<1>(sort_result);

    std::cout << "\n\nSorted key indexes:" << std::endl;
    std::cout << sort_index.index({torch::indexing::Slice(0, 100)}) << std::endl;

    // allocate tensor memory in accoding device
    auto sorted_in_tensor_opt = torch::TensorOptions().dtype(torch::kFloat32).device(cuda_is_avalible? torch::kCUDA: torch::kCPU);
    auto sorted_in = torch::zeros({sample_num * 200,}, sorted_in_tensor_opt);

    // sort feature vector based on given sort index
    for (int i = 0; i < sample_num; ++i) {
        for (int j = 0; j < key_num; ++j) {
            // check whether the atom is the padding one
            auto atom_val = sort_value.index({i, j}).item().to<float>();
            if (atom_val >= *padding_val) {
                break;
            }

            // write unnull atom grid (x, y, z, r) to the sorted feature array
            auto atom_ind = sort_index.index({i, j}).item().to<int>();
            for (int k = 0; k < 4; ++k) {
                sorted_in.index_put_({200 * i + 4 * j + k,}, in[200 * i + 4 * atom_ind + k]);
            }
        }
    }

    std::cout << "\n\nSorted input features:" << std::endl;
    std::cout << sorted_in.reshape({sample_num, 200}).index({torch::indexing::Slice(0, 100)}) << std::endl;

    return sorted_in;
}

/**
 * @brief MLSES Pipeline Wrapper
 *  computes the MLSES model prediction based on custom input features:
 *    (1) sort input features based on atom distances
 *    (2) call `mlses_predict_wrapper_`
 *    (3) call `mlses_gradient_wrapper_`
 *
 * @param isamples_nums the number of input samples
 * @param ibatch_size the number of samples for each batch
 * @param iworker_num the number of process for data loading
 * @param in the input feature pointer where all samples are flattened into a 1-D array.
 * @param key the input key (1-D) but 4 times as short as in() because it is the list of distances of atoms
 * @param out_feat the output allocated variable where the feature size of prediction output is 1
 */
extern "C" void mlses_preprocess_and_predict_wrapper_(int *isamples_nums, int *ibatch_size, int *iworker_num, float *in, float *key, float *out_feat){
    // sort feature
    float padding_value = 999.;
    auto sorted_in = mlses_preprocess_helper_(isamples_nums, in, key, &padding_value);

    // run forward
    mlses_predict_wrapper_(isamples_nums, ibatch_size, iworker_num, sorted_in, out_feat);
}