#pragma once

#include <cstddef>
#include <fstream>
#include <ios>
#include <iostream>
#include <string>
#include <vector>
#include "torch/data.h"


namespace mlses {
    // MLSES dataset with file input
    class MLSESDataset: public torch::data::Dataset<MLSESDataset> {
        public:
            MLSESDataset(std::string txt_file, size_t max_rows = std::numeric_limits<size_t>::max()): file_path(txt_file) {
                std::ifstream input_file(txt_file);
                if (!input_file) {
                    throw std::invalid_argument(txt_file + " is not a valid file path.");
                };

                // set max_rows
                max_rows = max_rows > 0? max_rows: std::numeric_limits<size_t>::max();

                // record line position for quick getting method
                std::string line;
                while (input_file && sample_positions.size() < max_rows) {
                    this->sample_positions.push_back(input_file.tellg());
                    std::getline(input_file, line);

                    // log out readding process
                    if (this->sample_positions.size() % 50000 == 0) {
                        std::cout << "Read " << this->sample_positions.size() << " lines\n";
                    }
                }
                // file reach the end therefore the last one is null
                if (!input_file) {
                    this->sample_positions.pop_back();
                }

                this->sample_nums = this->sample_positions.size();
                std::cout << "Build dataset with " << this->sample_nums << " samples\n";
            }

            // Override get() function to return tensor at location index
            torch::data::Example<> get(size_t index) override {
                // read line based on index
                auto pos = sample_positions.at(index);
                std::ifstream input_file(this->file_path);
                input_file.seekg(pos);

                // parse feat
                float x;
                std::string line;
                std::vector<float> feat;

                std::getline(input_file, line);
                std::istringstream ss(line);
                while (ss >> x) {
                    feat.push_back(x);
                }

                int feats_sz = feat.size();

                auto options = torch::TensorOptions().dtype(torch::kFloat);
                torch::Tensor feat_tensor = torch::from_blob(feat.data(), {feats_sz}, options);
                torch::Tensor label_tensor = torch::ones({1}); // padding label
                return {feat_tensor.clone(), label_tensor.clone()};
            }

            // Override size() function, return the length of data
            torch::optional<size_t> size() const override {
                return this->sample_nums;
            }

        private:
            size_t sample_nums = 0;
            std::string file_path = "";
            std::vector<std::streampos> sample_positions = {};
    };

    // MLSES dataset with array pointer input
    class MLSESArrayDataset: public torch::data::Dataset<MLSESArrayDataset> {
        public:
            MLSESArrayDataset(const torch::Tensor& features, int feature_size, size_t sample_nums): sample_nums(sample_nums), feat_sz(feature_size), feature_tensor(features) {
                std::cout << "Build dataset with " << this->sample_nums << " samples\n";

                // create dummpy label tensor during initialization
                auto label_tensor_opt = torch::TensorOptions().device(feature_tensor.device());
                label_tensor = torch::ones({1}, label_tensor_opt);
            }

            // Override get() function to return tensor at location index
            torch::data::Example<> get(size_t index) override {
                // get pointer location
                auto start_index = this->feat_sz * index;
                auto end_index = this->feat_sz * (index + 1);

                // slice tensor
                auto feat_tensor = feature_tensor.index({torch::indexing::Slice(start_index, end_index)});
                return {feat_tensor, label_tensor};
            }

            // Override size() function, return the length of data
            torch::optional<size_t> size() const override {
                return this->sample_nums;
            }

        private:
            size_t sample_nums = 0;
            const int feat_sz = 0;
            const torch::Tensor& feature_tensor;
            torch::Tensor label_tensor;
    };
}