#include <string>
#include <Initializer/tree/LTSTree.hpp>
#include <Initializer/tree/Layer.hpp>

#include "hdf5.h"
#include <iostream>
#include <fstream>

#ifndef OUTPUT_HELPER_H_
#define OUTPUT_HELPER_H_


#if REAL_SIZE==8
#   define S3_H5_TYPE H5T_NATIVE_DOUBLE
#elif REAL_SIZE==4
#   define S3_H5_TYPE H5T_NATIVE_FLOAT
#endif

void write_dofs_to_file(seissol::initializers::LTSTree &ltsTree,
                        const seissol::initializers::Variable<real[tensor::Q::size()]> dofs_variable,
                        const seissol::initializers::LayerMask layer_type,
                        const std::string file_name) {

    hid_t file_id = H5Fcreate(file_name.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);


    unsigned cluster_counter  = 0;
    for (auto it = ltsTree.beginLeaf(); it != ltsTree.endLeaf(); ++it) {
        if (it->isMasked(layer_type)) {
            std::string cluster_name("/cluster=" + std::to_string(cluster_counter));

            hsize_t dims[3] = {it->getNumberOfCells(), tensor::Q::Shape[1], tensor::Q::Shape[0]};
            hid_t dataspace_id = H5Screate_simple(3, dims, NULL);
            hid_t dataset_id = H5Dcreate2(file_id,
                                          cluster_name.c_str(),
                                          S3_H5_TYPE,
                                          dataspace_id,
                                          H5P_DEFAULT,
                                          H5P_DEFAULT,
                                          H5P_DEFAULT);

            herr_t status = H5Dwrite(dataset_id,
                                     S3_H5_TYPE,
                                     H5S_ALL,
                                     H5S_ALL,
                                     H5P_DEFAULT,
                                     it->var(dofs_variable)[0]);

            status = H5Dclose(dataset_id);
            status = H5Sclose(dataspace_id);
        }
        ++cluster_counter;
    }

    hid_t status = H5Fclose(file_id);
}

void compare_dofs_with_file(seissol::initializers::LTSTree &ltsTree,
                            const seissol::initializers::Variable<real[tensor::Q::size()]> dofs_variable,
                            const seissol::initializers::LayerMask layer_type,
                            const std::string file_name) {

    std::ifstream file(file_name);
    if ((bool)file) {
        std::cout << "File for comparison has been found" << std::endl;
        std::cout << "Begin a comparison" << std::endl;
        hid_t file_id = H5Fopen(file_name.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

        unsigned cluster_counter = 0;
        for (auto it = ltsTree.beginLeaf(); it != ltsTree.endLeaf(); ++it) {
            if (it->isMasked(layer_type)) {
                std::string cluster_name("/cluster=" + std::to_string(cluster_counter));


                real *values_from_file = new real[it->getNumberOfCells() * tensor::Q::Shape[1] * tensor::Q::Shape[0]];

                hid_t dataset_id = H5Dopen2(file_id, cluster_name.c_str(), H5P_DEFAULT);
                herr_t status = H5Dread(dataset_id,
                                        S3_H5_TYPE,
                                        H5S_ALL,
                                        H5S_ALL,
                                        H5P_DEFAULT,
                                        values_from_file);

                status = H5Dclose(dataset_id);

                for (unsigned element = 0; element < it->getNumberOfCells(); ++element) {
                    unsigned offset = element * tensor::Q::Shape[1] * tensor::Q::Shape[0];
                    for (unsigned index = 0; index < tensor::Q::Shape[1] * tensor::Q::Shape[0]; ++index) {
                        const real difference = values_from_file[index + offset]
                                             - it->var(dofs_variable)[element][index];

                        const real eps = 1e-12;
                        if (std::fabs(difference) > eps) {
                            std::cout << "element:" << element << "|"
                                      << "index:" << index
                                      << std::endl;

                            std::cout << "must be: " << values_from_file[offset + index] << "|"
                                      << "computed: " << it->var(dofs_variable)[element][index] << "|"
                                      << "difference: " << difference
                                      << std::endl << std::endl;

                            hid_t status = H5Fclose(file_id);
                            delete[] values_from_file;
                            throw std::string("ERROR: comparison failed");
                        }
                    }
                }
                delete[] values_from_file;
            }
            ++cluster_counter;
        }

        hid_t status = H5Fclose(file_id);
        std::cout << "Everything is correct" << std::endl;
    }
    else {
        std::cout << "ERROR: cannot find a file for comparison" << std::endl;
    }
}

#endif  //OUTPUT_HELPER_H_