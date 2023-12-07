#include "csr.h"
#include "file_name.h"
#include "libsvm.h"
#include "parse_param.h"

#include <chrono>
#include <iostream>
#include <random>
#include <set>
#include <vector>

namespace {

    struct data_t {
        std::size_t value;
        std::size_t index_one;
        std::size_t index_two;

        bool operator < (data_t const & rhs) const {
            return value < rhs.value; // ascending order!
        }
    };

    file_name_c * fn_data{nullptr};
    csr_c data;
    std::vector<double> model;

    std::size_t random_seed;
    std::size_t sample_size;
    std::size_t sets;

    bool no_refine_step{false};

    double rho{0.1};
    double tau{0.8};

    std::mt19937 mt_random;

    std::vector<std::size_t> transform_row_idx;
    std::vector<std::size_t> transform_col_idx;

    std::vector<std::vector<std::size_t>> clusters;
    std::vector<std::size_t> noise;

/*
 * psKC
 */

    void compute_psi_D(std::vector<double> & psi_D, std::set<std::size_t> const & D) {
        psi_D.clear();
        psi_D.resize(sets * sample_size, 0.0);

        for (std::size_t x : D) {
            std::size_t col_s = transform_row_idx[x];
            std::size_t col_e = transform_row_idx[x + 1];

            for (std::size_t idx2 = col_s; idx2 < col_e; ++idx2) {
                psi_D[transform_col_idx[idx2]]++;
            }
        }

        for (std::size_t i = 0; i < psi_D.size(); ++i) {
            psi_D[i] /= D.size();
        }
    }

    std::size_t find_max(std::vector<double> const & psi_D, std::set<std::size_t> const & D) {
        std::size_t x_p = 0;
        double best_sim = 0.0;

        for (std::size_t x : D) {
            std::size_t col_s = transform_row_idx[x];
            std::size_t col_e = transform_row_idx[x + 1];

            double sim = 0.0;

            for (std::size_t c_idx = col_s; c_idx < col_e; c_idx++) {
            sim += psi_D[transform_col_idx[c_idx]];
            }

            sim /= sets;

            if (sim > best_sim) {
                best_sim = sim;
                x_p = x;
            }
        }

        return x_p;
    }

    double find_max(std::size_t const x_p, std::size_t & x_q, std::set<std::size_t> const & D) {
        double best_sim = 0.0;
        x_q = 0;

        std::size_t p_col_s = transform_row_idx[x_p];
        std::size_t p_col_e = transform_row_idx[x_p + 1];

        for (std::size_t x : D) {
            std::size_t col_s = transform_row_idx[x];
            std::size_t col_e = transform_row_idx[x + 1];

            std::size_t p_col = p_col_s;
            double sim = 0.0;

            while ((p_col < p_col_e) && (col_s < col_e)) {
                if (transform_col_idx[p_col] == transform_col_idx[col_s]) {
                    sim++;
                    p_col++;
                    col_s++;
                } else if (transform_col_idx[p_col] < transform_col_idx[col_s]) {
                    p_col++;
                } else {
                    col_s++;
                }
            }

            if (sim > best_sim) {
                best_sim = sim;
                x_q = x;
            }
        }

        return best_sim / sets;
    }

    void compute_set(std::vector<std::size_t> const & set, std::vector<double> & pt) {
        pt.clear();
        pt.resize(sets * sample_size, 0.0);

        for (std::size_t idx : set) {
            std::size_t col_s = transform_row_idx[idx];
            std::size_t col_e = transform_row_idx[idx + 1];

            for (std::size_t idx2 = col_s; idx2 < col_e; idx2++) {
                pt[transform_col_idx[idx2]]++;
            }
        }

        for (std::size_t i = 0; i < pt.size(); i++) {
            pt[i] /= set.size();
        }
    }

    double compute_score(std::size_t idx, std::vector<double> const & pt) {
        std::size_t col_s = transform_row_idx[idx];
        std::size_t col_e = transform_row_idx[idx + 1];

        double sim = 0.0;

        for (std::size_t c_idx = col_s; c_idx < col_e; c_idx++) {
            sim += pt[transform_col_idx[c_idx]];
        }

        return sim / sets;
    }

    void find_min_score(std::size_t & min_idx, double & min_score, std::vector<std::size_t> const & G_k, std::vector<double> const & pt) {
        min_idx = 0;
        min_score = 10.0;

        for (std::size_t x : G_k) {
            double score = compute_score(x, pt);

            if (score < min_score) {
                min_score = score;
                min_idx = x;
            }
        }
    }

    void refine_step() {
        std::vector<std::vector<double>> G_pt;

        for (std::vector<std::size_t> & G_k : clusters) {
          std::vector<double> pt(sets * sample_size, 0.0);
          compute_set(G_k, pt);
          G_pt.push_back(pt);
        }

        for (std::size_t i = 0; i < clusters.size(); ++i) {
            std::vector<std::size_t> & G_k_i = clusters[i];
            std::vector<double> & pt_i = G_pt[i];

            bool finished = false;

            while (!finished) {
                finished = true;
                std::size_t min_idx;
                double min_score;

                find_min_score(min_idx, min_score, G_k_i, pt_i);

                std::size_t best_idx = 0;

                for (std::size_t j = 0; j < clusters.size(); ++j) {
                    if (i != j) {
                        double score = compute_score(min_idx, G_pt[j]);

                        if (score > min_score) {
                            min_score = score;
                            best_idx = j;
                            finished = false;
                        }
                    }
                }

                if (!finished) {
                    std::vector<std::size_t> & G_k_j = clusters[best_idx];
                    std::vector<double> & pt_j = G_pt[best_idx];

                    G_k_i.erase(std::remove(G_k_i.begin(), G_k_i.end(), min_idx), G_k_i.end());
                    G_k_j.push_back(min_idx);

                    compute_set(G_k_i, pt_i);
                    compute_set(G_k_j, pt_j);
                }
            }
        }
    }

    void assign() {
        std::vector<std::size_t> D_temp(data.rows());
        std::iota(D_temp.begin(), D_temp.end(), 0);
        std::set<std::size_t> D(D_temp.begin(), D_temp.end());

        std::vector<double> psi_D;

        while (D.size() > 1) {
            compute_psi_D(psi_D, D);

            std::size_t x_p = find_max(psi_D, D);
            D.erase(x_p);
            std::size_t x_q;
            double sim_pq = find_max(x_p, x_q, D);

            double gamma = (1.0 - rho) * sim_pq;

            if (gamma <= tau) {
                break;
            }

            std::vector<std::size_t> G_k;
            G_k.push_back(x_p);
            G_k.push_back(x_q);

            D.erase(x_q);

            while (gamma > tau) {
                std::vector<std::size_t> S;

                std::vector<double> pt(sets * sample_size, 0.0);

                for (std::size_t idx : G_k) {
                    std::size_t col_s = transform_row_idx[idx];
                    std::size_t col_e = transform_row_idx[idx + 1];

                    for (std::size_t idx2 = col_s; idx2 < col_e; idx2++) {
                        pt[transform_col_idx[idx2]]++;
                    }
                }

                for (std::size_t i = 0; i < pt.size(); ++i) {
                    pt[i] /= G_k.size();
                }

                for (std::size_t x : D) {
                    std::size_t col_s = transform_row_idx[x];
                    std::size_t col_e = transform_row_idx[x + 1];

                    double sim = 0.0;

                    for (std::size_t c_idx = col_s; c_idx < col_e; ++c_idx) {
                        sim += pt[transform_col_idx[c_idx]];
                    }

                    sim /= sets;

                    if (sim > gamma) {
                        S.push_back(x);
                    }
                }

                G_k.insert(G_k.end(), S.begin(), S.end());

                for (std::size_t idx : S) {
                    D.erase(idx);
                }

                gamma *= 1.0 - rho;
            }

            clusters.push_back(G_k);
        }

        noise.insert(noise.end(), D.begin(), D.end());

        if (!no_refine_step) {
            refine_step();
        }
    }

/*
 * Auxiliary codes
 */

    double aNNE_distance_pt(std::size_t idx_x, std::size_t idx_y) {
        double sum = 0.0;

        std::size_t col_s = data.r_idx(idx_x);
        std::size_t col_e = data.r_idx(idx_x + 1);

        std::size_t y = idx_y * data.cols();

        for (std::size_t i = 0, col_idx = col_s; i < data.cols(); ++i) {
            if ((col_idx < col_e) && (i == data.c_idx(col_idx))) {
#ifdef NDEBUG
                double val = data.value(col_idx) - model[y + i];
#else
                double val = data.value(col_idx) - model.at(y + i);
#endif
                sum += val * val;
                col_idx++;
            } else {
#ifdef NDEBUG
                sum += model[y + i] * model[y + i];
#else
                sum += model.at(y + i) * model.at(y + i);
#endif
            }
        }

        return sum;
    }

    void aNNE_transform_sparse() {
    // CSR format without storing the values
        transform_row_idx.push_back(0);

        for (std::size_t i = 0, r_idx = 0; i < data.rows(); ++i) {
            for (std::size_t j = 0; j < sets; ++j) {
                std::size_t idx_t = j * sample_size;

                double best_dist = std::numeric_limits<double>::max();
                std::size_t best_idx = 0;

                for (std::size_t k = 0; k < sample_size; ++k) {
                    std::size_t idx = idx_t + k;

                    double dist = aNNE_distance_pt(i, idx);

                    if (dist < best_dist) {
                        best_dist = dist;
                        best_idx = k;
                    }
                }

                transform_col_idx.push_back(idx_t + best_idx);
                r_idx++;
            }

            transform_row_idx.push_back(r_idx);
        }
    }

    void sub_sampling(std::vector<std::size_t> & list, std::uniform_int_distribution<> & uid, std::mt19937 & mt_random, std::set<std::size_t> & sampling_without_replacement, std::size_t max_rows, std::size_t psi) {
        list.clear();

        std::size_t min_psi = std::min(psi, max_rows);

        for (std::size_t i = 0; i < min_psi; ++i) {
            if (sampling_without_replacement.size() >= max_rows) {
                sampling_without_replacement.clear();
                sampling_without_replacement.insert(list.begin(), list.end());
            }

            std::size_t random_index = std::size_t(uid(mt_random));

            while (sampling_without_replacement.find(random_index) != sampling_without_replacement.end()) {
                random_index = std::size_t(uid(mt_random));
            }

            list.push_back(random_index);
            sampling_without_replacement.insert(random_index);
        }
    }

    void clean_up() {
        delete fn_data;
    }

    void init() {
        mt_random.seed(random_seed);
        std::uniform_int_distribution<> uid = std::uniform_int_distribution<>(0, int(data.rows() - 1));

        std::set<std::size_t> sampling_without_replacement;
        std::vector<std::size_t> new_idx;

        for (std::size_t i = 0; i < sets; ++i) {
            sub_sampling(new_idx, uid, mt_random, sampling_without_replacement, data.rows(), sample_size);

            for (uint64_t val : new_idx) {
                std::size_t col_s = data.r_idx(val);
                std::size_t col_e = data.r_idx(val + 1);

                for (std::size_t j = 0, col_idx = col_s; j < data.cols(); ++j) {
                    if ((col_idx < col_e) && (j == data.c_idx(col_idx))) {
                        model.push_back(data.value(col_idx));
                        col_idx++;
                    } else {
                        model.push_back(0.0);
                    }
                }
            }
        }

        aNNE_transform_sparse();

        std::cout << "psi  : " << sample_size << std::endl;
        std::cout << "sets : " << sets << std::endl;
        std::cout << "# idx: " << model.size() << std::endl;
        std::cout << std::endl;
    }

    void load_data() {
        std::cout << "Loading " << fn_data->get_file_name() << "..." << std::endl;

        std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();

        libsvm_c * libsvm_data;

        libsvm_data = new libsvm_c(*fn_data);

        libsvm_data->load(data);
        delete libsvm_data;

        sample_size = std::min(sample_size, data.rows());

        std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
        std::chrono::duration<double> elapsed = end - start;
        double time = elapsed.count();
        std::cout << "Loading time: " << time << " secs" << std::endl;
        std::cout << std::endl;

        std::cout << "Rows: " << data.rows() << std::endl;
        std::cout << "Cols: " << data.cols() << std::endl;
        std::cout << std::endl;
    }

    void parse_command_line(int argc, char ** argv) {
        {
            std::string data_file_name;
            parse_param(argc, argv, "--data", data_file_name, "No data file name given (--data).");
            fn_data = new file_name_c(data_file_name);
        }

        parse_param(argc, argv, "--random_seed", random_seed, "No random seed given (--random_seed).");
        parse_param(argc, argv, "--sample_size", sample_size, "No sample size given (--sample_size).");
        parse_param(argc, argv, "--sets", sets, "No sets given (--sets).");

        parse_param(argc, argv, "--growth_rate", rho, "No growth rate given (--growth_rate).");
        parse_param(argc, argv, "--threshold", tau, "No threshold given (--threshold).");

        parse_param(argc, argv, "--no_refine_step", no_refine_step);

        std::cout << "Random seed: " << random_seed << std::endl;
        std::cout << "Growth rate: " << rho << std::endl;
        std::cout << "Threshold  : " << tau << std::endl;
        std::cout << std::endl;
    }

/*
 * Computation of F1 measure
 */

    double compute_precision(std::vector<std::size_t> const & confusion_matrix, std::size_t idx_cls, std::size_t idx_clu, std::size_t num_classes, std::size_t num_clusters_detected) {
        double correct = 0.0;
        double total = 0.0;

        for (std::size_t i = 0; i < num_classes; ++i) {
            if (i == idx_cls) {
                correct += confusion_matrix[i * num_clusters_detected + idx_clu];
            }

            total += confusion_matrix[i * num_clusters_detected + idx_clu];
        }

        return correct / total;
    }

    double compute_recall(std::vector<std::size_t> const & confusion_matrix, std::size_t idx_cls, std::size_t idx_clu, std::size_t num_clusters_detected) {
        double correct = 0.0;
        double total = 0.0;

        for (std::size_t i = 0; i < num_clusters_detected; ++i) {
            if (i == idx_clu) {
                correct += confusion_matrix[idx_cls * num_clusters_detected + i];
            }

            total += confusion_matrix[idx_cls * num_clusters_detected + i];
        }


        return correct / total;
    }

    void quick_mapping(bool noise_found, std::size_t num_clusters, std::size_t num_classes, std::vector<std::size_t> & counts, std::vector<std::size_t> & total, std::vector<double> & best) {
        std::vector<data_t> list;

        std::size_t count = noise_found ? num_clusters - 1 : num_clusters;

        for (std::size_t i = 0; i < count; ++i) {
            best[i] = -1;

            for (std::size_t j = 0; j < num_classes; j++) {
                data_t data{counts[j * num_clusters + i], i, j};
                list.push_back(data);
            }
       }

       std::sort(list.begin(), list.end());
       std::reverse(list.begin(), list.end()); // descending order

       for (data_t & data : list) {
            if (data.value == 0) {
                break;
            }

            if (best[data.index_one] == -1) {
                bool found = false;

                for (std::size_t i = 0; i < num_clusters; i++) {
                    if (best[i] == data.index_two) {
                        found = true;
                        break;
                    }
                }

                if (!found) {
                    best[data.index_one] = data.index_two;
                }
            }
        }
    }

    void compute_F1_measure() {
        std::size_t num_clusters_detected = clusters.size();
        std::size_t num_classes = data.num_labels();

        if (clusters.size() == 0) {
            std::cout << "No clusters detected!" << std::endl;
            return;
        }

        std::cout << "# Classes       : " << num_classes << std::endl;
        std::cout << "# Clusters found: " << num_clusters_detected << std::endl;

        bool noise_found = false;

        if (noise.size() > 0) {
            num_clusters_detected++;
            noise_found = true;
        }

        std::cout << "Noise points    : " << (noise_found ? "yes" : "no") << std::endl;

        if (noise_found) {
            std::cout << "# Not assigned  : " << noise.size() << std::endl;
        }

        std::cout << std::endl;

    // confustion[actual][predicted]
        std::vector<std::size_t> confusion_matrix(num_clusters_detected * num_classes, 0);
        std::vector<std::size_t> cluster_total(num_clusters_detected, 0);

        for (std::size_t i = 0; i < num_clusters_detected; ++i) {
            std::vector<std::size_t> const & next = (noise_found && (i == (num_clusters_detected - 1))) ? noise : clusters[i];

            for (std::size_t idx : next) {
                std::size_t cls_idx = data.get_label_idx(idx);

                confusion_matrix[cls_idx * num_clusters_detected + i]++;
                cluster_total[i]++;
            }
        }

        std::cout << "Confusion Matrix:" << std::endl;
        int spc_wth = 8;

        std::cout << std::setw(spc_wth) << " ";

        for (std::size_t i = 0; i < num_clusters_detected; ++i) {
            if (noise_found && (i == (num_clusters_detected - 1))) {
                std::cout << std::setw(spc_wth) << "N";
            } else {
                std::cout << std::setw(spc_wth) << i;
            }
        }

        std::cout << std::setw(spc_wth) << "Tot" << std::endl;

        for (std::size_t i = 0; i < num_classes; ++i) {
            std::cout << std::setw(spc_wth) << data.get_label_from_idx(i);

            std::size_t total = 0;

            for (std::size_t j = 0; j < num_clusters_detected; ++j) {
                std::size_t val = confusion_matrix[i * num_clusters_detected + j];
                std::cout << std::setw(spc_wth) << val;
                total += val;
            }

            std::cout << std::setw(spc_wth) << total << std::endl;
        }

        {
            std::size_t total = 0;
            std::cout << std::setw(spc_wth) << "Tot";

            for (std::size_t i = 0; i < num_clusters_detected; ++i) {
                std::cout << std::setw(spc_wth) << cluster_total[i];
                total += cluster_total[i];
            }

            std::cout << std::setw(spc_wth) << total << std::endl;
            std::cout << std::endl;
        }

        std::vector<double> best(num_clusters_detected + 1, -1.0);
        quick_mapping(noise_found, num_clusters_detected, num_classes, confusion_matrix, cluster_total, best);

        for (std::size_t i = 0; i < num_clusters_detected; i++) {
            std::cout << "Cluster " << i << " <-- ";

            if (best[i] < 0) {
                std::cout << "no class";
            } else {
                std::cout << data.get_label_from_idx(std::size_t(best[i])) << " (" << best[i] << ")";
            }

            std::cout << std::endl;
        }

        std::cout << std::endl;

        double total_F1_measure = 0.0;

        for (std::size_t i = 0; i < num_classes; i++) {
            bool found = false;
            std::size_t idx = 0;

            for (std::size_t j = 0; j < num_clusters_detected; j++) {
                if (best[j] == i) {
                    idx = j;
                    found = true;
                    break;
                }
            }

            if (!found) {
                std::cout << "Class (" << data.get_label_from_idx(i) << ") F-Measure: 0.0 (not assigned)" << std::endl;
                continue;
            }


            double precision = compute_precision(confusion_matrix, std::size_t(best[idx]), idx, num_classes, num_clusters_detected);
            double recall = compute_recall(confusion_matrix, std::size_t(best[idx]), idx, num_clusters_detected);

            double f1_measure = 2.0 * precision * recall / (precision + recall);

            std::cout << "Class (" << data.get_label_from_idx(i)
                      << ") F-Measure: " << f1_measure
                      << " Precision: " << precision
                      << " Recall: " << recall
                      << std::endl;

            if (!std::isnan(f1_measure)) {
                total_F1_measure += f1_measure;
            }
        }

        std::cout << std::endl;

        std::cout << "Avg F1-Measure: " << (total_F1_measure / num_classes) << std::endl;
    }

}

int main(int argc, char ** argv) {
    parse_command_line(argc, argv);
    load_data();

    std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();

    init();

    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    double time = elapsed.count();
    std::cout << "Init time: " << time << " secs" << std::endl;
    std::cout << std::endl;

    start = std::chrono::steady_clock::now();

    assign();

    end = std::chrono::steady_clock::now();
    elapsed = end - start;
    time = elapsed.count();
    std::cout << "Alg time: " << time << " secs" << std::endl;
    std::cout << std::endl;

    compute_F1_measure();

    clean_up();

    return EXIT_SUCCESS;
}
