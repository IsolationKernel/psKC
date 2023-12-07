#include "libsvm.h"

#include "csr.h"
#include "file_name.h"

#include <fstream>
#include <iostream>

libsvm_c::libsvm_c(file_name_c const & fn) : file_name(fn) {
}

void libsvm_c::load(csr_c & data) {
    std::ifstream input;
    input.open(file_name.get_full_name());

    if (!input.is_open()) {
        std::cerr << "Error: unable to open " << file_name.get_file_name() << std::endl;
        exit(EXIT_SUCCESS);
    }

    data.clear();

    std::size_t r_idx = 0;
    data.add_row_idx(r_idx);

    for (std::string line; std::getline(input, line); ) {
        std::stringstream is(line);
        is >> std::skipws;

        std::string label;

        if (!(is >> label)) {
            break;
        }

        data.add_label(label);

        for (;; r_idx++) {
            std::size_t idx;
            double value;
            char token;

            if (!(is >> idx)) {
                break;
            }

            is >> token; // :
            is >> value;

            data.add_col_idx(idx - 1); // libsvm is one based not zero
            data.add_value(value);
        }

        data.add_row_idx(r_idx);
    }

    input.close();
}
