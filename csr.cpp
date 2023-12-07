#include "csr.h"

void csr_c::add_col_idx(std::size_t idx) {
    col_idx.push_back(idx);
    max_col = std::max(max_col, idx);
}

void csr_c::add_label(std::string const & label) {
    labels.push_back(label);

    label_map_t::iterator const & it = label_map.find(label);

    if (it == label_map.end()) {
        label_map.insert(label_pair_t(label, map_idx));
        map_idx++;
        idx_map.push_back(label);
    }
}

void csr_c::add_row_idx(std::size_t idx) {
    row_idx.push_back(idx);
}

void csr_c::add_value(double value) {
    values.push_back(value);
}

std::size_t csr_c::c_idx(std::size_t idx) const {
#ifdef NDEBUG
    return col_idx[idx];
#else
    return col_idx.at(idx);
#endif
}

void csr_c::clear() {
    row_idx.clear();
    col_idx.clear();

    labels.clear();
    values.clear();
}

std::size_t csr_c::cols() const {
    return max_col + 1;
}

std::size_t csr_c::get_label_idx(std::size_t idx) const {
    std::string lbl = label(idx);
    label_map_t::const_iterator const & it = label_map.find(lbl);

    return it->second;
}

std::string const & csr_c::get_label_from_idx(std::size_t idx) const {
#ifdef NDEBUG
    return idx_map[idx];
#else
    return idx_map.at(idx);
#endif
}

std::string const & csr_c::label(std::size_t idx) const {
#ifdef NDEBUG
    return labels[idx];
#else
    return labels.at(idx);
#endif
}

std::size_t csr_c::num_labels() const {
    return label_map.size();
}

std::size_t csr_c::r_idx(std::size_t idx) const {
#ifdef NDEBUG
    return row_idx[idx];
#else
    return row_idx.at(idx);
#endif
}

std::size_t csr_c::rows() const {
    return labels.size();
}

double csr_c::value(std::size_t idx) const {
#ifdef NDEBUG
    return values[idx];
#else
    return values.at(idx);
#endif
}
