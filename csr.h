#ifndef CSR_H
#define CSR_H

#include <map>
#include <string>
#include <vector>

class csr_c {
    typedef std::map<std::string, std::size_t> label_map_t;
    typedef std::pair<std::string, std::size_t> label_pair_t;

    public:
        void add_col_idx(std::size_t idx);
        void add_row_idx(std::size_t idx);

        void add_label(std::string const & label);
        void add_value(double value);

        void clear();

        std::size_t cols() const;
        std::size_t rows() const;

        std::size_t c_idx(std::size_t idx) const;
        std::size_t r_idx(std::size_t idx) const;

        std::string const & label(std::size_t idx) const;
        double value(std::size_t idx) const;

        std::size_t num_labels() const;
        std::size_t get_label_idx(std::size_t idx) const;
        std::string const & get_label_from_idx(std::size_t idx) const;

        std::vector<std::size_t> const & get_col_idx() const;

    private:
        std::vector<std::size_t> row_idx;
        std::vector<std::size_t> col_idx;

        std::vector<std::string> labels;
        std::vector<double> values;

        label_map_t label_map;
        std::vector<std::string> idx_map;

        std::size_t max_col{0};
        std::size_t map_idx{0};
};

#endif // CSR_H
