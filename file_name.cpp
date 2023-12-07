#include "file_name.h"

file_name_c::file_name_c(std::string const & file_name_) : file_name(file_name_) {
}

std::string file_name_c::get_file_name() const {
    return file_name.filename().string();
}

std::string file_name_c::get_full_name() const {
    return file_name.string();
}
