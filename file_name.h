#ifndef FILE_NAME_H
#define FILE_NAME_H

#include <filesystem>
#include <string>

class file_name_c {
    public:
        file_name_c(std::string const & file_name);
        virtual ~file_name_c() = default;

        std::string get_file_name() const;
        std::string get_file_name_without_ex() const;
        std::string get_full_name() const;

    private:
        std::filesystem::path file_name;
};

#endif // FILE_NAME_H
