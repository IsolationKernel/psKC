#ifndef PARSE_PARAM_H
#define PARSE_PARAM_H

#include <algorithm>
#include <iostream>
#include <sstream>

// https://stackoverflow.com/questions/865668/how-to-parse-command-line-arguments-in-c
template <typename T>
    void parse_param(int argc, char ** argv, std::string const & option, T & result, std::string const & err_msg) {
        char ** end = argv + argc;
        char ** itr = std::find(argv, end, option);

        if ((itr != end) && (++itr != end)) {
            std::stringstream is(*itr);
            is >> std::skipws;
            is >> result;
        } else {
            std::cerr << "Error: " << err_msg << std::endl;
            exit(EXIT_FAILURE);
        }
    }

void parse_param(int argc, char ** argv, std::string const & option, bool & result) {
    char ** end = argv + argc;
    char ** itr = std::find(argv, end, option);

    result = (itr != end);
}

#endif // PARSE_PARAM_H
