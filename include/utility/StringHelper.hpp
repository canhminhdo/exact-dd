//
// Created by Canh Minh Do on 2026/07/14.
//

#ifndef STRINGHELPER_HPP
#define STRINGHELPER_HPP

#include <string>

class StringHelper {
public:
    static std::string toLowerCase(const std::string &str);

    static std::string toUpperCase(const std::string &str);

    static bool startsWith(const std::string &str, const std::string &prefix);

    static bool endsWith(const std::string &str, const std::string &suffix);
};

#endif//STRINGHELPER_HPP
