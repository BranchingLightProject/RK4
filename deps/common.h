#ifndef __COMMON_H
#define __COMMON_H

#include <exception>

struct FileNotOpened : public std::exception {
    const char *what() const throw() {
        return "Unable to open file";
    }
};

#endif __COMMON_H