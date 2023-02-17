/**
 * Json Handler main module. @jpvanegasc/FileHandler
 *
 * @author: Juan Vanegas. git: jpvanegasc
 *
 * Last modification: Feb 16, 2023.
 */
#ifndef __JSON_HANDLER_H
#define __JSON_HANDLER_H

#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <sstream>
#include <string>

#include "common.h"

/**
 * Mixin to read and write a file
 *
 * @param filename: name or path of file.
 */
class FileHandler {
   public:
    std::string content;
    void read_file(std::string filename);
    void write_file(std::string filename);
};

/**
 * Read and write a json file
 *
 * @param filename: name or path of file.
 * @param delim_char : file data delimiter. Default is ','.
 * @param comm_char : file comment indicator. Default is '#'
 */
class JsonHandler : public FileHandler {
   public:
    std::map<std::string, std::string> json;
    void dumps(std::map<std::string, std::string> json);
    void dump(std::map<std::string, std::string> json, std::string filename);

    void loads(std::string string);
    void load(std::string filename);
};

#endif  // __JSON_HANDLER_H