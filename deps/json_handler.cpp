#include "json_handler.h"
void FileHandler::read_file(std::string filename) {
    std::ifstream file(filename);

    if (file.is_open()) {
        std::stringstream buffer;
        buffer << file.rdbuf();
        content = buffer.str();
    } else
        throw FileNotOpened();

    file.close();
}

void FileHandler::write_file(std::string filename) {
    std::ofstream file(filename);
    file << content;
    file.close();
}

void JsonHandler::dumps(std::map<std::string, std::string> json) {
    content = "{\n";
    for (const auto& [key, value] : json) {
        content += "    \"" + key + "\": \"" + value + "\",\n";
    }
    content[content.length() - 2] = ' ';
    content += "}";
}

void JsonHandler::dump(std::map<std::string, std::string> json, std::string filename) {
    dumps(json);
    write_file(filename);
}

void JsonHandler::loads(std::string string) {
    std::string key, value;
    short quote_position = 0;

    for (int i = 0; i < string.length(); i++) {
        char c = string[i];
        if (c == '"') {
            quote_position = (quote_position + 1) % 5;
            continue;
        }

        if (quote_position == 1)
            key += c;
        else if (quote_position == 2)
            // values between key and value (i.e. whitespace and '"')
            continue;
        else if (quote_position == 3)
            value += c;
        else if (quote_position == 4) {
            json[key] = value;
            quote_position = (quote_position + 1) % 5;
        }
    }
}
