#include "babelsim/config.h"

#include <cctype>
#include <fstream>
#include <sstream>
#include <stdexcept>

namespace babelsim {

std::vector<ConfigLine> readConfigLines(const std::filesystem::path& path) {
    std::ifstream input(path);
    if (!input) {
        throw std::runtime_error("cannot open BabelSim dictionary: " + path.string());
    }

    std::vector<ConfigLine> lines;
    std::string text;
    std::size_t number = 0;
    while (std::getline(input, text)) {
        ++number;
        const std::size_t comment = text.find('#');
        if (comment != std::string::npos) {
            text.resize(comment);
        }
        std::istringstream words(text);
        ConfigLine line;
        line.number = number;
        for (std::string word; words >> word;) {
            line.tokens.push_back(std::move(word));
        }
        if (!line.tokens.empty()) {
            lines.push_back(std::move(line));
        }
    }
    return lines;
}

std::vector<ConfigToken> readConfigTokens(const std::filesystem::path& path) {
    std::ifstream input(path);
    if (!input) {
        throw std::runtime_error("cannot open BabelSim dictionary: " + path.string());
    }

    std::vector<ConfigToken> tokens;
    std::string text;
    std::size_t number = 0;
    while (std::getline(input, text)) {
        ++number;
        const std::size_t comment = text.find('#');
        if (comment != std::string::npos) text.resize(comment);
        std::string word;
        const auto append = [&] {
            if (!word.empty()) {
                tokens.push_back({number, std::move(word)});
                word.clear();
            }
        };
        for (char character : text) {
            if (std::isspace(static_cast<unsigned char>(character))) {
                append();
            } else if (character == '{' || character == '}' ||
                       character == '(' || character == ')') {
                append();
                tokens.push_back({number, std::string(1, character)});
            } else {
                word += character;
            }
        }
        append();
    }
    return tokens;
}

}  // babelsim 命名空间
