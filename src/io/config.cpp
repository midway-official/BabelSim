#include "babelsim/config.h"

#include <cctype>
#include <cmath>
#include <limits>
#include <set>
#include <fstream>
#include <sstream>
#include <stdexcept>

namespace babelsim {

Parameters::Parameters(const std::filesystem::path& path)
    : m_path(path), m_lines(readConfigLines(path)), m_used(m_lines.size(), false)
{
    std::set<std::string> keys;
    for (const ConfigLine& line : m_lines) {
        if (!keys.insert(line.tokens.front()).second) invalid(line, "duplicate entry");
    }
}

void Parameters::invalid(const ConfigLine& line, const std::string& message) const {
    throw std::runtime_error(m_path.string() + ":" + std::to_string(line.number) +
                             ": " + line.tokens.front() + ": " + message);
}

bool Parameters::contains(const std::string& key) const {
    for (const ConfigLine& line : m_lines) if (line.tokens.front() == key) return true;
    return false;
}

const ConfigLine& Parameters::entry(const std::string& key) const {
    for (std::size_t i = 0; i < m_lines.size(); ++i) {
        if (m_lines[i].tokens.front() == key) {
            m_used[i] = true;
            return m_lines[i];
        }
    }
    throw std::runtime_error(m_path.string() + ": missing entry " + key);
}

double Parameters::number(const std::string& key) const {
    const ConfigLine& line = entry(key);
    if (line.tokens.size() != 2) invalid(line, "expected one number");
    try {
        std::size_t consumed = 0;
        const double value = std::stod(line.tokens[1], &consumed);
        if (consumed == line.tokens[1].size() && std::isfinite(value)) return value;
    } catch (const std::exception&) {
    }
    invalid(line, "expected a finite number");
}

double Parameters::number(const std::string& key, double fallback) const {
    return contains(key) ? number(key) : fallback;
}

double Parameters::positive(const std::string& key) const {
    const double value = number(key);
    if (!(value > 0.0)) invalid(entry(key), "must be positive");
    return value;
}

double Parameters::nonnegative(const std::string& key) const {
    const double value = number(key);
    if (value < 0.0) invalid(entry(key), "must be nonnegative");
    return value;
}

int Parameters::integer(const std::string& key, int fallback) const {
    if (!contains(key)) return fallback;
    const double value = number(key);
    if (std::trunc(value) != value || value < std::numeric_limits<int>::min() ||
        value > std::numeric_limits<int>::max()) invalid(entry(key), "expected an integer");
    return static_cast<int>(value);
}

void Parameters::requireAllUsed() const {
    for (std::size_t i = 0; i < m_lines.size(); ++i) {
        if (!m_used[i]) invalid(m_lines[i], "unused or unknown entry");
    }
}

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
