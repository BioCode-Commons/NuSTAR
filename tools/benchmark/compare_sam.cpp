#include <algorithm>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

struct Reader {
    std::ifstream input;
    std::vector<std::string> sq;
    std::string pending, name;
    unsigned long long records = 0;
    bool ready = false;
    explicit Reader(const char* path) : input(path) {
        if (!input) throw std::runtime_error("cannot open SAM");
        advance();
    }
    void advance() {
        std::string line;
        ready = false;
        while (std::getline(input, line)) {
            if (!line.empty() && line.back() == '\r') line.pop_back();
            if (line.starts_with("@")) {
                if (line.starts_with("@SQ\t")) sq.push_back(line);
                continue;
            }
            std::vector<std::string_view> fields;
            std::string_view tail(line);
            for (;;) {
                const auto tab = tail.find('\t');
                fields.push_back(tail.substr(0, tab));
                if (tab == tail.npos) break;
                tail.remove_prefix(tab + 1);
            }
            if (fields.size() < 11 || fields[0].empty()) throw std::runtime_error("invalid SAM record");
            name = fields[0];
            std::sort(fields.begin()+11, fields.end());
            pending.clear();
            for (size_t i=0; i<fields.size(); ++i) {
                if (i) pending += '\t';
                pending += fields[i];
            }
            ++records;
            ready = true;
            return;
        }
        if (!input.eof()) throw std::runtime_error("SAM read error");
    }
    std::vector<std::string> group() {
        const auto qname = name;
        std::vector<std::string> result;
        do {
            result.push_back(pending);
            advance();
        } while (ready && name == qname);
        std::sort(result.begin(), result.end());
        return result;
    }
};
int main(int argc, char** argv) {
    if (argc == 3 && std::string(argv[1]) == "--canonical") {
        try {
            std::ios::sync_with_stdio(false);
            Reader reader(argv[2]);
            auto emit = [](const std::string& value) {
                std::cout << value.size() << ':' << value;
            };
            std::cout << "SAM-canonical-v1\n";
            for (const auto& sq : reader.sq) emit(sq);
            std::cout << "\nBODY\n";
            unsigned long long groups = 0;
            while (reader.ready) {
                const auto rows = reader.group();
                std::cout << rows.size() << '[';
                for (const auto& row : rows) emit(row);
                std::cout << ']';
                ++groups;
            }
            std::cout.flush();
            if (!std::cout) throw std::runtime_error("canonical stream write error");
            std::cerr << "{\"records\":" << reader.records << ",\"groups\":" << groups << "}\n";
            return 0;
        } catch (const std::exception& e) { std::cerr << e.what() << '\n'; return 2; }
    }
    if (argc != 4) return 2;
    try {
        Reader left(argv[1]), right(argv[2]);
        std::ofstream diff(argv[3]);
        if (!diff) throw std::runtime_error("cannot create differences file");
        unsigned long long groups=0, different=0;
        while (left.ready || right.ready) {
            auto a = left.ready ? left.group() : std::vector<std::string>{};
            auto b = right.ready ? right.group() : std::vector<std::string>{};
            ++groups;
            if (a != b) {
                if (different++ < 10) {
                    diff << "Group " << groups << "\nLEFT\n";
                    for (const auto& row:a) diff << row << '\n';
                    diff << "RIGHT\n";
                    for (const auto& row:b) diff << row << '\n';
                }
            }
        }
        const bool sq_equal=left.sq==right.sq;
        std::cout << "{\"equal\":" << (different==0 && sq_equal ? "true":"false")
                  << ",\"left_records\":" << left.records << ",\"right_records\":" << right.records
                  << ",\"groups\":" << groups << ",\"different_groups\":" << different
                  << ",\"sq_equal\":" << (sq_equal?"true":"false") << "}\n";
        return different==0 && sq_equal ? 0:1;
    } catch (const std::exception& e) { std::cerr << e.what() << '\n'; return 2; }
}
