/*	MCM file compressor

  Copyright (C) 2014, Google Inc.
  Authors: Mathieu Chartier

  LICENSE

    This file is part of the MCM file compressor.

    MCM is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    MCM is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with MCM.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef ARCHIVE_HPP_
#define ARCHIVE_HPP_

#include <thread>
#include <cctype>

#include "CM.hpp"
#include "Compressor.hpp"
#include "File.hpp"
#include "Stream.hpp"
#include "Util.hpp"

class Phase1Profiler {
public:
  static constexpr size_t WINDOW_SIZE = 4096;
  static constexpr size_t STRIDE = 1024;
  static constexpr size_t NUM_CLASSES = 12;

  struct WindowRecord {
    uint64_t offset;
    float entropy;
    uint8_t class_histogram[NUM_CLASSES];
  };

  std::vector<WindowRecord> records;
  uint64_t current_offset = 0;
  double mean_ = 0, stdev_ = 0;

  void log(uint64_t offset, float entropy, const uint8_t* class_counts) {
    WindowRecord r;
    r.offset = offset;
    r.entropy = entropy;
    memcpy(r.class_histogram, class_counts, NUM_CLASSES);
    records.push_back(r);
  }

  static int classify(const std::string& token) {
    // Markup (10) - highest precedence
    if (!token.empty()) {
      if (token[0] == '<' || (token.size() >= 2 && token.substr(0, 2) == "</") || token.back() == '>' ||
          (token[0] == '&' && token.back() == ';')) {
        return 10;
      }
    }
    // Capitalized (2)
    if (!token.empty() && isupper(token[0])) {
      bool is_cap = true;
      for (size_t i = 1; i < token.size(); ++i) {
        if (!islower(token[i])) {
          is_cap = false;
          break;
        }
      }
      if (is_cap) return 2;
    }
    // Alphanumeric (4)
    bool has_letter = false, has_digit = false;
    for (char c : token) {
      if (isalpha(c)) has_letter = true;
      if (isdigit(c)) has_digit = true;
    }
    if (has_letter && has_digit) return 4;
    // Lowercase (0)
    bool is_lower = true;
    for (char c : token) {
      if (!islower(c)) {
        is_lower = false;
        break;
      }
    }
    if (is_lower) return 0;
    // Uppercase (1)
    bool is_upper = true;
    for (char c : token) {
      if (!isupper(c)) {
        is_upper = false;
        break;
      }
    }
    if (is_upper) return 1;
    // Digits (3)
    bool is_digit = true;
    for (char c : token) {
      if (!isdigit(c)) {
        is_digit = false;
        break;
      }
    }
    if (is_digit) return 3;
    // Whitespace (5)
    if (token.find_first_not_of(" \t\n\r") == std::string::npos) return 5;
    // Brackets (7)
    if (token.size() == 1 && (token[0] == '(' || token[0] == ')' || token[0] == '[' || token[0] == ']' || token[0] == '{' || token[0] == '}')) return 7;
    // Operators (8)
    if (token.size() == 1 && (token[0] == '=' || token[0] == '+' || token[0] == '-' || token[0] == '*' || token[0] == '/' || token[0] == '<' || token[0] == '>' || token[0] == '|' || token[0] == '&' || token[0] == '^' || token[0] == '%')) return 8;
    // Quotes (9)
    if (token.size() == 1 && (token[0] == '\'' || token[0] == '"')) return 9;
    // Punctuation (6)
    if (token.size() == 1 && (token[0] == '.' || token[0] == ',' || token[0] == ';' || token[0] == ':' || token[0] == '!' || token[0] == '?')) return 6;
    // Other (11)
    return 11;
  }

  void compute_global_stats() {
    if (records.empty()) return;
    std::vector<double> ents;
    for (auto& r : records) ents.push_back(r.entropy);
    std::sort(ents.begin(), ents.end());
    mean_ = std::accumulate(ents.begin(), ents.end(), 0.0) / ents.size();
    double var = 0;
    for (auto e : ents) var += (e - mean_) * (e - mean_);
    stdev_ = sqrt(var / ents.size());
  }

  std::vector<uint64_t> find_boundaries() {
    std::vector<uint64_t> boundaries;
    if (records.size() < 2) return boundaries;
    for (size_t i = 1; i < records.size(); i++) {
      double H = records[i].entropy;
      double prev_H = records[i-1].entropy;
      if (H > mean_ + 2.5 * stdev_ && H - prev_H > stdev_) {
        boundaries.push_back(records[i].offset);
      }
    }
    return boundaries;
  }

  void write_phase1_file(const std::string& filename, uint64_t total_bytes) {
    compute_global_stats();
    std::ofstream fout(filename, std::ios::binary);
    // magic
    uint32_t magic = 0x50483130; // "PH10"
    fout.write((char*)&magic, 4);
    // window_size
    uint32_t ws = WINDOW_SIZE;
    fout.write((char*)&ws, 4);
    // stride
    uint32_t st = STRIDE;
    fout.write((char*)&st, 4);
    // total_bytes
    fout.write((char*)&total_bytes, 8);
    // num_windows
    uint32_t nw = records.size();
    fout.write((char*)&nw, 4);
    // global stats
    double p75 = 0, p90 = 0, p95 = 0;
    if (!records.empty()) {
      std::vector<double> ents;
      for (auto& r : records) ents.push_back(r.entropy);
      std::sort(ents.begin(), ents.end());
      size_t idx75 = ents.size() * 0.75;
      size_t idx90 = ents.size() * 0.9;
      size_t idx95 = ents.size() * 0.95;
      p75 = ents[idx75];
      p90 = ents[idx90];
      p95 = ents[idx95];
    }
    fout.write((char*)&mean_, 8);
    fout.write((char*)&stdev_, 8);
    fout.write((char*)&p75, 8);
    fout.write((char*)&p90, 8);
    fout.write((char*)&p95, 8);
    // window records
    for (auto& r : records) {
      fout.write((char*)&r.offset, 8);
      fout.write((char*)&r.entropy, 4);
      fout.write((char*)r.class_histogram, NUM_CLASSES);
    }
    // boundaries
    auto boundaries = find_boundaries();
    uint32_t nb = boundaries.size();
    fout.write((char*)&nb, 4);
    for (auto b : boundaries) {
      fout.write((char*)&b, 8);
      float ent = 0;
      for (auto& r : records) {
        if (r.offset == b) {
          ent = r.entropy;
          break;
        }
      }
      fout.write((char*)&ent, 4);
    }
  }

  void dump_csv(std::ostream& os) {
    os << "offset,entropy\n";
    for (auto& r : records) {
      os << r.offset << "," << r.entropy << "\n";
    }
  }
};

// Force filter
enum FilterType {
  kFilterTypeNone,
  kFilterTypeDict,
  kFilterTypeX86,
  kFilterTypeAuto,
  kFilterTypeCount,
};

enum CompLevel {
  kCompLevelStore,
  kCompLevelTurbo,
  kCompLevelFast,
  kCompLevelMid,
  kCompLevelHigh,
  kCompLevelMax,
  kCompLevelSimple,
};
std::ostream& operator<<(std::ostream& os, CompLevel comp_level);

enum LZPType {
  kLZPTypeAuto,
  kLZPTypeEnable,
  kLZPTypeDisable,
};

class CompressionOptions {
public:
  static const size_t kDefaultMemUsage = 6;
  static const CompLevel kDefaultLevel = kCompLevelMid;
  static const FilterType kDefaultFilter = kFilterTypeAuto;
  static const LZPType kDefaultLZPType = kLZPTypeAuto;

public:
  size_t mem_usage_ = kDefaultMemUsage;
  CompLevel comp_level_ = kDefaultLevel;
  FilterType filter_type_ = kDefaultFilter;
  LZPType lzp_type_ = kDefaultLZPType;
  std::string dict_file_;
  std::string out_dict_file_;
};

// File headers are stored in a list of blocks spread out through data.
class Archive {
public:
  class Header {
  public:
    static const size_t kCurMajorVersion = 0;
    static const size_t kCurMinorVersion = 84;
    static const size_t kMagicStringLength = 10;

    static const char* getMagic() {
      return "MCMARCHIVE";
    }
    Header();
    void read(Stream* stream);
    void write(Stream* stream);
    bool isArchive() const;
    bool isSameVersion() const;
    uint16_t majorVersion() const {
      return major_version_;
    }
    uint16_t minorVersion() const {
      return minor_version_;
    }

  private:
    char magic_[10]; // MCMARCHIVE
    uint16_t major_version_ = kCurMajorVersion;
    uint16_t minor_version_ = kCurMinorVersion;
  };

  class Algorithm {
  public:
    Algorithm() {}
    Algorithm(const CompressionOptions& options, Detector::Profile profile);
    Algorithm(Stream* stream);
    // Freq is the approximate distribution of input frequencies for the compressor.
    Compressor* CreateCompressor(const FrequencyCounter<256>& freq);
    void read(Stream* stream);
    void write(Stream* stream);
    Filter* createFilter(Stream* stream, Analyzer* analyzer, Archive& archive, size_t opt_var = 0);
    Detector::Profile profile() const {
      return profile_;
    }

  private:
    uint8_t mem_usage_;
    Compressor::Type algorithm_;
    bool lzp_enabled_;
    FilterType filter_;
    Detector::Profile profile_;
  };

  class SolidBlock {
  public:
    Algorithm algorithm_;
    std::vector<FileSegmentStream::FileSegments> segments_;
    // Not stored, obtianed from segments.
    uint64_t total_size_ = 0u;

    SolidBlock() = default;
    SolidBlock(const Algorithm& algorithm) : algorithm_(algorithm) {}
    void write(Stream* stream);
    void read(Stream* stream);
  };

  class Blocks : public std::vector<std::unique_ptr<SolidBlock>> {
  public:
    void write(Stream* stream);
    void read(Stream* stream);
  };

  // Compression.
  Archive(Stream* stream, const CompressionOptions& options);

  // Decompression.
  Archive(Stream* stream);

  // Compression with profiling.
  Archive(Stream* stream, const CompressionOptions& options, Phase1Profiler* profiler);

  // Construct blocks from analyzer.
  void constructBlocks(Analyzer::Blocks* blocks_for_file);

  const Header& getHeader() const {
    return header_;
  }

  CompressionOptions& Options() {
    return options_;
  }

  bool setOpt(size_t var) {
    opt_var_ = var;
    return true;
  }

  bool setOpts(size_t* vars) {
    opt_vars_ = vars;
    setOpt(vars[0]);
    return true;
  }

  void writeBlocks();
  void readBlocks();

  // Analyze and compress. Returns how many bytes wre compressed.
  uint64_t compress(const std::vector<FileInfo>& in_files);

  // Decompress.
  void decompress(const std::string& out_dir, bool verify = false);

  // List files and info.
  void list();

  size_t* opt_vars_ = nullptr;
private:
  Stream* stream_;
  Header header_;
  CompressionOptions options_;
  size_t opt_var_;  
  FileList files_;  // File list.
  Blocks blocks_;
  Phase1Profiler* profiler_ = nullptr;

  void init();
  Compressor* createMetaDataCompressor();
};

#endif
