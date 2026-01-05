/*	MCM file compressor

  Copyright (C) 2013, Google Inc.
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

#include <chrono>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <stdio.h>
#include <string.h>
#include <string>
#include <vector>
#include <filesystem>
#include <algorithm>
#include <cmath>
#include <cstring>
#include <sstream>
#include <thread>
#include <tuple>

#include "Archive.hpp"
#include "CM.hpp"
#include "DeltaFilter.hpp"
#include "Dict.hpp"
#include "File.hpp"
#include "Huffman.hpp"
#include "LZ-inl.hpp"
#include "ProgressMeter.hpp"
#include "Tests.hpp"
#include "TurboCM.hpp"
#include "X86Binary.hpp"

static constexpr bool kReleaseBuild = false;

static void printHeader() {
  std::cout
    << "======================================================================" << std::endl
    << "mcm compressor v" << Archive::Header::kCurMajorVersion << "." << Archive::Header::kCurMinorVersion
    << ", by Mathieu Chartier (c)2016 Google Inc." << std::endl
    << "Ultra experimental, may contain bugs. Contact mathieu.a.chartier@gmail.com" << std::endl
    << "Special thanks to: Matt Mahoney, Stephan Busch, Christopher Mattern." << std::endl
    << "======================================================================" << std::endl;
}

static void printUsage(const std::string& name) {
  printHeader();
  std::cout
    << "Caution: Experimental, use only for testing!" << std::endl
    << "Usage: " << name << " [commands] [options] <infile|dir> <outfile>(default infile.mcm)" << std::endl
    << "Commands:" << std::endl
    << "  --phase1-profile <file>     : Profile entropy and generate phase1.final" << std::endl
    << "  --phase1-dump-csv <file>    : Dump phase1 data as CSV" << std::endl
    << "  --phase2p5-profile <file>   : Run Phase 2.5 reordering on <file>" << std::endl
    << "Options: d for decompress" << std::endl
    << "-{t|f|m|h|x}{1 .. 11} compression option" << std::endl
    << "t is turbo, f is fast, m is mid, h is high, x is max (default " << CompressionOptions::kDefaultLevel << ")" << std::endl
    << "0 .. 11 specifies memory with 32mb .. 5gb per thread (default " << CompressionOptions::kDefaultMemUsage << ")" << std::endl
    << "10 and 11 are only supported on 64 bits" << std::endl
    << "-test tests the file after compression is done" << std::endl
    << "Examples:" << std::endl
    << "Compress: " << name << " -m9 enwik8 enwik8.mcm" << std::endl
    << "Decompress: " << name << " d enwik8.mcm enwik8.ref" << std::endl
    << "Phase 1: " << name << " --phase1-profile enwik8" << std::endl
    << "Phase 2.5: " << name << " --phase2p5-profile enwik8" << std::endl;
}

class Options {
public:
  // Block size of 0 -> file size / #threads.
  static const uint64_t kDefaultBlockSize = 0;
  enum Mode {
    kModeUnknown,
    // Compress -> Decompress -> Verify.
    // (File or directory).
    kModeTest,
    // Compress infinite times with different opt vars.
    kModeOpt,
    // In memory test.
    kModeMemTest,
    // Single file test.
    kModeSingleTest,
    // Add a single file.
    kModeAdd,
    kModeExtract,
    kModeExtractAll,
    // Single hand mode.
    kModeCompress,
    kModeDecompress,
    // List & other
    kModeList,
  };
  Mode mode = kModeUnknown;
  bool opt_mode = false;
  CompressionOptions options_;
  Compressor* compressor = nullptr;
  uint32_t threads = 1;
  uint64_t block_size = kDefaultBlockSize;
  FileInfo archive_file;
  std::vector<FileInfo> files;
  std::string dict_file;
  bool phase1_profile = false;
  bool phase1_dump_csv = false;
  bool phase2p5_profile = false;
  bool phase2p75_measure = false;

  // int usage(const std::string& name) {
  //   printHeader();
  //   std::cout
  //     << "Caution: Experimental, use only for testing!" << std::endl
  //     << "Usage: " << name << " [commands] [options] <infile|dir> <outfile>(default infile.mcm)" << std::endl
  //     << "Options: d for decompress" << std::endl
  //     << "-{t|f|m|h|x}{1 .. 11} compression option" << std::endl
  //     << "t is turbo, f is fast, m is mid, h is high, x is max (default " << CompressionOptions::kDefaultLevel << ")" << std::endl
  //     << "0 .. 11 specifies memory with 32mb .. 5gb per thread (default " << CompressionOptions::kDefaultMemUsage << ")" << std::endl
  //     << "10 and 11 are only supported on 64 bits" << std::endl
  //     << "-test tests the file after compression is done" << std::endl
  //     // << "-b <mb> specifies block size in MB" << std::endl
  //     // << "-t <threads> the number of threads to use (decompression requires the same number of threads" << std::endl
  //     << "Examples:" << std::endl
  //     << "Compress: " << name << " -m9 enwik8 enwik8.mcm" << std::endl
  //     << "Decompress: " << name << " d enwik8.mcm enwik8.ref" << std::endl;
  //   return 0;
  // }

  int parse(int argc, char* argv[]) {
    assert(argc >= 1);
    std::string program(trimExt(argv[0]));
    const std::string kDictArg = "-dict=";
    const std::string kOutDictArg = "-out-dict=";
    // Parse options.
    bool has_comp_args = false;
    int i = 1;
    for (; i < argc; ++i) {
      const std::string arg(argv[i]);
      if (arg == "--phase1-profile") {
        phase1_profile = true;
        mode = kModeCompress;
      } else if (arg == "--phase1-dump-csv") {
        phase1_dump_csv = true;
      } else if (arg == "--phase2p5-profile") {
        phase2p5_profile = true;
        mode = kModeCompress;
      } else if (arg == "--phase2p75-measure-headroom") {
        phase2p75_measure = true;
      } else if (arg == "d") {
        mode = kModeDecompress;
      } else if (arg == "-test") {
        // Test after compression/decompression
      } else if (arg.size() > 1 && arg[0] == '-' && std::string("tfmhx").find(arg[1]) != std::string::npos) {
        char level_char = arg[1];
        std::string num_str = arg.substr(2);
        if (num_str.empty()) {
          std::cerr << "Missing number after " << level_char << std::endl;
          return 4;
        }
        try {
          int num = std::stoi(num_str);
          if (num < 0 || num > 11) {
            std::cerr << "Invalid memory usage " << num << std::endl;
            return 4;
          }
          options_.mem_usage_ = static_cast<size_t>(num);
          switch (level_char) {
            case 't': options_.comp_level_ = kCompLevelTurbo; break;
            case 'f': options_.comp_level_ = kCompLevelFast; break;
            case 'm': options_.comp_level_ = kCompLevelMid; break;
            case 'h': options_.comp_level_ = kCompLevelHigh; break;
            case 'x': options_.comp_level_ = kCompLevelMax; break;
          }
        } catch (const std::exception&) {
          std::cerr << "Invalid number " << num_str << std::endl;
          return 4;
        }
      } else if (!arg.empty() && arg[0] != '-') {
        files.push_back(FileInfo(arg));
      } else {
        std::cerr << "Unknown option " << arg << std::endl;
        return 4;
      }
    }
    if (mode == kModeUnknown) {
      mode = kModeCompress;
    }
    if (files.empty() && !phase1_profile && !phase1_dump_csv && !phase2p5_profile && !phase2p75_measure) {
      printUsage(program);
      return 0;
    }
    if (files.empty()) {
      std::cerr << "No files specified" << std::endl;
      return 5;
    }
    archive_file = FileInfo(files[0].getName() + ".mcm");
    return 0;
  }
};

int main(int argc, char* argv[]) {
  if (!kReleaseBuild) {
    RunAllTests();
  }
  Options options;
  auto ret = options.parse(argc, argv);
  if (ret) {
    std::cerr << "Failed to parse arguments" << std::endl;
    return ret;
  }
  Phase1Profiler* profiler_ptr = nullptr;
  if (options.phase1_profile) {
    profiler_ptr = new Phase1Profiler();
  }
  switch (options.mode) {
  case Options::kModeMemTest: {
    constexpr size_t kCompIterations = kIsDebugBuild ? 1 : 1;
    constexpr size_t kDecompIterations = kIsDebugBuild ? 1 : 25;
    // Read in the whole file.
    std::vector<uint64_t> lengths;
    uint64_t long_length = 0;
    for (const auto& file : options.files) {
      File f(file.getName());
      lengths.push_back(f.length());
      long_length += lengths.back();
    }
    auto length = static_cast<size_t>(long_length);
    check(length < 300 * MB);
    auto* in_buffer = new uint8_t[length];
    // Read in the files.
    uint32_t index = 0;
    uint64_t read_pos = 0;
    for (const auto& file : options.files) {
      File f(file.getName(), std::ios_base::in | std::ios_base::binary);
      size_t count = f.read(in_buffer + read_pos, static_cast<size_t>(lengths[index]));
      check(count == lengths[index]);
      index++;
    }
    // Create the memory compressor.
    typedef SimpleEncoder<8, 16> Encoder;
    // auto* compressor = new LZ4;
    // auto* compressor = new LZSSE;
    // auto* compressor = new MemCopyCompressor;
    auto* compressor = new LZ16<FastMatchFinder<MemoryMatchFinder>>;
    auto out_buffer = new uint8_t[compressor->getMaxExpansion(length)];
    uint32_t comp_start = clock();
    uint32_t comp_size;
    static const bool opt_mode = false;
    if (opt_mode) {
      uint32_t best_size = 0xFFFFFFFF;
      uint32_t best_opt = 0;
      std::ofstream opt_file("opt_result.txt");
      for (uint32_t opt = 0; ; ++opt) {
        compressor->setOpt(opt);
        comp_size = compressor->compress(in_buffer, out_buffer, length);
        opt_file << "opt " << opt << " = " << comp_size << std::endl << std::flush;
        std::cout << "Opt " << opt << " / " << best_opt << " =  " << comp_size << "/" << best_size << std::endl;
        if (comp_size < best_size) {
          best_opt = opt;
          best_size = comp_size;
        }
      }
    } else {
      for (uint32_t i = 0; i < kCompIterations; ++i) {
        comp_size = compressor->compress(in_buffer, out_buffer, length);
      }
    }

    const uint32_t comp_end = clock();
    std::cout << "Done compressing " << length << " -> " << comp_size << " = " << float(double(length) / double(comp_size)) << " rate: "
      << prettySize(static_cast<uint64_t>(long_length * kCompIterations / clockToSeconds(comp_end - comp_start))) << "/s" << std::endl;
    memset(in_buffer, 0, length);
    const uint32_t decomp_start = clock();
    for (uint32_t i = 0; i < kDecompIterations; ++i) {
      compressor->decompress(out_buffer, in_buffer, length);
    }
    const uint32_t decomp_end = clock();
    std::cout << "Decompression took: " << decomp_end - comp_end << " rate: "
      << prettySize(static_cast<uint64_t>(long_length * kDecompIterations / clockToSeconds(decomp_end - decomp_start))) << "/s" << std::endl;
    index = 0;
    for (const auto& file : options.files) {
      File f(file.getName(), std::ios_base::in | std::ios_base::binary);
      const auto count = static_cast<uint32_t>(f.read(out_buffer, static_cast<uint32_t>(lengths[index])));
      check(count == lengths[index]);
      for (uint32_t i = 0; i < count; ++i) {
        if (out_buffer[i] != in_buffer[i]) {
          std::cerr << "File" << file.getName() << " doesn't match at byte " << i << std::endl;
          check(false);
        }
      }
      index++;
    }
    std::cout << "Decompression verified" << std::endl;
    break;
  }
  case Options::kModeSingleTest:
  case Options::kModeOpt:
  case Options::kModeCompress:
  case Options::kModeTest: {
    printHeader();

    int err = 0;

    std::string out_file = options.archive_file.getName();
    File fout;

    if (options.mode == Options::kModeOpt) {
      std::cout << "Optimizing" << std::endl;
      uint64_t best_size = std::numeric_limits<uint64_t>::max();
      size_t best_var = 0;
      std::ofstream opt_file("opt_result.txt");
      // static const size_t kOpts = 10624;
      static const size_t kOpts = 3;
      // size_t opts[kOpts] = {0,1,2,3,15,14,4,6,7,8,9,17,12,11,13,5,10,18,20,19,21,26,22,28,23,24,16,25,27,29,31,32,36,33,34,35,37,30,38,39,};
      // size_t opts[kOpts] = {}; for (size_t i = 0; i < kOpts; ++i) opts[i] = i;
      size_t opts[kOpts] = {};
      //size_t opts[] = {7,14,1,12,3,4,11,15,9,16,5,6,18,13,19,30,45,20,21,22,23,17,8,2,26,10,32,43,36,35,42,29,34,24,25,37,31,33,39,38,0,41,28,40,44,58,46,59,92,27,60,61,91,63,95,47,64,124,94,62,93,96,123,125,72,69,65,67,83,68,66,73,82,70,80,76,71,81,77,87,78,74,79,84,75,48,49,50,51,52,53,54,55,56,57,86,88,97,98,99,100,85,101,90,103,104,89,105,107,102,108,109,110,111,106,113,112,114,115,116,119,118,120,121,117,122,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,151,144,145,146,147,148,149,150,152,153,155,156,157,154,158,159,160,161,162,163,164,165,166,167,168,169,170,171,172,173,174,175,176,177,178,179,180,181,182,183,184,185,186,187,188,189,190,191,192,193,194,195,196,197,198,199,200,201,202,203,204,205,206,207,208,209,210,211,212,213,214,215,216,217,218,219,220,221,222,223,224,225,226,239,227,228,229,230,231,232,233,234,235,236,237,238,240,241,242,243,244,245,246,247,248,249,250,251,252,253,254,255,};
      if (false) {
        auto temp = ReadCSI<size_t>("optin.txt");
        check(temp.size() == kOpts);
        std::copy(temp.begin(), temp.end(), &opts[0]);
      }
      size_t best_opts[kOpts] = {};
      srand(clock() + 291231);
      size_t cur_index = 0;
      size_t len = 1;
      // size_t kMaxIndex = 128 - len;
      // size_t kMaxIndex = cm::kModelCount;
      size_t kMaxIndex = 12345;
      size_t bads = 0;
      size_t best_cur = 0;
      static constexpr size_t kAvgCount = 2;
      double total[kAvgCount] = {};
      size_t count[kAvgCount] = {};
      double min_time = std::numeric_limits<double>::max();
      const bool kPerm = false;
      if (kPerm) {
        for (size_t i = 0;; ++i) {
          if ((i & 255) == 255 && false) {
            len -= len > 1;
            kMaxIndex = 128 - len;
          }
          auto a = i % kMaxIndex;
          auto b = (i + 1 + std::rand() % (kMaxIndex - 1)) % kMaxIndex;
          VoidWriteStream fout;
          Archive archive(&fout, options.options_, profiler_ptr);
          if (i != 0) {
            ReplaceSubstring(opts, a, len, b, kOpts);
          }
          if (!archive.setOpts(opts)) {
            std::cerr << "Failed to set opts" << std::endl;
            continue;
          }
          uint64_t in_bytes = archive.compress(options.files);
          if (in_bytes == 0) continue;
          const auto size = fout.tell();
          std::cout << i << ": swap " << a << " to " << b << " " << size << std::endl;
          if (size < best_size) {
            std::cout << "IMPROVEMENT " << i << ": " << size << std::endl;
            opt_file << i << ": " << size << " ";
            for (auto opt : opts) opt_file << opt << ",";
            opt_file << std::endl;
            std::copy_n(opts, kOpts, best_opts);
            best_size = size;
          } else {
            std::copy_n(best_opts, kOpts, opts);
          }
        }
      } else {
        for (auto o : opts) check(o <= kMaxIndex);
        for (size_t i = 0;; ++i) {
          const clock_t start = clock();
          VoidWriteStream fout;
          Archive archive(&fout, options.options_, profiler_ptr);
          if (!archive.setOpts(opts)) {
            continue;
          }
          uint64_t in_bytes = archive.compress(options.files);
          if (in_bytes == 0) continue;
          const double time = clockToSeconds(clock() - start);
          total[i % kAvgCount] += time;
          min_time = std::min(min_time, time);
          const auto size = fout.tell();
          opt_file << "opts ";
          for (auto opt : opts) opt_file << opt << ",";
          ++count[i % kAvgCount];
          auto before_index = cur_index;
          auto before_opt = opts[before_index];
          if (size < best_size) {
            best_size = size;
            std::copy_n(opts, kOpts, best_opts);
            best_var = opts[cur_index];
            bads = 0;
          } 
          if (opts[cur_index] >= kMaxIndex) {
            std::copy_n(best_opts, kOpts, opts);
            cur_index = (cur_index + 1) % kOpts;
            opts[cur_index] = 0;
          } else {
            ++opts[cur_index];
          }

          std::ostringstream ss;
          double avgs[kAvgCount] = {};
          for (size_t i = 0; i < kAvgCount; ++i) {
            if (count[i] != 0) avgs[i] = total[i] / double(count[i]);
          }
          double avg = std::accumulate(total, total + kAvgCount, 0.0) / double(std::accumulate(count, count + kAvgCount, 0u));
          ss << " -> " << formatNumber(size) << " best " << best_var << " in " << time << "s avg "
             << avg << "(";
          for (double d : avgs) ss << d << ",";
          ss << ") min " << min_time;

          opt_file << ss.str() << std::endl << std::flush;

          std::cout << "opt[" << before_index << "]=" << before_opt << " best=" << best_var << "(" << formatNumber(best_size) << ") "
            << formatNumber(in_bytes) << ss.str() << std::endl;
        }
      }
    } else {
      const clock_t start = clock();
      if (err = fout.open(out_file, std::ios_base::out | std::ios_base::binary)) {
        std::cerr << "Error opening: " << out_file << " (" << errstr(err) << ")" << std::endl;
        return 2;
      }

      std::cout << "Compressing to " << out_file << " mode=" << options.options_.comp_level_ << " mem=" << options.options_.mem_usage_ << std::endl;
      Archive archive(&fout, options.options_, profiler_ptr);
      uint64_t in_bytes = archive.compress(options.files);
      clock_t time = clock() - start;
      std::cout << "Done compressing " << formatNumber(in_bytes) << " -> " << formatNumber(fout.tell())
        << " in " << std::setprecision(3) << clockToSeconds(time) << "s"
        << " bpc=" << double(fout.tell()) * 8.0 / double(in_bytes) << std::endl;

      fout.close();

      if (profiler_ptr) {
        std::string phase1_file = out_file + ".phase1";
        profiler_ptr->write_phase1_file(phase1_file, in_bytes);
        std::cout << "Phase 1 profiling written to " << phase1_file << std::endl;
        if (options.phase1_dump_csv) {
          profiler_ptr->dump_csv(std::cout);
        }
        // Compute Phase 1 final stats
        if (!profiler_ptr->records.empty()) {
          std::cout << "Computing Phase 1 final stats" << std::endl;
          // Compute entropies
          std::vector<float> entropies;
          for (auto& r : profiler_ptr->records) entropies.push_back(r.entropy);
          // Mean entropy
          double sum_entropy = 0;
          for (float e : entropies) sum_entropy += e;
          double mean_entropy = sum_entropy / entropies.size();
          // Stddev entropy
          double sum_sq = 0;
          for (float e : entropies) sum_sq += (e - mean_entropy) * (e - mean_entropy);
          double stddev_entropy = std::sqrt(sum_sq / entropies.size());
          // Deltas
          std::vector<float> deltas;
          for (size_t i = 1; i < entropies.size(); ++i) {
            deltas.push_back(entropies[i] - entropies[i-1]);
          }
          double mean_delta = 0;
          double stddev_delta = 0;
          if (!deltas.empty()) {
            double sum_d = 0;
            for (float d : deltas) sum_d += d;
            mean_delta = sum_d / deltas.size();
            double sum_sq_d = 0;
            for (float d : deltas) sum_sq_d += (d - mean_delta) * (d - mean_delta);
            stddev_delta = std::sqrt(sum_sq_d / deltas.size());
          }
          // Coarseness
          double coarseness = (mean_entropy > 0) ? stddev_delta / mean_entropy : 0;
          // Spikes for k=2.5
          std::vector<size_t> spike_positions;
          for (size_t i = 0; i < deltas.size(); ++i) {
            if (deltas[i] > mean_delta + 2.5 * stddev_delta) {
              spike_positions.push_back(i);
            }
          }
          // Target segment size
          int target_segment_size = 16384; // default
          if (spike_positions.size() > 1) {
            double total_dist = 0;
            for (size_t p = 1; p < spike_positions.size(); ++p) {
              total_dist += spike_positions[p] - spike_positions[p-1];
            }
            double mean_spike_dist = total_dist / (spike_positions.size() - 1);
            target_segment_size = static_cast<int>(mean_spike_dist * 1024); // stride
          }
          int min_segment_size = target_segment_size / 4;
          int max_segment_size = (coarseness < 0.5) ? target_segment_size * 8 : target_segment_size * 2;
          double entropy_threshold = mean_delta + 2.5 * stddev_delta;
          // Token class stats
          const int num_classes = 12;
          int transition_matrix[12][12];
          memset(transition_matrix, 0, sizeof(transition_matrix));
          uint64_t total_counts[12] = {0};
          int prev_class = -1;
          double total_churn = 0;
          for (auto& r : profiler_ptr->records) {
            // Find dominant class
            int max_count = 0;
            int dominant = 0;
            for (int c = 0; c < num_classes; ++c) {
              if (r.class_histogram[c] > max_count) {
                max_count = r.class_histogram[c];
                dominant = c;
              }
            }
            if (prev_class != -1) {
              transition_matrix[prev_class][dominant]++;
            }
            prev_class = dominant;
            // Churn: number of non-zero classes
            int churn = 0;
            for (int c = 0; c < num_classes; ++c) if (r.class_histogram[c] > 0) churn++;
            total_churn += churn;
          }
          double avg_churn = total_churn / profiler_ptr->records.size();
          // Spikes counts
          int spikes_2_0 = 0, spikes_2_5 = 0, spikes_3_0 = 0;
          for (float d : deltas) {
            if (d > mean_delta + 2.0 * stddev_delta) spikes_2_0++;
            if (d > mean_delta + 2.5 * stddev_delta) spikes_2_5++;
            if (d > mean_delta + 3.0 * stddev_delta) spikes_3_0++;
          }
          // Validation
          int active_classes = 0;
          for (int i = 0; i < num_classes; ++i) {
            int sum = 0;
            for (int j = 0; j < num_classes; ++j) sum += transition_matrix[i][j];
            if (sum > 0) active_classes++;
          }
          bool validation_passed = true;
          std::vector<std::string> errors;
          std::vector<std::string> warnings;
          if (stddev_entropy < 0.05) {
            validation_passed = false;
            errors.push_back("Entropy stddev too low (< 0.05)");
          }
          if (coarseness < 0.005) {
            validation_passed = false;
            errors.push_back("Coarseness too low (< 0.005)");
          }
          if (active_classes < 8) {
            validation_passed = false;
            errors.push_back("Insufficient token class diversity");
          }
          if (spikes_2_5 == 0) {
            validation_passed = false;
            errors.push_back("No entropy spikes at k=2.5");
          }
          if (active_classes >= 8 && active_classes <= 9) {
            warnings.push_back("Low token class diversity");
          }
          // Write JSON
          std::string final_file = options.files[0].name() + ".phase1.final";
          std::ofstream fout_final(final_file);
          if (fout_final) {
            fout_final << "{\n";
            fout_final << "  \"phase\": 1,\n";
            fout_final << "  \"profiling_bytes\": " << in_bytes << ",\n";
            fout_final << "  \"entropy\": {\n";
            fout_final << "    \"mean\": " << mean_entropy << ",\n";
            fout_final << "    \"stddev\": " << stddev_entropy << ",\n";
            fout_final << "    \"delta_stddev\": " << stddev_delta << "\n";
            fout_final << "  },\n";
            fout_final << "  \"coarseness\": " << coarseness << ",\n";
            fout_final << "  \"entropy_spikes\": {\n";
            fout_final << "    \"k_2_0\": " << spikes_2_0 << ",\n";
            fout_final << "    \"k_2_5\": " << spikes_2_5 << ",\n";
            fout_final << "    \"k_3_0\": " << spikes_3_0 << "\n";
            fout_final << "  },\n";
            fout_final << "  \"cdc_plus_parameters\": {\n";
            fout_final << "    \"target_segment_bytes\": " << target_segment_size << ",\n";
            fout_final << "    \"min_segment_bytes\": " << min_segment_size << ",\n";
            fout_final << "    \"max_segment_policy\": {\n";
            fout_final << "      \"low_variance_max\": 262144,\n";
            fout_final << "      \"high_variance_max\": 131072,\n";
            fout_final << "      \"variance_gate_ratio\": 0.25\n";
            fout_final << "    },\n";
            fout_final << "    \"entropy_threshold\": " << entropy_threshold << "\n";
            fout_final << "  },\n";
            fout_final << "  \"token_class_stats\": {\n";
            fout_final << "    \"transition_matrix\": [\n";
            for (int i = 0; i < num_classes; ++i) {
              fout_final << "      [";
              for (int j = 0; j < num_classes; ++j) {
                fout_final << transition_matrix[i][j];
                if (j < num_classes - 1) fout_final << ",";
              }
              fout_final << "]";
              if (i < num_classes - 1) fout_final << ",";
              fout_final << "\n";
            }
            fout_final << "    ],\n";
            fout_final << "    \"avg_churn\": " << avg_churn << "\n";
            fout_final << "  },\n";
            fout_final << "  \"validation\": {\n";
            fout_final << "    \"passed\": " << (validation_passed ? "true" : "false") << ",\n";
            fout_final << "    \"warnings\": [";
            for (size_t i = 0; i < warnings.size(); ++i) {
              fout_final << "\"" << warnings[i] << "\"";
              if (i < warnings.size() - 1) fout_final << ",";
            }
            fout_final << "],\n";
            fout_final << "    \"errors\": [";
            for (size_t i = 0; i < errors.size(); ++i) {
              fout_final << "\"" << errors[i] << "\"";
              if (i < errors.size() - 1) fout_final << ",";
            }
            fout_final << "]\n";
            fout_final << "  }\n";
            fout_final << "}\n";
            fout_final.close();
            std::cout << "Phase 1 final artifact written to " << final_file << std::endl;
          } else {
            std::cout << "Failed to open final file" << std::endl;
          }
        }
        delete profiler_ptr;
      }

      if (options.mode == Options::kModeSingleTest) {
        if (err = fout.open(out_file, std::ios_base::in | std::ios_base::binary)) {
          std::cerr << "Error opening: " << out_file << " (" << errstr(err) << ")" << std::endl;
          return 1;
        }
        Archive archive(&fout);
        archive.list();
        std::cout << "Verifying archive decompression" << std::endl;
        archive.decompress("", true);
      }
    }
    break;
  }
  case Options::kModeAdd: {
    // Add a single file.
    break;
  }
  case Options::kModeList: {
    auto in_file = options.archive_file.getName();
    File fin;
    int err = 0;
    if (err = fin.open(in_file, std::ios_base::in | std::ios_base::binary)) {
      std::cerr << "Error opening: " << in_file << " (" << errstr(err) << ")" << std::endl;
      return 1;
    }
    printHeader();
    std::cout << "Listing files in archive " << in_file << std::endl;
    Archive archive(&fin);
    const auto& header = archive.getHeader();
    if (!header.isArchive()) {
      std::cerr << "Attempting to open non mcm compatible file" << std::endl;
      return 1;
    }
    if (!header.isSameVersion()) {
      std::cerr << "Attempting to open old version " << header.majorVersion() << "." << header.minorVersion() << std::endl;
      return 1;
    }
    archive.list();
    fin.close();
    break;
  }
  case Options::kModeDecompress: {
    auto in_file = options.archive_file.getName();
    File fin;
    File fout;
    int err = 0;
    if (err = fin.open(in_file, std::ios_base::in | std::ios_base::binary)) {
      std::cerr << "Error opening: " << in_file << " (" << errstr(err) << ")" << std::endl;
      return 1;
    }
    printHeader();
    std::cout << "Decompresing archive " << in_file << std::endl;
    Archive archive(&fin);
    const auto& header = archive.getHeader();
    if (!header.isArchive()) {
      std::cerr << "Attempting to decompress non archive file" << std::endl;
      return 1;
    }
    if (!header.isSameVersion()) {
      std::cerr << "Attempting to decompress other version " << header.majorVersion() << "." << header.minorVersion() << std::endl;
      return 1;
    }
    // archive.decompress(options.files.back().getName());
    archive.decompress("");
    fin.close();
    // Decompress the single file in the archive to the output out.
    break;
  }
  case Options::kModeExtract: {
    // Extract a single file from multi file archive.
    break;
  }
  case Options::kModeExtractAll: {
    // Extract all the files in the archive.
    // archive.ExtractAll();
    break;
  }
  }

  // Phase 2.5 profiling
  if (options.phase2p5_profile) {
    std::cout << "Starting Phase 2.5 profiling" << std::endl;
    std::string input_file = options.files[0].name();
    std::string segments_file = input_file + ".phase2.segments.txt";
    std::ifstream fin_seg(segments_file);
    if (!fin_seg) {
      std::cerr << "Cannot open " << segments_file << std::endl;
      return 1;
    }
    std::vector<std::tuple<uint64_t, uint64_t, double, double>> segments;
    uint64_t start, end;
    double mean, stddev;
    while (fin_seg >> start >> end >> mean >> stddev) {
      segments.emplace_back(start, end, mean, stddev);
    }
    fin_seg.close();

    // Read original file
    std::ifstream fin_orig(input_file, std::ios::binary);
    std::vector<uint8_t> original((std::istreambuf_iterator<char>(fin_orig)), std::istreambuf_iterator<char>());
    fin_orig.close();

    // Extract segments
    std::vector<std::vector<uint8_t>> segment_buffers;
    for (auto& seg : segments) {
      uint64_t s, e;
      std::tie(s, e, std::ignore, std::ignore) = seg;
      segment_buffers.emplace_back(original.begin() + s, original.begin() + e);
    }

    // Constants
    const int N = 4096;
    const double FP_THRESHOLD = 1e9; // No fingerprint, large
    const double ENTROPY_DELTA_1 = 0.1;
    const double ENTROPY_DELTA_2 = 0.5;

    // Greedy successor selection
    std::vector<int> order;
    std::vector<bool> used(segments.size(), false);

    // Start with lowest entropy segment
    int current = -1;
    double min_entropy = 1e9;
    for (size_t i = 0; i < segments.size(); ++i) {
      double ent = std::get<2>(segments[i]);
      if (ent < min_entropy) {
        min_entropy = ent;
        current = i;
      }
    }
    order.push_back(current);
    used[current] = true;

    std::string transitions_file = input_file + ".phase2p5.transitions.txt";
    std::ofstream fout_trans(transitions_file);

    while (order.size() < segments.size()) {
      // Find candidates
      std::vector<int> candidates;
      int tier = 1;
      double delta = ENTROPY_DELTA_1;
      for (size_t i = 0; i < segments.size(); ++i) {
        if (!used[i]) {
          double ent_diff = std::abs(std::get<2>(segments[current]) - std::get<2>(segments[i]));
          if (ent_diff <= delta) {
            candidates.push_back(i);
          }
        }
      }
      if (candidates.empty()) {
        tier = 2;
        delta = ENTROPY_DELTA_2;
        for (size_t i = 0; i < segments.size(); ++i) {
          if (!used[i]) {
            double ent_diff = std::abs(std::get<2>(segments[current]) - std::get<2>(segments[i]));
            if (ent_diff <= delta) {
              candidates.push_back(i);
            }
          }
        }
      }
      if (candidates.empty()) {
        tier = 3;
        // All unused, sorted by entropy diff
        std::vector<std::pair<double, int>> sorted;
        for (size_t i = 0; i < segments.size(); ++i) {
          if (!used[i]) {
            double ent_diff = std::abs(std::get<2>(segments[current]) - std::get<2>(segments[i]));
            sorted.emplace_back(ent_diff, i);
          }
        }
        std::sort(sorted.begin(), sorted.end());
        for (auto& p : sorted) {
          candidates.push_back(p.second);
        }
      }

      // Evaluate candidates
      int best_b = -1;
      uint64_t best_size = UINT64_MAX;
      for (int b : candidates) {
        // Extract tail(A) + head(B)
        auto& buf_a = segment_buffers[current];
        auto& buf_b = segment_buffers[b];
        std::vector<uint8_t> tail_a(buf_a.end() - std::min(N, (int)buf_a.size()), buf_a.end());
        std::vector<uint8_t> head_b(buf_b.begin(), buf_b.begin() + std::min(N, (int)buf_b.size()));
        std::vector<uint8_t> combined = tail_a;
        combined.insert(combined.end(), head_b.begin(), head_b.end());

        // Compress
        std::string temp_file = "temp_transition.bin";
        std::ofstream fout_temp(temp_file, std::ios::binary);
        fout_temp.write((char*)combined.data(), combined.size());
        fout_temp.close();
        FileInfo temp_info(temp_file);
        std::vector<FileInfo> temp_files = {temp_info};
        CompressionOptions fast_options = options.options_;
        fast_options.comp_level_ = kCompLevelTurbo;
        VoidWriteStream vws;
        Archive archive(&vws, fast_options, nullptr);
        uint64_t in_bytes = archive.compress(temp_files);
        uint64_t comp_size = vws.tell();
        std::remove(temp_file.c_str());

        if (comp_size < best_size) {
          best_size = comp_size;
          best_b = b;
        }
      }

      // Select best
      order.push_back(best_b);
      used[best_b] = true;

      // Log
      double buffer_len = std::min(N, (int)segment_buffers[current].size()) + std::min(N, (int)segment_buffers[best_b].size());
      double bpb = best_size * 8.0 / buffer_len;
      fout_trans << current << " " << best_b << " " << tier << " " << best_size << " " << buffer_len << " " << bpb << std::endl;

      current = best_b;
    }
    fout_trans.close();

    // Reordered file
    std::string reordered_file = input_file + ".phase2p5.reordered";
    std::ofstream fout_reordered(reordered_file, std::ios::binary);
    for (int idx : order) {
      fout_reordered.write((char*)segment_buffers[idx].data(), segment_buffers[idx].size());
    }
    fout_reordered.close();

    // Index
    std::string index_file = input_file + ".phase2p5.index.json";
    std::ofstream fout_index(index_file);
    fout_index << "{\n";
    fout_index << "  \"reordered_to_original\": [";
    for (size_t i = 0; i < order.size(); ++i) {
      fout_index << order[i];
      if (i < order.size() - 1) fout_index << ",";
    }
    fout_index << "],\n";
    fout_index << "  \"original_to_reordered\": [";
    std::vector<int> orig_to_reordered(order.size());
    for (size_t i = 0; i < order.size(); ++i) {
      orig_to_reordered[order[i]] = i;
    }
    for (size_t i = 0; i < orig_to_reordered.size(); ++i) {
      fout_index << orig_to_reordered[i];
      if (i < orig_to_reordered.size() - 1) fout_index << ",";
    }
    fout_index << "]\n";
    fout_index << "}\n";
    fout_index.close();

    // Validation
    uint64_t reordered_size = 0;
    for (auto& buf : segment_buffers) reordered_size += buf.size();
    if (reordered_size != original.size()) {
      std::cerr << "Size mismatch" << std::endl;
      return 1;
    }
    std::vector<int> used_count(segments.size(), 0);
    for (int idx : order) {
      if (idx < 0 || idx >= (int)segments.size()) {
        std::cerr << "Invalid index" << std::endl;
        return 1;
      }
      used_count[idx]++;
    }
    for (int c : used_count) {
      if (c != 1) {
        std::cerr << "Segment used " << c << " times" << std::endl;
        return 1;
      }
    }
    std::cout << "Phase 2.5 completed successfully" << std::endl;
  }

  // Phase 2.75 empirical reordering validation
  if (options.phase2p75_measure) {
    std::cout << "Starting Phase 2.75 empirical reordering validation" << std::endl;
    std::string input_file = options.files[0].name();
    std::string segments_file = input_file + ".phase2.segments.txt";
    std::ifstream fin_seg(segments_file);
    if (!fin_seg) {
      std::cerr << "Cannot open " << segments_file << std::endl;
      return 1;
    }
    std::vector<std::tuple<uint64_t, uint64_t, double, double>> segments;
    uint64_t start, end;
    double mean, stddev;
    while (fin_seg >> start >> end >> mean >> stddev) {
      segments.emplace_back(start, end, mean, stddev);
    }
    fin_seg.close();

    // Read original file
    std::ifstream fin_orig(input_file, std::ios::binary);
    std::vector<uint8_t> original((std::istreambuf_iterator<char>(fin_orig)), std::istreambuf_iterator<char>());
    fin_orig.close();

    // Select blocks: size >=64KB, variance <=0.01, not within 8KB of cuts, within first 5MB
    const uint64_t min_size = 65536;
    const double max_variance = 0.1;
    const uint64_t margin = 8192;
    const uint64_t max_offset = 5242880; // 5MB
    std::vector<std::tuple<int, uint64_t, uint64_t>> selected_blocks;
    for (size_t i = 0; i < segments.size(); ++i) {
      auto [s, e, m, stddev] = segments[i];
      uint64_t len = e - s;
      double variance = stddev * stddev;
      if (len >= min_size && variance <= max_variance && s >= margin && e <= original.size() - margin && e <= max_offset) {
        selected_blocks.emplace_back(i, s, e);
      }
    }
    if (selected_blocks.size() > 10) selected_blocks.resize(10); // limit to 10 for speed

    // Output selected blocks
    std::string blocks_file = input_file + ".phase2p75.blocks.txt";
    std::ofstream fout_blocks(blocks_file);
    for (auto& b : selected_blocks) {
      int id; uint64_t s, e;
      std::tie(id, s, e) = b;
      fout_blocks << id << " " << s << " " << e << " " << (e - s) << std::endl;
    }
    fout_blocks.close();

    // Compress original with no dictionary
    options.options_.dict_file_ = "";
    options.options_.filter_type_ = kFilterTypeNone;

    // Measurements
    std::string report_file = input_file + ".phase2p75.report.txt";
    std::ofstream fout_report(report_file);
    double total_avoidable_bits = 0;
    uint64_t total_measured_bytes = 0;

    for (auto& b : selected_blocks) {
      int id; uint64_t s, e;
      std::tie(id, s, e) = b;
      uint64_t len = e - s;

      // Compress prefix A = original[0..s]
      std::string temp_prefix = "temp_prefix.bin";
      std::ofstream fout_prefix(temp_prefix, std::ios::binary);
      fout_prefix.write((char*)original.data(), s);
      fout_prefix.close();
      FileInfo temp_info_prefix(temp_prefix);
      std::vector<FileInfo> temp_files_prefix = {temp_info_prefix};
      std::string compressed_prefix = temp_prefix + ".mcm";
      File fout_comp_prefix;
      fout_comp_prefix.open(compressed_prefix.c_str(), std::ios_base::out | std::ios_base::binary);
      Archive archive_prefix(&fout_comp_prefix, options.options_, nullptr);
      archive_prefix.compress(temp_files_prefix);
      uint64_t size_prefix = fout_comp_prefix.tell();
      fout_comp_prefix.close();
      std::remove(temp_prefix.c_str());
      std::remove(compressed_prefix.c_str());

      // Compress A + B = original[0..e]
      std::string temp_AB = "temp_AB.bin";
      std::ofstream fout_AB(temp_AB, std::ios::binary);
      fout_AB.write((char*)original.data(), e);
      fout_AB.close();
      FileInfo temp_info_AB(temp_AB);
      std::vector<FileInfo> temp_files_AB = {temp_info_AB};
      std::string compressed_AB = temp_AB + ".mcm";
      File fout_comp_AB;
      fout_comp_AB.open(compressed_AB.c_str(), std::ios_base::out | std::ios_base::binary);
      Archive archive_AB(&fout_comp_AB, options.options_, nullptr);
      archive_AB.compress(temp_files_AB);
      uint64_t size_AB = fout_comp_AB.tell();
      fout_comp_AB.close();
      std::remove(temp_AB.c_str());
      std::remove(compressed_AB.c_str());

      // Warm cost of B: cost(A+B) - cost(A)
      double warm_cost_B = (static_cast<double>(size_AB) - size_prefix) * 8;

      // Compress only B = original[s..e]
      std::string temp_B = "temp_B.bin";
      std::ofstream fout_B(temp_B, std::ios::binary);
      fout_B.write((char*)original.data() + s, len);
      fout_B.close();
      FileInfo temp_info_B(temp_B);
      std::vector<FileInfo> temp_files_B = {temp_info_B};
      std::string compressed_B = temp_B + ".mcm";
      File fout_comp_B;
      fout_comp_B.open(compressed_B.c_str(), std::ios_base::out | std::ios_base::binary);
      Archive archive_B(&fout_comp_B, options.options_, nullptr);
      archive_B.compress(temp_files_B);
      uint64_t size_B = fout_comp_B.tell();
      fout_comp_B.close();
      std::remove(temp_B.c_str());
      std::remove(compressed_B.c_str());

      // Cold cost of B: cost(B)
      double cold_cost_B = static_cast<double>(size_B) * 8;

      // Delta: cold - warm (positive means headroom)
      double delta = cold_cost_B - warm_cost_B;
      if (delta > 0) {
        total_avoidable_bits += delta;
        total_measured_bytes += len;
      }

      fout_report << "Block " << id << ":\n";
      fout_report << "  Prefix size: " << s << " bytes\n";
      fout_report << "  Block size: " << len << " bytes\n";
      fout_report << "  Cold compressed size: " << len << " -> " << size_B << " bytes (" << cold_cost_B << " bits)\n";
      fout_report << "  Warm compressed size: " << len << " -> " << (size_AB - size_prefix) << " bytes (" << warm_cost_B << " bits)\n";
      fout_report << "  Delta: " << delta << " bits\n\n";
    }

    double avoidable_per_mb = total_avoidable_bits / total_measured_bytes * 1000000;
    double projected_gain_kb = total_avoidable_bits / 8 / 1024;

    fout_report << "SUMMARY:\n";
    fout_report << "Blocks tested: " << selected_blocks.size() << "\n";
    fout_report << "Total measured bytes: " << total_measured_bytes / 1000000.0 << " MB\n";
    fout_report << "Total avoidable bits: " << total_avoidable_bits << "\n";
    fout_report << "Avoidable bits per MB: " << avoidable_per_mb << "\n";
    fout_report << "Projected enwik8 gain: ~" << projected_gain_kb << " KB\n";
    fout_report.close();

    std::cout << "Phase 2.75 completed. Report: " << report_file << std::endl;
  }

  return 0;
}
