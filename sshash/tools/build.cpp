#include <tuple>

using namespace sshash;

int build(int argc, char** argv) {
    cmd_line_parser::parser parser(argc, argv);

    /* Required arguments. */
    parser.add(
        "input_filename",
        "Must be a FASTA file (.fa/fasta extension) or a CUTTLEFISH file (.cfseg extension, "
        "produced by CUTTLEFISH with option -f 3) "
        "compressed with gzip (.gz) or not:\n"
        "\t- without duplicate nor invalid kmers\n"
        "\t- one DNA sequence per line.\n"
        "\tFor example, it could be the de Bruijn graph topology output by BCALM or CUTTLEFISH.",
        "-i", true);
    parser.add("k", "K-mer length (must be <= " + std::to_string(default_kmer_t::max_k) + ").",
               "-k", true);
    parser.add("m", "Minimizer length (must be < k).", "-m", true);

    /* Optional arguments. */
    parser.add("seed",
               "Seed for construction (default is " + std::to_string(constants::seed) + ").", "-s",
               false);
    parser.add("t", "Number of threads (default is 1).", "-t", false);
    parser.add("lambda",
               "A (floating point) constant that trades construction speed for space effectiveness "
               "of minimal perfect hashing. "
               "A reasonable value lies between 3.0 and 10.0 (default is " +
                   std::to_string(constants::lambda) + ").",
               "-a", false);
    parser.add("output_filename", "Output file name where the data structure will be serialized.",
               "-o", false);
    parser.add(
        "tmp_dirname",
        "Temporary directory used for construction in external memory. Default is directory '" +
            constants::default_tmp_dirname + "'.",
        "-d", false);
    parser.add("RAM",
               "RAM limit in GiB. Default value is " +
                   std::to_string(constants::default_ram_limit_in_GiB) + " GiB.",
               "-g", false);
    parser.add("weighted", "Also store the weights in compressed format.", "--weighted", false,
               true);
    parser.add("check", "Check correctness after construction.", "--check", false, true);
    parser.add("verbose", "Verbose output during construction.", "--verbose", false, true);
    parser.add("order",
               "Minimizer order: native (default), xor-hash, xor-pcg-hash, greedymini, or one "
               "of the 4 mod-minimizer orders (mod-minimizer-native-hash, "
               "mod-minimizer-xor-hash, mod-minimizer-pcg-xor-hash, mod-minimizer-greedymini). "
               "See orders/README.md.",
               "-O", false);
    parser.add("greedy_mini_filename",
               "GreedyMini .gm rank table (required with -O greedymini or "
               "-O mod-minimizer-greedymini).",
               "--greedy-mini-file", false);
    parser.add("gm_alphabet", "GreedyMini table alphabet: binary (default) or dna.",
               "--gm-alphabet", false);
    parser.add("t_length",
               "t-mer length for mod-minimizer's inner order O_t (required, and only "
               "meaningful, with a mod-minimizer-* order; must be < m).",
               "--t-length", false);
    parser.add("regular_parsing",
               "Force two-probe 'regular parsing' for a non-mod-minimizer order "
               "(default is single-probe canonical parsing). No-op for "
               "mod-minimizer-* orders, which are always regular parsing.",
               "--regular", false, true);

    if (!parser.parse()) return 0;

    auto input_filename = parser.get<std::string>("input_filename");
    auto k = parser.get<uint64_t>("k");
    auto m = parser.get<uint64_t>("m");

    build_configuration build_config;
    build_config.k = k;
    build_config.m = m;

    if (parser.parsed("seed")) build_config.seed = parser.get<uint64_t>("seed");
    if (parser.parsed("lambda")) build_config.lambda = parser.get<double>("lambda");
    build_config.weighted = parser.get<bool>("weighted");
    build_config.verbose = parser.get<bool>("verbose");
    if (parser.parsed("tmp_dirname")) {
        build_config.tmp_dirname = parser.get<std::string>("tmp_dirname");
        essentials::create_directory(build_config.tmp_dirname);
    }
    if (parser.get<uint64_t>("RAM")) {
        build_config.ram_limit_in_GiB = parser.get<uint64_t>("RAM");
    }
    if (parser.parsed("t")) build_config.num_threads = parser.get<uint64_t>("t");

    if (parser.parsed("order")) {
        std::tie(build_config.order, build_config.mod_minimizer) =
            parse_order_name(parser.get<std::string>("order"));
    }
    if (parser.parsed("greedy_mini_filename")) {
        build_config.greedy_mini_filename = parser.get<std::string>("greedy_mini_filename");
    }
    if (parser.parsed("gm_alphabet")) {
        build_config.gm_alphabet = parse_gm_alphabet(parser.get<std::string>("gm_alphabet"));
    }
    if (build_config.order == custom_order_mode::GREEDYMINI &&
        build_config.greedy_mini_filename.empty()) {
        std::cerr << "Error: -O greedymini requires --greedy-mini-file" << std::endl;
        return 1;
    }
    if (build_config.mod_minimizer) {
        if (!parser.parsed("t_length")) {
            std::cerr << "Error: a mod-minimizer-* order requires --t-length" << std::endl;
            return 1;
        }
        build_config.t = parser.get<uint64_t>("t_length");
        if (build_config.t == 0 or build_config.t >= m) {
            std::cerr << "Error: --t-length must be > 0 and < m" << std::endl;
            return 1;
        }
    } else if (parser.parsed("t_length")) {
        std::cerr << "Error: --t-length is only meaningful with a mod-minimizer-* order"
                  << std::endl;
        return 1;
    }
    build_config.regular_parsing = parser.get<bool>("regular_parsing");

    // build_config.print();

    essentials::logger("building data structure...");
    dictionary_type dict;
    dict.build(input_filename, build_config);

    bool check = parser.get<bool>("check");
    if (check) {
        check_correctness_lookup_access(dict, input_filename);
        check_correctness_navigational_kmer_query(dict, input_filename);
        check_correctness_navigational_string_query(dict);
        if (build_config.weighted) check_correctness_weights(dict, input_filename);
        check_correctness_kmer_iterator(dict);
        check_correctness_string_iterator(dict);
    }

    if (parser.parsed("output_filename")) {
        auto output_filename = parser.get<std::string>("output_filename");
        essentials::logger("saving data structure to disk...");
        essentials::save(dict, output_filename.c_str());
        essentials::logger("DONE");
    }

    return 0;
}
