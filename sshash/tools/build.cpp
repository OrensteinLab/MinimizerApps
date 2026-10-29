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
    parser.add("bucket_stats",
               "Also print the bucket-size and super-k-mer histograms after construction.",
               "--bucket-stats", false, true);
    parser.add("scheme",
               "Minimizer scheme: native (default), xor, xor-pcg-hash or greedymini; or the "
               "strand-symmetric mod-minimizer with t-mers ranked by xor (mod-minimizer-xor) or "
               "GreedyMini (mod-minimizer-greedymini), with --mod-t. See orders/README.md.",
               "--scheme", false);
    parser.add("greedy_mini_filename",
               "GreedyMini .gm rank table (required with --scheme greedymini and --scheme mod-minimizer-greedymini; "
               "for mod-minimizer-greedymini a table of t-mers).",
               "--greedy-mini-file", false);
    parser.add("mod_t",
               "t-mer length of the mod-minimizer schemes (0 < t < m; the usual choice is "
               "4 + (m - 4) mod (k - m + 1)).",
               "--mod-t", false);
    parser.add("gm_projection",
               "DNA->binary projection polarity for GreedyMini: "
               "at0-cg1 (default; A,T->0 and C,G->1) or at1-cg0 (flipped).",
               "--gm-projection", false);
    parser.add("regular_parsing",
               "Force two-probe 'regular parsing' (default is single-probe canonical parsing).",
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
    build_config.bucket_stats = parser.get<bool>("bucket_stats");
    if (parser.parsed("tmp_dirname")) {
        build_config.tmp_dirname = parser.get<std::string>("tmp_dirname");
        essentials::create_directory(build_config.tmp_dirname);
    }
    if (parser.get<uint64_t>("RAM")) {
        build_config.ram_limit_in_GiB = parser.get<uint64_t>("RAM");
    }
    if (parser.parsed("t")) build_config.num_threads = parser.get<uint64_t>("t");

    if (parser.parsed("scheme")) {
        build_config.scheme = parse_custom_scheme_mode(parser.get<std::string>("scheme"));
    }
    if (parser.parsed("greedy_mini_filename")) {
        build_config.greedy_mini_filename = parser.get<std::string>("greedy_mini_filename");
    }
    if (parser.parsed("gm_projection")) {
        build_config.gm_proj = parse_gm_projection(parser.get<std::string>("gm_projection"));
        if (build_config.gm_proj != gm_projection::AT0_CG1 and
            !is_gm_mode(build_config.scheme)) {
            std::cerr << "Error: --gm-projection at1-cg0 requires a GreedyMini scheme" << std::endl;
            return 1;
        }
    }
    if (is_gm_mode(build_config.scheme) && build_config.greedy_mini_filename.empty()) {
        std::cerr << "Error: --scheme " << custom_scheme_mode_name(build_config.scheme)
                  << " requires --greedy-mini-file" << std::endl;
        return 1;
    }
    if (parser.parsed("mod_t")) build_config.mod_t = parser.get<uint64_t>("mod_t");
    if (is_mod_mode(build_config.scheme) != (build_config.mod_t != 0)) {
        std::cerr << "Error: --mod-t is required by, and only valid with, the mod-minimizer "
                     "schemes (mod-minimizer-xor, mod-minimizer-greedymini)"
                  << std::endl;
        return 1;
    }
    build_config.regular_parsing = parser.get<bool>("regular_parsing");
    if (build_config.regular_parsing and is_mod_mode(build_config.scheme)) {
        std::cerr << "Error: --regular does not apply to the mod-minimizer schemes" << std::endl;
        return 1;
    }

    if (is_mod_mode(build_config.scheme) and
        !mod_rule::t_is_canonical(int(m), int(k - m + 1), int(build_config.mod_t))) {
        std::cerr << "Error: the mod-minimizer needs --mod-t == m (mod k - m + 1) (the centre tie "
                     "rule; e.g. t = 4 + (m - 4) mod (k - m + 1))"
                  << std::endl;
        return 1;
    }

    if (is_gm_mode(build_config.scheme) or is_mod_mode(build_config.scheme)) {
        /* fail early with a clear message (e.g. an odd m - gm_k,
           or a table that is not mirror-symmetric) instead of an uncaught exception */
        try {
            hasher_type probe(build_config.seed);
            probe.configure(build_config.scheme, m, build_config.greedy_mini_filename,
                            build_config.gm_proj, build_config.mod_t);
        } catch (std::exception const& e) {
            std::cerr << "Error: " << e.what() << std::endl;
            return 1;
        }
    }

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
