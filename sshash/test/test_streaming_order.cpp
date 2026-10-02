// streaming_query<Dict, StrandInvariant> must be instantiated with the iterator kind matching the
// dictionary's order: a mismatched instantiation throws std::logic_error, a matching one does not.
// Usage: test_streaming_order INPUT.fa[.gz] GREEDYMINI_TABLE.gm
#include <iostream>

#include "include/dictionary_types.hpp"
#include "include/streaming_query.hpp"

/* the library is compiled into each executable that needs it (see tools/sshash.cpp) */
#include "src/builder/build.cpp"
#include "src/dictionary.cpp"
#include "src/query.cpp"
#include "src/info.cpp"

using namespace sshash;

template <bool StrandInvariant>
static bool constructs(dictionary_type const& dict) {
    try {
        streaming_query<dictionary_type, StrandInvariant> q(&dict);
        return true;
    } catch (std::logic_error const&) {
        return false;
    }
}

int main(int argc, char** argv) {
    if (argc < 3) {
        std::cerr << "usage: " << argv[0] << " INPUT.fa[.gz] GREEDYMINI_TABLE.gm" << std::endl;
        return 1;
    }
    build_configuration config;
    config.k = 31;
    config.m = 21;
    config.verbose = false;
    config.num_threads = 2;

    dictionary_type native;
    config.order = custom_order_mode::NATIVE;
    native.build(argv[1], config);

    dictionary_type gm;
    config.order = custom_order_mode::GREEDYMINI;
    config.greedy_mini_filename = argv[2];
    gm.build(argv[1], config);

    int failures = 0;
    auto expect = [&](bool got, bool want, const char* what) {
        if (got != want) {
            std::cerr << "FAILED: " << what << std::endl;
            ++failures;
        }
    };
    expect(constructs<false>(native), true, "native + classical iterator constructs");
    expect(constructs<true>(native), false, "native + strand-invariant iterator throws");
    expect(constructs<true>(gm), true, "greedymini + strand-invariant iterator constructs");
    expect(constructs<false>(gm), false, "greedymini + classical iterator throws");
    if (failures) return 1;
    std::cout << "OK" << std::endl;
    return 0;
}
