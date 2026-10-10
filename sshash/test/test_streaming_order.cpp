// streaming_query<Dict, Mode> must be instantiated with the order of the dictionary: a mismatched
// instantiation throws std::logic_error, a matching one does not.
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

template <custom_order_mode Mode>
static bool constructs(dictionary_type const& dict) {
    try {
        streaming_query<dictionary_type, Mode> q(&dict);
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
    constexpr auto NATIVE = custom_order_mode::NATIVE;
    constexpr auto XOR_HASH = custom_order_mode::XOR_HASH;
    constexpr auto XOR_PCG_HASH = custom_order_mode::XOR_PCG_HASH;
    constexpr auto GREEDYMINI = custom_order_mode::GREEDYMINI;
    expect(constructs<NATIVE>(native), true, "native + native instantiation constructs");
    expect(constructs<XOR_HASH>(native), false, "native + xor instantiation throws");
    expect(constructs<XOR_PCG_HASH>(native), false, "native + xor-pcg-hash instantiation throws");
    expect(constructs<GREEDYMINI>(native), false, "native + greedymini instantiation throws");
    expect(constructs<GREEDYMINI>(gm), true, "greedymini + greedymini instantiation constructs");
    expect(constructs<NATIVE>(gm), false, "greedymini + native instantiation throws");
    expect(constructs<XOR_HASH>(gm), false, "greedymini + xor instantiation throws");
    expect(constructs<XOR_PCG_HASH>(gm), false, "greedymini + xor-pcg-hash instantiation throws");
    if (failures) return 1;
    std::cout << "OK" << std::endl;
    return 0;
}
