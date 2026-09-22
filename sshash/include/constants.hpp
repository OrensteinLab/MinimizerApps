#pragma once

namespace sshash::constants {

constexpr uint64_t invalid_uint64 = uint64_t(-1);
constexpr uint64_t default_ram_limit_in_GiB = 8;
constexpr uint64_t seed = 1;

/* for PTHash */
constexpr double lambda = 5.0;
constexpr uint64_t avg_partition_size = 3000000;

constexpr uint64_t min_l = 6;
/* min_l/max_l are NOT a bound on the minimizer length m (a prior
   misconception this fork briefly repeated -- see git history): they size
   sshash's skew-index bucket-size partitioning
   (src/builder/build_sparse_and_skew_index.cpp), which packs a
   `partition_id` into exactly 3 bits (`heavy_load_buckets_size << 3 |
   partition_id`), so num_partitions must never exceed 8. Keeping
   max_l - min_l + 1 <= 8 guarantees that arithmetically, for *any* dataset's
   bucket-size distribution -- not just ones tested so far -- so do not
   raise max_l without also widening that 3-bit field (and every place that
   assumes it) to match. */
constexpr uint64_t max_l = 13;
static_assert(min_l < max_l);
static_assert(max_l - min_l + 1 <= 8);

static const std::string default_tmp_dirname(".");
constexpr int forward_orientation = 1;
constexpr int backward_orientation = -1;

namespace current_version_number {
constexpr uint8_t x = 6;
constexpr uint8_t y = 0;
constexpr uint8_t z = 0;
}  // namespace current_version_number

}  // namespace sshash::constants
