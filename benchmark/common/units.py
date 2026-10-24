"""Memory, space and I/O units shared by all exporters: always decimal GB.

1 GB = 10^9 bytes (not GiB). /usr/bin/time -v reports peak RSS in KiB
(1024 bytes), so convert via bytes. Column names end in `_gb`.
"""

BYTES_PER_GB = 10 ** 9
BYTES_PER_KIB = 1024


def gb_from_bytes(n, digits=4):
    return None if n is None else round(float(n) / BYTES_PER_GB, digits)


def gb_from_kib(kib, digits=4):
    """Peak RSS as reported by GNU time (KiB) -> decimal GB."""
    return None if kib is None else round(float(kib) * BYTES_PER_KIB / BYTES_PER_GB, digits)


BYTES_PER_BLOCK = 512


def gb_from_blocks(blocks, digits=4):
    """GNU time's 'File system inputs/outputs' (512-byte blocks on Linux) -> decimal GB."""
    return None if blocks is None else round(float(blocks) * BYTES_PER_BLOCK / BYTES_PER_GB, digits)
