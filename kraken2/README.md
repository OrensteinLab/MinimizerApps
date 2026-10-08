[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/kraken2/README.html)
[![European Galaxy server](https://img.shields.io/badge/usegalaxy-.eu-brightgreen?logo=data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABgAAAASCAYAAABB7B6eAAAABGdBTUEAALGPC/xhBQAAACBjSFJNAAB6JgAAgIQAAPoAAACA6AAAdTAAAOpgAAA6mAAAF3CculE8AAAACXBIWXMAAAsTAAALEwEAmpwYAAACC2lUWHRYTUw6Y29tLmFkb2JlLnhtcAAAAAAAPHg6eG1wbWV0YSB4bWxuczp4PSJhZG9iZTpuczptZXRhLyIgeDp4bXB0az0iWE1QIENvcmUgNS40LjAiPgogICA8cmRmOlJERiB4bWxuczpyZGY9Imh0dHA6Ly93d3cudzMub3JnLzE5OTkvMDIvMjItcmRmLXN5bnRheC1ucyMiPgogICAgICA8cmRmOkRlc2NyaXB0aW9uIHJkZjphYm91dD0iIgogICAgICAgICAgICB4bWxuczp0aWZmPSJodHRwOi8vbnMuYWRvYmUuY29tL3RpZmYvMS4wLyI+CiAgICAgICAgIDx0aWZmOlJlc29sdXRpb25Vbml0PjI8L3RpZmY6UmVzb2x1dGlvblVuaXQ+CiAgICAgICAgIDx0aWZmOkNvbXByZXNzaW9uPjE8L3RpZmY6Q29tcHJlc3Npb24+CiAgICAgICAgIDx0aWZmOk9yaWVudGF0aW9uPjE8L3RpZmY6T3JpZW50YXRpb24+CiAgICAgICAgIDx0aWZmOlBob3RvbWV0cmljSW50ZXJwcmV0YXRpb24+MjwvdGlmZjpQaG90b21ldHJpY0ludGVycHJldGF0aW9uPgogICAgICA8L3JkZjpEZXNjcmlwdGlvbj4KICAgPC9yZGY6UkRGPgo8L3g6eG1wbWV0YT4KD0UqkwAAAn9JREFUOBGlVEuLE0EQruqZiftwDz4QYT1IYM8eFkHFw/4HYX+GB3/B4l/YP+CP8OBNTwpCwFMQXAQPKtnsg5nJZpKdni6/6kzHvAYDFtRUT71f3UwAEbkLch9ogQxcBwRKMfAnM1/CBwgrbxkgPAYqlBOy1jfovlaPsEiWPROZmqmZKKzOYCJb/AbdYLso9/9B6GppBRqCrjSYYaquZq20EUKAzVpjo1FzWRDVrNay6C/HDxT92wXrAVCH3ASqq5VqEtv1WZ13Mdwf8LFyyKECNbgHHAObWhScf4Wnj9CbQpPzWYU3UFoX3qkhlG8AY2BTQt5/EA7qaEPQsgGLWied0A8VKrHAsCC1eJ6EFoUd1v6GoPOaRAtDPViUr/wPzkIFV9AaAZGtYB568VyJfijV+ZBzlVZJ3W7XHB2RESGe4opXIGzRTdjcAupOK09RA6kzr1NTrTj7V1ugM4VgPGWEw+e39CxO6JUw5XhhKihmaDacU2GiR0Ohcc4cZ+Kq3AjlEnEeRSazLs6/9b/kh4eTC+hngE3QQD7Yyclxsrf3cpxsPXn+cFdenF9aqlBXMXaDiEyfyfawBz2RqC/O9WF1ysacOpytlUSoqNrtfbS642+4D4CS9V3xb4u8P/ACI4O810efRu6KsC0QnjHJGaq4IOGUjWTo/YDZDB3xSIxcGyNlWcTucb4T3in/3IaueNrZyX0lGOrWndstOr+w21UlVFokILjJLFhPukbVY8OmwNQ3nZgNJNmKDccusSb4UIe+gtkI+9/bSLJDjqn763f5CQ5TLApmICkqwR0QnUPKZFIUnoozWcQuRbC0Km02knj0tPYx63furGs3x/iPnz83zJDVNtdP3QAAAABJRU5ErkJggg==)](https://usegalaxy.eu/root?tool_id=kraken2)

Kraken 2
========

The second version of the Kraken taxonomic sequence classification system

Please refer to the Operating Manual (in `docs/MANUAL.html`) for details on
how to use Kraken 2.

For minimizer / mod-minimizer sketching benchmarks and library setup, see
[`../benchmark/kraken2/README.md`](../benchmark/kraken2/README.md). On a fresh machine, run
[`../benchmark/kraken2/setup.sh`](../benchmark/kraken2/setup.sh) (host build deps, Miniconda if needed, and a
conda env with [bioconda `kraken2`](https://anaconda.org/channels/bioconda/packages/kraken2/overview)).

Custom minimizer
----------------

This fork adds pluggable minimizer orders (`--order xor|xor-pcg|greedymini`)
and a sketching choice (`--sketching minimizer|mod-minimizer`; on DNA the
mod-minimizer is strand-symmetric, see `src/mod_minimizer.h`; with
`--greedymini-file FILE` for the greedymini order) on top of upstream
Kraken 2's own default order (`--sketching minimizer --order xor`), including a shared GreedyMini-based
low-density order also used by this monorepo's
[KMC fork](../kmc/README.md#custom-minimizer). See
[`../orders/README.md`](../orders/README.md) for the option/CLI reference
and shared implementation.

With `--order greedymini` and `--sketching minimizer`, a table with `gm_k < l` is
extended by its **centre** (rank of the middle `gm_k` bases, then tie-breaks),
which keeps the table's designed density (`w5_k15_canonical`: 1.466 at every
valid `l`).
`l - gm_k` must be even (`w5_k15_canonical`: `l` = 15, 17, …, 31). The
GreedyMini order is DNA only: protein databases (`--protein`) refuse it.
`--gm-projection at1-cg0` flips its DNA->binary projection (A,T -> 1, C,G -> 0;
default `at0-cg1`), as in the sshash and minimap2 forks; it is stored in the
database. Databases
by other Kraken 2 builds (any `db_version` other than 2, including upstream
databases, and databases built before the opts.k2d layout change: wrong size or
missing magic) are refused by `classify` and `dump_table` and must be rebuilt. Details, numbers and the commands to reproduce them:
[`minimizers/README.md`](minimizers/README.md).

Publications
------------

For additional implementation details and guidance on using Kraken 2, you can also refer to:
- The [paper describing the first version of Kraken](https://dx.doi.org/10.1186/gb-2014-15-3-r46).  Note that Kraken 2 is a rewrite of Kraken 1 and is not backwards compatible with Kraken 1.
- The [paper describing KrakenUniq](https://dx.doi.org/10.1186/s13059-018-1568-0).  KrakenUniq's HLL-based functionality is incorporated in Kraken 2 when using the `--report-minimizer-data` flag described in the Operating Manual.
- The [Kraken 2 paper](https://dx.doi.org/10.1186/s13059-019-1891-0).
- A [protocol paper describing use of Kraken 2](https://dx.doi.org/10.1038/s41596-022-00738-y).
