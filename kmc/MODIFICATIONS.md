# Modifications to KMC

This folder is KMC 3.2.4 ([refresh-bio/KMC](https://github.com/refresh-bio/KMC),
commit `751ef36a3c1ccc6dda664f529ad218dc51d76f55`, imported in commit `3e38ae98`),
modified by the Orenstein Lab in 2026. Like the original, it is distributed under
the GNU General Public License, version 3 ([LICENSE](LICENSE)).

The modifications add custom signature orders (`xor`, `xor-pcg-hash` and GreedyMini
tables, see [`../orders/README.md`](../orders/README.md)) and statistics of the
superwindows for the benchmark. The order code shared with the other tools of this
repository is in [`../common/minimizer-schemes/`](../common/minimizer-schemes/) (MIT,
GPL-compatible). `git diff 3e38ae98 -- kmc` shows every change.

| File | Change |
|---|---|
| `kmc_api/custom_order.cpp`, `kmc_api/custom_order.h` | added: the custom signature orders |
| `kmc_CLI/kmc.cpp` | command-line options `-or`, `-co`, `-os`, `-ol`, `-ode`, `-oat1` |
| `kmc_api/mmer.cpp`, `kmc_api/mmer.h` | ranking signatures by a custom order |
| `kmc_api/kmc_file.cpp`, `kmc_api/kmc_file.h` | reading the custom order stored with a database |
| `kmc_core/kb_completer.cpp`, `kmc_core/kb_completer.h` | storing the custom order with the database |
| `kmc_core/kmc.h`, `kmc_core/params.h`, `kmc_core/s_mapper.h` | passing the custom order through a run |
| `kmc_core/kmc_runner.cpp`, `kmc_core/kmc_runner.h` | API for custom signature orders |
| `kmc_core/splitter.cpp`, `kmc_core/splitter.h` | custom orders in the splitter; superwindow statistics |
| `Makefile` | dynamic link by default, Python autodetect, raduls flags; the new sources |
| `*.vcxproj` | the new sources in the Visual Studio projects |
| `README.md`, `.gitignore` | this notice; build outputs |
