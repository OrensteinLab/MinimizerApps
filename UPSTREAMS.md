# Tracking upstream

This repo folds four tools into one monorepo. Each lives in its own top-level
folder, starts from one **import commit** that is a byte-exact copy of an
upstream commit, and is followed by lab-specific commits (orders, modifications,
benchmark harness). The import commits are plain single-parent commits; no
upstream history is carried in this repo.

## Pinned versions

| Folder | Upstream | Version | Upstream commit | Date | Remote |
|---|---|---|---|---|---|
| `kraken2/` | [DerrickWood/kraken2](https://github.com/DerrickWood/kraken2) | tag `2.17.2` (`master`; `VERSION` file reads `2.17.1`) | `8c190b1b668825935dbf6dee5f969227dc8269bb` | 2026-09-08 | `upstream-kraken2` |
| `kmc/` | [refresh-bio/KMC](https://github.com/refresh-bio/KMC) | KMC 3.2.4 (`master`, untagged) | `751ef36a3c1ccc6dda664f529ad218dc51d76f55` | 2026-07-07 | `upstream-kmc` |
| `sshash/` | [jermp/sshash](https://github.com/jermp/sshash) | tag `v6.0.0` (latest release) | `094ec2f342859e7713ceb82cbd0d0426d71a66a1` | 2026-08-29 | `upstream-sshash` |
| `minimap2/` | [lh3/minimap2](https://github.com/lh3/minimap2) | release `v2.31` (r1302) | `3c28777e7e2dcc90f825de1b9f17a89cca7d4452` | 2026-05-19 | `upstream-minimap2` |

Submodules (gitlinks identical to the ones the pinned upstream commits use;
registered in the root `.gitmodules`, since git only honors that file at the
repo root):

| Path | URL | Commit |
|---|---|---|
| `kmc/3rd_party/cloudflare` | https://github.com/refresh-bio/dependencies-zlib | `3cccc9613d79781a8486f9bae1ff38103963318a` |
| `sshash/external/pthash` | https://github.com/jermp/pthash.git | `e04a1920ffeae9e7d876acd0362cab79605f7af3` (`v2.0.0-30-ge04a192`) |

The only difference between an import commit and its upstream tree is that
KMC's own root `.gitmodules` is not carried into `kmc/`; its entry lives in the
monorepo root `.gitmodules` instead.

Each `Import ...` commit message repeats the version, upstream URL and SHA.

## Checking an import

```
git remote add upstream-<app> <url>   # once; see the table
git fetch upstream-<app>
git rev-parse <import-commit>:<app>   # tree of the import
git rev-parse <upstream-sha>^{tree}   # must be equal (KMC: apart from .gitmodules)
```

## Pulling in an upstream update

1. `git fetch <remote> <branch-or-tag>`
2. Replace `<dir>/` with the new upstream tree in its own commit, so the lab
   changes stay separable:
   `git rm -r -q <dir>` and `git read-tree --prefix=<dir>/ -u <new-upstream-sha>`,
   then `git commit -m "Import <app> <version> @ <sha>"` and update the table above.
3. Re-apply the lab-specific changes on top (they are the later commits that
   touch `<dir>/`: orders wrappers, modifications, harness). Expect conflicts
   in the files marked as modifications; reconcile them as with any fork, scoped
   to one subfolder.
4. If a submodule pin moved (`sshash/external/pthash`,
   `kmc/3rd_party/cloudflare`), run
   `git submodule update --init --recursive -- <path>`; `read-tree` only
   updates the gitlink, not the submodule content.

`git subtree pull --prefix=<dir> <remote> <branch>` does the same in one step
where `git-subtree` is installed.

`minimap2/` is pinned to a release tag rather than a moving branch so that a
custom-order benchmark does not shift under an unrelated upstream commit;
bumping the pin is the recipe above against the new tag.
