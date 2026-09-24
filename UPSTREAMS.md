# Tracking upstream

This repo folds four forked tools into one monorepo, each living in its own
top-level folder with lab-specific work on top of an unmodified import of
the real upstream project. `origin` points at this monorepo; each folder
also has its own remote pointing at the real upstream project:

| Folder | Upstream | Remote |
|---|---|---|
| `kraken2/` | [DerrickWood/kraken2](https://github.com/DerrickWood/kraken2) | `upstream-kraken2` |
| `kmc/` | [refresh-bio/KMC](https://github.com/refresh-bio/KMC) | `upstream-kmc` |
| `sshash/` | [jermp/sshash](https://github.com/jermp/sshash) | `upstream-sshash` |
| `minimap2/` | [lh3/minimap2](https://github.com/lh3/minimap2), pinned to release tag `v2.31` | `upstream-minimap2` |

## Pulling in an upstream update

Each folder was merged in additively (`git subtree`-style: `read-tree
--prefix=<dir>/` plus an `-s ours` merge, since the `git-subtree` contrib
command isn't installed on every machine this repo is cloned on), so no
existing branch was rewritten. Updating from upstream uses the same
manual recipe. Replace `<dir>` with `kraken2`, `kmc`, `sshash`, or `minimap2`
and `<upstream>` with the matching remote from the table above:

1. `git fetch <upstream> <branch>` (e.g. `upstream-kraken2 master`)
2. `git merge -s ours --no-commit --allow-unrelated-histories <upstream>/<branch>`
3. `git read-tree --prefix=<dir>/ -u <upstream>/<branch>`
4. `git commit -m "Merge <dir> update from <upstream>/<branch>"`
5. Re-apply (rebase or manually reconcile) any lab-specific commits that touch files upstream also changed, same as keeping a normal fork in sync, just scoped to one subfolder of a monorepo instead of a whole repo.

`minimap2/` is pinned to a release tag rather than tracking a moving branch
(reproducibility — a custom-order benchmark shouldn't silently shift under
an unrelated upstream commit): step 1 is `git fetch upstream-minimap2 tag
vX.Y` and steps 2-3 use the tag (e.g. `v2.31`) directly in place of
`<upstream>/<branch>`. Bumping the pin later is the same recipe against the
new tag.

`sshash/` additionally vendors a real git submodule (`sshash/external/pthash`,
registered in the root `.gitmodules` since git only honors that file at the
repo root). If the pinned pthash commit moves upstream, also run
`git submodule update --init --recursive -- sshash/external/pthash` after
step 3 to pull in the new submodule content — the tree/blob content is what
`read-tree` updates, but a submodule reference is a gitlink pointing at a
commit, and that commit's content isn't fetched automatically.

Step 2's merge commit brings the new upstream commit in as a parent, so
`git log <merge-commit>^2` shows exactly what changed upstream, and step 3
updates `<dir>/`'s tree to match it.

If a machine has the real `git-subtree` command installed, `git subtree
pull --prefix=<dir> <upstream> <branch>` does the same thing in one step.
