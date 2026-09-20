# Tracking upstream

This repo folds two forked tools into one monorepo, each living in its own
top-level folder with lab-specific work on top of an unmodified import of
the real upstream project. `origin` points at this monorepo; each folder
also has its own remote pointing at the real upstream project:

| Folder | Upstream | Remote |
|---|---|---|
| `kraken2/` | [DerrickWood/kraken2](https://github.com/DerrickWood/kraken2) | `upstream-kraken2` |
| `kmc/` | [refresh-bio/KMC](https://github.com/refresh-bio/KMC) | `upstream-kmc` |

## Pulling in an upstream update

Both folders were merged in additively (`git subtree`-style: `read-tree
--prefix=<dir>/` plus an `-s ours` merge, since the `git-subtree` contrib
command isn't installed on every machine this repo is cloned on), so no
existing branch was rewritten. Updating from upstream uses the same
manual recipe. Replace `<dir>` with `kraken2` or `kmc` and `<upstream>`
with the matching remote from the table above:

1. `git fetch <upstream> <branch>` (e.g. `upstream-kraken2 master`)
2. `git merge -s ours --no-commit --allow-unrelated-histories <upstream>/<branch>`
3. `git read-tree --prefix=<dir>/ -u <upstream>/<branch>`
4. `git commit -m "Merge <dir> update from <upstream>/<branch>"`
5. Re-apply (rebase or manually reconcile) any lab-specific commits that touch files upstream also changed, same as keeping a normal fork in sync, just scoped to one subfolder of a monorepo instead of a whole repo.

Step 2's merge commit brings the new upstream commit in as a parent, so
`git log <merge-commit>^2` shows exactly what changed upstream, and step 3
updates `<dir>/`'s tree to match it.

If a machine has the real `git-subtree` command installed, `git subtree
pull --prefix=<dir> <upstream> <branch>` does the same thing in one step.
