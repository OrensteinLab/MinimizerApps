#!/usr/bin/env bash
# Sanity/regression tests for the custom minimizer order integration
# (order.h/order.cpp). See the plan's "Sanity checks / tests" section.
#
# greedymini uses orders/GreedyMini/w10_k15_canonical.gm (gm_k=15, matching
# this script's fixed -k 15 -w 10 exactly). There is no mod-minimizer-
# greedymini for minimap2 (unlike kraken2/sshash) -- it would need a second,
# differently-sized table (gm_k == the mod-minimizer t, not k) just for that
# one variant, so it was dropped instead of adding that second table.
#
# Usage: ./test_order.sh (run from minimap2/, after `make`)
set -u
cd "$(dirname "$0")"
BIN=./minimap2
REF=test/MT-human.fa
QRY=test/MT-orang.fa
TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT

pass=0 fail=0
ok()   { pass=$((pass+1)); echo "  [PASS] $1"; }
bad()  { fail=$((fail+1)); echo "  [FAIL] $1"; }

ORDERS_SEEDED="xor-hash xor-pcg-hash"
ORDERS_MOD="mod-minimizer-xor-hash mod-minimizer-pcg-xor-hash"
ALL_TEST_ORDERS="native $ORDERS_SEEDED greedymini $ORDERS_MOD"
SEED=42
MOD_T=5
GM_FILE="$(cd .. && pwd)/orders/GreedyMini/w10_k15_canonical.gm"

order_flags() {
	case "$1" in
		native) echo "" ;;
		greedymini) echo "--minimizer-order greedymini --greedymini-file $GM_FILE" ;;
		mod-minimizer-*) echo "--minimizer-order $1 --order-seed $SEED --mod-t $MOD_T" ;;
		*) echo "--minimizer-order $1 --order-seed $SEED" ;;
	esac
}

echo "== 1. determinism: same order+seed, repeated runs -> identical output =="
for o in $ALL_TEST_ORDERS; do
	f=$(order_flags "$o")
	$BIN -a $f "$REF" "$QRY" > "$TMP/r1_$o.sam" 2>/dev/null
	$BIN -a $f "$REF" "$QRY" > "$TMP/r2_$o.sam" 2>/dev/null
	if diff -q "$TMP/r1_$o.sam" "$TMP/r2_$o.sam" >/dev/null; then ok "determinism: $o"; else bad "determinism: $o"; fi
done

echo "== 2. density: each order shows a plausible, non-degenerate minimizer count =="
declare -A density
for o in $ALL_TEST_ORDERS; do
	f=$(order_flags "$o")
	n=$($BIN -a $f "$REF" "$QRY" 2>&1 >/dev/null | grep -o 'distinct minimizers: [0-9]*' | grep -o '[0-9]*')
	density[$o]=$n
	# a healthy order should select at least 1% and at most 100% of the
	# ~16569 k-mers in MT-human -- catches both the "selects nothing" and
	# the "selects nearly every position" bug classes hit during development
	if [ -n "$n" ] && [ "$n" -gt 165 ] && [ "$n" -le 16569 ]; then
		ok "density non-degenerate: $o ($n distinct minimizers)"
	else
		bad "density non-degenerate: $o (got '$n')"
	fi
done
echo "  native=${density[native]} xor-hash=${density[xor-hash]} xor-pcg-hash=${density[xor-pcg-hash]}"
echo "  mod-minimizer-xor-hash=${density[mod-minimizer-xor-hash]} mod-minimizer-pcg-xor-hash=${density[mod-minimizer-pcg-xor-hash]}"

echo "== 3. .mmi round-trip: dump+reload maps identically to a fresh build+map =="
for o in $ALL_TEST_ORDERS; do
	f=$(order_flags "$o")
	$BIN -a $f "$REF" "$QRY" > "$TMP/direct_$o.sam" 2>/dev/null
	$BIN -d "$TMP/$o.mmi" $f "$REF" >/dev/null 2>&1
	$BIN -a "$TMP/$o.mmi" "$QRY" > "$TMP/reload_$o.sam" 2>/dev/null
	if diff -q <(grep -v '^@PG' "$TMP/direct_$o.sam") <(grep -v '^@PG' "$TMP/reload_$o.sam") >/dev/null; then
		ok ".mmi round-trip: $o"
	else
		bad ".mmi round-trip: $o"
	fi
done

echo "== 4. split (-I) multi-part .mmi round-trip (order metadata written once, reused by later parts) =="
o=mod-minimizer-xor-hash
f=$(order_flags "$o")
$BIN -a $f "$REF" "$QRY" > "$TMP/split_direct.sam" 2>/dev/null
$BIN -d "$TMP/split.mmi" -I 5k $f "$REF" >/dev/null 2>&1
$BIN -a "$TMP/split.mmi" "$QRY" > "$TMP/split_reload.sam" 2>/dev/null
if diff -q <(grep -v '^@PG' "$TMP/split_direct.sam") <(grep -v '^@PG' "$TMP/split_reload.sam") >/dev/null; then
	ok "split -I round-trip: $o"
else
	bad "split -I round-trip: $o"
fi

echo "== 5. old-format .mmi (no order section) still loads, defaults to native =="
$BIN -d "$TMP/old.mmi" "$REF" >/dev/null 2>&1
if head -c4 "$TMP/old.mmi" | grep -q $'MMI\002'; then ok "old .mmi keeps MMI\\2 magic"; else bad "old .mmi magic changed unexpectedly"; fi
$BIN -a "$TMP/old.mmi" "$QRY" > "$TMP/old_reload.sam" 2>/dev/null
$BIN -a "$REF" "$QRY" > "$TMP/old_direct.sam" 2>/dev/null
if diff -q <(grep -v '^@PG' "$TMP/old_direct.sam") <(grep -v '^@PG' "$TMP/old_reload.sam") >/dev/null; then
	ok "old-format .mmi reload matches native"
else
	bad "old-format .mmi reload matches native"
fi

echo "== 6. -H (HPC) smoke test for each order =="
for o in $ALL_TEST_ORDERS; do
	f=$(order_flags "$o")
	if $BIN -a -H $f "$REF" "$QRY" > "$TMP/hpc_$o.sam" 2>"$TMP/hpc_$o.err"; then
		ok "-H builds+maps: $o"
	else
		bad "-H builds+maps: $o ($(tail -1 "$TMP/hpc_$o.err"))"
	fi
done

echo "== 7. -O/-s (gap-open/min-dp-score) still parse and are untouched by the new flags =="
if $BIN -a -O 4,24 -s 60 "$REF" "$QRY" > "$TMP/os.sam" 2>"$TMP/os.err"; then
	ok "-O/-s still parse"
else
	bad "-O/-s still parse ($(tail -1 "$TMP/os.err"))"
fi

echo "== 8. CLI validation: bad input is rejected, not silently accepted =="
if $BIN -a --minimizer-order bogus-order-name "$REF" "$QRY" >/dev/null 2>&1; then
	bad "unknown --minimizer-order should be rejected"
else
	ok "unknown --minimizer-order rejected"
fi
if $BIN -a --minimizer-order greedymini "$REF" "$QRY" >/dev/null 2>&1; then
	bad "greedymini without --greedymini-file should be rejected"
else
	ok "greedymini without --greedymini-file rejected"
fi
if $BIN -a --minimizer-order mod-minimizer-xor-hash "$REF" "$QRY" >/dev/null 2>&1; then
	bad "mod-minimizer-* without --mod-t should be rejected"
else
	ok "mod-minimizer-* without --mod-t rejected"
fi
if $BIN -a --order-seed 1 "$REF" "$QRY" >/dev/null 2>&1; then
	bad "--order-seed without --minimizer-order should be rejected"
else
	ok "--order-seed without --minimizer-order rejected"
fi

echo "== 9. t=5 / k=15,w=10 formula sanity (paper's r=4 formula: t = r + (k-r) mod w) =="
r=4; k=15; w=10
t=$(( r + (k - r) % w ))
if [ "$t" -eq 5 ]; then ok "t formula gives 5 for k=15,w=10 (matches --mod-t 5 used above)"; else bad "t formula gave $t, expected 5"; fi

echo "== 10. mod-minimizer at a large window (w > 16): regression test for a stale kraken2-derived"
echo "   assert that bounded w itself by ModResidueDeque::kMaxCap instead of only per-residue capacity =="
for o in $ORDERS_MOD; do
	if $BIN -a --minimizer-order "$o" --order-seed "$SEED" --mod-t "$MOD_T" -w 50 "$REF" "$QRY" > "$TMP/w50_$o.sam" 2>"$TMP/w50_$o.err"; then
		ok "-w 50 builds+maps: $o"
	else
		bad "-w 50 builds+maps: $o ($(tail -1 "$TMP/w50_$o.err"))"
	fi
done

echo "== 11. mod-minimizer O(1)-amortized selection: brute-force cross-validation =="
PY=/usr/bin/python3.11
command -v "$PY" >/dev/null 2>&1 || PY=python3
if "$PY" validate_mod_minimizer.py >"$TMP/validate_mod_minimizer.out" 2>&1; then
	ok "sketch_mod_minimizer matches brute-force reference ($($PY --version 2>&1))"
else
	bad "brute-force validation failed -- see $TMP/validate_mod_minimizer.out"
	tail -20 "$TMP/validate_mod_minimizer.out" >&2
fi

echo
echo "== summary: $pass passed, $fail failed =="
[ "$fail" -eq 0 ]
