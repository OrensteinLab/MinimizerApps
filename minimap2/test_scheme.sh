#!/usr/bin/env bash
# Sanity/regression tests for the custom minimizer order integration
# (scheme.h/scheme.cpp). See the plan's "Sanity checks / tests" section.
#
# greedymini uses orders/GreedyMini/w10_k15_canonical.gm (gm_k=15) at this
# script's -k 15 (no extension) and, centre-extended, at -k 19 (section 12),
# with both AT/CG projections (the `greedymini_at1-cg0` test label = greedymini with --gm-projection at1-cg0). The
# mod-minimizer schemes are strand-symmetric: mod-minimizer-xor (O_t = xor, sections 1-13)
# and mod-minimizer-greedymini (O_t = the w20_k9 table at t = 9, -k 19; section 15).
#
# Usage: ./test_scheme.sh (run from minimap2/, after `make all sketch_dump`)
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

SCHEMES_SEEDED="xor xor-pcg-hash"
SCHEMES_MOD="mod-minimizer-xor"
ALL_TEST_SCHEMES="native $SCHEMES_SEEDED greedymini greedymini_at1-cg0 $SCHEMES_MOD"
SEED=42
MOD_T=5
GM_FILE="$(cd .. && pwd)/orders/GreedyMini/w10_k15_canonical.gm"

scheme_flags() {
	case "$1" in
		native) echo "" ;;
		greedymini) echo "--scheme greedymini --greedymini-file $GM_FILE" ;;
		greedymini_at1-cg0) echo "--scheme greedymini --greedymini-file $GM_FILE --gm-projection at1-cg0" ;;
		mod-minimizer-*) echo "--scheme $1 --order-seed $SEED --mod-t $MOD_T" ;;
		*) echo "--scheme $1 --order-seed $SEED" ;;
	esac
}

echo "== 1. determinism: same order+seed, repeated runs -> identical output =="
for o in $ALL_TEST_SCHEMES; do
	f=$(scheme_flags "$o")
	$BIN -a $f "$REF" "$QRY" > "$TMP/r1_$o.sam" 2>/dev/null
	$BIN -a $f "$REF" "$QRY" > "$TMP/r2_$o.sam" 2>/dev/null
	if diff -q "$TMP/r1_$o.sam" "$TMP/r2_$o.sam" >/dev/null; then ok "determinism: $o"; else bad "determinism: $o"; fi
done

echo "== 2. density: each order shows a plausible, non-degenerate minimizer count =="
declare -A density
for o in $ALL_TEST_SCHEMES; do
	f=$(scheme_flags "$o")
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
echo "  native=${density[native]} xor=${density[xor]} xor-pcg-hash=${density[xor-pcg-hash]}"
echo "  greedymini=${density[greedymini]} greedymini_at1-cg0=${density[greedymini_at1-cg0]} mod-minimizer-xor=${density[mod-minimizer-xor]}"
# The flipped projection is a different order: same table, complemented keys.
if [ "${density[greedymini]}" != "${density[greedymini_at1-cg0]}" ]; then
	ok "--gm-projection at1-cg0 changes the selection"
else
	bad "--gm-projection at1-cg0 selected exactly the same minimizer count as at0-cg1"
fi

echo "== 3. .mmi round-trip: dump+reload maps identically to a fresh build+map =="
for o in $ALL_TEST_SCHEMES; do
	f=$(scheme_flags "$o")
	$BIN -a $f "$REF" "$QRY" > "$TMP/direct_$o.sam" 2>/dev/null
	$BIN -d "$TMP/$o.mmi" $f "$REF" >/dev/null 2>&1
	$BIN -a "$TMP/$o.mmi" "$QRY" > "$TMP/reload_$o.sam" 2>/dev/null
	if diff -q <(grep -v '^@PG' "$TMP/direct_$o.sam") <(grep -v '^@PG' "$TMP/reload_$o.sam") >/dev/null; then
		ok ".mmi round-trip: $o"
	else
		bad ".mmi round-trip: $o"
	fi
done

echo "== 4. split (-I) multi-part .mmi round-trip (scheme metadata written once, reused by later parts) =="
o=mod-minimizer-xor
f=$(scheme_flags "$o")
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
for o in $ALL_TEST_SCHEMES; do
	f=$(scheme_flags "$o")
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
if $BIN -a --scheme bogus-scheme-name "$REF" "$QRY" >/dev/null 2>&1; then
	bad "unknown --scheme should be rejected"
else
	ok "unknown --scheme rejected"
fi
if $BIN -a --scheme greedymini "$REF" "$QRY" >/dev/null 2>&1; then
	bad "greedymini without --greedymini-file should be rejected"
else
	ok "greedymini without --greedymini-file rejected"
fi
if $BIN -a --scheme mod-minimizer-xor "$REF" "$QRY" >/dev/null 2>&1; then
	bad "mod-minimizer-* without --mod-t should be rejected"
else
	ok "mod-minimizer-* without --mod-t rejected"
fi
if $BIN -a --scheme xor --gm-projection at1-cg0 "$REF" "$QRY" >/dev/null 2>&1; then
	bad "--gm-projection with a non-greedymini order should be rejected"
else
	ok "--gm-projection with a non-greedymini order rejected"
fi
if $BIN -a $(scheme_flags greedymini) --gm-projection bogus "$REF" "$QRY" >/dev/null 2>&1; then
	bad "unknown --gm-projection value should be rejected"
else
	ok "unknown --gm-projection value rejected"
fi
for old in xor-hash mod-minimizer-xor-hash mod-minimizer-pcg-xor-hash mod-minimizer-xor-canonical; do
	if $BIN -a --scheme $old --order-seed 1 --mod-t 5 "$REF" "$QRY" >/dev/null 2>&1; then
		bad "removed/renamed order name '$old' should be rejected"
	else
		ok "removed/renamed order name '$old' rejected"
	fi
done
if $BIN -a --order-seed 1 "$REF" "$QRY" >/dev/null 2>&1; then
	bad "--order-seed without --scheme should be rejected"
else
	ok "--order-seed without --scheme rejected"
fi

echo "== 9. t=5 / k=15,w=10 formula sanity (paper's r=4 formula: t = r + (k-r) mod w) =="
r=4; k=15; w=10
t=$(( r + (k - r) % w ))
if [ "$t" -eq 5 ]; then ok "t formula gives 5 for k=15,w=10 (matches --mod-t 5 used above)"; else bad "t formula gave $t, expected 5"; fi

echo "== 10. mod-minimizer window limits: t == k (mod w) (the centre rule), widest window, ring sized from w+k-t =="
for o in $SCHEMES_MOD; do
	if $BIN -a --scheme "$o" --order-seed "$SEED" --mod-t 1 -k 28 -w 27 "$REF" "$QRY" > "$TMP/w27_$o.sam" 2>"$TMP/w27_$o.err"; then
		ok "-k 28 -w 27 --mod-t 1 (L = 54) builds+maps: $o"
	else
		bad "-k 28 -w 27 --mod-t 1 builds+maps: $o ($(tail -1 "$TMP/w27_$o.err"))"
	fi
	for bad_tw in "6 10" "5 50"; do
		set -- $bad_tw
		if $BIN -a --scheme "$o" --order-seed "$SEED" --mod-t "$1" -w "$2" "$REF" "$QRY" >/dev/null 2>&1; then
			bad "--mod-t $1 -w $2 (t != k mod w) should be rejected: $o"
		else
			ok "--mod-t $1 -w $2 (t != k mod w) rejected: $o"
		fi
	done
done

echo "== 11. mod-minimizer selection: brute-force cross-validation (Python port + compiled sketch_dump) =="
PY=/usr/bin/python3.11
command -v "$PY" >/dev/null 2>&1 || PY=python3
if "$PY" ../benchmark/minimap2/tests/test_mod_minimizer_sketch.py >"$TMP/validate_mod_minimizer.out" 2>&1; then
	ok "sketch_mod_minimizer matches brute-force reference ($($PY --version 2>&1))"
else
	bad "brute-force validation failed -- see $TMP/validate_mod_minimizer.out"
	tail -20 "$TMP/validate_mod_minimizer.out" >&2
fi

echo "== 12. greedymini centre extension: k = gm_k + 2i (k=19 here), both projections =="
for proj in at0-cg1 at1-cg0; do
	f="--scheme greedymini --greedymini-file $GM_FILE --gm-projection $proj -k 19"
	$BIN -a $f "$REF" "$QRY" > "$TMP/k19_$proj.sam" 2>/dev/null
	$BIN -a $f "$REF" "$QRY" > "$TMP/k19b_$proj.sam" 2>/dev/null
	if [ -s "$TMP/k19_$proj.sam" ] && diff -q "$TMP/k19_$proj.sam" "$TMP/k19b_$proj.sam" >/dev/null; then ok "k=19 determinism: $proj"; else bad "k=19 determinism: $proj"; fi
	$BIN -d "$TMP/k19_$proj.mmi" $f "$REF" >/dev/null 2>&1
	$BIN -a "$TMP/k19_$proj.mmi" "$QRY" > "$TMP/k19r_$proj.sam" 2>/dev/null
	if diff -q <(grep -v '^@PG' "$TMP/k19_$proj.sam") <(grep -v '^@PG' "$TMP/k19r_$proj.sam") >/dev/null; then ok "k=19 .mmi round-trip: $proj"; else bad "k=19 .mmi round-trip: $proj"; fi
	if $BIN -a -H $f "$REF" "$QRY" >/dev/null 2>&1; then ok "k=19 -H builds+maps: $proj"; else bad "k=19 -H builds+maps: $proj"; fi
done
for badk in 16 25; do
	if $BIN -a --scheme greedymini --greedymini-file "$GM_FILE" -k $badk "$REF" "$QRY" >/dev/null 2>&1; then
		bad "greedymini -k $badk should be rejected (k - gm_k odd, or k > 23)"
	else
		ok "greedymini -k $badk rejected"
	fi
done
if "$PY" ../benchmark/minimap2/tests/test_greedymini_sketch.py >"$TMP/validate_greedymini.out" 2>&1; then
	ok "greedymini rank and sketch match the Python references (k=15..23, both projections)"
else
	bad "greedymini validation failed -- see $TMP/validate_greedymini.out"
	tail -20 "$TMP/validate_greedymini.out" >&2
fi

echo "== 13. robustness: bad k, unreadable stored orders, overridden options, rank sentinel =="
for kk in 29 32; do
	if $BIN -a $(scheme_flags mod-minimizer-xor) -k $kk "$REF" "$QRY" >/dev/null 2>&1; then
		bad "mod-minimizer-xor -k $kk should be rejected (k > 28)"
	else
		ok "mod-minimizer-xor -k $kk rejected"
	fi
done
# The order blob ends the file: [u64 len][u32 kind]...; patch the kind byte.
patch_kind() { # <in.mmi> <out.mmi> <kind>
	"$PY" - "$1" "$2" "$3" <<'EOF'
import struct, sys
b = bytearray(open(sys.argv[1], 'rb').read())
for off in range(len(b) - 12, 0, -1):  # find the len field that spans to EOF
	if struct.unpack_from('<Q', b, off)[0] == len(b) - off - 8:
		b[off + 8] = int(sys.argv[3]); break
open(sys.argv[2], 'wb').write(b)
EOF
}
for kind in 4 0 3; do # retired kind, native, greedymini without a table
	patch_kind "$TMP/xor.mmi" "$TMP/kind$kind.mmi" $kind
	if $BIN -a "$TMP/kind$kind.mmi" "$QRY" >/dev/null 2>"$TMP/kind$kind.err"; then
		bad "index with a patched order kind $kind should be refused"
	elif grep -q 'failed to load the minimizer order' "$TMP/kind$kind.err"; then
		ok "index with a patched order kind $kind refused"
	else
		bad "patched order kind $kind: unexpected error ($(tail -1 "$TMP/kind$kind.err"))"
	fi
done
# A mod-minimizer index written before the centre tie rule: its order section
# lacks the trailing tie-rule byte (drop it and shorten the section length).
"$PY" - "$TMP/mod-minimizer-xor.mmi" "$TMP/mod_gt.mmi" <<'EOF2'
import struct, sys
b = bytearray(open(sys.argv[1], 'rb').read())
for off in range(len(b) - 12, 0, -1):
	if struct.unpack_from('<Q', b, off)[0] == len(b) - off - 8:
		struct.pack_into('<Q', b, off, len(b) - off - 9); break
open(sys.argv[2], 'wb').write(b[:-1])
EOF2
if $BIN -a "$TMP/mod_gt.mmi" "$QRY" >/dev/null 2>"$TMP/mod_gt.err"; then
	bad "a mod-minimizer index with the retired G/T tie rule should be refused"
elif grep -q 'retired G/T-strand tie rule' "$TMP/mod_gt.err"; then
	ok "mod-minimizer index with the retired G/T tie rule refused"
else
	bad "retired tie rule: unexpected error ($(tail -1 "$TMP/mod_gt.err"))"
fi
head -c $(( $(stat -c %s "$TMP/xor.mmi") - 4 )) "$TMP/xor.mmi" > "$TMP/trunc.mmi"
if $BIN -a "$TMP/trunc.mmi" "$QRY" >/dev/null 2>&1; then bad "truncated order section should be refused"; else ok "truncated order section refused"; fi
$BIN -a --scheme xor-pcg-hash --order-seed $SEED "$TMP/xor.mmi" "$QRY" >/dev/null 2>"$TMP/override.err"
if grep -q 'overridden by the scheme of the prebuilt index (xor)' "$TMP/override.err"; then ok "scheme options overridden by a prebuilt index: warned"; else bad "no warning when --scheme differs from the prebuilt index"; fi
$BIN -a $(scheme_flags xor) "$TMP/xor.mmi" "$QRY" >/dev/null 2>"$TMP/same.err"
if grep -q 'Minimizer scheme options' "$TMP/same.err"; then bad "warning although the options match the prebuilt index"; else ok "no warning when the options match the prebuilt index"; fi
$BIN -a $(scheme_flags greedymini) "$TMP/greedymini.mmi" "$QRY" >/dev/null 2>"$TMP/same_gm.err"
if grep -q 'Minimizer scheme options' "$TMP/same_gm.err"; then bad "warning although the greedymini table matches the prebuilt index"; else ok "no warning when the greedymini table matches the prebuilt index"; fi
polya=$(printf '>a\n%0100d\n' 0 | tr 0 A)
for o in xor xor-pcg-hash; do
	n=$(echo "$polya" | ./sketch_dump $o 10 15 0xFFFFFFFFFFFFFFFF 0 | grep -c .)
	if [ "$n" -gt 0 ]; then ok "$o seed ~0: poly-A still selected ($n)"; else bad "$o seed ~0: poly-A never selected"; fi
done

echo "== 14. greedymini --order-seed: the Gamma tie-break toggle =="
fs="--scheme greedymini --greedymini-file $GM_FILE --order-seed $SEED -k 19"
$BIN -a $fs "$REF" "$QRY" > "$TMP/gms.sam" 2>/dev/null
$BIN -d "$TMP/gms.mmi" $fs "$REF" >/dev/null 2>&1
$BIN -a "$TMP/gms.mmi" "$QRY" > "$TMP/gms_r.sam" 2>/dev/null
if [ -s "$TMP/gms.sam" ] && diff -q <(grep -v '^@PG' "$TMP/gms.sam") <(grep -v '^@PG' "$TMP/gms_r.sam") >/dev/null; then ok "seeded greedymini .mmi round-trip"; else bad "seeded greedymini .mmi round-trip"; fi
grep -v '>' "$REF" | tr -d '\n' > "$TMP/ref.seq"
./sketch_dump greedymini 10 15 0 0 "$GM_FILE" < "$TMP/ref.seq" > "$TMP/gs0.txt" && ./sketch_dump greedymini 10 15 0x9E3779B97F4A7C15 0 "$GM_FILE" < "$TMP/ref.seq" > "$TMP/gs1.txt"; rc=$?
n0=$(md5sum < "$TMP/gs0.txt"); n1=$(md5sum < "$TMP/gs1.txt")
if [ "$rc" = 0 ] && [ -s "$TMP/gs0.txt" ] && [ -s "$TMP/gs1.txt" ] && [ "$n0" != "$n1" ]; then ok "--order-seed changes the greedymini selection"; else bad "--order-seed did not change the greedymini selection"; fi
$BIN -a --scheme greedymini --greedymini-file "$GM_FILE" --order-seed 7 "$TMP/gms.mmi" "$QRY" >/dev/null 2>"$TMP/gms_w.err"
if grep -q 'Minimizer scheme options' "$TMP/gms_w.err"; then ok "a different greedymini seed than the prebuilt index's: warned"; else bad "no warning for a different greedymini seed"; fi

echo "== 15. mod-minimizer-greedymini: O_t = w20_k9 at t = 9 (-k 19 -w 10), both projections =="
GM_MOD_FILE="$(cd .. && pwd)/orders/GreedyMini/w20_k9_canonical.gm"
declare -A mgd
for proj in at0-cg1 at1-cg0; do
	f="--scheme mod-minimizer-greedymini --greedymini-file $GM_MOD_FILE --gm-projection $proj --order-seed $SEED --mod-t 9 -k 19"
	$BIN -a $f "$REF" "$QRY" > "$TMP/mg_$proj.sam" 2>"$TMP/mg_$proj.err"
	$BIN -a $f "$REF" "$QRY" > "$TMP/mgb_$proj.sam" 2>/dev/null
	if [ -s "$TMP/mg_$proj.sam" ] && diff -q "$TMP/mg_$proj.sam" "$TMP/mgb_$proj.sam" >/dev/null; then ok "determinism: $proj"; else bad "determinism: $proj"; fi
	mgd[$proj]=$(grep -o 'distinct minimizers: [0-9]*' "$TMP/mg_$proj.err" | grep -o '[0-9]*')
	$BIN -d "$TMP/mg_$proj.mmi" $f "$REF" >/dev/null 2>&1
	$BIN -a "$TMP/mg_$proj.mmi" "$QRY" > "$TMP/mgr_$proj.sam" 2>/dev/null
	if diff -q <(grep -v '^@PG' "$TMP/mg_$proj.sam") <(grep -v '^@PG' "$TMP/mgr_$proj.sam") >/dev/null; then ok ".mmi round-trip: $proj"; else bad ".mmi round-trip: $proj"; fi
	$BIN -a $f "$TMP/mg_$proj.mmi" "$QRY" >/dev/null 2>"$TMP/mgw_$proj.err"
	if grep -q 'Minimizer scheme options' "$TMP/mgw_$proj.err"; then bad "warning although the options match the prebuilt index: $proj"; else ok "no warning when the options match the prebuilt index: $proj"; fi
	if $BIN -a -H $f "$REF" "$QRY" >/dev/null 2>&1; then ok "-H builds+maps: $proj"; else bad "-H builds+maps: $proj"; fi
done
if [ -n "${mgd[at0-cg1]}" ] && [ "${mgd[at0-cg1]}" != "${mgd[at1-cg0]}" ]; then ok "--gm-projection at1-cg0 changes the selection (${mgd[at0-cg1]} vs ${mgd[at1-cg0]})"; else bad "--gm-projection did not change the mod-minimizer-greedymini selection"; fi
grep -v '>' "$REF" | tr -d '\n' > "$TMP/ref.seq"
./sketch_dump mod-minimizer-greedymini 10 19 0 9 "$GM_MOD_FILE" < "$TMP/ref.seq" > "$TMP/ms0.txt" && ./sketch_dump mod-minimizer-greedymini 10 19 0x9E3779B97F4A7C15 9 "$GM_MOD_FILE" < "$TMP/ref.seq" > "$TMP/ms1.txt"; rc=$?
n0=$(md5sum < "$TMP/ms0.txt"); n1=$(md5sum < "$TMP/ms1.txt")
if [ "$rc" = 0 ] && [ -s "$TMP/ms0.txt" ] && [ -s "$TMP/ms1.txt" ] && [ "$n0" != "$n1" ]; then ok "--order-seed changes the selection"; else bad "--order-seed did not change the selection"; fi
if $BIN -a --scheme mod-minimizer-greedymini --mod-t 9 -k 19 "$REF" "$QRY" >/dev/null 2>&1; then bad "without --greedymini-file should be rejected"; else ok "without --greedymini-file rejected"; fi
if $BIN -a --scheme mod-minimizer-greedymini --greedymini-file "$GM_MOD_FILE" --mod-t 5 -k 15 "$REF" "$QRY" >/dev/null 2>&1; then bad "t = 5 with a gm_k = 9 table should be rejected"; else ok "t = 5 with a gm_k = 9 table rejected"; fi

echo
echo "== summary: $pass passed, $fail failed =="
[ "$fail" -eq 0 ]
