# Dataset CSV helpers + download (Zenodo eulertigs + ENA query FASTQ) + md5 verification.
#
# data/datasets.csv columns:
#   dataset,species,eulertigs_filename,eulertigs_size,eulertigs_md5,eulertigs_url,
#   query_filename,query_size,query_md5,query_url,local_path,local_query_path
#
# `local_path` (and optional `local_query_path`) let a dataset row point at an
# already-on-disk file instead of downloading anything -- used by the `smoke`
# row (the in-repo sshash sample dataset) so the whole pipeline (download ->
# run -> export) can be smoke-tested without any network access. A
# `--dataset-file PATH` CLI override does the same for an arbitrary ad hoc
# local FASTA, without needing a datasets.csv row at all.

csv_field() {
  local ds="$1" col="$2"
  awk -F, -v ds="$ds" -v col="$col" '
    NR == 1 {
      for (i = 1; i <= NF; i++) h[$i] = i
      next
    }
    $1 == ds && h[col] { print $h[col]; exit }
  ' "$DATASETS_CSV"
}

dataset_is_local() {
  local ds="$1" lp
  lp="$(csv_field "$ds" local_path)"
  [[ -n "$(trim "$lp")" ]]
}

# Absolute path to dataset `ds`'s eulertigs FASTA (downloaded, or local_path).
dataset_fasta_path() {
  local ds="$1" lp fn
  lp="$(trim "$(csv_field "$ds" local_path)")"
  if [[ -n "$lp" ]]; then
    case "$lp" in
      /*) echo "$lp" ;;
      *) echo "$(cd "$BENCHMARK_DIR" && cd "$(dirname "$lp")" && pwd)/$(basename "$lp")" ;;
    esac
    return
  fi
  fn="$(trim "$(csv_field "$ds" eulertigs_filename)")"
  [[ -n "$fn" ]] || die "dataset $ds: no eulertigs_filename or local_path in $DATASETS_CSV"
  echo "$DATA_DIR/data/$ds/$fn"
}

# Absolute path to dataset `ds`'s query FASTQ, or empty if none configured
# (e.g. the smoke dataset, which has no paired real-read query file).
dataset_query_path() {
  local ds="$1" lqp fn
  lqp="$(trim "$(csv_field "$ds" local_query_path)")"
  if [[ -n "$lqp" ]]; then
    case "$lqp" in
      /*) echo "$lqp" ;;
      *) echo "$(cd "$BENCHMARK_DIR" && cd "$(dirname "$lqp")" && pwd)/$(basename "$lqp")" ;;
    esac
    return
  fi
  fn="$(trim "$(csv_field "$ds" query_filename)")"
  [[ -n "$fn" ]] || return 0
  echo "$DATA_DIR/data/$ds/$fn"
}

# ENA fastq path convention: subfolder is empty for a 6-digit run number,
# else 00<last digit> / 0<last 2 digits> / <last 3 digits> for 7/8/9 digits.
ena_url() {
  local filename="$1" acc digits prefix6 len sub=""
  if [[ "$filename" =~ ^((SRR|ERR|DRR)[0-9]+) ]]; then
    acc="${BASH_REMATCH[1]}"
  else
    die "cannot parse accession from $filename"
  fi
  digits="${acc:3}"
  prefix6="${acc:0:6}"
  len=${#digits}
  if ((len == 7)); then
    sub="/00${digits: -1}"
  elif ((len == 8)); then
    sub="/0${digits: -2}"
  elif ((len == 9)); then
    sub="/${digits: -3}"
  fi
  echo "${ENA_FASTQ}/${prefix6}${sub}/${acc}/${filename}"
}

file_complete() {
  [[ -f "$1" && -s "$1" ]]
}

md5_of() {
  if command -v md5sum >/dev/null 2>&1; then
    md5sum "$1" | awk '{print $1}'
  else
    md5 -q "$1"
  fi
}

verify_md5() {
  local path="$1" expected="$2" actual
  [[ -n "$expected" ]] || return 0
  actual="$(md5_of "$path")"
  if [[ "$actual" != "$expected" ]]; then
    die "md5 mismatch for $path: expected $expected, got $actual"
  fi
  echo "md5 OK: $(basename "$path")"
}

curl_code_in_list() {
  local code="$1" c
  for c in "${CURL_RESUME_CODES[@]}"; do
    [[ "$code" -eq "$c" ]] && return 0
  done
  return 1
}

download_one() {
  local url="$1" dest="$2"
  if file_complete "$dest"; then
    echo "skip existing $dest"
    return 0
  fi
  command -v curl >/dev/null 2>&1 || die "curl is required for download"
  mkdir -p "$(dirname "$dest")"
  local tmp="${dest}.part"
  local attempts=40
  local i have wait last=1
  local connect_timeout=$(( DOWNLOAD_TIMEOUT < 60 ? DOWNLOAD_TIMEOUT : 60 ))
  for ((i = 1; i <= attempts; i++)); do
    have=0
    [[ -f "$tmp" ]] && have=$(stat -c%s "$tmp" 2>/dev/null || echo 0)
    echo "curl -fL -C - --connect-timeout $connect_timeout -o $tmp $url  # ${have} bytes already"
    if curl -fL -C - --connect-timeout "$connect_timeout" -o "$tmp" "$url"; then
      mv -f "$tmp" "$dest"
      return 0
    fi
    last=$?
    if ! curl_code_in_list "$last"; then
      die "curl failed for $url (exit $last)"
    fi
    have=0
    [[ -f "$tmp" ]] && have=$(stat -c%s "$tmp" 2>/dev/null || echo 0)
    wait=$(( 5 * i ))
    (( wait > 60 )) && wait=60
    echo "warning: curl exit $last after ${have} bytes of $(basename "$dest"); resume in ${wait}s ($i/$attempts)" >&2
    sleep "$wait"
  done
  die "curl failed for $url after $attempts attempts (exit $last)"
}

# A read file shared by several datasets (hprc is queried with human's SRR5833294.fastq.gz): when
# `dest` is missing but a non-empty file of the same name sits in another dataset's directory,
# symlink it instead of downloading it again.
link_existing_query() {
  local dest="$1" fn other
  fn="$(basename "$dest")"
  [[ -s "$dest" ]] && return 0
  for other in "$DATA_DIR"/data/*/"$fn"; do
    if [[ -s "$other" && "$(readlink -f "$other")" != "$(readlink -f "$dest")" ]]; then
      ln -sf "$(readlink -f "$other")" "$dest"
      echo "$fn: linked to $(readlink -f "$other")"
      return 0
    fi
  done
  return 1
}

download_dataset() {
  local ds="$1"
  if dataset_is_local "$ds"; then
    echo "$ds: local_path dataset, nothing to download"
    return 0
  fi
  local out_dir="$DATA_DIR/data/$ds"
  mkdir -p "$out_dir"

  local eu_fn eu_url eu_md5 eu_dest
  eu_fn="$(trim "$(csv_field "$ds" eulertigs_filename)")"
  eu_url="$(trim "$(csv_field "$ds" eulertigs_url)")"
  eu_md5="$(trim "$(csv_field "$ds" eulertigs_md5)")"
  [[ -n "$eu_fn" && -n "$eu_url" ]] || die "dataset $ds: missing eulertigs_filename/eulertigs_url in $DATASETS_CSV"
  eu_dest="$out_dir/$eu_fn"
  if download_one "$eu_url" "$eu_dest"; then
    verify_md5 "$eu_dest" "$eu_md5"
  elif [[ "$SKIP_MISSING" -ne 1 ]]; then
    die "download failed: $eu_url"
  fi

  local q_fn q_url q_md5 q_dest
  q_fn="$(trim "$(csv_field "$ds" query_filename)")"
  if [[ -n "$q_fn" ]]; then
    q_url="$(trim "$(csv_field "$ds" query_url)")"
    [[ -n "$q_url" ]] || q_url="$(ena_url "$q_fn")"
    q_md5="$(trim "$(csv_field "$ds" query_md5)")"
    q_dest="$out_dir/$q_fn"
    link_existing_query "$q_dest" || true
    if download_one "$q_url" "$q_dest"; then
      verify_md5 "$q_dest" "$q_md5"
    elif [[ "$SKIP_MISSING" -ne 1 ]]; then
      die "download failed: $q_url"
    fi
  fi
  echo "$ds: download complete -> $out_dir"
}

cmd_download() {
  local ds
  IFS=',' read -r -a ds_arr <<<"$DATASETS"
  for ds in "${ds_arr[@]}"; do
    ds="$(trim "$ds")"
    [[ -n "$ds" ]] || continue
    download_dataset "$ds"
  done
}
