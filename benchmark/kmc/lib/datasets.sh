# Dataset CSV helpers, ENA download, and FASTQ prepare

csv_rows_for_dataset() {
  local ds="$1"
  awk -F, -v ds="$ds" '
    NR == 1 { next }
    $1 == ds { print }
  ' "$DATASETS_CSV"
}

csv_field() {
  local ds="$1" col="$2"
  awk -F, -v ds="$ds" -v col="$col" '
    NR == 1 {
      for (i = 1; i <= NF; i++) if ($i == col) c = i
      next
    }
    $1 == ds && c { print $c; exit }
  ' "$DATASETS_CSV"
}

dataset_mem() {
  local ds="$1"
  if [[ -n "$MEM_OVERRIDE" ]]; then
    echo "$MEM_OVERRIDE"
  else
    csv_field "$ds" mem_gb
  fi
}

ena_url() {
  local filename="$1"
  if [[ "$filename" =~ ^((SRR|ERR|DRR)[0-9]+) ]]; then
    local acc="${BASH_REMATCH[1]}"
    echo "${ENA_FASTQ}/${acc:0:6}/${acc}/${filename}"
  else
    die "cannot parse accession from $filename"
  fi
}

file_complete() {
  [[ -f "$1" && -s "$1" ]]
}

write_list() {
  local path="$1"
  shift
  mkdir -p "$(dirname "$path")"
  : >"$path"
  local f
  for f in "$@"; do
    echo "$(cd "$(dirname "$f")" && pwd)/$(basename "$f")" >>"$path"
  done
}

curl_code_in_list() {
  local code="$1"
  local c
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

download_dataset() {
  local ds="$1"
  local out_dir="$DATA_DIR/data/$ds/fastq_gz"
  mkdir -p "$out_dir"
  local -a kept=()
  local line dataset filename url dest ok
  while IFS= read -r line || [[ -n "$line" ]]; do
    [[ -z "$line" ]] && continue
    IFS=',' read -r dataset _ _ filename url <<<"$line"
    dataset="$(trim "$dataset")"
    filename="$(trim "$filename")"
    url="$(trim "$url")"
    url="${url//$'\r'/}"
    [[ -n "$filename" ]] || continue
    if [[ -z "$url" ]]; then
      url="$(ena_url "$filename")"
    fi
    dest="$out_dir/$filename"
    ok=0
    if download_one "$url" "$dest"; then
      ok=1
    elif [[ "$SKIP_MISSING" -eq 1 ]]; then
      echo "warning: failed $url" >&2
    else
      die "download failed: $url"
    fi
    if [[ "$ok" -eq 1 ]] && file_complete "$dest"; then
      kept+=("$dest")
    fi
  done < <(csv_rows_for_dataset "$ds")

  if ((${#kept[@]} == 0)); then
    if [[ "$SKIP_MISSING" -eq 1 ]]; then
      echo "warning: no files downloaded for $ds; skipping dataset" >&2
      return 0
    fi
    die "no files downloaded for $ds"
  fi
  write_list "$DATA_DIR/data/$ds/files_gz.lst" "${kept[@]}"
  echo "$ds: ${#kept[@]} gzip FASTQ(s) -> $DATA_DIR/data/$ds/files_gz.lst"
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

prepare_dataset() {
  local ds="$1"
  local gz_list="$DATA_DIR/data/$ds/files_gz.lst"
  [[ -f "$gz_list" ]] || die "missing $gz_list; run download first"
  local plain_root="$DATA_DIR/data/$ds/fastq_plain"
  mkdir -p "$plain_root"
  local -a plains=()
  local gz name plain_name dest
  while IFS= read -r gz || [[ -n "$gz" ]]; do
    gz="$(trim "$gz")"
    [[ -z "$gz" ]] && continue
    [[ -f "$gz" ]] || die "missing gzip file $gz"
    name="$(basename "$gz")"
    if [[ "$name" == *.fastq.gz ]]; then
      plain_name="${name%.gz}"
    elif [[ "$name" == *.fq.gz ]]; then
      plain_name="${name%.gz}"
    else
      plain_name="${name}.fastq"
    fi
    dest="$plain_root/$plain_name"
    if file_complete "$dest" && [[ "$FORCE" -eq 0 ]]; then
      echo "skip existing $dest"
    else
      echo "gunzip $gz -> $dest"
      gunzip -c "$gz" >"${dest}.part"
      mv -f "${dest}.part" "$dest"
    fi
    plains+=("$dest")
  done <"$gz_list"
  write_list "$DATA_DIR/data/$ds/files_plain.lst" "${plains[@]}"
  echo "$ds: ${#plains[@]} uncompressed FASTQ(s)"
}

cmd_prepare() {
  local ds
  IFS=',' read -r -a ds_arr <<<"$DATASETS"
  for ds in "${ds_arr[@]}"; do
    ds="$(trim "$ds")"
    [[ -n "$ds" ]] || continue
    prepare_dataset "$ds"
  done
}

inputs_include_plain() {
  local part
  IFS=',' read -r -a parts <<<"$INPUTS"
  for part in "${parts[@]}"; do
    [[ "$(trim "$part")" == "plain" ]] && return 0
  done
  return 1
}
