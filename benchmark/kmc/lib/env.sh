# .env loading and small shell helpers

load_env_file() {
  local f="$1"
  [[ -f "$f" ]] || return 0
  while IFS= read -r line || [[ -n "$line" ]]; do
    line="${line%%#*}"
    line="$(echo "$line" | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')"
    [[ -z "$line" || "$line" != *=* ]] && continue
    local key="${line%%=*}"
    local val="${line#*=}"
    key="$(echo "$key" | sed 's/[[:space:]]*$//')"
    val="$(echo "$val" | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')"
    if [[ ${#val} -ge 2 && ( "$val" == \"*\" || "$val" == \'*\' ) ]]; then
      val="${val:1:${#val}-2}"
    fi
    if [[ -n "$key" && -z "${!key:-}" ]]; then
      export "$key=$val"
    fi
  done <"$f"
}

load_dotenv() {
  load_env_file "$HERE/.env"
  load_env_file "$REPO_ROOT/.env"
}

die() {
  echo "error: $*" >&2
  exit 1
}

trim() {
  echo "$1" | sed 's/^[[:space:]]*//;s/[[:space:]]*$//'
}
