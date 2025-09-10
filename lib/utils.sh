require_file()  { [ -f "$1" ] || die "Required file not found: $1"; }
require_nonempty() { [ -n "$1" ] || die "Missing required value: $2"; }
ensure_parent_dir() {
  # ensure parent directory exists for output files
  local p; p="$(dirname -- "$1")"
  [ -d "$p" ] || mkdir -p "$p"
}

