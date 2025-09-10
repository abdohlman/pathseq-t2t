log_level=${log_level:-1} # 0=quiet, 1=info, 2=debug

log()   { [ "${log_level}" -ge 1 ] && printf '[INFO] %s\n' "$*" >&2; }
debug() { [ "${log_level}" -ge 2 ] && printf '[DEBUG] %s\n' "$*" >&2; }
err()   { printf '[ERROR] %s\n' "$*" >&2; }
die()   { err "$*"; exit 2; }
