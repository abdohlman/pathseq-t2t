#!/usr/bin/env bash

# Logging utility with timestamp + severity.
# Format: 00:38:03.980 INFO  message

# Default log level: 1 = INFO
: "${log_level:=1}"

_log_timestamp() {
  date +"%H:%M:%S.%3N"
}

_log() {
  local level="$1"
  local min_level="$2"
  shift 2
  if [[ ${log_level} -ge ${min_level} ]]; then
    printf "%s %-5s %s\n" "$(_log_timestamp)" "${level}" "$*"
  fi
}

log() { _log INFO 1 "$@"; }
warn() { _log WARN 0 "$@"; }
err() { _log ERROR 0 "$@"; }
die() { err "$@"; exit 1; }
