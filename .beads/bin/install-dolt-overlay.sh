#!/bin/bash
# Idempotent installer that relocates dolt's storage off virtiofs to overlayfs
# so directory fsync (required by dolt's manifest-rename protocol) actually works.
#
# Background: this Docker sandbox bind-mounts /Users/bence/code/QuantumFurnace.jl
# as virtiofs. virtiofs returns EBADF on directory fsync after a manifest file
# rename, which corrupts dolt's noms store on every write. Symptom:
#   "fatal error: sync .../noms: bad file descriptor:
#    error fsyncing directory after manifest file rename"
# /home/agent is on overlayfs, where fsync works. Overlay storage is ephemeral
# (wipes on container rebuild), so we lazily rebuild from .beads/issues.jsonl
# whenever the overlay store is missing or empty.
#
# What this script does (idempotently):
#   1. Ensures the overlay parent dir exists with 0700 perms.
#   2. Drops a stub metadata.json there so bd's project-root discovery (which
#      follows symlinks) finds the right database name "QuantumFurnace".
#   3. Replaces .beads/embeddeddolt and .beads/dolt with symlinks to overlay
#      paths if they're real dirs (preserves any prior state as .corrupt.backup
#      so nothing is silently destroyed).
#   4. If the overlay db is missing, runs bd init --from-jsonl --force to
#      rebuild from the source-of-truth issues.jsonl.
#
# Safe to run multiple times.
set -e

REPO_ROOT="$(cd "$(dirname "$(readlink -f "$0")")"/../.. && pwd)"
BEADS_DIR="$REPO_ROOT/.beads"
OVERLAY_BASE="/home/agent/qf-beads-dolt"

[ -d "$BEADS_DIR" ] || exit 0  # not a beads project; nothing to do
[ -f "$BEADS_DIR/issues.jsonl" ] || exit 0  # no source of truth; bail

mkdir -p "$OVERLAY_BASE/embeddeddolt" "$OVERLAY_BASE/dolt"
chmod 700 "$OVERLAY_BASE"

# Stub metadata.json at the symlink-target parent so bd's project-discovery
# (which resolves symlinks before walking up) picks up the right db name.
META="$OVERLAY_BASE/metadata.json"
SRC_META="$BEADS_DIR/metadata.json"
if [ -f "$SRC_META" ]; then
  if ! cmp -s "$SRC_META" "$META" 2>/dev/null; then
    cp "$SRC_META" "$META"
  fi
fi

# Replace real dirs with symlinks (preserve any existing state as .corrupt.backup).
relink() {
  local name="$1"
  local target="$OVERLAY_BASE/$name"
  local link="$BEADS_DIR/$name"
  if [ -L "$link" ]; then
    [ "$(readlink "$link")" = "$target" ] && return 0
    rm "$link"
  elif [ -e "$link" ]; then
    local ts
    ts="$(date +%Y%m%d-%H%M%S)"
    mv "$link" "$link.corrupt.backup.$ts"
  fi
  ln -s "$target" "$link"
}
relink embeddeddolt
relink dolt

# Lazily rebuild the dolt store from JSONL if the overlay is empty.
# (Container rebuilds wipe overlay storage; this restores it.)
DB_DIR="$OVERLAY_BASE/embeddeddolt/QuantumFurnace"
if [ ! -d "$DB_DIR/.dolt" ]; then
  if command -v bd >/dev/null 2>&1; then
    (cd "$REPO_ROOT" && bd init --prefix qf --from-jsonl --force >/dev/null 2>&1) || true
  fi
fi
