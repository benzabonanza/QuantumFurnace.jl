#!/bin/bash
# Idempotent installer that ensures the bd binary at the npm location is
# wrapped with tini for zombie reaping. Runs from the SessionStart hook.
#
# Background: this Docker sandbox has `sleep infinity` as PID 1 with no
# subreaper. Without wrapping, every `bd` call leaves dolt subprocess zombies
# that never get reaped (PID 1 doesn't reap), eventually exhausting PIDs.
# Tini -s reaps most of them; the rest are cosmetic and clear on restart.
#
# This script:
#   1. Ensures the vendored tini binary exists and is executable.
#   2. Detects if the bd binary at the npm path is the wrapper or the original.
#   3. If original, swaps it: bd -> bd.real, installs wrapper as bd.
#
# Safe to run multiple times. Re-runs after container rebuild restore the wrap.
set -e

REPO_ROOT="$(cd "$(dirname "$(readlink -f "$0")")"/../.. && pwd)"
TINI="$REPO_ROOT/.beads/bin/tini"
BD_DIR="/usr/local/share/npm-global/lib/node_modules/@beads/bd/bin"
BD="$BD_DIR/bd"
BD_REAL="$BD_DIR/bd.real"

[ -x "$TINI" ] || { echo "install-bd-wrapper: tini missing at $TINI" >&2; exit 0; }
[ -d "$BD_DIR" ] || exit 0  # bd not installed in this environment

# If bd is already the wrapper (small shell script), nothing to do.
if [ -f "$BD" ] && head -1 "$BD" 2>/dev/null | grep -q "^#!/bin/bash"; then
  exit 0
fi

# bd is the original binary; wrap it.
[ -f "$BD_REAL" ] || sudo mv "$BD" "$BD_REAL"

sudo tee "$BD" >/dev/null <<'WRAPPER'
#!/bin/bash
# Wraps the real bd binary with tini as a subreaper (-s) so dolt subprocess
# zombies are reaped. PID 1 in this sandbox is `sleep infinity` (no reaper).
set -e
REAL_BD="$(dirname "$(readlink -f "$0")")/bd.real"
TINI=""
d="$PWD"
while [ "$d" != "/" ]; do
  if [ -x "$d/.beads/bin/tini" ]; then TINI="$d/.beads/bin/tini"; break; fi
  d="$(dirname "$d")"
done
[ -z "$TINI" ] && [ -x /usr/bin/tini ] && TINI=/usr/bin/tini
[ -z "$TINI" ] && [ -x /usr/local/bin/tini ] && TINI=/usr/local/bin/tini
if [ -n "$TINI" ]; then
  exec "$TINI" -s -- "$REAL_BD" "$@"
else
  exec "$REAL_BD" "$@"
fi
WRAPPER
sudo chmod +x "$BD"
