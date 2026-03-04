#!/usr/bin/env bash
# Ensure Claude Code's auto-memory directory points into the project repo,
# so memories survive Docker image rebuilds.
#
# Called automatically via .claude/settings.json PreToolUse hook.
# Safe to run repeatedly — all operations are idempotent.

set -euo pipefail

REPO_DIR="$(cd "$(dirname "$0")/.." && pwd)"
MEMORY_DIR="${REPO_DIR}/.planning/memory"
CLAUDE_MEMORY_DIR="${HOME}/.claude/projects/-Users-bence-code-QuantumFurnace-jl/memory"

# Ensure the project-local memory directory exists
mkdir -p "${MEMORY_DIR}"

# Ensure the parent of the Claude auto-memory path exists
mkdir -p "$(dirname "${CLAUDE_MEMORY_DIR}")"

# Create or update the symlink (idempotent)
if [ -L "${CLAUDE_MEMORY_DIR}" ]; then
    # Already a symlink — verify it points to the right place
    current="$(readlink "${CLAUDE_MEMORY_DIR}")"
    if [ "${current}" != "${MEMORY_DIR}" ]; then
        ln -sfn "${MEMORY_DIR}" "${CLAUDE_MEMORY_DIR}"
    fi
elif [ -d "${CLAUDE_MEMORY_DIR}" ]; then
    # Real directory exists — migrate any contents, then replace with symlink
    if [ "$(ls -A "${CLAUDE_MEMORY_DIR}" 2>/dev/null)" ]; then
        cp -rn "${CLAUDE_MEMORY_DIR}"/* "${MEMORY_DIR}/" 2>/dev/null || true
    fi
    rm -rf "${CLAUDE_MEMORY_DIR}"
    ln -sfn "${MEMORY_DIR}" "${CLAUDE_MEMORY_DIR}"
else
    # Nothing there yet — just create the symlink
    ln -sfn "${MEMORY_DIR}" "${CLAUDE_MEMORY_DIR}"
fi
