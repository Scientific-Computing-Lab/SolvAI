#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")"
umask 077

if [[ ! -f .env ]]; then
  penpot_secret=$(python3 -c 'import secrets; print(secrets.token_urlsafe(64))')
  db_password=$(python3 -c 'import secrets; print(secrets.token_urlsafe(36))')
  {
    printf 'PENPOT_VERSION=2.17.2\n'
    printf 'PENPOT_HTTP_PORT=9001\n'
    printf 'PENPOT_PUBLIC_URI=http://localhost:9001\n'
    printf 'PENPOT_SECRET_KEY=%s\n' "$penpot_secret"
    printf 'PENPOT_DB_PASSWORD=%s\n' "$db_password"
  } > .env
fi

docker compose pull
docker compose up -d
docker compose ps
