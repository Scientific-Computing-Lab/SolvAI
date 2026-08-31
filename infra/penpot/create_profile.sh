#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")"
umask 077

credentials_file=credentials.env
if [[ ! -f "$credentials_file" ]]; then
  profile_password=$(python3 -c 'import secrets; print(secrets.token_urlsafe(24))')
  {
    printf 'PENPOT_PROFILE_EMAIL=solvai@local.invalid\n'
    printf 'PENPOT_PROFILE_FULLNAME=SolvAI Figure Editor\n'
    printf 'PENPOT_PROFILE_PASSWORD=%s\n' "$profile_password"
  } > "$credentials_file"
fi

profile_email=$(sed -n 's/^PENPOT_PROFILE_EMAIL=//p' "$credentials_file")
profile_fullname=$(sed -n 's/^PENPOT_PROFILE_FULLNAME=//p' "$credentials_file")
profile_password=$(sed -n 's/^PENPOT_PROFILE_PASSWORD=//p' "$credentials_file")

docker compose exec -T penpot-backend python3 manage.py create-profile \
  --skip-tutorial \
  --skip-walkthrough \
  --email "$profile_email" \
  --fullname "$profile_fullname" \
  --password "$profile_password"
