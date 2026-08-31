# Self-hosted Penpot for SolvAI Figure 1

This directory pins the official Penpot application and first-party MCP service to
Penpot `2.17.2`. The browser application is bound to loopback only at
`http://127.0.0.1:9001`; no unauthenticated HTTP service is exposed publicly.

Persistent state is held in Docker volumes:

- `solvai_penpot_postgres` — accounts, projects, files and native design objects;
- `solvai_penpot_assets` — imported SVGs and exported design assets.

Local secrets are generated into `.env` and `credentials.env`, both ignored by Git.

```bash
cd infra/penpot
./bootstrap.sh
./create_profile.sh
```

For access from another workstation, create an SSH tunnel:

```bash
ssh -L 9001:127.0.0.1:9001 USER@SERVER
```

Then open `http://localhost:9001`. Read the local credentials on the server with
`sed -n '1,3p' infra/penpot/credentials.env`; do not commit or paste them into logs.

The public URI and loopback binding are intentionally HTTP-only because access is
through an encrypted SSH tunnel. A public deployment would require a real domain,
TLS reverse proxy, secure cookies and email verification.
