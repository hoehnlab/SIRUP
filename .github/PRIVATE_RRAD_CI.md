# Private rrad CI credential

The `Private rrad integration` workflow performs an uncached installation of
the exact `rrad` release pinned in `DESCRIPTION`. GitHub's workflow token for
SIRAP cannot read a separate private repository, so the workflow expects one
repository secret named `RRAD_READ_TOKEN`.

An administrator should:

1. Create a fine-grained GitHub token under an account authorized for
   `indianewok/rrad`.
2. Limit repository access to `indianewok/rrad` and grant only
   **Contents: read**.
3. Add it to `hjm2217/SIRAP` under **Settings > Secrets and variables >
   Actions** with the name `RRAD_READ_TOKEN`.
4. Run the workflow manually once and confirm that the pinned-version and
   native-core handshake steps pass.

Do not substitute a maintainer's broad classic personal token. The
secret-bearing workflow runs only after a push to `main` or a trusted manual
dispatch; it is not exposed to pull-request code.
