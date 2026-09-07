# Private rrad CI credential

The `Private rrad integration` workflow performs an uncached installation of
the exact `rrad` release pinned in `DESCRIPTION`. GitHub's workflow token for
SIRAP cannot read a separate private repository, so the workflow expects one
repository secret named `RRAD_READ_TOKEN`.

Because `rrad` is owned by the personal account `indianewok`, setup requires
the `rrad` owner and a SIRAP administrator:

1. The `indianewok` account owner creates a fine-grained GitHub token whose
   resource owner is `indianewok`.
2. Limit the token to the single `rrad` repository, grant only
   **Contents: read**, and choose an appropriate expiration.
3. Transfer the token once through an approved secret-sharing channel; do not
   send it in an issue, pull request, email, or chat history.
4. A `hjm2217/SIRAP` administrator adds it under **Settings > Secrets and
   variables > Actions** with the name `RRAD_READ_TOKEN`.
5. Run the workflow manually once and confirm that the pinned-version and
   native-core handshake steps pass.

Do not substitute a maintainer's broad classic personal token. The
secret-bearing workflow runs only after a push to `main` or a trusted manual
dispatch; it is not exposed to pull-request code.
