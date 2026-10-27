# Independent OpenCode review

The [review workflow](../../.github/workflows/opencode-review.yml) uses
`openrouter/z-ai/glm-5.3` to independently review pull requests. It runs on
GitHub-hosted Ubuntu runners, with a 25-minute job timeout. It does not replace
the project's numerical validation or human review.

## Repository setup

1. Add an OpenRouter API key as the repository Actions secret
   `OPENROUTER_API_KEY` in
   [Actions secrets](https://github.com/okolo/mcray/settings/secrets/actions).
   Keep the key out of source files, comments, logs, and chat messages.
   OpenRouter usage is billed to the account owning this key.
2. Merge the workflow into the default branch, `main`. Comment-triggered
   workflows must exist there before GitHub will dispatch them.
3. Open or update a non-draft PR from a branch in this repository. Check
   [Actions](https://github.com/okolo/mcray/actions) and the resulting PR comment.
   Confirm the reviewed head SHA matches the PR's current head.

This setup uses the automatic `GITHUB_TOKEN`, not a personal access token or
the OpenCode GitHub App. The token can read code and write PR comments and
reactions, including inline review comments. GitHub requires the PR write
permission for these operations; the token has no code-write permission.
No OIDC permission is requested.
Automatic session sharing is disabled.

## Triggers and scope

Review runs when a non-draft same-repository PR is opened, updated, reopened,
or marked ready for review. Repository owners, members, and collaborators can
request another review with `/opencode` or `/oc` in a PR comment or an inline
review comment. These commands always request a review, including when their
text asks for implementation. Issue comments outside PRs are ignored.

Automatic fork-PR runs are skipped because GitHub does not expose repository
API secrets to them. A trusted participant can request review through a PR
comment; this runs the workflow from the default branch with repository
credentials. Only request this for code you trust the agent to inspect with
those credentials. This workflow does not use `pull_request_target`.

One review job per PR is retained at a time; a newer eligible trigger cancels
an older review. Ordinary comments and other skipped events do not cancel it.
The review reports head/base SHAs, concrete findings, checks actually executed,
and remaining verification gaps. A successful Actions job means the agent ran,
not that the PR is correct. The workflow is advisory and does not configure a
required branch-protection check or approve/merge PRs.

## Maintenance and verification

The checkout and OpenCode Action references are pinned to commit SHAs. The
upstream OpenCode Action itself downloads the latest OpenCode CLI and uses
its own dependencies, so the complete runtime is not immutable. Recheck
upstream behavior when updating it. The job requires outbound access to
GitHub, OpenCode installation endpoints, and OpenRouter.

Local YAML/static checks do not verify provider credentials, model access,
GitHub comment permissions, or review quality. Verify those with an actual
PR after publication and secret setup; do not describe the integration as
operational until that run completes and its comment is checked.

References: [OpenCode GitHub integration](https://opencode.ai/docs/github/),
[OpenRouter models](https://openrouter.ai/models?q=z-ai).
