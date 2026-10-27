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
2. In Actions settings, enable "Allow GitHub Actions to create and approve
   pull requests". The workflow already grants `pull-requests: write`.
   For `main`, retain at least one required approval and enable dismissal of
   stale approving reviews on new changes. An approval by `github-actions[bot]`
   is distinct from the PR author and counts toward the ordinary required
   approval count. Separate code-owner requirements still require code owners;
   a GitHub App bot is not a user or team that can be assigned in CODEOWNERS.
3. Merge the workflow into the default branch, `main`. Comment-triggered
   workflows must exist there before GitHub will dispatch them.
4. Open or update a non-draft PR from a branch in this repository. Check
   [Actions](https://github.com/okolo/mcray/actions), the detailed PR comment,
   and the formal review. Confirm the reviewed head SHA matches the PR's head.

This setup uses the automatic `GITHUB_TOKEN`, not a personal access token or
the OpenCode GitHub App. The token can read code and write PR comments and
reactions, including inline review comments. GitHub requires the PR write
permission for these operations; the token has no code-write permission.
No OIDC permission is requested.
Automatic session sharing is disabled.

## Triggers and scope

Review runs when a non-draft same-repository PR is opened, updated, reopened,
or marked ready for review. Repository owners, members, and collaborators
with write or admin access can request another review in a PR comment or an
inline review comment. Start the comment with exactly `/opencode` or `/oc`,
optionally followed by a space and review context. Mentions inside URLs or
other prose do not trigger review. These commands always request a review,
including when their text asks for implementation. Issue comments outside
PRs are ignored. The upstream action verifies the requester's permissions.

Fork PRs are excluded from automatic review and from comment-triggered review.
For comment triggers, a metadata-only step checks the PR's source repository
before checkout or passing provider credentials to OpenCode. The agent can
execute checks with credentials in its environment, so same-repository code
and its instructions must still be trusted. This workflow does not use
`pull_request_target`. Draft and closed PRs are also excluded.

One review job per PR is retained at a time; a newer eligible trigger cancels
an older review. Ordinary comments and other skipped events do not cancel it.
The review reports head/base SHAs, concrete findings, checks actually executed,
and remaining verification gaps. The model also writes a structured verdict
outside the checkout. The publisher validates it and submits a formal review:

- `APPROVE`: completed independent review with no blocking findings.
- `REQUEST_CHANGES`: concrete blocking defects, listed in the review.
- `COMMENT`: incomplete or inconclusive analysis, without approval.

Passing build checks alone never establishes approval. Missing or malformed
verdicts, failed agent runs, and changed head/base SHAs prevent publication.
The publisher pins the review to the frozen head SHA and checks again after
approval, withdrawing it if publication raced a change or the final check fails.
An earlier marked approval from this bot is withdrawn before a new review
starts, so a failed/inconclusive rerun cannot retain its previous approval.
Human reviews and other automation's approvals are preserved. The workflow
never merges PRs or bypasses branch protection, and it does not add a required
status check. A successful job can also mean `COMMENT` or `REQUEST_CHANGES`;
inspect the formal verdict rather than treating the job result as approval.

The review must leave the checkout clean. Build directories, caches, logs,
and test outputs belong outside it. Otherwise the upstream action tries to
commit/push the generated files and fails with these restricted credentials,
potentially losing the review. The agent writes its JSON verdict only to
`OPENCODE_REVIEW_RESULT`; it must not submit formal reviews itself. Only the
publisher should do that after validation. The prompt also forbids merges
and changes to PR state or metadata. These are instructions to the agent,
not restrictions on the token, which needs PR write access for the integration.

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

Run the publication guard tests without creating checkout artifacts:

```sh
PYTHONDONTWRITEBYTECODE=1 python3 -m unittest discover -s .github/scripts -p 'test_*.py' -v
```

References: [OpenCode GitHub integration](https://opencode.ai/docs/github/),
[OpenRouter models](https://openrouter.ai/models?q=z-ai),
[GitHub review API](https://docs.github.com/en/rest/pulls/reviews#create-a-review-for-a-pull-request).
