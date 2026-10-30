# Independent OpenCode review

The [review workflow](../../.github/workflows/opencode-review.yml) uses
`openrouter/deepseek/deepseek-v4.1-flash` as its only model. It runs on
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

### Model routing

The `reviewer` performs all code and architectural analysis directly using
DeepSeek V4.1 Flash. Task/subagent delegation is denied globally and for the
reviewer; no architectural subagent or consultation result is configured.
The final verdict has only `head_sha`, `base_sha`, `review_complete`, `verdict`,
`blocking_findings`, and `summary`. The formal heading remains `AI review`.

The main `model` and auxiliary `small_model` both explicitly select
`openrouter/deepseek/deepseek-v4.1-flash`. The built-in `title` agent is disabled
because CI does not need generated conversation titles; otherwise it can use
an automatically selected small model and bill that request to the same API
key. The `compaction` and `summary` agents explicitly use DeepSeek as well.
Only the OpenRouter provider is enabled, and its model whitelist contains
only `deepseek/deepseek-v4.1-flash`, so other catalog models fail resolution
rather than silently receive API calls. These are OpenCode v1 settings;
revalidate their equivalents before a runtime migration.

This policy limits this OpenCode runtime, not every other application using
the API key. Verify actual `providerID`/`modelID`/`agent` stream records in
Actions logs after updates; the set of billable models should contain only
OpenRouter DeepSeek V4.1 Flash, with no `title` stream or subagent invocation.
No new secret or GitHub permission is required. Inline configuration survives
PR branch switches. Set both `agent: reviewer` and `default_agent: reviewer`:
OpenCode v1.19.5's GitHub handler selects the configured default.

The reviewer writes `OPENCODE_REVIEW_RESULT`; the existing publisher validates
SHAs and submits the single formal verdict. Failed/incomplete review cannot
create a new approval and does not preemptively withdraw an existing approval
for unchanged code.

### Events

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
Starting a new review never withdraws an existing approval or submits a
temporary negative review. If the head is unchanged, a failed/inconclusive
rerun preserves the previous approval. A new completed `REQUEST_CHANGES`
verdict supersedes an earlier approval. GitHub's stale-review dismissal rule
handles new code changes. Human reviews are never modified.

Only a newly published approval that fails the final revision check is
withdrawn. If GitHub denies that dismissal, the publisher supersedes its own
invalid approval with `REQUEST_CHANGES`, explaining that the reviewed revision
could not be confirmed. This is an exceptional publication-race safeguard,
not a preliminary step before review. The workflow
never merges PRs or bypasses branch protection, and it does not add a required
status check. A successful job can also mean `COMMENT` or `REQUEST_CHANGES`;
inspect the formal verdict rather than treating the job result as approval.

The publisher is copied outside the checkout before the agent starts. This
preserves the publication code even when OpenCode switches to an older PR
branch that does not yet contain the integration scripts.
Runtime permissions explicitly allow the runner's temporary directory and
`/tmp` for context, verdicts, and scratch work. Other external directories
and interactive question/loop prompts are denied so CI cannot hang waiting
for approval to read the frozen context or write its verdict.

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
[OpenCode agents](https://opencode.ai/docs/agents/) describes per-agent models
and task permissions.
