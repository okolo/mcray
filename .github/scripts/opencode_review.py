#!/usr/bin/env python3
"""Publish an OpenCode verdict only for the revision that was reviewed."""

import json
import os
from pathlib import Path
import subprocess
import sys


MARKER = "<!-- opencode-independent-review -->"


def api(endpoint, method="GET", payload=None):
    command = ["gh", "api", endpoint, "--method", method]
    if payload is not None:
        command += ["--input", "-"]
    result = subprocess.run(
        command, input=json.dumps(payload) if payload is not None else None,
        text=True, capture_output=True, check=True,
    )
    return json.loads(result.stdout)


def dismiss(endpoint, review_id, reason, head_sha):
    try:
        api(f"{endpoint}/reviews/{review_id}/dismissals", "PUT", {"message": reason})
    except subprocess.CalledProcessError as error:
        # Protected branches can reserve dismissal for administrators even
        # though this bot has permission to submit its own reviews. Supersede
        # its approval with a temporary negative verdict in that case.
        if "HTTP 403" not in (error.stderr or ""):
            raise
        api(f"{endpoint}/reviews", "POST", {
            "commit_id": head_sha, "event": "REQUEST_CHANGES",
            "body": f"{MARKER}\nApproval withdrawn pending independent re-review.\n\n{reason}",
        })


def prepare(repository, number, context_path, output_path):
    endpoint = f"repos/{repository}/pulls/{number}"
    pr = api(endpoint)
    allowed = (pr["state"] == "open" and not pr["draft"]
               and pr["head"]["repo"] is not None
               and pr["head"]["repo"]["full_name"] == repository)
    with output_path.open("a") as output:
        output.write(f"allowed={'true' if allowed else 'false'}\n")
    if not allowed:
        print("Review skipped: expected an open, non-draft, same-repository PR.")
        return
    # Preparing a rerun must not change any existing review. GitHub dismisses
    # stale approvals on new changes; a later decisive verdict supersedes the
    # previous verdict for the same head.
    context_path.write_text(json.dumps({
        "repository": repository, "number": number,
        "head_sha": pr["head"]["sha"], "base_sha": pr["base"]["sha"],
    }))
    # Step outputs are captured by the runner before the model starts. The
    # publisher uses them, not the model-accessible context file, as its anchor.
    with output_path.open("a") as output:
        output.write(f"head_sha={pr['head']['sha']}\nbase_sha={pr['base']['sha']}\n")


def review_payload(context, verdict, current):
    fields = {"head_sha", "base_sha", "review_complete", "verdict", "blocking_findings", "summary"}
    if not isinstance(verdict, dict) or set(verdict) != fields:
        raise ValueError("Missing or unexpected fields in the OpenCode verdict")
    if (verdict["head_sha"] != context["head_sha"]
            or verdict["base_sha"] != context["base_sha"]):
        raise ValueError("The model did not review the frozen head/base SHAs")
    if (current["state"] != "open" or current["draft"]
            or current["head"]["repo"] is None
            or current["head"]["repo"]["full_name"] != context["repository"]
            or current["head"]["sha"] != context["head_sha"]
            or current["base"]["sha"] != context["base_sha"]):
        raise ValueError("PR state or head/base SHA changed during review; rerun required")
    blockers = verdict["blocking_findings"]
    if (not isinstance(blockers, list)
            or any(not isinstance(item, str) or not item.strip() for item in blockers)):
        raise ValueError("blocking_findings must be a list of nonempty strings")
    if type(verdict["review_complete"]) is not bool:
        raise ValueError("review_complete must be a JSON boolean")
    summary = verdict["summary"]
    if not isinstance(summary, str) or not summary.strip() or len(summary) > 50000:
        raise ValueError("Expected a nonempty review summary of at most 50000 characters")
    decision = verdict["verdict"]
    if decision not in ("APPROVE", "REQUEST_CHANGES", "COMMENT"):
        raise ValueError("Unknown review verdict")
    if decision == "APPROVE" and (not verdict["review_complete"] or blockers):
        raise ValueError("Cannot approve an incomplete review or blocking findings")
    if decision == "REQUEST_CHANGES" and not blockers:
        raise ValueError("REQUEST_CHANGES requires concrete blocking findings")
    if blockers and decision != "REQUEST_CHANGES":
        raise ValueError("Blocking findings require REQUEST_CHANGES")
    body = (f"{MARKER}\n## Independent OpenCode review — GLM-5.3 via OpenRouter\n\n"
            f"Head: `{context['head_sha']}`\nBase: `{context['base_sha']}`\n\n{summary}")
    if blockers:
        body += "\n\n### Blocking findings\n\n" + "\n".join(f"- {item}" for item in blockers)
    if len(body) > 60000:
        raise ValueError("Review body exceeds the publication limit")
    return {"commit_id": context["head_sha"], "event": decision, "body": body}


def publish(context, verdict_path):
    verdict = json.loads(verdict_path.read_text())
    endpoint = f"repos/{context['repository']}/pulls/{context['number']}"
    payload = review_payload(context, verdict, api(endpoint))
    review = api(f"{endpoint}/reviews", "POST", payload)
    # The REST create-review API has no compare-and-swap guard. Detect a push
    # racing publication and withdraw this approval rather than approve new code.
    if payload["event"] == "APPROVE":
        try:
            current = api(endpoint)
            if (current["head"]["sha"] != context["head_sha"]
                    or current["base"]["sha"] != context["base_sha"]):
                raise ValueError("PR changed during approval publication")
        except Exception:
            dismiss(endpoint, review["id"], "Could not confirm the reviewed head/base after publication.",
                    context["head_sha"])
            raise
    print(f"Published {payload['event']}: {review['html_url']}")


def main():
    number = int(os.environ["PR_NUMBER"])
    if number <= 0:
        raise ValueError("Expected a positive PR number")
    if sys.argv[1:] == ["prepare"]:
        prepare(os.environ["GITHUB_REPOSITORY"], number,
                Path(os.environ["OPENCODE_REVIEW_CONTEXT"]), Path(os.environ["GITHUB_OUTPUT"]))
    elif sys.argv[1:] == ["publish"]:
        context = {
            "repository": os.environ["GITHUB_REPOSITORY"], "number": number,
            "head_sha": os.environ["OPENCODE_EXPECTED_HEAD_SHA"],
            "base_sha": os.environ["OPENCODE_EXPECTED_BASE_SHA"],
        }
        publish(context, Path(os.environ["OPENCODE_REVIEW_RESULT"]))
    else:
        raise ValueError("Expected prepare or publish")


if __name__ == "__main__":
    main()
