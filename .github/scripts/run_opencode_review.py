#!/usr/bin/env python3
"""Run analysis without the GitHub command's automatic issue comments."""

import os
import subprocess


def run_review(executable, prompt, expected_head):
    environment = os.environ.copy()
    for key in ("GH_TOKEN", "GITHUB_TOKEN", "GH_ENTERPRISE_TOKEN", "GITHUB_ENTERPRISE_TOKEN"):
        environment.pop(key, None)
    subprocess.run([
        executable, "run", "--model", "openrouter/deepseek/deepseek-v4.1-flash",
        "--agent", "reviewer", "--print-logs", "--log-level", "INFO", prompt,
    ], env=environment, check=True)
    head = subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip()
    if head != expected_head:
        raise ValueError("Review changed the frozen checkout revision")
    if subprocess.check_output(["git", "status", "--porcelain"], text=True).strip():
        raise ValueError("Review left tracked changes or untracked files in the checkout")


if __name__ == "__main__":
    run_review(os.environ["OPENCODE_REVIEW_CLI"], os.environ["OPENCODE_REVIEW_PROMPT"],
               os.environ["OPENCODE_EXPECTED_HEAD_SHA"])
