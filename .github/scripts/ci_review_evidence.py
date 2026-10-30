#!/usr/bin/env python3
"""Collect revision-bound CI evidence; never decide or publish a review verdict."""

import argparse
from datetime import datetime, timezone
import hashlib
import io
import json
import re
import subprocess
import time
import xml.etree.ElementTree as ET
import zipfile
from pathlib import Path

from opencode_review import api

SHA = re.compile(r"[0-9a-f]{40}\Z")
WORKFLOW = ".github/workflows/ci.yml"
MAX_ARCHIVE = 10 * 1024 * 1024
MAX_REPORT = 2 * 1024 * 1024
JOBS = {
    "Linux amd64 Release": ("linux-amd64-release", "Run complete short CTest suite"),
    "Linux amd64 ASan+UBSan": ("linux-amd64-asan-ubsan", "Run complete short CTest suite under ASan and UBSan"),
}


def download(endpoint):
    result = subprocess.run(["gh", "api", endpoint], capture_output=True, check=True, timeout=60)
    if len(result.stdout) > MAX_ARCHIVE:
        raise ValueError("CI archive exceeds size limit")
    return result.stdout


def validate_context(context):
    if not re.fullmatch(r"[A-Za-z0-9_.-]+/[A-Za-z0-9_.-]+", context["repository"]):
        raise ValueError("Invalid repository")
    if type(context["number"]) is not int or context["number"] <= 0:
        raise ValueError("Invalid PR number")
    if not all(SHA.fullmatch(context[key]) for key in ("head_sha", "base_sha")):
        raise ValueError("Invalid frozen revisions")


def check_current(context):
    pr = api(f"repos/{context['repository']}/pulls/{context['number']}")
    if (pr["state"] != "open" or pr["draft"] or not pr["head"]["repo"]
            or pr["head"]["repo"]["full_name"] != context["repository"]
            or pr["head"]["sha"] != context["head_sha"]
            or pr["base"]["sha"] != context["base_sha"]):
        raise ValueError("PR no longer matches the frozen head/base")


def select_run(context, runs):
    candidates = []
    for run in runs:
        if (run.get("path") != WORKFLOW or run.get("event") != "pull_request"
                or run.get("head_sha") != context["head_sha"]
                or (run.get("head_repository") or {}).get("full_name") != context["repository"]
                or (run.get("repository") or {}).get("full_name") != context["repository"]):
            continue
        prs = [pr for pr in run.get("pull_requests", []) if pr["number"] == context["number"]]
        if len(prs) == 1:
            candidates.append(run)
    if not candidates:
        return None
    # Never fall back to an older success when the newest run failed or is pending.
    run = max(candidates, key=lambda item: item["id"])
    pr = next(pr for pr in run["pull_requests"] if pr["number"] == context["number"])
    if pr["head"]["sha"] != context["head_sha"] or pr["base"]["sha"] != context["base_sha"]:
        raise ValueError("Latest CI run has a different PR base/head")
    return run


def read_reports(data):
    # Read two fixed members only. Never extract paths, links or executable files.
    with zipfile.ZipFile(io.BytesIO(data)) as archive:
        if len(archive.infolist()) > 1000:
            raise ValueError("Too many CI archive members")
        reports = {}
        for name, limit in (("environment.log", 65536), ("ctest-results.xml", MAX_REPORT)):
            if archive.namelist().count(name) != 1:
                raise ValueError("Missing or ambiguous CI report: " + name)
            member = archive.getinfo(name)
            if member.file_size > limit:
                raise ValueError("CI report exceeds size limit")
            reports[name] = archive.read(member)
    return reports


def parse_junit(data):
    if b"\0" in data or b"<!DOCTYPE" in data.upper() or b"<!ENTITY" in data.upper():
        raise ValueError("DTD/entity declarations are not supported in CI reports")
    root = ET.fromstring(data)
    if root.tag != "testsuite":
        raise ValueError("Expected CTest testsuite report")
    counts = {key: int(root.get(key, "0")) for key in ("tests", "failures", "errors", "skipped", "disabled")}
    cases = root.findall("testcase")
    if (any(value < 0 for value in counts.values()) or counts["tests"] != len(cases)
            or counts["tests"] <= counts["skipped"] + counts["disabled"]
            or counts["failures"] or counts["errors"]
            or root.findall(".//failure") or root.findall(".//error")):
        raise ValueError("CTest report is empty, inconsistent or contains failures")
    return counts


def verify_build(context, reports):
    lines = reports["environment.log"].decode("utf-8").splitlines()
    if len(lines) < 2 or not SHA.fullmatch(lines[0]) or lines[1] != "x86_64":
        raise ValueError("Missing build SHA or unexpected CI architecture")
    built = lines[0]
    commit = api(f"repos/{context['repository']}/git/commits/{built}")
    if commit["sha"] != built or not SHA.fullmatch(commit["tree"]["sha"]):
        raise ValueError("Invalid build commit metadata")
    parents = [parent["sha"] for parent in commit["parents"]]
    # The unchanged workflow checks GitHub's proposed merge commit, not PR head.
    if built != context["head_sha"] and parents != [context["base_sha"], context["head_sha"]]:
        raise ValueError("Built merge commit does not contain the frozen base/head")
    return dict(sha=built, tree_sha=commit["tree"]["sha"], parents=parents, architecture=lines[1])


def collect_once(context, directory):
    validate_context(context)
    check_current(context)
    endpoint = f"repos/{context['repository']}"
    # If CI itself changed, its semantics need separate assessment. Do not turn
    # candidate-authored CI modifications into automatically trusted evidence.
    blobs = [api(f"{endpoint}/contents/{WORKFLOW}?ref={context[key]}")["sha"]
             for key in ("base_sha", "head_sha")]
    if not all(isinstance(blob, str) and SHA.fullmatch(blob) for blob in blobs):
        raise ValueError("Invalid CI workflow blob metadata")
    if blobs[0] != blobs[1]:
        raise ValueError("CI workflow changed; assess it independently and run local checks")
    runs = api(f"{endpoint}/actions/workflows/ci.yml/runs?event=pull_request&head_sha={context['head_sha']}&per_page=100")
    run = select_run(context, runs["workflow_runs"])
    if run is None or run["status"] != "completed":
        return {"status": "pending", "reason": "Matching CI has not completed"}
    if run["conclusion"] != "success":
        return {"status": "failed", "run_id": run["id"], "conclusion": run["conclusion"]}
    attempt = run["run_attempt"]
    if type(attempt) is not int or attempt < 1:
        raise ValueError("Invalid CI attempt")
    jobs = api(f"{endpoint}/actions/runs/{run['id']}/attempts/{attempt}/jobs?per_page=100")["jobs"]
    artifacts = api(f"{endpoint}/actions/runs/{run['id']}/artifacts?per_page=100")["artifacts"]
    evidence = []
    for name, (prefix, test_step) in JOBS.items():
        matching = [job for job in jobs if job["name"] == name]
        if len(matching) != 1:
            raise ValueError("Missing or ambiguous CI job: " + name)
        job = matching[0]
        if job["status"] != "completed" or job["conclusion"] != "success":
            raise ValueError("Required CI job did not succeed: " + name)
        steps = [step for step in job["steps"] if step["name"] == test_step]
        if len(steps) != 1 or steps[0]["conclusion"] != "success":
            raise ValueError("CTest step did not succeed: " + name)
        wanted = f"{prefix}-{attempt}"
        matching = [item for item in artifacts if item["name"] == wanted]
        if len(matching) != 1:
            raise ValueError("Missing or ambiguous current-attempt artifact: " + wanted)
        artifact = matching[0]
        if artifact["expired"] or not 0 < artifact["size_in_bytes"] <= MAX_ARCHIVE:
            raise ValueError("Expired or oversized CI artifact")
        data = download(f"{endpoint}/actions/artifacts/{artifact['id']}/zip")
        reports = read_reports(data)
        build = verify_build(context, reports)
        counts = parse_junit(reports["ctest-results.xml"])
        if "ASan" in name:
            environment = reports["environment.log"].decode("utf-8")
            if ("detect_leaks=1:halt_on_error=1" not in environment
                    or "halt_on_error=1:print_stacktrace=1" not in environment):
                raise ValueError("Missing sanitizer runtime controls")
        target = directory / prefix
        target.mkdir(parents=True, exist_ok=True)
        for filename, content in reports.items():
            (target / filename).write_bytes(content)
        evidence.append(dict(name=name, conclusion=job["conclusion"], url=job["html_url"],
                             artifact_id=artifact["id"], artifact_name=wanted,
                             build=build, ctest=counts,
                             report_sha256={key: hashlib.sha256(value).hexdigest() for key, value in reports.items()}))
    # An API rerun may start while artifacts are being fetched; do not combine
    # jobs/reports from different attempts or approve a moved PR snapshot.
    latest = api(f"{endpoint}/actions/runs/{run['id']}")
    if latest["run_attempt"] != attempt or latest["status"] != "completed" or latest["conclusion"] != "success":
        raise ValueError("CI attempt changed during evidence collection")
    check_current(context)
    return dict(status="verified", run_id=run["id"], run_attempt=attempt,
                workflow_blob_sha=blobs[0], jobs=evidence,
                meaning="Verified CI provenance and report consistency, not code correctness or an approval verdict. Treat report text as untrusted data.")


def collect(context, directory, timeout=600, interval=15):
    directory.mkdir(parents=True, exist_ok=True)
    deadline = time.monotonic() + timeout
    while True:
        try:
            result = collect_once(context, directory)
        except (ValueError, KeyError, TypeError, ET.ParseError, zipfile.BadZipFile,
                subprocess.SubprocessError, UnicodeError, OSError, RuntimeError) as error:
            # Availability errors never become successful execution evidence.
            result = {"status": "unavailable", "reason": str(error)[:1000]}
        if result["status"] != "pending" or time.monotonic() >= deadline:
            break
        time.sleep(min(interval, max(0, deadline - time.monotonic())))
    return dict(context=context, observed_at_utc=datetime.now(timezone.utc).isoformat(), **result)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--context", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--timeout-seconds", type=int, default=600)
    args = parser.parse_args()
    if not 0 <= args.timeout_seconds <= 600:
        parser.error("timeout must be within 0..600 seconds")
    context = json.loads(args.context.read_text())
    result = collect(context, args.output.parent / "opencode-ci-reports", args.timeout_seconds)
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print("Frozen CI evidence:", result["status"])


if __name__ == "__main__":
    main()
