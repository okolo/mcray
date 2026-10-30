"""Revision/attempt/artifact guards for CI evidence; no live service writes."""

from copy import deepcopy
import io
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import zipfile

import ci_review_evidence as evidence

HEAD, BASE, MERGE, TREE = (character * 40 for character in "abcd")
CONTEXT = dict(repository="owner/repo", number=12, head_sha=HEAD, base_sha=BASE)


def run(identifier=20):
    return dict(id=identifier, path=evidence.WORKFLOW, event="pull_request",
                head_sha=HEAD, repository={"full_name": "owner/repo"},
                head_repository={"full_name": "owner/repo"}, run_attempt=2,
                status="completed", conclusion="success",
                pull_requests=[dict(number=12, head={"sha": HEAD}, base={"sha": BASE})])


def archive(sha=MERGE, xml=None, sanitizer=True):
    xml = xml or b'<testsuite tests="1"><testcase name="real-test"/></testsuite>'
    env = sha + "\nx86_64\ncompiler-version\n"
    if sanitizer:
        env += "detect_leaks=1:halt_on_error=1\nhalt_on_error=1:print_stacktrace=1\n"
    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, "w") as zipped:
        zipped.writestr("environment.log", env)
        zipped.writestr("ctest-results.xml", xml)
        zipped.writestr("../../not-extracted", "untrusted")
    return buffer.getvalue()


class EvidenceTests(unittest.TestCase):
    def setUp(self):
        self.pr = dict(state="open", draft=False, head=dict(sha=HEAD, repo={"full_name": "owner/repo"}), base=dict(sha=BASE))
        self.run = run()
        self.commit = dict(sha=MERGE, tree={"sha": TREE}, parents=[{"sha": BASE}, {"sha": HEAD}])
        self.artifacts = [dict(id=index+1, name=prefix+"-2", size_in_bytes=1000, expired=False)
                          for index, (prefix, _) in enumerate(evidence.JOBS.values())]
        self.jobs = [dict(name=name, status="completed", conclusion="success", html_url="https://example/job",
                          steps=[dict(name=step, conclusion="success")])
                     for name, (_, step) in evidence.JOBS.items()]

    def api(self, endpoint):
        if endpoint.endswith("/pulls/12"): return self.pr
        if "/contents/" in endpoint: return {"sha": "e"*40}
        if "/workflows/ci.yml/runs?" in endpoint: return {"workflow_runs": [self.run]}
        if "/attempts/2/jobs?" in endpoint: return {"jobs": self.jobs}
        if "/runs/20/artifacts?" in endpoint: return {"artifacts": self.artifacts}
        if "/git/commits/" in endpoint: return self.commit
        if endpoint.endswith("/runs/20"): return self.run
        raise AssertionError(endpoint)

    def collect(self):
        with tempfile.TemporaryDirectory() as directory:
            with patch.object(evidence, "api", side_effect=self.api), patch.object(evidence, "download", return_value=archive()):
                return evidence.collect(CONTEXT, Path(directory), timeout=0)

    def test_verified_ci_has_both_jobs_and_real_merge_revision(self):
        result = self.collect()
        self.assertEqual(result["status"], "verified")
        self.assertEqual(result["run_attempt"], 2)
        self.assertEqual(len(result["jobs"]), 2)
        self.assertEqual(result["jobs"][0]["build"]["sha"], MERGE)
        self.assertEqual(result["jobs"][0]["ctest"]["tests"], 1)
        self.assertNotIn("verdict", result)

    def test_selection_ignores_wrong_repo_workflow_event_and_head(self):
        for change in ({"path": ".github/workflows/other.yml"}, {"event": "push"},
                       {"head_sha": "f"*40}, {"head_repository": {"full_name": "fork/repo"}},
                       {"repository": {"full_name": "other/repo"}}, {"pull_requests": []}):
            with self.subTest(change=change):
                self.assertIsNone(evidence.select_run(CONTEXT, [run() | change]))

    def test_newest_pending_or_failed_run_does_not_reuse_old_success(self):
        newer = run(21) | {"status": "in_progress", "conclusion": None}
        self.assertEqual(evidence.select_run(CONTEXT, [run(), newer])["id"], 21)
        self.run["status"] = "in_progress"
        self.assertEqual(self.collect()["status"], "pending")
        self.run.update(status="completed", conclusion="failure")
        self.assertEqual(self.collect()["status"], "failed")

    def test_changed_base_or_live_pr_cannot_supply_evidence(self):
        self.run["pull_requests"][0]["base"]["sha"] = "f"*40
        self.assertEqual(self.collect()["status"], "unavailable")
        self.run = run()
        self.pr["head"]["sha"] = "f"*40
        self.assertEqual(self.collect()["status"], "unavailable")

    def test_ci_workflow_changes_are_not_automatically_trusted(self):
        def changed(endpoint):
            if "/contents/" in endpoint: return {"sha": endpoint[-40:]}
            return self.api(endpoint)
        with tempfile.TemporaryDirectory() as directory, patch.object(evidence, "api", side_effect=changed):
            self.assertEqual(evidence.collect(CONTEXT, Path(directory), 0)["status"], "unavailable")

    def test_merge_parents_and_tree_must_match(self):
        for parents in ([{"sha": BASE}, {"sha": "f"*40}], [{"sha": HEAD}, {"sha": BASE}], []):
            self.commit["parents"] = parents
            with self.subTest(parents=parents):
                self.assertEqual(self.collect()["status"], "unavailable")
        self.commit.update(sha=HEAD, parents=[])
        with patch.object(evidence, "api", return_value=self.commit):
            self.assertEqual(evidence.verify_build(CONTEXT, evidence.read_reports(archive(HEAD)))["sha"], HEAD)

    def test_missing_failed_or_ambiguous_jobs_and_test_steps_are_not_verified(self):
        variants = [[], self.jobs+self.jobs[:1], [self.jobs[0] | {"conclusion": "failure"}, self.jobs[1]],
                    [self.jobs[0] | {"steps": []}, self.jobs[1]]]
        for jobs in variants:
            with self.subTest(jobs=jobs):
                original = self.jobs; self.jobs = jobs
                self.assertEqual(self.collect()["status"], "unavailable")
                self.jobs = original

    def test_expired_missing_old_attempt_or_oversized_artifacts_are_not_verified(self):
        for change in ({"expired": True}, {"name": "linux-amd64-release-1"},
                       {"size_in_bytes": evidence.MAX_ARCHIVE+1}):
            original = deepcopy(self.artifacts)
            self.artifacts[0].update(change)
            self.assertEqual(self.collect()["status"], "unavailable")
            self.artifacts = original
        self.artifacts += self.artifacts[:1]
        self.assertEqual(self.collect()["status"], "unavailable")

    def test_junit_failures_errors_empty_counts_and_entities_fail_closed(self):
        for xml in (b'<testsuite tests="0"/>', b'<testsuite tests="2"><testcase/></testsuite>',
                    b'<testsuite tests="1" failures="1"><testcase/></testsuite>',
                    b'<testsuite tests="1"><testcase><failure/></testcase></testsuite>',
                    b'<testsuite tests="1" errors="1"><testcase/></testsuite>',
                    b'<testsuite tests="1" skipped="1"><testcase/></testsuite>',
                    b'<!DOCTYPE testsuite><testsuite tests="0"/>'):
            with self.subTest(xml=xml), self.assertRaises(ValueError):
                evidence.parse_junit(xml)

    def test_archive_reads_fixed_files_without_extracting_other_members(self):
        reports = evidence.read_reports(archive())
        self.assertEqual(set(reports), {"environment.log", "ctest-results.xml"})
        buffer = io.BytesIO()
        with zipfile.ZipFile(buffer, "w") as zipped:
            zipped.writestr("environment.log", "a" * 65537)
            zipped.writestr("ctest-results.xml", "<testsuite/>")
        with self.assertRaises(ValueError): evidence.read_reports(buffer.getvalue())

    def test_missing_sanitizer_controls_prevent_verified_status(self):
        with tempfile.TemporaryDirectory() as directory, patch.object(evidence, "api", side_effect=self.api), \
             patch.object(evidence, "download", return_value=archive(sanitizer=False)):
            result = evidence.collect(CONTEXT, Path(directory), 0)
        self.assertEqual(result["status"], "unavailable")

    def test_rerun_during_collection_invalidates_evidence(self):
        def changed(endpoint):
            result = self.api(endpoint)
            return result | {"run_attempt": 3} if endpoint.endswith("/runs/20") else result
        with tempfile.TemporaryDirectory() as directory, patch.object(evidence, "api", side_effect=changed), \
             patch.object(evidence, "download", return_value=archive()):
            self.assertEqual(evidence.collect(CONTEXT, Path(directory), 0)["status"], "unavailable")

    def test_pr_moving_after_download_invalidates_evidence(self):
        calls = 0
        def moved(endpoint):
            nonlocal calls
            result = self.api(endpoint)
            if endpoint.endswith("/pulls/12"):
                calls += 1
                if calls > 1:
                    result = deepcopy(result)
                    result["base"]["sha"] = "f"*40
            return result
        with tempfile.TemporaryDirectory() as directory, patch.object(evidence, "api", side_effect=moved), \
             patch.object(evidence, "download", return_value=archive()):
            self.assertEqual(evidence.collect(CONTEXT, Path(directory), 0)["status"], "unavailable")

    def test_api_failure_never_becomes_verified(self):
        with tempfile.TemporaryDirectory() as directory, patch.object(evidence, "api", side_effect=KeyError("unavailable")):
            self.assertEqual(evidence.collect(CONTEXT, Path(directory), 0)["status"], "unavailable")

    def test_workflow_keeps_analysis_guards_and_installs_dependencies(self):
        text = (Path(__file__).parents[1]/"workflows/opencode-review.yml").read_text()
        self.assertIn("runs-on: ubuntu-22.04", text)
        self.assertIn("actions: read", text)
        self.assertIn("libyaml-cpp-dev libssl-dev", text)
        self.assertIn("OPENCODE_CI_EVIDENCE:", text)
        self.assertIn('bash "$RUNNER_TEMP/prepare_review_toolchain.sh"', text)
        self.assertLess(text.index("Verify review publication guards"), text.index("Checkout frozen review head"))
        self.assertIn("never approve by default or because checks passed", text)
        self.assertIn("No full duplicate sanitizer build", text)
        self.assertLess(text.index("Verify review toolchain outside checkout"), text.index("Review with DeepSeek and GLM"))


if __name__ == "__main__":
    unittest.main()
