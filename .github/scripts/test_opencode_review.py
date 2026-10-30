"""Safety checks for formal approval publication; no live GitHub writes."""

from copy import deepcopy
import json
from pathlib import Path
import tempfile
import subprocess
import unittest
from unittest.mock import patch

import opencode_review as review


class ReviewTests(unittest.TestCase):
    def setUp(self):
        self.context = {"repository": "owner/repo", "number": 6,
                        "head_sha": "a" * 40, "base_sha": "b" * 40}
        self.pr = {"state": "open", "draft": False,
                   "head": {"sha": "a" * 40, "repo": {"full_name": "owner/repo"}},
                   "base": {"sha": "b" * 40}}
        self.verdict = {"head_sha": "a" * 40, "base_sha": "b" * 40,
                        "review_complete": True, "verdict": "APPROVE",
                        "blocking_findings": [], "summary": "Reviewed the diff; no blocking findings.",
                        "architecture_review": {"required": False, "status": "not_required",
                                                "summary": "Localized documentation change."}}

    def test_completed_review_approves_exact_commit(self):
        payload = review.review_payload(self.context, self.verdict, self.pr)
        self.assertEqual(payload["event"], "APPROVE")
        self.assertEqual(payload["commit_id"], "a" * 40)
        self.assertIn(review.MARKER, payload["body"])

    def test_incomplete_or_blocking_review_cannot_approve(self):
        for changes in ({"review_complete": False}, {"review_complete": "true"},
                        {"blocking_findings": ["Confirmed data corruption"]},
                        {"summary": ""}, {"verdict": "approve"}):
            with self.subTest(changes=changes), self.assertRaises(ValueError):
                review.review_payload(self.context, self.verdict | changes, self.pr)

    def test_stale_or_wrong_revision_cannot_approve(self):
        for part in ("head", "base"):
            changed = deepcopy(self.pr)
            changed[part]["sha"] = "c" * 40
            with self.subTest(part=part), self.assertRaises(ValueError):
                review.review_payload(self.context, self.verdict, changed)
        with self.assertRaises(ValueError):
            review.review_payload(self.context, self.verdict | {"head_sha": "c" * 40}, self.pr)

    def test_draft_closed_or_fork_pr_cannot_approve(self):
        for changes in ({"draft": True}, {"state": "closed"},
                        {"head": {"sha": "a" * 40, "repo": None}},
                        {"head": {"sha": "a" * 40, "repo": {"full_name": "fork/repo"}}}):
            with self.subTest(changes=changes), self.assertRaises(ValueError):
                review.review_payload(self.context, self.verdict, self.pr | changes)

    def test_blockers_request_changes_and_uncertainty_comments(self):
        blocked = self.verdict | {"verdict": "REQUEST_CHANGES", "blocking_findings": ["Crash on empty input"]}
        self.assertEqual(review.review_payload(self.context, blocked, self.pr)["event"], "REQUEST_CHANGES")
        incomplete = self.verdict | {"verdict": "COMMENT", "review_complete": False}
        self.assertEqual(review.review_payload(self.context, incomplete, self.pr)["event"], "COMMENT")
        with self.assertRaises(ValueError):
            review.review_payload(self.context, blocked | {"verdict": "COMMENT"}, self.pr)

    def test_missing_or_unknown_fields_fail_closed(self):
        missing = self.verdict.copy()
        del missing["review_complete"]
        for verdict in (None, [], missing, self.verdict | {"unexpected": True}):
            with self.subTest(verdict=verdict), self.assertRaises(ValueError):
                review.review_payload(self.context, verdict, self.pr)

    def test_required_consultation_must_complete_before_approval(self):
        for status in ("failed", "inconclusive"):
            verdict = self.verdict | {"architecture_review": {
                "required": True, "status": status, "summary": "GLM consultation did not complete."}}
            with self.subTest(status=status), self.assertRaises(ValueError):
                review.review_payload(self.context, verdict, self.pr)
            verdict.update(verdict="COMMENT", review_complete=False)
            self.assertEqual(review.review_payload(self.context, verdict, self.pr)["event"], "COMMENT")

    def test_completed_consultation_is_reported_in_formal_review(self):
        verdict = self.verdict | {"architecture_review": {
            "required": True, "status": "completed", "summary": "Task ses_example: contracts are consistent."}}
        payload = review.review_payload(self.context, verdict, self.pr)
        self.assertEqual(payload["event"], "APPROVE")
        self.assertIn("## AI review\n", payload["body"])
        self.assertIn("DeepSeek V4.1 Flash", payload["body"])
        self.assertIn("GLM-5.3: completed", payload["body"])

    def test_invalid_or_missing_consultation_evidence_cannot_approve(self):
        invalid = [None, {}, {"required": "true", "status": "completed", "summary": "Checked."},
                   {"required": True, "status": "not_required", "summary": "Skipped."},
                   {"required": False, "status": "completed", "summary": "Checked."},
                   {"required": True, "status": "completed", "summary": ""}]
        for architecture in invalid:
            with self.subTest(architecture=architecture), self.assertRaises(ValueError):
                review.review_payload(self.context, self.verdict | {"architecture_review": architecture}, self.pr)

    @patch.object(review, "api")
    def test_prepare_only_reads_metadata_and_preserves_existing_reviews(self, api):
        api.return_value = self.pr
        with tempfile.TemporaryDirectory() as directory:
            context = Path(directory) / "context.json"
            output = Path(directory) / "output"
            review.prepare("owner/repo", 6, context, output)
            self.assertEqual(json.loads(context.read_text()), self.context)
            self.assertIn("allowed=true\n", output.read_text())
            self.assertIn(f"head_sha={self.context['head_sha']}\n", output.read_text())
        api.assert_called_once_with("repos/owner/repo/pulls/6")

    @patch.object(review, "api")
    def test_push_racing_publication_withdraws_approval(self, api):
        changed = deepcopy(self.pr)
        changed["head"]["sha"] = "c" * 40
        api.side_effect = [self.pr, {"id": 99, "html_url": "https://example/review"}, changed, {}]
        with tempfile.TemporaryDirectory() as directory:
            verdict = Path(directory) / "verdict.json"
            verdict.write_text(json.dumps(self.verdict))
            with self.assertRaises(ValueError):
                review.publish(self.context, verdict)
        self.assertEqual(api.call_args_list[-1].args[:2],
                         ("repos/owner/repo/pulls/6/reviews/99/dismissals", "PUT"))

    @patch.object(review, "api")
    def test_missing_result_never_calls_github(self, api):
        with tempfile.TemporaryDirectory() as directory:
            with self.assertRaises(FileNotFoundError):
                review.publish(self.context, Path(directory) / "missing.json")
        api.assert_not_called()

    @patch.object(review, "api")
    def test_failed_final_revision_check_withdraws_approval(self, api):
        api.side_effect = [self.pr, {"id": 99, "html_url": "https://example/review"},
                           RuntimeError("GitHub unavailable"), {}]
        with tempfile.TemporaryDirectory() as directory:
            verdict = Path(directory) / "verdict.json"
            verdict.write_text(json.dumps(self.verdict))
            with self.assertRaises(RuntimeError):
                review.publish(self.context, verdict)
        self.assertEqual(api.call_args_list[-1].args[:2],
                         ("repos/owner/repo/pulls/6/reviews/99/dismissals", "PUT"))

    @patch.object(review, "api")
    def test_fork_prepare_stops_before_context_creation(self, api):
        fork = deepcopy(self.pr)
        fork["head"]["repo"]["full_name"] = "fork/repo"
        api.return_value = fork
        with tempfile.TemporaryDirectory() as directory:
            context = Path(directory) / "context.json"
            output = Path(directory) / "output"
            review.prepare("owner/repo", 6, context, output)
            self.assertFalse(context.exists())
            self.assertEqual(output.read_text(), "allowed=false\n")
        self.assertEqual(api.call_count, 1)

    @patch.object(review, "api")
    def test_protected_branch_dismissal_denied_supersedes_approval(self, api):
        api.side_effect = [subprocess.CalledProcessError(
            1, ["gh", "api"], stderr="gh: Resource not accessible (HTTP 403)"), {}]
        review.dismiss("repos/owner/repo/pulls/6", 99, "Re-review required", "a" * 40)
        endpoint, method, payload = api.call_args.args
        self.assertEqual((endpoint, method), ("repos/owner/repo/pulls/6/reviews", "POST"))
        self.assertEqual(payload["event"], "REQUEST_CHANGES")
        self.assertEqual(payload["commit_id"], "a" * 40)


if __name__ == "__main__":
    unittest.main()
