"""Check analysis credential isolation and publication prerequisites."""

import os
import subprocess
import unittest
from unittest.mock import patch

from run_opencode_review import run_review


class AnalysisTests(unittest.TestCase):
    @patch("run_opencode_review.subprocess.check_output", side_effect=["abc\n", ""])
    @patch("run_opencode_review.subprocess.run")
    def test_analysis_has_provider_key_but_no_github_credentials(self, run, output):
        credentials = {"GH_TOKEN": "test-gh", "GITHUB_TOKEN": "test-github",
                       "GH_ENTERPRISE_TOKEN": "test-enterprise",
                       "GITHUB_ENTERPRISE_TOKEN": "test-enterprise-github",
                       "OPENROUTER_API_KEY": "test-provider"}
        with patch.dict(os.environ, credentials):
            run_review("opencode", "Review the frozen diff", "abc")
            environment = run.call_args.kwargs["env"]
            self.assertEqual(environment["OPENROUTER_API_KEY"], "test-provider")
            for key in credentials:
                if key != "OPENROUTER_API_KEY":
                    self.assertNotIn(key, environment)
                    self.assertEqual(os.environ[key], credentials[key])
        self.assertEqual(run.call_args.args[0][1], "run")
        self.assertNotIn("github", run.call_args.args[0])

    @patch("run_opencode_review.subprocess.check_output")
    @patch("run_opencode_review.subprocess.run", side_effect=subprocess.CalledProcessError(1, "opencode"))
    def test_failed_analysis_stops_before_publication_prerequisites(self, run, output):
        with self.assertRaises(subprocess.CalledProcessError):
            run_review("opencode", "Review", "abc")
        output.assert_not_called()

    @patch("run_opencode_review.subprocess.check_output", return_value="different\n")
    @patch("run_opencode_review.subprocess.run")
    def test_branch_switch_cannot_reach_publication(self, run, output):
        with self.assertRaisesRegex(ValueError, "frozen checkout"):
            run_review("opencode", "Review", "abc")

    @patch("run_opencode_review.subprocess.check_output", side_effect=["abc\n", " M source.cpp\n"])
    @patch("run_opencode_review.subprocess.run")
    def test_dirty_checkout_cannot_reach_publication(self, run, output):
        with self.assertRaisesRegex(ValueError, "checkout"):
            run_review("opencode", "Review", "abc")


if __name__ == "__main__":
    unittest.main()
