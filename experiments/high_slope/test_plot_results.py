"""Check that lightweight figure rendering does not bypass raw-data validation."""
import contextlib
import io
import json
from pathlib import Path
import unittest
from unittest.mock import patch

import plot_results


class PlotCommandTests(unittest.TestCase):
    def test_summary_mode_only_renders_existing_statistics(self):
        summary = {"rows": [{"N": 10}], "two_point_slopes": [], "schema_version": 1}
        with patch.object(Path, "read_text", return_value=json.dumps(summary)), \
                patch.object(Path, "mkdir") as mkdir, \
                patch.object(plot_results, "plot_summary") as render, \
                patch.object(plot_results, "write_outputs") as raw, \
                contextlib.redirect_stdout(io.StringIO()) as output:
            plot_results.main(["--summary", "summary.json", "--output", "figures"])
        mkdir.assert_called_once_with(parents=True, exist_ok=True)
        render.assert_called_once_with(summary, Path("figures"))
        raw.assert_not_called()
        self.assertIn("not revalidated", output.getvalue())

    def test_raw_mode_keeps_validating_write_outputs_path(self):
        records = [{"id": "w7_baseline"}]
        with patch.object(Path, "read_text", return_value=json.dumps(records)), \
                patch.object(plot_results, "write_outputs",
                             return_value={"rows": [], "two_point_slopes": []}) as raw, \
                patch.object(plot_results, "plot_summary") as render, \
                contextlib.redirect_stdout(io.StringIO()):
            plot_results.main(["--output", "results", "--distance-records", "distance.json"])
        raw.assert_called_once_with(Path("results"), records)
        render.assert_not_called()

    def test_summary_rejects_distance_record_override(self):
        with contextlib.redirect_stderr(io.StringIO()), self.assertRaises(SystemExit) as error:
            plot_results.main(["--summary", "summary.json", "--distance-records", "distance.json"])
        self.assertEqual(error.exception.code, 2)


if __name__ == "__main__":
    unittest.main()
