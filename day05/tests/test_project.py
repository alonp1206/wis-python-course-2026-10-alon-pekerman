import csv
import io
import os
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import numpy as np
import pandas as pd


PROJECT_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_DIR))
os.environ.setdefault("MAVEN_SECRET_KEY", "maven-test-secret")

import app as maven_app
from pipeline import validate_count_matrix, validate_input_tables


class MavenProjectTests(unittest.TestCase):
    def setUp(self):
        self.temp_dir = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp_dir.cleanup)
        root = Path(self.temp_dir.name)
        self.upload_root = root / "uploads"
        self.static_root = root / "static"
        self.upload_root.mkdir()
        self.static_root.mkdir()

        for name, value in (
            ("UPLOAD_ROOT", str(self.upload_root)),
            ("STATIC_ROOT", str(self.static_root)),
        ):
            patcher = patch.object(maven_app, name, value)
            patcher.start()
            self.addCleanup(patcher.stop)

        self.original_static_folder = maven_app.app.static_folder
        self.addCleanup(
            setattr, maven_app.app, "static_folder", self.original_static_folder
        )
        maven_app.app.static_folder = str(self.static_root)
        maven_app.app.config.update(TESTING=True)
        self.client = maven_app.app.test_client()

    @staticmethod
    def _csv_upload(frame, index_label):
        buffer = io.StringIO()
        frame.to_csv(buffer, index_label=index_label)
        return io.BytesIO(buffer.getvalue().encode("utf-8")), "input.csv"

    def _synthetic_inputs(self):
        rng = np.random.default_rng(2026)
        sample_ids = ["C1", "C2", "C3", "T1", "T2", "T3"]
        rows = []
        for gene_number in range(120):
            if gene_number < 30:
                counts = np.concatenate(
                    (rng.poisson(25, size=3), rng.poisson(450, size=3))
                )
            else:
                baseline = int(rng.integers(30, 250))
                counts = rng.poisson(baseline, size=len(sample_ids))
            rows.append([f"gene_{gene_number:03}", *counts.tolist()])

        count_table = pd.DataFrame(
            rows, columns=["gene_id", *sample_ids]
        ).set_index("gene_id")
        duplicate = count_table.loc[["gene_000"]].copy()
        duplicate.iloc[:, :] = 1
        count_table = pd.concat([count_table, duplicate])
        metadata = pd.DataFrame(
            {
                "Condition": ["Control"] * 3 + ["Treated"] * 3,
                "Lineage": ["A", "A", "B", "A", "B", "B"],
            },
            index=sample_ids,
        )
        return count_table, metadata

    def _upload(self, count_table, metadata):
        counts_file = self._csv_upload(count_table, "gene_id")
        metadata_file = self._csv_upload(metadata, "sample_id")
        return self.client.post(
            "/upload",
            data={
                "counts_file": counts_file,
                "meta_file": metadata_file,
            },
            content_type="multipart/form-data",
            follow_redirects=True,
        )

    def test_complete_analysis_workflow(self):
        count_table, metadata = self._synthetic_inputs()
        response = self._upload(count_table, metadata)
        self.assertEqual(response.status_code, 200)
        self.assertIn(b"Merged 1 duplicate gene row", response.data)
        self.assertIn(b"Step 2: Explore Low-Count Filtering", response.data)

        preview = self.client.get("/filter_preview?min_count=10&min_samples=2")
        self.assertEqual(preview.status_code, 200)
        preview_data = preview.get_json()
        self.assertAlmostEqual(sum(item["y"] for item in preview_data["histogram"]), 100)
        self.assertGreaterEqual(preview_data["selected_percent"], 0)
        self.assertLessEqual(preview_data["selected_percent"], 100)

        response = self.client.post(
            "/filter",
            data={"min_count": "1", "min_samples": "2"},
            follow_redirects=True,
        )
        self.assertEqual(response.status_code, 200)
        self.assertIn(b"Analysis Dashboard", response.data)
        with self.client.session_transaction() as session:
            run_id = session["run_id"]

        filtered_counts = pd.read_csv(
            self.upload_root / run_id / "filtered_counts.csv", index_col=0
        )
        self.assertEqual(len(filtered_counts), 120)

        response = self.client.get("/run_pca?color_by=Condition")
        self.assertEqual(response.status_code, 200)
        self.assertIn(b"PCA Results", response.data)
        self.assertTrue((self.static_root / run_id / "pca_plot.png").is_file())

        response = self.client.post(
            "/run_deseq",
            data={
                "design_col": "Condition",
                "target_group": "Treated",
                "reference_group": "Control",
                "comparison_mode": "pairwise",
                "lfc_threshold": "1",
            },
        )
        self.assertEqual(response.status_code, 200)
        self.assertIn(b"Differential Expression (DESeq2)", response.data)
        self.assertTrue((self.static_root / run_id / "volcano_plot.png").is_file())
        results = pd.read_csv(
            self.upload_root / run_id / "deseq_results.csv", index_col="gene"
        )
        self.assertGreaterEqual(int(results["is_significant"].sum()), 2)
        download = self.client.get("/download/deseq")
        self.assertEqual(download.status_code, 200)
        download.close()

        response = self.client.post("/run_clustering", data={"k_clusters": "2"})
        self.assertEqual(response.status_code, 200)
        self.assertIn(b"K-Means Clustering", response.data)
        self.assertTrue(
            (self.static_root / run_id / "kmeans_heatmap.png").is_file()
        )
        assignments = pd.read_csv(
            self.upload_root / run_id / "cluster_assignments.csv"
        )
        self.assertGreaterEqual(assignments["cluster"].nunique(), 2)
        download = self.client.get("/download/clusters")
        self.assertEqual(download.status_code, 200)
        download.close()

        with patch.object(maven_app, "run_enrichment", return_value=(None, [])) as enrich:
            response = self.client.post(
                "/run_enrichment",
                data={"database": "KEGG_2021_Human", "direction": "up"},
            )
        self.assertEqual(response.status_code, 200)
        self.assertIn(b"Pathway Enrichment Analysis", response.data)
        enrich.assert_called_once()

    def test_input_validation_and_filter_guardrails(self):
        duplicate_counts = pd.DataFrame(
            [[2, 4], [3, 5]], index=["geneA", "geneA"], columns=["S1", "S2"]
        )
        merged = validate_count_matrix(duplicate_counts)
        self.assertEqual(merged.loc["geneA"].tolist(), [5, 9])

        with self.assertRaisesRegex(ValueError, "non-negative integer"):
            validate_count_matrix(
                pd.DataFrame([[1.5, 2]], index=["geneA"], columns=["S1", "S2"])
            )
        with self.assertRaisesRegex(ValueError, "non-negative"):
            validate_count_matrix(
                pd.DataFrame([[-1, 2]], index=["geneA"], columns=["S1", "S2"])
            )

        metadata = pd.DataFrame({"Condition": ["A", "B"]}, index=["S1", "S1"])
        with self.assertRaisesRegex(ValueError, "metadata must be unique"):
            validate_input_tables(
                pd.DataFrame([[1, 2]], index=["geneA"], columns=["S1", "S2"]),
                metadata,
            )

        invalid_counts = pd.DataFrame(
            [[1.5, 2]], index=["geneA"], columns=["S1", "S2"]
        )
        valid_metadata = pd.DataFrame({"Condition": ["A", "B"]}, index=["S1", "S2"])
        response = self._upload(invalid_counts, valid_metadata)
        self.assertEqual(response.status_code, 200)
        self.assertIn(b"non-negative integer raw counts", response.data)
        with self.client.session_transaction() as session:
            self.assertNotIn("run_id", session)

        response = self.client.get("/run_pca", follow_redirects=True)
        self.assertEqual(response.status_code, 200)
        self.assertIn(b"Upload RNA-seq Data", response.data)


if __name__ == "__main__":
    unittest.main()
