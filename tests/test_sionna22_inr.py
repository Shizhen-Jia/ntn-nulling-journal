"""Regression checks for the Sionna 2.2 aggregate INR extension.

Run in sionna22:
    MPLBACKEND=Agg python -m unittest discover -s tests -p test_sionna22_inr.py -v
"""

from contextlib import redirect_stdout
import importlib.metadata
import io
import unittest

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import torch
import sionna.rt as rt
import mitsuba as mi
import drjit as dr


def make_scene():
    scene = rt.load_scene()
    scene.tx_array = rt.PlanarArray(num_rows=1, num_cols=1,
                                   pattern="iso", polarization="V")
    scene.rx_array = rt.PlanarArray(num_rows=1, num_cols=1,
                                   pattern="iso", polarization="V")
    scene.add(rt.Transmitter("tx0", position=[-2, 0, 5], power_dbm=30))
    scene.add(rt.Transmitter("tx1", position=[2, 0, 5],
                             power_dbm=30 + 10 * np.log10(2)))
    scene.add(rt.Receiver("rx0", position=[0, 0, 1]))
    return scene


class AggregateINRTests(unittest.TestCase):
    def setUp(self):
        self.scene = make_scene()
        self.rm = rt.PlanarRadioMap(self.scene, cell_size=(1, 1),
                                    center=(0, 0, 0), orientation=(0, 0, 0),
                                    size=(2, 2))
        self.rm._thermal_noise_power = mi.Float(1e-6)
        self.rm._pathgain_map = mi.TensorXf(np.array([
            [[1, 2], [0, 3]], [[2, 1], [0, 4]],
        ], dtype=np.float32) * 1e-6)
        # 1 W * gain_tx0 + 2 W * gain_tx1, divided by 1 uW noise.
        self.expected = np.array([[5, 4], [0, 11]], dtype=np.float32)

    def tearDown(self):
        plt.close("all")

    def test_formula_shape_and_alias(self):
        self.assertEqual(self.rm.inr_ntn.shape, (1, 2, 2))
        np.testing.assert_allclose(self.rm.inr_ntn.numpy()[0], self.expected,
                                   rtol=1e-6)
        np.testing.assert_array_equal(self.rm.inr.numpy(), self.rm.inr_ntn.numpy())
        for metric in ("inr", "inr_ntn"):
            for tx in (None, 0, 1, "tx0", "tx1"):
                np.testing.assert_allclose(
                    self.rm.transmitter_radio_map(metric, tx).numpy(),
                    self.expected, rtol=1e-6)
            for tx in (-1, 2, "missing"):
                with self.assertRaises(ValueError):
                    self.rm.transmitter_radio_map(metric, tx)

    def test_show_and_cdf_use_db_not_dbm(self):
        for metric in ("inr", "inr_ntn"):
            fig = self.rm.show(metric=metric)
            plotted = fig.axes[0].images[0].get_array()
            mask = self.expected > 0
            np.testing.assert_allclose(plotted[mask],
                                       10 * np.log10(self.expected[mask]), rtol=1e-6)
            self.assertIn("INR", fig.axes[0].get_title())
            self.assertIn("[dB]", fig.axes[1].get_ylabel())
            fig, x, cdf = self.rm.cdf(metric=metric, bins=8)
            self.assertIn("INR", fig.axes[0].get_xlabel())
            self.assertAlmostEqual(float(cdf.numpy()[-1]), 1.0)
            self.assertAlmostEqual(float(x.numpy()[-1]),
                                   10 * np.log10(11), places=4)
            self.assertTrue(np.all(np.diff(cdf.numpy()) >= 0))

    def test_constant_and_empty_maps(self):
        self.rm._pathgain_map = dr.full(mi.TensorXf, 1e-6, (2, 2, 2))
        _, x, cdf = self.rm.cdf(metric="inr")
        np.testing.assert_allclose(x.numpy(), 10 * np.log10(3), rtol=1e-6)
        np.testing.assert_array_equal(cdf.numpy(), np.ones(200))
        self.rm._pathgain_map = dr.zeros(mi.TensorXf, (2, 2, 2))
        np.testing.assert_array_equal(self.rm.inr.numpy(), np.zeros((1, 2, 2)))
        with self.assertRaisesRegex(ValueError, "no coverage"):
            self.rm.cdf(metric="inr")

    def test_sampling_aggregate_metric(self):
        # Only the cell with INR=11 falls in [10, 11] dB.
        for metric in ("inr", "inr_ntn"):
            cells, count = self.rm.sample_cells(
                8, metric=metric, min_val_db=10, max_val_db=11,
                tx_association=False)
            np.testing.assert_array_equal(count.numpy(), [8, 8])
            np.testing.assert_array_equal(cells.numpy(), np.full((2, 8), 3))
            with self.assertRaisesRegex(ValueError, "tx_association=False"):
                self.rm.sample_cells(8, metric=metric)
            with self.assertRaisesRegex(ValueError, "no transmitter association"):
                self.rm.tx_association(metric)
        # Preserve upstream sampling and association for ordinary metrics.
        association = self.rm.tx_association("rss").numpy()
        self.assertEqual(association.reshape(-1)[0], 1)
        self.assertEqual(association.reshape(-1)[2], -1)

    def test_mesh_map(self):
        mesh = mi.Mesh("test_mesh", vertex_count=4, face_count=2)
        params = mi.traverse(mesh)
        params["vertex_positions"] = mi.Float([-1, -1, 0, 1, -1, 0,
                                               -1, 1, 0, 1, 1, 0])
        params["faces"] = mi.UInt([0, 1, 2, 1, 3, 2])
        params.update()
        rm = rt.MeshRadioMap(self.scene, mesh)
        rm._thermal_noise_power = mi.Float(1e-6)
        rm._pathgain_map = mi.TensorXf(np.array([[1, 2], [2, 1]],
                                               dtype=np.float32) * 1e-6)
        self.assertEqual(rm.inr.shape, (1, 2))
        np.testing.assert_allclose(rm.inr.numpy(), [[5, 4]], rtol=1e-6)
        self.assertAlmostEqual(float(rm.cdf(metric="inr")[2].numpy()[-1]), 1.0)

    def test_scene_render_and_preview(self):
        camera = rt.Camera(position=[0, -6, 6], look_at=[0, 0, 0])
        fig = self.scene.render(camera=camera, radio_map=self.rm,
                                rm_metric="inr_ntn", resolution=(32, 32),
                                num_samples=1)
        self.assertIsNotNone(fig)
        # preview displays its cached widget and intentionally returns None.
        with redirect_stdout(io.StringIO()):
            self.scene.preview(radio_map=self.rm, rm_metric="inr_ntn")
        self.assertIsNotNone(self.scene._preview_widget)


class GPUEnvironmentTests(unittest.TestCase):
    def test_pytorch_cuda(self):
        self.assertEqual(importlib.metadata.version("sionna"), "2.2.0")
        self.assertEqual(importlib.metadata.version("sionna-rt"), "2.2.0")
        self.assertTrue(torch.cuda.is_available())
        a = torch.arange(16, device="cuda", dtype=torch.float32).reshape(4, 4)
        np.testing.assert_allclose((a @ a.T).cpu().numpy(),
                                   np.arange(16).reshape(4, 4) @ np.arange(16).reshape(4, 4).T)

    def test_gpu_paths_and_radio_map(self):
        self.assertTrue(mi.variant().startswith("cuda"), mi.variant())
        scene = make_scene()
        paths = rt.PathSolver()(scene=scene, max_depth=1,
                                 specular_reflection=False, refraction=False)
        a, _ = paths.cir(out_type="numpy")
        self.assertTrue(np.isfinite(a).all())
        self.assertTrue(np.any(np.abs(a) > 0))
        rm = rt.RadioMapSolver()(scene=scene, center=(0, 0, 0),
                                  orientation=(0, 0, 0), size=(8, 8),
                                  cell_size=(2, 2), samples_per_tx=8192,
                                  max_depth=1, specular_reflection=False,
                                  refraction=False)
        inr = rm.inr_ntn.numpy()[0]
        self.assertTrue(np.isfinite(inr).all())
        self.assertTrue(np.any(inr > 0))
        expected = rm.rss.numpy().sum(axis=0) / float(rm._thermal_noise_power[0])
        np.testing.assert_allclose(inr, expected, rtol=1e-5)


if __name__ == "__main__":
    unittest.main(verbosity=2)
