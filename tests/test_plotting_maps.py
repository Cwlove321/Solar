from types import SimpleNamespace

import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import numpy as np

from solarutils.plotting.maps import mask_weak_hmi_field, plot_hmi_on_aia


def test_mask_weak_hmi_field_sets_weak_pixels_to_nan():
    data = np.array([[-40.0, -20.0], [29.0, 50.0]])

    masked = mask_weak_hmi_field(data, threshold=30.0)

    assert np.isfinite(masked[0, 0])
    assert np.isnan(masked[0, 1])
    assert np.isnan(masked[1, 0])
    assert np.isfinite(masked[1, 1])
    assert np.array_equal(data, np.array([[-40.0, -20.0], [29.0, 50.0]]))


def test_plot_hmi_on_aia_overlays_masked_hmi_data():
    aia_map = SimpleNamespace(
        data=np.arange(4, dtype=float).reshape(2, 2),
        plot_settings={"norm": Normalize()},
        cmap="gray",
        wavelength=SimpleNamespace(value=193),
        date=SimpleNamespace(value="2024-10-08T23:00:04.843"),
    )
    hmi_map = SimpleNamespace(data=np.array([[-50.0, -10.0], [20.0, 80.0]]))

    artists = plot_hmi_on_aia(
        aia_map,
        hmi_map,
        extent=(100.0, 200.0, 300.0, 400.0),
        hmi_threshold=30.0,
        colorbar=False,
        boxes=[(120.0, 150.0, 320.0, 350.0)],
    )

    hmi_array = np.asarray(artists["hmi"].get_array())
    assert np.isfinite(hmi_array[0, 0])
    assert np.isnan(hmi_array[0, 1])
    assert np.isnan(hmi_array[1, 0])
    assert np.isfinite(hmi_array[1, 1])
    assert len(artists["boxes"]) == 1
    assert "HMI (|B| > 30 G) overlaid on AIA 193" in artists["axis"].get_title()
    plt.close(artists["figure"])
