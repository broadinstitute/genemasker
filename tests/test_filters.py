import pandas as pd
import pytest

from genemasker import filters


def test_curated_mask_sets_cover_registered_baseline_masks():
    registered = {func.__name__ for func in filters.MASKS}

    assert len(filters.LOW_FREQUENCY_MASKS) == 9
    assert len(filters.RARE_MASKS) == 9
    assert len(filters.MASK_SETS["all"]) == 18
    assert len(set(filters.MASK_SETS["all"])) == 18
    assert set(filters.MASK_SETS["all"]) == registered


def test_removed_legacy_masks_are_not_registered():
    registered = {func.__name__ for func in filters.MASKS}

    assert "new_damaging_ic25" not in registered
    assert "new_damaging_og25" not in registered
    assert "new_damaging_og50" not in registered
    assert "x23633568_m1" not in registered
    assert "x30828346_m1" not in registered
    assert "x32141622_m4" not in registered
    assert "x32141622_m7" not in registered


def test_resolve_masks_preserves_requested_order_and_rejects_unknown_names():
    requested = ["LoF_HC_0_001", "x37348876_m8"]

    assert [func.__name__ for func in filters.resolve_masks(requested)] == requested
    with pytest.raises(ValueError, match="Unknown mask name"):
        filters.resolve_masks(["not_a_mask"])


def test_lof_masks_apply_their_frequency_thresholds():
    df = pd.DataFrame(
        {
            "LoF": ["HC", "HC", "LC"],
            "MAF": [0.00005, 0.005, 0.00005],
        }
    )

    assert filters.LoF_HC_0_01(df).tolist() == [True, True, False]
    assert filters.LoF_HC_0_001(df).tolist() == [True, False, False]
