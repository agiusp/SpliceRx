from app import palette


def test_categorical_fixed_order_and_other():
    cmap = palette.categorical_colors(["b", "a", "b", "c"] + [f"x{i}" for i in range(10)])
    assert cmap["b"] == palette.CATEGORICAL[0][0]
    assert cmap["a"] == palette.CATEGORICAL[1][0]
    assert cmap["c"] == palette.CATEGORICAL[2][0]
    assert cmap["x7"] == palette.OTHER[0]  # 9th distinct onward


def test_sequential_monotone_endpoints():
    assert palette.sequential_color(0) == palette.SEQUENTIAL[0]
    assert palette.sequential_color(1) == palette.SEQUENTIAL[-1]


def test_categorical_offset_gives_non_overlapping_hues():
    a = palette.categorical_colors(["LUAD", "LUSC"], offset=0)
    b = palette.categorical_colors(["CR", "PR"], offset=2)
    assert set(a.values()).isdisjoint(b.values())


def test_natural_key_orders_numeric_strings():
    assert sorted(["10", "2", "1"], key=palette.natural_key) == ["1", "2", "10"]


def test_sequential_ramps_are_distinct():
    ends = [r[-1] for r in palette.SEQUENTIAL_RAMPS]
    assert len(set(ends)) == len(ends)
    assert palette.sequential_color(1, ramp=1) == palette.SEQUENTIAL_RAMPS[1][-1]


def test_bivariate_corners():
    assert palette.bivariate_color(0, 0).lower() == palette.BIV_LO_LO.lower()
    assert palette.bivariate_color(1, 0).lower() == palette.BIV_HI_LO.lower()
    assert palette.bivariate_color(0, 1).lower() == palette.BIV_LO_HI.lower()
    assert palette.bivariate_color(1, 1).lower() == palette.BIV_HI_HI.lower()


def test_shape_set():
    smap = palette.shape_for(["p", "q", "r"])
    assert smap["p"] == "circle" and smap["q"] == "square"
    assert len(palette.SHAPES) == 8
