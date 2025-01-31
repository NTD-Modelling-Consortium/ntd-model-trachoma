import pytest

from ntdmc_trachoma.trachoma_functions import readPlatformData

EXPECTED = [
    [2020.0, 0, 100, 0.8, 0, 2],
    [2021.0, 0, 100, 0.8, 0, 2],
    [2022.0, 0, 100, 0.8, 0, 2],
    [2023.0, 0, 100, 0.8, 0, 2],
    [2024.0, 0, 100, 0.85, 0, 2],
    [2024.083, 1, 10, 0.85, 1, 2],
    [2024.17, 1, 10, 0.85, 1, 2],
    [2025.0, 0, 100, 0.85, 0, 2],
]

def test_read_platform_data():
    result = readPlatformData("test_coverage_data.csv", "MDA")
    assert result == EXPECTED

def test_read_platform_data_other_data_path():
    result = readPlatformData(
        "test_cov_data.csv",
        "MDA",
        data_path="data"
    )
    assert result == EXPECTED

def test_read_platform_data_wrong_data_path():
    with pytest.raises(FileNotFoundError):
        result = readPlatformData(
            "test_coverage_data.csv",
            "MDA",
            data_path="/does/not/exist"
        )
