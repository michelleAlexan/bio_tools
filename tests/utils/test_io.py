import pytest
from pathlib import Path
from bio_tools.utils.io import load_yaml

def test_load_yaml():
    """
    test the utils.io.load_yaml function
    """
    # test non existing file
    non_existing_yaml_path = Path(__file__).parents[1] / "data/filedoesnotexist.yaml"
    with pytest.raises(FileNotFoundError):
        load_yaml(non_existing_yaml_path)

    # test existing file
    yaml_path = Path(__file__).parents[1] / "data/example.yaml"
    result = load_yaml(yaml_path)
    assert isinstance(result, dict)
    assert result["Oryza sativ"] == "Oryza sativa"

