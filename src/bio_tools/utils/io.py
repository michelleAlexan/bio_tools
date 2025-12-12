#%%
import yaml
from pathlib import Path


def load_yaml(path: str | Path) -> dict:
    """
    Load yaml file and return content as dictionary.
    """
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"The yaml file {p} does not exist.")
    
    with p.open("r") as y:
        yaml_dict = yaml.safe_load(y)
    return yaml_dict


# %%
