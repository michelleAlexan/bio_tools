from typing import Iterable, TypeVar
import warnings

T = TypeVar("T")
#%%

def internet_on(timeout: float = 3.0) -> bool:
    """Check internet by resolving a hostname."""
    import socket
    try:
        socket.setdefaulttimeout(timeout)
        host = socket.gethostbyname("www.google.com")
        s = socket.create_connection((host, 80), timeout)
        s.close()
        return True
    except OSError:
        return False


def ensure_list(x: Iterable[T]) -> list[T]:
    """
    If species is a single value and not an iterable, convert to list. 
    """
    if isinstance(x, str):
        return [x]

    if isinstance(x, Iterable):
        return list(x)

    raise TypeError("Species must be a string or an iterable of strings.")


def ensure_no_duplicates(l: list[T]) -> list[T]:
    """
    Take a list, and in case of duplicated values, remove them and raise warning. 
    """
    if not isinstance(l, list):
        raise TypeError(f"Passed argument is of type {type(l)}, but must be a list.")
    
    seen = set()
    result = []
    duplicated_values = set()

    for v in l:
        if v in seen:
            duplicated_values.add(v)
            continue
        seen.add(v)
        result.append(v)


    if duplicated_values:
        warnings.warn(
            f"Duplicate values {duplicated_values} detected and removed.",
            UserWarning
        )

    return result

