import pytest
import socket

from bio_tools.utils.constraint import ensure_no_duplicates, ensure_list, internet_on


@pytest.mark.parametrize(
    "input_list, expected, warning",
    [
        ([1, 2, 3], [1, 2, 3], None),
        ([1, 2, 2, 3, 1], [1, 2, 3], UserWarning),
        (["a", "b", "a", "c", "b"], ["a", "b", "c"], UserWarning),
        ([], [], None),
        ([42], [42], None),
    ],
)
def test_ensure_no_duplicates(input_list, expected, warning):
    if warning:
        with pytest.warns(warning):
            result = ensure_no_duplicates(input_list)
        assert result == expected

    else:
        result = ensure_no_duplicates(input_list)
        assert result == expected


@pytest.mark.parametrize(
    "input_value, expected",
    [
        ("Arabidopsis thaliana", ["Arabidopsis thaliana"]),
        (3702, [3702]),
        (["A", "B"], ["A", "B"]),
        (("A", "B"), ["A", "B"]),
        ({"A", "B"}, ["A", "B"]),
    ],
)
def test_ensure_list(input_value, expected):
    result = ensure_list(input_value)
    assert sorted(result) == sorted(expected)



@pytest.mark.parametrize(
    "gethostbyname_ok, create_connection_ok, expected",
    [
        (True, True, True),
        (False, True, False),
        (True, False, False),
    ],
)
def test_internet_on(monkeypatch, gethostbyname_ok, create_connection_ok, expected):
    
    def mock_gethostbyname(host):
        if not gethostbyname_ok:
            raise socket.gaierror
        return "127.0.0.1"

    def mock_create_connection(address, timeout):
        if not create_connection_ok:
            raise OSError
        class DummySocket:
            def close(self): ...
        return DummySocket()

    monkeypatch.setattr(socket, "gethostbyname", mock_gethostbyname)
    monkeypatch.setattr(socket, "create_connection", mock_create_connection)

    assert internet_on() is expected