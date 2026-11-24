"""Run regression tests with outgoing socket connections blocked."""

from pathlib import Path
import socket
import unittest
from unittest.mock import patch


def main():
    root = Path(__file__).resolve().parent.parent
    suite = unittest.defaultTestLoader.discover(str(root / "tests"), top_level_dir=str(root.parents[1]))
    with patch.object(socket.socket, "connect", side_effect=AssertionError("Live network access in offline test")), \
            patch.object(socket.socket, "connect_ex", side_effect=AssertionError("Live network access in offline test")):
        result = unittest.TextTestRunner(verbosity=2).run(suite)
    return 0 if result.wasSuccessful() else 1


if __name__ == "__main__":
    raise SystemExit(main())
