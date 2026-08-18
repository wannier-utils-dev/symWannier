import numpy as np
import pytest

from symwannier.mmn import Mmn


class DummyNnkp:
    nk = 2
    nb = 2
    bvec = np.array([[1.0, 0.0, 0.0], [-1.0, 0.0, 0.0]])

    @staticmethod
    def calc_bvec(info):
        return np.asarray(info[2:5], dtype=float)


def write_mmn(path, blocks):
    lines = ["synthetic full-k MMN", "1 2 2"]
    for ik, ikb, gx, value in blocks:
        lines.append(f"{ik} {ikb} {gx} 0 0")
        lines.append(f"{value.real} {value.imag}")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def test_full_k_mmn_blocks_are_reordered_from_headers(tmp_path):
    mmn_file = tmp_path / "permuted.mmn"
    write_mmn(
        mmn_file,
        [
            (1, 1, -1, 20.0 + 2.0j),
            (1, 2, 1, 10.0 + 1.0j),
            (2, 2, -1, 40.0 + 4.0j),
            (2, 1, 1, 30.0 + 3.0j),
        ],
    )

    mmn = Mmn(str(mmn_file), nnkp=DummyNnkp())

    expected = np.array(
        [
            [10.0 + 1.0j, 20.0 + 2.0j],
            [30.0 + 3.0j, 40.0 + 4.0j],
        ]
    )
    assert np.allclose(mmn.mmn[:, :, 0, 0], expected)
    assert np.array_equal(mmn.kb2k, [[1, 0], [0, 1]])
    assert np.array_equal(mmn.kpb_info[:, :, 2], [[1, -1], [1, -1]])


def test_full_k_mmn_rejects_duplicate_b_vectors(tmp_path):
    mmn_file = tmp_path / "duplicate.mmn"
    write_mmn(
        mmn_file,
        [
            (1, 2, 1, 10.0 + 1.0j),
            (1, 2, 1, 20.0 + 2.0j),
            (2, 1, 1, 30.0 + 3.0j),
            (2, 2, -1, 40.0 + 4.0j),
        ],
    )

    with pytest.raises(ValueError, match="duplicate MMN b-vector 1 at k-point 1"):
        Mmn(str(mmn_file), nnkp=DummyNnkp())
