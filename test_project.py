from project import *
import pytest
import numpy as np


def main():
    test_get_random_matrix()
    test_check_matrix()
    test_compute_statistics()


def test_get_random_matrix():
    n = np.random.randint(0, 100, 1)[0]
    k = np.random.randint(0, n, 1)[0]
    k = int(k)
    n = int(n)
    mat = get_random_matrix(n, k)
    assert mat.shape[0] == mat.shape[1] == n
    assert all(mat[i, :].sum() == k + 1 for i in range(n))
    assert all(mat[:, j].sum() == k + 1 for j in range(n))

    with pytest.raises(ValueError):
        get_random_matrix(n, n + 1)
        get_random_matrix(n, "k")


def test_check_matrix():
    mat = [
        [1, 0, 0, 0, 1, 0],
        [0, 1, 0, 0, 0, 1],
        [1, 0, 1, 0, 0, 0],
        [0, 1, 0, 1, 0, 0],
        [0, 0, 1, 0, 1, 0],
        [0, 0, 0, 1, 0, 1],
    ]
    with pytest.raises(TypeError):
        check_matrix(mat, 1)
    mat = np.array(mat)
    assert check_matrix(mat, 1) == True
    assert check_matrix(mat, 0) == False
    assert check_matrix(mat, 5) == False
    with pytest.raises(TypeError):
        check_matrix(mat, "1")
    mat = [
        [1, 0, 0, 0, 1, 0],
        [0, 1, 0, 0, 0, 1],
        [1, 0, 1, 0, 0, 0],
    ]
    mat = np.array(mat)
    assert check_matrix(mat, 1) == False


def test_compute_statistics():
    fitness = [
        [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7],
        [0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65],
        [0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75],
    ]
    fitness = np.array(fitness)
    stats = compute_statistics(fitness, 0.99)
    test_values = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7]
    test_values = np.array(test_values)
    length = len(test_values)
    assert type(stats) == dict
    assert len(stats["mean"][:]) == length
    for i in range(length):
        assert round(stats["mean"][i], 5) == round(test_values[i], 5)
    with pytest.raises(AttributeError):
        compute_statistics("fitness", 0.99)


if __name__ == "__main__":
    main()
