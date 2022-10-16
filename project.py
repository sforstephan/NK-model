from NK import Landscape, Agent
from matrices import matrices
import argparse
import numpy as np
from numpy.random import default_rng
import scipy.stats as st
import matplotlib.pyplot as plt
import json


def main():

    parser = argparse.ArgumentParser(description="Get parameters for NK model")
    parser.add_argument(
        "-n",
        type=int,
        default=5,
        help="an integer for the number of genes in the genome",
    )
    parser.add_argument(
        "-k",
        type=int,
        default=2,
        help="an integer for the number of dependencies between genes",
    )
    parser.add_argument(
        "-matrix",
        type=str,
        default="random",
        help="a name of a predefined matrix (overrules n and k and takes predefined matrix from matrices.py)",
    )
    parser.add_argument(
        "-time",
        type=int,
        default=100,
        help="an integer for length of the observation period",
    )
    parser.add_argument(
        "-repeat",
        type=int,
        default=8,
        help="an integer for the number of repetitions of the process of evolution",
    )

    parser.add_argument(
        "-mean",
        type=float,
        default=0,
        help="mean of the error when predicting fitness",
    )

    parser.add_argument(
        "-std",
        type=float,
        default=0,
        help="standard deviation of the error when predicting fitness",
    )
    parser.add_argument(
        "-confidence",
        type=float,
        default=0.9,
        help="confidence interval for plot",
    )

    args = parser.parse_args()
    args_dict = vars(args)

    if args.matrix != "random" and args.matrix in matrices.keys():
        pattern = matrices[args.matrix]
    elif args.matrix == "random":
        pattern = get_random_matrix(args.n, args.k)
    else:
        raise ValueError("Invalid parameter for matrix")

    fitness = np.zeros((args.repeat, args.time), dtype=float)

    for simround in range(args.repeat):
        landscape = Landscape(pattern)
        agent = Agent(pattern, args.mean, args.std)
        for timestep in range(args.time):
            if timestep == 0:
                fitness[simround, timestep] = landscape.get_fitness_genome(
                    agent.position
                )
            if timestep > 0:
                agent.evolve(landscape)
                fitness[simround, timestep] = landscape.get_fitness_genome(
                    agent.position
                )

    stats = compute_statistics(fitness, args.confidence)

    # create and save figure
    x = np.arange(0, args.time)
    plt.plot(x, stats["mean"], label="Mean normalized fitness")
    plt.fill_between(
        x,
        stats["lower_conf"],
        stats["upper_conf"],
        alpha=0.2,
        label=f"Confidence interval, \u03B1 = {round(1-args.confidence,2)}",
    )
    plt.ylim(0.5, 1.1)
    plt.xlabel("Time")
    plt.ylabel("Normalized fitness")
    plt.title("Fitness")
    plt.grid(alpha=0.2)
    plt.legend(loc="lower right")
    plt.savefig("fitness.jpg", dpi=300, transparent=True)

    # format data for json and save
    with open("parameters.json", "w") as file:
        json.dump(args_dict, file)

    for key in stats.keys():
        stats[key] = stats[key][:].tolist()
    with open("statistics.json", "w") as file:
        json.dump(stats, file)

    fitness_dict = {}
    for i in range(args.repeat):
        fitness_dict[f"Round {i}"] = fitness[i, :].tolist()
    with open("fitness.json", "w") as file:
        json.dump(fitness_dict, file)


def get_random_matrix(n: int, k: int):
    if n > k:
        matrix = np.identity(n, dtype=int)
        matrix = np.array(matrix)

        randon_number_generator = default_rng()
        idxs: list = []
        for i in range(n):
            for j in range(n):
                while True:
                    idxs = randon_number_generator.choice(n, size=k, replace=False)
                    if i not in idxs:
                        break
            for idx in idxs:
                matrix[i, idx] = 1

        while not check_matrix(matrix, k):
            missing_ones = []
            excess_ones = []
            for j in range(n):
                if matrix[:, j].sum() < k + 1:
                    missing_ones.append(j)
                elif matrix[:, j].sum() > k + 1:
                    excess_ones.append(j)
            excess_col = randon_number_generator.choice(excess_ones, size=1)
            missing_col = randon_number_generator.choice(missing_ones, size=1)
            while True:
                row = np.random.randint(0, n)
                if (
                    matrix[row, missing_col] == 0
                    and matrix[row, excess_col] == 1
                    and row != excess_col
                    and row != missing_col
                ):
                    matrix[row, missing_col] = 1
                    matrix[row, excess_col] = 0
                    break
        return np.array(matrix)
    else:
        raise ValueError("Check values for parameters N and K")


def check_matrix(var: np.ndarray, k: int):
    if type(var) is np.ndarray and type(k) is int:
        if (
            var.shape[0] == var.shape[1]
            and all(var[i, :].sum() == k + 1 for i in range(var.shape[0]))
            and all(var[:, j].sum() == k + 1 for j in range(var.shape[0]))
        ):
            return True
        else:
            return False
    else:
        raise TypeError


def compute_statistics(fitness: np.ndarray, confidence):
    stats = {
        "mean": np.zeros(fitness.shape[1], dtype=float),
        "upper_conf": np.zeros(fitness.shape[1], dtype=float),
        "lower_conf": np.zeros(fitness.shape[1], dtype=float),
    }

    for i in range(fitness.shape[1]):
        if fitness.shape[0] <= 30:
            tmp = st.t.interval(
                confidence=confidence,
                df=len(fitness[:, i]) - 1,
                loc=np.mean(fitness[:, i]),
                scale=st.sem(fitness[:, i]),
            )
        else:
            tmp = st.norm.interval(
                confidence=confidence,
                loc=np.mean(fitness[:, i]),
                scale=st.sem(fitness[:, i]),
            )
        stats["lower_conf"][i] = tmp[0]
        stats["upper_conf"][i] = tmp[1]
        stats["mean"][i] = fitness[:, i].mean()

    return stats


if __name__ == "__main__":
    main()
