import numpy as np


class NK:
    def __init__(self, dependencymap: np.array):
        self.n = self.get_n(dependencymap)
        self.k = self.get_k(dependencymap)

    @property
    def n(self):
        return self._n

    @n.setter
    def n(self, var: int):
        self._n = var

    @property
    def k(self):
        return self._k

    @k.setter
    def k(self, var: int):
        self._k = var

    @staticmethod
    def get_n(dependencymap: np.array) -> int:
        """
        Returns the number of genes in the genome
        :param dependencymap: interdependence matrix
        :return: number of genes in genome
        """
        return dependencymap.shape[0]

    @staticmethod
    def get_k(dependencymap: np.array) -> int:
        """
        Returns the number of interdependencies based on the dependency matrix
        :param dependencymap: interdependence matrix
        :return: number of interdependencies per gene
        """
        return dependencymap[0].sum() - 1

    def convert_number_to_bin_nparray(self, var: int):
        # convert int into binary bitstring
        rd: str = bin(var)[2:].zfill(self.n)
        # unpack, convert to int and into numpy array
        rd_unpacked: list = [*rd]
        for i in range(len(rd_unpacked)):
            rd_unpacked[i] = int(rd_unpacked[i])
        return np.array(rd_unpacked)


class Landscape(NK):
    def __init__(self, dependencymap: np.array, alleles: int = 2):
        if self.check_dependencymap(dependencymap):
            super().__init__(dependencymap)
            # Options per gene, default setting is 2, i.e. binary problem
            self.alleles = alleles
            # interdependence matrix
            self.dependencymap = dependencymap
            # random fitness contributions
            self.contributions = np.random.uniform(
                0, 1, self.get_required_fitness_contributions()
            )
            # lookup table to compute indices for fitness contributions
            self.lookuptable = np.zeros([self.n, self.k + 1])
            for i in range(self.n):
                k = 0
                for j in range(self.n):
                    if self.dependencymap[i, j] == 1:
                        self._lookuptable[i, k] = j
                        k = k + 1
            # store the position of the global maximum in the landscape
            self.global_max_position, self.global_max_fitness = self.get_global_max()

        else:
            raise ValueError("Invalid interaction matrix")

    def __str__(self) -> str:
        return (
            f"NK landscape with N = {self.n} Genes and K = {self.k} interdependencies."
        )

    @property
    def global_max_position(self):
        return self._global_max_position

    @global_max_position.setter
    def global_max_position(self, var: np.ndarray):
        self._global_max_position = var

    @property
    def contributions(self):
        return self._contributions

    @contributions.setter
    def contributions(self, var: np.array):
        self._contributions = var

    @property
    def lookuptable(self):
        return self._lookuptable

    @lookuptable.setter
    def lookuptable(self, var: np.ndarray):
        self._lookuptable = var

    @property
    def alleles(self):
        return self._alleles

    @alleles.setter
    def alleles(self, var: int):
        self._alleles = var

    @property
    def dependencymap(self):
        return self._dependencymap

    @dependencymap.setter
    def dependencymap(self, dependencymap: np.ndarray):
        self._dependencymap = dependencymap

    @classmethod
    def get(cls, dependencymap: np.ndarray, alleles: int = 2) -> object:
        """
        Initializes the NK landscape
        :param dependencymap: square interaction matrix of datatype boolean
        :param alleles: options per gene, optimal parameter, set to two if no input (results in binary options for every gene)
        :return: NK landscape object
        """
        return cls(dependencymap, alleles)

    @staticmethod
    def check_dependencymap(dependencymap: np.ndarray) -> bool:
        """
        Checks the interaction matrix for correctness. The method checks for
        1. Two dimensions
        2. Square matrix
        3. Ones along the main diagonal
        4. Symmetric patterns, i.e., is the sum in all columns and all rows the same
        :param dependencymap: interaction pattern in the form of a numpy array
        :return: True if all conditions are met, False otherwise
        """
        # checks
        is_two_dimensional: bool = len(dependencymap.shape) == 2
        is_square: bool = dependencymap.shape[0] == dependencymap.shape[1]
        is_one_along_diagonal: bool = all(
            values == 1 for values in np.diagonal(dependencymap)
        )
        # check if sum of columns and rows in the array are all the same
        is_symmetric: bool = all(
            dependencymap[i, :].sum() == dependencymap[0, :].sum()
            and dependencymap[:, i].sum() == dependencymap[0, :].sum()
            for i in range(dependencymap.shape[0])
        )
        return (
            is_two_dimensional and is_square and is_one_along_diagonal and is_symmetric
        )

    def get_required_fitness_contributions(self) -> int:
        """
        Returns the number of fitness contributions required to compute the NK landscape
        :return: number of required fitness contributions
        """
        return self.n * self.alleles ** (self.k + 1)

    def check_genome(self, genome: np.ndarray) -> None:
        """
        Checks whether a genome is valid, i.e, of the right length (n) and has only binary values
        :rtype: None
        :param genome: numpy array containing 0 and 1 of length n
        :return: raises ValueError if genome is invalid
        """
        is_length_n = len(genome) == self.n
        is_binary = all(values == 1 or values == 0 for values in genome)
        if is_length_n and is_binary:
            pass
        else:
            raise ValueError("Invalid genome")

    def get_fitness_gene(
        self, idx: int, genome: np.ndarray, normalized: bool = True
    ) -> float:
        """
        Return the contribution of one gene to the fitness of the genome
        :param idx: index of the gene for which the fitness contribution should be returned
        :param genome: configuration of the genome
        :param normalized: indicates whether the absolute or the normalized fitness is returned by the function
        :return: fitness contribution
        """
        self.check_genome(genome)
        tmp_idx: int = 0
        for i in range(self.k + 1):
            if genome[int(self.lookuptable[idx, i])] == 1:
                tmp_idx = 2 * tmp_idx + 1
            else:
                tmp_idx = 2 * tmp_idx
        contribution_idx = self.n * tmp_idx + idx
        if normalized:
            return self.contributions[contribution_idx] / self.global_max_fitness
        else:
            return self.contributions[contribution_idx]

    def get_fitness_genome(self, genome: np.ndarray, normalized: bool = True) -> float:
        """
        Returns the fitness of a genome
        :param genome: genome for which the fitness should be computed
        :param normalized: indicates whether the function returns the absolute of the normalized fitness of the genome
        :return: fitness of the genome
        """
        self.check_genome(genome)
        perf: float = 0
        for i in range(self.n):
            perf += self.get_fitness_gene(i, genome, normalized)
        return perf / self.n

    def get_global_max(self):
        tmp_pos = []
        tmp_fitness: float = 0.0
        for i in range((2 ** self.n) - 1):
            if (
                self.get_fitness_genome(
                    self.convert_number_to_bin_nparray(i), normalized=False
                )
                > tmp_fitness
            ):
                tmp_fitness = self.get_fitness_genome(
                    self.convert_number_to_bin_nparray(i), normalized=False
                )
                tmp_pos = self.convert_number_to_bin_nparray(i)
        return tmp_pos, tmp_fitness


#


class Agent(NK):
    def __init__(
        self,
        dependencymap: np.ndarray,
        error_mean: float = 0.0,
        error_std: float = 0.0,
    ):
        super().__init__(dependencymap)
        self.error_mean = error_mean
        self.error_std = error_std
        self.position = self.get_random_position()

    def __str__(self):
        return f"Agent, operates on landscape with with N = {self.n} genes and K = {self.k} interdependencies, makes errors with mean {self.error_mean} and std {self.error_std}"

    @property
    def error_mean(self):
        return self._error_mean

    @error_mean.setter
    def error_mean(self, var: float):
        self._error_mean = var

    @property
    def error_std(self):
        return self._error_std

    @error_std.setter
    def error_std(self, var: float):
        self._error_std = var

    @property
    def position(self):
        return self._position

    @position.setter
    def position(self, var: np.ndarray):
        self._position = var

    @classmethod
    def get(
        cls, dependencymap: np.ndarray, error_mean: float = 0.0, error_std: float = 0.0,
    ) -> object:
        return cls(dependencymap, error_mean, error_std)

    @staticmethod
    def flip_bit(var: int) -> int:
        if var == 0 or var == 1:
            if var == 0:
                return 1
            else:
                return 0
        else:
            raise ValueError("Bit must be one or zero")

    def get_random_position(self):
        # get random number
        rd = np.random.randint(0, 2 ** self.n - 1)
        return self.convert_number_to_bin_nparray(rd)

    def get_alternative_position(self, discoveries: int = 1):
        # get random position of a bit to flip
        rd = np.random.choice(
            np.arange(0, self.n, dtype=int), size=discoveries, replace=True
        )
        options = np.empty((discoveries, self.n))
        for i in range(len(rd)):
            options[i, :] = self.position
            options[i, rd[i]] = self.flip_bit(options[i, rd[i]])
        return options

    def evolve(self, landscape: object):
        options = self.get_alternative_position(1)
        if landscape.get_fitness_genome(self.position) < landscape.get_fitness_genome(
            options[0]
        ) + np.random.normal(self.error_mean, self.error_std):
            self.update_position(options[0])

    def update_position(self, var: np.ndarray):
        self._position = var


def main():
    print("This is an extended version of Kauffman's NK model!")


if __name__ == "__main__":
    main()
