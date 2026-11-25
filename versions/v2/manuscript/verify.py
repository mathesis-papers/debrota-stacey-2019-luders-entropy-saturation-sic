import itertools
import json
import sys

import sympy as symbolic


def normalize(matrix):
    return matrix.applyfunc(symbolic.simplify)


def require_equal(actual, expected, label):
    difference = actual - expected
    zero = normalize(difference) if isinstance(difference, symbolic.MatrixBase) else symbolic.simplify(difference)
    if zero != (symbolic.zeros(*zero.shape) if isinstance(zero, symbolic.MatrixBase) else 0):
        raise AssertionError(label + ': ' + str(zero))


def examples(dimension):
    identity = symbolic.eye(dimension)
    if dimension == 1:
        return [identity]
    if dimension == 2:
        pauli = [symbolic.Matrix([[0, 1], [1, 0]]),
                 symbolic.Matrix([[0, -symbolic.I], [symbolic.I, 0]]),
                 symbolic.diag(1, -1)]
        signs = [(1, 1, 1), (1, -1, -1), (-1, 1, -1), (-1, -1, 1)]
        return [normalize((identity + sum((sign * matrix / symbolic.sqrt(3)
                 for sign, matrix in zip(row, pauli)), symbolic.zeros(2))) / 2)
                for row in signs]
    omega = (-1 + symbolic.I * symbolic.sqrt(3)) / 2
    projectors = []
    for exponent in range(3):
        power = symbolic.simplify(omega ** exponent)
        for entries in [(0, 1, -power), (-power, 0, 1), (1, -power, 0)]:
            vector = symbolic.Matrix(entries)
            projectors.append(normalize(vector * vector.conjugate().T / 2))
    return projectors


def check_example(dimension):
    identity = symbolic.eye(dimension)
    projectors = examples(dimension)
    effects = [projector / dimension for projector in projectors]
    require_equal(sum(effects, symbolic.zeros(dimension)), identity, 'POVM normalization')
    for projector in projectors:
        require_equal(projector.conjugate().T, projector, 'Hermitian projector')
        require_equal(projector * projector, projector, 'idempotent projector')
        require_equal(symbolic.trace(projector), 1, 'rank-one trace')
    gram = symbolic.zeros(dimension ** 2)
    for row, first in enumerate(effects):
        for column, second in enumerate(effects):
            gram[row, column] = symbolic.simplify(symbolic.trace(first * second))
            target = symbolic.Rational(dimension * (row == column) + 1,
                                       dimension ** 2 * (dimension + 1))
            require_equal(gram[row, column], target, 'SIC overlap')
    if gram.det() == 0:
        raise AssertionError('MIC linear independence')
    for row in range(dimension):
        for column in range(dimension):
            matrix = symbolic.zeros(dimension)
            matrix[row, column] = 1
            channel = sum((dimension * effect * matrix * effect for effect in effects),
                          symbolic.zeros(dimension))
            require_equal(channel, (matrix + symbolic.trace(matrix) * identity) / (dimension + 1),
                          'depolarizing action on matrix units')
    swap = symbolic.zeros(dimension ** 2)
    for row in range(dimension):
        for column in range(dimension):
            swap[row * dimension + column, column * dimension + row] = 1
    moment = sum((symbolic.kronecker_product(projector, projector) for projector in projectors),
                 symbolic.zeros(dimension ** 2)) / dimension ** 2
    require_equal(moment, (symbolic.eye(dimension ** 2) + swap) / (dimension * (dimension + 1)),
                  'finite projective two-design')
    entropy = symbolic.log(dimension + 1) - symbolic.Rational(2, dimension + 1) * symbolic.log(2)
    return {'dimension': dimension, 'projectors': len(projectors),
            'gram_entries': dimension ** 4, 'matrix_units': dimension ** 2,
            'two_design': True, 'entropy_expression': str(entropy),
            'arithmetic': 'exact in Q(i,sqrt(3))'}


def check_spectral_gaps():
    cases = 0
    for dimension in range(1, 7):
        for spectrum in itertools.product(range(3), repeat=dimension):
            gap = sum(spectrum) ** 2 - sum(value ** 2 for value in spectrum)
            pair_sum = 2 * sum(spectrum[first] * spectrum[second]
                               for first in range(dimension) for second in range(first + 1, dimension))
            if gap != pair_sum or gap < 0 or ((gap == 0) != (sum(value != 0 for value in spectrum) <= 1)):
                raise AssertionError('positive spectral rank defect')
            cases += 1
    return {'dimensions': [1, 6], 'eigenvalue_alphabet': [0, 1, 2], 'spectra': cases,
            'rank_at_most_one_iff_zero_gap': True, 'arithmetic': 'integer'}


def check_strict_example():
    dimension = 2
    identity = symbolic.eye(dimension)
    projectors = examples(dimension)
    beta = symbolic.Rational(1, 2)
    effects = [((1 - beta) * identity + 2 * beta * projector) / 4 for projector in projectors]
    coefficient = (symbolic.sqrt(1 + beta) + symbolic.sqrt(1 - beta)) / 4
    deviation = (symbolic.sqrt(1 + beta) - symbolic.sqrt(1 - beta)) / 4
    roots = [coefficient * identity + deviation * (2 * projector - identity) for projector in projectors]
    alpha = (1 + symbolic.sqrt(3)) / 3
    for effect, root in zip(effects, roots):
        require_equal(root * root, effect, 'positive square root of full-rank effect')
    gram = symbolic.Matrix([[symbolic.simplify(symbolic.trace(first * second))
                            for second in effects] for first in effects])
    if gram.det() == 0:
        raise AssertionError('full-rank MIC independence')
    for row in range(2):
        for column in range(2):
            matrix = symbolic.zeros(2)
            matrix[row, column] = 1
            channel = sum((root * matrix * root for root in roots), symbolic.zeros(2))
            require_equal(channel, alpha * matrix + (1 - alpha) * symbolic.trace(matrix) * identity / 2,
                          'strict-example channel')
    rank_sum = symbolic.simplify(sum(symbolic.trace(root) ** 2 for root in roots))
    require_equal(rank_sum, 2 + symbolic.sqrt(3), 'strict-example rank sum')
    return {'dimension': 2, 'beta': '1/2', 'effect_eigenvalues': ['1/8', '3/8'],
            'depolarizing_parameter': str(alpha), 'T_minus_d': str(rank_sum - 2),
            'entropy_deficit_lower_bound': 'sqrt(3)*log(2)/6 > 0', 'arithmetic': 'exact radicals'}


def check_overcomplete_design():
    identity = symbolic.eye(2)
    pauli = [symbolic.Matrix([[0, 1], [1, 0]]),
             symbolic.Matrix([[0, -symbolic.I], [symbolic.I, 0]]), symbolic.diag(1, -1)]
    projectors = [(identity + sign * matrix) / 2 for matrix in pauli for sign in [-1, 1]]
    require_equal(sum(projectors, symbolic.zeros(2)) / 3, identity, 'six-effect normalization')
    swap = symbolic.Matrix([[1, 0, 0, 0], [0, 0, 1, 0],
                            [0, 1, 0, 0], [0, 0, 0, 1]])
    moment = sum((symbolic.kronecker_product(projector, projector) / 6
                  for projector in projectors), symbolic.zeros(4))
    require_equal(moment, (symbolic.eye(4) + swap) / 6, 'six-state weighted design')
    maximally_entangled = symbolic.Matrix([1, 0, 0, 1]) / symbolic.sqrt(2)
    entangled_projector = maximally_entangled * maximally_entangled.conjugate().T
    choi = symbolic.zeros(4)
    for row in range(2):
        for column in range(2):
            matrix = symbolic.zeros(2)
            matrix[row, column] = 1
            channel = sum((projector * matrix * projector / 3 for projector in projectors),
                          symbolic.zeros(2))
            require_equal(channel, (matrix + symbolic.trace(matrix) * identity) / 3,
                          'six-effect channel')
            choi += symbolic.kronecker_product(matrix, channel) / 2
    require_equal(choi, (symbolic.eye(4) + 2 * entangled_projector) / 6,
                  'entropy-exchange target state')
    require_equal(symbolic.trace(entangled_projector * choi), symbolic.Rational(1, 2),
                  'entropy-exchange target overlap')
    return {'dimension': 2, 'effects': 6, 'matrix_units': 4, 'weighted_two_design': True,
            'choi_identity': True, 'choi_eigenvalues': ['1/2', '1/6', '1/6', '1/6'],
            'arithmetic': 'exact in Q(i,sqrt(2))'}


def check_spectral_average_example():
    identity = symbolic.eye(2)
    pauli = [symbolic.Matrix([[0, 1], [1, 0]]),
             symbolic.Matrix([[0, -symbolic.I], [symbolic.I, 0]]), symbolic.diag(1, -1)]
    signs = [(1, 1, 1), (1, -1, -1), (-1, 1, -1), (-1, -1, 1)]
    axes = [symbolic.sqrt(2) / 2, symbolic.Rational(1, 2), symbolic.Rational(1, 2)]
    vectors = [symbolic.Matrix([sign * axis for sign, axis in zip(row, axes)]) for row in signs]
    projectors = [(identity + sum((entry * matrix for entry, matrix in zip(vector, pauli)),
                                  symbolic.zeros(2))) / 2 for vector in vectors]
    effects = [projector / 2 for projector in projectors]
    require_equal(sum(effects, symbolic.zeros(2)), identity, 'anisotropic MIC normalization')
    for projector in projectors:
        require_equal(projector * projector, projector, 'anisotropic rank-one effect')
        require_equal(projector.conjugate().T, projector, 'anisotropic Hermitian effect')
        require_equal(symbolic.trace(projector), 1, 'anisotropic projector trace')
    gram = symbolic.Matrix([[symbolic.trace(first * second) for second in effects] for first in effects])
    if gram.det() == 0:
        raise AssertionError('anisotropic MIC independence')
    bloch = sum((vector * vector.T / 4 for vector in vectors), symbolic.zeros(3))
    require_equal(bloch, symbolic.diag(symbolic.Rational(1, 2), symbolic.Rational(1, 4),
                                      symbolic.Rational(1, 4)), 'anisotropic Bloch matrix')
    for index, matrix in enumerate(pauli):
        channel = sum((projector * matrix * projector / 2 for projector in projectors),
                      symbolic.zeros(2))
        require_equal(channel, bloch[index, index] * matrix, 'anisotropic Pauli action')
    for row in signs:
        vector = symbolic.Matrix(row) / symbolic.sqrt(3)
        require_equal((vector.T * bloch.T * bloch * vector)[0], symbolic.Rational(1, 8),
                      'tetrahedral squared output radius')
    require_equal(symbolic.trace(bloch.T * bloch) / 3, symbolic.Rational(1, 8),
                  'Haar mean squared output radius')
    if bloch[0, 0] ** 2 == bloch[1, 1] ** 2:
        raise AssertionError('nonconstant Haar output radius')
    return {'dimension': 2, 'effects': 4, 'gram_determinant': str(symbolic.factor(gram.det())),
            'bloch_eigenvalues': ['1/2', '1/4', '1/4'],
            'tetrahedral_squared_output_radius': '1/8', 'haar_mean_squared_output_radius': '1/8',
            'haar_squared_radius_nonconstant': True,
            'spectral_mean_comparison': 'Strict inequality follows from the written strict-concavity argument.',
            'arithmetic': 'exact in Q(i,sqrt(2),sqrt(3))'}


def main():
    result = {'schema_version': 1, 'status': 'pass',
              'scope': 'Finite exact SIC identities in dimensions 1, 2, 3, spectral rank defects, a full-rank MIC, an overcomplete design and an anisotropic qubit MIC.',
              'software': {'python': '.'.join(map(str, sys.version_info[:3])), 'sympy': symbolic.__version__},
              'sic_examples': [check_example(dimension) for dimension in [1, 2, 3]],
              'spectral_gaps': check_spectral_gaps(), 'strict_example': check_strict_example(),
              'overcomplete_design': check_overcomplete_design(),
              'spectral_average_example': check_spectral_average_example()}
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
