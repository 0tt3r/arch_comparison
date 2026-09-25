"""Pinnacle architecture ledger (Webster et al., arXiv:2602.11457v2; refs/2602.11457v2.pdf).

Generalized-bicycle qLDPC processing blocks and magic engines, with the paper's logical-error ansatz.
Values are the paper's (Table I, Table II "Log. Meas." row, Eqs. 4-13, 18-20) and the magic-engine
operating points vendored in CUDA-Q Logical (cudaq.logical.architectures.pinnacle, which cites the same
paper). The paper assumes uniform circuit-level depolarizing noise and most-likely-error decoding.

A sampler is one processing unit: ceil(N_logical/k) processing blocks plus one magic engine (Eq. 18).
It performs one logical measurement per logical cycle (d + 2 code cycles), which covers each T-state
injection, so T gates are sequential within a unit.
"""
from dataclasses import dataclass
from math import ceil


@dataclass(frozen=True)
class Code:
    name: str
    n: int
    k: int
    d: int
    block_qubits: int      # processing block: code block + 4 gadgets + 4 bridges (Table I)

    @property
    def cycle_rounds(self):
        return self.d + 2


CODES = [Code('gb30', 30, 8, 4, 140), Code('gb62', 62, 10, 6, 244), Code('gb126', 126, 12, 10, 452),
         Code('gb254', 254, 14, 16, 860), Code('gb510', 510, 16, 24, 1620)]

# Table II, logical-measurement fit: p_L per logical qubit per logical cycle = (A/k) (p/B)^(d/2 + C).
FIT = dict(A=6.2, B=0.0158, C=0.47)


@dataclass(frozen=True)
class Engine:
    name: str
    p: float
    output_error: float
    qubits: int
    cycles_per_attempt: int
    acceptance: float


# CUDA-Q Logical _MAGIC_ENGINES (output_infidelity, physical_qubits, cycles_per_attempt, acceptance).
ENGINES = [Engine('p1e-4_1e-9', 1e-4, 9.6e-10, 592, 14, .9981),
           Engine('p1e-4_1e-11', 1e-4, 8.4e-12, 1807, 18, .98),
           Engine('p1e-3_1e-9', 1e-3, 6.2e-10, 4410, 26, .90),
           Engine('p1e-3_1e-11', 1e-3, 1.2e-11, 5430, 26, .90)]


def logical_error(code, p):
    """Per logical qubit per logical cycle (Eq. 13)."""
    return FIT['A']/code.k*(p/FIT['B'])**(code.d/2 + FIT['C'])


def unit_qubits(code, n_logical, engine):
    """Eq. 18: one processing unit plus its magic engine."""
    return code.block_qubits*ceil(n_logical/code.k) + engine.qubits


def table_iv():
    """Replay of the paper's Fermi-Hubbard Table IV: N = 2L^2 + 2 logical qubits, d = 24 (p = 1e-3) or 10 (1e-4),
    p_out = 1e-9 engines. Eq. 20 prints n_me = 1807 (the 1e-11 engine), but Table IV matches the
    stated p_out = 1e-9 engine at p = 1e-4 (592 qubits), which we use."""
    rows = []
    codes = {c.d: c for c in CODES}
    for L in range(8, 33, 2):
        n1 = codes[24].block_qubits*ceil((L*L+1)/8) + 4410
        n2 = codes[10].block_qubits*ceil((L*L+1)/6) + 592
        rows.append(dict(L=L, qubits_p1e3=n1, qubits_p1e4=n2))
    return rows


if __name__ == '__main__':
    published = {8: (19e3, 5.6e3), 16: (58e3, 20e3), 24: (123e3, 44e3), 32: (213e3, 78e3)}
    for r in table_iv():
        if r['L'] in published:
            a, b = published[r['L']]
            print(r['L'], r['qubits_p1e3'], f'(paper {a:.0f})', r['qubits_p1e4'], f'(paper {b:.0f})')
    for c in CODES:
        print(c.name, [f'{logical_error(c, p):.1g}' for p in (1e-3, 1e-4)])
