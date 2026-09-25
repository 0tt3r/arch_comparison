"""Explain the physical circuit: Jordan-Wigner ordering, long Pauli strings, and the fSWAP network.

(a) The row snake of hubbard_circuits.snake_sites on a 4x4 lattice: horizontal bonds are adjacent
    in the Jordan-Wigner order; vertical bonds can be 2L-1 positions apart.
(b) Compiling one long hop directly: a CNOT parity ladder across every qubit in between.
(c) A fermionic-swap network on a two-row strip (3 sites per row): odd-even layers of fSWAPs
    reverse the mode order, every pair of modes is adjacent exactly once, and each vertical hop
    is applied (merged with its fSWAP) at that moment. All gates act on neighbouring qubits.
"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

from hubbard_circuits import snake_sites
from project_paths import ROOT

INK, MUTED, GRID = '#1f2933', '#6b7785', '#c9d1db'
SHORT, LONG, HOP = '#2a6fb0', '#c2552e', '#3a9a5b'


def snake_panel(ax, L=4):
    sites = snake_sites(L)
    index = {s: i for i, s in enumerate(sites)}
    for (x, y) in sites:
        for nb in ((x+1, y), (x, y+1)):
            if nb in index:
                ax.plot([x, nb[0]], [-y, -nb[1]], color=GRID, lw=1.2, zorder=1)
    xs, ys = zip(*sites)
    ax.plot(xs, [-y for y in ys], color=MUTED, lw=1.6, ls=(0, (1, 2)), zorder=2)
    short = ((1, 0), (2, 0))
    long_ = ((1, 0), (1, 1))
    for (a, b), color in ((short, SHORT), (long_, LONG)):
        ax.plot([a[0], b[0]], [-a[1], -b[1]], color=color, lw=4, zorder=3, solid_capstyle='round')
    for s, i in index.items():
        ax.add_patch(plt.Circle((s[0], -s[1]), .27, facecolor='white', edgecolor=INK, lw=1.2, zorder=4))
        ax.text(s[0], -s[1], str(i), ha='center', va='center', fontsize=10, zorder=5)
    ax.text(1.5, .62, 'horizontal bond: sites 1–2\nadjacent in JW order → $X_1X_2$', color=SHORT,
            ha='center', fontsize=8.5)
    ax.text(1.5, -L+.35, 'vertical bond: sites 1–6, five apart\n$X_1Z_2Z_3Z_4Z_5X_6$, drawn in (b)', color=LONG,
            ha='center', va='top', fontsize=8.5)
    ax.set_xlim(-.8, L-.2)
    ax.set_ylim(-L-.8, 1.3)
    ax.set_aspect('equal')
    ax.axis('off')
    ax.set_title('(a) Jordan–Wigner snake, 4×4', loc='left', fontsize=11)


def wires(ax, n, x0, x1, labels):
    for q in range(n):
        ax.plot([x0, x1], [-q, -q], color=INK, lw=.9, zorder=1)
        ax.text(x0-.15, -q, labels[q], ha='right', va='center', fontsize=9)


def cnot(ax, x, c, t):
    ax.plot([x, x], [-c, -t], color=INK, lw=1, zorder=2)
    ax.add_patch(plt.Circle((x, -c), .09, color=INK, zorder=3))
    ax.add_patch(plt.Circle((x, -t), .2, facecolor='white', edgecolor=INK, lw=1, zorder=3))
    ax.plot([x-.2, x+.2], [-t, -t], color=INK, lw=1, zorder=4)
    ax.plot([x, x], [-t-.2, -t+.2], color=INK, lw=1, zorder=4)


def box(ax, x, y0, y1, text, color, width=.6):
    ax.add_patch(FancyBboxPatch((x-width/2, -y1-.3), width, (y1-y0)+.6, boxstyle='round,pad=0.02',
                                facecolor=color, edgecolor=INK, lw=.8, alpha=.9, zorder=3))
    ax.text(x, -(y0+y1)/2, text, ha='center', va='center', fontsize=8, zorder=4, color='white' if color != 'white' else INK)


def direct_panel(ax):
    n = 6
    wires(ax, n, 0, 12.4, [f'$q_{{{i}}}$' for i in range(1, 7)])
    box(ax, .6, 0, 0, 'B', '#8894a3', .5)
    box(ax, .6, 5, 5, 'B', '#8894a3', .5)
    for k in range(n-1):
        cnot(ax, 1.5+k, k, k+1)
    box(ax, 6.7, 5, 5, '$R_Z$', LONG, .7)
    for k in reversed(range(n-1)):
        cnot(ax, 7.9+(n-2-k), k, k+1)
    box(ax, 12.0, 0, 0, 'B', '#8894a3', .5)
    box(ax, 12.0, 5, 5, 'B', '#8894a3', .5)
    ax.text(6.2, -6.1, 'one Pauli term of the long hop: parity ladder over every qubit between\n'
            '2(k−1) CNOTs in depth ∝ k, for modes k apart (B: basis change)', ha='center', fontsize=8.5, color=MUTED)
    ax.set_xlim(-1, 12.8)
    ax.set_ylim(-6.9, .8)
    ax.axis('off')
    ax.set_title('(b) Direct: $e^{-i\\theta X_1Z_2\\cdots Z_5X_6}$', loc='left', fontsize=11)


def network_panel(ax):
    # Two-row strip, 3 sites per row, JW snake: positions 0,1,2 = (0,0),(1,0),(2,0); 3,4,5 = (2,1),(1,1),(0,1).
    coords = [(0, 0), (1, 0), (2, 0), (2, 1), (1, 1), (0, 1)]
    vertical = {frozenset(((x, 0), (x, 1))) for x in range(3)}
    horizontal = {frozenset((a, b)) for a in coords for b in coords if a[1] == b[1] and abs(a[0]-b[0]) == 1}
    order = list(range(6))
    n = len(order)
    name = {i: f'{c[0]},{c[1]}' for i, c in enumerate(coords)}
    wires(ax, n, 0, n+2.2, [f'({name[m]})' for m in order])
    # Horizontal hops first: adjacent in the initial order.
    for q in range(n-1):
        if frozenset((coords[order[q]], coords[order[q+1]])) in horizontal:
            box(ax, .45 + .75*(q % 2), q, q+1, 'hop', SHORT, .5)
    hops = 0
    for layer in range(n):
        x = 2.0 + layer
        for q in range(layer % 2, n-1, 2):
            a, b = order[q], order[q+1]
            is_hop = frozenset((coords[a], coords[b])) in vertical
            box(ax, x, q, q+1, 'hop\n+fS' if is_hop else 'fS', HOP if is_hop else '#b8c1cc', .62)
            hops += is_hop
            order[q], order[q+1] = b, a
    for q in range(n):
        ax.text(n+2.35, -q, f'({name[order[q]]})', ha='left', va='center', fontsize=9)
    ax.text((n+2)/2, -6.25, f'odd–even fSWAP layers reverse the order: every pair of modes is adjacent once,\n'
            f'so all {hops} vertical hops act on neighbouring qubits; no long Z strings.', ha='center',
            fontsize=8.5, color=MUTED)
    ax.set_xlim(-1.2, n+3.4)
    ax.set_ylim(-7.0, .8)
    ax.axis('off')
    ax.set_title('(c) Fermionic-swap network, 2×3 strip', loc='left', fontsize=11)


def main():
    fig = plt.figure(figsize=(15.5, 5.2))
    grid = fig.add_gridspec(1, 3, width_ratios=[.8, 1.25, 1.25], wspace=.12)
    snake_panel(fig.add_subplot(grid[0]))
    direct_panel(fig.add_subplot(grid[1]))
    network_panel(fig.add_subplot(grid[2]))
    fig.text(.5, .01, 'fS: fermionic SWAP (SWAP with a −1 phase on |11⟩), which exchanges two adjacent modes and keeps the '
             'Jordan–Wigner order consistent. The physical arm uses such networks on the full lattice (Kivlichan et al.),\n'
             'so each Trotter step needs O(N·L) nearest-neighbour two-qubit gates; hence the physical arm\'s G ∝ rL³ gate count.',
             ha='center', fontsize=8.5, color=MUTED)
    for ext in ('pdf', 'png'):
        fig.savefig(ROOT/f'figures/fswap_network.{ext}', dpi=180, bbox_inches='tight')
    plt.close(fig)


if __name__ == '__main__':
    main()
