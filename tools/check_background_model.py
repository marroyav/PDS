#!/usr/bin/env python3
"""Reproduce the note's benchmarks; requires NumPy and SciPy.

Geometry: lossless direct light, centered cathode aperture, equal source
boxes on both faces. Lengths are metres; activity is atmospheric Ar-39.
These checks verify the arithmetic/model implementation, not detector inputs.
"""

from math import atan2, isclose, pi, sqrt

import numpy as np
from scipy.integrate import dblquad, quad


def rectangle_solid_angle(xmin, xmax, ymin, ymax, distance):
    """Solid angle of an unobstructed rectangle for positive normal distance."""
    if distance <= 0 or xmin >= xmax or ymin >= ymax:
        raise ValueError("Require positive distance and ordered bounds")
    return sum(
        sx * sy * atan2(x * y, distance * sqrt(distance**2 + x**2 + y**2))
        for sx, x in ((-1, xmin), (1, xmax))
        for sy, y in ((-1, ymin), (1, ymax))
    )


def geometric_length(width, height, drift, x=0.0, y=0.0):
    """G_DS/A for a differential cathode aperture at (x,y)."""
    value, error = quad(
        lambda z: rectangle_solid_angle(
            -width / 2 - x, width / 2 - x,
            -height / 2 - y, height / 2 - y, z
        ),
        0.0, drift, epsabs=1e-10, epsrel=1e-10,
    )
    assert error < 1e-7
    return value / (2 * pi)


def finite_aperture_length(width, height, drift, side=0.6, order=12):
    """Average the differential result over a square aperture by quadrature.

    Interchanging the source-volume and aperture-area integrals gives the
    exact finite-aperture direct acceptance, subject to quadrature accuracy.
    The 1/4 normalizes Gauss-Legendre weights on the two [-1,1] intervals.
    """
    nodes, weights = np.polynomial.legendre.leggauss(order)
    return sum(
        wx * wy * geometric_length(width, height, drift, side*x/2, side*y/2)
        for x, wx in zip(nodes, weights)
        for y, wy in zip(nodes, weights)
    ) / 4


def check_solid_angle():
    # Independent direct surface integration checks signs, including off-axis.
    for bounds in ((-0.3, 0.3, -0.3, 0.3), (0.2, 0.8, -0.1, 0.5)):
        xmin, xmax, ymin, ymax = bounds
        distance = 0.7
        numeric, error = dblquad(
            lambda y, x: distance / (x*x + y*y + distance*distance)**1.5,
            xmin, xmax, lambda x: ymin, lambda x: ymax,
            epsabs=1e-10, epsrel=1e-10,
        )
        analytic = rectangle_solid_angle(*bounds, distance)
        assert error < 1e-7 and isclose(analytic, numeric, rel_tol=1e-10)
        assert 0 < analytic < 2*pi
    assert isclose(geometric_length(1e8, 1e8, 3.5), 3.5, rel_tol=1e-6)


def main():
    check_solid_angle()
    activity_density = 0.964 * 1395.0
    area = 0.6**2
    results = []
    for name, dimensions in (("ProtoDUNE-VD", (6.8, 3.0, 3.5)),
                             ("FD-VD", (13.5, 60.0, 6.5))):
        point = geometric_length(*dimensions)
        finite = finite_aperture_length(*dimensions)
        finer = finite_aperture_length(*dimensions, order=20)
        assert 0 < finite <= point <= dimensions[2]
        assert isclose(finite, finer, rel_tol=1e-9)
        assert isclose(finite_aperture_length(*dimensions, side=1e-4),
                       point, rel_tol=1e-8)
        results.append((point, finite))
        print(f"{name}: ell = {point:.9f} m; G = {area*point:.9f} m^3; "
              f"aG = {activity_density*area*point:.3f} /s")
        print(f"  Finite aperture: ell = {finite:.9f} m; "
              f"change = {100*(finite/point-1):.4f}%")
    print(f"FD/PD ratio: point = {results[1][0]/results[0][0]:.8f}; "
          f"finite aperture = {results[1][1]/results[0][1]:.8f}")

    mass_per_tpc = 2600.0 / 8
    decays_per_tpc = 0.964 * mass_per_tpc
    print(f"ND: active mass/TPC = {mass_per_tpc:g} kg; "
          f"active volume/TPC = {mass_per_tpc/1395:.9f} m^3; "
          f"decays/TPC = {decays_per_tpc:g} /s")
    print(f"Old envelope interpreted as fully active: "
          f"{4*0.67*0.67*1.8*1395/1000:.6f} t")
    for coverage in (0.25, 0.30):
        # Assumptions only: uniform first-hit flux, equal technology areas,
        # 4400 photons/decay, and the representative 0.2%/0.6% device PDEs.
        toy = decays_per_tpc * 4400 * coverage * (0.002 + 0.006)/2
        print(f"ND illustration, coverage={coverage:g}: "
              f"{toy:.3f} PE/s/TPC; {8*toy:.3f} PE/s total")
    print("PASS: solid-angle surface integration, infinite-plane limit, "
          "small-aperture limit, and finite-aperture convergence")


if __name__ == "__main__":
    main()
