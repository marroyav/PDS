#!/usr/bin/env python3
"""Reproduce the note's benchmarks; requires NumPy and SciPy.

Geometry: lossless direct light, cathode and representative wall apertures;
ND uses uniform first-hit wall flux. Lengths are metres; masses kilograms.
Activity is atmospheric Ar-39. Cosmic estimates include MIP tracks only.
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


def per_channel_estimates():
    """Print all four detectors in the same primary PE/s/channel units.

    Provenance: background_model.tex bibliography; per_channel_estimates.tex,
    cosmic_estimates.tex and nd_2x2_background.tex state the assumptions.
    XA PDE is whole-device; ND PDE is whole-panel, summed over its SiPMs.
    """
    specific_activity, density, photons_per_decay = 0.964, 1395.0, 4400.0
    xa_area, xa_pde, xa_channels = 0.36, 0.045, 2
    track_loss, muon_yield = 210.5, 20000.0  # MeV/m and photons/MeV
    surface_track_density = (4/3) * (1e4/60)  # cos^2 intensity; horizontal flux
    surf_track_density = 5.31e-5  # scalar flux approximation; see note
    vd_groups = (
        ("ProtoDUNE-VD cathode", 16, geometric_length(6.8, 3, 3.5)),
        ("ProtoDUNE-VD wall", 16, geometric_length(6.8, 3.5, 3)/2),
        ("FD-VD cathode", 640, geometric_length(13.5, 60, 6.5)),
        ("FD-VD long wall", 640,
         geometric_length(60, 6.5, 13.5, y=1.25)/2),
        ("FD-VD short wall", 64,
         geometric_length(13.5, 6.5, 60, y=1.25)/2),
    )
    print("\nPrimary scintillation PE/s per electronic channel:")
    print(f"{'Group':28s} {'Channels':>8s} {'Ar-39':>13s} {'Cosmic MIP':>13s}")
    for name, count, length in vd_groups:
        response_volume = xa_area * length * xa_pde / xa_channels
        ar_rate = specific_activity * density * photons_per_decay * response_volume
        track_density = (surface_track_density if name.startswith("Proto")
                         else surf_track_density)
        cosmic_rate = track_density * track_loss * muon_yield * response_volume
        print(f"{name:28s} {count:8d} {ar_rate:13.3f} {cosmic_rate:13.3f}")

    # Full-size CDR: 35 modules, 60 LCM and 20 ACL panels/module.
    full_acl_channels, full_lcm_channels = 35*20*6, 35*60*2
    # Two 0.5 x 1 x 3 m optical boxes/module, counting both cathode faces.
    interior_area = 2 * 2 * (0.5*1 + 0.5*3 + 1*3)
    full_h_acl = 20 * (0.30*0.50) / interior_area
    full_h_lcm = 60 * (0.10*0.50) / interior_area
    # Prototype: 32 ACL and 96 LCM, measured design coverage 29%.
    proto_acl_area, proto_lcm_area = 32*(0.28*0.30), 96*(0.10*0.30)
    proto_area = proto_acl_area + proto_lcm_area
    cases = (
        ("Full ND-LAr", 35*3*density, 70,
         full_acl_channels, full_lcm_channels, full_h_acl, full_h_lcm, 100, 3),
        ("ND 2x2", 2400.0, 8, 32*6, 96*2,
         0.29*proto_acl_area/proto_area, 0.29*proto_lcm_area/proto_area, 2, 1),
    )
    for name, mass, compartments, n_acl, n_lcm, h_acl, h_lcm, mu_rate, chord in cases:
        ar_photons = specific_activity * mass * photons_per_decay
        mu_photons = mu_rate * chord * track_loss * muon_yield
        ar_total, cosmic_total = 0.0, 0.0
        for technology, count, hit_fraction, pde in (
            ("ArCLight", n_acl, h_acl, 0.002), ("LCM", n_lcm, h_lcm, 0.006)
        ):
            response = hit_fraction * pde
            ar_per_channel = ar_photons * response / count
            # Independent compartment normalization must give the same result.
            compartment_rate = (specific_activity * (mass/compartments)
                                * photons_per_decay * response
                                / (count/compartments))
            assert isclose(ar_per_channel, compartment_rate, rel_tol=1e-12)
            cosmic_per_channel = mu_photons * response / count
            ar_total += ar_per_channel * count
            cosmic_total += cosmic_per_channel * count
            print(f"{name + ' ' + technology:28s} {count:8d} "
                  f"{ar_per_channel:13.3f} {cosmic_per_channel:13.3f}")
        channels = n_acl + n_lcm
        print(f"  {name} mean: Ar-39={ar_total/channels:.3f}; "
              f"cosmic={cosmic_total/channels:.3f}; "
              f"Ar-39 detector total={ar_total:.3f} PE/s")

    print("\nActivity bookkeeping (not a detected event rate):")
    for name, volume, channels in (
        ("ProtoDUNE-VD", 6.8*3*7, 32), ("FD-VD", 13.5*60*13, 1344),
        ("Full ND-LAr", 35*3, full_acl_channels+full_lcm_channels),
        ("ND 2x2", 2400/density, 384),
    ):
        mass = volume * density
        decays = specific_activity * mass
        print(f"{name}: mass={mass/1000:.3f} t; total={decays:.3f} Bq; "
              f"total/channels={decays/channels:.6f} Bq/channel")


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

    per_channel_estimates()
    print("PASS: solid-angle surface integration, infinite-plane limit, "
          "small-aperture limit, finite-aperture convergence, "
          "and equivalent ND compartment/channel normalization")


if __name__ == "__main__":
    main()
