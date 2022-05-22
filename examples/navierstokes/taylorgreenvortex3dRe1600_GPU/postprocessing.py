#! /usr/bin/env python3
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2022 Pi-Yueh Chuang <pychuang@pm.me>
#
# Distributed under terms of the BSD 3-Clause license.

"""Post-processing.
"""
import multiprocessing
import pathlib
import re
import struct
import yaml
import numpy
import h5py
from matplotlib import pyplot


def interpolate_vorticity_to_center(wx, wy, wz):
    """Interpolate edge-centerd vorticit fields to cell center.
    """

    assert wx.shape[2] == wy.shape[2] - 1
    assert wx.shape[2] == wz.shape[2] - 1
    assert wy.shape[1] == wx.shape[1] - 1
    assert wy.shape[1] == wz.shape[1] - 1
    assert wz.shape[0] == wx.shape[0] - 1
    assert wz.shape[0] == wy.shape[0] - 1

    wxc = (wx[:-1, :-1, :] + wx[:-1, 1:, :] + wx[1:, :-1, :] + wx[1:, 1:, :]) / 4.
    wyc = (wy[:-1, :, :-1] + wy[:-1, :, 1:] + wy[1:, :, :-1] + wy[1:, :, 1:]) / 4.
    wzc = (wz[:, :-1, :-1] + wz[:, :-1, 1:] + wz[:, 1:, :-1] + wz[:, 1:, 1:]) / 4.

    assert wxc.shape == wyc.shape
    assert wyc.shape == wzc.shape

    return wxc, wyc, wzc


def plot_energy_properties(workdir, snapshots, nu, rho):
    """Plot kinetic enerty v.s. time
    """

    # read in baseline data
    baseline = numpy.loadtxt(workdir.joinpath("resources", "spectral_Re1600_512.gdiag"))

    simtimes = []
    simke = []
    simkedr = []
    for time, vals in snapshots.items():
        simtimes.append(time)
        simke.append(vals["ke"])
        simkedr.append(vals["kedr"])

    workdir.joinpath("figures").mkdir(exist_ok=True)

    pyplot.figure(figsize=(8, 6), dpi=166, constrained_layout=True)
    pyplot.plot(baseline[:, 0], baseline[:, 1], "k-", lw=2, label="Baseline")
    pyplot.scatter(simtimes, simke, 80, "tab:red", marker="x", label="Simulation")
    pyplot.xlabel("Time (seconds)")
    pyplot.ylabel(r"Volumetric kinetic energy ($kg / s^2 - m$)")
    pyplot.title("Mean (volumetric) kinetic energy v.s. time")
    pyplot.legend(loc=0)
    pyplot.savefig(workdir.joinpath("figures", f"mean-kinetic-energy.png"), dpi=166)
    pyplot.close("all")

    pyplot.figure(figsize=(8, 6), dpi=166, constrained_layout=True)
    pyplot.plot(baseline[:, 0], baseline[:, 2], "k-", lw=2, label="Baseline")
    pyplot.scatter(simtimes, simkedr, 80, "tab:red", marker="x", label="Simulation")
    pyplot.xlabel("Time (seconds)")
    pyplot.ylabel(r"Volumetric kinetic energy dissipation rate ($kg / s^3 - m$)")
    pyplot.title("Mean (volumetric) kinetic energy dissipation rate v.s. time")
    pyplot.legend(loc=0)
    pyplot.savefig(workdir.joinpath("figures", f"mean-kinetic-energy-dissipation-rate.png"), dpi=166)
    pyplot.close("all")


def plot_slice_contour(workdir):
    """Plot vorticity magnitude at t=8 on x=-pi, y=[0, pi/2], z=[pi/2, pi].
    """

    # read-in baseline data
    with open(root.joinpath("resources", "wn_slice_x0_08000.out"), "rb") as fobj:
        baseline = fobj.read()

    # make sure the byte size is correct
    assert len(baseline) == 8 * 512 * 512  # baseline data is 512x512

    # read all elements as a flatten 1D list (row dominate; as in Fortran)
    baseline = struct.unpack(f"{512*512}d", baseline)

    # convert to 2D ndarray and use C-style ordering
    baseline = numpy.array(baseline, dtype=float).reshape(512, 512).transpose()

    # add the missing nodes due to periodic BCs
    baseline = numpy.concatenate((baseline, baseline[:, :1]), axis=1)
    baseline = numpy.concatenate((baseline, baseline[:1, :]), axis=0)

    # gridlines for baseline (seems using different domain)
    yb = numpy.linspace(0, 2*numpy.pi, 513)
    zb = numpy.linspace(-numpy.pi, numpy.pi, 513)

    # read-in our simulation data
    with h5py.File(root.joinpath("output", "0000800.h5"), "r") as dset:
        time = dset["p"].attrs["time"]
        pred = {key: dset[key][...] for key in ["wx", "wy", "wz"]}

    # using periodic BC to interpolate to cell faces' centers on x=-pi
    wx, wy, wz = pred["wx"], pred["wy"], pred["wz"]
    wx0 = (
        ( # wx has no values at the face of x=-pi, interpolated from neighbors (using periodic BC)
            wx[:-1, :-1, 0] + wx[:-1, 1:, 0] + wx[1:, :-1, 0] + wx[1:, 1:, 0] +
            wx[:-1, :-1, -1] + wx[:-1, 1:, -1] + wx[1:, :-1, -1] + wx[1:, 1:, -1]
        )/ 8 +
        (wy[:-1, :, 0] + wy[1:, :, 0]) / 2. +  # wy has values on x=-pi, just they are on the edges
        (wz[:, :-1, 0] + wz[:, 1:, 0]) / 2  # wz has values on x=-pi, just they are on the edges
    )

    # shapes
    nz, ny = wx0.shape

    # gridlines
    y = numpy.linspace(-numpy.pi, numpy.pi, ny+1)
    z = numpy.linspace(-numpy.pi, numpy.pi, nz+1)
    y = (y[1:] + y[:-1]) / 2.
    z = (z[1:] + z[:-1]) / 2.

    pyplot.figure(figsize=(8, 6), dpi=166, constrained_layout=True)
    pyplot.contour(
        yb[:513//4+1], zb[3*513//4:], baseline[3*513//4:, :513//4+1],
        levels=[1, 5, 10, 20, 30], colors="k"
    )
    pyplot.contour(
        y[ny//2:3*ny//4], z[-nz//4:], wx0[-nz//4:, ny//2:(3*ny)//4], levels=[1, 5, 10, 20, 30],
        colors="tab:red", alpha=0.75
    )
    pyplot.xlim(0, numpy.pi/2.)
    pyplot.ylim(numpy.pi/2., 2.6)
    pyplot.xlabel("y")
    pyplot.ylabel("z")
    pyplot.title(r"Vorticity norm at $x=-\pi L$ 2")
    pyplot.savefig(root.joinpath("figures", f"vorticity-x=-piL.png"), dpi=166)
    pyplot.close("all")


def calculate_mean_kinetic_energy(u, v, w):
    """Simle linear interpolation + summation to calculate volumetric kinetic energy.
    """
    uintrp = (u[:, :, 1:] + u[:, :, :-1]) / 2.
    vintrp = (v[:, 1:, :] + v[:, :-1, :]) / 2.
    wintrp = (w[1:, :, :] + w[:-1, :, :]) / 2.

    nz, ny, nx = uintrp.shape
    mke = numpy.sum(uintrp**2 + vintrp**2 + wintrp**2) / ny / nx / nz / 2.
    return mke


def calculate_mean_enstropy(wx, wy, wz):
    """Simle linear interpolation + summation to calculate mean enstropy.
    """

    wxc, wyc, wzc = interpolate_vorticity_to_center(wx, wy, wz)
    nz, ny, nx = wxc.shape
    ens = numpy.sum(wxc**2 + wyc**2 + wzc**2) / ny / nx / nz / 2.
    return ens


def worker(rank, inputs, outputs, nu, fields):
    """A worker that processes snatshot data.
    """

    while True:
        try:
            filename = inputs.get(True, 2)
        except multiprocessing.queues.Empty:
            print(f"[Rank {rank}] no more work to do, going home")
            inputs.close()
            outputs.close()
            break

        print(f"[Rank {rank}] processing {filename.name}.")

        with h5py.File(filename, "r") as dset:
            time = dset["p"].attrs["time"]
            pred = {key: dset[key][...] for key in fields}

        # padding the ignored nodes due to periodic BCs (arrays' orders are k, j, i)
        pred["u"] = numpy.concatenate((pred["u"][:, :, -1:], pred["u"]), axis=2)
        pred["v"] = numpy.concatenate((pred["v"][:, -1:, :], pred["v"]), axis=1)
        pred["w"] = numpy.concatenate((pred["w"][-1:, :, :], pred["w"]), axis=0)

        # mean (volumetric) kinetic energy
        ke = calculate_mean_kinetic_energy(pred["u"], pred["v"], pred["w"])

        # mean (volumetric) enstropy
        ens = calculate_mean_enstropy(pred["wx"], pred["wy"], pred["wz"])

        # kinetic energy dissipation rate
        kedr = 2. * ens * nu

        outputs.put({"time": time, "pred": pred, "ke": ke, "ens": ens, "kedr": kedr})
        print(f"[Rank {rank}] done processing {filename.name}.")

        inputs.task_done()  # required by JoinableQueue


if __name__ == "__main__":

    root = pathlib.Path(__file__).resolve().parent

    # target fields
    fields = ["u", "v", "w", "wx", "wy", "wz", "p"]

    # read in configuration
    with open(root.joinpath("config.yaml"), "r") as fobj:
        cfg = yaml.load(fobj, Loader=yaml.CLoader)

    # gather all existing solution files
    filenames = [
        filename
        for filename in root.joinpath("output").glob("*.h5")
        if re.search(r"\d+?\.h5", filename.name) is not None
    ]

    # prepare work loads
    works = multiprocessing.JoinableQueue()
    for filename in filenames:
        works.put(filename)

    # start dynamic work distribution
    outcomes = multiprocessing.Queue()
    procs = []
    for rank in range(multiprocessing.cpu_count()//2):
        proc = multiprocessing.Process(
            target=worker, args=(rank, works, outcomes, cfg["flow"]["nu"], fields)
        )
        proc.start()
        procs.append(proc)

    # block until all works in `works` are done (i.e., `works` is empty)
    works.join()

    # get result from outcome queue to EMPTY the queue
    # workers won't terminate unitil involved queues (works & outcomes) are empty
    results = {}
    while len(results) != len(filenames):
        result = outcomes.get(block=True)  # won't return until it actually gets something
        time = result.pop("time")
        results[time] = result
    else:
        assert outcomes.empty()

    # plot kinetic energy
    plot_energy_properties(root, results, cfg["flow"]["nu"], 1.0)

    # plot contourf
    plot_slice_contour(root)
