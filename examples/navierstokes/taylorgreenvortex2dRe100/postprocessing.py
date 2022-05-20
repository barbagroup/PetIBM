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
import yaml
import numpy
import h5py
from matplotlib import pyplot


def get_exact_solution(x, y, t, nu, rho, field):
    """Exact solution of 2D Taylor-Green vortex.
    """
    if field == "u":
        return numpy.cos(x) * numpy.sin(y) * numpy.exp(-2.*nu*t)
    elif field == "v":
        return - numpy.sin(x) * numpy.cos(y) * numpy.exp(-2.*nu*t)
    elif field == "p":
        return - rho * numpy.exp(-4.*nu*t) * (numpy.cos(2.*x) + numpy.cos(2.*y)) / 4.
    elif field == "wz":
        return - 2. * numpy.cos(x) * numpy.cos(y) * numpy.exp(-2.*nu*t)
    elif field == "KE":  # kinetic energy
        return rho * numpy.pi * numpy.pi * numpy.exp(-4.*nu*t)
    elif field == "KEDR":  # kinetic energy dissipation rate
        return - 4. * nu * rho * numpy.pi * numpy.pi * numpy.exp(-4.*nu*t)
    else:
        raise ValueError


def calculate_ke(u, v):
    """Simle linear interpolation + summation to calculate total kinetic energy.
    """
    uintrp = (u[:, 1:] + u[:, :-1]) / 2.
    vintrp = (v[1:, :] + v[:-1, :]) / 2.

    ny, nx = uintrp.shape
    area = 4 * numpy.pi * numpy.pi
    ke = area * numpy.sum(uintrp**2 + vintrp**2) / ny / nx / 2.

    return ke


def plot_kinetic_energy(workdir, snapshots, nu, rho):
    """Plot kinetic enerty v.s. time
    """

    simtimes = []
    simke = []
    for time, vals in snapshots.items():
        simtimes.append(time)
        simke.append(calculate_ke(vals["pred"]["u"], vals["pred"]["v"]))

    times = numpy.linspace(0., 100., 101)
    ke = get_exact_solution(0, 0, times, nu, rho, "KE")  # x, y not important for calc. KE
    kedr = get_exact_solution(0, 0, times, nu, rho, "KEDR")  # x, y not important for calc. KEDR

    workdir.joinpath("figures").mkdir(exist_ok=True)

    pyplot.figure(figsize=(8, 6), dpi=166, constrained_layout=True)
    pyplot.plot(times, ke, "k-", lw=2, label="Analytical solution")
    pyplot.scatter(simtimes, simke, 80, "tab:red", marker="x", label="Simulation")
    pyplot.xlabel("Time (seconds)")
    pyplot.ylabel(r"Kinetic energy ($kg-m^2 / s^2$)")
    pyplot.title("Total kinetic energy v.s. time")
    pyplot.legend(loc=0)
    pyplot.savefig(workdir.joinpath("figures", f"kinetic-energy.png"), dpi=166)
    pyplot.close("all")

    pyplot.figure(figsize=(8, 6), dpi=166, constrained_layout=True)
    pyplot.plot(times, kedr, "k-", lw=3, label="Analytical solution")
    pyplot.xlabel("Time (seconds)")
    pyplot.ylabel(r"Kinetic energy dissipation rate ($kg-m^2 / s^3$)")
    pyplot.title("Kinetic energy dissipation rate v.s. time")
    pyplot.legend(loc=0)
    pyplot.savefig(workdir.joinpath("figures", f"kinetic-energy-dissipation-rate.png"), dpi=166)
    pyplot.close("all")


def contourf(workdir, x, y, vals, title, fname):
    """Generate a contourf
    """

    assert x.shape == vals.shape, f"{x.shape}, {vals.shape}"
    assert y.shape == vals.shape, f"{x.shape}, {vals.shape}"

    workdir.joinpath("figures").mkdir(exist_ok=True)

    pyplot.figure(figsize=(8, 6), dpi=166, constrained_layout=True)
    pyplot.contourf(x, y, vals, 128)
    pyplot.colorbar()
    pyplot.xlabel("x")
    pyplot.ylabel("y")
    pyplot.title(title)
    pyplot.savefig(workdir.joinpath("figures", f"{fname}.png"), dpi=166)
    pyplot.close("all")


def process_single_file(workdir, filename, x, y, nu, fields):
    """Process a single file.
    """

    with h5py.File(filename, "r") as dset:
        time = dset["p"].attrs["time"]
        pred = {key: dset[key][...] for key in fields}

    # padding the ignored nodes due to periodic BCs
    pred["u"] = numpy.concatenate((pred["u"][:, -1:], pred["u"]), axis=1)
    pred["v"] = numpy.concatenate((pred["v"][-1:, :], pred["v"]), axis=0)

    ans = {
        key: get_exact_solution(x[key], y[key], time, cfg["flow"]["nu"], 1.0, key)
        for key in fields
    }

    err = {key: abs(pred[key]-ans[key]) for key in fields}
    l1norm = {key: err[key].sum() / err[key].size for key in fields}
    l2norm = {key: numpy.sqrt((err[key]**2).sum()/err[key].size) for key in fields}

    for key in fields:
        contourf(
            workdir, x[key], y[key], ans[key],
            f"analytical solution; t={round(time, 2)}; field {key}",
            f"analyticalsolution-t{round(time, 2)}-{key}"
        )
        contourf(
            workdir, x[key], y[key], pred[key],
            f"simulation result; t={round(time, 2)}; field {key}",
            f"simulationresult-t{round(time, 2)}-{key}"
        )
        contourf(
            workdir, x[key], y[key], err[key],
            f"absolute error; t={round(time, 2)}; field {key}",
            f"absoluteerror-t{round(time, 2)}-{key}"
        )

    return {
        "time": time, "ans": ans, "pred": pred,
        "err": err, "l1norm": l1norm, "l2norm": l2norm
    }


def worker(rank, inputs, outputs, kwargs):
    """A worker that processes snatshot data.
    """

    while True:
        try:
            filename = inputs.get(False)  # w/ False, it returns even not getting anything
        except multiprocessing.queues.Empty:
            print(f"[Rank {rank}] has no more work to do.")
            inputs.close()
            outputs.close()
            return

        print(f"[Rank {rank}] processing {filename.name}.")
        outcome = process_single_file(filename=filename, **kwargs)
        outputs.put(outcome)
        print(f"[Rank {rank}] done processing {filename.name}.")


if __name__ == "__main__":

    root = pathlib.Path(__file__).resolve().parent

    # target fields
    fields = ["u", "v", "p"]

    # read in configuration
    with open(root.joinpath("config.yaml"), "r") as fobj:
        cfg = yaml.load(fobj, Loader=yaml.CLoader)

    # read in gridlines
    with h5py.File(root.joinpath("output", "grid.h5"), "r") as dset:
        x = {key: dset[f"{key}/x"][...] for key in fields}
        y = {key: dset[f"{key}/y"][...] for key in fields}

    # padding the nodes not used in simulation due to periodic BC
    x["u"] = numpy.concatenate((numpy.zeros(1, float), x["u"]))
    y["v"] = numpy.concatenate((numpy.zeros(1, float), y["v"]))

    # convert gridlines to mesh grids
    xy = {key: numpy.meshgrid(x[key], y[key]) for key in fields}
    x = {key: val[0] for key, val in xy.items()}
    y = {key: val[1] for key, val in xy.items()}

    # gather all existing solution files
    filenames = [
        filename
        for filename in root.joinpath("output").glob("*.h5")
        if re.search(r"\d+?\.h5", filename.name) is not None
    ]

    # prepare work loads
    works = multiprocessing.Queue()
    for filename in filenames:
        works.put(filename)

    # start dynamic work distribution
    outcomes = multiprocessing.Queue()
    procs = []
    for rank in range(multiprocessing.cpu_count()//2):
        proc = multiprocessing.Process(
            target=worker,
            args=(
                rank, works, outcomes,
                {"workdir": root, "x": x, "y": y, "nu": cfg["flow"]["nu"], "fields": fields}
            )
        )
        proc.start()
        procs.append(proc)

    # get result from outcome queue to EMPTY the queue
    # workers won't terminate unitil involved queues (works & outcomes) are empty
    results = {}
    while len(results) != len(filenames):
        result = outcomes.get(block=True)  # won't return until it actually gets something
        time = result.pop("time")
        results[time] = result
    else:
        assert outcomes.empty()

    # block the main program until all workser are done and terminated
    for proc in procs:
        proc.join()

    # plot kinetic energy
    plot_kinetic_energy(root, results, cfg["flow"]["nu"], 1.0)
