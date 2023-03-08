# SPDX-FileCopyrightText: : 2023 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT

def memory(w):
    factor = 3.0
    for o in w.opts.split("-"):
        m = re.match(r"^(\d+)h$", o, re.IGNORECASE)
        if m is not None:
            factor /= int(m.group(1))
            break
    for o in w.opts.split("-"):
        m = re.match(r"^(\d+)seg$", o, re.IGNORECASE)
        if m is not None:
            factor *= int(m.group(1)) / 8760
            break
    if w.clusters.endswith("m"):
        return int(factor * (18000 + 180 * int(w.clusters[:-1])))
    elif w.clusters == "all":
        return int(factor * (18000 + 180 * 4000))
    else:
        return int(factor * (10000 + 195 * int(w.clusters)))


rule solve_network:
    input:
        RESOURCES + "networks/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}.nc",
    output:
        RESULTS + "networks/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}.nc",
    log:
        solver=normpath(
            LOGS + "solve_network/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_solver.log"
        ),
        python=LOGS
        + "solve_network/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_python.log",
        memory=LOGS
        + "solve_network/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_memory.log",
    benchmark:
        BENCHMARKS + "solve_network/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}"
    threads: 4
    resources:
        mem_mb=memory,
    shadow:
        "minimal"
    script:
        "scripts/solve_network.py"


rule solve_operations_network:
    input:
        unprepared=RESOURCES + "networks/elec_s{simpl}_{clusters}_ec.nc",
        optimized=RESULTS + "networks/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}.nc",
    output:
        RESULTS + "networks/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_op.nc",
    log:
        solver=normpath(
            LOGS
            + "solve_operations_network/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_op_solver.log"
        ),
        python=LOGS
        + "solve_operations_network/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_op_python.log",
        memory=LOGS
        + "solve_operations_network/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_op_memory.log",
    benchmark:
        (
            BENCHMARKS
            + "solve_operations_network/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}"
        )
    threads: 4
    resources:
        mem_mb=(lambda w: 5000 + 372 * int(w.clusters)),
    shadow:
        "minimal"
    script:
        "scripts/solve_operations_network.py"