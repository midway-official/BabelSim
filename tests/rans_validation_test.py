"""Model-equation verification, homogeneous decay, and complete SIMPLE integration.

The decay oracle is derived from the published transport equations independently
of the C++ discretization. These are equation verification tests, not a claim of
general wall-flow validation for high-Re models without a wall-function workflow.
"""
import csv
import json
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[1]
BASE = Path(tempfile.mkdtemp(prefix="babelsim-rans-validation-"))
ENV = dict(os.environ, TMPDIR="/tmp")


def run(case, ranks, label, expected=0, executable="babelsim-solve"):
    args = ["mpirun", "-np", str(ranks), str(ROOT / "build" / executable)]
    args += [str(case)] if executable == "rans_equations_test" else ["-case", str(case), "-time", label]
    result = subprocess.run(args, cwd=ROOT, env=ENV, text=True,
                            stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=240)
    (BASE / (case.name + "-" + label + ".log")).write_text(result.stdout)
    assert result.returncode == expected, (args, result.returncode, result.stdout[-4000:])
    return result.stdout


def field(name, kind, value, boundary="zeroGradient"):
    patches = ["hot", "cold", "lower", "upper", "front", "back"]
    entries = "\n".join(f"    {p} {{ type {boundary}" +
        (f" value ({value})" if boundary == "fixedValue" else "") + " }" for p in patches)
    return (f"field {name}\n{{\ntype {kind}\nlocation cell\ninternal uniform ({value})\n"
            f"boundary\n{{\n{entries}\n}}\n}}\n")


def case(label, model, method="euler", dt=0.01, end=0.1, fixed=False, clipping=False):
    target = BASE / label
    shutil.copytree(ROOT / "cases/heat", target, ignore=shutil.ignore_patterns("results", "post"))
    path = target / "case.bs"
    path.write_text(path.read_text().replace("solver heat", "solver " +
        ("simple" if method == "steady" else "transientSimple")) + "\nghostLayers 3\n")
    (target / "numerics/methods.bs").write_text(
        f"interpolation linear\ngradient leastSquares\nconvection upwind\ndiffusion orthogonal\ntime {method}\n")
    (target / "control.bs").write_text(f"startTime 0\nendTime {end}\ndeltaT {dt}\n")
    (target / "output.bs").write_text("directory results\ntimeName final\nwriteInterval 100000\n")
    (target / "numerics/solution.bs").write_text(
        "scalarSolver bicgstab ilut 1e-14 1e-12 2000\n"
        "vectorSolver bicgstab ilut 1e-14 1e-12 2000\n"
        f"maxIterations {40 if clipping else 3000}\n"
        "velocityTolerance 1e-9\ncontinuityTolerance 1e-9\n"
        "pressureCorrectionTolerance 1e-9\nmomentumTolerance 1e-9\n")
    physics = f"density 1\ndynamicViscosity 0.01\nturbulenceModel {model}\n"
    if model != "none": physics += "turbulenceRelaxation 0.7\nturbulenceTolerance 1e-10\n"
    if clipping: physics += "kMin 1\nepsilonMin 1\n"
    (target / "physics/thermal.bs").write_text(physics)
    for name, kind, value in [("U","vector","0 0 0"),("p","scalar","1"),
        ("k","scalar","1"),("omega","scalar","1"),("epsilon","scalar","1"),
        ("nuTilda","scalar","0.3"),("wallDistance","scalar","0.5")]:
        bc = "fixedValue" if name == "U" or (fixed and name not in ("p","wallDistance")) else "zeroGradient"
        (target / "fields/initial" / (name + ".field")).write_text(field(name,kind,value,bc))
    return target


def values(case, label, name):
    result = {}
    for path in sorted((case / "results" / label).glob("rank-*/" + name + ".csv")):
        for row in csv.DictReader(path.open()):
            key = int(row["global_id"])
            assert key not in result
            result[key] = float(row["value0"])
    assert len(result) == 32, (case, label, name, len(result))
    return result


def sa_rhs(nu):
    chi = nu/0.01
    fv1 = chi**3/(chi**3+7.1**3)
    fv2 = 1-chi/(1+chi*fv1)
    ft2 = 1.2*math.exp(-0.5*chi**2)
    st = max(nu*fv2/(0.41**2*0.5**2), 1e-30)
    r = max(0, min(nu/(st*0.41**2*0.5**2), 10))
    g = r+0.3*(r**6-r)
    fw = g*(65/(g**6+64))**(1/6)
    cw1 = 0.1355/0.41**2+(1+0.622)/(2/3)
    return 0.1355*(1-ft2)*st*nu - (cw1*fw-0.1355*ft2/0.41**2)*nu**2/0.5**2


def exact(model, end):
    if model == "kOmega":
        return {"omega": 1/(1+0.075*end), "k": (1+0.075*end)**(-0.09/0.075)}
    if model == "kEpsilon":
        k = (1+0.92*end)**(-1/0.92)
        return {"k": k, "epsilon": k**1.92}
    nu, steps = 0.3, 20000
    h = end/steps
    for _ in range(steps):
        a = sa_rhs(nu); b = sa_rhs(nu+h*a/2)
        c = sa_rhs(nu+h*b/2); d = sa_rhs(nu+h*c)
        nu += h*(a+2*b+2*c+d)/6
    return {"nuTilda": nu}


summary = {"evidence": str(BASE), "decay": [], "steady": []}
print("RANS evidence:", BASE, flush=True)
for model in ("SA", "kOmega", "kEpsilon"):
    fixture = case("coefficients-" + model, model, dt=0.001, end=0.001)
    for ranks in (1,2,4):
        run(fixture, ranks, f"coefficients-{model}-{ranks}", executable="rans_equations_test")
    reference = exact(model, 0.1)
    for method in ("euler", "bdf2"):
        errors = []
        for dt in (0.01, 0.005, 0.0025):
            fixture = case(f"decay-{model}-{method}-{dt}", model, method, dt)
            run(fixture, 1, "serial")
            error = max(abs(v-reference[name]) for name in reference
                        for v in values(fixture,"serial",name).values())
            errors.append(error)
            if dt == 0.005:
                for ranks in (2,4):
                    run(fixture,ranks,f"np{ranks}")
                    for name in reference:
                        a, b = values(fixture,"serial",name), values(fixture,f"np{ranks}",name)
                        assert max(abs(a[i]-b[i]) for i in a) < 1e-9
        ratios = [errors[i]/errors[i+1] for i in (0,1)]
        assert min(ratios) > (1.8 if method == "euler" else 3.2), (model,method,errors,ratios)
        row = dict(model=model,method=method,errors=errors,ratios=ratios)
        summary["decay"].append(row)
        print(json.dumps(row), flush=True)

for model in ("none","SA","kOmega","kEpsilon"):
    fixture = case("steady-"+model,model,"steady",fixed=True)
    for ranks in (1,2,4):
        text = run(fixture,ranks,f"np{ranks}")
        assert "converged=true" in text
        if ranks > 1:
            for name in ({"SA":["nuTilda","mut"],"kOmega":["k","omega","mut"],
                          "kEpsilon":["k","epsilon","mut"],"none":["p"]}[model]):
                a,b = values(fixture,"np1",name),values(fixture,f"np{ranks}",name)
                assert max(abs(a[i]-b[i]) for i in a) < 1e-8
    summary["steady"].append(model)

fixture = case("clipping-rejection","kEpsilon",dt=1,end=1,clipping=True)
for ranks in (1,2,4):
    text = run(fixture,ranks,f"np{ranks}",expected=2)
    assert "converged=true" not in text
    assert max(map(float,re.findall(r"rTurb=([0-9.eE+-]+)",text))) > 1e-4
    assert not (fixture / "results" / f"np{ranks}").exists()

(BASE / "summary.json").write_text(json.dumps(summary,indent=2))
print("rans_validation_test: published terms, temporal order, full SIMPLE 1/2/4 ranks and clipping rejection passed",flush=True)
