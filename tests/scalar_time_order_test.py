"""Production heat/transport time order with nonzero deferred spatial corrections."""
import argparse
import csv
import math
import os
from pathlib import Path
import subprocess
import tempfile

parser = argparse.ArgumentParser()
parser.add_argument("--solver", type=Path, required=True)
SOLVER = parser.parse_args().solver.resolve()
BASE = Path(tempfile.mkdtemp(prefix="babelsim-scalar-time-order-"))
ENV = dict(os.environ, TMPDIR="/tmp", OMPI_ALLOW_RUN_AS_ROOT="1", OMPI_ALLOW_RUN_AS_ROOT_CONFIRM="1")
PATCHES = ("left", "right", "bottom", "top", "front", "back")


def mesh(path, nx, ny, length, shear):
    def vertex(i, j, k):
        return i + (nx+1)*(j+(ny+1)*k)
    vertices = [(length*i/nx+shear*length*j/ny, length*j/ny, k)
                for k in range(2) for j in range(ny+1) for i in range(nx+1)]
    layout = ((0,4,7,3), (1,2,6,5), (0,1,5,4), (3,7,6,2), (0,3,2,1), (4,5,6,7))
    faces, lookup, centres = [], {}, []
    for j in range(ny):
        for i in range(nx):
            cell = i+nx*j
            corners = [vertex(i,j,0), vertex(i+1,j,0), vertex(i+1,j+1,0), vertex(i,j+1,0),
                       vertex(i,j,1), vertex(i+1,j,1), vertex(i+1,j+1,1), vertex(i,j+1,1)]
            centres.append(tuple(sum(vertices[v][d] for v in corners)/8 for d in range(3)))
            boundary = (i==0, i==nx-1, j==0, j==ny-1, True, True)
            for side, order in enumerate(layout):
                ring = [corners[q] for q in order]
                key = tuple(sorted(ring))
                if key in lookup:
                    faces[lookup[key]][1] = cell
                else:
                    lookup[key] = len(faces)
                    faces.append([cell, -1, side if boundary[side] else -1, ring])
    text = ["BABELSIM_MESH 3", "vertices", str(len(vertices))]
    text += [" ".join(map(str, point)) for point in vertices]
    text += ["faces", str(len(faces))]
    text += [f"face 4 {o} {n} {p} " + " ".join(map(str, ring)) for o,n,p,ring in faces]
    text += ["patches", "6"]
    text += [f"patch\n{name}\ngeneric" for name in PATCHES]
    path.write_text("\n".join(text)+"\nend\n")
    return centres


def fixture(kind, steps, shortened):
    target = BASE/f"{kind}-{steps}-{shortened}"
    target.mkdir()
    (target/"fields").mkdir()
    heat = kind == "heat"
    name, equation = ("T", "temperature") if heat else ("C", "transport")
    centres = mesh(target/"mesh.mesh", 3 if heat else 8, 3 if heat else 1,
                   3 if heat else 1, 0.7 if heat else 0)
    (target/"case.bs").write_text(f"solver {kind}\nmesh mesh.mesh\nfields fields\nphysics physics.bs\n"
        "methods methods.bs\nsolution solution.bs\ncontrol control.bs\noutput output.bs\nghostLayers 3\n")
    (target/"physics.bs").write_text("density 1\nheatCapacity 1\nconductivity 1\nsource 0\n" if heat
        else "storage 1\ndiffusivity 0\nsource 0\n")
    (target/"methods.bs").write_text("time bdf2\n" + "\n".join(f"equation.{equation}.{key} {value}" for key,value in
        (("interpolation","corrected"),("gradient","leastSquares"),
         ("convection","upwind" if heat else "linearUpwind"),("diffusion","corrected" if heat else "orthogonal")))+"\n")
    (target/"solution.bs").write_text("\n".join(f"equation.{equation}.{key} {value}" for key,value in
        (("kspType","bcgs"),("pcType","bjacobi"),("absoluteTolerance","1e-14"),
         ("relativeTolerance","1e-13"),("maxIterations","2000")))+"\n")
    dt = 0.2/(steps+(0.5 if shortened else 0))
    (target/"control.bs").write_text(f"startTime 0\nendTime 0.2\ndeltaT {dt:.17g}\n")
    (target/"output.bs").write_text("directory results\ntimeName final\nwriteInterval 100000\n")
    boundary = "\n".join(f"{patch} {{ type " + ("fixedValue value (0)" if heat or index==0 else "zeroGradient") + " }"
                         for index,patch in enumerate(PATCHES))
    (target/f"fields/{name}.field").write_text(f"field {name}\n{{\ntype scalar\nlocation cell\ninternal file initial.dat\n"
        f"boundary\n{{\n{boundary}\n}}\n}}\n")
    (target/"fields/initial.dat").write_text("\n".join(f"{i} {1+0.1*x*y if heat else math.sin(math.pi*x):.17g}"
        for i,(x,y,z) in enumerate(centres))+"\n")
    if not heat:
        boundaries = "\n".join(f"{p} {{ type zeroGradient }}" for p in PATCHES)
        (target/"fields/U.field").write_text("field U\n{\ntype vector\nlocation cell\ninternal uniform (1 0 0)\n"
            f"boundary\n{{\n{boundaries}\n}}\n}}\n")
    return target, name


def run(target, field, ranks):
    label = f"np{ranks}"
    result = subprocess.run(["mpirun", "-np", str(ranks), str(SOLVER), "-case", str(target), "-time", label],
        env=ENV, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=90)
    (target/f"{label}.log").write_text(result.stdout)
    assert result.returncode == 0, (target, ranks, result.stdout[-3000:])
    values = {}
    for path in (target/"results"/label).glob(f"rank-*/{field}.csv"):
        for row in csv.DictReader(path.open()):
            index = int(row["global_id"])
            assert index not in values
            values[index] = float(row["value0"])
    assert values, target
    return [values[i] for i in sorted(values)]


print("Scalar time-order evidence:", BASE, flush=True)
for kind in ("heat", "transport"):
    for shortened in (False, True):
        answers = []
        for steps in (20,40,80,160):
            target, field = fixture(kind, steps, shortened)
            answer = run(target, field, 1)
            answers.append(answer)
            if steps == 40:
                for ranks in (2,4):
                    parallel = run(target, field, ranks)
                    assert len(answer) == len(parallel)
                    assert max(abs(a-b) for a,b in zip(answer,parallel)) < 1e-9
        differences = [math.sqrt(sum((a-b)**2 for a,b in zip(answers[i],answers[i+1]))) for i in range(3)]
        ratios = [differences[i]/differences[i+1] for i in range(2)]
        assert min(ratios) > 3.4, (kind,shortened,differences,ratios)
        print(kind, "shortened="+str(shortened), "differences=", differences, "ratios=", ratios, flush=True)
print("scalar_time_order_test: production BDF2 deferred corrections, shortened last step, 1/2/4 ranks passed")
