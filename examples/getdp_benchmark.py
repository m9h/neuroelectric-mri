"""Cross-backend CEM verification release."""

def write_mesh(data, path):
    raise NotImplementedError("Export tagged linear geometry without changing it")

def compare_solvers(mesh_path, sigma, *, degree=2, output_dir=None):
    raise NotImplementedError("Compare independently assembled CEM solutions")

def write_benchmark(output, *, plot=False):
    raise NotImplementedError("Archive independent verification evidence")
