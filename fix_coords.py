import json, os

nb_dir = "C:/Users/peter/OneDrive/Documents/Python/ThermalLandscapes/src/2_process_sites"
nbs = [f for f in os.listdir(nb_dir) if f.endswith(".ipynb")]

OLD = "x = pixels.x\n                    y = pixels.y"
NEW = (
    'yy, xx = np.meshgrid(pixels.y.values, pixels.x.values, indexing="ij")\n'
    '                    x = xr.DataArray(xx, dims=["y", "x"], coords={"x": pixels.x, "y": pixels.y})\n'
    '                    y = xr.DataArray(yy, dims=["y", "x"], coords={"x": pixels.x, "y": pixels.y})'
)

patched, skipped = [], []

for nb_name in nbs:
    nb_path = f"{nb_dir}/{nb_name}"
    with open(nb_path) as f:
        nb = json.load(f)

    found = False
    for cell in nb["cells"]:
        if cell["cell_type"] == "code":
            src = "".join(cell["source"])
            if OLD in src:
                cell["source"] = [src.replace(OLD, NEW)]
                found = True

    if found:
        with open(nb_path, "w") as f:
            json.dump(nb, f, indent=1)
        patched.append(nb_name)
    else:
        skipped.append(nb_name)

print("Patched:", patched)
print("Skipped (pattern not found):", skipped)
