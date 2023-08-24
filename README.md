Links:
- https://github.com/scrgiorgio/pyplasm-v2

# Examples

Dependencies:
- numpy for MatrixNd
- glfw, OpenGL for the viewer
- scipy.spatial for qhull

```
set PYTHONPATH=src
python src/pyplasm/hpc.py
python src/pyplasm/fenvs.py
python src/pyplasm/viewer.py


python examples/python/temple.py
python examples/python/pisa.py
```

# Julia

```
set PATH=%PATH%;c:\Julia-1.9.2\bin

julia

using Pkg
Pkg.add("ModernGL")
Pkg.add("GLFW")
Pkg.add("StaticArrays")
exit()

julia src/julia/Viewer.jl

``````