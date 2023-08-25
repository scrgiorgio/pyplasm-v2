Links:
- https://github.com/scrgiorgio/pyplasm-v2
- https://docs.google.com/document/d/1e_wwyUngGepHJ_7Xf-VKPADNJEOtebo2Y-mP-YXp0B0/edit?hl=it#heading=h.o51y9rmfsi0r

Translators:
- https://github.com/GunnarFarneback/CrudePythonTranslator.jl  (NOT GOOD)
- https://www.codeconvert.ai/python-to-julia-converter (VERY GOOD, spent 30$)


# Examples

Dependencies:
- numpy for MatrixNd
- glfw, OpenGL for the viewer
- scipy.spatial for qhull

```
set PYTHONPATH=src
python src/pyplasm/viewer.py
python src/pyplasm/hpc.py
python src/pyplasm/fenvs.py



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