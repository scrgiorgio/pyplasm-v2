Links:
- https://github.com/scrgiorgio/pyplasm-v2

# Examples

```
set PYTHONPATH=src
python src/pyplasm/fenvs.py
python src/pyplasm/viewer.py


python examples/python/temple.py
python examples/python/pisa.py
```

Dependencies:
- numpy for some matrix things (can be eliminated?)
- glfw, OpenGL for the viewer
- scipy.spatial (alternatives?)
- 
- 
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