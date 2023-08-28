Links:
- https://github.com/scrgiorgio/pyplasm-v2
- https://github.com/scrgiorgio/Plasm.jl
- https://docs.google.com/document/d/1e_wwyUngGepHJ_7Xf-VKPADNJEOtebo2Y-mP-YXp0B0/edit?hl=it#heading=h.o51y9rmfsi0r
- - https://julialang.org/contribute/developing_package/

Translators:
- https://github.com/GunnarFarneback/CrudePythonTranslator.jl  (NOT GOOD)
- https://www.codeconvert.ai/python-to-julia-converter (VERY GOOD, spent 30$)



# Plasm.jl

```
set PATH=%PATH%;c:\Julia-1.9.2\bin

julia

using Pkg
Pkg.add("Combinatorics)
Pkg.add("ModernGL")
Pkg.add("GLFW")
Pkg.add("StaticArrays")
Pkg.add(PackageSpec(name="PyCall", rev="master"))
Pkg.build("PyCall")
exit()

julia src/viewer.jl
julia src/hpc.jl
julia src/fenvs.jl

``````

# PyPlasm


```
set PYTHONPATH=.
python pyplasm/viewer.py
python pyplasm/hpc.py
python pyplasm/fenvs.py
python pyplasm/temple.py
python pyplasm/pisa.py
```