import os
from setuptools import setup, Extension
from Cython.Build import cythonize
import platform

# Windows-specific GSL paths
if platform.system() == "Windows":
    vcpkg_root = os.environ.get('VCPKG_ROOT', r'C:\vcpkg')

    gsl_include_dirs = [
        os.path.join(vcpkg_root, 'installed', 'x64-windows', 'include')
    ]
    gsl_library_dirs = [
        os.path.join(vcpkg_root, 'installed', 'x64-windows', 'lib')
    ]

    # Find existing GSL installation
    include_dirs = ["nlbayes", "core/include"]
    library_dirs = []
    
    for path in gsl_include_dirs:
        if os.path.exists(path):
            include_dirs.append(path)
            break
    
    for path in gsl_library_dirs:
        if os.path.exists(path):
            library_dirs.append(path)
            break

    libraries = ["gsl", "gslcblas"]
    extra_compile_args = ["/std:c++14"]
else:
    include_dirs = ["core/include"]
    library_dirs = []
    libraries = ["gsl", "gslcblas", "m"]
    extra_compile_args = ["-std=c++17"]
    


ext_modules = [
    Extension(
        "nlbayes.ModelORNOR",
        sources=[
            "nlbayes/ModelORNOR.pyx",
            "core/src/Beta.cpp",
            "core/src/Dirichlet.cpp",
            "core/src/GraphBase.cpp",
            "core/src/GraphORNOR.cpp",
            "core/src/HNode.cpp",
            "core/src/HNodeORNOR.cpp",
            "core/src/HParentNode.cpp",
            "core/src/ModelBase.cpp",
            "core/src/ModelORNOR.cpp",
            "core/src/Multinomial.cpp",
            "core/src/NodeDictionary.cpp",
            "core/src/RVNode.cpp",
            "core/src/SNode.cpp",
            "core/src/TNode.cpp",
            "core/src/XNode.cpp",
            "core/src/YDataNode.cpp",
            "core/src/YNoiseNode.cpp"
        ],
        include_dirs=include_dirs,
        library_dirs=library_dirs,
        libraries=libraries,
        language="c++",
        extra_compile_args=extra_compile_args,
        depends=["core/include/*.h"],
    )
]

setup(
    ext_modules=cythonize(ext_modules),
)
