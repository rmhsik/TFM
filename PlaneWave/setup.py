from setuptools import setup, Extension
from Cython.Build import cythonize

ext_modules = [
    Extension(
        "MatMult",
        ["MatMult.pyx"],
        extra_compile_args=['-fopenmp'],
        extra_link_args=['-fopenmp'],
    )
]
setup(
    name='TFM utils',
    ext_modules=cythonize(ext_modules),
    zip_safe=False,
)

