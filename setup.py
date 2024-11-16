from setuptools import setup, Extension

module = Extension("symnmfmodule", sources=["symnmf.c", "symnmfmodule.c"])

setup(
    name="symnmfmodule",
    version="1.0",
    description="Python interface for the symnmf C extension",
    ext_modules=[module],
)
