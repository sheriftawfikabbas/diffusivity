from setuptools import setup, find_packages
from pathlib import Path

here = Path(__file__).resolve().parent
README = (here / "README.md").read_text(encoding="utf-8")
VERSION = (here / 'diffusivity' / "VERSION").read_text(encoding="utf-8").strip()

setup(
    name='diffusivity',
    packages=['diffusivity',
              ] + find_packages(exclude=['tests', 'tests.*']),
    include_package_data=True,
    entry_points={
        "console_scripts": ["diffusivity=diffusivity.cli:execute_cli"],
    },
    version=VERSION,
    license='mit',
    description='diffusivity calculates the diffusion properties by processing trajectory files',
    long_description=README,
    long_description_content_type='text/markdown',
    author='Sherif Abdulkader Tawfik Abbas',
    author_email='sherif.tawfic@gmail.com',
    url='https://github.com/sheriftawfikabbas/diffusivity',
    keywords=['material science', 'simulation'],
    install_requires=['ase',
                      'pandas',
                      'numpy'],

)
