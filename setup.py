import setuptools

if __name__ == "__main__":
    setuptools.setup(
        name='janus',
        version="0.1.1",
        description='A Python library for adaptive QM/MM methods',
        author='Boyi Zhang',
        author_email='boyi.zhang25@uga.edu',
        url="https://github.com/bzhang/janus",
        license='BSD-3-Clause',
        packages=setuptools.find_packages(),
        install_requires=[
            'numpy>=1.7',
        ],
        extras_require={
            'docs': [
                'sphinx==1.2.3',  # autodoc was broken in 1.3.1
                'sphinxcontrib-napoleon',
                'sphinx_rtd_theme',
                'numpydoc',
            ],
            'tests': [
                'pytest',
                'pytest-cov',
                'pytest-pep8',
                'tox',
            ],
        },

        tests_require=[
            'pytest',
            'pytest-cov',
            'pytest-pep8',
            'tox',
        ],
      #  data_files=[('janus/default_input/', ['janus/default_input/qmmm.json',
      #                                         'janus/default_input/aqmmm.json',
      #                                         'janus/default_input/psi4.json',
      #                                         'janus/default_input/md.json',
      #                                         'janus/default_input/system.json',
      #                                         'janus/default_input/openmm.json'])],

        entry_points={
            'console_scripts': [
                'janus = janus.janus:main',
            ],
        },

        classifiers=[
            'Development Status :: 4 - Beta',
            'Intended Audience :: Science/Research',
            'Programming Language :: Python :: 2.7',
            'Programming Language :: Python :: 3',
        ],
    )
