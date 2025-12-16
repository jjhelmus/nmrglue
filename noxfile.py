import nox

PYTHON_NUMPY_MATRIX = {
    "3.9": ["1.24", "1.25", "1.26", "2.0"],
    "3.10": ["1.24", "1.25", "1.26", "2.0", "2.1", "2.2"],
    "3.11": ["1.24", "1.25", "1.26", "2.0", "2.1", "2.2", "2.3"],
    "3.12": ["1.26", "2.0", "2.1", "2.2", "2.3"],
    "3.13": ["1.26", "2.0", "2.1", "2.2", "2.3"],
    "3.14": ["2.0", "2.1", "2.2", "2.3"],
}

ALL_NUMPY_VERSIONS = tuple(set(sum(PYTHON_NUMPY_MATRIX.values(), [])))

DEPENDENCIES = ["pytest", "pluggy", "scipy"]

@nox.session(python=list(PYTHON_NUMPY_MATRIX.keys()))
@nox.parametrize("numpy_version", ALL_NUMPY_VERSIONS)
def tests(session: nox.Session, numpy_version):
    """
    Main test function

    Parameters
    ----------
    session : nox.Session
        current nox session
    numpy_version : str
        numpy version to test
        
    """
    py_version = session.python
    if numpy_version not in PYTHON_NUMPY_MATRIX[py_version]:
        session.skip(f"NumPy version '{numpy_version}' is not supported on Python {py_version}")
        return

    numpy_dep = f"numpy=={numpy_version}.*"
    session.install(".", numpy_dep, *DEPENDENCIES)
    session.run(
        "pytest",
        "--ignore=nmrglue/util/xcpy_test.py",
        *session.posargs
    )