"""
Test Program

"""
import sys


def main():

    python = sys.version.split("\n")[0]

    try:
        name, curdir, curexpno, curprocno = sys.argv
        message = """
        Welcome to Python {}.
        The current directory is set to {}
        EXPNO {} is currently open at PROCNO {}
        """.format(python, curdir, curexpno, curprocno)

    except ValueError:
        message = """
        Welcome to Python {}.
        No directory is currently open, 
        or none was passed on to this script.
        """.format(python)

    print(message)


if __name__ == "__main__":
    main()
