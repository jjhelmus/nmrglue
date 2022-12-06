"""
Test Program

"""
import sys


def main():

    python = sys.version.split("\n")[0]

    try:
        name, curdir, curexpno, curprocno = sys.argv
        message = f"""
        Welcome to Python {python}.
        The current directory is set to {curdir}
        EXPNO {curexpno} is currently open at PROCNO {curprocno}
        """

    except ValueError:
        message = f"""
        Welcome to Python {python}.
        No directory is currently open,
        or none was passed on to this script.
        """

    print(message)


if __name__ == "__main__":
    main()
