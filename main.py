#main.py

from core import App
import sys

def main():
    args = sys.argv[1:]

    config_file = args[0]

    App(config_file).run()

if __name__ == "__main__":
    main()