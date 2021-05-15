#main.py

from core import app

if __name__ == "__main__":

    app = app("files/config.info")
    app.run()