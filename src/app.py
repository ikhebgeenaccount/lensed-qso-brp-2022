import configparser


class App:

    __conf = None

    @staticmethod
    def config():
        if App.__conf is None:  # Read only once, lazy.
            App.__conf = configparser.ConfigParser(interpolation=configparser.ExtendedInterpolation())
            App.__conf.read('settings.ini')
        return App.__conf
